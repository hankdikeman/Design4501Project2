"""
Implements an inheritance-based component substructure for reactor network,
couples components to inlets and outlets to other components. Determines outlet
compositions as a function of inlet compositions
"""
import numpy as np
from abc import ABCMeta, abstractmethod
from thermoutils import get_psat, get_tsat
from scipy.optimize import newton


# baseclass for component
class Component(metaclass=ABCMeta):
    def __init__(self, inlets, outlets, fixed=False):
        self.inlets = inlets
        self.outlets = outlets
        self.fixed = fixed

    # return inlet dictionary
    def get_inlets(self):
        return self.inlets

    # return calculated outlets (must be implemented in each component!)
    @abstractmethod
    def calc_outlets(self):
        pass

    # check whether stream is fixed
    def is_fixed(self):
        return self.fixed

    # check whether solution has converged for this component
    def check_solution(self, next_inlets):
        # check all inlet flows to be the same within tolerance
        for inflow in self.inlets.keys():
            previous = self.inlets[inflow]
            iterated = next_inlets[inflow]
            distance = np.linalg.norm(previous - iterated)/np.count_nonzero(previous)
            print('inflow:',inflow,', lin alg dist =',distance)
            if distance > 1E-4:
                return False
        return True


class Feed(Component):
    def __init__(self, outlets):
        super(Feed, self).__init__({}, outlets, fixed=True)

    def calc_outlets(self):
        return self.outlets


class Removal(Component):
    def __init__(self, inlets):
        super(Removal, self).__init__(inlets, {})

    def calc_outlets(self):
        return {}


class ProductRemoval(Removal):
    def __init__(self, inlets, fixed=True):
        super(Removal, self).__init__(inlets, {})


class Mixer(Component):
    def __init__(self, inlets, outlets):
        super(Mixer, self).__init__(inlets, outlets)

    def calc_outlets(self):
        mixed = np.zeros(len(list(self.inlets.values())[0]))
        for keys in self.inlets.keys():
            mixed += self.inlets[keys]
        self.outlets[next(iter(self.outlets))] = mixed
        return self.outlets


class Splitter(Component):
    def __init__(self, inlets, outlets, split_fraction, split_stream_name):
        self.split = split_fraction
        self.splitstream = split_stream_name
        super(Splitter, self).__init__(inlets, outlets)

    def calc_outlets(self):
        for stream in self.outlets.keys():
            if stream == self.splitstream:
                self.outlets[stream] = list(self.inlets.values())[0] \
                    * self.split
            else:
                self.outlets[stream] = list(self.inlets.values())[0] \
                    * (1 - self.split)
        return self.outlets


class HeatExchanger(Component):
    def __init__(self, inlettemp, outlettemp, inlets, outlets):
        self.in_temp = inlettemp
        self.out_temp = outlettemp
        # run super constructor
        super(HeatExchanger, self).__init__(inlets, outlets)

    def calc_outlets(self):
        return 0


class Compressor(Component):
    def __init__(self, in_t, out_t, in_p, out_p, inlets, outlets):
        self.in_temp = in_t
        self.out_temp = out_t
        self.in_press = in_p
        self.out_temp = out_p
        # run super constructor
        super(Compressor, self).__init__(inlets, outlets)

    def calc_outlets(self):
        return 0


class Turbine(Component):
    def __init__(self, in_t, out_t, in_p, out_p, inlets, outlets):
        self.in_temp = in_t
        self.out_temp = out_t
        self.in_press = in_p
        self.out_temp = out_p
        # run super constructor
        super(Turbine, self).__init__(inlets, outlets)

    def calc_outlets(self):
        return 0


# reactor baseclass
class Reactor(Component):
    def __init__(self, temperature, pressure, inlets, outlets, conv_function):
        self.temperature = temperature
        self.pressure = pressure
        self.conversion_function = conv_function
        # run super constructor
        super(Reactor, self).__init__(inlets, outlets)

    def calc_outlets(self):
        # calculate new outlet values and return
        self.outlets[list(self.outlets.keys())[0]] = self.conversion_function(
            list(self.inlets.values())[0], self.pressure, self.temperature)
        return self.outlets


class FlashTank(Component):
    def __init__(self, pressure, recov, bub_key, inlets, v_out, l_out):
        self.pressure = pressure
        self.xi_key, self.ind_key = recov['Key']
        self.vapor_out = list(v_out.keys())[0]
        self.liquid_out = list(l_out.keys())[0]
        self.bub_key = bub_key
        # run super constructor
        super(FlashTank, self).__init__(inlets, {**v_out, **l_out})

    def calc_outlets(self):
        if (np.sum(list(self.inlets.values())[0])):
            temp = newton(self.temperature_iteration, 300)
            alpha = get_psat(temp) / get_psat(temp)[self.ind_key]
            # Calculate recoveries
            xi = np.divide(alpha * self.xi_key, 1 + (alpha - 1) * self.xi_key)
            # Mass Balance
            self.outlets[self.vapor_out] = v = xi * list(self.inlets.values())[0]
            self.outlets[self.liquid_out] = l = (
                1 - xi) * list(self.inlets.values())[0]
            inletstuff = list(self.inlets.values())[0]
            # manually remove gas from outlet liquid stream
            gasinliq = np.copy(self.outlets[self.liquid_out][0:3])
            gasingas = np.copy(self.outlets[self.vapor_out][0:3])
            # add gas stream to the vapor out stream and set in liq to 0
            self.outlets[self.liquid_out][0:3] = 0
            self.outlets[self.vapor_out][0:3] = gasinliq + gasingas
            # manually remove liquid from gas outlet stream
            liqinliq = np.copy(self.outlets[self.liquid_out][3:5])
            liqingas = np.copy(self.outlets[self.vapor_out][3:5])
            # add liquid to liquid out stream and set in vap to 0
            self.outlets[self.vapor_out][3:5] = liqingas*0.005
            self.outlets[self.liquid_out][3:5] = liqinliq + liqingas*0.995
            return self.outlets
        return self.outlets

    def temperature_iteration(self, temp_guess):
        alpha = get_psat(temp_guess) / get_psat(temp_guess)[self.ind_key]
        # Calculate recoveries
        xi = np.divide(alpha * self.xi_key, 1 + (alpha - 1) * self.xi_key)
        # Mass Balance
        self.outlets[self.vapor_out] = v = xi * list(self.inlets.values())[0]
        self.outlets[self.liquid_out] = l = (
            1 - xi) * list(self.inlets.values())[0]
        x = l / sum(l)
        # Bubble point
        temp_comp = (get_tsat(self.pressure *
                              alpha[self.bub_key] / np.sum(alpha * x)))[self.bub_key]
        return temp_guess - temp_comp


class Adsorber(Component):
    def __init__(self, recov, recov_index, inlets, outlets, ads_out):
        self.outkey = list(outlets.keys())[0]
        self.inkey = list(inlets.keys())[0]
        self.adskey = list(ads_out.keys())[0]
        self.recov = recov
        self.recov_index = recov_index
        # run super constructor
        super(Adsorber, self).__init__(inlets, {**ads_out, **outlets})

    def calc_outlets(self):
        flows_in = self.inlets[self.inkey]
        # initialize streams to all zero
        self.outlets[self.adskey] = np.zeros_like(self.inlets[self.inkey], dtype=np.float64)
        self.outlets[self.outkey] = np.zeros_like(self.inlets[self.inkey], dtype=np.float64)
        # split non-water flows
        self.outlets[self.outkey] = 0.995 * (flows_in)
        self.outlets[self.adskey] = 0.005 * (flows_in)
        # split water flow
        self.outlets[self.outkey][self.recov_index] = self.inlets[self.inkey][self.recov_index]*0.01
        self.outlets[self.adskey][self.recov_index] = self.inlets[self.inkey][self.recov_index]*0.99
        # print out stream compositions
        return self.outlets

class DistillationColumn(Component):
    def __init__(self, temp_in, pressure, recov, inlets, v_out, l_out):
        self.temperature_feed = temp_in
        self.pressure = pressure
        self.xi_hk, self.ind_hk = recov['HK']
        self.xi_lk, self.ind_lk = recov['LK']
        # self.temperature_bottom
        # self.temperature_top
        self.vapor_out = list(v_out.keys())[0]
        self.liquid_out = list(l_out.keys())[0]
        self.alpha = get_psat(temp_in) / get_psat(temp_in)[self.ind_hk]
        # run super constructor
        super(DistillationColumn, self).__init__(inlets, {**v_out, **l_out})

    def calc_outlets(self):
        if (np.sum(list(self.inlets.values())[0])):
            # Compute minimum number of trays
            Nmin = np.log(self.xi_lk * (1 - self.xi_hk) / ((1 - self.xi_lk)
                                                           * self.xi_hk)) / np.log(self.alpha[self.ind_lk])
            # Compute split fractions of all components
            xi = np.power(self.alpha, Nmin) * self.xi_hk / (
                1 + (np.power(self.alpha, Nmin) - 1) * self.xi_hk)
            # compute distillate and vapor flowrates
            self.outlets[self.vapor_out] = d = xi * list(self.inlets.values())[0]
            self.outlets[self.liquid_out] = b = (
                1 - xi) * list(self.inlets.values())[0]
            B_hk = 1 - xi[self.ind_hk]
            # calculate and print flows through reboiler and condenser
            B_lk = xi[self.ind_lk]
            R_lk = 1.38 / (np.power(self.alpha[self.ind_lk] - 1,
                                    0.9) * np.power(1 - B_lk, 0.1))
            R_hk = 1.38 / (np.power(self.alpha[self.ind_lk] - 1,
                                    0.9) * np.power(1 - B_hk, 0.1))
            Ractual = (0.8 * np.amax([R_lk, R_hk]) + 0.2 * np.amin([R_lk, R_hk]))
            D = self.outlets[self.vapor_out]
            L = np.sum(D * Ractual + list(self.inlets.values())[0] - self.outlets[self.liquid_out])
            V = np.sum(D * (1 + Ractual))
            return self.outlets
        return self.outlets

    # Calculate key temperatures
    def calc_bot_temp(self):
        # Compute mole fraction
        x_B = self.outlets[self.liquid_out] / \
            np.sum(self.outlets[self.liquid_out])
        return get_tsat((self.pressure) * np.sum(x_B / self.alpha))

    def calc_top_temp(self):
        # Compute mole fraction
        x_D = self.outlets[self.vapor_out] / \
            np.sum(self.outlets[self.vapor_out])
        return get_tsat((self.pressure) * self.alpha[self.ind_lk] / np.sum(x_D * self.alpha))


class Absorber(Component):
    def __init__(self, solvent_t_in, pressure, key_recovery, key_index, solvent_index, v_in, v_out, l_in, l_out):
        self.solvent_temp_in = solvent_t_in
        self.pressure = pressure
        self.key_recovery = key_recovery
        self.key_index = key_index
        self.solvent_index = solvent_index
        self.vapor_in = list(v_in.keys())[0]
        self.vapor_out = list(v_out.keys())[0]
        self.liquid_in = list(l_in.keys())[0]
        self.liquid_out = list(l_out.keys())[0]
        # run super constructor
        super(Absorber, self).__init__({**v_in, **l_in}, {**v_out, **l_out})

    def calc_outlets(self):
        self.outlets[self.liquid_out][3] = self.inlets[self.vapor_in][3] + self.inlets[self.liquid_in][3]
        self.outlets[self.liquid_out][4] = self.inlets[self.vapor_in][4]
        self.outlets[self.liquid_out][5] = self.inlets[self.vapor_in][5]
        self.outlets[self.vapor_out][0] = self.inlets[self.vapor_in][0]
        self.outlets[self.vapor_out][1] = self.inlets[self.vapor_in][1]
        self.outlets[self.vapor_out][2] = self.inlets[self.vapor_in][2]
        self.outlets[self.vapor_out][6] = self.inlets[self.vapor_in][6]
        self.outlets[self.vapor_out][7] = self.inlets[self.vapor_in][7]
        return self.outlets


# generate streams with components specified with kwargs
def StreamGen(H2=0, CO2=0, CO=0, H2O=0, MEOH=0, FA=0, N2=0, O2=0, OME1=0, OME2=0, OME3=0, OME4=0, OME5=0, OME6=0):
    return 6.47*np.array([H2, CO2, CO, H2O, MEOH, FA, N2, O2, OME1, OME2, OME3, OME4, OME5, OME6], dtype=np.float64)


# unit tests for components
if __name__ == "__main__":
    print(__doc__)
    dict1 = {'mu1': np.array([1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0])}
    dict2 = {'mu2': np.array([0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13])}
    dict3 = {'mu3': np.array([0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13])}
    dict4 = {'mu4': np.array([0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13])}
    absorber1 = Absorber(300, 2, 0.97, 5, 3, dict1, dict2, dict3, dict4)
    recov = {'LK': (0.99, 1), 'HK': (0.001, 2)}
    distillationColumn = DistillationColumn(300, 2, recov, dict1, dict2, dict3)
    # print(distillationColumn.calc_outlets())

    recov1 = {'Key': (0.99, 4)}
    flash = FlashTank(760, recov1, 3, dict1, dict2, dict3)
    print(flash.calc_outlets())

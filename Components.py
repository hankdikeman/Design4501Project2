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

    def is_fixed(self):
        return self.fixed

    def check_solution(self, next_inlets):
        # check all inlet flows to be the same within tolerance
        for inflow in self.inlets.keys():
            previous = self.inlets[inflow]
            iterated = next_inlets[inflow]
            if np.linalg.norm(previous - iterated) > 1E-9:
                return False
        return True


class Feed(Component):
    def __init__(self, outlets):
        super(Feed, self).__init__({}, outlets)

    def calc_outlets(self):
        return self.outlets


class Removal(Component):
    def __init__(self, inlets):
        super(Removal, self).__init__(inlets, {})

    def calc_outlets(self):
        return {}


class ProductRemoval(Removal):
    def __init__(self, inlets, fixed=True):
        super(Removal, self).__init__(inlets, {}, fixed=fixed)


class Mixer(Component):
    def __init__(self, inlets, outlets):
        super(Mixer, self).__init__(inlets, outlets)

    def calc_outlets(self):
        mixed = np.zeros(len(list(self.inlets.values())[0]))
        for keys in self.inlets.keys():
            mixed += self.inlets[keys]
        self.outlets[next(iter(self.outlets))] = mixed


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
        list(self.outlets.values())[0] = self.conversion_function(
            list(self.inlets.values)[0], self.pressure, self.temperature)
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
        temp = newton(self.temperature_iteration, 300)
        alpha = get_psat(temp) / get_psat(temp)[self.ind_key]
        # Calculate recoveries
        xi = np.divide(alpha * self.xi_key, 1 + (alpha - 1) * self.xi_key)
        # Mass Balance
        self.outlets[self.vapor_out] = v = xi * list(self.inlets.values())[0]
        self.outlets[self.liquid_out] = l = (
            1 - xi) * list(self.inlets.values())[0]
        print("vout", self.outlets[self.vapor_out])
        print("lout", self.outlets[self.liquid_out])
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
        # Variables to track
        l_0 = np.zeros(len(list(self.inlets.values())[0]))
        AF = np.zeros(len(list(self.inlets.values())[0]))
        v_Np1 = self.inlets[self.vapor_in]
        n = self.key_index
        while (True):
            # Vapor pressure
            P0 = get_psat(self.solvent_temp_in)
            # Solvent flowrate
            AF[n] = 1.4  # Key component effective absorption factor
            # Also need to set air to zero (dont know what component)
            l_0[self.solvent_index] = AF[n] * \
                sum(v_Np1) * P0[n] / self.pressure  # solvent
            # Number of stages (deleted l_0[n]):
            N = np.log(((self.key_recovery - AF[n]) * v_Np1[n]) / (
                AF[n] * (1 - self.key_recovery) * v_Np1[n])) / np.log(AF[n])
            # Alpha
            alpha = P0 / P0[n]
            # Check to see if AF of solvent has been adjusted at the bottom yet
            if (AF[self.solvent_index] == 0):
                # Need someway to set absorption factors of remaining species (so minus key and air)
                AF[self.solvent_index] = AF[n] / alpha[self.solvent_index]
            # Beta values
            beta_N = (1 - AF**(N + 1)) / (1 - AF)
            beta_Nm1 = (1 - AF**N) / (1 - AF)
            # Mass balance
            v_1 = 1 / beta_N * v_Np1 + beta_Nm1 / beta_N * l_0
            l_N = v_Np1 + l_0 - v_1
            x_N = l_N / sum(l_N)
            # Temperature of solvent out (bubble point)
            alpha_avg = sum(x_N * alpha)
            # Need to fix Antoine coeff.
            T_N = get_tsat(np.log(self.pressure / alpha_avg)
                           )[self.solvent_index]
            # Check top vs bottom solvent temp
            if (abs(self.solvent_temp_in - T_N) > 100):
                # Absortion factors with bottom temp
                P0_N = get_psat(T_N)
                alpha_N_wat = P0_N[self.solvent_index] / P0_N[n]
                AF_N_wat = AF[n] / alpha_N_wat
                # Update absorption factors (Edmister equations)
                AF[self.solvent_index] = (
                    AF_N_wat * (1 + AF[self.solvent_index]) + 0.25)**0.5 - 0.5
            else:
                self.inlets[self.liquid_in] = l_0
                self.outlets[self.vapor_out] = v_1
                self.outlets[self.liquid_out] = l_N
                break
        # The liquid component flowrates Solvent (water), FORM and MeOH
        return self.outlets


def StreamGen(H2=0, CO2=0, CO=0, H2O=0, MEOH=0, FA=0, N2=0, O2=0, OME1=0, OME2=0, OME3=0, OME4=0, OME5=0, OME6=0):
    return np.array([H2, CO2, CO, H2O, MEOH, FA, N2, O2, OME1, OME2, OME3, OME4, OME5, OME6])


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

"""
Implements an inheritance-based component substructure for reactor network,
couples components to inlets and outlets to other components. Determines outlet
compositions as a function of inlet compositions
"""
import numpy as np
from abc import ABCMeta, abstractmethod
from thermoutils import get_psat, get_tsat


# baseclass for component
class Component(metaclass=ABCMeta):
    next_inlets = {}

    def __init__(self, inlets, outlets):
        self.inlets = inlets
        self.outlets = outlets

    # return inlet dictionary
    def get_inlets(self):
        return self.inlets

    # return calculated outlets (must be implemented in each component!)
    @abstractmethod
    def calc_outlets(self):
        pass

    def check_solution(self):
        # check all inlet flows to be the same within tolerance
        for inflow in self.inlets.keys():
            previous = self.inlets[inflow]
            iterated = self.next_inlets[inflow]
            if np.linalg.norm(previous - iterated) > 1E-9:
                return False
        return True


class HeatExchanger(Component):
    def __init__(self, inlettemp, outlettemp, inlets, outlets):
        self.in_temp = inlettemp
        self.out_temp = outlettemp
        # run super constructor
        super(HeatExchanger, self).__init__(inlets, outlets)


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


# reactor baseclass
class Reactor(Component):
    def __init__(self, temperature, pressure, inlets, outlets):
        self.temperature = temperature
        self.pressure = pressure
        # run super constructor
        super(Reactor, self).__init__(inlets, outlets)


class DistillationColumn(Component):
    def __init__(self, temp_in, pressure, recov, inlets, v_out, l_out):
        self.temperature_feed = temp_in
        self.pressure = pressure
        self.xi_hk, self.ind_hk = recov['HK']
        self.xi_lk, self.ind_lk = recov['LK']
        # self.temperature_bottom
        # self.temperature_top
        self.vapor_out = v_out.keys()[0]
        self.liquid_out = l_out.keys()[0]
        self.alpha = get_psat(temp_in) / get_psat(temp_in)[self.ind_hk]
        # run super constructor
        super(DistillationColumn, self).__init__(inlets, {**v_out, **l_out})

        def get_outlet():
            # Compute minimum number of trays
            Nmin = np.log(self.xi_lk * (1 - self.xi_hk) / ((1 - self.xi_lk)
                                                           * self.xi_hk)) / np.log(self.alpha[self.ind_lk])
            # Compute split fractions of all components
            xi = np.power(self.alpha, Nmin) * self.xi_hk / (
                1 + (np.power(self.alpha, Nmin) - 1) * self.xi_hk)
            # compute distillate and vapor flowrates
            self.outlets[self.vapor_out] = d = xi * self.inlets.values()[0]
            self.outlets[self.liquid_out] = b = (
                1 - xi) * self.inlets.values()[0]
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
        self.vapor_in = v_in.keys()[0]
        self.vapor_out = v_out.keys()[0]
        self.liquid_in = l_in.keys()[0]
        self.liquid_out = l_out.keys()[0]
        # run super constructor
        super(Absorber, self).__init__({**v_in, **l_in}, {**v_out, **l_out})

        def calc_outlets():
            # Variables to track
            l_0 = np.zeros(len(self.inlets.values()[0]))
            AF = np.zeros(len(self.inlets.values()[0]))
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
                T_N = get_tsat(np.log(self.pressure / alpha_avg))
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
                    self.outlets[self.vapor_out] = v_1
                    self.outlets[self.liquid_out] = l_N
                    break
            # The liquid component flowrates Solvent (water), FORM and MeOH
            return self.outlets

if __name__ == "__main__":
    print(__doc__)
    dict = {'mu1': np.array([0, 1, 2, 3, 4, 5, 6 ,7 ,8 , 9, 10, 11, 12, 13]), 'mu2': np.array([0, 1, 2, 3, 4, 5, 6 ,7 ,8 , 9, 10, 11, 12, 13]), 'mu3': np.array([0, 1, 2, 3, 4, 5, 6 ,7 ,8 , 9, 10, 11, 12, 13]), 'mu4': np.array([0, 1, 2, 3, 4, 5, 6 ,7 ,8 , 9, 10, 11, 12, 13])}
    absorber1 = Absorber(300, 2, 0.97, 5, 3, dict['mu1'], dict['mu2'], dict['mu3'], dict['mu4'])
    recov = {'LK':0.99, 'HK': 0.001}
    distillationColumn = DistillationColumn(300, 2, recov, dict['mu1'], dict['mu2'], dict['mu3'])

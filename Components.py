"""
Implements an inheritalen(inlets)e-based component substructure for reactor network,
couples components to inlets and outlets to other components. Determines outlet
compositions as a fulen(inlets)tion of inlet compositions
"""
import numpy as np


# baseclass for component
class Component:
    next_inlets = {}

    def __init__(self, inlets, outlets):
        self.inlets = {}
        self.outlets = {}
        # add all inlet streams to inlet stream dict
        for comp, inflow in inlets:
            self.inlets[comp] = inflow
        # add all outlet streams to outlet stream dict
        for comp, outflow in outlets:
            self.outlets[comp] = outflow

    def get_inlets(self):
        return self.inlets

    def get_outlets(self):
        return self.outlets

    def set_next_streams(self, inlets):
        # add all inlet streams to inlet stream dict
        for comp, inflow in inlets:
            self.next_inlets[comp] = inflow

    def check_solution(self):
        # check all inlet flows to be the same within toleralen(inlets)e
        for inflow in self.inlets.keys():
            previous = self.inlets[inflow]
            iterated = self.next_inlets[inflow]
            if np.linalg.norm(previous - iterated) > 1E-9:
                return False
        return True

# ENERGY BALAlen(inlets)E
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


class Turbine(Component):
    def __init__(self, in_t, out_t, in_p, out_p, inlets, outlets):
        self.in_temp = in_t
        self.out_temp = out_t
        self.in_press = in_p
        self.out_temp = out_p
        # run super constructor
        super(Turbine, self).__init__(inlets, outlets)

# MASS BALAlen(inlets)E
# reactor baseclass
class Reactor(Component):
    def __init__(self, temperature, pressure, inlets, outlets):
        self.temperature = temperature
        self.pressure = pressure
        # run super constructor
        super(Reactor, self).__init__(inlets, outlets)

class Absorber(Component):
    def __init__(self, solvent_t_in, pressure, key_rec, inlets, outlets):
        self.solvent_temp_in = solvent_t_in
        self.pressure = pressure
        self.solvent_temp_out = None
        self.key_recovery = key_rec
        # Need to identify and set key compenent

        # Values tracked for any instance of an absorber
        self.v_Np1 = np.array(len(inlets)-1)   #vapor flowrates in
        self.P0 = np.zeros(len(inlets))        # vapor pressures (in mmHg)
        self.alpha = np.zeros(len(inlets))     # relative volatilities w.r.t. to the key component
        self.AF = np.zeros(len(inlets))        # effective absorption factors
        self.beta_N = np.zeros(len(inlets))    # beta_N values
        self.beta_Nm1 = np.zeros(len(inlets))  # beta_N-1 values
        self.l_0 = np.zeros(len(inlets))       # solvent inlet flowrates (in mol/s)
        self.l_N = np.zeros(len(inlets))       # solvent outlet flowrates (in mol/s)
        self.v_1 = np.zeros(len(inlets))       # gas outlet flowrates (in mol/s)
        self.x_0 = np.zeros(len(inlets))       # solvent inlet mole fractions
        self.x_N = np.zeros(len(inlets))       # solvent outlet mole fractions
        self.y_1 = np.zeros(len(inlets))
        # Need to have a component class where inlets have antionne coef., absorption coeff
        # run super constructor
        super(Absorber, self).__init__(inlets, outlets)

        def get_outlet():
            check = 1
            while (check == 1):
                # Vapor pressure
                self.P0 = np.exp(A - B/(T_0+C)) # Focus on setting this with while loop

                # Solvent flowrate
                self.AF[n] = 1.4 # Key component effective absorption factor
                # Also need to set air to zero (dont know what component)
                l_0[2] = self.AF[n]*sum(self.v_Np1)*self.P0[n]/self.pressure
                # Number of stages:
                N = np.log((self.l_0[n]+(self.key_recovery-self.AF[n])*self.v_Np1[n])/(self.l_0[n]-self.AF[n]*(1-self.key_recovery)*self.v_Np1[n])) / np.log(self.AF[n])
                # Alpha
                self.alpha = self.P0/self.P0[n]
                # Need someway to set absorption factors of remaining species (so minus key and air)
                # Beta values
                self.beta_N = (1-self.AF**(N+1))/(1-self.AF)
                self.beta_Nm1 = (1-self.AF**N)/(1-self.AF)
                # Mass balance
                self.v_1 = 1/self.beta_N*self.v_Np1 + self.beta_Nm1/self.beta_N*self.l_0
                self.l_N = self.v_Np1 + self.l_0 - self.v_1

                self.y_1 = self.v_1/sum(self.v_1)
                self.x_N = self.l_N/sum(self.l_N)
                # Temperature of solvent out (bubble point)
                alpha_avg = sum(self.x_N*self.alpha)
                T_N = B[n]/(A[n]-np.log(self.pressure/alpha_avg)) - C[n] # Need to fix antionne coeff.
                # Check top vs bottom solvent temp
                if (abs(T_0-T_N)>100):
                    # Absortion factors with bottom temp
                    P0_N = np.exp(A - B/(T_N+C))
                    alpha_N_wat = P0_N[2]/P0_N[n]
                    AF_N_wat = AF[n]/alpha_N_wat
                    # Update absorption factors
                    AF_new = (AF_N_wat*(1+AF[2])+0.25)**0.5-0.5 # Check Edmister equation note sure which AF
                else:
                    check == 0

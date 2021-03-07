"""
This file holds the implementation of reactor conversion files for the
formaldehyde reactor, OME reactor, and methanol reactor for the equilibration
function to determine steady state solution of reactor network.
"""
import numpy as np
from scipy.optimize import newton, fsolve


def MethanolReactor(inlets, temperature, pressure):
    R = 8.314 # Check units
    # Kps Van't Hoff EQ
    Kp1 = -1/R*((get_GRxn(temperature, 'R1') - get_HRxn(temperature, 'R1'))/298 + get_HRxn(temperature, 'R1')/temperature)
    Kp2 = -1/R*((get_GRxn(temperature, 'R2') - get_HRxn(temperature, 'R2'))/298 + get_HRxn(temperature, 'R2')/temperature)
    # Total moles in
    n_total = np.sum(list(inlets.values())[0])
    # Solve for extent of reaction
    extent = newton(FormRxnFunc, [0.5, 0.5], (n_total, list(inlets.values())[0], pressure, Kp1, Kp2))

def MeOHRxnFunc(extent, n_total, inlets, pressure, Kp1, Kp2):
    return [((inlets[4]+extent[0])*(inlets[3]+extent[0]+extent[1])*n_total**2)/((inlets[1]-extent[0]-extent[1])*(inlets[0]-3*extent[0]-extent[1])**3)*(pressure/760)**(-2)-Kp1, \
    ((inlets[2]+extent[1])*(inlets[3]+extent[0]+extent[1]))/((inlets[1]-extent[0]-extent[1])*(inlets[0]-3*extent[0]-extent[1]))-Kp2]


def FormaldehydeReactor(inlets, temperature, pressure):
    # Kp with memo2 equation
    Kp = np.power(10, ((-1.4722E-8)*temperature**3 + (5.2525E-5)*temperature**2 + (-6.889E-2)*temperature + 43))
    # Total moles in
    n_total = np.sum(list(inlets.values())[0])
    # Solve for extent of reaction
    extent = newton(FormRxnFunc, 0.5, (n_total, list(inlets.values())[0], pressure, Kp))
    # Caclulate outlet flow rates of reacting species
    list(inlets.values())[0][3] += extent # H2O
    list(inlets.values())[0][4] -= extent # MeOH
    list(inlets.values())[0][7] -= 0.5*extent # O2
    list(inlets.values())[0][5] += extent # FA
    return inlets

# Helper function to specify nonlinear equation to find extent of reaction
def FormRxnFunc(extent, n_total, inlets, pressure, Kp):
    Kp_guess = (extent*(inlets[3]-extent)*n_total**(-0.5))/((inlets[4]-extent)*(inlets[7]-0.5*extent))**(0.5)*(pressure/760)**(0.5)
    return Kp - Kp_guess


def OMEReactor(inlets, temperature, pressure):
    pass

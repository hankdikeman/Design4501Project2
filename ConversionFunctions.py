"""
This file holds the implementation of reactor conversion files for the
formaldehyde reactor, OME reactor, and methanol reactor for the equilibration
function to determine steady state solution of reactor network.
"""
import numpy as np
from scipy.optimize import newton, fsolve
from thermoutils import get_psat, get_tsat, get_HRxn, get_GRxn

# Methanol reactor model
def MethanolReactor(inlets, temperature, pressure):
    R = 8.314 # Check units
    # Kps Van't Hoff EQ
    Kp1 = -1/R*((get_GRxn(temperature, 'R1') - get_HRxn(temperature, 'R1'))/298 + get_HRxn(temperature, 'R1')/temperature)
    Kp2 = -1/R*((get_GRxn(temperature, 'R2') - get_HRxn(temperature, 'R2'))/298 + get_HRxn(temperature, 'R2')/temperature)
    # Total moles in
    n_total = np.sum(list(inlets.values())[0])
    # Solve for extent of reaction
    extent = fsolve(FormRxnFunc, [0.5, 0.5], (n_total, list(inlets.values())[0], pressure, Kp1, Kp2))
    # Calculate outlet flow rates of reacting species
    list(inlets.values())[0][0] -= 3*extent[0] + extent[1]  # H2
    list(inlets.values())[0][1] -= extent[0] + extent[1]    # CO2
    list(inlets.values())[0][2] += extent[1]                # CO
    list(inlets.values())[0][3] += extent[0] + extent[1]    # H2O
    list(inlets.values())[0][4] -= extent[0]                # MeOH
    return inlets
# Helper function to specify nonlinear equations from EQ constants (Use flsove)
def MeOHRxnFunc(extent, n_total, inlets, pressure, Kp1, Kp2):
    return 
    [
    ((inlets[4]+extent[0])*(inlets[3]+extent[0]+extent[1])*n_total**2)/((inlets[1]-extent[0]-extent[1])*(inlets[0]-3*extent[0]-extent[1])**3)*(pressure/760)**(-2)-Kp1, 
    ((inlets[2]+extent[1])*(inlets[3]+extent[0]+extent[1]))/((inlets[1]-extent[0]-extent[1])*(inlets[0]-3*extent[0]-extent[1]))-Kp2                                    
    ]

# Formaldehyde reactor model
def FormaldehydeReactor(inlets, temperature, pressure):
    # Kp with memo2 equation
    Kp = np.power(10, ((-1.4722E-8)*temperature**3 + (5.2525E-5)*temperature**2 + (-6.889E-2)*temperature + 43))
    # Total moles in
    n_total = np.sum(list(inlets.values())[0])
    # Solve for extent of reaction
    extent = newton(FormRxnFunc, 0.5, (n_total, list(inlets.values())[0], pressure, Kp))
    # Calculate outlet flow rates of reacting species
    list(inlets.values())[0][3] += extent       # H2O
    list(inlets.values())[0][4] -= extent       # MeOH
    list(inlets.values())[0][5] += extent       # FA
    list(inlets.values())[0][7] -= 0.5*extent   # O2
    return inlets
# Helper function to specify nonlinear equation from EQ constants (Use newton)
def FormRxnFunc(extent, n_total, inlets, pressure, Kp):
    Kp_guess = (extent*(inlets[3]-extent)*n_total**(-0.5))/((inlets[4]-extent)*(inlets[7]-0.5*extent))**(0.5)*(pressure/760)**(0.5)
    return Kp - Kp_guess

# OME reactor model
def OMEReactor(inlets, temperature, pressure):
    # Constants for equilibrium constant correlations (Rxns 1-7)
    rxnConstA = np.array([-1.9020, 0.8147, -2.454, -2.454, -2.454, -2.454, -2.454])
    rxnConstB = np.array([3512, 240.25, 3029.6, 3029.6, 3029.6, 3029.6, 3029.6])
    # Initial EQ constants
    K = rxnConstA + rxnConstB/temperature
    # Adjust K for combination of rxns 1 + 2
    Ka = K[0]*K[1]
    Keq = K[2,8] # Slicing Ks 3-7 out of initial K array
    np.insert(Keq, 0, Ka, axis=0) # Double check correct syntax
    # Total moles in
    n_total = np.sum(list(inlets.values())[0])
    # Solve for extent of reaction
    extent = fsolve(OMERxnFunc, [0.5, 0.5, 0.5, 0.5, 0.5], (n_total, list(inlets.values())[0], Keq))
    # Calculate outlet flow rates of reacting species
    list(inlets.values())[0][3] += extent[0]                # H2O
    list(inlets.values())[0][4] -= 2*extent[0]              # MeOH
    list(inlets.values())[0][5] -= extent[0]                # form
    list(inlets.values())[0][8] += extent[0] - extent[1]    # OME1, etc
    list(inlets.values())[0][9] += extent[1] - extent[2]
    list(inlets.values())[0][10] += extent[2] - extent[3]
    list(inlets.values())[0][11] += extent[3] - extent[4]
    list(inlets.values())[0][12] += extent[4] - extent[5]
    list(inlets.values())[0][13] += extent[6]
    return inlets
# Helper function to specify nonlinear equations from EQ constants (Use flsove)
def OMERxnFunc(extent, n_total, inlets, Keq):
    # Nomenclature: extents[0:5]~[A:F], Keq[0:5]~[A:F]
    return [ 
    ((inlets[8]+extent[0]-extent[1])*(inlets[3]+extent[0])*n_total)/((inlets[5]-extent[0]-extent[1]-extent[2]-extent[3]-extent[4])*(inlets[4]-2*extent[0])**2)-Keq[0], 
    ((inlets[9]+extent[1]-extent[2])*n_total)/((inlets[8]+extent[0]-extent[1])*(inlets[5]-extent[0]-extent[1]-extent[2]-extent[3]-extent[4]))-Keq[1], 
    ((inlets[10]+extent[2]-extent[3])*n_total)/((inlets[9]+extent[1]-extent[2])*(inlets[5]-extent[0]-extent[1]-extent[2]-extent[3]-extent[4]))-Keq[2], 
    ((inlets[11]+extent[3]-extent[4])*n_total)/((inlets[10]+extent[2]-extent[3])*(inlets[5]-extent[0]-extent[1]-extent[2]-extent[3]-extent[4]))-Keq[3], 
    ((inlets[12]+extent[4]-extent[3])*n_total)/((inlets[11]+extent[3]-extent[4])*(inlets[5]-extent[0]-extent[1]-extent[2]-extent[3]-extent[4]))-Keq[4], 
    ((inlets[13]+extent[5])*n_total)/((inlets[12]+extent[4]-extent[5])*(inlets[5]-extent[0]-extent[1]-extent[2]-extent[3]-extent[4]))-Keq[5] 
    ]

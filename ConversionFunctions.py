"""
This file holds the implementation of reactor conversion files for the
formaldehyde reactor, OME reactor, and methanol reactor for the equilibration
function to determine steady state solution of reactor network.
"""
import numpy as np
from scipy.optimize import newton, fsolve
import math

# Methanol reactor model
def MethanolReactor(inlets, temperature, pressure):
    # Make outlet array
    new_outlets = list(inlets.values())[0]
    # Rxn 1
    react_ind1 = [0, 1]
    react_coeff1 = [3, 1]
    prod_ind1 = [3, 4]
    prod_coeff1 = [1, 1]
    Grxn1 = 0
    # Rxn 2
    react_ind2 = [0, 1]
    react_coeff2 = [1, 1]
    prod_ind2 = [2, 3]
    prod_coeff2 = [1, 1]
    Grxn2 = 0
    R = 8.314 # Check units
    # Kps Van't Hoff EQ
    Kp1 = math.exp(-1/R*((Grxn1 - get_HRxn(react_ind1, react_coeff1, prod_ind1, prod_coeff1))/298 + get_HRxn(react_ind1, react_coeff1, prod_ind1, prod_coeff1)/temperature))
    Kp2 = math.exp(-1/R*((Grxn2 - get_HRxn(react_ind2, react_coeff2, prod_ind2, prod_coeff2))/298 + get_HRxn(react_ind2, react_coeff2, prod_ind2, prod_coeff2)/temperature))
    # Total moles in
    n_total = np.sum(list(inlets.values())[0])
    # Solve for extent of reaction
    extent = fsolve(MeOHRxnFunc, [0.5, 0.5], (n_total, list(inlets.values())[0], pressure, Kp1, Kp2))
    # Calculate outlet flow rates of reacting species
    new_outlets[0] -= 3*extent[0] + extent[1]  # H2
    new_outlets[1] -= extent[0] + extent[1]    # CO2
    new_outlets[2] += extent[1]                # CO
    new_outlets[3] += extent[0] + extent[1]    # H2O
    new_outlets[4] -= extent[0]                # MeOH
    return new_outlets
# Helper function to specify nonlinear equations from EQ constants (Use flsove)
def MeOHRxnFunc(extent, n_total, inlets, pressure, Kp1, Kp2):
    return
    [
    ((inlets[4]+extent[0])*(inlets[3]+extent[0]+extent[1])*n_total**2)/((inlets[1]-extent[0]-extent[1])*(inlets[0]-3*extent[0]-extent[1])**3)*(pressure/760)**(-2)-Kp1,
    ((inlets[2]+extent[1])*(inlets[3]+extent[0]+extent[1]))/((inlets[1]-extent[0]-extent[1])*(inlets[0]-3*extent[0]-extent[1]))-Kp2
    ]

# Formaldehyde reactor model
def FormaldehydeReactor(inlets, temperature, pressure):
    # Make outlet array
    new_outlets = list(inlets.values())[0]
    # Kp with memo2 equation
    Kp = np.power(10, ((-1.4722E-8)*temperature**3 + (5.2525E-5)*temperature**2 + (-6.889E-2)*temperature + 43))
    # Total moles in
    n_total = np.sum(list(inlets.values())[0])
    # Solve for extent of reaction
    extent = newton(FormRxnFunc, 0.5, (n_total, list(inlets.values())[0], pressure, Kp))
    # Calculate outlet flow rates of reacting species
    new_outlets[3] += extent       # H2O
    new_outlets[4] -= extent       # MeOH
    new_outlets[5] += extent       # FA
    new_outlets[7] -= 0.5*extent   # O2
    return new_outlets
# Helper function to specify nonlinear equation from EQ constants (Use newton)
def FormRxnFunc(extent, n_total, inlets, pressure, Kp):
    Kp_guess = (extent*(inlets[3]-extent)*n_total**(-0.5))/((inlets[4]-extent)*(inlets[7]-0.5*extent))**(0.5)*(pressure/760)**(0.5)
    return Kp - Kp_guess

# OME reactor model
def OMEReactor(inlets, temperature, pressure):
    # Make outlet array
    new_outlets = list(inlets.values())[0]
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
    new_outlets[3] += extent[0]                # H2O
    new_outlets[4] -= 2*extent[0]              # MeOH
    new_outlets[5] -= extent[0]                # form
    new_outlets[8] += extent[0] - extent[1]    # OME1, etc
    new_outlets[9] += extent[1] - extent[2]
    new_outlets[10] += extent[2] - extent[3]
    new_outlets[11] += extent[3] - extent[4]
    new_outlets[12] += extent[4] - extent[5]
    new_outlets[13] += extent[6]
    return new_outlets
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

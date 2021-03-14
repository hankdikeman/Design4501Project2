"""
This file holds the implementation of reactor conversion files for the
formaldehyde reactor, OME reactor, and methanol reactor for the equilibration
function to determine steady state solution of reactor network.
"""
import numpy as np
from scipy.optimize import newton, minimize, Bounds, fsolve
import math
from thermoutils import get_HRxn

# Methanol reactor model
def MethanolReactor(inlets, temperature, pressure):
    # Make outlet array
    new_outlets = np.copy(inlets)
    # Rxn 1
    react_ind1 = [0, 1]
    react_coeff1 = [3, 1]
    prod_ind1 = [3, 4]
    prod_coeff1 = [1, 1]
    Grxn1 = 17.3 # kJ/mol
    # Rxn 2
    react_ind2 = [0, 1]
    react_coeff2 = [1, 1]
    prod_ind2 = [2, 3]
    prod_coeff2 = [1, 1]
    Grxn2 = 28.6 # kJ/mol

    R = 8.314/1000 # kJ/(K*mol)
    # Kps Van't Hoff EQ
    # K1 = math.exp(-1/R*(((Grxn1 - get_HRxn(react_ind1, react_coeff1, prod_ind1, prod_coeff1))/298 + get_HRxn(react_ind1, react_coeff1, prod_ind1, prod_coeff1)/temperature)))
    # K2 = math.exp(-1/R*(((Grxn2 - get_HRxn(react_ind2, react_coeff2, prod_ind2, prod_coeff2))/298 + get_HRxn(react_ind2, react_coeff2, prod_ind2, prod_coeff2)/temperature)))
    K1 = math.exp(-Grxn1/(R*temperature))
    K2 = math.exp(-Grxn2/(R*temperature))
    # Total moles in
    n_total = np.sum(new_outlets)
    # Solve for extent of reaction
    opt_result = minimize(MeOHRxnFunc, (0.5, 0.5), (n_total, new_outlets, pressure, K1, K2))
    # if not opt_result['success']:
    #    raise ValueError('MeOH minimization function did not converge')
    extent = opt_result['x']
    # Calculate outlet flow rates of reacting species
    new_outlets[0] -= (3*extent[0] + extent[1])  # H2
    new_outlets[1] -= (extent[0] + extent[1])    # CO2
    new_outlets[2] += extent[1]                # CO
    new_outlets[3] += (extent[0] + extent[1])    # H2O
    new_outlets[4] += extent[0]                # MeOH
    return new_outlets
# Helper function to specicfy nonlinear equations from EQ constants (Use flsove)
def MeOHRxnFunc(extent_tup, n_total, inlets, pressure, K1, K2):
    extent = [extent_tup[0], extent_tup[1]]
    extent_calc = np.array([
    ((inlets[4]+extent[0])*(inlets[3]+extent[0]+extent[1])*n_total**2)/((inlets[1]-extent[0]-extent[1])*(inlets[0]-3*extent[0]-extent[1])**3)*(pressure/760)**(-2)-K1,
    ((inlets[2]+extent[1])*(inlets[3]+extent[0]+extent[1]))/((inlets[1]-extent[0]-extent[1])*(inlets[0]-3*extent[0]-extent[1]))-K2
    ])
    return np.sum(np.power(extent_calc,2))

# Formaldehyde reactor model
def FormaldehydeReactor(inlets, temperature, pressure):
    # Make outlet array
    new_outlets = np.copy(inlets)
    # Kp with memo2 equation
    # Kp = np.power(10, ((-1.4722E-8)*temperature**3 + (5.2525E-5)*temperature**2 + (-6.889E-2)*temperature + 43))
    # print("K val: ", Kp)
    # # Total moles in
    # n_total = np.sum(new_outlets)
    # Solve for extent of reaction
    # opt_result = minimize(FormRxnFunc, 5, (n_total, new_outlets, pressure, Kp))
    # extent = opt_result['x']
    # print(opt_result)
    extent = new_outlets[4]
    # Calculate outlet flow rates of reacting species
    new_outlets[3] += extent       # H2O
    new_outlets[4] -= extent       # MeOH
    new_outlets[5] += extent       # FA
    new_outlets[7] -= 0.5*extent   # O2
    return new_outlets
# Helper function to specify nonlinear equation from EQ constants (Use newton)
# def FormRxnFunc(extent, n_total, inlets, pressure, Kp):
#     return np.power((extent*(inlets[3]-extent)*n_total**(-0.5))/((inlets[4]-extent)*(inlets[7]-0.5*extent))**(0.5)*(pressure/760)**(0.5)-Kp, 2)

# OME reactor model
def OMEReactor(inlets, temperature, pressure):
    # Make outlet array
    new_outlets = np.copy(inlets)
    # Constants for equilibrium constant correlations (Rxns 1-7)
    rxnConstA = np.array([-1.9020, 0.8147, -2.454, -2.454, -2.454, -2.454, -2.454])
    rxnConstB = np.array([3512, 240.25, 3029.6, 3029.6, 3029.6, 3029.6, 3029.6])
    # Initial EQ constants
    K = np.exp(rxnConstA + rxnConstB/temperature)
    # Adjust K for combination of rxns 1 + 2
    Ka = K[0]*K[1]
    Keq = K[2:7] # Slicing Ks 3-7 out of initial K array
    Keq = np.insert(Keq, 0, Ka, axis=0) # Double check correct syntax
    # Total moles in
    n_total = np.sum(new_outlets)
    print('\n\n\n\nnew outlets:',new_outlets)
    print('Keq',Keq)
    # Solve for extent of reaction
    opt_result = minimize(OMERxnFunc, np.array([2.9,2.3,1.6,1.3,1.3,1.2]), (n_total, new_outlets, Keq), method='SLSQP', bounds=Bounds([0,0,0,0,0,0],[n_total,n_total,n_total,n_total,n_total,n_total]))
    # if not opt_result['success']:
    #     raise ValueError('OME minimization function did not converge')
    extent = opt_result['x']
    print('DID OME CONVERGE??')
    print('Should be zero')
    print(OMERxnFunc(extent, n_total, new_outlets, Keq))
    # Calculate outlet flow rates of reacting species
    new_outlets[3] += extent[0]                # H2O
    new_outlets[4] -= 2*extent[0]              # MeOH
    new_outlets[5] -= extent[0]                # form
    new_outlets[8] += extent[0] - extent[1]    # OME1, etc
    new_outlets[9] += extent[1] - extent[2]
    new_outlets[10] += extent[2] - extent[3]
    new_outlets[11] += extent[3] - extent[4]
    new_outlets[12] += extent[4] - extent[5]
    new_outlets[13] += extent[5]
    return new_outlets
# Helper function to specify nonlinear equations from EQ constants (Use flsove)
def OMERxnFunc(extent_tup, n_total, inlets, Keq):
    # Nomenclature: extents[0:5]~[A:F], Keq[0:5]~[A:F]
    extent = [extent_tup[0], extent_tup[1], extent_tup[2], extent_tup[3], extent_tup[4], extent_tup[5]]
    extent_calc = np.array([
    ((inlets[8]+extent[0]-extent[1]) * (inlets[3]+extent[0]) * n_total) / ((inlets[5]-extent[0]-extent[1]-extent[2]-extent[3]-extent[4]-extent[5]) * (inlets[4]-2*extent[0])**2)-Keq[0],
    ((inlets[9]+extent[1]-extent[2]) * n_total) / ((inlets[8]+extent[0]-extent[1]) * (inlets[5]-extent[0]-extent[1]-extent[2]-extent[3]-extent[4]-extent[5]))-Keq[1],
    ((inlets[10]+extent[2]-extent[3]) * n_total) / ((inlets[9]+extent[1]-extent[2]) * (inlets[5]-extent[0]-extent[1]-extent[2]-extent[3]-extent[4]-extent[5]))-Keq[2],
    ((inlets[11]+extent[3]-extent[4]) * n_total) / ((inlets[10]+extent[2]-extent[3]) * (inlets[5]-extent[0]-extent[1]-extent[2]-extent[3]-extent[4]-extent[5]))-Keq[3],
    ((inlets[12]+extent[4]-extent[5]) * n_total) / ((inlets[11]+extent[3]-extent[4]) * (inlets[5]-extent[0]-extent[1]-extent[2]-extent[3]-extent[4]-extent[5]))-Keq[4],
    ((inlets[13]+extent[5]) * n_total) / ((inlets[12]+extent[4]-extent[5]) * (inlets[5]-extent[0]-extent[1]-extent[2]-extent[3]-extent[4]-extent[5]))-Keq[5]
    ])
    return np.sum(np.power(extent_calc,2))

# # def OMERxnDer(extent, n_total, inlets, Keq):
#     # calculate deriv of extent 1
#     x = extent[0]
#     a = inlets[8] - extent[1]
#     b = inlets[3]
#     c = inlets[5]  - np.sum(extent(
#     return Jacobian(lambda x: OMERxnFunc(x, n_total, inlets, Keq))(x).ravel()
# 
# def OMERxnHess(x, n_total, inlets, Keq):
#     return Hessian(lambda x: OMERxnFunc(x, n_total, inlets, Keq))(x)
# 
# 
# 
# labels = ['H2-0', 'CO2-1', 'CO-2', 'H2O-3', 'MEOH-4', 'FA-5', 'N2-6',
#           'O2-7', 'OME1-8', 'OME2-9', 'OME3-10', 'OME4-11', 'OME5-12', 'OME6-13']

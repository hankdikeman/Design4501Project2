"""
This file holds the implementation of reactor conversion files for the
formaldehyde reactor, OME reactor, and methanol reactor for the equilibration
function to determine steady state solution of reactor network.
"""
import numpy as np
from numpy.random import rand
from scipy.optimize import newton, minimize, root, Bounds, fsolve, brute, differential_evolution, basinhopping
import math
from thermoutils import get_HRxn
import sympy

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
    K1 = math.exp(-Grxn1/(R*temperature))
    K2 = math.exp(-Grxn2/(R*temperature))
    # Total moles in
    n_total = np.sum(new_outlets)
    # Solve for extent of reaction
    opt_result = minimize(MeOHRxnFunc, (6.47*0.5, 6.47*0.5), (n_total, new_outlets, pressure, K1, K2))
    extent = opt_result['x']
    # Calculate outlet flow rates of reacting species
    new_outlets[0] -= (3*extent[0] + extent[1])  # H2
    new_outlets[1] -= (extent[0] + extent[1])    # CO2
    new_outlets[2] += extent[1]                # CO
    new_outlets[3] += (extent[0] + extent[1])    # H2O
    new_outlets[4] += extent[0]                # MeOH
    return new_outlets
# Helper function to specicfy nonlinear equations from EQ constants (Use flsove)
def MeOHRxnFunc(extent_tup, n_tot, inlets, pressure, K1, K2):
    extent = [extent_tup[0], extent_tup[1]]
    # adjust n_total
    n_total = n_tot - np.sum(extent)
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
    # # Total moles in
    # n_total = np.sum(new_outlets)
    # Solve for extent of reaction
    # opt_result = minimize(FormRxnFunc, 5, (n_total, new_outlets, pressure, Kp))
    # extent = opt_result['x']
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
    global OME_inflows
    OME_inflows = inlets
    rxnConstA = np.array([-1.9020, 0.8147, -2.454, -2.454, -2.454, -2.454, -2.454])
    rxnConstB = np.array([3512, 240.25, 3029.6, 3029.6, 3029.6, 3029.6, 3029.6])
    # Initial EQ constants
    K = np.exp(rxnConstA + rxnConstB/temperature)
    # Adjust K for combination of rxns 1 + 2
    Ka = K[0]*K[1]
    Keq = K[2:7] # Slicing Ks 3-7 out of initial K array
    Keq = np.append(Ka, Keq) # Double check correct syntax
    # Total moles in
    n_total = np.sum(new_outlets)
    global OME_Keq
    global OME_lastextent
    global OME_ntotal
    # set initial guess to validated IC if it is not already set
    if not 'OME_lastextent' in globals():
        OME_lastextent = 6.47*np.array([5.51479, 4.17782, 3.35746, 2.12493, 1.19877, 0.5085])
    OME_Keq = Keq
    OME_ntotal = n_total
    # Solve for extent of reaction with basinhopping algorithm and BFGS
    # minimizer_kwargs = {"method": "BFGS",'args': (n_total, new_outlets, Keq)}
    opt_result = root(OMERxnFunc, OME_lastextent, (n_total, new_outlets, Keq), method='lm', jac=OMEJacobian) # , bounds=((0,n_total),(0,n_total),(0,n_total),(0,n_total),(0,n_total),(0,n_total)))
    # basinhopping(OMERxnFunc, np.array([20,20,20,20,20,20], dtype=np.float64), niter=100,T=2,minimizer_kwargs=minimizer_kwargs, accept_test=OMEAccept)
    extent = opt_result['x']
    # write new extent to global IC variable
    OME_lastextent = extent
    # Calculate outlet flow rates of reacting species
    new_outlets[3] += extent[0]                # H2O
    new_outlets[4] -= 2*extent[0]              # MeOH
    new_outlets[5] -= np.sum(extent)               # form
    new_outlets[8] += extent[0] - extent[1]    # OME1, etc
    new_outlets[9] += extent[1] - extent[2]
    new_outlets[10] += extent[2] - extent[3]
    new_outlets[11] += extent[3] - extent[4]
    new_outlets[12] += extent[4] - extent[5]
    new_outlets[13] += extent[5]
    if np.any(new_outlets < 0):
        raise ValueError('the OME equilibrium converged to an invalid solution!!')
    return new_outlets

# Helper function to specify nonlinear equations from EQ constants (Use flsove)
def OMERxnFunc(extent_tup, n_tot, inlets, Keq):
    # Nomenclature: extents[0:5]~[A:F], Keq[0:5]~[A:F]
    extent = [extent_tup[0], extent_tup[1], extent_tup[2], extent_tup[3], extent_tup[4], extent_tup[5]]
    # adjust n_total
    n_total = n_tot - np.sum(extent)
    FA = inlets[5]-extent[0]-extent[1]-extent[2]-extent[3]-extent[4]-extent[5]
    MEOH = inlets[4]-2*extent[0]
    H2O = inlets[3]+extent[0]
    OME1 = inlets[8]+extent[0]-extent[1]
    OME2 = inlets[9]+extent[1]-extent[2]
    OME3 = inlets[10]+extent[2]-extent[3]
    OME4 = inlets[11]+extent[3]-extent[4]
    OME5 = inlets[12]+extent[4]-extent[5]
    OME6 = inlets[13]+extent[5]
    # calculate extents of reaction
    extent_calc = np.array([
    (OME1 * H2O * n_total) / (FA * np.power(MEOH,2))-Keq[0],
    (OME2 * n_total) / (OME1 * FA)-Keq[1],
    (OME3 * n_total) / (OME2 * FA)-Keq[2],
    (OME4 * n_total) / (OME3 * FA)-Keq[3],
    (OME5 * n_total) / (OME4 * FA)-Keq[4],
    (OME6 * n_total) / (OME5 * FA)-Keq[5]
    ])
    # return the mean squared error
    return extent_calc

def OMEAccept(**kwargs):
    global OME_inflows
    extent = kwargs['x_new']
    rejected = np.any(
        np.array([
        OME_inflows[3] + extent[0],
        OME_inflows[4] - 2*extent[0],              # MeOH
        OME_inflows[5] - np.sum(extent),               # form
        OME_inflows[8] + extent[0] - extent[1],    # OME1, etc
        OME_inflows[9] + extent[1] - extent[2],
        OME_inflows[10] + extent[2] - extent[3],
        OME_inflows[11] + extent[3] - extent[4],
        OME_inflows[12] + extent[4] - extent[5],
        OME_inflows[13] + extent[5]
        ]) < 0
        )
    return not rejected

def OMEJacobian(extent, n_total, new_outlets, Keq): # Pass in np.array of extents
    x0, x1, x2, x3, x4, x5 = sympy.symbols('x0,x1,x2,x3,x4,x5', real=True)
    # use global parameters from above
    global OME_inflows
    global OME_Keq
    global OME_ntotal
    inlets = OME_inflows
    Keq = OME_Keq
    n_tot = OME_ntotal
    n_total = n_tot - np.sum(extent)

    f1 = ((inlets[8]+x0-x1) * (inlets[3]+x0) * n_total) / ((inlets[5]-x0-x1-x2-x3-x4-x5) * np.power(inlets[4]-2*x0,2))-Keq[0]
    f2 = ((inlets[9]+x1-x2) * n_total) / ((inlets[8]+x0-x1) * (inlets[5]-x0-x1-x2-x3-x4-x5))-Keq[1]
    f3 = ((inlets[10]+x2-x3) * n_total) / ((inlets[9]+x1-x2) * (inlets[5]-x0-x1-x2-x3-x4-x5))-Keq[2]
    f4 = ((inlets[11]+x3-x4) * n_total) / ((inlets[10]+x2-x3) * (inlets[5]-x0-x1-x2-x3-x4-x5))-Keq[3]
    f5 = ((inlets[12]+x4-x5) * n_total) / ((inlets[11]+x3-x4) * (inlets[5]-x0-x1-x2-x3-x4-x5))-Keq[4]
    f6 = ((inlets[13]+x5) * n_total) / ((inlets[12]+x4-x5) * (inlets[5]-x0-x1-x2-x3-x4-x5))-Keq[5]

    F = sympy.Matrix([f1,f2, f3, f4, f5, f6])
    J = F.jacobian([x0,x1,x2,x3,x4,x5]).subs([(x0,extent[0]), (x1,extent[1]), (x2,extent[2]), (x3,extent[3]), (x4,extent[4]), (x5,extent[5])])
    
    return np.array(J).astype(np.float64)


# labels = ['H2-0', 'CO2-1', 'CO-2', 'H2O-3', 'MEOH-4', 'FA-5', 'N2-6',
#           'O2-7', 'OME1-8', 'OME2-9', 'OME3-10', 'OME4-11', 'OME5-12', 'OME6-13']

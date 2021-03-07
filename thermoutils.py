"""
This file contains all the thermodynamics utility files, such as Antoine,
bubble point, enthalpies, etc
"""
import numpy as np

labels = ['H2', 'CO2', 'CO', 'H2O', 'MEOH', 'FA', 'N2',
          'O2', 'OME1', 'OME2', 'OME3', 'OME4', 'OME5', 'OME6']

A = np.array([-100, -100, -100, 18.3036, 18.5875, 16.4775, -100,
     -100, 15.8714, 14.6956, 15.9303, 15.728, 16.2705, 18.511])
B = np.array([-100, -100, -100, 3816.44, 3626.55, 2204.13, -100,
     -100, 2428.14, 2313.66, 3283.16, 3448.25, 4067.99, 5642.97])
C = np.array([-100, -100, -100, -46.13, -34.29, -30.15, -100,
     -100, -52.35, -91.73, -76.75, -96.99, -95.59, -65.12])


# get saturation pressure from temperature
def get_psat(temperature):
    return np.exp(A - B / (temperature + C))


# get saturation temperature from saturation pressure
def get_tsat(psat):
    return B / (A - np.log(psat)) - C

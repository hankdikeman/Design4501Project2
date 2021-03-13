"""
This file contains all the thermodynamics utility files, such as Antoine,
bubble point, enthalpies, etc
"""
import numpy as np

labels = ['H2-0', 'CO2-1', 'CO-2', 'H2O-3', 'MEOH-4', 'FA-5', 'N2-6',
          'O2-7', 'OME1-8', 'OME2-9', 'OME3-10', 'OME4-11', 'OME5-12', 'OME6-13']

# Antoinne Coeff.
A = np.array([13.6333, 22.5898, 14.3686, 18.3036, 18.5875, 16.4775, 14.9342,
              15.4075, 15.8714, 14.6956, 15.9303, 15.728, 16.2705, 18.511])
B = np.array([164.90, 3103.39, 530.22, 3816.44, 3626.55, 2204.13, 588.72,
              734.55, 2428.14, 2313.66, 3283.16, 3448.25, 4067.99, 5642.97])
C = np.array([3.19, -0.16, -13.15, -46.13, -34.29, -30.15, -6.6,
              -6.45, -52.35, -91.73, -76.75, -96.99, -95.59, -65.12])


# Heat of formations (kJ/mol)
Hf = np.array([0, -393.5, -110.5, -241.8, -201.3, -108.6, 0,
              0, -369.447, -522.205, -674.963, -827.721, -980.479, -1285.994])
# Heat vaporization (kJ/mol)
Hvap = np.array([0.44936, 15.326, 6.0, 40.7, 35.21, 23.3, 5.56,
                3.4099, 27.058, 31.685, 36.317, 40.949, 45.580, 50.212, 54.844])
# Specific heat capacity
Cp = np.array([])

# get saturation pressure from temperature
def get_psat(temperature):
    return np.exp(A - B / (temperature + C))


# get saturation temperature from saturation pressure
def get_tsat(psat):
    return B / (A - np.log(psat)) - C

# get heat of reaction (kj/mol)
def get_HRxn(reactant_indexs, reactant_coeff, product_indexs, product_coeff):
    reactant_sum = 0
    product_sum = 0
    # Sum reactant heat of formations with coeff
    for i in range(len(reactant_indexs)):
        reactant_sum += reactant_coeff[i]*Hf[reactant_index[i]]
    # Sum product heat of formations with coeff
    for j in range(len(product_indexs)):
        product_sum += product_coeff[j]*Hf[product_index[j]]
    return (product_sum - reactant_sum)

if __name__ == "__main__":
    rel_vol_h20 = get_psat(298)/get_psat(298)[3]

    for i in range(len(labels)):
        print(labels[i],rel_vol_h20[i])

# Delta G of rxn calculations for methanol reactor (kJ/mol)
#  CO2:-394.4 H2:0 MeOH:-148.5 H2O:-228.6 CO:-137.2
# RXN1: 17.3
# RXN2: 28.6

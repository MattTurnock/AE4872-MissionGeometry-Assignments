#Purpose is to convert from cartesian --> Kepler and vice-versa

import Kep2Cart_utils
import kep_orbit_utils
import numpy as np

def do_tests(state_cart, state_kep, mu):
    SMA, ECC, INC, Omega, omega, theta = state_kep
    state_cart_calculated = Kep2Cart_utils.Kep2Cart(SMA, ECC, INC, Omega, omega, theta, mu)

    rvector =np.vstack(state_cart_cryo[0:3])
    Vvector = np.vstack(state_cart_cryo[3:6])



#Set variables to test against. All numbers in km and rad
# for keplerian state should be [SMA, ECC, INC, Omega, omega, theta]
state_cart_iss = [-2700816.14E-3, -3314092.8E-3, 5266346.42E-3, 5168.606550E-3, -5597.546618E-3, -868.878445E-3]
state_kep_iss = [6787746.891E-3, 0.000731104, np.deg2rad(51.68714486), np.deg2rad(127.5486706), np.deg2rad(74.21987137), np.deg2rad(24.10027677)]

state_cart_cryo = [3126974.99E-3, -6374445.74E-3, 28673.59E-3, -254.91197E-3, -83.30107E-3, 7485.70674E-3]
state_kep_cryo = [7096137.00E-3, 0.0011219, np.deg2rad(92.0316), np.deg2rad(296.1384), np.deg2rad(120.6878), np.deg2rad(239.5991)]

print(np.vstack(state_cart_cryo[3:6]))




# mu = kep_orbit_utils.mu_Earth
# SMA, ECC, INC, Omega, omega, theta = state_kep_iss
# ISS_state = Kep2Cart_utils.Kep2Cart(SMA, ECC, INC, Omega, omega, theta, mu)
#
# print('State vector of ISS in km and km/s \n' + str(ISS_state))
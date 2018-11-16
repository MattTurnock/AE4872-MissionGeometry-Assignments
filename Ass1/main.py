#Purpose is to convert from cartesian --> Kepler and vice-versa

import Kep2Cart_utils
import kep_orbit_utils
import numpy as np

SMA = 6787746.891E-3
ECC = 0.000731104
INC = np.deg2rad(51.68714486)
Omega = np.deg2rad(127.5486706)
omega = np.deg2rad(74.21987137)
theta = np.deg2rad(24.10027677)
mu = kep_orbit_utils.mu_Earth
ISS_state = Kep2Cart_utils.Kep2Cart(Omega, omega, INC, SMA, ECC, theta, mu)

print('State vector of ISS in km and km/s \n' + str(ISS_state))
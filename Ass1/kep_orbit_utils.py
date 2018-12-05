#Defines a selection of relevant orbit-related conversions/ functions and constants
import numpy as np
from Ass1 import misc_utils
import astropy.units as u
cos = np.cos
sin = np.sin
pi = np.pi

mu_Earth = 398600.441E9
mu_Sun = 1.327178E11
mu_Jup = 1.2669E8
mu_Mars = 4.2832E4

#Function to take the Keplerian state and convert the angles to degrees
def Kep_Rad2Deg(state_kep_in):
    state_kep = state_kep_in.copy()
    for i in range(len(state_kep)):
        if i!= 0 and i!= 1:
            state_kep[i] = np.rad2deg(state_kep[i])

    return state_kep

#reverse of Kep_Rad2Deg
def Kep_Deg2Rad(state_kep_in):
    state_kep = state_kep_in.copy()
    for i in range(len(state_kep)):
        if i!= 0 and i!= 1:
            state_kep[i] = np.deg2rad(state_kep[i])

    return state_kep

def get_SMA(orbital_radius, mu, V):
    return 1/(2/orbital_radius - V**2/mu)

def get_ECCvector(Vvector, hvector, rvector, mu):
    r = np.linalg.norm(rvector)
    return np.cross(Vvector, hvector, axis=0)/mu - rvector/r

def get_INC(hvector):
    hz = hvector[2][0]
    hmag = np.linalg.norm(hvector)
    return np.arccos(hz/hmag)

def get_hvector(rvector, Vvector):
    return np.cross(rvector, Vvector, axis=0)

def get_N(hvector):
    unitz = np.vstack((0,0,1))
    return np.cross(unitz, hvector, axis=0)

def get_orbit_radius(SMA, ECC, theta):
    return SMA*(1-ECC**2)/(1 + ECC*np.cos(theta))

def get_orbit_radius_2(p, ECC, theta):
    return p / (1 + ECC * np.cos(theta))

def get_radius_Pe(SMA, ECC):
    return SMA * (1-ECC)

def get_radius_Ap(SMA, ECC):
    return SMA * (1+ECC)

def get_orbital_velocity(mu, orbital_radius, SMA):
    return mu * (2/orbital_radius - 1/SMA)

def get_circular_velocity(mu, orbital_radius):
    return np.sqrt(mu/orbital_radius)

def get_escape_velocity(mu, orbital_radius):
    return np.sqrt(2*mu/orbital_radius)

def get_orbit_period(SMA, mu, astropy_units=None):
    T = 2 * pi * np.sqrt(SMA ** 3 / mu)
    if astropy_units!=None:
        return T.to(u.Unit(astropy_units))
    else:
        return T

def get_H(mu, SMA, ECC):
    return np.sqrt(mu * SMA * (1 - ECC**2))

def E2theta(E, ECC):
    root = np.sqrt((1+ECC)/(1-ECC))
    inside = root * np.tan(E/2)
    theta = 2*np.arctan(inside)
    if theta < 0:
        theta = misc_utils.angle2positive(theta)
    return theta

def M2E(M, ECC, verbose=False):
    precision_in = misc_utils.get_precision(M)
    precision = precision_in + 1
    E_i = M
    E_i_rounded = round(M, precision)
    running = True

    iterations = 0
    while running:
        E_ip1 = E_i + (M - E_i + ECC*np.sin(E_i))/(1-ECC*np.cos(E_i))
        E_ip1_rounded = round(E_ip1, precision)

        iterations = iterations + 1
        if E_i_rounded == E_ip1_rounded:
            running = False
        E_i_rounded = E_ip1_rounded
        E_i = E_ip1

        if iterations == 100:
            print('RUNTIME ERROR IN FUNCTION M2E')
            running=False
        #print(np.rad2deg(np.array([E_i, E_ip1])))

    if E_i < 0:
        E_i = misc_utils.angle2positive(E_i)

    if verbose:
        return E_i, iterations
    else: return E_i

def theta2E(theta, ECC):
    root = np.sqrt((1-ECC)/(1+ECC))
    theta = 2*np.arctan(root*np.tan(theta/2))
    if theta <0:
        theta = misc_utils.angle2positive(theta)
    return theta

def E2M(E, ECC):
    M = E - ECC*np.sin(E)
    if M <0:
        M = misc_utils.angle2positive(M)
    return M







#Defines a few basic keplerian orbit functions and constants
import numpy as np
cos = np.cos
sin = np.sin
pi = np.pi

mu_Earth = 398600.441E9
mu_Sun = 1.327178E11
mu_Jup = 1.2669E8
mu_Mars = 4.2832E4

def Kep_Rad2Deg(state_kep_in):
    state_kep = state_kep_in.copy()
    for i in range(len(state_kep)):
        if i!= 0 and i!= 1:
            state_kep[i] = np.rad2deg(state_kep[i])

    return state_kep

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

def get_orbit_period(SMA, mu):
    return 2*pi*np.sqrt(SMA**3 / mu)

def get_H(mu, SMA, ECC):
    return np.sqrt(mu * SMA * (1 - ECC**2))




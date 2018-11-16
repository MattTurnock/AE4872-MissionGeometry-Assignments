#Defines a few basic keplerian orbit functions and constants
import numpy as np
cos = np.cos
sin = np.sin
pi = np.pi

mu_Earth = 398600.441
mu_Sun = 1.327178E11
mu_Jup = 1.2669E8
mu_Mars = 4.2832E4

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
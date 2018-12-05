from json_to_dict import constants
import numpy as np
from astropy import units as u
from Ass1.kep_orbit_utils import get_orbit_period
from Ass1.misc_utils import get_precision
pi = np.pi
u.sday = u.def_unit('sday', 86164.1004*u.s)
#print((1*u.sday).to(u.s))

#Functions for less accurate measurement (only Omega)
# def get_H0(j_k, mu=constants["muEarth"], DS=constants["DS"], R=constants["RE"]):
#     r0 = mu**(1/3)*((2*pi)/(DS)*j_k)**(-2/3)
#     H0 = r0 - R
#     return H0

# #needs astropy.unit inputs, outputs in rad (but is actually rad/ revolution)
# def get_DL1(TE, T, astropy_units='rad'):
#     DL1 = -2*pi*T/TE
#     DL1 = DL1.decompose()
#     DL1 = DL1 * u.rad
#     return DL1.to(u.Unit(astropy_units))
#
# def get_i_from_a(a, j_k, J2=constants["J2"], R=constants["RE"], mu=constants["muEarth"], TE=u.sday):
#     T = get_orbit_period(a, mu, astropy_units='s')
#     DL1 = get_DL1(TE, T)
#     jk_term = (j_k**-1 * 2*pi)*u.rad
#     DL2_noi = (a**2/(3*pi*J2*R**2))/u.rad
#
#     inner_bracket_p = DL1 + jk_term
#     inner_bracket_n = DL1 - jk_term
#
#     print(DL2_noi*inner_bracket_p)
#
#     i_p = np.arccos(DL2_noi*inner_bracket_p)
#     i_n = np.arccos(DL2_noi*inner_bracket_n)
#
#     return i_p.to(u.deg), i_n
#
# a = 7258.69*u.km
# jk = 14
#
# T = get_orbit_period(a, constants["muEarth"], astropy_units='minute')
# print(T)
# DL1 = get_DL1(u.sday, T, astropy_units='deg')
# print(DL1, '\n')
#
# i_p_i_n = get_i_from_a(a, jk)
# print(i_p_i_n)
###############################################################################################
#Functions for more accurate things
def get_a0(j_k, mu=constants["muEarth"], DS=constants["DS"], return_unit=(u.km)):
    frac = (DS/(2*pi)*j_k**-1)**(2/3)
    a0 = frac * mu**(1/3)
    return a0.to(return_unit)

#Functions for more accurate estimate
def get_Omegadot(a, i, e=0, J2=constants["J2"], mu=constants["muEarth"], R=constants["RE"], DS=constants["DS"], return_unit=(u.deg/u.sday)):
    #this k2 doesn't have the right value because of fancy stuff the book does with units
    k2 = (0.75 * J2 * mu**0.5 * R**2)
    #calulate a working version of Omegadot
    Omegadot_temp = -2* k2 * a**(-7/2) * np.cos(i) * (1 - e**2)**(-2)
    #convert to proper units, and take value
    Omegadot_temp = (Omegadot_temp.to(u.s**-1)).value
    #convert to deg/sday
    Omegadot = (Omegadot_temp * 180/pi * DS.value)*(u.deg/u.sday)
    return Omegadot.to(return_unit)

def get_omegadot(a, i, e=0, J2=constants["J2"], mu=constants["muEarth"], R=constants["RE"], DS=constants["DS"], return_unit=(u.deg/u.sday)):
    #this k2 doesn't have the right value because of fancy stuff the book does with units
    k2 = (0.75 * J2 * mu**0.5 * R**2)

    #calulate a working version of Omegadot
    omegadot_temp = k2 * a**(-7/2) * (5*(np.cos(i))**2 - 1)*(1-e**2)**(-2)
    #convert to proper units, and take value
    omegadot_temp = (omegadot_temp.to(u.s**-1)).value
    #convert to deg/sday
    omegadot = (omegadot_temp * 180/pi * DS.value)*(u.deg/u.sday)
    return omegadot.to(return_unit)

def get_Mdot(a, i, e=0, J2=constants["J2"], mu=constants["muEarth"], R=constants["RE"], DS=constants["DS"], return_unit=(u.deg/u.sday)):
    #this k2 doesn't have the right value because of fancy stuff the book does with units
    k2 = (0.75 * J2 * mu**0.5 * R**2)
    #calulate a working version of Omegadot
    Mdot_temp = k2 * a**(-7/2) * (3*(np.cos(i))**2 - 1)*(1-e**2)**(-3/2)
    #convert to proper units, and take value
    Mdot_temp = (Mdot_temp.to(u.s**-1)).value
    #convert to deg/sday
    Mdot = (Mdot_temp * 180/pi * DS.value)*(u.deg/u.sday)
    return Mdot.to(return_unit)

def get_n(j_k, Omegadot, omegadot=None, Mdot=None, Ldot=constants["Ldot"], return_unit=(u.deg/u.sday), approach='1'):
    if approach == '1':
        n = j_k * (Ldot - Omegadot)
    elif approach == '2':
        n = j_k*(Ldot -Omegadot) - (omegadot + Mdot)
    return n.to(return_unit)

def get_a(n, mu=constants["muEarth"], DS=constants["DS"], return_unit=(u.km)):
    #make sure everything is in units we want, and take value
    mu_val = (mu.to(u.km**3/(u.s**2))).value
    n_val = (n.to(u.deg/u.sday)).value
    DS_val = (DS.to(u.s)).value
    #perform calculation
    bottom = ((n_val * pi) / (DS_val * 180)) ** 2
    a_val = (mu_val/bottom)**(1/3)
    a = a_val*(u.km)
    return a.to(return_unit)

def do_a_iterations(j_k, i, e=0, da_precision=5, return_its=False, approach='1'):
    a0 = get_a0(j_k)
    iterations=0
    running=True
    while running:
        #calc the components
        Omegadot = get_Omegadot(a0,i)
        if approach=='1':
            n = get_n(j_k, Omegadot, approach=approach)
        elif approach=='2':
            omegadot = get_omegadot(a0,i)
            Mdot = get_Mdot(a0,i)
            n = get_n(j_k, Omegadot, omegadot=omegadot, Mdot=Mdot, approach=approach)
        a_new = get_a(n)
        #find difference
        da = a_new - a0

        a0 = a_new
        iterations += 1
        #set convergence condition as the precision level given
        da_rounded = round(da.value, da_precision)
        if da_rounded == 0:
            running = False

    if return_its: return a_new, iterations
    else: return a_new


# j_k = 14
# a=7258.69*1000*u.m
# i=28*u.deg
#
# #print(do_a_iterations(j_k, i, return_its=True, da_precision=1))
# #print(get_H0(j_k))
# Omegadot = get_Omegadot(a, i)
# print(Omegadot)
# omegadot = get_omegadot(a, i)
# #print(omegadot)
# Mdot = get_Mdot(a, i)
# #print(Mdot)
#
# a0 = get_a0(j_k)
# #print(a0)
# n = get_n(j_k, Omegadot, omegadot, Mdot)
# #print(n)
# a = get_a(n)
# #print(a)
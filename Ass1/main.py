#Purpose is to convert from cartesian --> Kepler and vice-versa

import Kep2Cart_utils
import Cart2Kep_utils
import kep_orbit_utils
import numpy as np

def do_conversion(state_cart, state_kep, mu=kep_orbit_utils.mu_Earth):
    #SMA, ECC, INC, Omega, omega, theta = state_kep
    state_cart_calculated = Kep2Cart_utils.Kep2Cart(state_kep)
    state_kep_calculated = Cart2Kep_utils.Cart2Kep(state_cart)
    return state_cart_calculated, state_kep_calculated

#state_precision is the number of decimal places present in the true state
def round_state(state, state_calculated):#, state_precision):
    state_str = str(state).replace('[','').replace(']','').split(',')
    state_precision=[]
    for i in range(len(state)):
        number = state_str[i]
        decimals = number[::-1].find('.')
        state_precision.append(decimals)

    state_calculated_rounded = []
    for i in range(len(state_calculated)):
        rounded = round(state_calculated[i], state_precision[i])
        state_calculated_rounded.append(rounded)
    return state_calculated_rounded


#Set variables to test against. All numbers in km and rad
# for keplerian state should be [SMA, ECC, INC, Omega, omega, theta]

#ISS given state
state_cart_iss = [-2700816.14, -3314092.8, 5266346.42, 5168.606550, -5597.546618, -868.878445]
state_kep_iss_degrees = [6787746.891, 0.000731104, 51.68714486, 127.5486706, 74.21987137, 24.10027677]

#convert to radians
state_kep_iss = kep_orbit_utils.Kep_Deg2Rad(state_kep_iss_degrees)
#perform conversion, rounding and conversion to degrees
state_cart_iss_calculated, state_kep_iss_calculated = do_conversion(state_cart_iss, state_kep_iss)
state_cart_iss_calculated_rounded = round_state(state_cart_iss, state_cart_iss_calculated)
state_kep_iss_calculated_degrees = kep_orbit_utils.Kep_Rad2Deg(state_kep_iss_calculated)
state_kep_iss_calculated_degrees_rounded = round_state(state_kep_iss_degrees, state_kep_iss_calculated_degrees)



state_cart_cryo = [3126974.99, -6374445.74, 28673.59, -254.91197, -83.30107, 7485.70674]
state_kep_cryo_degrees = [7096137.00, 0.0011219, 92.0316, 296.1384, 120.6878, 239.5437]

#convert to radians
state_kep_cryo = kep_orbit_utils.Kep_Deg2Rad(state_kep_cryo_degrees)
#perform conversion, rounding and conversion to degrees
state_cart_cryo_calculated, state_kep_cryo_calculated = do_conversion(state_cart_cryo, state_kep_cryo)
state_cart_cryo_calculated_rounded = round_state(state_cart_cryo, state_cart_cryo_calculated)
state_kep_cryo_calculated_degrees = kep_orbit_utils.Kep_Rad2Deg(state_kep_cryo_calculated)
state_kep_cryo_calculated_degrees_rounded = round_state(state_kep_cryo_degrees, state_kep_cryo_calculated_degrees)



state_cart_cryo_calculated, state_kep_cryo_calculated = do_conversion(state_cart_cryo, state_kep_cryo)


print('\n=========================================ISS_TEST=====================================================\n')
print('ISS true cartesian coordinates:\n%s \n\nISS converted cartesian coordinates:\n%s\n' %(state_cart_iss, state_cart_iss_calculated_rounded))
print('ISS true keplerian coordinates:\n%s \n\nISS converted keplerian coordinaes:\n%s\n' %(state_kep_iss_degrees, state_kep_iss_calculated_degrees_rounded))
print('\n=========================================CRYO_TEST=====================================================\n')
print('Cryo true cartesian coordinates:\n%s \n\nCryo converted cartesian coordinates:\n%s\n' %(state_cart_cryo, state_cart_cryo_calculated_rounded))
print('Cryo true keplerian coordinates:\n%s \n\nCryo converted keplerian coordinaes:\n%s\n' %(state_kep_cryo_degrees, state_kep_cryo_calculated_degrees_rounded))

#Purpose is to convert from cartesian --> Kepler and vice-versa


from Ass1 import kep_orbit_utils
from Ass1 import Kep2Cart_utils
from Ass1 import Cart2Kep_utils
from Ass1 import latex_utils, misc_utils
import json
import os
import numpy as np
AE4878_path = os.path.dirname(os.path.realpath(__file__))
with open(os.path.join(AE4878_path,'constants.json')) as handle:
    course_constants = json.loads(handle.read())

def do_conversion(state_cart, state_kep, mu=kep_orbit_utils.mu_Earth):
    #SMA, ECC, INC, Omega, omega, theta = state_kep
    state_cart_calculated = Kep2Cart_utils.Kep2Cart(state_kep, given_angle='theta')
    state_kep_calculated = Cart2Kep_utils.Cart2Kep(state_cart)
    return state_cart_calculated, state_kep_calculated

#state_precision is the number of decimal places present in the true state
def round_state(state, state_calculated):
    state_str = str(state).replace('[','').replace(']','').split(',')
    state_precision=[]
    for i in range(len(state)):
        number = state_str[i]
        decimals = misc_utils.get_precision(number)
        state_precision.append(decimals)

    state_calculated_rounded = []
    for i in range(len(state_calculated)):
        rounded = round(state_calculated[i], state_precision[i])
        state_calculated_rounded.append(rounded)
    return state_calculated_rounded


#Set variables to test against. All numbers in km and rad
# for keplerian state should be [SMA, ECC, INC, Omega, omega, theta]

#ISS given state
state_cart_iss = [-2700816.14, -3314092.8, 5266346.42, 5168.606550, -5597.546618, -868.878445]          #in m and m/s
state_kep_iss_degrees = [6787746.891, 0.000731104, 51.68714486, 127.5486706, 74.21987137, 24.10027677]  #in m and degrees, given theta

#convert to radians
state_kep_iss = kep_orbit_utils.Kep_Deg2Rad(state_kep_iss_degrees)
#perform conversion, rounding and conversion to degrees
state_cart_iss_calculated, state_kep_iss_calculated = do_conversion(state_cart_iss, state_kep_iss)
state_cart_iss_calculated_rounded = round_state(state_cart_iss, state_cart_iss_calculated)
state_kep_iss_calculated_degrees = kep_orbit_utils.Kep_Rad2Deg(state_kep_iss_calculated)
state_kep_iss_calculated_degrees_rounded = round_state(state_kep_iss_degrees, state_kep_iss_calculated_degrees)

state_cart_cryo = [3126974.99, -6374445.74, 28673.59, -254.91197, -83.30107, 7485.70674]                # in m and m/s
state_kep_cryo_degrees = [7096137.00, 0.0011219, 92.0316, 296.1384, 120.6878, 239.5437]                 # in m and degrees, given theta

#convert to radians
state_kep_cryo = kep_orbit_utils.Kep_Deg2Rad(state_kep_cryo_degrees)
#perform conversion, rounding and conversion to degrees
state_cart_cryo_calculated, state_kep_cryo_calculated = do_conversion(state_cart_cryo, state_kep_cryo)
state_cart_cryo_calculated_rounded = round_state(state_cart_cryo, state_cart_cryo_calculated)
state_kep_cryo_calculated_degrees = kep_orbit_utils.Kep_Rad2Deg(state_kep_cryo_calculated)
state_kep_cryo_calculated_degrees_rounded = round_state(state_kep_cryo_degrees, state_kep_cryo_calculated_degrees)

state_cart_cryo_calculated, state_kep_cryo_calculated = do_conversion(state_cart_cryo, state_kep_cryo)

ISS_theta = state_kep_iss_degrees[-1]
ISS_E = 24.08317766  #in degrees
ISS_M = 24.06608426 #in degrees
ISS_ECC = state_kep_iss_degrees[1]
ISS_theta_calculated = np.rad2deg(kep_orbit_utils.E2theta(np.deg2rad(ISS_E), ISS_ECC)) #in degrees, from E
ISS_theta_calculated_rounded = round(ISS_theta_calculated, misc_utils.get_precision(ISS_theta))
ISS_E_calculated_rads, ISS_iterations = kep_orbit_utils.M2E(np.deg2rad(ISS_M), ISS_ECC, verbose=True)
ISS_E_calculated = np.rad2deg(ISS_E_calculated_rads)
ISS_E_calculated_rounded = round(ISS_E_calculated, misc_utils.get_precision(ISS_E))

cryo_theta = state_kep_cryo_degrees[-1]
cryo_E = 239.5991  #in degrees
cryo_M = 239.6546 #in degrees
cryo_ECC = state_kep_cryo_degrees[1]
cryo_theta_calculated = np.rad2deg(kep_orbit_utils.E2theta(np.deg2rad(cryo_E), cryo_ECC)) #in degrees, from E
cryo_theta_calculated_rounded = round(cryo_theta_calculated, misc_utils.get_precision(cryo_theta))
cryo_E_calculated_rads, cryo_iterations = kep_orbit_utils.M2E(np.deg2rad(cryo_M), cryo_ECC, verbose=True)
cryo_E_calculated = np.rad2deg(cryo_E_calculated_rads)
cryo_E_calculated_rounded = round(cryo_E_calculated, misc_utils.get_precision(cryo_E))


print('\n=========================================theta2E MODULE TEST=====================================================\n')
print('Cryo true value for [theta, E] = [%s, %s]' %(cryo_theta, cryo_E))
print('Cryo calculated value for E (given theta): E=%s\n' %(round(np.rad2deg(kep_orbit_utils.theta2E(np.deg2rad(cryo_theta), cryo_ECC)), misc_utils.get_precision(cryo_E))))

print('\n=========================================E2M MODULE TEST=====================================================\n')
print('Cryo true value for [E,M] = [%s, %s]' %(cryo_E, cryo_M))
print('Cryo calculated value for M (given E): M=%s\n' %(round(np.rad2deg(kep_orbit_utils.E2M(np.deg2rad(cryo_E), cryo_ECC)), misc_utils.get_precision(cryo_M)+1)))


print('\n=========================================E2theta MODULE TEST=====================================================\n')
doISS = False
if doISS:
    print('ISS true value for [theta, E] = [%s, %s]' %(ISS_theta, ISS_E))
    print('ISS calculated value for theta (given E): theta=%s\n' %(ISS_theta_calculated_rounded))
print('Cryo true value for [theta, E] = [%s, %s]' %(cryo_theta, cryo_E))
print('Cryo calculated value for theta (given E): theta=%s\n' %(cryo_theta_calculated_rounded))


print('\n=========================================M2E MODULE TEST=====================================================\n')
if doISS:
    print('ISS true value for [E, M] = [%s, %s]' %(ISS_E, ISS_M))
    print('ISS calculated value for E (given M): E=%s, in %s iterations\n' %(ISS_E_calculated_rounded, ISS_iterations))
print('cryo true value for [E, M] = [%s, %s]' %(cryo_E, cryo_M))
print('cryo calculated value for E (given M): E=%s, in %s iterations\n' %(cryo_E_calculated_rounded, cryo_iterations))



print('\n=========================================ISS TEST FULL=====================================================\n')
print('ISS true cartesian coordinates:\n%s \nISS converted cartesian coordinates:\n%s\n' %(state_cart_iss, state_cart_iss_calculated_rounded))
print('ISS true keplerian coordinates:\n%s \nISS converted keplerian coordinaes:\n%s\n' %(state_kep_iss_degrees, state_kep_iss_calculated_degrees_rounded))
print('\n=========================================CRYO TEST FULL=====================================================\n')
print('Cryo true cartesian coordinates:\n%s \nCryo converted cartesian coordinates:\n%s\n' %(state_cart_cryo, state_cart_cryo_calculated_rounded))
print('Cryo true keplerian coordinates:\n%s \nCryo converted keplerian coordinaes:\n%s\n' %(state_kep_cryo_degrees, state_kep_cryo_calculated_degrees_rounded))


state_cart_iss_latex = latex_utils.Cart2Latex(state_cart_iss)
state_cart_iss_calculated_rounded_latex = latex_utils.Cart2Latex(state_cart_iss_calculated_rounded)

state_kep_iss_degrees_latex = latex_utils.Kep2Latex(state_kep_iss_degrees)
state_kep_iss_calculated_degrees_rounded_latex = latex_utils.Kep2Latex(state_kep_iss_calculated_degrees_rounded)

state_cart_cryo_latex = latex_utils.Cart2Latex(state_cart_cryo)
state_cart_cryo_calculated_rounded_latex = latex_utils.Cart2Latex(state_cart_cryo_calculated_rounded)

state_kep_cryo_degrees_latex = latex_utils.Kep2Latex(state_kep_cryo_degrees)
state_kep_cryo_calculated_degrees_rounded_latex = latex_utils.Kep2Latex(state_kep_cryo_calculated_degrees_rounded)

dolatex=False
if dolatex:
    print('ISS states, in order [true cartesian, converted cartesian, true keplerian, converted keplerian]',
          state_cart_iss_latex,
          state_cart_iss_calculated_rounded_latex,
          state_kep_iss_degrees_latex,
          state_kep_iss_calculated_degrees_rounded_latex)
    print('cryo states, in order [true cartesian, converted cartesian, true keplerian, converted keplerian]',
          state_cart_cryo_latex,
          state_cart_cryo_calculated_rounded_latex,
          state_kep_cryo_degrees_latex,
          state_kep_cryo_calculated_degrees_rounded_latex)



# for i in state_cart_iss:
#     print(i)
#     print(misc_utils.get_sigfigs(i))
# print('\n')
# for i in state_kep_iss_degrees:
#     print(i)
#     print(misc_utils.get_sigfigs(i))
# print('\n')
# for i in state_cart_cryo:
#     print(i)
#     print(misc_utils.get_sigfigs(i))
# print('\n')
# for i in state_kep_cryo_degrees:
#     print(i)
#     print(misc_utils.get_sigfigs(i))



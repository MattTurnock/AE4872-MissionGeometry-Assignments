from Ass1 import Cart2Kep_utils, Kep2Cart_utils, kep_orbit_utils, misc_utils, latex_utils
import numpy as np

state_1_cartesian = [10157768.1264, -6475997.0091, 2421205.9518, 1099.2953996, 3455.1059240, 4355.0978095]        #in m and m/s
state_2_keplerian = [12269687.5912, 0.004932091570, 109.823277603, 134.625563565, 106.380426142, 301.149932402]   #in m and degrees, with M given

state_1_keplerian = Cart2Kep_utils.Cart2Kep(state_1_cartesian)
state_2_cartesian = Kep2Cart_utils.Kep2Cart(state_2_keplerian, given_angle='M')


sigfigs = 10
state_1_keplerian_rounded = misc_utils.round_list_sigfigs(kep_orbit_utils.Kep_Rad2Deg(state_1_keplerian), 10)
state_2_cartesian_rounded = misc_utils.round_list_sigfigs(state_2_cartesian, 10)

print('State 1 Cartesian (given): %s' %state_1_cartesian)
print('State 1 Keplerian (calculated): %s\n' %state_1_keplerian_rounded)

print('State 2 Keplerian (given): %s' %state_2_keplerian)
print('State 2 Cartesian (calculated): %s' %state_2_cartesian_rounded)

print('\nState 1 [Cartesian, Keplerian]:\n',
      latex_utils.Cart2Latex(state_1_cartesian),
      latex_utils.Kep2Latex(state_1_keplerian_rounded), '\n')
print('State 2 [Keplerian, cartesian]: \n',
      latex_utils.Kep2Latex(state_2_keplerian),
      latex_utils.Cart2Latex(state_2_cartesian_rounded))
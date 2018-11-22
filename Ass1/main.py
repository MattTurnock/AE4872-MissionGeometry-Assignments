#Script performing the wanted transformation

from Ass1 import Cart2Kep_utils, Kep2Cart_utils, kep_orbit_utils, misc_utils, latex_utils
import numpy as np

#Define given states
state_1_cartesian = [10157768.1264, -6475997.0091, 2421205.9518, 1099.2953996, 3455.1059240, 4355.0978095]        #in m and m/s
state_2_keplerian_degrees = [12269687.5912, 0.004932091570, 109.823277603, 134.625563565, 106.380426142, 301.149932402]   #in m and degrees, with M given
state_2_keplerian = kep_orbit_utils.Kep_Deg2Rad(state_2_keplerian_degrees)

#Transform the states
state_1_keplerian = Cart2Kep_utils.Cart2Kep(state_1_cartesian)
state_2_cartesian = Kep2Cart_utils.Kep2Cart(state_2_keplerian, given_angle='M')

#Round the states
sigfigs = 10
state_1_keplerian_rounded = misc_utils.round_list_sigfigs(kep_orbit_utils.Kep_Rad2Deg(state_1_keplerian), 10)
state_2_cartesian_rounded = misc_utils.round_list_sigfigs(state_2_cartesian, 10)

#Print an output for convenience
print('State 1 Cartesian (given): %s' %state_1_cartesian)
print('State 1 Keplerian (calculated): %s\n' %state_1_keplerian_rounded)

print('State 2 Keplerian (given): %s' %state_2_keplerian)
print('State 2 Cartesian (calculated): %s\n' %state_2_cartesian_rounded)

#Print latex matrices if wanted
dolatex = False
if dolatex:
      print('\nState 1 [Cartesian, Keplerian]:\n',
            latex_utils.Cart2Latex(state_1_cartesian),
            latex_utils.Kep2Latex(state_1_keplerian_rounded), '\n')
      print('State 2 [Keplerian, cartesian]: \n',
            latex_utils.Kep2Latex(state_2_keplerian),
            latex_utils.Cart2Latex(state_2_cartesian_rounded))

#Also calculate theta, E and M for both states, for completeness
state_1_theta = state_1_keplerian[-1]
state_1_E = kep_orbit_utils.theta2E(state_1_theta, state_1_keplerian[1])
state_1_M = kep_orbit_utils.E2M(state_1_E, state_1_keplerian[1])

state_2_M = state_2_keplerian[-1]
state_2_E = kep_orbit_utils.M2E(state_2_M, state_2_keplerian[1])
state_2_theta = kep_orbit_utils.E2theta(state_2_E, state_2_keplerian[1])

print('State 1 [theta, E, M] = [%s, %s, %s]\n' %(misc_utils.round2sigfig(np.rad2deg(state_1_theta), 10), misc_utils.round2sigfig(np.rad2deg(state_1_E), 10), misc_utils.round2sigfig(np.rad2deg(state_1_M),10)))
print('State 2 [theta, E, M] = [%s, %s, %s]\n' %(misc_utils.round2sigfig(np.rad2deg(state_2_theta), 10), misc_utils.round2sigfig(np.rad2deg(state_2_E), 10), misc_utils.round2sigfig(np.rad2deg(state_2_M), 10)))

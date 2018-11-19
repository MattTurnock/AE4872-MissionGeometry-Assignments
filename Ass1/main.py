from Ass1 import Cart2Kep_utils, Kep2Cart_utils, kep_orbit_utils

state_1_cartesian = [10157768.1264, -6475997.0091, 2421205.9518, 1099.2953996, 3455.1059240, 4355.0978095]        #in m and m/s
state_2_keplerian = [12269687.5912, 0.004932091570, 109.823277603, 134.625563565, 106.380426142, 301.149932402]   #in m and degrees, with M given

state_1_keplerian = Cart2Kep_utils.Cart2Kep(state_1_cartesian)
state_2_cartesian = Kep2Cart_utils.Kep2Cart(state_2_keplerian, given_angle='M')

print(kep_orbit_utils.Kep_Rad2Deg(state_1_keplerian))
print(state_2_cartesian)
import numpy as np
from astropy import units as u
from poliastro.bodies import Earth
from poliastro.twobody import Orbit
from Ass3.integrator_utils import do_integration
from json_to_dict import constants
from Ass1.Cart2Kep_utils import Cart2Kep
from Ass1.Kep2Cart_utils import Kep2Cart
pi = np.pi

a = 7500*u.km
e = 0.1 *u.one
i = 1E-100*u.deg
Omega = 0*u.deg
omega = 0*u.deg
theta = 0*u.deg

kepstate = [a.to(u.m).value, e.to(u.one).value, i.to(u.rad).value, Omega.to(u.rad).value, omega.to(u.rad).value, theta.to(u.rad).value]
cartsate = Kep2Cart(kepstate)

X0 = cartsate
t0 = 0
t_end = ((1*u.wk).to(u.s)).value
ts = np.linspace(t0, t_end, 10000)

Xall, Xdots = do_integration(X0, ts, method="RK4")
print(Xall, Xdots)
cartstate_end = Xall[-1, 1:]
kepstate_end = Cart2Kep(cartstate_end, do_units_out=True)
print(cartstate_end)
print(kepstate_end)




init_orbit = Orbit.from_classical(Earth, a, e, i, Omega, omega, theta)
init_pos = init_orbit.state.r
init_vel = init_orbit.state.v
print(init_pos, init_vel)
print(init_orbit)












# a = 7096137.00*u.m
# e = 0.0011219*u.one
# i = 92.0316*u.deg
# Omega = 296.1384*u.deg
# omega = 120.6878*u.deg
# theta = 239.5437*u.deg
# M = 239.6546*u.deg
# theta = pl.twobody.angles.M_to_nu(M, e)
# print(theta + 360*u.deg)
#
#
# # a = 7500*u.km
# # e = 0.1 *u.one
# # i = 0*u.deg
# # Omega = 0*u.deg
# # omega = 0*u.deg
# # nu = 0*u.deg
#

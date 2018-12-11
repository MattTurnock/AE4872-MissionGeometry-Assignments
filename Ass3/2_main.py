#Does calculations for INTEG-2
####################################################################################
import numpy as np
from astropy import units as u
from Ass3.integrator_utils import get_integrator_table_orbit
from matplotlib import pyplot as plt
pi = np.pi
np.set_printoptions(suppress=True)

a = 7500*u.km
a_true = 7500*1000
e = 0.1 *u.one
e_true = 0.1
i = 1E-100*u.deg
Omega = 0*u.deg
omega = 0*u.deg
theta = 0*u.deg
kepstate_0 = [a.to(u.m).value, e.to(u.one).value, i.to(u.rad).value, Omega.to(u.rad).value, omega.to(u.rad).value, theta.to(u.rad).value]

t0 = 0
t_end = ((1*u.wk).to(u.s)).value
dts = [0.1, 1, 10, 100, 1000, 10000, 100000]
###############################################################################################################################################
# doing methods
calc=True
printing=True
if calc:
    euler_table = get_integrator_table_orbit(kepstate_0, t_end, dts, t0=0, method="euler", printing=printing)
    np.save("orbit_error_table_euler", euler_table)

    RK4_table = get_integrator_table_orbit(kepstate_0, t_end, dts, t0=0, method="RK4", printing=printing)
    np.save("orbit_error_table_RK4", RK4_table)

euler_table = np.load("orbit_error_table_euler.npy")
euler_table_tosave = np.copy(euler_table)
euler_table_tosave[:,1] = euler_table_tosave[:,1]/10**3
euler_table_tosave[:,3] = euler_table_tosave[:,3]/10**6
np.savetxt('orbit_error_table_euler.txt', euler_table_tosave, delimiter=' & ', newline=' \\\\\n \hline \n', fmt='%1.E & %1.f & %1.2f & %1.2f & %1.3E & %1.3E')
RK4_table = np.load("orbit_error_table_RK4.npy")
RK4_table_tosave = np.copy(RK4_table)
RK4_table_tosave[:,1] = RK4_table_tosave[:,1]/10**6
RK4_table_tosave[:,3] = RK4_table_tosave[:,3]/10**3
np.savetxt('orbit_error_table_RK4.txt', RK4_table_tosave, delimiter=' & ', newline=' \\\\\n \hline \n', fmt='%1.E & %1.3f & %1.2f & %1.3E & %1.3E & %1.3E')

print(RK4_table, '\n')
print(euler_table)

euler_dist_error =abs(euler_table[:,3])/1000
euler_calc_no = euler_table[:,-1]

RK4_dist_error = abs(RK4_table[:,3])/1000
RK4_calc_no = RK4_table[:,-1]

plot=True
plt.figure()
plt.loglog(RK4_calc_no, RK4_dist_error)
plt.loglog(euler_calc_no, euler_dist_error)
plt.legend(["RK4", "Euler"])
plt.ylabel("Final semi-major axis error [km]")
plt.xlabel("Number of derivative evaluations, N")
plt.grid()
plt.savefig("orbit_errors.pdf")
if plot: plt.show()

from Ass3.integrator_utils import do_integration, do_err_plot, get_fulldata_car, get_integrator_table_car, get_A_car
import numpy as np
from astropy import units as u

from matplotlib import pyplot as plt
pi = np.pi
np.set_printoptions(suppress=True)

x_end_true = 3600
v_end_true = 120

###############################################################################
#Euler
print('Basic euler (with dt=1s) end:')
V0 = 0.0*u.m/u.s
D0 = 0.0*u.m
t0 = 0*u.s
t_end = 60*u.s
a0 = 2*u.m/u.s**2
dt = 1*u.s
a0 = a0.to(u.m/u.s**2).value
V0 = V0.to(u.m/u.s).value
D0 = D0.to(u.m).value
t0 = t0.to(u.s).value
t_end = t_end.to(u.s).value
dt = dt.to(u.s).value

X0 = np.array([D0,
               V0])


B = np.array([0,
              a0])
ts = np.linspace(t0, t_end, 100)


Xall_e, Xdots_e = do_integration(X0, ts, get_A=get_A_car, B=B, method="euler")
fulldata_e = get_fulldata_car(Xall_e, Xdots_e)


#has cols: [t, x, v, a]
print(fulldata_e[-1], '\n')
print('Doing a full euler dt test. Following array has columns [dt, x_end, v_end, error x_end, error v_end, number of evaluations]')

dts = [0.001, 0.01, 0.1, 1, 10, 20, 21]
tablen=6

calc=False
if calc:
    euler_table = get_integrator_table_car(X0, dts, t_end, tablen, B, x_end_true, v_end_true, method="euler")
    np.save('car_error_table_euler', euler_table)

euler_table = np.load('car_error_table_euler.npy')
np.savetxt('car_error_table_euler.txt', euler_table, delimiter=' & ', newline=' \\\\\n \hline \n', fmt='%1.1E & %1.2f & %1.f & %1.2E& %1.f & %1.3E')

print(euler_table, '\n')

err_e = euler_table[:, 3]




title = "Euler car error"
scale="linear"
#do_err_plot(dts, err_e, show=False, save=False, method="euler", figure=True, scale=scale)


###############################################################################
# RK4
print('Basic RK4 (with dt=1s) end:')
V0 = 0.0*u.m/u.s
D0 = 0.0*u.m
t0 = 0*u.s
t_end = 60*u.s
a0 = 2*u.m/u.s**2
dt = 1*u.s
a0 = a0.to(u.m/u.s**2).value
V0 = V0.to(u.m/u.s).value
D0 = D0.to(u.m).value
t0 = t0.to(u.s).value
t_end = t_end.to(u.s).value
dt = dt.to(u.s).value

X0 = np.array([D0,
               V0])


B = np.array([0,
              a0])


Xall_rk, Xdots_rk = do_integration(X0, ts, get_A=get_A_car, B=B, method="RK4")
fulldata_rk = get_fulldata_car(Xall_rk, Xdots_rk)


#has cols: [t, x, v, a]
print(fulldata_rk[-1], '\n')
print('Doing a full RK4 dt test. Following array has columns [dt, x_end, v_end, error x_end, error v_end, number of evaluations]')

dts = [0.001, 0.01, 0.1, 1, 10, 20, 21]
tablen=6

if calc:
    RK_table = get_integrator_table_car(X0, dts, t_end, tablen, B, x_end_true, v_end_true, method="RK4")
    np.save('car_error_table_RK4', RK_table)



RK_table = np.load('car_error_table_RK4.npy')
np.savetxt('car_error_table_RK4.txt', RK_table, delimiter=' & ', newline=' \\\\\n \hline \n', fmt='%1.1E & %1.2f & %1.f & %1.2E& %1.f & %1.3E')
#np.savetxt('test.txt', RK_table)
print(RK_table)

err_rk = RK_table[:, 3]
title = None
#do_err_plot(dts, err_rk, show=True, save=True, method="both", figure=False, scale=scale, legend=["Euler", "RK4"])


euler_dist_error = abs(euler_table[:,3])
euler_calc_no = euler_table[:,-1]

RK4_dist_error = abs(RK_table[:,3])
RK4_calc_no = RK_table[:,-1]

plt.figure()
plt.plot(RK4_calc_no, RK4_dist_error, linewidth=3)
plt.plot(euler_calc_no, euler_dist_error)

plt.ylim([0, 10])
plt.legend(["RK4", "Euler"])
plt.ylabel("Final distance eror [m]")
plt.xlabel("Number of derivative evaluations, N")
plt.grid()
#plt.savefig("car_errors.pdf")
#plt.show()


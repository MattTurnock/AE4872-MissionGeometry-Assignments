#Verification script for integrator functions
####################################################################################3
from Ass3.integrator_utils import do_integration
import numpy as np
from astropy import units as u
from Ass1.Kep2Cart_utils import Kep2Cart
pi = np.pi
np.set_printoptions(suppress=True)

a = 7378.137*u.km
e = 0.1 *u.one
i = 1E-100*u.deg
Omega = 0*u.deg
omega = 0*u.deg
theta = 0*u.deg

kepstate = [a.to(u.m).value, e.to(u.one).value, i.to(u.rad).value, Omega.to(u.rad).value, omega.to(u.rad).value, theta.to(u.rad).value]
cartsate = Kep2Cart(kepstate)

X0 = cartsate
dt=10
ts = np.arange(0,30,dt)



##############################################################################################################################
#euler verification

Xall, Xdots = do_integration(X0, ts, method="euler")
print('==========================================================Euler===========================================================\n')
print('To mimic what is in the slides, from 0s to 10s:')

X0 = Xall[0][1:]
X0dot = Xdots[0][1:]
X0new = Xall[1][1:]

DATA_0to10 = np.zeros([6,4])
DATA_0to10[:,0] = np.arange(1,7)
DATA_0to10[:,1] = X0
DATA_0to10[:,2] = X0dot
DATA_0to10[:,3] = X0new
np.savetxt('verification_euler_step1.txt', DATA_0to10, delimiter=' & ', newline=' \\\\\n \hline \n', fmt='%i & %1.2f & %1.2f & %1.2f')
print(DATA_0to10, '\n')

print('And from 10s to 20s:')
X10 = Xall[1][1:]
X10dot = Xdots[1][1:]
X10new = Xall[2][1:]

DATA_10to20 = np.zeros([6,4])
DATA_10to20[:,0] = np.arange(1,7)
DATA_10to20[:,1] = X10
DATA_10to20[:,2] = X10dot
DATA_10to20[:,3] = X10new
np.savetxt('verification_euler_step2.txt', DATA_10to20, delimiter=' & ', newline=' \\\\\n \hline \n', fmt='%i & %1.2f & %1.2f & %1.2f')
print(DATA_10to20, '\n')
print('==========================================================RK4===========================================================\n')

print('Again mimicing the slides, from 0s to 10s:')
Xall, Xdots, kall = do_integration(X0, ts, method="RK4", return_ks=True)
ks = kall[0]

X0 = Xall[0][1:]
k1 = ks[0]
X0pk1 = X0 + 0.5 * dt * k1
k2 = ks[1]
X0pk2 = X0 + 0.5 * dt * k2
k3 = ks[2]
X0pk3 = X0 + dt * k3
k4 = ks[3]
X0new = Xall[1][1:]
to_add = np.array([
                    X0,
                    k1,
                    X0pk1,
                    k2,
                    X0pk2,
                    k3,
                    X0pk3,
                    k4,
                    X0new
                    ])

DATA_0to10 = np.zeros([6,10])
DATA_0to10[:,0] = np.arange(1,7)
DATA_0to10[:,1] = X0
DATA_0to10[:,2] = k1
DATA_0to10[:,3] = X0pk1
DATA_0to10[:,4] = k2
DATA_0to10[:,5] = X0pk2
DATA_0to10[:,6] = k3
DATA_0to10[:,7] = X0pk3
DATA_0to10[:,8] = k4
DATA_0to10[:,9] = X0new
print(ks)
print(DATA_0to10, '\n')
np.savetxt('verification_RK4_step1.txt', DATA_0to10, delimiter=' & ', newline=' \\\\\n \hline \n', fmt='%i & %1.2f & %1.2f & %1.2f & %1.2f & %1.2f & %1.2f & %1.2f & %1.2f & %1.2f')


print('And from 10s to 20s:')
Xall, Xdots, kall = do_integration(X0, ts, method="RK4", return_ks=True)
ks = kall[1]

X10 = Xall[1][1:]
k1 = ks[0]
X10pk1 = X10 + 0.5 * dt * k1
k2 = ks[1]
X10pk2 = X10 + 0.5 * dt * k2
k3 = ks[2]
X10pk3 = X10 + dt * k3
k4 = ks[3]
X10new = Xall[2][1:]
to_add = np.array([
                    X10,
                    k1,
                    X10pk1,
                    k2,
                    X10pk2,
                    k3,
                    X10pk3,
                    k4,
                    X10new
                    ])

DATA_10to20 = np.zeros([6,10])
DATA_10to20[:,0] = np.arange(1,7)
DATA_10to20[:,1] = X10
DATA_10to20[:,2] = k1
DATA_10to20[:,3] = X10pk1
DATA_10to20[:,4] = k2
DATA_10to20[:,5] = X10pk2
DATA_10to20[:,6] = k3
DATA_10to20[:,7] = X10pk3
DATA_10to20[:,8] = k4
DATA_10to20[:,9] = X10new
print(ks)
print(DATA_10to20)
np.savetxt('verification_RK4_step2.txt', DATA_10to20, delimiter=' & ', newline=' \\\\\n \hline \n', fmt='%i & %1.2f & %1.2f & %1.2f & %1.2f & %1.2f & %1.2f & %1.2f & %1.2f & %1.2f')





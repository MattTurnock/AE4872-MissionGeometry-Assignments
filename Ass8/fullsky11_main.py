#!/usr/bin/env python

"""
fullsky11_main.py: Main script for the FULLSKY-11 assignment
"""
__author__      = "Matthew Turnock"
__email__ = "matthew.turnock@protonmail.com"
__version__ = "1.0"

########################################################################################################################

from Ass8.fullsky11_utils import *
from Ass6.design7_utils import *
from Ass1.kep_orbit_utils import *
from json_to_dict import constants

# Delfi constants
TDelfi = constants["TDelfi-n3Xt"]
iDelfi = constants["iDelfi-n3Xt"]
nDelfi = getn(TDelfi, outputUnits=u.rad/u.s)

# Set problem parameters
rho1 = iDelfi.to(u.deg)
rho2 = 90*u.deg
omega1 = -((2*pi*u.rad)/(24*u.hour)).to(u.rad/u.s, equivalencies=u.dimensionless_angles())
omega2 = nDelfi
phi10 = 0*u.rad
phi20 = 0*u.rad

ts = np.linspace(0.01*u.s,TDelfi - 0.01*u.s, 100)

# Open list to put data, has order: delta, Deltaalpha, alpha, deltaEPrime, rhoE, omegaE, V, DeltaPsi, Psi, phi1, phi2
delfiData = np.zeros( (len(ts), 11), dtype=object )

for i in range(len(ts)):
    t = (ts[i]).to(u.s)
    # print(rho1, rho2, omega1, omega2, phi10, phi20)
    # print(t)

    # Gives values of chosen units
    delfiData[i, :] = doDualAxisSpiral(rho1, rho2, omega1, omega2, phi10, phi20, t, angleUnits=u.deg, omegaUnits=u.rad/u.s, returnValue=True)

params = ["delta", "Deltaalpha", "alpha", "deltaEPrime", "rhoE", "omegaE", "V", "DeltaPsi", "Psi", "phi1", "phi2"]
ylabels = [r"$\delta$ [deg]", r"$\Delta \alpha$ [deg]", r"$\alpha$ [deg]", r"$\delta_E'$ [deg]", r"$\rho_E$ [deg]", r"$\omega_E$ [rad/s]", r"$V$ [rad/s]", r"$\Delta \Psi$ [deg]", r"$\Psi$ [deg]", r"$\phi_1$ [deg]", r"$\phi_2$ [deg]"]

saveBase = "plots/%sPlot.png"

for i in range(len(delfiData[0,:])):
    plt.figure()
    yData = delfiData[:, i]
    plt.plot(ts, yData)
    # plt.ylim((min(yData), max(yData)))

    plt.ylabel(ylabels[i])
    plt.xlabel("Time [mins]")
    plt.grid()
    plt.savefig(saveBase %(params[i]))
plt.show()


# thi = doDualAxisSpiral(rho1, rho2, omega1, omega2, phi10, phi20, 1.2*u.s, angleUnits=u.deg, omegaUnits=u.rad/u.s)
#
# print(thi)




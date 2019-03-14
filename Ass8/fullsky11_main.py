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
rho1 = iDelfi
rho2 = 90*u.deg
omega1 = -((2*pi*u.rad)/(24*u.hour)).to(u.rad/u.s, equivalencies=u.dimensionless_angles())
print("omega_1 = %s" %omega1)
omega2 = nDelfi
print("omega_2 = %s" %omega2)
phi10 = 0*u.rad
phi20 = 0*u.rad
N = 1000
ts = np.linspace(0.01*u.s, TDelfi - 10*TDelfi/N, N)
dt = ts[3]-ts[2]
print("timestep dt = %s" %dt.to(u.s))

# Set some general things
showing = True

# Open list to put data, has order: delta, Deltaalpha, alpha, deltaEPrime, rhoE, omegaE, V, DeltaPsi, Psi, phi1, phi2
delfiData = np.zeros( (len(ts), 11), dtype=object )

# Simulate the orbit
for i in range(len(ts)):
    t = (ts[i]).to(u.s)
    # Gives values of chosen units
    delfiData[i, :] = doDualAxisSpiral(rho1, rho2, omega1, omega2, phi10, phi20, t, angleUnits=u.deg, omegaUnits=u.rad/u.s, returnValue=True)

# Do relevant plots

params = ["delta", "Deltaalpha", "alpha", "deltaEPrime", "rhoE", "omegaE", "V", "DeltaPsi", "Psi", "phi1", "phi2"]
ylabels = [r"$\delta$ [deg]", r"$\Delta \alpha$ [deg]", r"$\alpha$ [deg]", r"$\delta_E'$ [deg]", r"$\rho_E$ [deg]", r"$\omega_E$ [rad/s]", r"$V$ [rad/s]", r"$\Delta \Psi$ [deg]", r"$\Psi$ [deg]", r"$\phi_1$ [deg]", r"$\phi_2$ [deg]"]

saveBase = "plots/%sPlot.%s"

for i in range(len(delfiData[0,:])):
    fig = plt.figure()
    yData = delfiData[:, i]
    plt.plot(ts, yData)
    plt.ylabel(ylabels[i])
    plt.xlabel("Time [mins]")
    plt.grid()
    plt.savefig(saveBase %(params[i], "png"), bbox_inches="tight")
    plt.savefig(saveBase % (params[i], "pdf"), bbox_inches="tight")


# Plot alpha (x) against delta (y)
plt.figure()
plt.plot(delfiData[:, 2], delfiData[:, 0])
plt.xlabel(ylabels[2])
plt.ylabel(ylabels[0])
plt.grid()
plt.savefig(saveBase % ("deltaVSalpha", "png"), bbox_inches="tight")
plt.savefig(saveBase % ("deltaVSalpha", "pdf"), bbox_inches="tight")

if showing: plt.show()

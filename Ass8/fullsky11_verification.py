#!/usr/bin/env python

"""
fullsky11_verification.py: Benchmarking and verification for the FULLSKY-11 assignment
"""
__author__      = "Matthew Turnock"
__email__ = "matthew.turnock@protonmail.com"
__version__ = "1.0"

########################################################################################################################
from Ass8.fullsky11_utils import *
from Ass6.design7_utils import *

rho1 = 40*u.deg
rho2 = 20*u.deg

phi1 = 5*u.deg
phi2s = [90*u.deg, 90*u.deg, 100*u.deg, 300*u.deg]

omega1s = [0*u.rad/u.s, 1*u.rad/u.s, 1*u.rad/u.s, 1*u.rad/u.s]
omega2 = 3*u.rad/u.s

veriArray = np.zeros( (9, 4), dtype=object)

for i in range(len(phi2s)):
    phi2 = phi2s[i]
    omega1 = omega1s[i]
    outputs = doPartialDualAxisSpiral(rho1, rho2, omega1, omega2, phi1, phi2,
                                      angleUnits=u.deg, omegaUnits=u.rad/u.s)
    veriArray[:, i] = outputs

veriArrayToSave = np.zeros( (9,5), dtype=float)

for i in range(len(veriArray[:, 0])):
    for j in range(len(veriArray[0, :])):

        val = veriArray[i, j]
        val = val.value

        if (i == 5) or (i == 6):
            decimals = 3
        else:
            decimals = 2

        val = str(np.round(val, decimals))
        veriArrayToSave[i, j+1] = val

print("Verification table saved as verificationLatex.txt")
npArray2LatexTable(veriArrayToSave, "verificationLatex.txt")













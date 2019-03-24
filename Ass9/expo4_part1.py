#!/usr/bin/env python

"""
expo4_part1.py: Script to do part 1 of the EXPO-4 assignment
"""
__author__      = "Matthew Turnock"
__email__ = "matthew.turnock@protonmail.com"
__version__ = "1.0"

########################################################################################################################
import numpy as np
from matplotlib import pyplot as plt

from Ass6.design7_utils import npArray2LatexTable
from Ass9.expo4_utils import *
from json_to_dict import constants

calculating = False
showing = True
printing = True



##########################################################################################################
# Set problem  parameters

muSun = constants["muSun"].decompose().value
k2 = 1/12
r1 = (1 * u.AU).decompose().value
r2 = (1.5 * u.AU).decompose().value
Psi = (90*u.deg).decompose().value
TOFReq = 2.0


# Set number of gammas to calculate, and number of sub-intervals to use and N list
numberOfGammas = 1000
numberOfSteps = [5000]
Ns = [0, 1, 2, 3]

TOFGammaDataSavenameBase = "Q1Data/TOFGammaData_%s.npy"
TOFGammaDatas = []
labels=[]
gamma1sArr = np.zeros(( numberOfGammas, len(Ns) ))
for i in range(len(Ns)):
    N = Ns[i]
    print("\nDoing N : %s\n" %N)
    # Calculate gamma max and min
    thetaBar = getthetaBar(Psi, N)
    Delta = getDelta(r1, r2, k2, thetaBar)
    gamma1Min = getgamma1Min(r1, r2, k2, thetaBar, Delta, outputUnits=None)
    gamma1Max = getgamma1Max(r1, r2, k2, thetaBar, Delta)

    # Make list of gammas
    gamma1s = np.linspace(gamma1Min, gamma1Max, numberOfGammas)
    gamma1sArr[:, i] = gamma1s


    TOFGammaDataSavename = TOFGammaDataSavenameBase %N
    if calculating:
        doTOFCalcMain(gamma1s, numberOfSteps, k2, r1, r2, Psi, N, muSun, saveName=TOFGammaDataSavename)

    TOFGammaData = np.load(TOFGammaDataSavename)
    TOFGammaDatas.append(TOFGammaData)
    ClosestIndex = find_nearest_index(TOFGammaData, TOFReq)
    print("Closest TOF at gamma = %s, TOF = %s" %(gamma1s[ClosestIndex], TOFGammaData[ClosestIndex]))


    labels.append(r"N=%s: $\gamma_m$=%s, $\gamma_M$=%s" %(N, np.round(gamma1Min, decimals=4 ), np.round(gamma1Max, decimals=4 )))

xlabel = r"$\gamma_1$ [deg]"
ylabel = "TOF [yrs]"
xlims = [-1.5, 1]
ylims= [0, 4]

plt.figure()
for i in range(len(Ns)):
    plt.plot(gamma1sArr[:, i], TOFGammaDatas[i])
plt.plot(xlims, [TOFReq, TOFReq], linestyle="--")
plt.xlim(xlims)
plt.ylim(ylims)
plt.legend(labels)
plt.xlabel(xlabel)
plt.ylabel(ylabel)
plt.grid()
plt.savefig("Q1Data/TOFPlot.pdf")
if showing: plt.show()






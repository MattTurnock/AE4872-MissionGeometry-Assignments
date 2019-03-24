#!/usr/bin/env python

"""
expo4_verification.py: Verification script for the EXPO-4 assignment
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

calculating = True
showing = True
printing = True

##########################################################################################################
# Do example case, with N=2 etc

muSun = constants["muSun"].decompose().value

r1 = (1 * u.AU).decompose().value
r2 = (1.5 * u.AU).decompose().value
k2 = 1/12
Psi = (90*u.deg).decompose().value
N = 2

gamma1Example = (-80*u.deg).decompose().value

dataExample = doFullExposin(k2, r1, r2, Psi, N, gamma1Example, muSun)

dataToSave1 = np.zeros( (len(dataExample), 6) )
dataToSave1[:,0] = dataExample[:,0]
dataToSave2 = np.copy(dataToSave1)

dataToSave1[:, 1:] = dataExample[:, 1:6]
dataToSave2[:, 1:] = dataExample[:, 6:]

formatList1 = ["%.7f", "%.3f", "%.3f", "%.3f", "%.8f", "%.7f"]
formatList2 = ["%.7f", "%.6E", "%.6E", "%.9E", "%.8E", "%.9E"]
formatListMaster = ["%.7f", "%.3f", "%.3f", "%.3f", "%.8f", "%.7f", "%.6E", "%.6E", "%.9E", "%.8E", "%.9E"]

npArray2LatexTable(dataExample, "verificationData/dataAll.txt", fmt=formatListMaster)
npArray2LatexTable(dataToSave1, "verificationData/dataPart1.txt", fmt=formatList1)
npArray2LatexTable(dataToSave2, "verificationData/dataPart2.txt", fmt=formatList2)

########## Get TOF for example case with different numbers of evals ################

stepNumbers = [10, 100, 1000, 2000, 10000]
TOFs = np.zeros(len(stepNumbers))

if calculating:

    for i in range(len(stepNumbers)):
        stepNumber = stepNumbers[i]

        dataStep = doFullExposin(k2, r1, r2, Psi, N, gamma1Example, muSun, thetaSplitNumber=stepNumber)
        TOF = (dataStep[-1, -1] * u.s).to(u.year).value
        TOFs[i] = TOF

    np.save("verificationData/TOFData.npy", TOFs)

else:
    TOFs = np.load("verificationData/TOFData.npy")

for i in range(len(TOFs)):
    TOF = TOFs[i]
    stepNumber = stepNumbers[i]
    print("Steps %s  : %s" % (stepNumber, TOF))
########## Do TOF plots for many gammas ################
# Get gamma1min and gamma1max
thetaBar = getthetaBar(Psi, N)
Delta = getDelta(r1, r2, k2, thetaBar)
gamma1Min = getgamma1Min(r1, r2, k2, thetaBar, Delta, outputUnits=None)
gamma1Max = getgamma1Max(r1, r2, k2, thetaBar, Delta)
gamma1Max = np.deg2rad(5)

gamma1s = np.linspace(gamma1Min, gamma1Max, 100)



# stepNumbers = [10, 100]


if calculating:
    doTOFCalcMain(gamma1s, stepNumbers, k2, r1, r2, Psi, N, muSun, saveName="verificationData/TOFGammaData.npy")

TOFGammaData = np.load("verificationData/TOFGammaData.npy")



gamma1s_degrees = np.rad2deg(gamma1s)

labels = ["10 steps", "100 steps", "1000 steps", "2000 steps", "10000 steps"]
xlabel = r"$\gamma_1$ [deg]"
ylabel = "TOF [yrs]"
savenameBase = "verificationData/verificationPlot_%s.pdf"

plt.figure()
for i in range(len(stepNumbers)):
    plt.plot(gamma1s_degrees, TOFGammaData[:, i])
plt.xlim([-90, 90])
plt.ylim([0, 5])
plt.legend(labels)
plt.xlabel(xlabel)
plt.ylabel(ylabel)
plt.grid()
plt.savefig(savenameBase %"main")

plt.figure()
for i in range(len(stepNumbers)):
    plt.plot(gamma1s_degrees, TOFGammaData[:, i])

plt.xlim([-90, -30])
plt.ylim([0, 0.5])
plt.legend(labels)
plt.xlabel(xlabel)
plt.ylabel(ylabel)
plt.grid()
plt.savefig(savenameBase %"cropped")
if showing: plt.show()




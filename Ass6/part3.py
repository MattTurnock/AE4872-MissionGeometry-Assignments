#!/usr/bin/env python

"""
part3.py: Script to perform part 3 of the DESIGN-7 assignment
"""
__author__      = "Matthew Turnock"
__email__ = "matthew.turnock@protonmail.com"
__version__ = "1.0"

########################################################################################################################

import numpy as np
from Ass6.design7_utils import *
from Ass6.part1 import *
from Ass6.part2 import *
from matplotlib import pyplot as plt

# np.random.seed(12)

iterations = 100000
calculating = False
fullLoad = False
finalLoad = True

monteCarloResults = np.zeros((2,2), dtype=object)
# monteCarloResults

print("\n==========================================Monte Carlo===================================================================")

#################### Do a run for the case with no systematic error #########################################

# Define random parameter indexes and values
randIndices = np.array([1, 4, 7, 9, 11, 13, 15, 17, 19, 21])
randSigmas = np.zeros(len(randIndices), dtype=object)
randSigmas[:] = [problemParameters.orbitDetermination,
              problemParameters.timing,
              problemParameters.starSensorMeasurement,
              problemParameters.starSensorMounting,
              problemParameters.starCatalog,
              problemParameters.attitudeComp,
              problemParameters.payloadSensorMeasurement,
              problemParameters.targetCentroiding,
              problemParameters.payloadSensorMounting,
              problemParameters.coordTransformation]

# Define systematic parameter indexes and values
systIndices = np.array([])
systValues = np.array([])

savename = "monteCarloArrayRand.npy"
savenameLastCol = "monteCarloArrayRandFinalColOnly.npy"
if calculating:
    monteCarloArrayRand = doMonteCarlo(iterations, randIndices, randSigmas, Xtarg_nom=Xtarg_nom, Vsc=Vsc, eta=eta, h=h, systIndices=systIndices, systValues=systValues)
    monteCarloArrayRandFinalCol = monteCarloArrayRand[:, 26]
    np.save(savename, monteCarloArrayRand)
    np.save(savenameLastCol, monteCarloArrayRandFinalCol)

if fullLoad:
    monteCarloArrayRand = np.load(savename)
    monteCarloArrayRandFinalCol = monteCarloArrayRand[:, 26]
if finalLoad:
    monteCarloArrayRandFinalCol = np.load(savenameLastCol)

RMS = np.sqrt(np.mean(np.square(monteCarloArrayRandFinalCol)))
mean = np.mean(monteCarloArrayRandFinalCol)
monteCarloResults[0, :] = mean, RMS
if printing:
    print("\nFor %s case:" % "random only")
    print("RMS = %s \nMean = %s" % (RMS, mean))

#################### Do a run for the case with systematic error #########################################

# Define random parameter indexes and values
randIndices = np.array([1, 4, 7, 9, 11, 13, 15, 17, 19, 21])
randSigmas = np.zeros(len(randIndices), dtype=object)
randSigmas[:] = [problemParameters.orbitDetermination,
              problemParameters.timing,
              problemParameters.starSensorMeasurement,
              problemParameters.starSensorMounting,
              problemParameters.starCatalog,
              problemParameters.attitudeComp,
              problemParameters.payloadSensorMeasurement,
              problemParameters.targetCentroiding,
              problemParameters.payloadSensorMounting,
              problemParameters.coordTransformation]

# Define systematic parameter indexes and values
systIndices = np.array([18])
systValues = np.array([problemParameters.payloadSensorMountingSystematic], dtype=object)

savename = "monteCarloArrayBoth.npy"
savenameLastCol = "monteCarloArrayBothFinalColOnly.npy"
if calculating:
    monteCarloArrayBoth = doMonteCarlo(iterations, randIndices, randSigmas, Xtarg_nom=Xtarg_nom, Vsc=Vsc, eta=eta, h=h, systIndices=systIndices, systValues=systValues)
    monteCarloArrayBothFinalCol = monteCarloArrayBoth[:, 26]
    np.save(savename, monteCarloArrayBoth)
    np.save(savenameLastCol, monteCarloArrayBothFinalCol)

if fullLoad:
    monteCarloArrayBoth = np.load(savename)
    monteCarloArrayBothFinalCol = monteCarloArrayBoth[:, 26]

if finalLoad:
    monteCarloArrayBothFinalCol = np.load(savenameLastCol)


RMS = np.sqrt(np.mean(np.square(monteCarloArrayBothFinalCol)))
mean = np.mean(monteCarloArrayBothFinalCol)
monteCarloResults[1, :] = mean, RMS
if printing:
    print("\nFor %s case:" % "random + systematic")
    print("RMS = %s \nMean = %s" % (RMS, mean))



np.savetxt("monteCarloResults.txt", monteCarloResults, fmt="%s", delimiter="      ", header="row0 = rand only.    row1 = rand+sys.    col0 = mean.    col1 = RMS")


monteCarloArrayRandFinalColValues = []
for i in range(len(monteCarloArrayRandFinalCol)):
    monteCarloArrayRandFinalColValues.append(monteCarloArrayRandFinalCol[i].value)

# plt.hist(monteCarloArrayRandFinalColValues, bins=1000, density=True, histtype="step")
# plt.show()k
#!/usr/bin/env python

"""
part2.py: Script to perform part 2 of the DESIGN-7 assignment
"""
__author__      = "Matthew Turnock"
__email__ = "matthew.turnock@protonmail.com"
__version__ = "1.0"

########################################################################################################################
from Ass6.part1 import *

nadirErrorSystematic = [[problemParameters.payloadSensorMountingSystematic,
           "PAYLOAD SENSOR MOUNTING (SYSTEMATIC):",
           nadMapError2]]

payloadSensorMountingSystematicMapping = nadMapError(nadirErrorSystematic[0][0], D0, epsilon, units=u.m)

nadirErrorSystematic[0].append(payloadSensorMountingSystematicMapping)
nadirErrorSystematic = np.array(nadirErrorSystematic)

nadirErrorSystematicToSave = np.array(nadirErrorSystematic)
nadirErrorSystematicToSave[0][2] = nadirErrorSystematicToSave[0][2].__name__

npArray2LatexTable(nadirErrorSystematicToSave, "nadirSystematicErrors.txt")

nadirRSSSystematic = doRSS(nadirErrorSystematic[:, 3])
# print(nadirRSSSystematic)

nadirRSSTotal = doRSS([nadirRSSRandom, nadirRSSSystematic])

if printing:
    print("Nadir RSS based on (my) calculated values for systematic = %s \n"
          "Nadir RSS for both random and systematic                 = %s" %(nadirRSSSystematic, nadirRSSTotal))


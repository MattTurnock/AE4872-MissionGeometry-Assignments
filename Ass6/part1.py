#!/usr/bin/env python

"""
part1.py: Script to perform part 1 of the DESIGN-7 assignment
"""
__author__      = "Matthew Turnock"
__email__ = "matthew.turnock@protonmail.com"
__version__ = "1.0"

########################################################################################################################
from json_to_dict import constants
import numpy as np
from astropy import units as u
from Ass6.design7_utils import *
import sys

pi=np.pi

saveOutput = True
printing = True

if saveOutput:
    orig_stdout = sys.stdout
    f = open('analyticalPrintout.txt', 'w')
    sys.stdout = f

########################################################################################################################
# Find common parameters

problemParameters = ProblemParameters()
D0 = problemParameters.D0
epsilon = problemParameters.epsilon
eta = problemParameters.eta
h = problemParameters.h
Vsc = problemParameters.Vsc
Xtarg_nom = problemParameters.Xtarg_nom

decimals = 1

########################################################################################################################
# Determine analytical relations
basePrint = "Error                               = %s \n" \
            "Known impact on mapping budget      = %s \n" \
            "Using the relation:                   %s \n" \
            "Calculated impact on mapping budget = %s\n" \
            "Calculated impact (rounded)         = %s\n" \
            "The error is:                         %s\n\n" \
            "========================================================================"


##################################


# List of things to loop through when printing mcgubbins
toLoop = [[problemParameters.starSensorMeasurement,
           problemParameters.starSensorMeasurementMapping,
           "STAR SENSOR MEASUREMENT:",
           nadMapError2],
          [problemParameters.starSensorMounting,
           problemParameters.starSensorMountingMapping,
           "STAR SENSOR MOUNTING:",
           nadMapError2],
          [problemParameters.starCatalog,
           problemParameters.starCatalogMapping,
           "STAR CATALOG ACCURACY:",
           nadMapError2],
          [problemParameters.attitudeComp,
           problemParameters.attitudeCompMapping,
           "ATTITUDE COMPUTATION:",
           nadMapError2],
          [problemParameters.payloadSensorMeasurement,
           problemParameters.payloadSensorMeasurementMapping,
           "PAYLOAD SENSOR MEASUREMENT:",
           nadMapError2],
          [problemParameters.targetCentroiding,
           problemParameters.targetCentroidingMapping,
           "TARGET CENTROIDING:",
           nadMapError2],
          [problemParameters.payloadSensorMounting,
           problemParameters.payloadSensorMountingMapping,
           "PAYLOAD SENSOR MOUNTING:",
           nadMapError2],
          [problemParameters.coordTransformation,
           problemParameters.coordTransformationMapping,
           "TRANSFORMATION OF TARGET LOCATION TO INERTIAL COORDINATES:",
           nadMapError2],
          [problemParameters.orbitDetermination,
           problemParameters.orbitDeterminationMapping,
           "ORBIT DETERMINATION:",
           None],
          [problemParameters.timing,
          problemParameters.timingMapping,
           "TIMING:",
           timingMapError2],
          [problemParameters.projection,
           problemParameters.projectionMapping,
           "PROJECTION ERROR DUE TO ALTITUDE:",
           None],
          [problemParameters.subsatellitePoint,
           problemParameters.subsatellitePointMapping,
           "SUBSATELLITE POINT:",
           None]]


# For loop to print out analytical mcgubbins
for i in range(len(toLoop)):
    # Define some parameters to feed into the functions
    problemParametersError, problemParametersMappingError, name, errorFunction = toLoop[i]

    # Determine the type of function, and asign inputs and names
    if errorFunction == nadMapError2:
        errorFunctionInputs = [problemParametersError, D0, epsilon, u.m]
        errorFunctionName = "Nadir Angle"
    elif errorFunction == None:
        errorFunctionInputs = []
        errorFunctionName = "Direct"
    elif errorFunction == timingMapError2:
        errorFunctionInputs = [problemParametersError, Vsc, u.m]
        errorFunctionName = "Timing Map Error (Basic)"


    # Run printout function with given parameters
    errorMappingCalc = mapErrorPrintout(problemParametersError,
                     problemParametersMappingError,
                     name=name,
                     errorFunction=errorFunction,
                     errorFunctionInputs=errorFunctionInputs,
                     errorFunctionName=errorFunctionName,
                     decimals=decimals,
                     basePrint=basePrint, printing=printing)

    toLoop[i].append(errorMappingCalc)

# Get total pointing error wrt Nadir

errorInfo = np.array(toLoop)

nadirErrorInfo = []
for i in range(len(errorInfo)):
    problemParametersError, problemParametersMappingError, name, errorFunction, errorMappingCalc = toLoop[i]
    if (errorFunction == nadMapError2) or (name == "ORBIT DETERMINATION:") or (name == "TIMING:"):
        nadirErrorInfo.append(toLoop[i])

nadirErrorRandom = np.array(nadirErrorInfo)

nadirRSSRandomCalc = doRSS(nadirErrorRandom[:, -1])
nadirRSSRandom = doRSS(nadirErrorRandom[:, 1])



if printing:
    print("\nNADIR RSS CALCULATIONS:")
    print("Nadir RSS based on (my) calculated values for random     = %s \n"
          "Nadir RSS based on given values for random               = %s" %(nadirRSSRandomCalc, nadirRSSRandom))

nadirErrorRandomToSave = np.copy(nadirErrorRandom)

for i in range(len(nadirErrorRandomToSave)):
    errorFunction = nadirErrorRandom[i, -2]
    if errorFunction != None:
        errorFunctionName = errorFunction.__name__
    nadirErrorRandomToSave[i, -2] = errorFunctionName

npArray2LatexTable(nadirErrorRandomToSave, "nadirRandomErrors.txt")


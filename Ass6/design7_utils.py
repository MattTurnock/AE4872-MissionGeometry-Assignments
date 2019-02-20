#!/usr/bin/env python

"""
design7_utils.py: A bunch of useful functions for the DESIGN-7 assignment
"""

from json_to_dict import constants
import numpy as np
from astropy import units as u
from Ass1.kep_orbit_utils import get_circular_velocity_units

pi=np.pi

def azMapError(dphi, D, eta, units=None):
    """
    Function to determine azimuthal mapping error.
    :param dphi: Error in azimuth angle phi
    :param D: Nominal distance from spacecraft to satellite
    :param eta: Nominal nadir angle eta
    :param units: The units to output as. Will give regular number if None. Angles must be radians and distances must be
    km in that case, output is in km
    :return error: The mapping error (in km, or whatever is given)
    """

    # Convert to radians if astropy units used
    if units != None:
        dphi = dphi.to(u.rad)

    # Calculate error
    error = dphi*D*np.sin(eta)

    # Convert to given output units
    if units != None:
        error = error.to(units, equivalencies=u.dimensionless_angles())

    return error

def nadMapError(deta, D, epsilon, units=None):
    """
    Function to determine nadir mapping error.
    :param deta: Error in nadir angle eta
    :param D: Nominal distance from spacecraft to satellite
    :param epsilon: Nominal observation angle epsilon
    :param units:The units to output as. Will give regular number if None. Angles must be radians and distances must be
    km in that case, output is in km
    :return error: The mapping error (in km, or whatever is given)
    """
    # Convert to radians if astropy units used
    if units != None:
        deta = deta.to(u.rad)

    # Calculate error
    error = (deta * D) / np.sin(epsilon)

    # Convert to given output units
    if units != None:
        error = error.to(units, equivalencies=u.dimensionless_angles())

    return error

def nadMapError2(paramList):
    """
    Function to allow other functions in part1.py work
    :param paramList: List of parameters, filled with identical ones to nadMapError
    :return error: Error identical to nadMapError
    """
    deta, D, epsilon, units = paramList
    error = nadMapError(deta, D, epsilon, units=units)
    return error

def clockMapError(dT, lat, Ve=464.0, units=None):
    # Calculate error
    error = dT * Ve * np.cos(lat)

    # Convert to given output units
    if units != None:
        error = error.to(units)

    return error

def timingMapError(dT, Vsc, units=None):
    """
    Function to find out (basic) timing error. Purely based on speed of s/c
    :param dT: Timing error (in s if no units given)
    :param Vsc: Speed of spacecraft (in m/s if no units given)
    :param units: The units to output as. Will give regular number if None.
    :param extra: Literally just there to make the function in part1.py work
    :return error:  The mapping error (in m, or whatever is given)
    """
    # Calculate error
    error = dT * Vsc

    # Convert to given output units
    if units != None:
        error = error.to(units)

    return error

def timingMapError2(paramList):
    """
    Function to allow other functions in part1.py work. Parameters identical to timingMapError
    :param paramList: List of parameters, filled with identical ones to timingMapError
    :return error:  Error identical to timingMapError
    """
    dT, Vsc, units = paramList
    error = timingMapError(dT, Vsc, units=units)
    return error

def doRSS(errorArray):
    """
    Function that will take an array (or list) and perform RSS
    :param errorArray: Array of values to perform RSS on
    :return RSS: The Root Sum Square of the listed values
    """
    SS = 0
    for i in range(len(errorArray)):
        error = errorArray[i]
        SS += error**2
    RSS = np.sqrt(SS)

    return RSS

def mapErrorPrintout(problemParametersError, problemParametersMappingError, name="", errorFunction=nadMapError,
                     errorFunctionInputs=[], errorFunctionName="", decimals=1, basePrint="%s %s %s %s %s %s",
                     printing=True):
    """
    Function to printout map errors TODO: finish this description and parameters description
    :param problemParametersError:
    :param problemParametersMappingError:
    :param name:
    :param errorFunction:
    :param errorFunctionInputs:
    :param errorFunctionName:
    :param decimals:
    :param basePrint:
    :return:
    """
    if printing: print("\n%s" %name)
    if errorFunction == None:
        errorMappingCalc = problemParametersError
    else:
        errorMappingCalc = errorFunction(errorFunctionInputs)

    errorMappingCalcRound = np.around(errorMappingCalc, decimals=decimals)
    if errorMappingCalcRound != problemParametersMappingError:
        outcome = "UNACCEPTABRU"
    else:
        outcome = "OK"
    if printing:
        print(basePrint % (problemParametersError,
                           problemParametersMappingError,
                           errorFunctionName,
                           errorMappingCalc,
                           errorMappingCalcRound,
                           outcome))

    return errorMappingCalc


def npArray2LatexTable(array, savename):
    np.savetxt(savename, array, fmt="%s", delimiter="\t&\t", newline="      // \n")


class ProblemParameters:

    # General parameters
    def __init__(self, h=1000*u.km, epsilon=21.6052*u.deg, lat=0*u.deg, Xssp=0*u.m):
        self.h = h
        self.epsilon = epsilon
        self.eta = 90*u.deg - epsilon
        self.D0 = (h/np.sin(epsilon)).to(u.km)
        self.lat = lat
        self.mu = constants["muEarth"]
        self.Re = constants["RE"]
        self.Rs = self.Re + self.h
        self.Vsc = get_circular_velocity_units(self.mu, self.Rs, units=u.m/u.s)
        self.Xtarg_nom = Xssp + self.h * np.tan(self.eta)

    # Individual errors
    starSensorMeasurement    = 0.0015*u.deg
    starSensorMounting       = 0.0020*u.deg
    starCatalog              = 0.0001*u.deg
    attitudeComp             = 0.0001*u.deg
    payloadSensorMeasurement = 0.0010*u.deg
    targetCentroiding        = 0.0020*u.deg
    payloadSensorMounting    = 0.0010*u.deg
    coordTransformation      = 0.0001*u.deg
    orbitDetermination       = 100*u.m
    timing                   = 50*u.ms
    projection               = 700*u.m
    subsatellitePoint        = 450*u.m

    payloadSensorMountingSystematic = 0.005*u.deg

    # Individual effects of mapping (known)
    starSensorMeasurementMapping    = 193.1*u.m
    starSensorMountingMapping       = 257.5*u.m
    starCatalogMapping              = 12.9*u.m
    attitudeCompMapping             = 12.9*u.m
    payloadSensorMeasurementMapping = 128.8*u.m
    targetCentroidingMapping        = 257.5*u.m
    payloadSensorMountingMapping    = 128.8*u.m
    coordTransformationMapping      = 12.9*u.m
    orbitDeterminationMapping       = 100.0*u.m
    timingMapping                   = 367.5*u.m
    projectionMapping               = 700.0*u.m
    subsatellitePointMapping        = 450.0*u.m



def doMonteCarlo(iterations, randIndices, randSigmas, Xtarg_nom=0, Vsc=0, eta=0, h=0, systIndices=np.array([]), systValues=np.array([])):
    # TODO: add neat documentation
    arrayLength = 27


    monteArray = np.zeros((iterations, arrayLength), dtype=object)

    # Plonk random numbers into the random slots and the systematic ones

    for i in range(len(monteArray)):
        row = monteArray[i]
        for j in range(len(randIndices)):
            randIndex = randIndices[j]
            row[randIndex] = np.random.randn() * randSigmas[j]
        for j in range(len(systIndices)):
            systIndex = systIndices[j]
            row[systIndex] = systValues[j]

    # Sum systematic and random errors for dX
    monteArray[:, 2] = np.sum(monteArray[:, 0:2], axis=1)

    # Sum systematic and random errors for dT
    monteArray[:, 5] = np.sum(monteArray[:, 3:5], axis=1)

    # Sum all errors in dTheta
    monteArray[:, 22] = np.sum(monteArray[:, 6:22], axis=1)

    # Calculate dX_SSP using dX_ode and dT
    dTs = np.copy(monteArray[:, 5])
    dXs = np.zeros(len(dTs), dtype=object)
    for i in range(len(dTs)):
        dT = dTs[i]
        dXs[i] = (dT * Vsc).to(u.m)

    monteArray[:, 23] = monteArray[:, 2] + dXs

    # Calculate new theta from theta_nom and dtheta
    toAdd = np.ones(len(monteArray[:, 24]), dtype=object)
    for i in range(len(toAdd)):
        toAdd[i] = eta

    monteArray[:, 24] = monteArray[:, 22] + toAdd

    # Calculate x_targ from equation, and also dXtargs
    Xtargs = np.ones(len(monteArray[:, 25]), dtype=object)
    dXtargs = np.copy(Xtargs)
    for i in range(len(Xtargs)):
        Xssp = monteArray[i, 23]
        theta = monteArray[i, 24]
        Xtarg = Xssp + h * np.tan(theta)

        Xtargs[i] = Xtarg
        dXtargs[i] = Xtarg - Xtarg_nom

    monteArray[:, 25] = Xtargs
    monteArray[:, 26] = dXtargs


    dXtargs = monteArray[:, 26]


    return monteArray
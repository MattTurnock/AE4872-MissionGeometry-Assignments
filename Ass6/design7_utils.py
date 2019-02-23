#!/usr/bin/env python

"""
design7_utils.py: A bunch of useful functions for the DESIGN-7 assignment
"""
from matplotlib import pyplot as plt
import numpy as np
from astropy import units as u
from Ass1.kep_orbit_utils import get_circular_velocity_units
from json_to_dict import constants
pi=np.pi

class ProblemParameters:
    """
    Class that contains all parameters for the problem (errors, altitudes etc)
    """
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
    :param units: The units to output as. Will give regular number if None. Angles must be radians and distances must be
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
    """
    Function to plotting clock mapping error, based on equation in OCDM
    :param dT: Timing error
    :param lat: Latitude
    :param Ve: Earth rotation velocity at equator
    :param units: The units to output as. Will give regular number if None.
    :return error: Mapping error
    """
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

    Function to return an array of mapping errors for a given parameter, and print if necessary
    :param problemParametersError: The representative error given in ProblemParameters class
    :param problemParametersMappingError: The mapping error given in ProblemParameters class that has been calculated
    by OCDM
    :param name: Name given to the parameter
    :param errorFunction: The error function that is to be uesed to analytically calculate error
    :param errorFunctionInputs: List of inputs to pass to the errorFunction
    :param errorFunctionName: String name for the error function, more descriptive
    :param decimals: Number of decimals to round to
    :param basePrint: Just a string to print into
    :param printing: Boolean of if printing should happen
    :return errorMappingCalc: The value of the mapping error based on representative value
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
    """
    Simple function to convert numpy array to a latex table component
    :param array: Array to convert
    :param savename: Name to give to generated txt file
    :return None:
    """
    np.savetxt(savename, array, fmt="%s", delimiter="\t&\t", newline="      \\\ \n\hline \n")

def doMonteCarlo(iterations, randIndices, randSigmas, Xtarg_nom=0, Vsc=0, eta=0, h=0, systIndices=np.array([]),
                 systValues=np.array([])):
    """
    Function to fully do the Monte-Carlo analysis
    :param iterations: Number of Monte-Carlo iterations to perform
    :param randIndices: np array containing indices in monteArray of random errors
    :param randSigmas: Values to place at locations defined by randIndices
    :param Xtarg_nom: Nominal value of the target X coordinate
    :param Vsc: Spacecraft velocity
    :param eta: observation angle
    :param h: Spacecraftv altitude
    :param systIndices: np array containing indices in monteArray of systematic errors
    :param systValues: Values to place at locations defined by systIndices
    :return monteArray: An array containing all Monte-Carlo pertinent info
    """
    # Array length specific to this problem
    arrayLength = 27

    # Create base array to be filled
    monteArray = np.zeros((iterations, arrayLength), dtype=object)

    # Plonk random numbers into the proper random slots and the systematic ones
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

    return monteArray

def doPlot(monteCarloArrayFinalCol, case="nocase", show=True, save=True, binno=100, linewidth=0.5, fontsz=20):
    """
    Function to plot a histogram and fit a normal distribution to it
    :param monteCarloArrayFinalCol: Numpy array listing all Monte-Carlo mapping errors
    :param case: A string to define the case for saving
    :param show: Choose to display plot
    :param save: Choose to save plot
    :param binno: Number of bins for the historgram
    :param linewidth: Width of histogram outline lines
    :return None:
    """

    # Define the plotting values (since given array has quantities)
    monteCarloArrayFinalColValues = []
    for i in range(len(monteCarloArrayFinalCol)):
        monteCarloArrayFinalColValues.append(monteCarloArrayFinalCol[i].value)

    # Define the gaussian distribution to overlay
    var = np.var(monteCarloArrayFinalColValues)
    avg = np.mean(monteCarloArrayFinalColValues)
    pdf_x = np.linspace(min(monteCarloArrayFinalColValues), max(monteCarloArrayFinalColValues), 100)
    pdf_y = 1.0/np.sqrt(2*np.pi*var)*np.exp(-0.5*(pdf_x-avg)**2/var)

    # Plot the things
    plt.figure()
    plt.hist(monteCarloArrayFinalColValues, bins=binno, density=True, histtype="bar", edgecolor="black", linewidth=linewidth)
    plt.plot(pdf_x, pdf_y)
    plt.xlabel("Mapping Error Magnitude [m]")
    plt.ylabel("Probability Density")
    ax = plt.gca()
    for item in ([ax.title, ax.xaxis.label, ax.yaxis.label, ax.yaxis.get_offset_text()] +
                 ax.get_xticklabels() + ax.get_yticklabels()):
        item.set_fontsize(fontsz)
    if save: plt.savefig("histogram%s.pdf" %case, bbox_inches="tight")
    if show: plt.show()

def swapCols(array, swapIndices):
    """
    Simple function to swap 2 columns in a numpy array
    :param array: Array to swap columns in
    :param swapIndices: length-2 list of indices to swap
    :return array: The new array with swapped columns
    """
    array[:, swapIndices[0]], array[:, swapIndices[1]] = array[:, swapIndices[1]], array[:, swapIndices[0]].copy()
    return array

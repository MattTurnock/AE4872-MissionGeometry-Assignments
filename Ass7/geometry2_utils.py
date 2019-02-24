#!/usr/bin/env python

"""
geometry2_utils.py: A bunch of useful functions for the GEOMETRY-2 assignment
"""
__author__      = "Matthew Turnock"
__email__ = "matthew.turnock@protonmail.com"
__version__ = "1.0"

########################################################################################################################

import numpy as np
from astropy import units as u
from numpy import pi, sin, cos, tan

def getAngularDistance(a1, a2, d1, d2, units=u.rad):
    """
    Function to get angular distance between 2 points P1 and P2
    :param a1: alpha of P1
    :param d1: delta of P1
    :param a2: alpha of P2
    :param d2: delta of P2
    :param units: astropy units to output. If None will input and output whatever is put in (deg or rad)
    :return theta12: The angular distance between P1 and P2
    """
    if units is None:
        a1 = a1 * u.rad
        a2 = a2 * u.rad
        d1 = d1 * u.rad
        d2 = d2 * u.rad

    theta12 = np.arccos( sin(d1)*sin(d2) + cos(d1)*cos(d2)*cos(a1-a2))

    if units is None:
        theta12 = theta12.value

    else:
        theta12 = theta12.to(units)

    return theta12

def getTriangleAngularDistances(inputList, units=u.rad):
    """
    Function to take 3 points of a spherical triangle, P1 P2 and P3, with coordinates ai and di, and find all 3 angular
    distances theta
    :param inputList: List of coords of form [a1, a2, a3, d1, d2, d3]
    :param units: astropy units to output. If None will input and output whatever is put in (deg or rad)
    :return theta12, theta13, theta23: The 3 angular distances between spherical triangle points
    """

    a1, a2, a3, d1, d2, d3 = inputList

    theta12 = getAngularDistance(a1, a2, d1, d2, units=units)
    theta13 = getAngularDistance(a1, a3, d1, d3, units=units)
    theta23 = getAngularDistance(a2, a3, d2, d3, units=units)

    if units is not None:
        theta12 = theta12.to(units)
        theta13 = theta13.to(units)
        theta23 = theta23.to(units)

    return theta12, theta13, theta23

def getRotationAngle(inputList, getFrom="coords", type="inside", units=u.rad):
    """
    Function to find rotation angle about point P3, from P1 to P2
    :param inputList: List of possible inputs.
    If getFrom is "coords" list is as [a1, a2, a3, d1, d2, d3] where ai and
    di are the alpha and delta of point Pi.
    If getFrom is "theta" list is as [theta12, theta13, theta23] where thetaij is the angular distance between points
    Pi and Pj
    :param getFrom: Boolean to determine what inputList format to use
    :param units: astropy units to output. If None will input and output radians (MUST BE RADIANS if NONE USED)
    :return Phi12_3: The rotation angle from P1 to P2 about P3
    """
    if getFrom == "theta":
        theta12, theta13, theta23 = inputList
    elif getFrom == "coords":
        theta12, theta13, theta23 = getTriangleAngularDistances(inputList, units=units)

    Phi12_3 = np.arccos( (cos(theta12) - cos(theta13)*cos(theta23)) / (sin(theta13)*sin(theta23)))

    if type == "outside":
        if units is not None:
            Phi12_3.to(u.rad)
        Phi12_3 = 2*pi*u.rad - Phi12_3

    if units is not None:
        Phi12_3 = Phi12_3.to(units)

    return Phi12_3

def getTotalArea(inputList, getFrom="coords", type="inside", R=1*u.m, units=u.m**2, angleUnits=u.rad):
    """
    Function to find the total area of a spherical triangle
    :param inputList: List of possible inputs.
    If getFrom is "coords" list is as [a1, a2, a3, d1, d2, d3] where ai and
    di are the alpha and delta of point Pi.
    If getFrom is "Phi" list is as [Phi23_1, Phi31_2, Phi12_3] where Phiij_k is the rotation angle
    :param getFrom: Boolean to define what type of inputList is used
    :param type: Boolean to determine if inside or outside triangle area is to be found
    :param R: Radius of the sphere. Default is unit sphere R=1*u.m. Must be input without astropy units if None used
    for other units
    :param units: Output units. Set to None if input units are also None
    :param angleUnits: Units to use for angles. Should be u.rad or None
    :return Omega_t: The spherical triangle area
    """
    # if (units is None) or (angleUnits is None):
    #     R = R.value
    if getFrom == "Phi":
        Phi23_1, Phi31_2, Phi12_3 = inputList
    elif getFrom == "coords":
        a1, a2, a3, d1, d2, d3 = inputList
        Phi23_1 = getRotationAngle([a2, a3, a1, d2, d3, d1], getFrom="coords", type=type, units=angleUnits)
        Phi31_2 = getRotationAngle([a3, a1, a2, d3, d1, d2], getFrom="coords", type=type, units=angleUnits)
        Phi12_3 = getRotationAngle([a1, a2, a3, d1, d2, d3], getFrom="coords", type=type, units=angleUnits)
        # Phi23_1, Phi31_2, Phi12_3 = getRotationAngle(inputList, getFrom="coords", type=type, units=angleUnits)
        # print(np.rad2deg(np.array([Phi23_1, Phi31_2, Phi12_3])))

    if angleUnits is not None:
        piUse = pi * u.rad
        Phi23_1 = Phi23_1.to(u.rad)
        Phi31_2 = Phi31_2.to(u.rad)
        Phi12_3 = Phi12_3.to(u.rad)
    else:
        piUse = pi

    Omega_ti = R**2 * (Phi23_1 + Phi31_2 + Phi12_3 - piUse)
    Omega_t = Omega_ti
    # if type == "inside":
    #     Omega_t = Omega_ti
    # elif type == "outside":
    #     Omega_t = 4 * piUse * R**2 - Omega_ti

    if units is not None:
        Omega_t = Omega_t.to(units, equivalencies=u.dimensionless_angles())

    return Omega_t

def findAllValues(inputArray):
    """
    Function to find all the values requested for the assignment
    :param inputArray: 2 dimenstional array with values to convert. Of form [[a1, a2, a3, d1, d2, d3, areaType],...].
    as and ds refer to alphas and deltas of P1 P2 and P3 making up triangle. areaType refers to calculation of inside or
    outside triangle area
    :return outputArray: 2 dimensional array containing output values.
    Of form [theta12, theta13, theta23, Phi23_1, Phi31_2, Phi12_3, Omega_t]
    """
    outputArray = np.zeros((len(inputArray), 7), dtype=object)
    for i in range(len(inputArray)):
        a1, a2, a3, d1, d2, d3, areaType = inputArray[i, :]

        # Find angular distances
        theta12 = getAngularDistance(a1, a2, d1, d2, units=u.deg)
        theta13 = getAngularDistance(a1, a3, d1, d3, units=u.deg)
        theta23 = getAngularDistance(a2, a3, d2, d3, units=u.deg)

        # Find rotation angles
        Phi23_1 = getRotationAngle([a2, a3, a1, d2, d3, d1], getFrom="coords", type=areaType, units=u.deg)
        Phi31_2 = getRotationAngle([a3, a1, a2, d3, d1, d2], getFrom="coords", type=areaType, units=u.deg)
        Phi12_3 = getRotationAngle([a1, a2, a3, d1, d2, d3], getFrom="coords", type=areaType, units=u.deg)

        # Find area
        Omega_t = getTotalArea([a1, a2, a3, d1, d2, d3], getFrom="coords", type=areaType, R=1*u.m, units=u.m**2, angleUnits=u.rad)

        outputArray[i, :] = [theta12, theta13, theta23, Phi23_1, Phi31_2, Phi12_3, Omega_t]

    return outputArray

#!/usr/bin/env python

"""
fullsky11_utils.py: A bunch of useful functions for the FULLSKY-11 assignment
"""
__author__      = "Matthew Turnock"
__email__ = "matthew.turnock@protonmail.com"
__version__ = "1.0"

########################################################################################################################
import numpy as np
from astropy import units as u

sin = np.sin
cos = np.cos
pi = np.pi

######################################################################################################################
# Generic functions

def mod360(angle, outputUnit=None):
    """
    Function to perform mod360 on numbers. if no unit given assumes radians
    :param angle: Angle to perform mod on
    :param outputUnit: Units to output in.
    :return angleMod: modulated angle
    """
    if (outputUnit is not None):
        angleVal = (angle.to(u.rad, equivalencies=u.dimensionless_angles())).value
    else:
        angleVal = angle

    angleMod = angleVal % (2*pi)

    if outputUnit is not None:
        angleMod = angleMod * u.rad
        angleMod = angleMod.to(outputUnit, equivalencies=u.dimensionless_angles())

    return angleMod

def H(phi, useAstroUnit=False):
    """
    Function to apply modulo360 to an angle phi, then set H as 1 or -1. If no units are given, assume radians are used
    :param phi: Angle to apply H to
    :param useAstroUnit: Boolean to choose if astrounits are used
    :return H: 1 or -1 depending on outcome
    """
    if (useAstroUnit is True):
        phiVal = (phi.to(u.rad, equivalencies=u.dimensionless_angles())).value
    else:
        phiVal = phi

    phiMod = phiVal % (2*pi)
    if 0 <= phiMod <= pi:
        H = 1
    elif pi <= phiMod <= 2*pi:
        H = -1
    else:
        return "Unkown stoopid error"

    return H

def acos2(alpha, Hfun, outputUnits=None):
    """
    version of acos function with H function applied to see if the output put should be negative. also applies mod360
    :param alpha: (dimensionless) number to perform acos on
    :param Hfun: Relevant hemisphere function to determine sign
    :param outputUnits: Astropy units to output as. None returns radians
    :return angelOut: Returns the acos angle
    """
    angleOut = np.arccos(alpha) * Hfun
    angleOut = mod360(angleOut, outputUnit=outputUnits)
    return angleOut

def getn(T, outputUnits=None):
    """
    Get the mean motion of a satellite with given period
    :param T: Period of satellite orbit
    :param outputUnits: Units to output n as. If None is rad/s
    :return n: mean motion of satellite
    """
    n = 2*pi/T
    if outputUnits is not None:
        n = n.to(outputUnits, equivalencies=u.dimensionless_angles())

    return n

######################################################################################################################
# Intermediate functions

def getphi1_2(phi10, omega1, t, outputUnits=u.deg):
    """
    Function to get phi1 or phi2 (the function is the same so we can combine them)
    :param phi10: phi10 or phi20. The initiaL FINISH
    :param omega1: omega1 or omega2. The rotation rate of S around C or P around S respectively
    :param t: elapsed time
    :param outputUnits: units to output phi12 as
    :return phi12: the wanted phi1 or phi2
    """
    phi12 = phi10 + omega1*t
    if outputUnits is not None:
        phi12 = phi12.to(outputUnits)

    return phi12

def getDeltaalpha(rho1, rho2, delta, phi2, outputUnits=None):
    """
    Function to get change in alpha
    :param rho1: angular distance between C and S
    :param rho2: angular distance between P and S
    :param delta: elevation angle
    :param phi2: azimuth angle about S
    :param outputUnits: Units to output as. If None is in radians
    :return returnAngle: returns DeltaAlpha, the change in azimuth angle
    """
    if outputUnits is None:
        useAstroUnit = False
    else:
        useAstroUnit = True

    top = cos(rho2) - cos(rho1)*sin(delta)
    bottom = sin(rho1)*cos(delta)
    returnAngle = acos2(top/bottom, -H(phi2, useAstroUnit=useAstroUnit), outputUnits=outputUnits)
    returnAngle = mod360(returnAngle, outputUnit=outputUnits)

    if useAstroUnit is True:
        returnAngle = returnAngle.to(outputUnits, equivalencies=u.dimensionless_angles())

    return returnAngle

def getrhoE(deltaEPrime, delta, Deltaalpha, outputUnits=None):
    """
    Gets the angle to euler axis E
    :param deltaEPrime: co-elevation of E
    :param delta: elevation angle
    :param Deltaalpha: change in azimuth angle
    :param outputUnits: Units to output as. If None returns radians
    :return rhoE: angle to E
    """
    rhoE = np.arccos(cos(deltaEPrime)*sin(delta) + sin(deltaEPrime)*cos(delta)*cos(Deltaalpha))

    if outputUnits is not None:
        rhoE = rhoE.to(outputUnits, equivalencies=u.dimensionless_angles())

    return rhoE

def getomegaE(omega2, rho1, deltaEPrime, outputUnits=None):
    """
    Intermediate function. Is a rotation rate so use outputunits=u.Unit("rad/s") for example. Otherwise assumes radians/s
    :param omega2: roatation rate of P around S
    :param rho1: angle between C and S
    :param deltaEPrime: co-elevation of E
    :param outputUnits: Units to output as. If None returns radians/s
    :return omegaE: rotation rate about E
    """
    top = omega2*sin(rho1)
    bottom = sin(deltaEPrime)
    omegaE = top/bottom

    if outputUnits is not None:
        omegaE = omegaE.to(outputUnits, equivalencies=u.dimensionless_angles())

    return omegaE

def getomegaE2(omega1, omega2, rho1, outputUnits=None):
    """
    Intermediate function. Is a rotation rate so use outputunits=u.Unit("rad/s") for example. Otherwise assumes radians/s
    :param omega1:
    :param omega2:
    :param rho1:
    :param outputUnits:
    :return:
    """
    omegaE = np.sqrt(omega1**2 + omega2**2 + 2*omega1*omega2*cos(rho1))

    if outputUnits is not None:
        omegaE = omegaE.to(outputUnits, equivalencies=u.dimensionless_angles())

    return omegaE

def getdeltaEPrime(omega1, omega2, rho1, outputUnits=None):
    """
    Intermediate function. Assumes rad/s and rad if no units given
    :param omega1:
    :param omega2:
    :param rho1:
    :param outputUnits:
    :return:
    """
    top = omega2 * sin(rho1)
    bottom = omega1 + omega2*cos(rho1)
    if outputUnits is not None:
        top = top.to(u.rad/u.s, equivalencies=u.dimensionless_angles())
        bottom = bottom.to(u.rad/u.s, equivalencies=u.dimensionless_angles())
    deltaEPrime = np.arctan(top/bottom)

    if outputUnits is not None:
        deltaEPrime = deltaEPrime.to(outputUnits, equivalencies=u.dimensionless_angles())

    return deltaEPrime

def getDeltaPsi(deltaEPrime, rhoE, delta, Deltaalpha, outputUnits=None):
    """
    intermediate function. Assumes radians if no units given
    :param deltaEPrime:
    :param rhoE:
    :param delta:
    :param Deltaalpha:
    :param outputUnits:
    :return:
    """
    top = cos(deltaEPrime) - cos(rhoE)*sin(delta)
    bottom = sin(rhoE)*cos(delta)

    if outputUnits is None:
        useAstroUnit = False
    else:
        useAstroUnit = True

    HDalpha = H(Deltaalpha, useAstroUnit=useAstroUnit)
    outAngle = acos2((top/bottom), HDalpha, outputUnits=outputUnits)

    if useAstroUnit:
        outAngle = outAngle.to(outputUnits, equivalencies=u.dimensionless_angles())

    outAngle = mod360(outAngle, outputUnit=outputUnits)

    return outAngle


######################################################################################################################
# Results functions

def getalpha(phi1, Deltaalpha, outputUnit=None):
    """
    Final coordinat alpha. If no units given, assume radians
    :param phi1:
    :param Deltaalpha:
    :param outputUnit:
    :return:
    """
    alpha = mod360(phi1 + Deltaalpha, outputUnit=outputUnit)

    return alpha

def getdelta(rho1, rho2, phi2, outputUnit=None):
    """
    Function to get delta coord. Uses radians if outputUnit=None
    :param rho1:
    :param rho2:
    :param phi2:
    :param outputUnit:
    :return:
    """
    minAngle = np.arccos(cos(rho1)*cos(rho2) + sin(rho1)*sin(rho2)*cos(phi2))

    if outputUnit is None:
        delta = pi/2 - minAngle
    else:
        minAngle = minAngle.to(u.rad, equivalencies=u.dimensionless_angles())
        delta = (pi/2)*u.rad - minAngle
        delta = delta.to(outputUnit, equivalencies=u.dimensionless_angles())

    return delta

def getPsi(DeltaPsi, outputUnit=None):
    """
    Function to get Psi. Uses radians if outputUnit=None
    :param DeltaPsi:
    :param outputUnit:
    :return:
    """
    if outputUnit is None:
        PsiTemp = DeltaPsi - pi/2
    else:
        DeltaPsi = DeltaPsi.to(u.rad, equivalencies=u.dimensionless_angles())
        PsiTemp = DeltaPsi - (pi/2)*u.rad

    Psi = mod360(PsiTemp, outputUnit=outputUnit)

    if outputUnit is not None:
        Psi = Psi.to(outputUnit, equivalencies=u.dimensionless_angles())

    return Psi

def getV(omegaE, rhoE, outputUnits=None):
    """
    Function to get final velocity. Assumes rad and s if no units given
    :param omegaE:
    :param rhoE:
    :param outputUnits:
    :return:
    """
    V = omegaE * sin(rhoE)

    if outputUnits is not None:
        V = V.to(outputUnits, equivalencies=u.dimensionless_angles())

    return V

#######################################################################################################################
# Full program to combine all, for a given time and initial conditions

def doPartialDualAxisSpiral(rho1, rho2, omega1, omega2, phi1, phi2,
                            angleUnits=None, omegaUnits=None):
    """
    Function to do half of full needs. Takes major parameters and give an output. Assumes standard SI units when not specified
    :param rho1:
    :param rho2:
    :param omega1:
    :param omega2:
    :param phi1:
    :param phi2:
    :param angleUnits:
    :param omegaUnits:
    :param velocityUnits:
    :return:
    """
    delta = getdelta(rho1, rho2, phi2, outputUnit=angleUnits)

    Deltaalpha = getDeltaalpha( rho1, rho2, delta, phi2, outputUnits=angleUnits )

    deltaEPrime = getdeltaEPrime(omega1, omega2, rho1, outputUnits=angleUnits)

    rhoE = getrhoE(deltaEPrime, delta, Deltaalpha, outputUnits=angleUnits)

    omegaE = getomegaE2(omega1, omega2, rho1, outputUnits=omegaUnits)

    DeltaPsi = getDeltaPsi(deltaEPrime, rhoE, delta, Deltaalpha, outputUnits=angleUnits)

    alpha = getalpha(phi1, Deltaalpha, outputUnit=angleUnits)

    Psi = getPsi(DeltaPsi, outputUnit=angleUnits)

    V = getV(omegaE, rhoE, outputUnits=omegaUnits)

    return delta, Deltaalpha, alpha, deltaEPrime, rhoE, omegaE, V, DeltaPsi, Psi

def doDualAxisSpiral(rho1, rho2, omega1, omega2, phi10, phi20, t, angleUnits=None, omegaUnits=None, returnValue=False):
    """
    Is like the PArtial dual axis funtion, but does it from initial conditions and t instead of known phi1 and phi2
    :param rho1:
    :param rho2:
    :param omega1:
    :param omega2:
    :param phi10:
    :param phi20:
    :param t:
    :param angleUnits:
    :param omegaUnits:
    :param velocityUnits:
    :return:
    """
    phi1 = getphi1_2(phi10, omega1, t, outputUnits=angleUnits)
    phi2 = getphi1_2(phi20, omega2, t, outputUnits=angleUnits)

    outputList = list(doPartialDualAxisSpiral(rho1, rho2, omega1, omega2, phi1, phi2,
                                     angleUnits=angleUnits, omegaUnits=omegaUnits))
    outputList.append(phi1)
    outputList.append(phi2)
    if returnValue is True:
        for i in range(len(outputList)):
            outputList[i] = outputList[i].value
    return tuple(outputList)

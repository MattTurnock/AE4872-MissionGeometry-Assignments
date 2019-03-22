#!/usr/bin/env python

"""
expo4_utils.py: A bunch of useful functions for the EXPO-4 assignment
"""
__author__      = "Matthew Turnock"
__email__ = "matthew.turnock@protonmail.com"
__version__ = "1.0"

########################################################################################################################
import numpy as np
from astropy import units as u
from json_to_dict import constants

sin = np.sin
cos = np.cos
pi = np.pi

######################################################################################################################

muSun = constants["muSun"]

r1 = 1 * u.AU
r2 = 1.5 * u.AU
k2 = 1/12
Psi = 90*u.deg
N = 2
thetaExample = 0.0*u.rad
gamma1Example = -80*u.deg

def getgamma1Min(r1, r2, k2, thetaBar, Delta, outputUnits=None):
    """
    Get minimum allowed value for flight path angle gamma_1 (uses constraints given in slides). If no units are None,
    use radians.
    :param r1: Departure radius
    :param r2: Arrival radius
    :param k2: k2 as given in slides
    :param thetaBar: thetaBar as given in slides (function defined below)
    :param Delta: Delta as given in slides (function defined below)
    :param outputUnits: Units to output as in astropy
    :return gamma1Min:
    """
    inner1 = -np.log(r1/r2) * (np.tan(0.5*k2*thetaBar))**-1

    inTan = 0.5*k2 * (inner1 - np.sqrt(Delta))

    gamma1Min = np.arctan(inTan)

    if outputUnits is not None:
        gamma1Min = gamma1Min.to(outputUnits)

    return gamma1Min


def getgamma1Max(r1, r2, k2, thetaBar, Delta, outputUnits=None):
    """
    Get maximum allowed value for flight path angle gamma_1 (uses constraints given in slides). If no units are None,
    use radians.
    :param r1: Departure radius
    :param r2: Arrival radius
    :param k2: k2 as given in slides
    :param thetaBar: thetaBar as given in slides (function defined below)
    :param Delta: Delta as given in slides (function defined below)
    :param outputUnits: Units to output as in astropy
    :return gamma1Max:
    """
    inner1 = -np.log(r1/r2) * (np.tan(0.5*k2*thetaBar))**-1

    inTan = 0.5*k2 * (inner1 + np.sqrt(Delta))

    gamma1Max = np.arctan(inTan)

    if outputUnits is not None:
        gamma1Max = gamma1Max.to(outputUnits)

    return gamma1Max



def getDelta(r1, r2, k2, thetaBar):
    """
    Get intermediate parameter Delta for previous functions. Dimensionless therefore can use astropy units or not.
    :param r1: Departure radius
    :param r2: Arrival radius
    :param k2: k2 as given in slides
    :param thetaBar: thetaBar as given in slides (function defined below)
    :param Delta: Delta as given in slides (function defined below)
    :param outputUnits: Units to output as in astropy
    :return Delta:
    """

    term1 = ( 2*(1 - cos(k2*thetaBar)))/k2**4
    term2 = (np.log(r1/r2))**2

    Delta = term1 - term2

    return Delta



def getthetaBar(Psi, N, outputUnits=None):
    """
    Get value for thetaBar as given in slides. If not units given use radians
    :param Psi: Psi from slides
    :param N: N from slides
    :param outputUnits: Units to output as
    :return thetaBar:
    """
    twoPiN = 2 * np.pi * N

    if outputUnits is not None:
         twoPiN = twoPiN * u.rad

    thetaBar = Psi + twoPiN

    if outputUnits is not None:
        thetaBar = thetaBar.to(outputUnits)

    return thetaBar



thetaBar = getthetaBar(Psi, N, outputUnits=u.rad)
print("thetaBar : %s" %thetaBar)

Delta = getDelta(r1, r2, k2, thetaBar)
print("Delta : %s" %Delta)

gamma1Min = getgamma1Min(r1, r2, k2, thetaBar, Delta, outputUnits=u.deg)
print("gamma1Min : %s" %gamma1Min)

gamma1Max = getgamma1Max(r1, r2, k2, thetaBar, Delta, outputUnits=u.deg)
print("gamma1Max : %s" %gamma1Max)





def getk1(r1, r2, k2, thetaBar, gamma1):
    """
    Gets dynamic range parameter k1 magnitude as given in slides. Dimensionless therefore use astropy units or dont
    :param r1: Departure radius
    :param r2: Arrival radius
    :param k2: k2 as given in slides
    :param thetaBar: thetaBar as given in slides
    :param gamma1: Chosen value for gamma1
    :return k1:
    """
    print(r1, r2, k2)
    top = np.log(r1/r2) + np.tan(gamma1)/k2 * sin(k2*thetaBar)
    bottom = 1 - cos(k2*thetaBar)


    k1Mag = np.sqrt( (top/bottom)**2 + (np.tan(gamma1))**2 / k2**2 )

    k1Sign = np.sign(top)
    k1 = k1Mag * k1Sign

    return k1

k1 = getk1(r1, r2, k2, thetaBar, gamma1Example)
print("k1 : %s" %k1)



def getphi(k1, k2, gamma1, outputUnits=None):
    """
    Gets phase angle phi as shown in slides. If no unit given use radians
    :param k1:
    :param k2:
    :param gamma1:
    :param outputUnits:
    :return phi:
    """
    phi = np.arccos( np.tan(gamma1) / (k1*k2) )

    if outputUnits is not None:
        phi = phi.to(outputUnits)

    return phi



phi = getphi(k1, k2, gamma1Example, outputUnits=u.rad)
print("phi : %s" %phi)


def getk0(r1, k1, phi, outputUnits=None):
    """
    Gets scaling factor as shown in slides. If no unit given use SI
    :param r1:
    :param k1:
    :param phi:
    :param outputUnits:
    :return k0:
    """
    k0 = r1/np.exp(k1*sin(phi))

    if outputUnits is not None:
        k0 = k0.to(outputUnits)

    return k0

k0 = getk0(r1, k1, phi, outputUnits=u.km)
print("k0 : %s" %k0)


def getS(k2, theta, phi):
    """
    Gets intermediate parameter s as shown in slides
    :param k2:
    :param theta:
    :param phi:
    :return:
    """
    s = sin(k2*theta + phi)

    return s

s = getS(k2, thetaExample, phi)
print("s : %s" %s)


def getC(k2, theta, phi):
    """
    Gets intermediate parameter C as shown in slides
    :param k2:
    :param theta:
    :param phi:
    :return:
    """
    c = cos(k2*theta + phi)

    return c

c = getC(k2, thetaExample, phi)
print("c : %s" %c)


def gettangamma(k1, k2, theta, phi):

    tangamma = k1 * k2 * cos(k2*theta + phi)

    return tangamma

tangamma = gettangamma(k1, k2, thetaExample, phi)
print("tan gamma : %s" %tangamma)

gamma = np.arctan(tangamma)
print("gamma : %s" %gamma)



def getaccelLocal(r, mu, outputUnits=None):

    aLocal = mu/r**2

    if outputUnits is not None:
        aLocal = aLocal.to(outputUnits)

    return aLocal

aLocal = getaccelLocal(r1, muSun, outputUnits=u.m/u.s**2)
print("Local accel : %s" %aLocal)

def getaccel(k1, k2, gamma, s, r, mu, outputUnits=None):

    out = np.tan(gamma) / 2*cos(gamma)
    in1 = ( (np.tan(gamma))**2 + k1*k2**2*s + 1 )**-1

    in2_top = k2**2 * (1 - 2*k1*s)
    in2_bottom = ( np.tan(gamma)**2 + k1*k2**2*s + 1 )**2

    #accelN is normalised acceleration wrt solar accel
    accelNTemp = out * (in1 - in2_top/in2_bottom)

    if np.sign(accelNTemp) == 1:
        n = 0
    else:
        n=1

    accelN = (-1)**n * accelNTemp

    accel = accelN * getaccelLocal(r, mu)

    if outputUnits is not None:
        accel = accel.to(outputUnits)

    return accel, n

accel, n = getaccel(k1, k2, gamma, s, r1, muSun, outputUnits=u.km/u.s**2)
print("Acceleration : %s In %s direction" %(accel, n))


def getdtdtheta(k1, k2, gamma, r, s, mu, outputUnits=None):

    top = r**3 * ( np.tan(gamma)**2 + k1*k2**2*s + 1)

    dtdtheta = np.sqrt(top/mu)

    if outputUnits is not None:
        dtdtheta = dtdtheta.to(outputUnits, equivalencies=u.dimensionless_angles())

    return dtdtheta

dtdtheta2 = getdtdtheta(k1, k2, gamma, r1, s, muSun, outputUnits=u.s/u.rad)**2
print("dt/dtheta^2 : %s" %dtdtheta2)


def getr(k0, k1, k2, theta, phi, outputUnits=None):

    r = k0 * np.exp(k1 * np.sin(k2*theta + phi))

    if outputUnits is not None:
        r = r.to(outputUnits)

    return r

r = getr(k0, k1, k2, thetaExample, phi, outputUnits=u.km)
print("r : %s" %r)







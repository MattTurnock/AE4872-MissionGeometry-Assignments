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

# muSun = constants["muSun"]
#
# r1 = (1 * u.AU)
# r2 = 1.5 * u.AU
# k2 = 1/12
# Psi = 90*u.deg
# N = 2
# thetaExample = 0.0*u.rad
# gamma1Example = -80*u.deg






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









# thetaBar = getthetaBar(Psi, N, outputUnits=u.rad)
# print("thetaBar : %s" %thetaBar)
#
# Delta = getDelta(r1, r2, k2, thetaBar)
# print("Delta : %s" %Delta)
#
# gamma1Min = getgamma1Min(r1, r2, k2, thetaBar, Delta, outputUnits=u.deg)
# print("gamma1Min : %s" %gamma1Min)
#
# gamma1Max = getgamma1Max(r1, r2, k2, thetaBar, Delta, outputUnits=u.deg)
# print("gamma1Max : %s" %gamma1Max)











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
    # print(r1, r2, k2)
    top = np.log(r1/r2) + np.tan(gamma1)/k2 * sin(k2*thetaBar)
    bottom = 1 - cos(k2*thetaBar)


    k1Mag = np.sqrt( (top/bottom)**2 + (np.tan(gamma1))**2 / k2**2 )

    k1Sign = np.sign(top)
    k1 = k1Mag * k1Sign

    return k1








#
# k1 = getk1(r1, r2, k2, thetaBar, gamma1Example)
# print("k1 : %s" %k1)
#
#








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










# phi = getphi(k1, k2, gamma1Example, outputUnits=u.rad)
# print("phi : %s" %phi)
#









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









#
# k0 = getk0(r1, k1, phi, outputUnits=u.km)
# print("k0 : %s" %k0)










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








#
# s = getS(k2, thetaExample, phi)
# print("s : %s" %s)










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







#
# c = getC(k2, thetaExample, phi)
# print("c : %s" %c)










def gettangamma(k1, k2, theta, phi):
    """
    Finds tan(gamma) using equation from slides
    :param k1:
    :param k2:
    :param theta:
    :param phi:
    :return tangamma:
    """

    tangamma = k1 * k2 * cos(k2*theta + phi)

    return tangamma








# tangamma = gettangamma(k1, k2, thetaExample, phi)
# print("tan gamma : %s" %tangamma)
#
# gamma = np.arctan(tangamma)
# print("gamma : %s" %gamma)










def getaccelLocal(r, mu, outputUnits=None):
    """
    Finds local acceleration at a given point from central body. If no units given use SI
    :param r:
    :param mu:
    :param outputUnits:
    :return:
    """
    aLocal = mu/r**2

    if outputUnits is not None:
        aLocal = aLocal.to(outputUnits)

    return aLocal









# aLocal = getaccelLocal(r1, muSun, outputUnits=u.m/u.s**2)
# print("Local accel : %s" %aLocal)










def getaccel(k1, k2, gamma, s, r, mu, outputUnits=None, alwaysPositive=False):
    """
    Gets the thrusting acceleration using equations given in slides. If no units given use SI
    :param k1:
    :param k2:
    :param gamma:
    :param s:
    :param r:
    :param mu:
    :param outputUnits:
    :param alwaysPositive: Determines whether to do as in Petropoulos (with a value for n to make accel always positive)
    or to leave as usual, and a negative acceleration indicates against the spacecraft velocity
    :return accel, n: When alwaysPositive is False, n is uneccessary and is equal to 0
    """
    out = np.tan(gamma) / (2*cos(gamma))
    in1 = ( (np.tan(gamma))**2 + k1*k2**2*s + 1 )**-1

    in2_top = k2**2 * (1 - 2*k1*s)
    in2_bottom = ( np.tan(gamma)**2 + k1*k2**2*s + 1 )**2

    #accelN is normalised acceleration wrt solar accel
    accelNTemp = out * (in1 - in2_top/in2_bottom)

    if alwaysPositive is True:
        if np.sign(accelNTemp) == 1:
            n = 0
        else:
            n=1
    else:
        n=0

    accelN = (-1)**n * accelNTemp

    accel = accelN * getaccelLocal(r, mu)


    if outputUnits is not None:
        accel = accel.to(outputUnits)

    return accel, n











# accel, n = getaccel(k1, k2, gamma, s, r1, muSun, outputUnits=u.km/u.s**2, alwaysPositive=False)
# print("Acceleration : %s In %s direction" %(accel, n))












def getdtdtheta(k1, k2, gamma, r, s, mu, outputUnits=None):
    """
    Gets inverse rate of change of theta. If no units given use SI
    :param k1:
    :param k2:
    :param gamma:
    :param r:
    :param s:
    :param mu:
    :param outputUnits:
    :return dtdtheta:
    """
    top = r**3 * ( np.tan(gamma)**2 + k1*k2**2*s + 1)

    dtdtheta = np.sqrt(top/mu)

    if outputUnits is not None:
        dtdtheta = dtdtheta.to(outputUnits, equivalencies=u.dimensionless_angles())

    return dtdtheta









# dtdtheta = getdtdtheta(k1, k2, gamma, r1, s, muSun, outputUnits=u.s/u.rad)
# dtdtheta2 = dtdtheta**2
# print("dt/dtheta^2 : %s" %dtdtheta2)










def getr(k0, k1, k2, theta, phi, outputUnits=None):
    """
    Gets radius from central body with equations as shown in slides. If no units given use SI
    :param k0:
    :param k1:
    :param k2:
    :param theta:
    :param phi:
    :param outputUnits:
    :return:
    """
    r = k0 * np.exp(k1 * np.sin(k2*theta + phi))

    if outputUnits is not None:
        r = r.to(outputUnits)

    return r











# r = getr(k0, k1, k2, thetaExample, phi, outputUnits=u.km)
# print("r : %s" %r)





def getdt(dtheta, dtdtheta, outputUnits=None):
    #TODO: might ditch this function (see TOF)
    """
    Gets value for t given dtdtheta. If no units given use SI
    :param dtheta: find from table with theta_(i+1) - theta_(i)
    :param dtdtheta: from previous equations
    :param outputUnits:
    :return dt:
    """
    dt = dtdtheta * dtheta

    if outputUnits is not None:
        dt.to(outputUnits)

    return dt





# theta_2 = 1.4137167*u.rad
# dtheta = theta_2 - thetaExample
# dt = getdt(dtheta, dtdtheta, outputUnits=u.s)
# print("dt : %s" %dt)
#
# tf = 4.045019805E7 * u.s
# t0 = 4.044862610E7 * u.s
# print((tf).to(u.year))









def doTOFFunc(k1, k2, gamma, r, s, mu, dtheta, outputUnits=None):
    """
    Performs function inside the sum for TOF calcs (can be used for rectangle or simpson) If no units given use SI
    :param k1:
    :param k2:
    :param gamma:
    :param r:
    :param s:
    :param mu:
    :param dtheta:
    :param outputUnits:
    :return fx:
    """
    top = r**3 * (np.tan(gamma)**2 + k1*k2**2*s + 1)

    fx = dtheta * np.sqrt(top/mu)

    if outputUnits is not None:
        fx = fx.to(outputUnits, equivalencies=u.dimensionless_angles())

    return fx



# fx0 = doTOFFunc(k1, k2, gamma, r, s, muSun, dtheta, outputUnits=u.s)
# print("fx0 : %s" %fx0)



def doExposinStep(theta, dtheta, k0, k1, k2, phi, mu, rUnits=None, accelUnits=None, dtdthetaUnits=None, timeUnits=None):

    s = getS(k2, theta, phi)
    k1s = k1*s
    c = getC(k2, theta, phi)
    tangamma = gettangamma(k1, k2, theta, phi)
    gamma = np.arctan(tangamma)
    r = getr(k0, k1, k2, theta, phi, outputUnits=rUnits)
    accel, n = getaccel(k1, k2, gamma, s, r, mu, outputUnits=accelUnits, alwaysPositive=False)
    dtdtheta = getdtdtheta(k1, k2, gamma, r, s, mu, outputUnits=dtdthetaUnits)
    dtdtheta2 = dtdtheta**2

    fx = doTOFFunc(k1, k2, gamma, r, s, mu, dtheta, outputUnits=timeUnits)

    return [k1s, s, c, tangamma, gamma, r, accel, dtdtheta2, fx]

# dtheta = 1.4137167 * u.rad
#
# theta0 = 0 * dtheta
# point0List = doExposinStep(theta0, dtheta, k0, k1, k2, phi, muSun, rUnits=u.km, accelUnits=u.km/u.s**2, dtdthetaUnits=u.s/u.rad, timeUnits=u.s)
# fx0 = point0List[-1]
# print("Point 0 things : %s" %point0List)











def doFullExposin(k2, r1, r2, Psi, N, gamma1, mu, method="rectangle_central", thetaSplitNumber=10,
                  angleUnits=None, distanceUnits=None, accelUnits=None, timeUnits=None):
    """
    Do the exposin simulation for a specific set of parameters. To run for many N's etc, use this function repeatedly
    :param k2:
    :param r1:
    :param r2:
    :param Psi:
    :param N:
    :param gamma1:
    :param mu:
    :param method: Whether to use central rectangle method or start of section method. Simpson not yet implemented
    :param thetaSplitNumber: Number of points to use when splitting thetaBar
    :param angleUnits:
    :param distanceUnits:
    :param accelUnits:
    :param timeUnits:
    :return:
    """
    # Find thetaBar and split the theta search space using central rectangle theorem
    thetaBar = getthetaBar(Psi, N, outputUnits=angleUnits)
    thetaArrBase = np.linspace(0, thetaBar, thetaSplitNumber + 1)[0:-1]
    dtheta = thetaArrBase[1] - thetaArrBase[0]
    if method == "rectangle_lower":
        thetaArr = thetaArrBase

    elif method == "rectangle_central":
        thetaArr = thetaArrBase + dtheta/2

    # Big array to put all data in (like the one on slides)
    if angleUnits is not None:
        dataFinal = np.zeros((len(thetaArr), 11), dtype=object)
    else:
        dataFinal = np.zeros((len(thetaArr), 11), dtype=float)

    dataFinal[:, 0] = thetaArr
    # print(dataFinal)

    # Calculate/ set common parameters
    k1 = getk1(r1, r2, k2, thetaBar, gamma1)
    phi = getphi(k1, k2, gamma1, outputUnits=angleUnits)
    k0 = getk0(r1, k1, phi, outputUnits=distanceUnits)

    fxList = []
    # print(thetaArr)
    for i in range(len(thetaArr)):
        theta = thetaArr[i]

        if (timeUnits is None) or (angleUnits is None):
            dtdthetaUnits = None
        else:
            dtdthetaUnits = timeUnits/angleUnits

        k1s, s, c, tangamma, gamma, r, accel, dtdtheta2, fx = doExposinStep(theta, dtheta, k0, k1, k2, phi, mu,
                      rUnits=distanceUnits, accelUnits=accelUnits, dtdthetaUnits=dtdthetaUnits, timeUnits=timeUnits)

        # Make time the sum
        fxList.append(fx)
        t = sum(fxList)
        dt = fx
        # print(k1s, s, c, tangamma, gamma, r, accel, dtdtheta2, dt, t)
        dataFinal[i, 1:] = k1s, s, c, tangamma, gamma, r, accel, dtdtheta2, dt, t

    return dataFinal

# print(doFullExposin(k2, r1, r2, Psi, N, gamma1Example, muSun, method="rectangle_central", thetaSplitNumber=10,
#                   angleUnits=u.rad, distanceUnits=u.km, accelUnits=u.km/u.s**2, timeUnits=u.s))


def doTOFCalcMain(gamma1s, stepNumbers, k2, r1, r2, Psi, N, mu, saveName="test.npy", TOFArrUnits=u.year, printing=True):
    """
    Does the calculations for plotting
    :param gamma1s:
    :param stepNumbers:
    :param k2:
    :param r1:
    :param r2:
    :param Psi:
    :param N:
    :param mu:
    :param saveName:
    :param TOFArrUnits:
    :param printing:
    :return:
    """
    TOFGammaData = np.zeros((len(gamma1s), len(stepNumbers)))
    for i in range(len(gamma1s)):
        for j in range(len(stepNumbers)):
            gamma1 = gamma1s[i]
            stepNumber = stepNumbers[j]
            yRange = abs(gamma1s[-1] - gamma1s[0])
            diff = abs((gamma1 - gamma1s[0]))
            percent = diff/yRange * 100
            if printing: print("gamma %s, step %s. percent completion: %s" % (np.round(gamma1, decimals=4), int(stepNumber), int(percent)))

            dataStep = doFullExposin(k2, r1, r2, Psi, N, gamma1, mu, thetaSplitNumber=stepNumber)
            TOF = (dataStep[-1, -1] * u.s).to(TOFArrUnits).value
            TOFGammaData[i, j] = TOF

    np.save(saveName, TOFGammaData)

def find_nearest_index(array, value):
    """
    Finds index nearest to desired value from array
    :param array:
    :param value:
    :return:
    """
    array = np.asarray(array)
    idx = (np.abs(array - value)).argmin()
    return idx




# theta1 = 1 * dtheta
# point1List = doExposinStep(theta1, dtheta, k0, k1, k2, phi, muSun, rUnits=u.km, accelUnits=u.km/u.s**2, dtdthetaUnits=u.s/u.rad, timeUnits=u.s)
# fx1 = point1List[-1]
# print("Point 1 things : %s" %point1List)
#
# theta2 = 2 * dtheta
# point2List = doExposinStep(theta2, dtheta, k0, k1, k2, phi, muSun, rUnits=u.km, accelUnits=u.km/u.s**2, dtdthetaUnits=u.s/u.rad, timeUnits=u.s)
# fx2 = point2List[-1]
# print("Point 2 things : %s" %point2List)
#
# TOF = (1/3 * dtheta * (fx0 + 4*fx1 + fx2)).to(u.s, equivalencies=u.dimensionless_angles())
# print("Simpson TOF : %s" %TOF)







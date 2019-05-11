#!/usr/bin/env python
"""
orbdesign14_utils.py: Script for utility functions of ORBDEGSIGN-14 assignment
"""
__author__      = "Matthew Turnock"
__email__ = "matthew.turnock@protonmail.com"
__version__ = "1.0"
###########################################################################################################

import numpy as np
from astropy import units as u

from json_to_dict import constants

###############################################################################################################
# DV budget functions
RE = constants["RE"]

def geta(Pe, Ap, outputUnits=None):
    """
    Get SMA from Pe and Ap
    :param Pe:
    :param Ap:
    :param outputUnits:
    :return:
    """
    a = 0.5*(Pe + Ap)

    if outputUnits is not None:
        a = a.to(outputUnits)

    return a

def getV(r, a, mu=constants["muEarth"], outputUnits=None):
    """
    To get orbital velocity given radius and SMA
    :param r:
    :param a:
    :param mu:
    :param outputUnits:
    :return:
    """
    V = np.sqrt( mu*( 2/r - 1/a ) )

    if outputUnits is not None:
        V = V.to(outputUnits)

    return V

def getDV_combined(V1, V2, Di, outputUnits=None):
    """
    Used for getting circularisation DV
    :param V1: Initial velocity
    :param V2: Final velocity
    :param Di: Inclination change
    :param outputUnits:
    :return:
    """
    DV = np.sqrt(V1**2 + V2**2 - 2*V1*V2*np.cos(Di))

    if outputUnits is not None:
        DV.to(outputUnits)

    return DV

def get_GEO_DV(Pe_GTO, i_GTO, Ap_GTO, i_GEO=0*u.deg, printing=True, outputUnits=None):
    """
    Function to get GTO--> GEO DV, with inclination change
    :param Pe_GTO:
    :param i_GTO:
    :param Ap_GTO:
    :param i_GEO:
    :param printing:
    :param outputUnits:
    :return:
    """
    # Get GTO parameters
    Ap_GTO = 35786*u.km + RE
    a_GTO = geta(Pe_GTO, Ap_GTO)
    # Get GEO parameters
    a_GEO = Ap_GTO
    # Get inc change
    Di = i_GTO - i_GEO
    # Get Apogee velocities
    V1 = getV(Ap_GTO, a_GTO)
    V2 = getV(a_GEO, a_GEO)

    # Get total DV
    DV_tot = getDV_combined(V1, V2, Di)

    #Print some info
    if printing:
        print("Initial V : %s.  Final V : %s"%(V1, V2))
        print("DV with no inc change : %s" %(abs(V1-V2)))
        print("Delta V (circ + inc) : %s" %DV_tot)

    if outputUnits is not None:
        DV_tot = DV_tot.to(outputUnits)

    return DV_tot

def getDV_circ2circ(a1, a2, printing=True, outputUnits=None):
    """
    Circular orbit to circular orbit DV. Assumes lower to higher
    :param a1:
    :param a2:
    :return:
    """
    # Get transfer orbit SMA
    aT = geta(a1, a2)
    # Get circular velocities
    Vc1 = getV(a1, a1)
    Vc2 = getV(a2,a2)
    # Get Transfer velocities
    VT_Pe = getV(a1, aT)
    VT_Ap = getV(a2, aT)
    # Get trnasfer DVs
    DV1 = abs(VT_Pe - Vc1)
    DV2 = abs(VT_Ap - Vc2)

    if outputUnits is not None:
        DV1 = DV1.to(outputUnits)
        DV2 = DV2.to(outputUnits)

    # Do some printing
    if printing:
        print("From orbit a1=%s to a2=%s:" %(a1, a2))
        print("DV1 : %s" %DV1)
        print("DV2 : %s" %DV2)

    return (DV1, DV2)

def get_DV_J22(lam, inputUnits=None, outputUnits=None):
    """
    Input longitude lam as DEGREES, or specify input units. Outputs as m/s/yr if nothing put in
    :param lam:
    :param outputUnits:
    :return:
    """

    if inputUnits is not None:
        lam = lam.to(u.deg)
    else:
        lam = lam*u.deg

    DV_J22 = abs(1.7*np.sin(2*(lam - 75*u.deg)))
    DV_J22 = DV_J22*u.m/u.s/u.year

    if outputUnits is not None:
        DV_J22 = DV_J22.to(outputUnits)

    return DV_J22

def get_DV_NS(V_GEO, Deltai_NS, outputUnits=None):
    """
    Get delta v due to polar wander
    :param V_GEO:
    :param Deltai_NS:
    :param outputUnits:
    :return:
    """
    DV_NS = np.sqrt(2 * V_GEO**2 * (1 - np.cos(Deltai_NS)))

    if outputUnits is not None:
        DV_NS = DV_NS.to(outputUnits)

    return DV_NS

###############################################################################################################
# Cost functions

def get_mProp(mDry, DV, Isp, g0=9.81*u.m/u.s**2, outputUnits=None):
    """
    Use tsiolkovsky to get propellant mass
    :param mDry:
    :param DV:
    :param Isp:
    :param g0:
    :param outputUnits:
    :return:
    """
    mProp = mDry * ( np.exp(DV/(Isp*g0)) - 1 )

    if outputUnits is not None:
        mProp = mProp.to(outputUnits)

    return mProp

def get_OCF(K, DV, Isp, g0=9.81*u.m/u.s**2):
    """
    Gets orbit cost function OCF with equation from Wertz
    :param K:
    :param DV:
    :param Isp:
    :param g0:
    :return:
    """
    OCF = (1+K) * np.exp( DV/(Isp*g0) ) - K
    OCF = (OCF.to(u.dimensionless_unscaled)).value

    return OCF

###############################################################################################################
# FOV stuff

def getW(R, beta, outputUnits=None):
    """
    Finds witdth of area for satellite to cover
    :param R: Radius of planet
    :param beta: Angle between most extreme points to cover
    :param outputUnits:
    :return:
    """
    W = 2*R*np.sin(0.5*beta)

    if outputUnits is not None:
        W = W.to(outputUnits)

    return W

def getFOV(W, r, outputUnits=None):
    """
    Finds required FOV of satellite
    :param W: width of area to cover
    :param r: orbital radius
    :param outputUnits:
    :return:
    """
    FOV = 2*np.arctan(W/(2*r))

    if outputUnits is not None:
        FOV =FOV.to(outputUnits)

    return FOV

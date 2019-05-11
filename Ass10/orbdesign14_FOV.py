#!/usr/bin/env python
"""
orbdesign14_FOV.py: Script for finding satellite FOV value of ORBDEGSIGN-14 assignment
"""
__author__      = "Matthew Turnock"
__email__ = "matthew.turnock@protonmail.com"
__version__ = "1.0"
###########################################################################################################
from Ass10.orbdesign14_utils import *
from Ass10.orbdesign14_OCF import a_GEO

beta = 126*u.deg
r = a_GEO

W = getW(RE, beta, outputUnits=u.km)
print("W : %s (=%s RE)" %(W, W/RE))

FOV = getFOV(W, r, outputUnits=u.deg)
print("FOV : %s" %FOV)

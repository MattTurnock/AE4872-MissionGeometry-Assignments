#!/usr/bin/env python

"""
geometry2_verification.py: Verification script for the GEOMETRY-2 assignment. Compares with values in slides
"""
__author__      = "Matthew Turnock"
__email__ = "matthew.turnock@protonmail.com"
__version__ = "1.0"

########################################################################################################################
from Ass6.design7_utils import npArray2LatexTable
from Ass7.geometry2_utils import *


a1 = 10*u.deg
d1 = 15*u.deg
a2 = 45*u.deg
d2 = 70*u.deg
a3 = 110*u.deg
d3 = 32*u.deg
inputArray = np.array([[a1, a2, a3, d1, d2, d3, "inside"],
                       [a1, a2, a3, d1, d2, d3, "outside"]], dtype=object)

outputArray = findAllValues(inputArray)
npArray2LatexTable(inputArray, "verification_input.txt")
npArray2LatexTable(outputArray, "verification_output.txt")
print("Given the input parameters:\n%s\nWe get output parameters:\n%s" %(inputArray, outputArray))



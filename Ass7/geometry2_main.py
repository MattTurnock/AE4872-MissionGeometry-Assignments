#!/usr/bin/env python

"""
geometry2_main.py: Main calculation script for the GEOMETRY-2 assignment.
"""
__author__      = "Matthew Turnock"
__email__ = "matthew.turnock@protonmail.com"
__version__ = "1.0"

########################################################################################################################
from Ass6.design7_utils import npArray2LatexTable
from Ass7.geometry2_utils import *

inputArray = np.array([[0*u.deg, 90*u.deg, 0*u.deg, 0*u.deg, 0*u.deg, 90*u.deg, "inside"],
                       [0*u.deg, 30*u.deg, 0*u.deg, 0*u.deg, 0*u.deg, 90*u.deg, "inside"],
                       [0*u.deg, 10*u.deg, 5*u.deg, 20*u.deg, 25*u.deg, 30*u.deg, "inside"],
                       [5*u.deg, 30*u.deg, 5*u.deg, 0*u.deg, 30*u.deg, 90*u.deg, "inside"],
                       [0*u.deg, 45*u.deg, 0*u.deg, -20*u.deg, -20*u.deg, 90*u.deg, "inside"],
                       [0*u.deg, 10*u.deg, 5*u.deg, 20*u.deg, 25*u.deg, 30, "outside"]], dtype=object)

outputArray = findAllValues(inputArray)
npArray2LatexTable(inputArray, "main_input.txt")
npArray2LatexTable(outputArray, "main_output.txt")
print(outputArray)



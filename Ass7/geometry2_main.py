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

decimals = 4

# Define input parameters
inputArray = np.array([[0*u.deg, 90*u.deg, 0*u.deg, 0*u.deg, 0*u.deg, 90*u.deg, "inside"],
                       [0*u.deg, 30*u.deg, 0*u.deg, 0*u.deg, 0*u.deg, 90*u.deg, "inside"],
                       [0*u.deg, 10*u.deg, 5*u.deg, 20*u.deg, 25*u.deg, 30*u.deg, "inside"],
                       [5*u.deg, 30*u.deg, 5*u.deg, 0*u.deg, 30*u.deg, 90*u.deg, "inside"],
                       [0*u.deg, 45*u.deg, 0*u.deg, -20*u.deg, -20*u.deg, 90*u.deg, "inside"],
                       [0*u.deg, 10*u.deg, 5*u.deg, 20*u.deg, 25*u.deg, 30*u.deg, "outside"]], dtype=object)

# Calculate output
outputArray = findAllValues(inputArray)

# Round off numbers to proper precision
outputArrayRounded = np.copy(outputArray)
for i in range(len(outputArrayRounded[:,0])):
    for j in range(len(outputArrayRounded[0,:])):
        outputArrayRounded[i, j] = outputArray[i, j].value
outputArrayRounded = np.array(outputArrayRounded, dtype=float)
outputArrayRounded = np.around(outputArrayRounded, decimals=decimals)

# Make an array that will be saved, for latex input
outputArrayToSave = np.zeros((len(outputArray), 11), dtype="<U10")
outputArrayToSave[:, 4:] = outputArrayRounded

# Add a column to the save array for angles of the problem
latexAnglesArray = np.zeros((len(outputArray), 3), dtype="<U20")
inputAngles = inputArray[:, 0:6]

for i in range(len(inputAngles[:,0])):
    for j in range(len(inputAngles[0,:])):
        inputAngles[i, j] = inputAngles[i, j].value
inputAngles = np.array(inputAngles, dtype=int)

baseAngleString = "(%s, %s)"

for i in range(len(latexAnglesArray)):
    a1, a2, a3, d1, d2, d3 = inputAngles[i]
    latexAnglesArray[i, 0] = baseAngleString %(a1, d1)
    latexAnglesArray[i, 1] = baseAngleString %(a2, d2)
    latexAnglesArray[i, 2] = baseAngleString %(a3, d3)

outputArrayToSave[:, 0:3] = latexAnglesArray

# Add a column for triangle types
triangleTypes = ["Inner"]*5 + ["Outer"]
outputArrayToSave[:,3] = triangleTypes

# Save arrays
npArray2LatexTable(inputArray, "main_input.txt")
npArray2LatexTable(outputArray, "main_output.txt")
npArray2LatexTable(outputArrayToSave, "main_output_to_latex.txt")

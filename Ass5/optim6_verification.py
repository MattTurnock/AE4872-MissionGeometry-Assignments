#!/usr/bin/env python

"""test.py: Runs verification for the Himmelblau GA, by comparing against "theoretical" minimum
for the OPTIM-6 assignment"""
__author__      = "Matthew Turnock"
__email__ = "matthew.turnock1@gmail.com"
__version__ = "1.0"

#######################################################################################################################
import numpy as np
from Ass5.runGACases import bitLengths, UBBoth, LBBoth
from Ass5.optim6_utils import himmelblau2
from Ass1.misc_utils import list_perms

# Use given grid sizes and initialise a grids list
gridSizes = [32, 128, 512]
grids = []

# Add the grid points to the grids list
for gridSize in gridSizes:
    grids.append(np.linspace(LBBoth, UBBoth, gridSize))

# Find the theoretical optimum that can be found by the GA given the grid sizes
for i in range(len(grids)):
    grid = grids[i]
    minimum = 999
    optimalPair = [0,0]
    permutations = list_perms(grid, grid)

    for permutation in permutations:
        evaluation = himmelblau2(permutation)
        if evaluation < minimum:
            minimum = evaluation
            optimalPair = permutation

    print("For bitLength: %s \nEvaluation = %s \nxns = %s\n" %(bitLengths[i], minimum, optimalPair) )

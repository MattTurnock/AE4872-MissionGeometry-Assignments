#!/usr/bin/env python

"""
part3.py: Script to perform part 3 of the DESIGN-7 assignment
"""
__author__      = "Matthew Turnock"
__email__ = "matthew.turnock@protonmail.com"
__version__ = "1.0"

########################################################################################################################

import numpy as np

iterations = 3
arrayLength = 26

monteArrayBase = np.zeros((iterations, arrayLength))
randIndices = [1, 4, 7, 9, 11, 13, 15, 17, 19, 21]

########################################################################################################################
# Do the random only case:

monteRandArray = np.copy(monteArrayBase)

# Plonk random numbers into the random slots
monteRandArray[:, randIndices] = np.random.rand(iterations, len(randIndices))
# monteRandArray[[0,1], :] = np.ones((2, 26))
print(monteRandArray) #TODO: implement the multiplication of created randoms with the representative value (use an array to do it?)


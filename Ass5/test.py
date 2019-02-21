#!/usr/bin/env python

"""test.py: Test file used for the OPTIM-6 assignment"""
__author__      = "Matthew Turnock"
__email__ = "matthew.turnock1@gmail.com"
__version__ = "1.0"

########################################################################################################################

from Ass4.optim1_utils import himmelblau

import numpy as np
import random


array = np.array([[1,2,3],[4,5,6]])
print(array)
newarray = swapCols(array, [0,-1])
print(newarray)

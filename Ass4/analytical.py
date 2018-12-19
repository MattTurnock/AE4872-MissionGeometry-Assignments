#!/usr/bin/env python

"""analytical.py: Finds minimum of Himmelblau using the analytical technique in conjunction with Newton-Raphson"""
__author__      = "Matthew Turnock"
__email__ = "matthew.turnock1@gmail.com"
__version__ = "1.0"

########################################################################################################################
from Ass4.optim1_utils import get_x1_x2
from time import clock

x1_init_1 = 1.0
x1_init_2 = 2.0

x2_limit_low = x1_limit_low = 0.0
x2_limit_high = x1_limit_high = 5.0
step_size = 1e-10

print("Following results are of form [x1, x2, h]\n")

start_time = clock()
X_atmin_1 = get_x1_x2(x1_init_1, step_size=step_size)
time_taken = clock() - start_time
print("Analytical evalutaion with x1=1: %s\nWith computation time: %s s\n" %(X_atmin_1, time_taken))

start_time = clock()
X_atmin_2 = get_x1_x2(x1_init_2, step_size=step_size)
time_taken = clock() - start_time
print("Analytical evalutaion with x1=2: %s\nWith computation time: %s s\n" %(X_atmin_2, time_taken))






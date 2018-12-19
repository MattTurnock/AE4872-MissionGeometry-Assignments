#!/usr/bin/env python

"""numerical.py: Finds minimum of Himmelblau using Monte-Carlo and Grid-Search techniques"""
__author__      = "Matthew Turnock"
__email__ = "matthew.turnock1@gmail.com"
__version__ = "1.0"

########################################################################################################################
from Ass4.optim1_utils import do_monte_himmel, get_himmel_grid, himmelblau
from random import seed
import numpy as np
import time

###########################################################
# Monte Carlo

mini = 0.0
maxi = 5.0
no_runs = 3
no_samples = [50,100,500]

seed(a=100)
print("===============================MonteCarlo Runs=========================")
Xs, Xs_min = do_monte_himmel(10, mini=mini, maxi=maxi)
time_start=1.0
time_end=1.0
time_taken=1.0
for sample in no_samples:
    time_start = time.clock()
    for run in range(no_runs):
        Xs, Xs_min = do_monte_himmel(sample, mini=mini, maxi=maxi)
        print("For N = %s and run %s, X at minimum is: %s" %(sample, run, Xs_min))
    time_end = time.clock()
    time_taken = time_end-time_start
    print("With computation time for all runs = %s s\nAnd average time per run = %s s\n" %(time_taken, time_taken/no_runs))

###########################################################
# Grid search
print("\n===========================Grid Runs=================================")

for sample in no_samples:
    time_start = time.clock()

    grid = get_himmel_grid(sample, mini=mini, maxi=maxi)
    no_evals = len(grid)
    Xs = np.zeros((no_evals, 3))
    Xs[:,0:2] = grid
    Xs[:,2] = himmelblau(Xs[:,0], Xs[:,1])

    index_min = np.argmin(Xs[:, -1])
    Xs_min = Xs[index_min]

    time_end = time.clock()
    time_taken = time_end - time_start
    print("For N = %s, X at minimum is: %s" %(no_evals, Xs_min))
    print("With computation time = %s s\n" %time_taken)











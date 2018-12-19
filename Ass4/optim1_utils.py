#!/usr/bin/env python

"""optim1_utils.py: General utilities used for the OPTIM-1 assignment"""
__author__      = "Matthew Turnock"
__email__ = "matthew.turnock1@gmail.com"
__version__ = "1.0"

########################################################################################################################
from Ass1.misc_utils import list_perms
import numpy as np
from random import uniform
np.set_printoptions(suppress=True)
###############################################################################################
# Analytical himmelblau utils
# Define g(x1)numerical
def g(x1):
    return x1**4 - 22*x1**2 + x1 + 114

# Derivative of g(x1)
def g_prime(x1):
    return 4*x1**3 - 44*x1 + 1

# Define himmelblau
def himmelblau(x1, x2):
    return (x1**2 + x2 - 11)**2 + (x1 + x2**2 - 7)**2

# x2 as a function of x1
def x2(x1):
    return 11 - x1**2

# gets x1, x2 and h for a given initial x1
def get_x1_x2(x1_init, g=g, g_prime=g_prime, step_size=1e-5):
    x1_at_min, der_g_at_min = newtons_method(g, g_prime, x1_init, step_size)
    x2_at_min = x2(x1_at_min)
    h_at_min = himmelblau(x1_at_min, x2_at_min)
    X_atmin = [x1_at_min, x2_at_min, h_at_min]
    return X_atmin

###############################################################################################
# Monte carlo utils

def do_monte_himmel(no_samples, mini=0.0, maxi=5.0):
    Xs = np.zeros([no_samples, 3])
    for i in range(len(Xs)):
        Xs[i][0] = uniform(mini, maxi)
        Xs[i][1] = uniform(mini, maxi)
        himmel_X = himmelblau(Xs[i][0], Xs[i][1])
        Xs[i][2] = himmel_X

    index_min = np.argmin(Xs[:,-1])
    X_min = Xs[index_min]

    return Xs, X_min

###############################################################################################
# Grid search utils

def get_himmel_grid(combos, mini=0.0, maxi=5.0):

    number = int(np.sqrt(combos))
    x1s = x2s = np.linspace(mini,maxi, number)
    combos = np.array(list_perms(x1s, x2s))
    return combos

#################################################################################################
#misc utils
# A couple of functions for doing newtion-raphson
def dx(f, x):
    return abs(0-f(x))
def newtons_method(f, df, x0, e, printing=False):
    delta = dx(f, x0)
    while delta > e:
        x0 = x0 - f(x0)/df(x0)
        delta = dx(f, x0)
    if printing:
        print('Root is at: ', x0)
        print('f(x) at root is: ', f(x0))

    return x0, f(x0)
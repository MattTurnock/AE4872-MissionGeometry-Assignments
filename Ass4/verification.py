#!/usr/bin/env python

"""verification.py: Performs verification of Newton-Raphson"""
__author__      = "Matthew Turnock"
__email__ = "matthew.turnock1@gmail.com"
__version__ = "1.0"

########################################################################################################################
from Ass4.optim1_utils import newtons_method
import numpy as np

########################################################################################
# Newton-Raphson verification

def f(x):
    return (x+4)*(x-7)

def f_prime(x):
    return 2*x - 3

e = 1e-10
x0s = np.linspace(-10,10, 100)

all = []
for x0 in x0s:
    all.append(newtons_method(f, f_prime, x0, e, printing=False)[0])
all = np.unique(np.around(np.array(all)))

print("Found roots for polynomial are: %s" %all)


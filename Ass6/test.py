#!/usr/bin/env python

"""
test.py:
"""
from matplotlib import pyplot as plt
import numpy as np
from astropy import units as u

x = [1*u.m, 2*u.m, 3*u.m]
for i in range(len(x)):
    x[i] = x[i].value
print(x)

plt.hist(x)
# print()
plt.show()
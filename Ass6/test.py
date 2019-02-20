#!/usr/bin/env python

"""
test.py:
"""
from __future__ import division

from matplotlib import pyplot as plt
import numpy as np
from astropy import units as u

# x = [1*u.m, 2*u.m, 3*u.m]
# for i in range(len(x)):
#     x[i] = x[i].value
# print(x)
#
# plt.hist(x)
# # print()
# plt.show()



import numpy as np
from matplotlib import pyplot as plt

# For the explanation, I simulate the data :
N=1000
data = np.random.randn(N)
# But in reality, you would read data from file, for example with :
#data = np.loadtxt("data.txt")

# Empirical average and variance are computed
avg = np.mean(data)
var = np.var(data)
# From that, we know the shape of the fitted Gaussian.
pdf_x = np.linspace(np.min(data),np.max(data),100)
pdf_y = 1.0/np.sqrt(2*np.pi*var)*np.exp(-0.5*(pdf_x-avg)**2/var)

# Then we plot :
plt.figure()
plt.hist(data,100,normed=True)
plt.plot(pdf_x,pdf_y,'k--')
plt.legend(("Fit","Data"),"best")
plt.show()

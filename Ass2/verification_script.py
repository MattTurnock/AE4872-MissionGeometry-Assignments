from json_to_dict import constants
import numpy as np
from astropy import units as u
from astropy.table import Table
pi = np.pi

from Ass2.repeating_orbit_utils import do_a_iterations

inc = 28*u.deg
i_jk_all = [[inc, 14],
            [inc, 43/3],
            [inc, 29/2],
            [inc, 59/4],
            [inc, 59/4],
            [inc, 74/5],
            [inc, 15]]


for pair in i_jk_all:
    i, j_k = pair
    a, iters = do_a_iterations(j_k, i, da_precision=100, return_its=True)
    print(a - constants["RE"], iters)






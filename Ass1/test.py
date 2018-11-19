import json
import os
import numpy as np
import pylatex
AE4878_path = os.path.dirname(os.path.realpath(__file__))
with open(os.path.join(AE4878_path,'constants.json')) as handle:
    course_constants = json.loads(handle.read())


# for i in course_constants:
#     print("Name: %s\nUnit: %s\nValue:%s" %(i, course_constants[i]['unit'][2:], course_constants[i]['val']))

def thing(E_i, e, M):
    # top = M - E_i + e*np.sin(E_i)
    # bottom = 1 - e*np.cos(E_i)
    #
    # E_ip1 = E_i + top/bottom
    ECC = e
    return E_i + (M - E_i + ECC*np.sin(E_i)/(1-ECC*np.cos(E_i)))

M = np.deg2rad(24.06608426)
e = 0.000731104
E_i = np.deg2rad(22.)

for i in range(0,10):
    #print(np.rad2deg(E_i))
    E_i = thing(E_i, e, M)

from Ass1.kep_orbit_utils import get_precision, M2E

print(M2E(0.3886468432, 0.3))




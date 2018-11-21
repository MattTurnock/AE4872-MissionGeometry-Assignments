import json
import os
import numpy as np
# import pylatex
# AE4878_path = os.path.dirname(os.path.realpath(__file__))
# with open(os.path.join(AE4878_path,'constants.json')) as handle:
#     course_constants = json.loads(handle.read())


from Ass1 import kep_orbit_utils
def get_sigfigs(number):
    '''Return the number of significant figures of the input digit string'''
    digits = str(number)
    integral, _, fractional = digits.partition(".")

    if fractional:
        return len((integral + fractional).lstrip('0'))
    else:
        return len(integral.strip('0'))

print(get_sigfigs('0.0005467576000'))
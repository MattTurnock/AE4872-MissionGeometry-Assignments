import json
import os
import numpy as np
# import pylatex
# AE4878_path = os.path.dirname(os.path.realpath(__file__))
# with open(os.path.join(AE4878_path,'constants.json')) as handle:
#     course_constants = json.loads(handle.read())


from Ass1 import kep_orbit_utils, misc_utils
a = np.rad2deg(kep_orbit_utils.E2M(np.deg2rad(239.5991), 0.0011219))
print(a)

print(np.deg2rad(239.5991) - 0.0011219*np.sin(np.deg2rad(239.5991)))
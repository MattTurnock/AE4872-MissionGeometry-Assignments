import json
import os
import numpy as np
import pylatex
AE4878_path = os.path.dirname(os.path.realpath(__file__))
with open(os.path.join(AE4878_path,'constants.json')) as handle:
    course_constants = json.loads(handle.read())


# for i in course_constants:
#     print("Name: %s\nUnit: %s\nValue:%s" %(i, course_constants[i]['unit'][2:], course_constants[i]['val']))


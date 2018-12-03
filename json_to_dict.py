from astropy import units as u
import json
import os
AE4878_path = os.path.dirname(os.path.realpath(__file__))
with open(os.path.join(AE4878_path, 'constants.json')) as  handle:
    course_constants = json.loads(handle.read())

constants = {}
for constant in course_constants:
    #print(course_constants[constant]['val'], type(course_constants[constant]['val']))
    #print(u.Unit(course_constants[constant]['unit']), type(u.Unit(course_constants[constant]['unit'])))
    constants[constant] = course_constants[constant]['val'] * u.Unit(course_constants[constant]['unit'])
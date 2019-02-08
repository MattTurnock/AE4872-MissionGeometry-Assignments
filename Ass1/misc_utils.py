#!/usr/bin/env python

"""test.py: Define a number of generic mathematical utilities"""
__author__      = "Matthew Turnock"
__email__ = "matthew.turnock1@gmail.com"
__version__ = "1.1"

#######################################################################################################################

import numpy as np
import itertools
#lists all permutations of 2 lists
def list_perms(list1, list2):
    perms_temp = list(itertools.product(list1, list2))
    perms=[]
    for perm in perms_temp:
        perms.append(list(perm))
    return perms

#Returns integer of the number of significant figures of a number
def get_sigfigs(number):
    '''Return the number of significant figures of the input digit string'''
    digits = str(number)
    if '-' in digits:
        digits = digits.replace('-','')
    integral, _, fractional = digits.partition(".")

    if fractional:
        return len((integral + fractional).lstrip('0'))
    else:
        return len(integral.strip('0'))

#returns a negative angle to the [0,360] domain
def angle2positive(angle):
    return 2 * np.pi + angle

#Finds the number of digits after the decimal point
def get_precision(number):
    number_str = str(number)
    return number_str[::-1].find('.')

#Rounds number to given significant figures
def round2sigfig(number, sigfigs, return_float=False):
    toround = '%.' + str(sigfigs) + 'g'
    rounded_str =  '%s' % float(toround % number)
    if '-' in rounded_str:
        addon = 2
    else:
        addon =1

    if ('.' in rounded_str) and (len(rounded_str)!=sigfigs+addon):
        rounded_str = rounded_str + '0'*(sigfigs - len(rounded_str)+addon)

    if return_float:
        return float(rounded_str)
    else: return rounded_str

#Returns the largest number of significant figures in a list
def get_large_sigfig_in_list(lst):
    sigfigs = []
    for i in lst:
        sigfigs.append(get_sigfigs(str(i)))

    return max(sigfigs)

#rounds a list to a certain number of significant figures
def round_list_sigfigs(lst, sigfigs):
    newlst = []
    for i in lst:
        newlst.append(round2sigfig(i, sigfigs))

    return newlst


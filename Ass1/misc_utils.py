import numpy as np

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

def angle2positive(angle):
    return 2 * np.pi + angle

def get_precision(number):
    number_str = str(number)
    return number_str[::-1].find('.')

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


def get_large_sigfig_in_list(lst):
    sigfigs = []
    for i in lst:
        sigfigs.append(get_sigfigs(str(i)))

    return max(sigfigs)

def round_list_sigfigs(lst, sigfigs):
    newlst = []
    for i in lst:
        newlst.append(round2sigfig(i, sigfigs))

    return newlst


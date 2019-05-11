from Ass10.orbdesign14_utils import *
from Ass10.orbdesign14_OCF import a_GEO

beta = 126*u.deg
r = a_GEO

W = getW(RE, beta, outputUnits=u.km)
print("W : %s (=%s RE)" %(W, W/RE))

FOV = getFOV(W, r, outputUnits=u.deg)
print("FOV : %s" %FOV)
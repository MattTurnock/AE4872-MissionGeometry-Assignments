#!/usr/bin/env python
"""
orbdesign14_DVBudget.py: Script for total DV budget of ORBDEGSIGN-14 assignment
"""
__author__      = "Matthew Turnock"
__email__ = "matthew.turnock@protonmail.com"
__version__ = "1.0"
###########################################################################################################
from Ass10.orbdesign14_utils import *

printing = False

if printing: print("\n===========TRANSFER DV1=============")
Pe_GTO = 185*u.km + RE
Ap_GTO = 35786*u.km + RE
a_GTO = geta(Pe_GTO, Ap_GTO, outputUnits=u.km)
V_Pe = getV(Pe_GTO, a_GTO)
V_Ap = getV(Ap_GTO, a_GTO)
DV_T_1 = abs(V_Ap - V_Pe)
if printing:
    print(V_Pe, V_Ap)
    print("DV_T1 Transfer : %s" %DV_T_1)

    print("\n=========Circ_US===========")
# Get US DV
Pe_GTO_US = Pe_GTO
i_GTO_US = 28.7*u.deg
DV_T_US = get_GEO_DV(Pe_GTO_US, i_GTO_US, Ap_GTO, printing=printing)

if printing: print("\n=========Circ_EU===========")
# Get Euro DV
Pe_GTO_EU = Pe_GTO
i_GTO_EU = 5.22*u.deg
DV_T_EU = get_GEO_DV(Pe_GTO_EU, i_GTO_EU, Ap_GTO,printing=printing)

if printing: print("\n=========Circ_IND===========")
# Get Euro DV
Pe_GTO_IND = Pe_GTO
i_GTO_IND = 13.7199*u.deg
DV_T_IND = get_GEO_DV(Pe_GTO_IND, i_GTO_IND, Ap_GTO,printing=printing)

if printing: print("\n=========GRAVEYARD============")
a1 = 35786*u.km + RE
a2 = a1 + 300*u.km
DVE_1, DVE_2 = getDV_circ2circ(a1, a2, outputUnits=u.m/u.s, printing=printing)
DVE = DVE_1 + DVE_2
if printing:
    print("Total DVE (EoL) : %s" %DVE)

    print("\n=============Perturbations================")
lam = 180*u.deg
DV_EW = get_DV_J22(lam, inputUnits=u.deg)

V_GEO = getV(a1, a1, outputUnits=u.km/u.s)
Deltai_NS_peryear = 0.8*u.deg
DV_NS = (get_DV_NS(V_GEO, Deltai_NS_peryear, outputUnits=u.m/u.s))/u.year

DV_M_yearly = DV_EW  + DV_NS

years = 10*u.year
DV_M = DV_M_yearly * years

if printing:
    print("DV for J22 term (E-W) at %s longitude : %s" %(lam, DV_EW))
    print("DV for 3rd body (N-S) : %s" %DV_NS)
    print("DV for all maintenance components : %s" %DV_M_yearly)
    print("DV for all maintenance over %s years : %s" %(years, DV_M))
    print(V_GEO)

    print("\n============TOTALS=============")

DVtot_US = DV_T_US + DV_M + DVE
DVtot_EU = DV_T_EU + DV_M + DVE
DVtot_IND = DV_T_IND + DV_M + DVE

if printing:
    print("Total DV for US launch  : %s" %DVtot_US)
    print("Total DV for EU launch  : %s" %DVtot_EU)
    print("Total DV for IND launch : %s" % DVtot_IND)

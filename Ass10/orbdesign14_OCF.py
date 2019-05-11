from Ass10.orbdesign14_utils import *
from Ass10.orbdesign14_DVBudget import DVtot_US, DVtot_EU, Pe_GTO, Ap_GTO, DVE, DV_M, DV_T_US

printing=False

##############################################################################################################
# Total spacecraft mass budgeting - first estimate using SMAD table

if printing: print("\n========SC_MASS=============")
# Set mass values and fractions
m_pl = 250*u.kg #Instrument mass
PL_frac = 0.32 # Uses high Earth
# Calculate total spacecraft dry mass
m_sc_dry = m_pl / PL_frac
if printing: print("Dry mass : %s" %m_sc_dry)

# Calculate prop and wet mass
Isp_sc = 300*u.s # Uses approximate value of liquid fueled apogee engine (Wertz)

mProp_US = get_mProp(m_sc_dry, DVtot_US, Isp_sc)
mProp_EU = get_mProp(m_sc_dry, DVtot_EU, Isp_sc)

m_sc_wet_US = m_sc_dry + mProp_US
m_sc_wet_EU = m_sc_dry + mProp_EU

if printing:
    print("US propellant mass : %s. Results in total mass %s" %(mProp_US, m_sc_wet_US))
    print("EU propellant mass : %s. Results in total mass %s" %(mProp_EU, m_sc_wet_EU))

    print("\n============LAUNCHER_DV===============")

a_LEO = Pe_GTO
a_GEO = Ap_GTO

DV_launcher, waste = getDV_circ2circ(a_LEO, a_GEO, outputUnits=u.km/u.s, printing=printing)

if printing:
    print("DV1 above is the launcher DV")

    print("\n============OCF===============")

# print(DV_launcher, DV_T_US, DV_M, DVE)

K = 0.1
Isp_launcher = 348*u.s

OCF_launcher = get_OCF(K, DV_launcher, Isp_launcher)
OCF_T_US = get_OCF(K, DV_T_US, Isp_sc)
OCF_M = get_OCF(K, DV_M, Isp_sc)
OCF_E = get_OCF(K, DVE, Isp_sc)
OCF_tot = OCF_launcher * OCF_T_US * OCF_M * OCF_E

printBase = "%s OCF (with %s Isp, %s DV) : %s"
if printing:
    print(printBase %("Launcher LEO to GTO", Isp_launcher, DV_launcher, OCF_launcher ))
    print(printBase %("SC circularise", Isp_sc, DV_T_US, OCF_T_US ))
    print(printBase %("Maintenance", Isp_sc, DV_M, OCF_M))
    print(printBase %("EoL", Isp_sc, DVE, OCF_E))
    print("Total OCF : %s " %(OCF_tot))




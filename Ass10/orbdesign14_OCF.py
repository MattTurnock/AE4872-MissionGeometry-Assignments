from Ass10.orbdesign14_utils import *
from Ass10.orbdesign14_DVBudget import DVtot_US, DVtot_EU, Pe_GTO, Ap_GTO

##############################################################################################################
# Total spacecraft mass budgeting - first estimate using SMAD table

print("\n========SC_MASS=============")
# Set mass values and fractions
m_pl = 250*u.kg #Instrument mass
PL_frac = 0.32 # Uses high Earth
# Calculate total spacecraft dry mass
m_sc_dry = m_pl / PL_frac
print("Dry mass : %s" %m_sc_dry)

# Calculate prop and wet mass
Isp_sc = 300*u.s # Uses approximate value of liquid fueled apogee engine (Wertz)

mProp_US = get_mProp(m_sc_dry, DVtot_US, Isp_sc)
mProp_EU = get_mProp(m_sc_dry, DVtot_EU, Isp_sc)

m_sc_wet_US = m_sc_dry + mProp_US
m_sc_wet_EU = m_sc_dry + mProp_EU

print("US propellant mass : %s. Results in total mass %s" %(mProp_US, m_sc_wet_US))
print("EU propellant mass : %s. Results in total mass %s" %(mProp_EU, m_sc_wet_EU))

print("\n============LAUNCHER_DV===============")

a_LEO = Pe_GTO
a_GEO = Ap_GTO

DV_launcher, waste = getDV_circ2circ(a_LEO, a_GEO, outputUnits=u.km/u.s, printing=True)

print("DV1 above is the launcher DV")







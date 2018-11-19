#Function definitions for Keplerian --> Cartesian conversion
import numpy as np
from Ass1 import kep_orbit_utils

#Currently uses vertical numpy vectors, but they can be a bit fiddly, so inputs/outputs are simple lists
def Kep2Cart(state_kep, mu=kep_orbit_utils.mu_Earth, given_angle='theta'):




    state_kep = list(np.float128(state_kep))
    SMA, ECC, INC, Omega, omega, angle = state_kep

    # if given_angle != theta, then convert it for further calcs
    if given_angle == 'theta':
        theta = angle
    elif given_angle =='E':
        theta = kep_orbit_utils.E2theta(angle, ECC)
    elif given_angle == 'M':
        theta = kep_orbit_utils.E2theta(kep_orbit_utils.M2E(angle, ECC), ECC)
    else: print('INPUT ERROR ON given_angle VARIABLE, SHOULD BE theta, E or M')

    sin = np.sin
    cos = np.cos
    sqrt = np.sqrt

    r = kep_orbit_utils.get_orbit_radius(SMA,ECC,theta)

    #Define intermediary terms
    l1 = cos(Omega)*cos(omega) - sin(Omega)*sin(omega)*cos(INC)
    l2 = -cos(Omega)*sin(omega) - sin(Omega)*cos(omega)*cos(INC)

    m1 = sin(Omega)*cos(omega) + cos(Omega)*sin(omega)*cos(INC)
    m2 = -sin(Omega)*sin(omega) + cos(Omega)*cos(omega)*cos(INC)

    n1 = sin(omega)*sin(INC)
    n2 = cos(omega)*sin(INC)

    # =============CALCULATE POSITION==============
    zeta_eta = np.vstack([r*cos(theta), r*sin(theta)])
    T = np.array([[l1,  l2],
                  [m1,  m2],
                  [n1,  n2]])

    position = np.dot(T, zeta_eta)

    # =============CALCULATE VELOCITY==============
    H = kep_orbit_utils.get_H(mu, SMA, ECC)

    xdot = mu/H*( -l1*sin(theta) + l2*(ECC + cos(theta)))
    ydot = mu/H*(-m1*sin(theta) + m2*(ECC + cos(theta)))
    zdot = mu/H*(-n1*sin(theta) + n2*(ECC + cos(theta)))

    velocity = np.vstack((xdot, ydot, zdot))

    state =np.vstack((position,velocity))
    state=list(state.T[0])
    return state








# Do an example if main
if __name__ == '__main__':


    SMA = 6787746.891E-3
    ECC = 0.000731104
    INC = np.deg2rad(51.68714486)
    Omega = np.deg2rad(127.5486706)
    omega = np.deg2rad(74.21987137)
    theta = np.deg2rad(24.10027677)
    state_kep = SMA, ECC, INC, Omega, omega, theta
    mu = kep_orbit_utils.mu_Earth

    Kep = Kep2Cart(state_kep, mu)
    print(Kep)


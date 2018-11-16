#Function definitions for converting Cartesian --> Keplerian
import numpy as np
import kep_orbit_utils

def get_signs_hats(N, ECCvector, rvector, hvector):
    Nhat = N / np.linalg.norm(N)
    ECChat = ECCvector / np.linalg.norm(ECCvector)
    rhat = rvector/ np.linalg.norm(rvector)

    if np.dot(np.cross(Nhat, ECCvector, axis=0).T[0], hvector.T[0]) > 0:
        sign_omega = 1
    else:
        sign_omega = -1

    if np.dot(np.cross(ECCvector, rvector, axis=0).T[0], hvector.T[0]) > 0:
        sign_theta = 1
    else:
        sign_theta = -1

    return sign_omega, sign_theta, Nhat, ECChat, rhat

def Cart2Kep(rvector, Vvector, mu):
    #basic/ derivative parameters
    r = np.linalg.norm(rvector)
    V = np.linalg.norm(Vvector)
    hvector = kep_orbit_utils.get_hvector(rvector, Vvector)
    N = kep_orbit_utils.get_N(hvector)

    #some easy orbital elements to calculate
    SMA = kep_orbit_utils.get_SMA(r, mu, V)
    ECCvector = kep_orbit_utils.get_ECCvector(Vvector, hvector, rvector, mu)
    ECC = np.linalg.norm(ECCvector)
    INC = kep_orbit_utils.get_INC(hvector)

    #Omega, omega and theta
    Nx, Ny, Nz = N[0][0], N[1][0], N[2][0]
    Nxy = np.sqrt(Nx**2 + Ny**2)

    Omega = np.arctan2(Ny/Nxy, Nx/Nxy)

    sign_omega, sign_theta, Nhat, ECChat, rhat = get_signs_hats(N, ECCvector, rvector, hvector)

    omega = sign_omega * np.arccos(np.dot(ECChat.T[0], Nhat.T[0]))
    theta = sign_theta * np.arccos(np.dot(rhat.T[0], ECChat.T[0]))

    state = SMA, ECC, INC, Omega, omega, theta
    return state

# Do an example if main
if __name__ == '__main__':


    x = -2700816.14E-3
    y = -3314092.8E-3
    z = 5266346.42E-3
    rvector = np.vstack((x,y,z))
    xdot = 5168.606550E-3
    ydot = -5597.546618E-3
    zdot = -868.878445E-3
    Vvector = np.vstack((xdot, ydot, zdot))
    mu = kep_orbit_utils.mu_Earth

    Cart = Cart2Kep(rvector, Vvector, mu)
    print(Cart)
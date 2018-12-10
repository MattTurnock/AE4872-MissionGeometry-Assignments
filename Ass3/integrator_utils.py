import numpy as np
from json_to_dict import constants
from matplotlib import pyplot as plt
from Ass1.Kep2Cart_utils import Kep2Cart
from Ass1.Cart2Kep_utils import Cart2Kep
pi = np.pi


#############################################################################################################################
#Misc things


def get_A_car(X):
    A = np.array([[0, 1],
                  [0, 0]])

    return A

def get_fulldata_car(Xall, Xdots):
    fulldata = np.zeros([len(Xall[:, 0]), 4])
    fulldata[:, 0] = Xall[:, 0]
    fulldata[:, 1:3] = Xall[:, 1:3]
    fulldata[:, 3] = Xdots[:, 2]

    return fulldata

def get_integrator_table_car(X0, dts, t_end, tablen, B,x_end_true, v_end_true, t0=0,method="RK4"):


    euler_table = np.zeros([len(dts),tablen])
    euler_table[:, 0] = dts
    #print(euler_table)
    for q in range(len(dts)):
        dt = dts[q]
        ts = np.arange(t0, t_end+dt, dt)
        Xall, Xdots = do_integration(X0, ts, get_A=get_A_car, B=B, method=method)
        X_end = get_fulldata_car(Xall, Xdots)[-1]

        n = len(ts)

        if method == "euler": factor = 1
        elif method == "RK4": factor = 4

        euler_table[q,1:3] = X_end[1:3]
        euler_table[q, 3] = x_end_true - X_end[1]
        euler_table[q, 4] = v_end_true - X_end[2]
        euler_table[q, 5] = factor*n
    return euler_table

def do_err_plot(dts, err, show=True, save=False, title=None, method="TEST", figure=True, scale="linear", legend=[]):
    if figure: plt.figure()
    if scale=="linear":
        plt.plot(dts, err)
    elif scale=="log":
        plt.loglog(dts, err)
    plt.xlabel("Timestep dt [s]")
    plt.ylabel("Final distance error [m]")
    if title != None:
        plt.title(title)

    plt.legend(legend)
    if save: plt.savefig("car_%s_error.pdf" %method)
    if show: plt.show()
    plt.grid()


def get_integrator_table_orbit(kepstate_0, t_end, dts, t0=0, method="RK4", printing=False):

    a_true = kepstate_0[0]
    e_true = kepstate_0[1]
    cartsate_0 = Kep2Cart(kepstate_0)
    X0 = cartsate_0

    fulldata = np.zeros([len(dts), len(kepstate_0)+1])
    fulldata[:,0] = dts
    fulldata[0, 1:] = kepstate_0

    for j in range(len(dts)):

        dt = dts[j]

        ts = np.arange(t0, t_end, dt)

        Xall, Xdots = do_integration(X0, ts, method=method)
        cartstate_end = Xall[-1, 1:]
        kepstate_end = Cart2Kep(cartstate_end, do_units_out=True)

        for el in range(len(kepstate_end)):
            kepstate_end[el] = kepstate_end[el].value

        fulldata[j, 1:] = kepstate_end

        if printing: print(kepstate_end)

    if method=="RK4": factor=4
    elif method=="euler": factor=1

    table = np.zeros([len(dts), 6])
    table[:,0] = dts
    table[:,1] = fulldata[:,1]
    table[:,2] = fulldata[:,2]
    table[:,3] = a_true - fulldata[:,1]
    table[:,4] = e_true - fulldata[:,2]
    table[:,5] = factor*(t_end - t0)/fulldata[:,0]
    return table



############################################################################################################################
#Euler integrator thignsa



mu = 100.
r = 12.
dt = 1

#vector to find derivative from the initial vectorino
def get_A_orbit(X, mu=constants["muEarth"].to("m ** 3 / (s * s)").value):
    r = np.linalg.norm(X[0:3])
    #print(mu)
    m = mu/r**3
    #print(m)
    A = np.array([[0, 0, 0, 1, 0, 0],
                  [0, 0, 0, 0, 1, 0],
                  [0, 0, 0, 0, 0, 1],
                  [-m, 0, 0, 0, 0, 0],
                  [0, -m, 0, 0, 0, 0],
                  [0, 0, -m, 0, 0, 0]])

    return A

def do_euler_step(Xn, dt, get_A, B):
    A = get_A(Xn)
    psi = np.matmul(A, Xn) + B
    Xnp1 = Xn + dt*psi

    return Xnp1, psi


###################################################################################################################################
# RK4 integrator things

def f_RK4(A, X, B):
    Xdot = np.matmul(A, X) + B
    return Xdot

# get k1-4 given the vector
def get_ks(X, dt, get_A, B, return_As=False):

    #print(X)
    A1 = get_A(X)
    #print(A1)
    k1 = f_RK4(A1, X, B)
    #print(k1)
    #print(A1, k1)

    Xk2 = X + 0.5*dt*k1
    A2 = get_A(Xk2)
    k2 = f_RK4(A2, Xk2, B)
    #print(A2, k2)

    Xk3 = X + 0.5 * dt * k2
    A3 = get_A(Xk3)
    k3 = f_RK4(A3, Xk3, B)
    #print(A3, k3)

    Xk4 = X + dt*k3
    A4 = get_A(Xk4)
    k4 = f_RK4(A4, Xk4, B)
    #print(A4, k4)

    if return_As:
        return k1,k2,k3,k4, A1,A2,A3,A4
    else:
        return k1,k2,k3,k4

def get_RK4_psi(k1,k2,k3,k4):
    psi = 1/6 * (k1 + 2*k2 + 2*k3 + k4)
    #print(psi)
    return psi

def do_RK4_step(Xn, dt, get_A, B, return_ks=False):
    k1,k2,k3,k4 = get_ks(Xn, dt, get_A, B)
    psi = get_RK4_psi(k1,k2,k3,k4)
    Xnp1 = Xn + dt*psi

    if return_ks:
        ks = [k1,k2,k3,k4]
        return Xnp1, psi, ks
    else:
        return Xnp1, psi


############################################################################################################################
#combined data getter

def do_integration(X0, ts, get_A=get_A_orbit, B=None, method="RK4", return_ks=False):
    try:
        if B == None:
            B = np.zeros(np.shape(X0))
    except ValueError:
        print("There was a value error, but if your input was an array don't worry about it")


    if method == "euler":
        integrator = do_euler_step
    elif method == "RK4":
        integrator = do_RK4_step

    Xall = np.zeros([len(ts), len(X0) + 1])
    Xall[:, 0] = ts

    Xdots = np.zeros(np.shape(Xall))
    Xdots[:, 0] = ts

    kall = []

    Xn = X0
    for q in range(len(ts) - 1):
        Xall[q, 1:] = Xn

        dt = ts[q + 1] - ts[q]
        if return_ks:
            Xnp1, Xdot_n, ks = integrator(Xn, dt, get_A, B, return_ks=return_ks)
            kall.append(ks)
        else:
            Xnp1, Xdot_n = integrator(Xn, dt, get_A, B)
        #print(Xnp1)
        Xdots[q, 1:] = Xdot_n

        Xn = Xnp1

    # one more append needed for the final timestamp
    Xall[-1, 1:] = Xn
    if return_ks:
        Xdots[-1, 1:], ks = integrator(Xn, 1, get_A, B, return_ks=return_ks)[1:3]
        kall.append(ks)
    else:
        Xdots[-1, 1:]= integrator(Xn, 1, get_A, B)[1]

    if return_ks:
        return Xall, Xdots, kall
    else:
        return Xall, Xdots


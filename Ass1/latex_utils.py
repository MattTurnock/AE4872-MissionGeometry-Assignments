#Defines a few functions for creating LaTeX matrices for state vectors

def Cart2Latex(state_cart):
    x, y, z, xdot, ydot, zdot = state_cart

    lat = """
\u005Cbegin{equation}
    \label{eq:myeq}
    \u005Cbegin{bmatrix}
        x\u005C\u005C
        y\u005C\u005C
        z\u005C\u005C
        \Dot{x}\u005C\u005C
        \Dot{y}\u005C\u005C
        \Dot{z}\u005C\u005C
    \end{bmatrix}
    =
    \u005Cbegin{bmatrix}
        %s\u005C\u005C
        %s\u005C\u005C
        %s\u005C\u005C
        %s\u005C\u005C
        %s\u005C\u005C
        %s
    \end{bmatrix}
\end{equation}
    """ %(x,y,z,xdot,ydot,zdot)
    return lat

def Kep2Latex(state_kep, E='Space for E, remove', M='Space for M, remove'):
    SMA, ECC, INC, Omega, omega, theta = state_kep

    lat = """
\u005Cbegin{equation}
    \label{eq:myeq}
    \u005Cbegin{bmatrix}
        a\u005C\u005C
        e\u005C\u005C
        i\u005C\u005C
        \Omega\u005C\u005C
        \omega\u005C\u005C
        \u005Ctheta\u005C\u005C
        E\u005C\u005C
        M\u005C\u005C
    \end{bmatrix}
    =
    \u005Cbegin{bmatrix}
        %s\u005C\u005C
        %s\u005C\u005C
        %s\u005C\u005C
        %s\u005C\u005C
        %s\u005C\u005C
        %s\u005C\u005C
        %s\u005C\u005C
        %s\u005C\u005C
    \end{bmatrix}
\end{equation}
    """ %(SMA, ECC, INC, Omega, omega, theta, E, M)

    return lat

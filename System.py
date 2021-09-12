import numpy as np


def host_virus(r):
    '''
    Calculates the rates of each variable at a particular timestep for
    the (nondimensionalised form of) the host-virus system.

    Inputs:
        r - a numpy array of positions in phase space
        parameters - an array of parameter values


    '''

    # Assign variables
    x1 = r[0]
    x2 = r[1]
    ys1 = r[2]
    y1 = r[3]
    y2 = r[4]
    zs = r[5]
    z = r[6]

    # Assign parameters
    alpha = 0.5
    alphas = 2.0
    beta1 = 1.50
    beta2 = 2.00
    mu = 0.10
    gammas1 = 0.25
    gamma1 = 0.25
    gamma2 = 0.25
    nu = 0.50
    nus = 0.50
    zeta = 0.22
    zetas = 0.22
    kappa1 = 1.00
    kappa2 = 1.00

    # Calculate rates
    x1_dot = (x1*(1-x1-x2))-(x1*alpha*z)-(x1*alphas*zs)
    x2_dot = (x2*(beta1-(beta2*(x1+x2))))-(x2*alpha*z)
    ys1_dot = (alphas*zs*x1) + (mu*y1)-(gammas1*ys1)
    y1_dot = (alpha*z*x1) - (mu*y1)-(gamma1*y1)
    y2_dot = (alpha*z*x2) - (gamma2*y2)
    zs_dot = (gammas1*ys1)-(nus*alphas*zs*x1)-(zetas*zs)
    z_dot = (kappa1*gamma1*y1)+(kappa2*gamma2*y2)-(nu*alpha*z*(x1+x2))-(zeta*z)

    state = np.array([x1_dot, x2_dot, ys1_dot, y1_dot, y2_dot, zs_dot, z_dot])

    return state

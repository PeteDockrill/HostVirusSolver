'''
Main.py - Main file for the Host-Virus system simulator

Peter Dockrill, 02/03/2021

'''

import numpy as np
import System as system
import RKFunctions as rk
import Plotter as pl
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from mpl_toolkits.mplot3d import Axes3D
from Model import Model

basic_params = {'alpha': 0.5,
                'alpha_s': 2.0,
                'beta_1': 1.50,
                'beta_2': 2.00,
                'mu': 0.10,
                'gamma_1_s': 0.25,
                'gamma_1': 0.25,
                'gamma_2': 0.25,
                'nu': 0.50,
                'nu_s': 0.50,
                'zeta': 0.22,
                'zeta_s': 0.22,
                'kappa_1': 1.00,
                'kappa_2': 1.00}


def main(params: dict) -> None:
    model = Model(params=params)
    sim = model.run_model(verbose=True)
    model.plot_general_cells()


if __name__ == '__main__':
    main(basic_params)

'''
Main.py - Main file for the Host-Virus system simulator

Peter Dockrill, 02/03/2021 (zcempjd@ucl.ac.uk)

'''

import numpy as np
import System as system
import RKFunctions as rk
import Plotter as pl
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from mpl_toolkits.mplot3d import Axes3D

# Set initial conditions
start_time = 0
x10 = 1
x20 = 0.05
y1s0 = 0
y10 = 0
y20 = 0
zs0 = 0.01
z0 = 0.5

initial_conds = np.array([start_time, x10, x20, y1s0, y10, y20, zs0, z0])

# Set simulation settings
time_step = 0.1
end_time = 16000
file_path = "alpha0.5_alphas_2.0_"+str(end_time/time_step)+"s_x1"

# Instantiate simuation
sim = rk.system_solver(initial_conds, end_time, time_step, file_path)

# Plot results
special_components = ['Time', 'x1', 'ys1', 'zs']
general_components = ['Time', 'x2', 'ys1', 'zs']
pl.plot_3d(sim, special_components, file_path)

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

#file_path = "alpha0.5_alphas_2.0_"+str(end_time/time_step)+"s_x1"

# Plot results
special_components = ['Time', 'x1', 'ys1', 'zs']
general_components = ['Time', 'x2', 'ys1', 'zs']
pl.plot_3d(sim, special_components, file_path)

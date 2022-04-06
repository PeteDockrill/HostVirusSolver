import numpy as np
import RKFunctions as rk
import Plotter as pl
import pandas as pd
from typing import List


class Model():

    def __init__(self, params: dict) -> None:
        self.params = params

    def run(self, verbose: bool = False, do_save: bool = False, filepath: str = '') -> pd.DataFrame:
        '''
        Runs an instance of the model

        Inputs:
            do_save - saves a csv of model results if set to True (boolean)
            filepath - file name for model results (string)

        Outputs:
            None
        '''

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
        end_time = 8000

        # Run simulation
        sim = rk.system_solver(initial_conds, end_time, time_step, self.params, verbose)

        self.sim = sim

        if do_save:
            self.save_model(filepath)

        return self.sim

    def save_model(self, filepath: str) -> None:
        '''
        Saves the results of the model as a csv

        Input:
            filepath - filepath - file name for model results (string)

        Output:
            None
        '''

        file_name = "Simulations/"+filepath+".csv"

        try:
            self.sim.to_csv(file_name)
            print('Model results saved to file: '+file_name)
        except:
            print("Simulation results do not exist")

    def plot_specialist_cells(self, xmax: int = 8000) -> None:
        '''
        Plots time series of the 'specialist' cell populations
        '''

        components = ['x1', 'y1', 'ys1']
        labels = ['Susceptible Cells', 'Infected Cells', 'Infected Specialist Cells']

        pl.plot_components(self.sim, components, labels, xmax=xmax)

    def plot_general_cells(self, xmax: int = 8000) -> None:
        '''
        Plots time series of the 'general' cell populations
        '''

        components = ['x1', 'x2', 'y2']
        labels = ['Susceptible Specialist Cells', 'Susceptible General Cells', 'Infected General Cells']

        pl.plot_components(self.sim, components, labels, xmax=xmax)

    def plot_components(self, components: List[str], labels: List[str], xmax=8000) -> None:
        '''
        Plots three components of the system

        Inputs:
            components - components to plot (list of strings)
            labels - labels for plot titles (list of strings)
        '''

        pl.plot_components(self.sim, components, labels, xmax=xmax)


class MultiModel():

    def __init__(self, models, params):

        self.n_models = len(models)
        self.models = models
        self.params = params

    def run(self) -> List[pd.DataFrame]:

        sims = []
        for model in self.models:
            sim = model.run()
            sims.append(sim)

        self.sims = sims

        return

    def plot(self, component: str, parameter: str, parameter_values: List[float], xmax: int = 8000) -> None:

        pl.plot_models(self.sims, component, parameter, parameter_values, xmax=xmax, label='')

import numpy as np
import RKFunctions as rk
import Plotter as pl
import pandas as pd
from typing import List


class Model():

    def __init__(self, params: dict, immune = False, perturb_conds = None) -> None:
        self.params = params
        self.immune = immune
        self.perturb_conds = perturb_conds

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

        if self.perturb_conds != None:
            initial_conds = self.perturb_conds
        else:
            initial_conds = np.array([start_time, x10, x20, y1s0, y10, y20, zs0, z0])

        # Set simulation settings
        self.time_step = 0.1
        self.end_time = 2000

        # Run simulation
        if self.immune:
            sim = rk.system_solver(initial_conds, self.end_time, self.time_step, self.params, verbose, True)
        else:
            sim = rk.system_solver(initial_conds, self.end_time, self.time_step, self.params, verbose)

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

    def calc_lyapunov_exponents(self, df, variables, max_tau=50, min_dist='auto'):
        """
        Estimates the largest Lyapunov exponents for multiple variables from a time series in a pandas DataFrame.

        :param df: Pandas DataFrame containing the time series data.
        :param variables: List of column names of the variables to analyze.
        :param dt: Sampling time interval of the data.
        :param max_tau: Maximum number of time steps over which to compute the divergence.
        :param min_dist: Minimum initial distance between points to consider for each variable. If 'auto', uses 1% of the standard deviation.
        :return: Dictionary with variable names as keys and their estimated largest Lyapunov exponents as values.
        """
        
        lyapunov_exponents = {}
        for variable in variables:
            time_series = df[variable].values
            exponent = self.estimate_lyapunov_exponent(time_series, self.time_step, max_tau, min_dist)
            lyapunov_exponents[variable] = exponent

        return lyapunov_exponents
    
    def estimate_lyapunov_exponent(self, time_series, max_tau, min_dist):
            N = len(time_series)
            if min_dist == 'auto':
                min_dist = 0.01 * np.std(time_series)

            distances = np.abs(np.subtract.outer(time_series, time_series))
            np.fill_diagonal(distances, np.inf)
            neighbors = distances.argmin(axis=1)

            divergences = []
            for tau in range(1, max_tau+1):
                # Avoid division by zero or log of zero by ensuring distances are above a small threshold
                dists = time_series[tau:] - time_series[neighbors[:-tau]]
                valid_dists = dists[np.abs(dists) > 1e-10]  # Exclude very small or zero distances
                if len(valid_dists) > 0:  # Ensure there are valid distances to prevent log(0)
                    divergence = np.mean(np.log(np.abs(valid_dists)))
                    divergences.append(divergence)
                else:
                    divergences.append(0)  # Append zero divergence if no valid distances are found

            if len(divergences) > 1:  # Prevent fitting to a single point or no points
                slopes, _ = np.polyfit(np.arange(1, len(divergences)+1) * self.time_step, divergences, 1)
            else:
                slopes = 0  # Return zero if unable to compute a slope
            return slopes

    
    def estimate_lyapunov_spectrum(self):
        variables = ["x1", "x2", "ys1", "y1", "y2", "zs", "z"]  # The variables for which to calculate Lyapunov exponents
        self.lyapunov_exponents = self.estimate_lyapunov_spectrum(self.sim, variables)


    def plot_specialist_cells(self) -> None:
        '''
        Plots time series of the 'specialist' cell populations
        '''

        components = ['x1', 'y1', 'ys1']
        labels = ['Susceptible Cells', 'Infected Cells', 'Infected Specialist Cells']

        pl.plot_components(self.sim, components, labels, xmax=self.end_time)

    def plot_general_cells(self) -> None:
        '''
        Plots time series of the 'general' cell populations
        '''

        components = ['x1', 'x2', 'y2']
        labels = ['Susceptible Specialist Cells', 'Susceptible General Cells', 'Infected General Cells']

        pl.plot_components(self.sim, components, labels, xmax=self.end_time)

    def plot_components(self, components: List[str], labels: List[str]) -> None:
        '''
        Plots three components of the system

        Inputs:
            components - components to plot (list of strings)
            labels - labels for plot titles (list of strings)
        '''

        pl.plot_components(self.sim, components, labels, xmax=self.end_time)


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

    def plot(self, component: str, parameter: str, parameter_values: List[float], xmax: int = 2000) -> None:

        pl.plot_models(self.sims, component, parameter, parameter_values, xmax=xmax, label='')

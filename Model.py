import numpy as np
import RKFunctions as rk
import Plotter as pl


class Model():

    def __init__(self, params):

        self.params = params

    def run(self, verbose=False, do_save=False, filepath=''):
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

    def save_model(self, filepath):
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

        return

    def plot_specialist_cells(self, xmax=8000):
        '''
        Plots time series of the 'specialist' cell populations
        '''

        components = ['x1', 'y1', 'ys1']
        labels = ['Susceptible Cells', 'Infected Cells', 'Infected Specialist Cells']

        pl.plot_components(self.sim, components, labels, xmax=xmax)

        return

    def plot_general_cells(self, xmax=8000):
        '''
        Plots time series of the 'general' cell populations
        '''

        components = ['x1', 'x2', 'y2']
        labels = ['Susceptible Specialist Cells', 'Susceptible General Cells', 'Infected General Cells']

        pl.plot_components(self.sim, components, labels, xmax=xmax)

        return

    def plot_components(self, components, labels, xmax=8000):
        '''
        Plots three components of the system

        Inputs:
            components - components to plot (list of strings)
            labels - labels for plot titles (list of strings)
        '''

        pl.plot_components(self.sim, components, labels, xmax=xmax)

        return


class MultiModel():

    def __init__(self, models, params):

        self.n_models = models.size()
        self.models = models
        self.params = params

    def run(self):

        for i, model in enumerate(self.models):
            model.run()

        return

    def plot(self, component, parameter, parameter_values):

        pl.plot_models(self.models, component, parameter, parameter_values, label)

        return

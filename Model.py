import numpy as np
import RKFunctions as rk
import Plotter as pl


class Model():

    def __init__(self, params):

        self.params = params

    def run_model(self, do_save=False, filepath=''):
        '''
        Runs an instance of the model for the
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
        end_time = 16000

        # Run simulation
        sim = rk.system_solver(initial_conds, end_time, time_step, self.params)

        self.sim = sim

        if do_save:
            self.save_model(filepath)

        return self.sim

    def save_model(self, filepath):

        file_name = "Simulations/"+filepath+".csv"

        try:
            self.sim.to_csv(file_name)
            print('Model results saved to file: '+file_name)
        except:
            print("Simulation results do not exist")

        return

import numpy as np
import pandas as pd
import System as system


def runge_kutta_step(t_0, r_0, delta, params):
    '''
    Performs one time step of a fourth-order Runge-Kutta method.

    Arguments:
            t_0 - initial time (numpy array)
            r_0 - initial spatial coordinates (numpy array)
            delta - the time step (float)
            params - model parameters (dict)

    Returns:
            r_n - the updated position (numpy array)
    '''

    k_1 = system.host_virus(t_0, r_0, params)
    k_2 = system.host_virus(t_0 + (delta/2), r_0 + (delta*(k_1/2)), params)
    k_3 = system.host_virus(t_0 + (delta/2), r_0 + (delta*(k_2/2)), params)
    k_4 = system.host_virus(t_0 + delta, r_0 + (delta*k_3), params)
    r_n = r_0 + ((delta/6)*(k_1 + (2*k_2) + (2*k_3) + k_4))
    t_n = t_0 + delta

    return r_n


def system_solver(initial_conditions, end, delta, params, verbose, do_save=False, filepath=''):
    '''
    Produces a simulated dataset for the system

    Arguments:
            Compulsory:
                initial_conditions - an array of the initial time (index 0) and positions (indices 1-3) (numpy array)
                end - the end point of the simultion (float)
                delta - the time step of the simulation (float)
                params - parameters for the system (dict)

            Optional:
                do_save - model output saved as csv if set to True
                filepath - filepath for simulation csv file

    Returns:
            solutions - the time and spatial coordinates for the solved system (pandas dataframe)
    '''

    times = np.arange(initial_conditions[0]+delta, end, step=delta)
    dimension = len(initial_conditions)-1
    solutions = np.array([initial_conditions])

    print("Simulation beginning for "+str(len(times))+" time steps.")
    print("Number of dimensions: "+str(dimension))

    for i, time in enumerate(times):
        new_position = runge_kutta_step(time, np.array(solutions)[i, 1:], delta, params)
        new_state = np.array([time, new_position[0], new_position[1], new_position[2], new_position[3], new_position[4], new_position[5], new_position[6]])
        solutions = np.vstack((solutions, new_state))

        if verbose:
            progress_percentage = round(i/len(times)*100, 2)
            if progress_percentage % 5 == 0:
                print("Simulation "+str(progress_percentage)+"\% complete")

    solutions_dict = {'Time': solutions[:, 0], 'x1': solutions[:, 1], 'x2': solutions[:, 2], 'ys1': solutions[:, 3], 'y1': solutions[:, 4], 'y2': solutions[:, 5], 'zs': solutions[:, 6], 'z': solutions[:, 7]}
    solutions = pd.DataFrame(solutions_dict)

    print("Simulation complete!")

    return solutions

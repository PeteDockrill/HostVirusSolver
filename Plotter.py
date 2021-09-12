import matplotlib.pyplot as plt
import matplotlib.animation as animation
from mpl_toolkits.mplot3d import Axes3D
import numpy as np


def plot_components(input):
    '''
    Plots the spatial coordinates as a function of time

    Arguments:
            input - the time and spatial coordinates of the system (panda dataframe)
    '''

    fig = plt.figure(figsize=(18, 12))

    ax1 = fig.add_subplot(311)
    ax1.plot(input.Time, input.X)
    ax1.grid(True)
    ax1.set_ylabel('x', fontsize=12)

    ax2 = fig.add_subplot(312)
    ax2.plot(input.Time, input.Y)
    ax2.grid(True)
    ax2.set_ylabel('y', fontsize=12)

    ax3 = fig.add_subplot(313)
    ax3.plot(input.Time, input.Z)
    ax3.grid(True)
    ax3.set_ylabel('z', fontsize=12)
    ax3.set_xlabel('Time', fontsize=12)

    plt.setp(ax1.get_xticklabels(), visible=False)
    plt.setp(ax2.get_xticklabels(), visible=False)

    plt.show()


def plot_3d(input, filepath):
    """
    Plots the system solution in 3D

    Arguments:
            input - the time and spatial coordinates of the system (panda dataframe)
            filepath - the main part of the simulation name (string)
    """

    special_components = input[['Time', 'x1', 'ys1', 'zs']]

    fig = plt.figure(figsize=(18, 18))
    ax = fig.add_subplot(111, projection='3d')
    ax.plot(input.ys1, input.x1, input.zs, lw=0.5)
    ax.set_xlabel(r"$y^s_1$", fontsize=12)
    ax.set_ylabel(r"$x_2$", fontsize=12)
    ax.set_xlim(1.2, 0)
    ax.set_zlabel(r"$z^s$", fontsize=12)
    # ax.set_title(system_name)
    plt.axis('off')
    fig = plt.gcf()
    plt.draw()
    plt.savefig(filepath+'_3D_plot.jpg', dpi=400)
    plt.show()
    #plt.savefig('Rossler_no_background.pdf', dpi = 500)


def plot_all_3d(input, components, filepath):
    '''
    Plots the system solution in 3D alongside the timeseries of the specified components

    Arguments:
            input - the time and spatial coordinates of the system (panda dataframe)
            components - the components of the system to plot (list of strings)
            filepath - the main part of the simulation name (string)
    '''

    #special_components = input[['Time','x2','ys1','zs']]

    fig = plt.figure(figsize=(14, 8))

    # Plot components
    axX = plt.subplot2grid(shape=(3, 2), loc=(0, 1), rowspan=1, colspan=1)
    axX.plot(input[components[0]], input[components[1]], lw=0.75)
    axX.grid(True)
    axX.set_ylabel(r"$x_1$", fontsize=12)

    axY = plt.subplot2grid(shape=(3, 2), loc=(1, 1), rowspan=1, colspan=1)
    axY.plot(input[components[0]], input[components[2]], lw=0.75)
    axY.grid(True)
    axY.set_ylabel(r"$y^s_1$", fontsize=12)

    axZ = plt.subplot2grid(shape=(3, 2), loc=(2, 1), rowspan=1, colspan=1)
    axZ.plot(input[components[0]], input[components[3]], lw=0.75)
    axZ.grid(True)
    axZ.set_ylabel(r"$z^s$", fontsize=12)
    axZ.set_xlabel('Time', fontsize=12)

    plt.setp(axX.get_xticklabels(), visible=False)
    plt.setp(axY.get_xticklabels(), visible=False)
    #plt.setp(axZ.get_xticklabels(), visible=False)

    # Plot phase space
    axXYZ = plt.subplot2grid(shape=(3, 2), loc=(0, 0), rowspan=3, colspan=1, projection='3d')
    axXYZ.plot(input[components[2]], input[components[1]], input[components[3]], lw=0.75)
    axXYZ.set_xlabel(r"$y^s_1$", fontsize=12)
    axXYZ.set_ylabel(r"$x_1$", fontsize=12)
    axXYZ.set_xlim(1.2, 0)
    axXYZ.set_zlabel(r"$z^s$", fontsize=12)
    # axXYZ.axis('off')

    fig.tight_layout()

    # plt.axis('off')
    fig = plt.gcf()
    plt.draw()
    plt.savefig('Simulations/'+filepath+'.pdf', dpi=500)
    plt.show()
    #plt.savefig('Rossler_no_background.pdf', dpi = 500)
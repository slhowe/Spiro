#!/bin/bash

# Add my extensions to path
import sys
sys.path.insert(0, '/home/sarah/Documents/Spirometry/python/extensions')

# Import my extensions
from data_struct import DataStore
from breath_analysis import BreathData, split_breaths, calc_flow_delay
from calculus import integral, derivative
from filters import hamming
from resistance_plots import ResistanceTable

#Import built-ins
import matplotlib.pyplot as plt
from math import sin, cos, pi, exp, log
from numpy import isnan, array, real, nan
from numpy.linalg import lstsq

''' Using data files:
    find E = (resistance*flow)/(volume)

    Plots are saved if plot_path given
    Data in kPa and L and s'''

# Definitions
INSP = 0
EXP = 1
sampling_frequency = 125

# Data path
path = '/home/sarah/Documents/Spirometry/data/'
# Place to save plots
plot_path = '/home/sarah/Documents/Spirometry/lungPlots/'

#files = ['Normal_1.csv']
files = ['Loops_1.csv']
files = [path + name for name in files]
filename = files[0]

plot=True

# Remove any old data sets hanging around
# (Needed when using ipython)
reset = DataStore()
reset.clear_dataset_registry()

# Create data storage
# Store pressure and flow data from files
dataset = DataStore(Fs=sampling_frequency)
dataset.store_data(filename)

# Filter data
# Use hamming filter
# data, Fc, Fs, bw
flow = hamming(dataset.flow, 50, 125, 20, plot=False)
flow = real(flow).tolist()
pressure = [-p for p in dataset.pressure]

# Split flow into breaths
# Find start, middle and end of each breath
flow_splits = split_breaths(dataset.flow)
flow_starts = flow_splits[0]
flow_middles = flow_splits[1]
flow_stops = flow_splits[2]

# Set up storage for breath data
# Corrects for pressure delay
delay = calc_flow_delay(pressure, flow)
data = BreathData(delay)

# Storage for full results
RC_all = [0]*(2*len(flow_starts))
vol_all = [0]*len(flow)

# Print starting info
print('Delay: {}'.format(delay))
print('Num breaths: {}'.format(len(flow_starts)))

# Iterate through breaths
print('\nRemember:\nRC should be NEGATIVE\n*** V = PC - RCQ ***\nand we find RC assuming positive')
for breath in range(1):#len(flow_starts)):
    # Extract breath data
    data.get_data(pressure,
                  flow,
                  flow_starts[breath],
                  flow_middles[breath],
                  flow_stops[breath],
                  dataset.sampling_frequency
                  )

    # Save arrays of full volume for breath
    vol_all[flow_starts[breath]:flow_middles[breath]] = data.insp_volume
    vol_all[flow_middles[breath]:flow_stops[breath]] = data.exp_volume

    # Go through each breath one half at a time
    # First: inspiration
    # Second: expiration
    for section in range(1):

        # Set data for insp/exp here
        if(section == INSP):
            print('\nInsp parameters:')
            breath_length = data.insp_length
            vol = data.insp_volume
            flw = data.insp_flow
            pres = data.insp_pressure
            time = data.time[:breath_length]
        elif(section == EXP):
            print('\nExp parameters:')
            breath_length = data.exp_length
            vol = data.exp_volume
            flw = data.exp_flow
            pres = data.exp_pressure
            time = data.time[-breath_length:]


        #---------------------------------------------#
        # Modelling driving pressure as an elastance  #
        # Trying elastance shape as exponential decay #
        #---------------------------------------------#

        # Take a guess at the windpipe resistance
        Rwp = 0.5

        # Get the spirometer resistance
        dependent = array([pres])
        independent = array([flw])
        res = lstsq(independent.T, dependent.T)
        print('Rsp: {}'.format(res[0][0][0]))

        # Get the shape of elastance
        shape = [nan]*breath_length
        for i in range(breath_length):
            if vol[i] == 0:
                vol[i] = 1e-8
            shape[i] = flw[i] / vol[i]

        crop = breath_length/10
        pres_crop = pres[crop:]
        flw_crop = flw[crop:]
        vol_crop = vol[crop:]
        shape_crop = shape[crop:]
        breath_length_crop = breath_length - crop

        # Setup arrays for least squares
        ones = [1]*breath_length

        # Make the left hand side of the equation
        Pin = [pres_crop[i] + Rwp*flw_crop[i] for i in range(breath_length_crop)]

        # Least squares
        dependent = array([Pin])
        independent = array([shape_crop])#, flw_crop])
        res = lstsq(independent.T, dependent.T)

        magnitude= res[0][0][0]
        #Raw = exp(res[0][1][0])
        #print('Raw: {}'.format(Raw))
        print('Magnitude: {}'.format(magnitude))
        print('Residuals:{}'.format(res[1]))


        # Remake data used to get parameters
        E = [magnitude * shape[i]
            for i in range(breath_length)]
        Pin_Sim = [E[i]*vol[i]# + (Raw * flw[i])
            for i in range(breath_length)]

        if(plot):
            f, (ax2, ax3, ax4, ax5) = plt.subplots(4, sharex=True)
            ax2.set_title('Flow')

            # Plot flow
            ax2.plot(flw)
            ax2.set_ylabel('Flow (L/s)')
            ax2.grid()

            # Plot volume
            ax3.plot(vol)
            ax3.set_ylabel('Volume (L)')
            ax3.grid()

            # Plot elastance and resistance over breath
            ax4.plot(E)
            #ax4.plot([Raw]*breath_length)
            ax4.legend(['Elastance', 'Resistance'])
            ax4.grid()

            # Plot real and simulated data
            ax5.plot(Pin)
            ax5.plot(Pin_Sim)
            ax5.set_ylabel('Pressure - Rwp*Q')
            ax5.legend(['Orig', 'Sim'])
            ax5.grid()

            plt.show()

            # Plot saved if plot path given
            if(plot_path != '\0'):
                fig_name = filename + '_iteration_{}'.format(iteration)
                print(fig_name)
                plt.savefig(fig_name)

            plt.close()


# plot parameters for each breath
# The plot the breaths
# Plot the volumes underneath
f, (ax1, ax2, ax3) = plt.subplots(3, sharex=True)
ax1.plot(flow_middles, RC_all[1::2], 'oc')
ax1.plot(flow_middles, RC_all[0::2], 'xr')
ax1.set_title('RC values')
ax1.legend(['exh', 'inh'])
ax1.grid()
ax2.plot(flow)
ax2.plot(flow_starts, [0]*len(flow_starts), 'og')
ax2.plot(flow_stops, [0]*len(flow_stops), 'or')
ax2.legend(['Data', 'Start point', 'End point'])
ax2.grid()
ax3.plot(vol_all)
ax3.plot([0]*len(vol_all), 'r')
ax3.set_ylabel('volume')
ax3.grid()
plt.show()

# Final clean up
dataset.clear_dataset_registry()



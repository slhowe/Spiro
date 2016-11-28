#!/bin/bash

# Add my extensions to path
import sys
sys.path.insert(0, '/home/sarah/Documents/Spirometry/python/extensions')

# Import my extensions
from data_struct import DataStore
from breath_analysis import BreathData, split_breaths, calc_flow_delay
from calculus import integral, derivative
from filters import hamming

#Import built-ins
import matplotlib.pyplot as plt
from math import sin, cos, pi, exp, log
from numpy import isnan, array, real, nan
from scipy import io
from numpy.linalg import lstsq

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# function definitions

def inflection_points(data, Fs, plot=False):
    """
    Inflection points are at zero crossings
    of second derivative.

    Returns indices of crossings, with Indices at the
    point before the crossing occurs

    Note: last 2 dp are missing from second derivative.
    Any inflections in the last points will be missed.
    """
    # Find second derivative of data
    der = derivative(data, Fs)
    derder = derivative(der, Fs)

    def signum(num):
        """
        Returns sign of a number
        -1 if negative
        1 if positive
        0 if zero
        """
        if(num < 0): return -1;
        elif(num > 0): return 1;
        else: return 0;

    # Find crossings
    crossings = [0]*len(derder)
    num_crossings = 0
    i = 0

    while(i < len(derder) - 1):
        # If the sign of the number changes

        if(signum(derder[i]) != signum(derder[i+1])):
            # Record current index
            # Increment number of crossings found
            crossings[num_crossings] = i
            num_crossings += 1
            i += 2
        # increment counter
        i += 1

    # Chop off unused indices
    if(crossings[0] < 4):
        crossings = crossings[1:]

    final_crossings = crossings[:num_crossings + 1]
    final_crossings[-2] = len(derder) - 1
    # Last inflection at end of data
    final_crossings[-1] = len(data) - 1

#    # Check for peaks
#    # Don't care about small noise, only big changes
#    MIN_PEAK = 1
#    final_crossings = [0]*(num_crossings + 1)
#    index = 0
#
#    # For every crossing found
#    for i in range(num_crossings):
#        peak_not_found = True
#        # Look for large peak in between this and next crossing
#        for j in range(crossings[i], crossings[i+1]):
#            # If found a peak, record crossing
#            if(abs(derder[j]) >= MIN_PEAK and peak_not_found):
#                if(index == 0):
#                    if(crossings[i] > 4):
#                        final_crossings[index] = crossings[i]
#                        peak_not_found = False
#                        index += 1
#                elif(j - final_crossings[index-1] > 4):
#                    final_crossings[index] = crossings[i]
#                    peak_not_found = False
#                    index += 1
#
#    index += 1
#
#    # Chop off unused indices
#    final_crossings = final_crossings[0:index]

    if(plot):
        # Plot data, derivative and second derivative
        # Data plot shows inflection points
        l, (axa, axb, axc) = plt.subplots(3, sharex=True)
        axa.plot(data, 'd-')
        axa.plot(crossings, [data[c] for c in crossings], 'rd')
        axa.plot(final_crossings, [data[c] for c in final_crossings], 'yd')
        axa.set_ylabel("Data")
        axa.grid()

        axb.plot(der, 'd-')
        axb.set_ylabel("Derivative")
        axb.grid()

        axc.plot(derder, 'd-')
        axc.set_ylabel("Second Derivative")
        axc.grid()

        plt.show()

    return final_crossings

def dirty_model_pressure(start, end, flow):
    pressure_estimation = [0]*(end-start)
    pressure_offset = 0

    # Find inflection points
    inflections = inflection_points(flow, Fs=50, plot=False)

    # Get index of the last inflection
    last_inflection = -1
    for inflection in inflections:
        if(inflection < end):
            last_inflection += 1

    start_point = start
    for index in inflections[:last_inflection]:
        if(index <= end):
            flow_section = flow[start_point:index]

            # Integrate if flow increasing
            if(flow[index] > flow[start_point]):
                pressure_section = integral(flow_section, 50)
                pressure_section = [p + pressure_offset for p in pressure_section]
            # Hold constant if flow decreasing
            else:
                pressure_section = [pressure_offset]*len(flow_section)

            # Update pressure estimate and pressure offset
            pressure_estimation[start_point-start:index-start] = pressure_section
            pressure_offset = pressure_section[-1]
            start_point = index

    # Carry over last points if not at inflection point
    if(inflections[last_inflection] < end):
        flow_section = flow[index:end]

        # Integrate if flow increasing
        if(flow[end] > flow[index]):
            pressure_section = integral(flow_section, 50)
            pressure_section = [p + pressure_offset for p in pressure_section]
        # Hold constant if flow decreasing
        else:
            pressure_section = [pressure_offset]*len(flow_section)

        # Update pressure estimate and pressure offset
        pressure_estimation[index-start:end] = pressure_section

    P_max = max(pressure_estimation)
    pressure_estimation = [p/P_max for p in pressure_estimation]

    return pressure_estimation

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Definitions
INSP = 0
EXP = 1
sampling_frequency = 125

# Data path
path = '/home/sarah/Documents/Spirometry/data/'

# Place to save plots
plot_path = '/home/sarah/Documents/Spirometry/lungPlots/'

files = ['Normal_1.csv']
#files = ['Loops_1.csv']
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
all_flow = hamming(dataset.flow, 50, 125, 20, plot=False)
all_flow = real(all_flow).tolist()

# Split flow into breaths
# Find start, middle and end of each breath
flow_splits = split_breaths(dataset.flow, plot=False)
flow_starts = flow_splits[0]
flow_middles = flow_splits[1]
flow_stops = flow_splits[2]

# Set up storage for breath data
data = BreathData(0)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Iterate through breaths
first_breath = 0
last_breath = len(flow_starts)

Ea = [nan]*last_breath
Ra = [nan]*last_breath
Es = [nan]*last_breath
Rs = [nan]*last_breath
for breath in range(first_breath,last_breath): # good 84-134

    print('\nBreath number {}\n~~~~~~~~~~~~~~~~'.format(breath))

    # Extract breath data
    data.get_data(dataset.pressure,
                  all_flow,
                  flow_starts[breath],
                  flow_middles[breath],
                  flow_stops[breath],
                  dataset.sampling_frequency
                  )

    pressure = data.insp_pressure + data.exp_pressure
    flow = data.insp_flow + data.exp_flow
    volume = data.insp_volume + data.exp_volume
    start_insp = 0
    end_insp = data.insp_length

    print('start_insp: {}'.format(start_insp))
    print('end_insp: {}'.format(end_insp))

    # Crop data to insp range
    flw = flow[start_insp:end_insp]
    pres = pressure[start_insp:end_insp]
    vol = volume[start_insp:end_insp]

    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    # Model the pressure
    pressure_estimation = dirty_model_pressure(start_insp, end_insp, flow)
    #exp_pressure_estimation = dirty_model_pressure(end_insp, len(flow), flow)

    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    # Get Eop and Rop from estimation
    dependent = array([pressure_estimation])
    independent = array([flw, vol])
    res = lstsq(independent.T, dependent.T)
    EoP = res[0][1][0]
    RoP = res[0][0][0]


    pressure_estimation_orig = pressure_estimation
    Q_orig = [(pressure_estimation[i] - vol[i]*EoP)/RoP for i in range(len(vol))]

    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    MAX_ITERATIONS = 195
    MAX_FITTING_ERROR = 0.4
    MAX_PERCENT_ERROR = 0.9

    iteration = 0
    flow_fitting_error = 1
    while(iteration < MAX_ITERATIONS and flow_fitting_error > MAX_FITTING_ERROR):

        #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~``

        # Forward simulate flow from pressure estimation
        Q = [(pressure_estimation[i] - vol[i]*EoP)/RoP for i in range(len(vol))]
        Q_error = [flw[i] - Q[i] for i in range(len(flw))]
        flow_fitting_error = [q**2 for q in Q_error]
        flow_fitting_error = sum(flow_fitting_error)
       # print('error: {}'.format(flow_fitting_error))

        max_flow = max(flow)
        max_flow += max_flow*(iteration/10.0)

        P_error = [0]*len(pressure_estimation)

        for i in range(len(Q)):
            # Work out the error
            if(flow[i] != 0):
                percent_error = (Q_error[i])/max_flow
                if(percent_error > MAX_PERCENT_ERROR):
                    percent_error = MAX_PERCENT_ERROR
                elif (percent_error < -MAX_PERCENT_ERROR):
                    percent_error = -MAX_PERCENT_ERROR
            else:
                percent_error = 0

            P_error[i] = (1 + percent_error)*pressure_estimation[i]

        pressure_estimation= P_error
        Q2 = [(P_error[i] - vol[i]*EoP)/RoP for i in range(len(vol))]

        iteration += 1

        #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~``

    # Get Eop and Rop after iteration
    dependent = array([pressure_estimation])
    independent = array([flw, vol])
    res = lstsq(independent.T, dependent.T)
    E_s = res[0][1][0]
    R_s = res[0][0][0]

    Ea[breath] = (EoP)
    Ra[breath] = (RoP)
    Es[breath] = (E_s)
    Rs[breath] = (R_s)

    print('')
    print('Eop: {}'.format(EoP))
    print('Rop: {}'.format(RoP))
    print('E_s: {}'.format(E_s))
    print('R_s: {}'.format(R_s))

    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # plot stuff
    print(iteration)
    if(0):
        f, (ax1, ax2, ax3) = plt.subplots(3, sharex=True)
        ax1.plot(range(start_insp,end_insp), pressure_estimation_orig, 'g*-', linewidth=3)
        ax1.plot(range(start_insp,end_insp), pressure_estimation, 'r*-', linewidth=3)

        ax2.plot(flow[0:end_insp], 'mx-')
        ax2.plot(range(start_insp,end_insp), flw, '^-', color='#ffddf4', linewidth=3)
        ax2.plot(range(start_insp,end_insp), Q_orig, 'g*-')
        ax2.plot(range(start_insp,end_insp), Q_error, 'r*-')
        ax2.plot(range(start_insp,end_insp), Q2, 'c*-')

        ax3.plot(volume[0:end_insp],'yx-')

        ax1.grid()
        ax2.grid()
        ax3.grid()
        plt.show()

if(1):
    plt.plot(all_flow)
    plt.grid()

    f, (ax1, ax2) = plt.subplots(2, sharex=True)
    ax1.plot(Ea, 'or')
    ax1.plot(Es, '^b')
    ax1.legend(['Before iteration', 'After iteration'])
    ax1.set_title('Elastance/scaling factor')
    ax1.grid()

    ax2.plot(Ra, 'or')
    ax2.plot(Rs, '^b')
    ax2.legend(['Before iteration', 'After iteration'])
    ax2.set_title('Resistance/scaling factor')
    ax2.grid()

    plt.show()

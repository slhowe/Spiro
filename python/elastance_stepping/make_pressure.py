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

# Definitions
INSP = 0
EXP = 1
sampling_frequency = 125

# Data path
path = '/home/sarah/Documents/Spirometry/data/ventilation/'

# Place to save plots
plot_path = '/home/sarah/Documents/Spirometry/lungPlots/'

files = ['ManualDetection_Patient4_PM.mat']
files = [path + name for name in files]
filename = files[0]

plot=True

# Load matlab file
mat = io.loadmat(filename)
full_data = mat['ManualDetection']

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
        # increment counter
        i += 1

    # Chop off unused indices
    crossings = crossings[:num_crossings]

    if(plot):
        # Plot data, derivative and second derivative
        # Second derivative plot also has inflection points
        l, (axa, axb, axc) = plt.subplots(3, sharex=True)
        axa.plot(data, 'd-')
        axa.set_ylabel("Data")
        axa.grid()

        axb.plot(der, 'd-')
        axb.set_ylabel("Derivative")
        axb.grid()

        axc.plot(derder, 'd-')
        axc.plot(crossings, [derder[c] for c in crossings], 'rd')
        axc.set_ylabel("Second Derivative")
        axc.grid()

        plt.show()

    return crossings

# Iterate through breaths
for breath in range(89,90): # good 84-134

    print('\n {}'.format(breath))
    # Extract breath data
    pressure = full_data['Pressure'][0][breath]
    flow = full_data['Flow'][0][breath]

    # pressure and flow in weird format
    # change into array
    pressure = pressure.tolist()
    pressure = [p[0] for p in pressure] # every item is own list
    flow = flow.tolist()
    flow = [f[0] for f in flow]

    # filter
    flow = hamming(flow, 20, 50, 10)
    flow = real(flow).tolist()
    print(type(flow))
    volume = integral(flow, 50)

    # only care about inspiration at the moment
    start_insp = 0
    end_insp = 0
    max_flow = 0
    i = 0
    while i < (len(flow)-5):
        if(flow[i] > max_flow):
            start_insp = i + 2
            max_flow = flow[i]
        if(flow[i]<= 0 and pressure[i+10]<(pressure[i]-5)):
            end_insp = i - 5
            i = len(flow)
        i += 1

    end_insp = start_insp
    start_insp = 0
    print('start_insp: {}'.format(start_insp))
    print('end_insp: {}'.format(end_insp))

    # Get peep
    peep_data = pressure[-30:-20]
    peep = sum(peep_data)/len(peep_data)
    print('peep: {}'.format(peep))

    # Crop data to insp range
    flw = flow[start_insp:end_insp]
    pres = pressure[start_insp:end_insp]
    pres = [p - peep for p in pres]
    vol = volume[start_insp:end_insp]

    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # Find R and E directly from pressure
    dependent = array([pres])
    independent = array([flw, vol])
    res = lstsq(independent.T, dependent.T)

    E = res[0][1][0]
    R = res[0][0][0]
    print('E: {}'.format(E))
    print('R: {}'.format(R))

    pressure_drop = [p - peep for p in pressure]
    remade_pres = [E*volume[i] + R*flow[i] for i in range(len(flow))]
    vol_sized_pres = [p/pres[-1]*vol[-1] for p in remade_pres]
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
    # Integral method for finding pressure from flow #
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
    # Estimate the pressure from flow
    # Try:
    #   integral(flow) = k*pressure + p0
    # over small ranges between inflection points
    # fitting linear and exponential curves between
    # sections and taking best fit. Note: exponential
    # flow fit gives results from constant pressure.

    inflections = inflection_points(flow, Fs=50, plot=True)
    PR = [q*1 for q in flw]
    int_PR = integral(PR, 50)

    # Fitting normalised estimate to known flow data
    # To see how it compares to direct calculation
    ones = [1]*len(flw)
    dependent = array([pres])
    independent = array([int_PR])#, ones])
    res = lstsq(independent.T, dependent.T)
    #offset = res[0][1][0]
    magnitude = res[0][0][0]

    #estimate = [(magnitude*ipr + offset) for ipr in int_PR]
    estimate = [(magnitude*ipr) for ipr in int_PR]
    estimate2 = [ipr for ipr in int_PR]

        # Find R and E
    dependent = array([estimate])
    independent = array([flw, vol])
    res = lstsq(independent.T, dependent.T)

    E_est = res[0][1][0]
    R_est = res[0][0][0]
    print('E_est: {}'.format(E_est))
    print('R_est: {}'.format(R_est))
    print('E/R_est: {}'.format(E_est/R_est))
    remade_pres_est = [E_est*volume[i] + R_est*flow[i] for i in range(len(flow))]
    flow_est = derivative(int_PR, 50)
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    # plot stuff
    f, (ax1, ax2, ax3) = plt.subplots(3, sharex=True)
    ax1.plot(pressure_drop, 'bx-')
    ax1.plot(range(start_insp,end_insp), pres, 'cx-')
    ax1.plot(range(start_insp,end_insp), estimate, 'o-', color='#b0e0e6')
    ax1.plot(remade_pres, 'kx-')
    ax1.plot(remade_pres_est, 'k-')
    ax2.plot(flow, 'mx-')
    ax2.plot(range(start_insp,end_insp-1), flow_est, 'r*-')
    ax2.plot(range(start_insp,end_insp), PR, 'ro-')
    ax2.plot(range(start_insp,end_insp), flw, '^-', color='#ffddf4')
    ax3.plot(volume,'yx-')
    ax3.plot(range(start_insp,end_insp),
            [p+vol_sized_pres[start_insp] for p in int_PR],'x-', color='#978248')
    ax3.plot(vol_sized_pres,'x-', color = '#f08080')
    ax1.grid()
    ax2.grid()
    ax3.grid()
    plt.show()

#    # Make a pressure curve
#    # Here a filtered square wave
#    pres = [(0.004*t)+1 for t in range(data.insp_length)]
#    offset = 6
#    pres[0:offset] = [0]*offset
#    pres[-offset:] = [0]*offset
#    pres[-2*offset:-offset] = [0.25]*offset
#    pres = hamming(pres, 2, 125, 50, plot=False)
#
#    # Make elastance shape
#    # Set any 0 volume measurement to be very small
#    for i in range(data.insp_length):
#        if(data.insp_volume[i] <= 0):
#            data.insp_volume[i] = 1e-5
#
#    # E = P/V
#    elas = [pres[i]/data.insp_volume[i]
#            for i in range(data.insp_length)]
#    elas = [e*0.001 for e in elas]
#
#    # Plot it
#    plt.plot(data.insp_flow)
#    plt.plot(data.insp_volume)
#    plt.plot(pres)
#    plt.plot(elas)
#    plt.show()
#
#    # Just try it out
#    dependent = array([pres])
#    independent = array([data.insp_volume, data.insp_flow])
#    res = lstsq(independent.T, dependent.T)
#
#    E = res[0][0][0]
#    R = res[0][1][0]
#    print('E = {}'.format(E))
#    print('R = {}'.format(R))
#    print('Residuals = {}'.format(res[1][0]))
#
#    # Remake pressure
#    sim_pres = [e*data.insp_volume[i] + R*data.insp_flow[i]
#                for i in range(data.insp_length)]
#
#    # Make varying elastance
#    Edrs = [(pres[i] - R*data.insp_flow[i])/data.insp_volume[i]
#            for i in range(data.insp_length)]
#
#    # Remake pressure
#    sim_pres_dyn = [Edrs[i]*data.insp_volume[i] + R*data.insp_flow[i]
#                for i in range(data.insp_length)]
#
#    plt.plot(pres)
#    plt.plot(sim_pres)
#    plt.plot(Edrs)
#    plt.plot(sim_pres_dyn)
#    plt.ylim([-5,5])
#    plt.grid()
#    plt.show()


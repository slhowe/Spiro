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
        sign_derder = signum(derder[i])
        # Record current index
        if(sign_derder != signum(derder[i+1])):
            # When going from negative to positive
            # The inflection is one point forward
            if(sign_derder == -1):
                crossings[num_crossings] = i - 1
            else:
                crossings[num_crossings] = i
            # Increment number of crossings found
            num_crossings += 1
        # increment counter
        i += 1

    # Chop off unused indices
    crossings = crossings[:num_crossings]

    # Check for peaks
    # Don't care about small noise, only big changes
    MIN_PEAK = 50
    final_crossings = [0]*num_crossings
    index = 0

    # For every crossing found
    for i in range(num_crossings - 1):
        peak_not_found = True
        # Look for large peak in between this and next crossing
        for j in range(crossings[i], crossings[i+1]):
            # If found a peak, record crossing
            if(abs(derder[j]) >= MIN_PEAK and peak_not_found):
                final_crossings[index] = crossings[i]
                peak_not_found = False
                index += 1

    # Last inflection at end of data
    final_crossings[index] = len(data) - 1
    index += 1

    # Chop off unused indices
    final_crossings = final_crossings[0:index]

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

# Iterate through breaths
for breath in range(89,90): # good 84-134

    print('\nBreath number {}'.format(breath))
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
    volume = integral(flow, 50)

    # only care about inspiration at the moment
    start_insp = 0
    end_insp = 0
    start_not_found = True
    i = 0
    while i < (len(flow)):
        # Find first positive index for start of insp
        if(flow[i] > 0 and start_not_found):
            start_insp = i
            start_not_found = False

        # Find negative flow for end of insp
        if(flow[i]<= 0):
            end_insp = i - 1
            i = len(flow)

        i += 1

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
    #   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # over small ranges between inflection points
    # fitting linear and exponential curves between
    # sections and taking best fit. Note: exponential
    # flow fit gives results from constant pressure.

    Pressure_estimation = [0]*len(pres)

    # Find inflection points
    inflections = inflection_points(flow, Fs=50, plot=True)

    # Go through flow data
    # Decide on shape between inflection points and integrate for the shape
    start_point = 0
    for index in inflections:
        # Only looking at inspiration at the moment
        if(index <= end_insp):
            flow_section = flw[start_point:index]

            # Test whether exponential or linear gives better fit to data
            time = range(len(flow_section))
            ones = [1]*len(time)

            #   Linear:
            dependent = array([flow_section])
            independent = array([time, ones])
            linear_res = lstsq(independent.T, dependent.T)

            #   Exponential:
            ln_flow = [log(f) for f in flow_section]
            dependent = array([ln_flow])
            independent = array([time, ones])
            exp_res = lstsq(independent.T, dependent.T)

            # Compare fits
            linear_error = linear_res[0][1]
            exp_error = exp_res[0][1]
            print('Linear_error: {}'.format(linear_error))
            print('exp_error: {}'.format(exp_error))

            # If the linear fit is best:
            # Integrate flow over the section with the offset
            # as the previous value. If no previous value start
            # at zero pressure offset
            if(abs(linear_error) > abs(exp_error)):
                print('--> exponential fit wins')
                section_pressure = integral(section_flow, 50)


            # If the exponential fit is best:
            # Approximate with a constant pressure
            # Hold pressure constant at previous offset
            # Hope there is no exponential from zero offset
            else:
                print('--> linear fit wins')

            # Integrate data between inflection points
            PR = [q*1 for q in flw]
            int_PR = integral(PR, 50)

            start_point = index

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

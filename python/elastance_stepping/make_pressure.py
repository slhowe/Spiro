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

    # Check for peaks
    # Don't care about small noise, only big changes
    MIN_PEAK = 40
    final_crossings = [0]*num_crossings
    index = 0

    # For every crossing found
    for i in range(num_crossings - 1):
        peak_not_found = True
        # Look for large peak in between this and next crossing
        for j in range(crossings[i], crossings[i+1]):
            # If found a peak, record crossing
            if(abs(derder[j]) >= MIN_PEAK and peak_not_found):
                if(index == 0):
                    final_crossings[index] = crossings[i]
                    peak_not_found = False
                    index += 1
                elif(j - final_crossings[index-1] > 4):
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


def model_pressure(pres, flow):
    """
    Integral method for finding pressure from flow
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Estimate the pressure from flow
    Try:
      integral(flow) = k*pressure + p0
      ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    over small ranges between inflection points
    fitting linear and exponential curves between
    sections and taking best fit. Note: exponential
    flow fit gives results from constant pressure.
    """

    pressure_estimation = [0]*(len(pres))
    pressure_offset = 0

    # Find inflection points
    inflections = inflection_points(flow, Fs=50, plot=True)

    # Go through flow data
    # Decide on shape between inflection points and integrate for the shape
    start_point = 0
    for index in inflections:
        # Only looking at inspiration at the moment
        if(index <= len(pres)):
            flow_section = flow[start_point:index]

            # Test whether exponential or linear gives better fit to data
            time = [t/50.0 for t in range(len(flow_section))]
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
            linear_error = abs(linear_res[1][0])
            exp_error = abs(exp_res[1][0])

           # Find the errors as a ratio of each other
            # Use this to define how much emphasis to put
            # on the fit
            total_error  = linear_error + exp_error
            linear_percent_error = 1 - (linear_error / total_error)
            exp_percent_error = 1 - (exp_error / total_error)
            #print('\nLinear_error: {} = {}%'.format(linear_error, linear_percent_error*100))
            #print('exp_error: {} = {}%'.format(exp_error, exp_percent_error*100))

            if(0):
                lin_line = [(t*linear_res[0][0][0] + linear_res[0][1][0]) for t in time]
                exp_line = [(exp(exp_res[0][0][0]*t) * exp(exp_res[0][1][0])) for t in time]
                plt.plot(flow_section, 'd')
                plt.plot(lin_line, 'd')
                plt.plot(exp_line, 'd')
                plt.grid()
                plt.show()

            if(exp_percent_error < 0.2):
                # Linear fit:
                # Integrate flow over the section with the offset
                # as the previous value. If no previous value start
                # at zero pressure offset
                pressure_section = integral(flow_section, 50)
                pressure_section = [p + pressure_offset for p in pressure_section]

            else:
                # If the exponential fit is best:
                # Approximate with a constant pressure
                # Hold pressure constant at previous offset
                # Hope there is no exponential from zero offset
                pressure_section = [pressure_offset]*len(flow_section)

            # Update pressure estimate and pressure offset
            pressure_estimation[start_point:index] = pressure_section
                    #                exp_percent_error*exp_pressure_section[i]
                    #              + linear_percent_error*linear_pressure_section[i]
                    #                for i in range(len(time))]
            pressure_offset = pressure_estimation[index - 1]

            start_point = index

    return pressure_estimation

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Iterate through breaths
for breath in range(89,92): # good 84-134

    print('\nBreath number {}\n~~~~~~~~~~~~~~~~'.format(breath))
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
    peep = min(sum(peep_data)/len(peep_data), pressure[0])
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

    # Model the pressure
    pressure_estimation = model_pressure(pres, flow)

    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # Make the estimation big to see how it overlaps
    max_P = max(pressure_estimation)
    pressure_estimation = [p/max_P for p in pressure_estimation]

    # Scale the pressure up
    dependent = array([pres])
    independent = array([pressure_estimation])
    res = lstsq(independent.T, dependent.T)
    pressure_estimation_scaled = [p*res[0][0][0] for p in pressure_estimation]

    # Use scaled data to directly predict
    dependent = array([pressure_estimation_scaled])
    independent = array([flw, vol])
    res = lstsq(independent.T, dependent.T)

    E_s = res[0][1][0]
    R_s = res[0][0][0]
    print('E scaled: {}'.format(E_s))
    print('R scaled: {}'.format(R_s))

    # Look at how good the prediction is for E/R
    dependent = array([pressure_estimation])
    independent = array([flw, vol])
    res = lstsq(independent.T, dependent.T)

    EoP = res[0][1][0]
    RoP = res[0][0][0]
    print('E/R actual: {}'.format(E/R))
    print('E/R estimate: {}'.format(EoP/RoP))

    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # Find flow back from estimation
    Q = [(pressure_estimation[i] - vol[i]*EoP)/RoP for i in range(len(vol))]
    Q_error = [flw[i] - Q[i] for i in range(len(flw))]


    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # plot stuff
    f, (ax1, ax2, ax3) = plt.subplots(3, sharex=True)
    ax1.plot(pressure_drop, 'bd-')
    ax1.plot(range(start_insp,end_insp), pres, 'cx-')
    ax1.plot(remade_pres, 'kx-')
    #ax1.plot(pressure_estimation, 'rx-')
    ax1.plot(pressure_estimation_scaled, 'md-')

    ax2.plot(flow, 'mx-')
    ax2.plot(range(start_insp,end_insp), flw, '^-', color='#ffddf4')
    ax2.plot(Q, 'g*-')
    ax2.plot(Q_error, 'r*-')

    ax3.plot(volume,'yx-')
    ax3.plot(vol_sized_pres,'x-', color = '#f08080')

    ax1.grid()
    ax2.grid()
    ax3.grid()
    plt.show()

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

files = ['ManualDetection_Patient8_PM.mat']
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
    i = 3

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
    crossings = crossings[:num_crossings]

    # Check for peaks
    # Don't care about small noise, only big changes
    MIN_PEAK = 1
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
                    if(crossings[i] > 4):
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

def dirty_model_pressure(start, end, flow):
    pressure_estimation = [0]*(end-start)
    pressure_offset = 0

    # Find inflection points
    inflections = inflection_points(flow, Fs=50, plot=True)

    # Go through flow data
    # Decide on shape between inflection points and integrate for the shape
    start_point = start
    for index in inflections:
        # Only looking at inspiration at the moment
        if(index <= end):
            flow_section = flow[start_point:index]

            if(flow[index] > flow[start_point]):
                pressure_section = integral(flow_section, 50)
                pressure_section = [p + pressure_offset for p in pressure_section]
            else:
                pressure_section = [pressure_offset]*len(flow_section)

            # Update pressure estimate and pressure offset
            pressure_estimation[start_point-start:index-start] = pressure_section
            pressure_offset = pressure_section[-1]
            start_point = index

    pressure_estimation[-1] = pressure_estimation[-2]
    return pressure_estimation

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Iterate through breaths
for breath in range(80,99): # good 84-134

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
            # Skip forwards. Don't want end close to start
            # if there is a small drop straight away.
            # 10 dp is about 200 ms so not too long
            i += 10

        # Find negative flow for end of insp
        if(flow[i]<= 0 and not start_not_found):
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

    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    # Find R and E directly from pressure
    dependent = array([pres])
    independent = array([flw, vol])
    res = lstsq(independent.T, dependent.T)

    E = res[0][1][0]
    R = res[0][0][0]
    print('E: {}'.format(E))
    print('R: {}'.format(R))
    print('E/R actual: {}'.format(E/R))
    print('')

    pressure_drop = [p - peep for p in pressure]
    remade_pres = [E*volume[i] + R*flow[i] for i in range(len(flow))]
    vol_sized_pres = [p/pres[-1]*vol[-1] for p in remade_pres]

    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    # Model the pressure
    pressure_estimation = dirty_model_pressure(start_insp, end_insp, flow)

    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    # Scale the pressure up
    dependent = array([pres])
    independent = array([pressure_estimation])
    res = lstsq(independent.T, dependent.T)
    scaling = res[0][0][0]
    pressure_estimation_scaled = [p*scaling for p in pressure_estimation]

    # Use scaled data to directly predict
    dependent = array([pressure_estimation_scaled])
    independent = array([flw, vol])
    res = lstsq(independent.T, dependent.T)

    E_s = res[0][1][0]
    R_s = res[0][0][0]
    print('E scaled: {}'.format(E_s))
    print('R scaled: {}'.format(R_s))
    print('E/R estimate: {}'.format(E_s/R_s))

    # Get Eop and Rop from estimation
    dependent = array([pressure_estimation])
    independent = array([flw, vol])
    res = lstsq(independent.T, dependent.T)
    EoP = res[0][1][0]
    RoP = res[0][0][0]

    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    for iteration in range(3):
        print('')

        #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~``

        # Forward simulate flow from pressure estimation
        Q = [(pressure_estimation[i] - vol[i]*EoP)/RoP for i in range(len(vol))]
        Q_error = [flw[i] - Q[i] for i in range(len(flw))]
        max_flow = max(flow)

        #P_error = model_pressure(end_insp, Q_error)
        #P_error = [(vol[i]*EoP) + (flw[i]*RoP) for i in range(len(vol))]
        P_error = [0]*len(pressure_estimation)
        aa = []

        for i in range(len(Q)):
            # Work out the error
            if(flow[i] != 0):
                percent_error = (Q_error[i])/max_flow
            else:
                percent_error = 0

            aa.append(percent_error)
            P_error[i] = (1 + percent_error)*pressure_estimation[i]

        #new_P = [pressure_estimation[i] + P_error[i] for i in range(len(flw))]
        #new_P_scaled = [p*scaling for p in new_P]
        P_error_scaled = [p*scaling for p in P_error]

        pressure_estimation= P_error
        Q2 = [(P_error[i] - vol[i]*EoP)/RoP for i in range(len(vol))]

        #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~``

        # Use scaled data to directly predict
        dependent = array([pres])
        independent = array([pressure_estimation])
        res = lstsq(independent.T, dependent.T)
        scaling = res[0][0][0]
        pressure_estimation_scaled = [p*scaling for p in pressure_estimation]

        dependent = array([pressure_estimation_scaled])
        independent = array([flw, vol])
        res = lstsq(independent.T, dependent.T)

        E_s = res[0][1][0]
        R_s = res[0][0][0]
        print('E scaled: {}'.format(E_s))
        print('R scaled: {}'.format(R_s))
        print('E/R estimate: {}'.format(E_s/R_s))

        #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        # plot stuff
        if(1):
            f, (ax1, ax2, ax3) = plt.subplots(3, sharex=True)
            ax1.plot(pressure_drop[0:90], 'bd-', linewidth=3)
            #ax1.plot(range(start_insp,end_insp), pres, 'cx-')
            #ax1.plot(remade_pres[0:90], 'kx-')
            #ax1.plot(pressure_estimation, 'rx-')
            ax1.plot(pressure_estimation_scaled, 'md-')
            #ax1.plot(new_P_scaled, 'g*-')
            ax1.plot(P_error_scaled, 'r*-', linewidth=3)
            ax1.plot(aa)

            ax2.plot(flow[0:90], 'mx-')
            ax2.plot(range(start_insp,end_insp), flw, '^-', color='#ffddf4')
            ax2.plot(Q, 'g*-')
            ax2.plot(Q_error, 'r*-')
            ax2.plot(Q2, 'c*-')

            ax3.plot(volume[0:90],'yx-')
            ax3.plot(vol_sized_pres[0:90],'x-', color = '#f08080')

            ax1.grid()
            ax2.grid()
            ax3.grid()
            plt.show()

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

def peak_points(data, Fs, plot=False):
    """
    Local minima and maxima are at zero crossings
    of derivative. Function returns indices of crossings

    Note: last 1 data points are missing from derivative.
    A crossing in the last point will be missed.
    """

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

    # Find second derivative of data
    der = derivative(data, Fs)

    # Find crossings
    crossings = [0]*len(der)
    num_crossings = 0
    i = 3

    while(i < len(der)):
        # If the sign of the number changes

        if(signum(flow[i-1]) != signum(flow[i])):
            crossings[num_crossings] = i
            num_crossings += 1
            i += 1
        elif(signum(der[i-1]) != signum(der[i])):
            # Record current index
            # Increment number of crossings found
            crossings[num_crossings] = i
            num_crossings += 1
            i += 1
        # increment counter
        i += 1

    # Chop off unused indices
    final_crossings = crossings[:num_crossings + 1]

    # Last point at end of data
    final_crossings[-1] = len(data) - 1

    if(plot):
        # Plot data, derivative and second derivative
        # Data plot shows point points
        l, (axa, axb) = plt.subplots(2, sharex=True)
        axa.plot(data, 'd-')
        axa.plot(crossings, [data[c] for c in crossings], 'rd')
        axa.plot(final_crossings, [data[c] for c in final_crossings], 'yd')
        axa.set_ylabel("Data")
        axa.grid()

        axb.plot(der, 'd-')
        axb.set_ylabel("Derivative")
        axb.grid()

        plt.show()

    return final_crossings


def dirty_model_pressure(start, end, flow, volume):

    def estimate_pressure(flow, volume, start, end, pressure_offset, factor):
        flow_section = flow[start:end]
        vol_section = integral(flow_section, 50)
        flow_gradient = (flow_section[-1] - flow_section[0])/(len(flow_section)/50.0)
        vol_gradient = (vol_section[-1] - vol_section[0])/(len(vol_section)/50.0)

        factor = 1 - (factor - 1)/2
        if(factor < 0):
            factor = 0

        # Want to flip the gradient if data is negative
        if(flow_section[0] < 0):
            flow_gradient = -flow_gradient
            vol_gradient = -vol_gradient

        # Flow increasing
        if(flow_gradient > 0):
            pressure_section = integral(flow_section, 50)
            pressure_section = [pressure_offset + p*flow_gradient
                                for p in pressure_section]

        # flow decreasing
        else:
            pressure_section = integral(flow_section, 50)
            pressure_section = [pressure_offset + p*factor
                                for p in pressure_section]

        return pressure_section

    # Setup
    pressure_estimation = [0]*(end-start)
    pressure_offset = 0

    Q_max = max(flow)
    V_max = max(volume)

    factor = (Q_max)/(V_max)
    print('factor: {}\n'.format(factor))

    # Find points
    # Pressure is estimated in small sections
    # between these points. Points are at minima,
    # maxima, and zero crossings in flow
    points = peak_points(flow, Fs=50, plot=True)

    # Get index of the last point
    last_point = 0
    for point in points:
        if(point < end):
            last_point += 1

    # Set starting points
    start_point = start

    # Looking at data in between points
    for index in points[:last_point]:
        if(index > start):

            pressure_section = estimate_pressure(flow,
                                                 volume,
                                                 start_point,
                                                 index,
                                                 pressure_offset,
                                                 factor
                                                 )
            # Update pressure estimate and pressure offset
            pressure_estimation[start_point-start:index-start] = pressure_section
            pressure_offset = pressure_section[-1]
            start_point = index

    # Carry over last points
    start_point = points[last_point - 1]
    if(start_point < end):

        pressure_section = estimate_pressure(flow,
                                             volume,
                                             start_point,
                                             end,
                                             pressure_offset,
                                             factor
                                             )
        # Update pressure estimate
        pressure_estimation[start_point-start:] = pressure_section

    return pressure_estimation

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Definitions
INSP = 0
EXP = 1
sampling_frequency = 50

# Place to save plots
plot_path = '/home/sarah/Documents/Spirometry/lungPlots/'

# Specify matlab files to use
files = ['ManualDetection_Patient4_PM.mat']
#files = ['ManualDetection_Patient14_FM.mat']
#files = ['ManualDetection_Patient8_PM.mat']
#files = ['ManualDetection_Patient17_PM.mat']

# Make the data path
path = '/home/sarah/Documents/Spirometry/data/ventilation/'
files = [path + name for name in files]
filename = files[0]

# Load matlab file
mat = io.loadmat(filename)
full_data = mat['ManualDetection']

# Specify breaths to iterate through
first_breath = 0
last_breath = 466

# Space to save results
scaling_factors = [nan]*last_breath
Ea = [nan]*last_breath
Ra = [nan]*last_breath
Es = [nan]*last_breath
Rs = [nan]*last_breath

# Iterate through every breath in range specified
for breath in range(first_breath,last_breath):
    print('\nBreath number {}\n~~~~~~~~~~~~~~~~'.format(breath))

    # Extract pressure and flow data
    pressure = full_data['Pressure'][0][breath]
    flow = full_data['Flow'][0][breath]

    # Pressure and flow are in a weird format
    # of array of lists of single values.
    # Convert array to list and pull values out
    # into main list.
    pressure = pressure.tolist()
    pressure = [p[0] for p in pressure] # every item is own list
    flow = flow.tolist()
    flow = [f[0] for f in flow]

    # filter data
    flow = hamming(flow, 20, 50, 10)
    flow = real(flow).tolist()

    # Get the volume
    volume = integral(flow, 50)

    # Get peep
    peep_data = pressure[-len(pressure)/3:]
    peep = sum(peep_data)/len(peep_data)
    print('peep: {}'.format(peep))

    # Remove peep from pressure
    # Offset for all pressure data is now 0
    # Offset of estimated pressure will be 0
    pressure = [p - peep for p in pressure]

    # Start of inspiration:
    # This is at first crossing from negative flow to
    # positive flow, or at the first index. Stop looking
    # after max flow value. That point is definitely in
    # inspiration. Work backwards from max flow, because
    # start data can be noisy around zero crossing.
    start_insp = 0
    Q_max = max(flow)
    Q_max_index = flow.index(Q_max)
    i = Q_max_index
    while(i > 0):
        if(flow[i] >= 0 and flow[i-1] < 0):
            start_insp = i
            i = 0
        i -= 1

    # End of inspiration:
    # This is at the point where flow goes negative
    # followed by a large drop in pressure
    # If that's not a thing, take the last index
    end_insp = len(flow) - 1
    i = start_insp
    while i < (len(flow)-15):
        # Find negative flow for end of insp
        if(flow[i] < 0 and (pressure[i+15] < peep+4)):
            end_insp = i - 1
            i = len(flow)
        i += 1

    print('start_insp: {}'.format(start_insp))
    print('end_insp: {}'.format(end_insp))

    # Check there are more than the minimum data points
    # needed for least squares in inspiration and that
    # flow starts close to zero.
    # If either condition isn't met, skip the dataset; It sucks.
    if(end_insp - start_insp <= 3
    or flow[start_insp] > 0.1):
        print('Bad data, ignoring')

    else:
        # Crop data to insp range
        flw = flow[start_insp:end_insp]
        pres = pressure[start_insp:end_insp]
        vol = volume[start_insp:end_insp]

        #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        ### THIS SECTION USES PRESSURE DATA ###

        # Find R and E directly from pressure and flow
        # and volume in the single compartment model
        dependent = array([pres])
        independent = array([flw, vol])
        res = lstsq(independent.T, dependent.T)

        E = res[0][1][0]
        R = res[0][0][0]

        print('E: {}'.format(E))
        print('R: {}'.format(R))
        print('R/E actual: {}'.format(R/E))
        print('')

        # Remake pressure and flow from parameters
        remade_pres = [E*volume[i] + R*flow[i] for i in range(len(flow))]
        remade_flow = [(pressure[i] - E*volume[i])/R for i in range(len(flow))]

        #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

        # Model the pressure
        pressure_estimation = dirty_model_pressure(start_insp, end_insp, flow, volume)

        # Get parameters from estimated pressure
        dependent = array([pressure_estimation])
        independent = array([flw, vol])
        res = lstsq(independent.T, dependent.T)
        E_est = res[0][1][0]
        R_est = res[0][0][0]

        print('E_est: {}'.format(E_est))
        print('R_est: {}'.format(R_est))
        print('R_est/E_est: {}'.format(R_est/E_est))
        print('')

        # Using the estimated pressure, E should = 1.
        # If not, there was an error in the magnitude of estimate.
        # Divide pressure estimate by E_est to correct for error.
        pressure_estimation_updated = [p/E_est for p in pressure_estimation]

        # Update changes to parameters
        E_est /= E_est
        R_est /= E_est

        print('E_updated: {}'.format(E_est))
        print('R_updated: {}'.format(R_est))
        print('R_updated/E_updated: {}'.format(R_est/E_est))
        print('')

        # Forward simulate flow from pressure estimate and parameters
        Q_orig = [(pressure_estimation_updated[i] - E_est*vol[i])/R_est
                 for i in range(len(vol))]

        #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        ### THIS SECTION USES PRESSURE DATA ###

        # Scale the pressure up
        scaling = E
        pressure_estimation_scaled_orig = [p * scaling
                                           for p in pressure_estimation]
        pressure_estimation_scaled_updated = [p * scaling
                                              for p in pressure_estimation_updated]

        # Use scaling to directly predict parameters
        E_s = E_est * E
        R_s = R_est * E
        print('E scaled: {}'.format(E_s))
        print('R scaled: {}'.format(R_s))
        print('R/E estimate: {}'.format(R_s/E_s))

        #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

        MAX_ITERATIONS = 35
        MAX_FITTING_ERROR = 0.4

        iteration = 0
        flow_fitting_error = 1
        while(iteration < MAX_ITERATIONS and flow_fitting_error > MAX_FITTING_ERROR):
            #print('')

            # Forward simulate flow from pressure estimation
            Q = [(pressure_estimation[i] - vol[i]*EoP)/RoP for i in range(len(vol))]
            Q_error = [flw[i] - Q[i] for i in range(len(flw))]
            flow_fitting_error = [q**2 for q in Q_error]
            flow_fitting_error = sum(flow_fitting_error)

            max_flow = max(flow)
            max_flow += max_flow*(iteration/10.0)

            P_error = [0]*len(pressure_estimation)

            for i in range(len(Q)):
                # Work out the error
                if(flow[i] != 0):
                    percent_error = (Q_error[i])/max_flow
                    if(percent_error > 1):
                        percent_error = 0.99
                    elif (percent_error < -1):
                        percent_error = -0.99
                else:
                    percent_error = 0

                P_error[i] = (1 + percent_error)*pressure_estimation[i]

            P_error_scaled = [p*scaling for p in P_error]

            pressure_estimation= P_error
            Q2 = [(P_error[i] - vol[i]*EoP)/RoP for i in range(len(vol))]

            #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~``
            #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~``
            #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~``
            #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~``

            # Use scaled data to directly predict
            dependent = array([pres])
            independent = array([pressure_estimation])
            res = lstsq(independent.T, dependent.T)
            scaling = res[0][0][0]
            scaling = E
            pressure_estimation_scaled = [p*scaling for p in pressure_estimation]

            dependent = array([pressure_estimation_scaled])
            independent = array([flw, vol])
            res = lstsq(independent.T, dependent.T)

            E_s = res[0][1][0]
            R_s = res[0][0][0]
            #print('E scaled: {}'.format(E_s))
            #print('R scaled: {}'.format(R_s))
            #print('E/R estimate: {}'.format(E_s/R_s))

            iteration += 1

        # Get Eop and Rop from estimation
        dependent = array([pressure_estimation])
        independent = array([flw, vol])
        res = lstsq(independent.T, dependent.T)
        EoP = res[0][1][0]
        RoP = res[0][0][0]
        print('EoP: {}'.format(EoP))
        print('RoP: {}'.format(RoP))
        print('R/E : {}'.format(RoP/EoP))

        scaling_factors[breath] = (scaling)
        Ea[breath] = (E)
        Ra[breath] = (R)
        Es[breath] = (E_s)
        Rs[breath] = (R_s)

        comb_flow = [(pressure_estimation_scaled[i] - vol[i]*E)/R for i in range(len(vol))]
        comb_pressure = [R*flw[i] +  vol[i]*E for i in range(len(vol))]

        #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        # plot stuff
        if(1):
            f, (ax1, ax2, ax3) = plt.subplots(3, sharex=True)

            ax1.plot(pressure[0:end_insp], 'b-', linewidth=3)
            ax1.plot(range(start_insp,end_insp), remade_pres, 'b-')
            ax1.plot(range(start_insp,end_insp), pressure_estimation_scaled_orig, 'm:')
            ax1.plot(range(start_insp,end_insp), pressure_estimation_scaled_updated, 'm-')
            ax1.plot(range(start_insp,end_insp), comb_pressure, 'b^-')

            ax2.plot(flow[0:end_insp], 'm-')
            ax2.plot(range(start_insp,end_insp), flw, '^-', color='#ffddf4', linewidth=3)
            ax2.plot(range(start_insp,end_insp), remade_flow, 'b-')
            ax2.plot(range(start_insp,end_insp), Q_orig, 'c*-')

            ax3.plot(volume[0:end_insp],'yx-')

            ax1.grid()
            ax2.grid()
            ax3.grid()
            plt.show()

if(1):
    f, (ax1, ax2, ax3) = plt.subplots(3, sharex=True)
    ax1.plot(Ea, 'or')
    ax1.plot(Es, '^b')
    ax1.legend(['Direct', 'Modelled and scaled'])
    ax1.set_title('Elastance')
    ax1.grid()

    ax2.plot(Ra, 'or')
    ax2.plot(Rs, '^b')
    ax2.legend(['Direct', 'Modelled and scaled'])
    ax2.set_title('Resistance')
    ax2.grid()

    ax3.plot(scaling_factors, 'dg')
    ax3.set_title('Scaling factor')
    ax3.grid()

    plt.show()

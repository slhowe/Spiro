#!/bin/bash

# Add my extensions to path
import sys
sys.path.insert(0, '/home/sarah/Documents/Spirometry/python/extensions')

# Import my extensions
from data_struct import DataStore
from breath_analysis import BreathData, split_breaths, calc_flow_delay
from calculus import integral, derivative
from filters import hamming
import data_extraction as extr

#Import built-ins
import matplotlib.pyplot as plt
from math import sin, cos, pi, exp, log
from numpy import isnan, array, real, nan
from scipy import io
from numpy.linalg import lstsq

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

#######################################
# Indices to split data into sections #
#######################################
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

    # Find derivative of data
    der = derivative(data, Fs)

    # Data is split into small chunks
    # Define minimum chunk size here
    MIN_SECTION_SIZE = 3

    # Setup
    crossings = [0]*len(der)
    num_crossings = 0
    i = MIN_SECTION_SIZE

    while(i < len(der)):
        # If the sign of the number changes

        # Indices are recorded at points where
        # flow crosses zero (inh <-> exh)
        # or at minima/maxima. Flow zero crossings
        # are priority.
        if(signum(flow[i-1]) != signum(flow[i])):
            crossings[num_crossings] = i
            num_crossings += 1
            i += (MIN_SECTION_SIZE - 2)
        elif(signum(der[i-1]) != signum(der[i])):
            crossings[num_crossings] = i
            num_crossings += 1
            i += (MIN_SECTION_SIZE - 2)

        # increment counter
        i += 1

    # Chop off unused indices
    final_crossings = crossings[:num_crossings + 1]

    # Last point at end of data
    final_crossings[-1] = len(data) - 1

    # Optional plotting
    if(plot):
        # Plot data, derivative
        # Data plot shows points
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

######################
# Pressure modelling #
######################
def model_pressure(start, end, flow, volume, pressure_offset):

    def estimate_pressure(flow, volume, start, end, pressure_offset, factor):
        # Looking at a small section of flow
        flow_section = flow[start:end]

        # Gradient of flow in section
        flow_gradient = (flow_section[-1] - flow_section[0])/(len(flow_section)/50.0)

       # Want to flip the gradient if flow is negative
        # This will subtract a slope from the pressure
        # offset so pressure decreases for negative flow
        if(flow_section[0] < 0):
            flow_gradient = -flow_gradient

        # Flow increasing
        # Models large jump at start of insp/exp
        # Jump size is very dependent on the rate of change of flow
        # Quicker rise time gives larger pressure jump
        if(flow_gradient > 0):
            print('flow_gradient: {}'.format(flow_gradient))
            pressure_section = integral(flow_section, 50)
            pressure_section = [pressure_offset + p*flow_gradient*12
                                for p in pressure_section]

        # Flow decreasing
        # Modelled as a constant pressure with a slight ramp
        # Slope of ramp depends on overall trend of breath
        # Trend found from (max flow/max volume) factor
        # If factor is large, pressure is quite flat
        # If factor is small, pressure has slope
        # Lower factor gives a steeper slope
        else:
            # Define factor used to control slope
            # of pressure when flow decreasing
            # 0.6 is a fudge-factor. Examples give
            # a range 0.5-0.7 so 0.6 was chosen
            new_factor = (1 - (factor - 1)/2) * 0.6
            if(new_factor < 0):
                new_factor = 0
            print('new factor: {}'.format(new_factor))

            pressure_section = integral(flow_section, 50)
            pressure_section = [pressure_offset + p*new_factor
                                for p in pressure_section]

        return pressure_section

    # Setup
    pressure_estimation = [0]*(end-start)

    Q_max = max(flow)
    V_max = max(volume)

    factor = (Q_max)/(V_max)
    print('\nfactor: {}\n'.format(factor))

    # Find points
    # Pressure is estimated in small sections
    # between these points. Points are at minima,
    # maxima, and zero crossings in flow
    points = peak_points(flow, Fs=50, plot=False)

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

    # Include any data after final point
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

#################
# Main function #
#################
if(__name__ == '__main__'):
    # Definitions
    INSP = 0
    EXP = 1
    sampling_frequency = 50

    # Place to save plots
    plot_path = '/home/sarah/Documents/Spirometry/lungPlots/'

    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    # There are different data files with different
    # data structures. This section extracts data
    # from different file types.
    using_ManualDetection_files = 1
    using_PS_vs_NAVA_invasive_files = 1

    files = []
    file_types = []

    # Declare file names for ManualDetection MV data
    if(using_ManualDetection_files):
        path = '/home/sarah/Documents/Spirometry/data/ventilation/ManualDetection/'
        MD_filenames = [
                        'ManualDetection_Patient4_PM.mat',
                        'ManualDetection_Patient14_FM.mat',
                        'ManualDetection_Patient17_FM.mat',
                       ]
        # Make full path names
        # Add names to list of all data
        # Add data typoe to data type list
        files += [path + name for name in MD_filenames]
        file_types += ['MD' for name in MD_filenames]

    # Declare file names for PS/Nava invasive ventilation data
    if(using_PS_vs_NAVA_invasive_files):
        path = '/home/sarah/Documents/Spirometry/data/ventilation/PS_vs_NAVA_invasive/'
        PNI_filenames = [
                         'BRU1-PS.mat',
                        ]
        # Make full path names
        # Add names to list of all data
        # Add data typoe to data type list
        files += [path + name for name in PNI_filenames]
        file_types += ['PNI' for name in PNI_filenames]

    # Go through every file declared
    file_index = 0
    for filename in files:
        # Load data from file
        mat_data = extr.load_mat_file(filename)

        # Determine which file type the data came from
        file_type = file_types[file_index]
        file_index += 1

        # Data extraction for ManualDetection type data,
        if(file_type == 'MD'):
            last_breath = 480
        # Data extraction for PS/NAVA invasive ventilation data,
        elif(file_type == 'PNI'):
            full_data = extr.PS_vs_NAVA_invasive_data(mat_data)
            PNI_pressure = full_data[0]
            PNI_flow = full_data[1]
            last_breath = len(PNI_flow)

        # Specify breaths to iterate through
        first_breath = 0

        # Make space to save results
        ER_actual = [nan]*last_breath
        ER_simulated = [nan]*last_breath

        # Iterate through every breath in range specified
        for breath in range(first_breath,last_breath):
            print('\nBreath number {}\n~~~~~~~~~~~~~~~~'.format(breath))

            # Data for each breath is extracted here depending on file type
            if(file_type == 'MD'):
                breath_data = extr.ManualDetection_data(mat_data, breath)
                pressure = breath_data[0]
                flow = breath_data[1]
            elif(file_type == 'PNI'):
                pressure = PNI_pressure[breath]
                flow = PNI_flow[breath]

        #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

            # filter data
            flow = hamming(flow, 20, sampling_frequency, 10)
            flow = real(flow).tolist()

            # Get the volume
            volume = integral(flow, sampling_frequency)

            # Get peep
            peep_data = pressure[-len(pressure)/3:]
            peep = sum(peep_data)/len(peep_data)
            print('peep: {}'.format(peep))

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

            # Remove peep from pressure
            # Offset for all pressure data is now 0
            # Offset of estimated pressure will be 0
            pressure = [p - peep for p in pressure]

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

                # Estimate the driving pressure
                # Guessing P = E * integral(Q) {approx}
                pressure_offset = 0
                pressure_estimation = model_pressure(start_insp,
                                                     end_insp,
                                                     flow[:end_insp],
                                                     volume[:end_insp],
                                                     pressure_offset)

                #pressure_offset = pressure_estimation[-1]
                #exp_pressure_estimation = model_pressure(end_insp,
                                                         #len(flow),
                                                         #flow[end_insp:],
                                                         #volume[:end_insp],
                                                         #pressure_offset)

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
                # This totally assumes the shape is correct
                pressure_estimation_updated = [p/E_est for p in pressure_estimation]

                # Update changes to parameters
                R_est /= E_est
                E_est /= E_est

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
                remade_pres = [E*vol[i] + R*flw[i] for i in range(len(flw))]
                remade_flow = [(pres[i] - E*vol[i])/R for i in range(len(flw))]

                # Scale the pressure up
                scaling = E
                pressure_estimation_scaled_orig = [p * scaling
                                                   for p in pressure_estimation]
                pressure_estimation_scaled_updated = [p * scaling
                                                      for p in pressure_estimation_updated]

                ER_actual[breath] = (R/E)
                ER_simulated[breath] = (R_est/E_est)

                # plot stuff
                if(1):
                    f, (ax1, ax2, ax3) = plt.subplots(3, sharex=True)

                    ax1.plot(pressure[0:end_insp], 'b-', linewidth=3)
                    ax1.plot(range(start_insp,end_insp), remade_pres, 'c-')
                    ax1.plot(range(start_insp,end_insp), pressure_estimation_scaled_orig, 'r:')
                    ax1.plot(range(start_insp,end_insp), pressure_estimation_scaled_updated, 'm-')
                    ax1.legend([
                                'Pressure',
                                'Forward sim from data',
                                'Original estimate (scaled)',
                                'Updated estimate (scaled)',
                                'Estimate after iteration (scaled)'
                                ], loc=4)

                    ax2.plot(flow[0:end_insp], 'r-', linewidth=3)
                    ax2.plot(range(start_insp,end_insp), remade_flow, 'b-')
                    ax2.plot(range(start_insp,end_insp), Q_orig, 'm*-')

                    ax2.legend([
                                'Flow',
                                'Forward sim from data',
                                'Forward sim from estimate',
                                'Forward sim after iteration',
                                'Error in flow'
                                ])

                    ax3.plot(volume[0:end_insp],'yx-')
                    ax3.legend([
                                'Volume'
                                ])

                    ax1.grid()
                    ax2.grid()
                    ax3.grid()
                    plt.show()

        if(1):
            f, (ax3) = plt.subplots(1, sharex=True)

            ax3.plot(ER_actual, 'or')
            ax3.plot(ER_simulated, '^b')
            ax3.legend(['Directly calculated from data', 'Estimated'])
            ax3.set_ylabel('Resistance/Elastance')
            ax3.set_xlabel('Breath')
            ax3.grid()

            plt.show()

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

##########
# SIGNUM #
##########
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

    def estimate_pressure(flow, volume, start, end, pressure_offset):
        # In typical breathing:
        #   - Breathing rate is low
        #   - Maximum flow rate is low
        # Because of low flow rate, pressure drop due to
        # airways is low (unless very restricted?).
        # Because breathing rate is low, lungs have time
        # to equalise to input pressure. This means the
        # lung pressure will closely follow input pressure.
        # So input pressure can be estimated as P = EV + P0.
        # Pressure estimation here is P^ = P/E so use
        # P^ = V + P^0 as our pressure estimate.

        # Volume is added to offset value so want first
        # value not equal to zero so the curve is smooth
        vol_section = [volume[i] - volume[start-1] for i in range(start,end)]

        # If flow is positive, pressure must be increasing.
        # If flow is negative, pressure must be decreasing.
        # Flow is split into sections between flow zero crossings,
        # minima and maxima
        if(flow[start] < 0):
            pressure_section = [pressure_offset - p
                                for p in vol_section]
        else:
            pressure_section = [pressure_offset + p
                                for p in vol_section]

        return pressure_section

    # Setup
    pressure_estimation = [nan]*(end-start)

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
    fac = []
    pnt = []
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
    using_dated_files = 1
    using_PS_vs_NAVA_invasive_files = 1
    using_PS_vs_NAVA_non_invasive_files = 1

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
        # Add data type to data type list
        files += [path + name for name in MD_filenames]
        file_types += ['MD' for name in MD_filenames]

    # Declare file names for dated data
    if(using_dated_files):
        path = '/home/sarah/Documents/Spirometry/data/ventilation/ManualDetection/'
        DF_filenames = [
                        '12_11_08.mat',
                        '13_11_21.mat',
                        '9_04_08.mat',
                        '9_04_09.mat',
                        '9_04_10A.mat',
                       ]
        # Make full path names
        # Add names to list of all data
        # Add data type to data type list
        files += [path + name for name in DF_filenames]
        file_types += ['DF' for name in DF_filenames]

   # Declare file names for PS/Nava invasive ventilation data
    if(using_PS_vs_NAVA_invasive_files):
        path = '/home/sarah/Documents/Spirometry/data/ventilation/PS_vs_NAVA_invasive/'
        PNI_filenames = [
                         'BRU1-PS.mat',
                         'BRU2-PS.mat',
                         'BRU4-PS.mat',
                         'BRU6-PS.mat',
                         'BRU14-PS.mat',
                         'GE04-PS.mat',
                         'GE05-PS.mat',
                         'GE11-PS.mat',
                         'GE21-PS.mat',
                        ]
        # Make full path names
        # Add names to list of all data
        # Add data type to data type list
        files += [path + name for name in PNI_filenames]
        file_types += ['PNI' for name in PNI_filenames]

    # Declare file names for PS/NAVA non-invasive ventilation data
    if(using_PS_vs_NAVA_non_invasive_files):
        path = '/home/sarah/Documents/Spirometry/data/ventilation/PS_vs_NAVA_non_invasive/'
        PNN_filenames = [
                        'NIV_BRU01.mat',
                        'NIV_BRU02.mat',
                        'NIV_BRU03.mat',
                        'NIV_BRU04.mat',
                        'NIV_BRU05.mat',
                        'NIV_BRU06.mat',
                        'NIV_BRU07.mat',
                        'NIV_BRU08.mat',
                        'NIV_BRU09.mat',
                        'NIV_BRU10.mat',
                        'NIV_LIE01.mat',
                        'NIV_LIE02.mat',
                        'NIV_LIE03.mat',
                       ]
        # Make full path names
        # Add names to list of all data
        # Add data type to data type list
        files += [path + name for name in PNN_filenames]
        file_types += ['PNN' for name in PNN_filenames]


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
        # Data extraction for dated ventilation data,
        elif(file_type == 'DF'):
            last_breath = 480
        # Data extraction for PS/NAVA invasive ventilation data,
        elif(file_type == 'PNI'):
            full_data = extr.PS_vs_NAVA_invasive_data(mat_data)
            full_pressure = full_data[0]
            full_flow = full_data[1]
            last_breath = len(full_flow)
        # Data extraction for PS/NAVA non-invasive ventilation data,
        elif(file_type == 'PNN'):
            full_data = extr.PS_vs_NAVA_noninvasive_data(mat_data)
            full_pressure = full_data[0]
            full_flow = full_data[1]
            last_breath = len(full_flow)

        # Specify breaths to iterate through
        first_breath = 0
        last_breath = 3

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
            elif(file_type == 'DF'):
                breath_data = extr.dated_data(mat_data, breath)
                pressure = breath_data[0]
                flow = breath_data[1]
            elif(file_type == 'PNI'):
                pressure = full_pressure[breath]
                flow = full_flow[breath]
            elif(file_type == 'PNN'):
                pressure = full_pressure[breath]
                flow = full_flow[breath]

            #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            # filter data
            flow = hamming(flow, 20, sampling_frequency, 10)
            flow = real(flow).tolist()
            pressure = hamming(pressure, 20, sampling_frequency, 10)
            pressure = real(pressure).tolist()

            # Get the volume
            volume = integral(flow, sampling_frequency)

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

            # Get peep
            #peep_data = pressure[-len(pressure)/3:]
            #peep = sum(peep_data)/len(peep_data)
            #peep -= 0.1
            peep_real = pressure[start_insp]
            print('peep real: {}'.format(peep_real))
            pressure_real = [p-peep_real for p in pressure]

            # End of inspiration:
            # Working backwards, find the last
            # point of positive flow in the data
            end_insp = len(flow)
            i = len(flow)*3/4
            while i > 0:
                if(flow[i] < 0 and flow[i-1] > 0):
                    end_insp = i - 1
                    i = 0
                i -= 1

           # Find first flow shoulder
            shoulder = Q_max
            shoulder_index = Q_max_index
            half_pres= pressure[start_insp:len(pressure)/8]
            if(len(half_pres) > 2):
                half_grad = (half_pres[-1] - half_pres[0])/(len(half_pres))
                flat_pres = [half_pres[i] - half_grad*i for i in range(len(half_pres))]
                shoulder = max(flat_pres)
                shoulder_index = flat_pres.index(shoulder)
                shoulder_index += start_insp + 2
            start_insp = shoulder_index
            print('shoulder: {}'.format(shoulder_index))

            start = start_insp + (end_insp - start_insp)*1/5
            end = end_insp - end_insp*1/18
            peep = pressure[start]

            print('start_insp: {}'.format(start_insp))
            print('end_insp: {}'.format(end_insp))
            print('start: {}'.format(start))
            print('end: {}'.format(end))

            # Remove peep from pressure
            # Offset for all pressure data is now 0
            # Offset of estimated pressure will be 0
            pressure = [p - peep for p in pressure]

            #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

            # Check there are more than the minimum data points
            # needed for least squares in inspiration and range
            # for estimating data
            if(end_insp - start_insp <= 3
            or end - start <= 3
            or isnan(flow[0])):
                print('Bad data, ignoring')

            else:
                # Crop data to insp range
                flw = flow[start:end]
                pres = pressure[start:end]
                vol = volume[start:end]

                # Estimate the driving pressure
                # Guessing P = E * integral(Q) {approx}
                pressure_offset = 0
                pressure_estimation = model_pressure(start,
                                                     end,
                                                     flow[:end],
                                                     volume[:end],
                                                     pressure_offset,
                                                     )

                #pressure_offset = pressure_estimation[-1]
                #exp_pressure_estimation = model_pressure(end_insp,
                                                         #len(flow),
                                                         #flow[end_insp:],
                                                         #volume[:end_insp],
                                                         #pressure_offset)

                # Get parameters from estimated pressure
                # E_est should = 1
                # R_est should = R/E
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
                # Divide estimates by E_est to correct for error.
                pressure_estimation_updated = [p/E_est for p in pressure_estimation]
                R_est /= E_est
                E_est /= E_est

                # Forward simulate flow from pressure estimate and parameters
                Q_orig = [(pressure_estimation_updated[i] - E_est*vol[i])/R_est
                         for i in range(len(vol))]
                P_est = [volume[i]*E_est + flow[i]*R_est
                        for i in range(len(flow)/2)]

                #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                # Now have an estimate for shape of pressure input over top
                # of inspiration curve. Don't have offset.
                # Assuming flat point in pressure at the end of inspiration
                # (goes from increasing/flat to decreasing at this point)
                # can take derivative to remove offset, with known value of
                # zero at the end.

                # Remake full pressure estimation of inspiration
                der_end = end - (end-start)/2
                der_offset = 3

                pressure_estimation_top = [E_est*volume[i] + R_est*flow[i]
                                         for i in range(start, end_insp + der_offset)]

                # Take the derivative
                pressure_estimation_der = derivative(pressure_estimation_top, 50)

                drop = 0#min(pressure_estimation_der)

                # Drop derivative so end at zero
                pressure_estimation_der = [d - drop
                            for d in pressure_estimation_der[:end-start-1]]

                # Find params and correct for error in V
                vol_der = (flow[start + 1:der_end])
                flow_der = derivative(flow[start:der_end], 50)
                dependent = array([pressure_estimation_der[:der_end - start - 1]])
                independent = array([flow_der, vol_der])
                res = lstsq(independent.T, dependent.T)
                E_est = res[0][1][0]
                R_est = res[0][0][0]

                pressure_estimation_der = [p/E_est for p in pressure_estimation_der]
                R_est /= E_est
                E_est = 1

                # Get parameters from estimated pressure
                dependent = array([pressure_estimation_der[:der_end - start - 1]])
                independent = array([flow_der, vol_der])
                res = lstsq(independent.T, dependent.T)
                E_est = res[0][1][0]
                R_est = res[0][0][0]

                print('E_est from estimated derivative: {}'.format(E_est))
                print('R_est from estimated derivative: {}'.format(R_est))
                print('R_est/E_est from estimated derivative: {}'.format(R_est/E_est))
                print('')

                #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                ### THIS SECTION USES PRESSURE DATA ###

                # Params from real pressure
                dependent = array([pressure_real[start_insp:end_insp]])
                independent = array([flow[start_insp:end_insp], volume[start_insp:end_insp]])
                res = lstsq(independent.T, dependent.T)
                Ea = res[0][1][0]
                Ra = res[0][0][0]

                print('E from actual pressure: {}'.format(Ea))
                print('R from actual pressure: {}'.format(Ra))
                print('R/E from actual pressure: {}'.format(Ra/Ea))
                print('')

                # From dropped pressure
                dependent = array([pressure[start_insp:end_insp]])
                independent = array([flow[start_insp:end_insp], volume[start_insp:end_insp]])
                res = lstsq(independent.T, dependent.T)

                Edrop = res[0][1][0]
                Rdrop = res[0][0][0]

                print('E from dropped pressure: {}'.format(Edrop))
                print('R from dropped pressure: {}'.format(Rdrop))
                print('R/E from dropped pressure: {}'.format(Rdrop/Edrop))
                print('')

                # Remake pressure and flow from parameters
                remade_pres = [Edrop*volume[i] + Rdrop*flow[i] for i in range(len(flow)/2)]
                remade_flow = [(pressure[i] - Edrop*volume[i])/Rdrop for i in range(len(flow)/2)]

                # Scale the pressure up
                scaling = 1/Edrop
                pressure_real_scaled = [p * scaling
                                   for p in pressure_real]
                pressure_scaled = [p * scaling
                                   for p in pressure]
                pressure_fwd_sim_scaled = [p * scaling
                                           for p in remade_pres]
                P_est_scaled = [p * scaling
                                for p in P_est]

                # Params from real scaled pressure derivative
                dependent = array([derivative(pressure_real_scaled[start:der_end], 50)])
                independent = array([flow_der, vol_der])
                res = lstsq(independent.T, dependent.T)

                Erd = res[0][1][0]
                Rrd = res[0][0][0]

                print('E from scaled pressure derivative: {}'.format(Erd))
                print('R from scaled pressure derivative: {}'.format(Rrd))
                print('R/E from scaled pressure derivative: {}'.format(Rrd/Erd))
                print('')

                # remake derivatives
                remade_est_der = [vol_der[i] + 0.35*flow_der[i]
                              for i in range(len(flow_der))]
                remade_der = [vol_der[i] + Rrd/Erd*flow_der[i]
                              for i in range(len(flow_der))]

                # Remake P/E
                remade_pres_act= [volume[i] + Rrd/Erd*flow[i]
                                      for i in range(start, end)]
                remade_pres_der= [volume[i] + R_est/E_est*flow[i]
                                      for i in range(start, end)]

                ER_actual[breath] = (Rdrop/Edrop)
                ER_simulated[breath] = (R_est/E_est)

                #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

                if(1):
                    f, (ax1, ax2, ax3) = plt.subplots(3, sharex=True)
                    plot_end = len(flow)/2

                    ax1.plot(range(start + 1,end), derivative(pressure_real_scaled[start:end], 50), 'b-', linewidth=3)
                    ax1.plot(range(start + 1,start+len(P_est[start:end])), derivative(P_est[start:end], 50), 'gs-', linewidth=3)
                    ax1.plot(range(start + 1,der_end), pressure_estimation_der[:der_end-start-1], 'cd-', linewidth=3)
                    ax1.plot(range(start + 1,der_end), remade_der, 'r-', linewidth=3)
                    ax1.plot(range(start + 1,der_end), remade_est_der, 'k*-', linewidth=3)
                    ax1.legend([
                                'Derivative pressure',
                                'Derivative pressure fwd sim',
                                'Derivative pressure estimate',
                                'Derivative from real R/E',
                                'Derivative from estimate R/E',
                                ], loc=4)

                    ax2.plot(range(start + 1,der_end), flow_der, 'ks-', linewidth=1)
                    ax2.legend([
                                'Derivative flow',
                                ])

                    ax3.plot(range(start + 1,der_end), vol_der, 'ys-', linewidth=3)
                    ax3.legend([
                                'Derivative Volume',
                                ])

                    ax1.grid()
                    ax2.grid()
                    ax3.grid()
                    plt.show()


                plt.plot([p/Ea for p in pressure_real[0:start]])
                plt.plot([(flow[i])*Ra/Ea*1.0 for i in range(start)])
                plt.show()

                # plot stuff
                if(1):
                    f, (ax1, ax2, ax3) = plt.subplots(3, sharex=True)
                    plot_end = len(flow)/2


                    ax1.plot(pressure_scaled[0:plot_end], 'b-', linewidth=3)
                    ax1.plot(range(plot_end), pressure_fwd_sim_scaled, 'k-', linewidth=2)
                    ax1.plot(range(plot_end), P_est, 'c-')
                    ax1.plot(range(start,end), pressure_estimation, 'r.')
                    ax1.plot(range(start,end), pressure_estimation_updated, 'm-', linewidth=3)
                    ax1.plot(range(start,end), remade_pres_act, 'r*', linewidth=3)
                    ax1.plot(range(start,end), remade_pres_der, 'm.', linewidth=3)
                    ax1.plot(start_insp, pressure_scaled[start_insp], 'go')
                    ax1.plot(end_insp, pressure_scaled[end_insp], 'ro')
                    ax1.legend([
                                'Pressure (scaled to P/E)',
                                'Forward sim from data (scaled to P/E)',
                                'Forward sim from estimate',
                                'Original estimate',
                                'Updated estimate',
                                'Forward sim from actual der',
                                'Forward sim from estimate der',
                                ], loc=4)

                    ax2.plot(flow[0:plot_end], 'r-', linewidth=3)
                    ax2.plot(range(plot_end), remade_flow, 'b-')
                    ax2.plot(range(start,end), Q_orig, 'm*-')
                    ax2.plot(start_insp, flow[start_insp], 'go')
                    ax2.plot(end_insp, flow[end_insp], 'ro')
                    ax2.legend([
                                'Flow',
                                'Forward sim from data',
                                'Forward sim from estimate',
                                ])

                    ax3.plot(volume[0:end_insp],'yx-')
                    ax3.legend([
                                'Volume',
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


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
from numpy import nan
import numpy as np
from scipy import io
from scipy.stats.stats import pearsonr
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
    derder = derivative(der, Fs)

    # Data is split into small chunks
    # Define minimum chunk size here
    MIN_SECTION_SIZE = 5

    # Setup
    crossings = [0]*len(der)
    num_crossings = 0
    i = MIN_SECTION_SIZE

    while(i < len(derder)):
        # If the sign of the number changes

        # Indices are recorded at points where
        # flow crosses zero (inh <-> exh)
        # or at minima/maxima. Flow zero crossings
        # are priority.
        if(signum(flow[i-1]) != signum(flow[i])):
            crossings[num_crossings] = i
            num_crossings += 1
            i += (MIN_SECTION_SIZE - 2)
        elif(signum(derder[i-1]) != signum(derder[i])):
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
def model_pressure(start, end, flow, volume, pressure_offset, points):

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

    # Get index of the last point
    last_point = 0
    for point in points:
        if(point < end):
            last_point += 1

    # Set starting points
    start_point = start

    # Looking at data in between points
    for index in points[:last_point]:
        if(index > start + 2):
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
    using_dated_files = 0
    using_PS_vs_NAVA_invasive_files = 0
    using_PS_vs_NAVA_non_invasive_files = 0

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
            sampling_frequency = 50.0
        # Data extraction for dated ventilation data,
        elif(file_type == 'DF'):
            last_breath = 480
            sampling_frequency = 50.0
        # Data extraction for PS/NAVA invasive ventilation data,
        elif(file_type == 'PNI'):
            full_data = extr.PS_vs_NAVA_invasive_data(mat_data)
            full_pressure = full_data[0]
            full_flow = full_data[1]
            last_breath = len(full_flow)
            sampling_frequency = 100.0
        # Data extraction for PS/NAVA non-invasive ventilation data,
        elif(file_type == 'PNN'):
            full_data = extr.PS_vs_NAVA_noninvasive_data(mat_data)
            full_pressure = full_data[0]
            full_flow = full_data[1]
            last_breath = len(full_flow)
            sampling_frequency = 100.0

        # Specify breaths to iterate through
        first_breath = 0
        last_breath = 3

        # Make space to save results
        ER_actual = []
        ER_simulated = []

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
            #flow = hamming(flow, 30, sampling_frequency, 10)
            #flow = np.real(flow).tolist()
            #pressure = hamming(pressure, 20, sampling_frequency, 10)
            #pressure = np.real(pressure).tolist()

            # Get the volume
            volume = integral(flow, sampling_frequency)

            # Start of inspiration:
            # This is at first crossing from negative flow to
            # positive flow, or at the first index. Stop looking
            # after max flow value. That point is definitely in
            # inspiration. Work backwards from max flow, because
            # start data can be noisy around zero crossing.
            start_insp = 0
            act_start_insp = 0
            Q_max = max(flow)
            Q_max_index = flow.index(Q_max)
            i = Q_max_index
            while(i > 0):
                if(flow[i] >= 0 and flow[i-1] < 0):
                    act_start_insp = i + 1
                    i = 0
                i -= 1

            # Get peep
            peep_data = pressure[-len(pressure)/3:]
            peep = sum(peep_data)/len(peep_data)

            # End of inspiration:
            # Working backwards, find the last
            # point of positive flow in the data
            end_insp = len(flow) - 1
            i = len(flow)*3/5
            while i > 0:
                if(flow[i] < 0 and flow[i-1] > 0):
                    end_insp = i - 1
                    i = 0
                i -= 1

           # Find first flow shoulder
            shoulder = Q_max
            shoulder_index = Q_max_index
            half_flow= flow[start_insp:len(pressure)/8]
            if(len(half_flow) > 2):
                half_grad = (half_flow[-1] - half_flow[0])/(len(half_flow))
                flat_flow = [half_flow[i] - half_grad*i for i in range(len(half_flow))]
                shoulder = max(flat_flow)
                shoulder_index = flat_flow.index(shoulder)
                shoulder_index += start_insp + 2
            start_insp = max(shoulder_index, Q_max_index)
            print('shoulder: {}'.format(shoulder_index))

            start = start_insp# + (end_insp - start_insp)*1/5
            end = end_insp - end_insp*1/3

            print('start_insp: {}'.format(start_insp))
            print('end_insp: {}'.format(end_insp))
            print('start: {}'.format(start))
            print('end: {}'.format(end))
            print('peep: {}'.format(peep))

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
            or np.isnan(flow[0])
            ):
                print('Bad data, ignoring')

            else:
                #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                ### THIS SECTION USES PRESSURE DATA ###

                # Params from real pressure
                dependent = np.array([pressure[start_insp:end_insp]])
                independent = np.array([flow[start_insp:end_insp], volume[start_insp:end_insp]])
                res = lstsq(independent.T, dependent.T)
                Ea = res[0][1][0]
                Ra = res[0][0][0]

                print('E from actual pressure: {}'.format(Ea))
                print('R from actual pressure: {}'.format(Ra))
                print('R/E from actual pressure: {}'.format(Ra/Ea))
                print('')

                # Remake pressure and flow from parameters
                remade_pres = [Ea*volume[i] + Ra*flow[i] for i in range(len(flow)/2)]
                remade_flow = [(pressure[i] - Ea*volume[i])/Ra for i in range(len(flow)/2)]

                # Scale the pressure
                scalinga = 1/Ea#max(pressure)
                pressure_scaled = [p * scalinga
                                   for p in pressure]
                pressure_fwd_sim_scaled = [p * scalinga
                                           for p in remade_pres]
                #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

                # This bit finds tau from pressure and flow data
                tau = Ra/Ea
                print('Tau: {}'.format(tau))
                V_max = max(volume)
                Q_max = max(flow)
                V_max_index = volume.index(V_max) - 3
                Q_max_index = flow.index(Q_max)
                V_est = [V_max*(1 - exp(-(i/sampling_frequency)/tau))
                        for i in range(V_max_index - act_start_insp)]
                Q_est = [Q_max*exp(-(i/sampling_frequency)/tau)
                        for i in range(end_insp - Q_max_index)]

                ln_flow = [0]*(V_max_index - act_start_insp)
                for i in range(V_max_index - act_start_insp):
                    if(flow[i + act_start_insp] > 0):
                        ln_flow[i] = [log(flow[i + act_start_insp])]
                    else:
                        ln_flow[i] = [log(1e-3)]

                # This bit guesses at dynamic tau
                #plt.plot(flow[:end_insp])
                #plt.plot([tau]*end_insp)
                time_varying_tau = [nan]*end_insp
                for i in range(end_insp):
                    if(volume[i] != 0):
                        if(flow[i] != 0):
                            time_varying_tau[i] = volume[i]/flow[i]
                    else:
                        time_varying_tau[0] = 1e-3

                #plt.plot(time_varying_tau)
                #plt.plot([np.mean(time_varying_tau)]*end_insp, ':')
                #plt.legend(['flow',
                 #           'tau',
                 #           'sec = 50',
                 #           'mn = 50',
                 #           ])
                #plt.show()

                # This bit guesses at tau from slope
                Q_max = max(flow[:end_insp])
                Q_max_index = flow.index(Q_max)
                t = 3/sampling_frequency
                tau_slp = -t/(log(flow[Q_max_index + 3]/flow[Q_max_index]))
                print(Q_max_index)
                print(tau_slp)

                # This bit does line stepping
                # sum of impulses signal thingy
                def fit_tau_to_flow(flow, end, tau_array, jump, altering_value):
                    # This bit does line stepping
                    # sum of impulses signal thingy
                    # using time varying tau
                    sp = 0
                    total_resp_tvt = [0] * end
                    total_jumps_tvt = [0] * end
                    add = 1
                    while(sp < end):
                        A = jump
                        time = range(end- sp)

                        resp = [A*(exp(-(time[i]/sampling_frequency)/tau_array[sp]))
                                for i in range(end- sp)]

                        for i in range(len(resp)):
                            total_resp_tvt[sp + i] += resp[i]

                        total_jumps_tvt[sp] += 1

                        i = 0
                        inc = False
                        add = 1
                        while(i < len(resp)):
                            if(total_resp_tvt[i+sp] < flow[i+sp]):
                                sp += i
                                inc = True
                                i = len(resp)
                            i += 1

                        if(not inc):
                            if(sp < end*0.5 and altering_value):
                                tau_array = [tau_array[0]*0.92]*end
                                print('New: {}'.format(0.92*tau_array[0]))
                                total_resp_tvt = [0] * end
                                total_jumps_tvt = [0] * end
                                sp = 0
                            else:
                                sp = end

                        #plt.plot(flow[:end_insp])
                        #plt.plot(total_resp_tvt)
                        #plt.show()

                    return(total_resp_tvt, total_jumps_tvt, tau_array[0])

                jump = 0.01
                res = fit_tau_to_flow(flow, end_insp, [tau]*end_insp, jump, False)
                total_resp = res[0]
                total_jumps = res[1]
                print(time_varying_tau)
                res = fit_tau_to_flow(flow, end_insp, [tau_slp]*end_insp, jump, True)
                total_resp_tvt = res[0]
                total_jumps_tvt = res[1]
                new_tau_slp= res[2]

                plt.plot(flow[:end_insp], 'b', linewidth=2)
                plt.plot(total_resp, 'g:')
                plt.plot(total_resp_tvt, 'r:')
                plt.plot(range(Q_max_index, end_insp),
                        [Q_max*exp(-t/sampling_frequency/tau)
                        for t in range(end_insp-Q_max_index)])
                plt.plot(range(Q_max_index, end_insp),
                        [Q_max*exp(-t/sampling_frequency/tau_slp)
                        for t in range(end_insp-Q_max_index)])
                plt.plot(range(Q_max_index, end_insp),
                        [Q_max*exp(-t/sampling_frequency/new_tau_slp)
                        for t in range(end_insp-Q_max_index)])
                plt.legend([
                            'flow',
                            'fit to real tau',
                            'fit to guessed tau',
                            'Decay from data',
                            'Decay estimate',
                            'Updated decay estimate',
                            ])
                plt.show()

                print(total_jumps)
                print(total_jumps_tvt)

                or i in range(len(total_jumps)):
                    total_jumps[i] *= jump*tau
                for i in range(len(total_jumps_tvt)):
                    total_jumps_tvt[i] *= jump*new_tau_slp#(time_varying_tau[0]/time_varying_tau[i])

                # Pressure guess
                # Guesses at P/R, I think...????
                # Is integral of all the jumps recorded scaled by the jump size
                guessy_pressure = [0]*len(total_jumps)
                for i in range(1, len(guessy_pressure)):
                    guessy_pressure[i] = guessy_pressure[i-1] + total_jumps[i-1]

                guessy_pressure_tvt = [0]*len(total_jumps)
                for i in range(1, len(guessy_pressure_tvt)):
                    guessy_pressure_tvt[i] = guessy_pressure_tvt[i-1] + total_jumps_tvt[i-1]
                guessy_pressure_scaled = [g/max(guessy_pressure_tvt)*max(guessy_pressure) for g in guessy_pressure_tvt]

                dependent = np.array([guessy_pressure])
                independent = np.array([volume[:end_insp], flow[:end_insp]])
                res = lstsq(independent.T, dependent.T)
                print(res)
                print(res[0][1][0]/res[0][0][0])

                if(1):
                    f, (ax1, ax2, ax3) = plt.subplots(3, sharex=True)
                    plot_end = len(flow)/2
                    ax1.set_title(file_type)

                    ax1.plot(pressure_scaled[:plot_end], 'b-', linewidth=3)
                    ax1.plot(range(plot_end), pressure_fwd_sim_scaled, 'k-', linewidth=2)
                    ax1.plot(guessy_pressure, 'g')
                    ax1.plot(guessy_pressure_tvt, 'r')
                    ax1.plot(guessy_pressure_scaled)
                    ax1.legend([
                                'Pressure (scaled to P/E)',
                                'Forward sim from data (scaled to P/E)',
                                'fit to real tau',
                                'fit to guessed tau',
                                ], loc=4)

                    ax2.plot(flow[0:plot_end], 'r-', linewidth=3)
                    ax2.plot(range(plot_end), remade_flow, 'b-')
                    ax2.plot(range(Q_max_index, end_insp), Q_est)
                    ax2.plot(range(Q_max_index, end_insp), [Q_max*exp(-t/sampling_frequency/tau) for t in range(end_insp-Q_max_index)])
                    ax2.plot(range(Q_max_index, end_insp), [Q_max*exp(-t/sampling_frequency/tau_slp) for t in range(end_insp-Q_max_index)])
                    ax2.legend([
                                'Flow',
                                'Flow, fwd sim',
                                'Flow from decay rate',
                                'Underlying decay rate',
                                'Guess at underlying decay rate',
                                ])

                    ax3.plot(volume[0:end_insp],'yx-')
                    ax3.plot(range(act_start_insp, V_max_index), V_est)
                    ax3.legend([
                                'Volume',
                                'Volume from decay rate',
                                ])

                    ax1.grid()
                    ax2.grid()
                    ax3.grid()
                    plt.show()

                #ER_actual.append(Ra/Ea)
                #ER_simulated.append(REs)

        #try:
            #corr = np.correlate(ER_actual, ER_simulated)
            #corr = pearsonr(ER_actual, ER_simulated)
            #print('Correlation: {}'.format(corr))
        #except:
            #pass

        #if(1):
            #f, (ax3) = plt.subplots(1, sharex=True)
            #ax3.set_title(file_type)

            #ax3.plot(ER_actual, 'or')
            #ax3.plot(ER_simulated, '^b')
            #ax3.legend(['Directly calculated from data', 'Estimated'])
            #ax3.set_ylabel('Resistance/Elastance')
            #ax3.set_xlabel('Breath')
            #ax3.grid()

            #plt.show()


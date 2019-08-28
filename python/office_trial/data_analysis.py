#!/bin/bash

# Import my extensions
import sys
sys.path.insert(0, '/home/sarah/Documents/Spirometry/python/extensions')

import csv
import matplotlib.pyplot as plt
import numpy as np
from data_struct import DataStore as Data
from calculus import integral
from breath_analysis import split_breaths
from numpy import array
from numpy.linalg import lstsq
from scipy import stats

path = './trial_data/'
subject = [
         # 'F01',
         # 'F02',
         # 'F03',
         # 'F04',
         # 'F05',
         # 'F09',
         # 'M01',
         # 'M02',
         # 'M03',
         # 'M04',
         # 'M05',
         # 'M06',
         # 'M07',
         # 'M08',
         # 'M09',
          'M10',
          ]
test_type = [
        #[
        #'_mp_peak.csv',
        #'_mp_peakRx.csv',
        #'_mp_peak.csv',
        #'_mp_peakRx.csv',
        # ],
        [
        '_mp_noRx.csv',
        '_mp_Rx.csv',
         ],
#        [
#        '_mp_noRx.csv',
#        '_mp_Rx.csv',
#         ],
        ]

# Coefficients for Ps -> Pm
# Pm = Ps*c[X][0] + c[X][1]
# X depends on pressure range [low, high]
# c = [crossing_point, [mp_RX, mp_woRX]]
conversion = [
        [-13, 1.3],
        [-7, 1.5]
        #[0.027, [1.1, 0], [2.9, 0.05]],
        #[0.023, [1.1, 0], [2.5, 0.03]],
        ]

sampling_frequency = 300


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Extract pressure two data arrays (assumed pressure and flow,
# but could be pressure & pressure etc) and time stamp from csv
# NOTE This function appends. Make sure empty arrays are passed in
def extract_data_arrays_from_file(filename, pressure, flow, time):
    # Store pressure and flow data
    with open(filename, 'rb') as csvfile:
        reader = csv.reader(csvfile, delimiter = ',')

        # skip header line
        header = reader.next()

        for row in reader:
            pressure_point = 1000 * float(row[0])
            flow_point = 1000 * float(row[1])
            time_point = float(row[2])

            pressure.append(pressure_point)
            flow.append(flow_point)
            time.append(time_point)

    return(1)
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Calibration for flow sensor is in
# ~/Documents/Spirometry/python/calibration
# This may need to be re-checked
# NOTE Only valid up to 250Pa across spirometer
def convert_spirometer_pressure_to_flow(pressure):
    flow = [p*0.0092 for p in pressure]
  #if (value >= 0) //inspiration
  #  flow = 0.0088*value;
  #else
  #  //flow = 0.0075*value;
  #  flow = 0.0092*value;
  #return flow;
    return flow
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Calibration for this conversion is in
# ~/Documents/Spirometry/python/maskLPF/Ps_Pm_calibration
# This conversion rate only applies to expiration
def convert_spirometer_pressure_to_mouth_pressure(spir_pressure, coeffs):
    mouth_pressure = [coeffs[0]*p**2 + coeffs[1]*p for p in spir_pressure]
#    mouth_pressure = [0]*len(spir_pressure)
#    for m in range(len(spir_pressure)):
#        if(abs(spir_pressure[m]) < coeffs[0]):
#            P = coeffs[1][0]*spir_pressure[m] + coeffs[1][1]
#        else:
#            P = coeffs[2][0]*spir_pressure[m] + coeffs[2][1]
#        mouth_pressure[m] = P

    return mouth_pressure
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
def root_mean_squared_error(data, fit, peak):
    mse = 0
    N = len(data)
    for m in range(N):
        mse += (data[m]/peak - fit[m]/peak)**2
    mse /= N
    rmse = np.sqrt(mse)
    return rmse

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#
def calculate_decay_rates(flow, starts, stops):
    PRINTING = False
    PLOTTING = False


    decay = [[np.nan, np.nan] for i in range(len(starts))]
    grad = [[np.nan, np.nan] for i in range(len(starts))]
    rvalues = [[np.nan] for i in range(len(starts))]
    decay_start = [0 for i in range(len(starts))]
    decay_end = [0 for i in range(len(starts))]
    peaks = [0 for i in range(len(starts))]

    # Fitting to range Ps = -50 <-> -70 Pa
    # Corresponds to Q = -20 <-> -27 L/min

    for breath in range(len(starts)):
        # Fit decay to each expiration
        expiration = flow[starts[breath]:stops[breath]]
        expiration_length = len(expiration)

        # Start from end and work way back up curve
        end_index = np.nan
        start_index = np.nan
        i = 0
        while i < expiration_length:
            # Working from end expiration back to start expiration
            # So the first points in range needed can be found
            index = expiration_length - i - 1
            if np.isnan(end_index):
                if expiration[index] <= -16:
                    end_index = index
            else:
                if expiration[index] <= -480:
                    start_index = index
                    i = expiration_length
            i += 1
        # Pick the minimum flow point as start if
        # expiration was too small to get start index
        if(PRINTING):
            print(expiration.index(min(expiration)))
            print(len(expiration))

        start_index = expiration.index(min(expiration)) + (len(expiration) - expiration.index(min(expiration)))/4
        end_index = len(expiration) - (len(expiration) - expiration.index(min(expiration)))/3

#        if np.isnan(start_index):
#            start_index = expiration.index(min(expiration)) + (len(expiration) - expiration.index(min(expiration)))/4
#            print('WARNING: Max flow less than start threshold (480 L/min)\n setting start to 1/3 range after flow minimum')
#        # Pick the end point as end of expiration
#        # if the breath was too small to get values
#        if np.isnan(end_index):
#            end_index = len(expiration) - (len(expiration) - expiration.index(min(expiration)))/3
#            print('WARNING: Max flow less than end threshold (16 L/min)\n setting end to 2/3 range after flow minimum')

        #halfway_point = (first_zero)/2
        if(PRINTING):
            print('start: {}'.format(start_index))
            print('end: {}'.format(end_index))

        decay_start[breath] = start_index
        decay_end[breath] = end_index

        # Grab the curve of interest and make it all positive
        # (assume no positive section in exp range)
        curve = expiration[start_index:end_index]
        curve = [np.abs(c) for c in curve]

        # Fit curve to exponential model
        # ln(Y) = At + ln(b)
        # for y(t) = e^(At) + b
        if(len(curve) > 3):
            Fs = 300.0
            measurements = array([np.log(curve)])

            # set up array of data multiplied by constants
            one_array = [1 for val in range(len(measurements.T))]
            times = [j/Fs for j in range(len(measurements.T))]
            independents = array([one_array,times])

            # Least squares fit to data
            result = lstsq(independents.T, measurements.T)

            # Squared Euclidean 2-norm for each col in (b - a*x)
            residual = result[1]
            if(PRINTING):
                print('       residual on line fit is {}'.format(residual[0]/len(curve)))

            # Parameters - offset and decay constant
            constants = result[0]
            constants[0] = np.exp(constants[0])
            if(PRINTING):
                print('       offset is {}, decay is {}'.format(constants[0], constants[1]))

            decay[breath] = constants

            # LINEAR FIT TO FLOW VS VOLUME :)
            volume = integral(expiration, sampling_frequency)
            volume_curve = volume[start_index:end_index]

            gradient, intercept, r_value, p_value, std_err = stats.linregress(volume_curve, expiration[start_index:end_index])


            #line2 = [gradient*s + intercept for s in steps]
            if(PRINTING):
                print('vvvvvvvvvvvv')
                print('       gradient is {}'.format(gradient))
                print('       offset is {}'.format(intercept))
                print('       R2 is {}'.format(r_value))
                print('       pValue is {}'.format(p_value))

                if abs(r_value) > 0.97:
                    print('***\nNICE AND LINEAR <3\n***')

            grad[breath] = [gradient, intercept]
            rvalues[breath] = [r_value]


            if(1):
                # find peak flow and peak flow index
                peak = min(expiration)
                peak_index = expiration.index(peak) #index relative end insp (flow middle)
                peaks[breath] = peak_index
                if(PRINTING):
                    print("peak:{}".format(peak))
                    print("peak_index:{}".format(peak_index))


                # Remake curve from best fit
                # Length of curve is from decay_start to flow_stop point
                Fs = 300.0
                times = [x/Fs for x in range(peak_index-start_index, end_index-start_index)]
                A = decay[breath][0][0]
                k = decay[breath][1][0]
                dec_curve = [A*np.exp(t*k) for t in times]
                #dec_curve = [-c for c in dec_curve]

                if(peak != 0):
                    rmse = root_mean_squared_error(expiration[peak_index:end_index], dec_curve, peak)
                    if(PRINTING):
                        print("rmse:{}".format(rmse))

                if(PLOTTING):
                    #times = [flow_middles[i]+peaks[i]+t for t in times]
                    plt.plot(range(peak_index, end_index), dec_curve, '--k', linewidth=2)
                    plt.plot([-e for e in expiration])
                    plt.plot(range(start_index, end_index), curve)
                    plt.show()

                    lin_fit = [grad[breath][0]*v + grad[breath][1] for v in volume_curve]
                    plt.plot(volume, expiration)
                    plt.plot(volume_curve, expiration[start_index:end_index], 'orange')
                    plt.plot(volume_curve, lin_fit, 'r')

                    plt.xlabel('volume')
                    plt.ylabel('flow')
                    plt.show()

        else:
            print('Less than 3 points to find decay. Can\'t do it')
            decay[breath] = np.nan


    return([grad, rvalues])


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


external_resistance = 0.3 # Dummy until value determined


for j in range(len(subject)):
    print('\nNEW subject\n--------\n')

    for i in range(len(test_type)):
        decays = []

        for k in range(2):

            print('\nNEW FILE\n--------\n')

            filename = path + subject[j] + '/normalised_' + subject[j] + test_type[i][k]
            print(filename)
            print('\n')

            try:
                # Extract data to empty arrays.
                # This will fail if file does not exist
                pressure_lowRes = []
                pressure_highRes = []
                time = []
                extract_data_arrays_from_file(filename, pressure_lowRes, pressure_highRes, time)

                # Flow and mouth pressure both calculated from spirometer pressure
                flow = convert_spirometer_pressure_to_flow(pressure_highRes)
                #mouth_pressure = convert_spirometer_pressure_to_mouth_pressure(pressure_highRes, conversion[k])

                plt.plot(pressure_highRes, 'b')
                plt.plot(pressure_lowRes, 'g')
                plt.plot(flow, 'r')
                plt.show()


                # Find the start and end of expiration
                flow_splits = split_breaths(flow, peak_height=0.025, Fs=300, plot=False)
                flow_middles = flow_splits[1]
                flow_stops = flow_splits[2]

                # return all decay rates of expiration
                # Some will be bull, because not held breaths
                results = calculate_decay_rates(flow, flow_middles, flow_stops)
                decay_rates = results[0]
                decay_rvalue = results[1]

            except IOError:
                print("File {} does not exist".format(filename))


            SF = 0.97 #sensitivity factor

            #Basic table print of all results (SF highlighted with *)
            print('\n')
            print('____ results for ____')
            print(filename)
            print('________')
            for n in range(len(decay_rates)):
                if(abs(decay_rvalue[n][0]) > SF):
                    print(decay_rates[n][0], decay_rates[n][1], decay_rvalue[n], '*')
                else:
                    print(decay_rates[n][0], decay_rates[n][1], decay_rvalue[n])

            used_grad = []
            count = 0
            for n in range(len(decay_rates)):
                if(abs(decay_rvalue[n][0]) > SF):
                    count += 1
                    used_grad.append(decay_rates[n][0])
            if(count > 0):
                mean_grad = sum(used_grad)/count
                print(' mean = {}'.format(mean_grad))
                SD = np.std(used_grad)
                print(' standard dev = {}'.format(SD))


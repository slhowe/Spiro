#!/bin/bash

# Import my extensions
import sys
sys.path.insert(0, '/home/sarah/Documents/Spirometry/python/extensions')

import csv
import matplotlib.pyplot as plt
from sympy import *
import numpy as np
from numpy import array
from numpy.linalg import lstsq
from breath_analysis import split_breaths

path = './'
files = [
        'data_recording.csv',
        'mask_example_4wRX_then_4woRX.csv',
        'mask_example.csv',
        'mask_example_Rx.csv'
         ]

relationship = [
        [-6.4, 1],
        [-6.4, 1],
        [-1.5, 1.25]
        ]

# Create data classes
spir_pressure = []
mask_pressure = []
time = []

for filename in range(0, 2):
    # Store pressure and flow data
    fname = path + files[filename]
    with open(fname, 'rb') as csvfile:
        reader = csv.reader(csvfile, delimiter = ',')

        # skip header line
        header = reader.next()

        for row in reader:
            spir_point = float(row[0])
            mask_point = float(row[1])
            time_point = float(row[2])

            # Convert spir_pressure to estimated mouth pressure
           # if(spir_point < 0):
           #     spir_point = relationship[filename][0]*spir_point**2 + relationship[filename][1]*spir_point

            spir_pressure.append(spir_point)
            mask_pressure.append(mask_point)
            time.append(time_point)

    plt.plot(spir_pressure)
    plt.plot(mask_pressure)
    plt.show()

    # Have the data not analysis
    # Split breaths
    # Find starts, middles and ends of flow data
    flow_splits = split_breaths(spir_pressure, peak_height=0.015, Fs=300, plot=False)
    flow_starts = flow_splits[0]
    flow_middles = flow_splits[1]
    flow_stops = flow_splits[2]

    decays = [[] for i in range(len(flow_middles))]
    peaks = [0 for i in range(len(flow_middles))]
    decay_end = [0 for i in range(len(flow_middles))]
    decay_start = [0 for i in range(len(flow_middles))]

    def get_decay_rate(pressure, i, flow_starts, flow_middles, flow_stops):
        #Find peak expiratory pressure
        # get min value between middle and end
        peak = min(pressure[flow_middles[i]:flow_stops[i]])
        peak_index = pressure[flow_middles[i]:flow_stops[i]].index(peak) #index relative end insp (flow middle)
        peaks[i] = peak_index
        print("peak:{}".format(peak))
        print("peak_index:{}".format(peak_index))

        # Find first index of value 0 relative to peak
        # If there isn't one, find the length from peak to end point
#        try:
#            first_zero = pressure[flow_middles[i]+peak_index:flow_stops[i]].index(0)
#        except(ValueError):
#            first_zero = flow_stops[i] - (flow_middles[i] + peak_index) - 2

        #linear relationship
        decay = [0, 0]
        decay_start = -1
        decay_end = -1
        j = peak_index
        # k_indexes len of expiration
        num_datapoints = flow_stops[i] - flow_middles[i]
        while j < num_datapoints:
            if(decay_start < 0):
                if(pressure[flow_middles[i] + j] > -0.10):
                    decay_start = j
                    j = j + 4
            else:
                if(pressure[flow_middles[i] + j] >= -0.045):
                    decay_end = j - 1
                    j = num_datapoints
            j+=1

        print("decay_start:{}".format(decay_start))
        print("decay_end:{}".format(decay_end))

#        # Find half way point to zero or end, relative to peak
#        drop_percentage = 63
#        thresh = (peak*(100-drop_percentage))/100
#        halfway_point = 0
#        search_index = 0
#        while search_index < first_zero:
#            value = pressure[flow_middles[i]+peak_index+search_index]
#            if value > thresh:
#                halfway_point = search_index
#                search_index = first_zero
#            search_index += 1
#
#        #halfway_point = (first_zero)/2
#        print("first_zero:{}".format(first_zero+peak_index))
#        print("halfway_point:{}".format(halfway_point+peak_index))
#
#        # Define range to determine decay rate here
#        # Relative to flow_middle
#        #decay_start = halfway_point + peak_index
#        decay_start = peak_index
#        decay_end = first_zero + peak_index
#        decay_end = halfway_point + peak_index

        # Grab the curve of interest and make it all positive
        # (assume no positive section in exp range)
        curve = pressure[flow_middles[i]+decay_start:flow_middles[i]+decay_end]
        curve = [np.abs(c) for c in curve]

        # Fit curve to exponential model
        # ln(Y) = At + ln(b)
        # for y(t) = e^(At) + b
        if(len(curve) > 3):
            Fs = 90.0
            measurements = array([np.log(curve)])

            # set up array of data multiplied by constants
            one_array = [1 for val in range(len(measurements.T))]
            times = [j/Fs for j in range(len(measurements.T))]
            independents = array([one_array,times])

            # Least squares fit to data
            result = lstsq(independents.T, measurements.T)

            # Squared Euclidean 2-norm for each col in (b - a*x)
            residual = result[1]
            print('       residual on line fit is {}'.format(residual[0]))

            # Parameters - offset and decay constant
            constants = result[0]
            constants[0] = np.exp(constants[0])
            print('       offset is {}, decay is {}'.format(constants[0], constants[1]))

            decay = constants

            # Optional plot range and decay for each breath
            if(1):
                # Remake curve from best fit
                # Length of curve is from decay_start to flow_stop point
                Fs = float(90)
                times = [x/Fs for x in range(flow_stops[i]-(flow_middles[i]+decay_start))]
                A = decay[0]
                k = decay[1]
                dec_curve = [A*exp(t*k) for t in times]
                dec_curve = [-c for c in dec_curve]

                #times = [flow_middles[i]+peaks[i]+t for t in times]
                plt.plot(range(decay_start,(flow_stops[i]-flow_middles[i])), dec_curve, '--k', linewidth=2)

                plt.plot(pressure[flow_middles[i]:flow_stops[i]], '-r')
                plt.plot(range(decay_start,decay_end), pressure[flow_middles[i]+decay_start:flow_middles[i]+decay_end], '-g', linewidth=2)
                plt.show()

        return (decay_start, decay_end, decay)



    # Fit decay to each expiration
    for i in range(len(flow_middles)):
        params = get_decay_rate(spir_pressure, i, flow_starts, flow_middles, flow_stops)
        decay_start[i] = params[0]
        decay_end[i] = params[1]
        decays[i] = params[2]

    #plot
    # Plot all data and show where breath middles & ends are located
    plt.plot(flow)
    dots = [flow[i] for i in flow_middles]
    plt.plot(flow_middles, dots, 'ro')
    dots = [flow[i] for i in flow_stops]
    plt.plot(flow_stops, dots, 'ko')
    plt.legend(['data', 'middles', 'ends'])
    plt.show()

    # Plot all data with decays of interest shown
    plt.plot(#pressure, 'b',
             flow, 'r',
             )

    # Define breaths of interest
    boi = [
           [0, 1, 2, 3, 4]
           ]

    # decays calculated from decay_start
    # decay start is relative to flow_middle
    for b in boi:
        print('\n')
        for i in b:
            Fs = float(90)
            times = [x/Fs for x in range(flow_stops[i]-(flow_middles[i]+decay_start[i]))]
            A = decay[i][0][0]
            k = decay[i][1][0]
            print("decay = {}".format(k))
            curve = [A*exp(t*k) for t in times]
            curve = [-c for c in curve]

            # Plot the decay overtop of the original data
            times = range(flow_middles[i]+decay_start[i], flow_stops[i])
            plt.plot(times, curve, '--k', linewidth=2)

    plt.xticks(fontsize=20)
    plt.yticks(fontsize=20)
    #plt.xlim(0, 2900)
    plt.ylim(-1, 1)
    plt.grid()
    plt.legend([#"Pressure",
                "Flow",
                ])
    plt.xlabel("Data Point", fontsize=32)
    plt.ylabel("Flow (L/s)", fontsize=32)
    plt.show()

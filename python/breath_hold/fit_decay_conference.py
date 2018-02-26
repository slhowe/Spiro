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
        'data_recording.csv'
         ]

# Create data classes
flow = []
pressure = []
time = []

# Store pressure and flow data
filename = path + files[0]
with open(filename, 'rb') as csvfile:
    reader = csv.reader(csvfile, delimiter = ',')

    # skip header line
    header = reader.next()

    for row in reader:
        pressure_point = float(row[0])
        flow_point = float(row[1])
        time_point = float(row[2])

        pressure.append(pressure_point)
        flow.append(flow_point)
        time.append(time_point)

# Have the data not analysis
# Split breaths
# Find starts, middles and ends of flow data
flow_splits = split_breaths(flow, plot=False)
flow_starts = flow_splits[0]
flow_middles = flow_splits[1]
flow_stops = flow_splits[2]

decay = [[] for i in range(len(flow_middles))]
peaks = [0 for i in range(len(flow_middles))]

# Fit decay to each expiration
for i in range(len(flow_middles)):
    #Find peak expiratory pressure
    # get min value between middle and end
    peak = min(flow[flow_middles[i]:flow_stops[i]])
    peak_index = flow[flow_middles[i]:flow_stops[i]].index(peak) #index relative end insp (flow middle)
    peaks[i] = peak_index
    print("peak:{}".format(peak))
    print("peak_index:{}".format(peak_index))

    # Find first index of value 0 relative to peak
    # If there isn't one, find the length from peak to end point
    try:
        first_zero = flow[flow_middles[i]+peak_index:flow_stops[i]].index(0)
    except(ValueError):
        first_zero = flow_stops[i] - (flow_middles[i] + peak_index)

    # Find half way point to zero or end, relative to peak
    halfway_point = (first_zero)/2
    print("first_zero:{}".format(first_zero+peak_index))
    print("halfway_point:{}".format(halfway_point+peak_index))

    # Define range to determine decay rate here
    # Relative to flow_middle
    decay_start = halfway_point + peak_index
    decay_end = first_zero + peak_index

    # Grab the curve of interest and make it all positive
    # (assume no positive section in exp range)
    curve = flow[flow_middles[i]+decay_start:flow_middles[i]+decay_end]
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

        decay[i] = constants

    if(1):
        # Remake curve from best fit
        # Length of curve is from decay_start to flow_stop point
        Fs = float(90)
        times = [x/Fs for x in range(flow_stops[i]-(flow_middles[i]+decay_start))]
        A = decay[i][0][0]
        k = decay[i][1][0]
        dec_curve = [A*exp(t*k) for t in times]
        dec_curve = [-c for c in dec_curve]

        #times = [flow_middles[i]+peaks[i]+t for t in times]
        plt.plot(dec_curve, '--k', linewidth=2)

        plt.plot(flow[flow_middles[i]:flow_stops[i]], '-r')
        plt.plot(range(decay_start,decay_end), flow[flow_middles[i]+decay_start:flow_middles[i]+decay_end], '-m')
        plt.show()


#plot
dots = [flow[i] for i in flow_middles]
plt.plot(flow)
plt.plot(flow_middles, dots, 'ro')
dots = [flow[i] for i in flow_stops]
plt.plot(flow_stops, dots, 'ko')
plt.show()

# Display only decaying breaths
plt.plot(#pressure, 'b',
         flow, 'r',
         )

boi = [1,4,6,8,14,17]
boi = [2,4,7,10,13,15,17]
for i in boi:
    Fs = float(90)
    times = [x/Fs for x in range(flow_stops[i]-(flow_middles[i]+peaks[i]))]
    A = decay[i][0][0]
    k = decay[i][1][0]
    print("decay = {}".format(k))
    curve = [A*exp(t*k) for t in times]
    curve = [-c for c in curve]

    #times = [flow_middles[i]+peaks[i]+t for t in times]
    times = range(flow_middles[i]+peaks[i],flow_stops[i])
    plt.plot(times, curve, '--k', linewidth=2)

plt.xticks(fontsize=20)
plt.yticks(fontsize=20)
plt.xlim(0, 2900)
plt.ylim(-1, 1)
plt.grid()
plt.legend([#"Pressure",
            "Flow",
            ])
plt.xlabel("Data Point", fontsize=32)
plt.ylabel("Flow (L/s)", fontsize=32)
plt.show()

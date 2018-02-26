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
        'example_1.csv'
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

# Fit decay to each expiration
decay = [[] for i in range(len(flow_middles))]
peaks = [0 for i in range(len(flow_middles))]
for i in range(len(flow_middles)):
    #Find peak expiratory pressure
    # get min value between middle and end
    peak = min(flow[flow_middles[i]:flow_stops[i]])
    print(peak)
    peak_index = flow[flow_middles[i]:flow_stops[i]].index(peak) #index relative end insp (flow middle)
    print(peak_index)
    peaks[i] = peak_index

#    end = max(flow[flow_middles[i]+peak_index:flow_stops[i]-5])
#    print(end)
#    end_index = flow[flow_middles[i]+peak_index:flow_stops[i]-5].index(end)
#    end_index -= peak_index
#    print(end_index)

    # End expiration range at half way through expiration
    # So only looking at first half of expiratory curve
    try:
        first_zero = flow[flow_middles[i]+peak_index:flow_stops[i]].index(0)
    except(ValueError):
        first_zero = flow_stops[i] - (flow_middles[i] + peak_index)
    end_index = peak_index + (first_zero)/2
    print(first_zero)
    #plt.plot(flow[flow_middles[i]:flow_stops[i]])
    #plt.show()

    # Grab the curve of interest and make it all positive (assume no positive section in exp range)
    curve = flow[flow_middles[i]+peak_index:flow_middles[i]+end_index]
    #plt.plot(curve)
    curve = [np.abs(c) for c in curve]
    #plt.plot(curve)
    #plt.show()

    if(len(curve) > 3):
        Fs = 90.0

        # ln(Y) = At + ln(b)
        # for y(t) = e^(At) + b
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

    if(0):
        # Remake curve from best fit
        Fs = float(90)
        times = [x/Fs for x in range(flow_stops[i]-(flow_middles[i]+peaks[i]))]
        A = decay[i][0][0]
        k = decay[i][1][0]
        dec_curve = [A*exp(t*k) for t in times]
        dec_curve = [-c for c in dec_curve]

        #times = [flow_middles[i]+peaks[i]+t for t in times]
        plt.plot(times, dec_curve, '--k', linewidth=2)
        plt.plot(dec_curve, '--k', linewidth=2)

        plt.plot(flow[flow_middles[i]+peak_index:flow_stops[i]], ':r')
        plt.plot(curve)
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

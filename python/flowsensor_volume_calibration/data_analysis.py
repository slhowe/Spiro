#!/bin/bash

# Import my extensions
import sys
sys.path.insert(0, '/home/sarah/Documents/Spirometry/python/extensions')

import csv
import matplotlib.pyplot as plt
import numpy as np
from numpy import array
from filters import hamming
from calculus import integral
from breath_analysis import split_breaths
from numpy.linalg import lstsq
from scipy import stats

path = './'
training_files = ['blue_flowsensor_calibration.csv',
                  'blue_flowsensor_calibration_3.csv']

test_file = ['blue_flowsensor_calibration_2.csv']


# Create data classes
pressure_high = []
flow = []
new_flow = []
pressure_low= []
time = []

for i in range(len(training_files)):
     # Create data classes
    pressure_high = []
    flow = []
    new_flow = []
    pressure_low= []
    time = []

   # Store pressure and flow data
    filename = path + training_files[i]
    with open(filename, 'rb') as csvfile:
        reader = csv.reader(csvfile, delimiter = ',')

        # skip header line
        header = reader.next()

        for row in reader:
            dp_pressure_high = -float(row[0])
            dp_time = float(row[3])

            pressure_high.append(dp_pressure_high)
            time.append(dp_time)


#    plt.plot(dataset.time, dataset.pressure, 'b',
#             dataset.time, dataset.flow, 'r',
#             )

    # Filter high resolution pressure data
    fc = 5#filter cutoff
    bw = 6#filter bandwidth (taps?)
    pressure = hamming(pressure_high, fc, 300, bw, plot=False)

    pressure_squared = [p**2 for p in pressure]

    # Find the end point of the 3L sections
    # Start with insp 3L, then exp 3L
    # Each insp and exp needs different fit
    # Ignore any data with pressure greater than 290

    # if pressure == 0:
    # at a flat point
    # start of new section
    # find end of section - next P=0 at different volume
    P_integral = integral(pressure, 300)
    P_integral = np.real(P_integral).tolist()
    P2_integral = integral(pressure**2, 300)
    P2_integral = np.real(P2_integral).tolist()

    starts = []
    min_stop_len = 100 #1/3 of a second must be flat
    for index in range(len(P_integral)-min_stop_len):
        if (abs(P_integral[index] - P_integral[index+min_stop_len]) < 0.1): #If flat for a while
            if(len(starts) > 0):
                if(abs(P_integral[index] - P_integral[starts[-1]]) > 200): #If there has been a big change in cumulative pressure
                    starts.append(index)
            else:
                starts.append(index)
        if (abs(pressure[index]) > 294):
            break

    print(starts)
    start_times = [time[s] for s in starts]
    start_pressures = [P_integral[s] for s in starts]


    # Work out volumes now
    pressures = []
    for start_index in range(len(starts)-1):
        # difference between start and end P_integral to hive total integral over 3L period
        start = starts[start_index]
        end = starts[start_index+1]
        P = P_integral[end] - P_integral[start]
        P2 = P2_integral[end] - P2_integral[start]
        pressures.append([P, P2])

    # Combine all the positive parts and negative parts to get datapoints for LSTSQ fit
    pos_V = [0]
    neg_V = [0]
    pos = [0]
    pos2 = [0]
    neg = [0]
    neg2 = [0]
    for index in range(len(pressures)):
        if pressures[index][0] > 0:
            res = pressures[index][0] + pos[-1]
            res2 = pressures[index][1] + pos2[-1]
            pos.append(res)
            pos2.append(res2)
            pos_V.append(pos_V[-1] + 3)
        else:
            res = pressures[index][0] + neg[-1]
            res2 = pressures[index][1] + neg2[-1]
            neg.append(res)
            neg2.append(res2)
            neg_V.append(neg_V[-1] - 3)

    print(pressures)
    print(pos)
    print(pos2)
    print(neg)
    print(neg2)
    print(pos_V)
    print(neg_V)

    # No offset, because Q = 0 at passive state
    # V = A*int(P^2) + B*int(P)

    # Least squares fit to data
    independents = array([pos2,pos])
    dependents = array([pos_V])
    result = lstsq(independents.T, dependents.T)
    print(result)
    pos_A = result[0][0][0]
    pos_B = result[0][1][0]
    pos_line = [pos_A*pos2[index] + pos_B*pos[index] for index in range(len(pos_V))]

    independents = array([neg2,neg])
    dependents = array([neg_V])
    result = lstsq(independents.T, dependents.T)
    print(result)
    neg_A = result[0][0][0]
    neg_B = result[0][1][0]
    neg_line = [neg_A*neg2[index] + neg_B*neg[index] for index in range(len(neg_V))]


    # REPEAT FOR LINEAR ASSUMPTION
    independents = array([pos])
    dependents = array([pos_V])
    result = lstsq(independents.T, dependents.T)
    print(result)
    pos_A_lin = result[0][0][0]

    #independents = array([neg2,neg])
    independents = array([neg])
    dependents = array([neg_V])
    result = lstsq(independents.T, dependents.T)
    print(result)
    neg_A_lin = result[0][0][0]


    # remake volume completely:
    flow = [0]*len(pressure)
    flow_lin = [0]*len(pressure)
    for index in range(len(flow)):
        if (pressure[index] >= 0):
            flow[index] = pos_A*(pressure[index])**2 + pos_B*pressure[index]
            flow_lin[index] = pos_A_lin*pressure[index]
        else:
            flow[index] = neg_A*(pressure[index])**2 + neg_B*pressure[index]
            flow_lin[index] = neg_A_lin*pressure[index]

    volume = integral(flow, 300)
    volume_lin = integral(flow_lin, 300)

# PLOTSSSSSSSS ################################################################
    plt.plot(time, pressure_high, 'b',)
    plt.plot(time, pressure, 'orange', linewidth=2,)
    plt.plot(time, P_integral, 'm')
    #plt.plot(time, P2_integral, 'c')
    plt.plot(start_times, start_pressures, 'or')
    plt.xlabel('Time')
    plt.grid()
    plt.legend(["Pressure high",
                "filtered pressure",
                "integral of pressure",
                ])
    plt.show()

    plt.plot(pos_V, pos, 'ro')
    plt.plot(neg_V, neg, 'co')
    plt.xlabel('volume')
    plt.ylabel('pressure')
    plt.grid()
    plt.show()

    plt.plot(time, pressure, 'orange')
    plt.plot(time, flow, 'r')
    plt.plot(time, volume, 'g')
    plt.plot(time, volume_lin, 'c')
    plt.xlabel('Time')
    plt.grid()
    plt.show()
# PLOTSSSSSSSS ################################################################




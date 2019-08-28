#!/bin/bash

# Import my extensions
import sys
sys.path.insert(0, '/home/sarah/Documents/Spirometry/python/extensions')

import csv
import matplotlib.pyplot as plt
from filters import hamming
from calculus import integral
from breath_analysis import split_breaths
import numpy as np
from numpy.linalg import lstsq


path = './'
files = [
        #'data_recording.csv',
        'nozzle_0095_adj.csv',
        'nozzle_0105_adj.csv',
        'nozzle_0125_adj.csv',
        'nozzle_0095.csv',
        'nozzle_0105.csv',
        'nozzle_0125.csv',
        'RX_human_washer11.csv',
        'RX_human_washer10.csv',
        'RX_human_washer9.csv',
        'RX_human_washer8.csv',
        'RX_human_washer7.csv',
        'RX_human_washer6.csv',
        'RX_human_washer5.csv',
        'RX_human_washer4.csv',
        'RX_human_washer3.csv',
        'RX_human_washer2.csv',
        'RX_human_washer1.csv',
        'RX_human_8_5mm8.csv',
        'RX_human_8_5mm1.csv',
        'RX_human_8_5mm2.csv',
        'RX_human_8_5mm3.csv',
        'RX_human_8_5mm4.csv',
        'RX_human_8_5mm5.csv',
        'RX_human_8_5mm6.csv',
        'RX_human_8_5mm7.csv',
        'RX_human_8_5mm9.csv',
        #'lung_pump_cap.csv'
         ]

# Create data classes
pressure_high = []
flow = []
new_flow = []
pressure_low= []
time = []

for i in range(len(files)):
    pressure_high = []
    flow = []
    new_flow = []
    pressure_low= []
    time = []

    # Store pressure and flow data
    filename = path + files[i]
    with open(filename, 'rb') as csvfile:
        reader = csv.reader(csvfile, delimiter = ',')

        # skip header line
        header = reader.next()

        for row in reader:
            dp_pressure_high = float(row[0]) # Rsp
            dp_flow = float(row[1])
            dp_pressure_low= float(row[2]) # Rsp + Rx
            dp_time = float(row[3])

            pressure_high.append(dp_pressure_high)
            flow.append(dp_flow)
            pressure_low.append(-dp_pressure_low)
            time.append(dp_time)


    # Filter high resolution pressure data
    fc = 15#filter cutoff
    bw = 10#filter bandwidth (taps?)
#    filt_pressure = hamming(pressure_high, fc, 300, bw, plot=False)
#    filt_pressure_Rx = hamming(pressure_low, fc, 300, bw, plot=False)
#
#    for p in filt_pressure:
#        if p >= 0:
#            new_flow.append(0.0071*p)
#        else:
#            new_flow.append(0.0075*p)
    filt_pressure = pressure_high
    filt_pressure_Rx = pressure_low
    new_flow = flow

    Rtot = []
    Rspir = []
    for index in range(len(flow)):
        if new_flow[index] != 0:
            Rspir.append(filt_pressure[index]/new_flow[index])
            Rtot.append(filt_pressure_Rx[index]/new_flow[index])
        else:
            Rspir.append(1)
            Rtot.append(1)


    volume = integral(flow, 300)
    new_volume = integral(new_flow, 300)

    plt.plot(time, pressure_high, 'b',)
    plt.plot(time, filt_pressure, 'orange', linewidth=2,)
    plt.plot(time, flow, 'r',
             time, new_flow, 'pink',
             time, pressure_low, 'c',
             time, volume, 'm',
             time, new_volume, 'purple',
             )

    plt.xlabel('time')
    plt.grid()
    plt.legend(["Pressure high",
                "filtered pressure",
                "Flow",
                "filtered flow",
                "pressure low",
                "volume",
                "integral filtered flow",
                ])
    plt.show()

#    plt.plot(time, filt_pressure, 'b')
#    plt.plot(time, filt_pressure_Rx, 'k')
#    plt.plot(time, new_flow, 'r')
#    plt.plot(time, Rspir, 'm')
#    plt.plot(time, Rtot, 'purple')
#    plt.ylim(1, -1)
#    plt.grid()
#    plt.legend([
#                "filtered pressure (spir)",
#                "filtered pressure (Rx + spir)",
#                "Flow (from filtered data)",
#                "resistance (spir)",
#                "resistance (Rx + spir)",
#                ])
#    plt.show()

    # grab the top bit of pressure and flow
    top_pressure = []
    top_flow = []
    MIN_PRESSURE = 0
    MAX_PRESSURE = 1000
    MIN_FLOW = 0.2
    MAX_FLOW = 0.7
    for index in range(len(filt_pressure_Rx)):
        if filt_pressure_Rx[index] > MIN_PRESSURE and filt_pressure_Rx[index] < MAX_PRESSURE:
            if new_flow[index] > MIN_FLOW and new_flow[index] < MAX_FLOW:
                top_pressure.append(filt_pressure_Rx[index])
                top_flow.append(new_flow[index])

    ones = [1]*len(top_pressure)
    top_flow2 = [q for q in top_flow]
    measurements = np.array([top_pressure])
    independents = np.array([top_flow2, ones])
    res = lstsq(independents.T, measurements.T)
    print('LOOK HEREvvvvvv')
    print(res)
    line = [res[0][0][0]*q + res[0][1][0] for q in top_flow2]
    #line = [80*q + res[0][1][0] for q in top_flow2]

    top_R = [top_pressure[i]/top_flow[i] for i in range(len(top_pressure))]
    print('\nNOTE:\nThis is only valid for expiration for pressures greater than 100 Pa')
    print('\nRESISTANCE RX + RSPIR IS:')
    print(np.mean(top_R))
    print('\nRSPIR is :\n1/0.0075 = 133.3')
    print('\nso RSPIR is:')
    print(np.real(np.mean(top_R)) - 1/0.0075)

    plt.plot(new_flow, filt_pressure, 'bo')
    plt.plot(new_flow, filt_pressure_Rx, 'ro')
    plt.plot(top_flow, top_pressure, 'mo')
    plt.plot(top_flow, line, 'g')
    plt.xlabel('flow')
    plt.ylabel('pressure')
    plt.grid()

    plt.figure(2)
    plt.plot(top_flow, top_R, 'oc')
    plt.xlabel('R')
    plt.grid()

    plt.show()



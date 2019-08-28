#!/bin/bash

# Import my extensions
import sys
sys.path.insert(0, '/home/sarah/Documents/Spirometry/python/extensions')

import csv
import matplotlib.pyplot as plt
from filters import hamming
from calculus import integral
from breath_analysis import split_breaths

path = './'
files = [
        'data_recording.csv'
         ]

# Create data classes
pressure_high = []
flow = []
new_flow = []
pressure_low= []
time = []

for i in range(len(files)):
    # Store pressure and flow data
    filename = path + files[i]
    with open(filename, 'rb') as csvfile:
        reader = csv.reader(csvfile, delimiter = ',')

        # skip header line
        header = reader.next()

        for row in reader:
            dp_pressure_high = -float(row[0])
            dp_flow = float(row[1])
            dp_pressure_low= float(row[2])
            dp_time = float(row[3])

            pressure_high.append(dp_pressure_high)
            flow.append(dp_flow)
            pressure_low.append(dp_pressure_low)
            time.append(dp_time)


#    plt.plot(dataset.time, dataset.pressure, 'b',
#             dataset.time, dataset.flow, 'r',
#             )

    # Filter high resolution pressure data
    fc = 15#filter cutoff
    bw = 10#filter bandwidth (taps?)
    filt_pressure = hamming(pressure_high, fc, 300, bw, plot=False)


    volume = integral(flow, 300)
    volume = integral(filt_pressure, 300)

    f, ax = plt.subplots(2, sharex=True)

    ax[0].plot(time, pressure_low, 'blue', linewidth=2,)
    ax[1].plot(time, flow, 'g',linewidth=2,)
    #         #time, pressure_low, 'c',
    #         time, volume, 'm',
    #         )

    ax[0].grid()
    ax[1].grid()
    ax[1].set_xlabel('Time (ms)', fontsize=32)
    ax[0].set_ylabel('Pressure (Pa)', fontsize=32)
    ax[1].set_ylabel('Flow (L/s)', fontsize=32)
#    plt.legend(["Pressure high",
#                "filtered pressure",
#                "Flow",
#                "pressure low",
#                "integral of flow",
#                ])
    plt.show()

    print("This is the plot of the High pressure measurement (low res one) vs volume")
    plt.plot(volume, pressure_low)

    plt.grid()
    plt.xlabel("volume")
    plt.ylabel("pressure (lung estimate when occluded)")
    plt.show()

    print("This is the plot of the High pressure measurement")
    plt.plot(time, pressure_low)

    plt.grid()
    plt.xlabel("time")
    plt.ylabel("pressure (lung estimate when occluded)")
    plt.show()


    print("This is the plot of the High pressure measurement (low res one) vs flow")
    plt.plot(volume, flow)

    plt.grid()
    plt.ylabel("flow")
    plt.xlabel("volume (lung estimate when occluded)")
    plt.show()

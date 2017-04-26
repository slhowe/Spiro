#!/bin/bash

# Import my extensions
import sys
sys.path.insert(0, '/home/sarah/Documents/Spirometry/python/extensions')
from breath_analysis import calc_flow_delay, split_breaths
from filters import hamming
from calculus import integral
from data_struct import DataStore as Data
import csv
import matplotlib.pyplot as plt
from numpy import array, sign, log
from numpy.linalg import lstsq

#path = '/home/sarah/Documents/Spirometry/python/flowComparisons/'
#files = ['Calibration_set']
path = '/home/sarah/Documents/Spirometry/data/'
files = ['Loops_1.csv', 'Loops_3.csv']

# Create data classes
reset = Data('reset')
reset.clear_dataset_registry()
dataset_1 = Data('dataset_1')
dataset_2 = Data('dataset_2')

# Points I prepared earlier =D
dataset_2_start_indices = [2465, 5737, 8939]
dataset_2_end_indices = [2865, 6150, 9340]

# Store pressure and flow data
file_index = 0
for dataset in Data:
    filename = path + files[file_index]
    with open(filename, 'rb') as csvfile:
        reader = csv.reader(csvfile, delimiter = ',')

        # skip header line
        header = reader.next()

        for row in reader:
            pressure = float(row[0])
            flow = float(row[1])

            dataset.pressure.append(pressure)
            dataset.flow.append(flow)
    file_index += 1

    # measured
    delay = calc_flow_delay(dataset.pressure[:len(dataset.flow)/10],
                            dataset.flow[:len(dataset.flow)/10],
                            125,
                            plot=False
                            )

    pressure = dataset.pressure
    flow = dataset.flow
    fpressure = hamming(pressure, 5, 125, 2, plot=True)
    fflow= hamming(flow, 5, 125, 2, plot=False)

    # Uncomment if splitting data into breaths
    if(1):
        flow_splits = split_breaths(dataset.flow)
        breaths = len(flow_splits[0])

        for breath in range(9):#(breaths - 1):
            # Get pressure and filter
            pres = pressure[flow_splits[0][breath]-delay:flow_splits[0][breath+1]-delay]
            fpres = fpressure[flow_splits[0][breath]-delay:flow_splits[0][breath+1]-delay]

            #Get flow and filter
            flw = flow[flow_splits[0][breath]:flow_splits[0][breath+1]]
            fflw = fflow[flow_splits[0][breath]:flow_splits[0][breath+1]]

            # Flip pressure
            pres = [-f for f in pres]
            fpres = [-f for f in fpres]

            # Get volume
            vol = integral(flw, 125)

            f, (ax1, ax2, ax3) = plt.subplots(3, sharex=True)
            ax1.plot(pres, 'b',
                    fpres, 'g'
                    )
            ax1.legend(["Original",
                        "Filtered"
                        ])
            ax1.set_title("Pressure")
            ax1.grid()

            ax2.plot(flw, 'r.-')
            ax2.plot(fflw, 'k')
            ax2.legend(["Original",
                        "Filtered"
                        ])
            ax2.set_title("Flow")
            ax2.grid()

            ax3.plot(vol, 'g.-')
            ax3.set_title("Volume")
            ax3.grid()

            plt.show()

    volume = integral(flow, 125)
    f, (ax1, ax2, ax3) = plt.subplots(3, sharex=True)
    ax1.plot(pressure, 'b',
            fpressure, 'g'
            )
    ax1.set_title("Pressure")
    ax1.legend(["Original",
                "Filtered"
                ])
    ax1.grid()

    ax2.plot(flow, 'r')
    ax2.set_title("Flow")
    ax2.grid()

    ax3.plot(volume, 'g')
    ax3.set_title("Volume")
    ax3.grid()

    plt.show()

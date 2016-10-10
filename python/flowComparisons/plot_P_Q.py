#!/bin/bash

# Import my extensions
import sys
sys.path.insert(0, '/home/sarah/Documents/Spirometry/python/extensions')

import csv
import matplotlib.pyplot as plt
from filters import hamming
from data_struct import DataStore as Data
from numpy import array, sign
from numpy.linalg import lstsq

#path = '/home/sarah/Documents/Spirometry/python/flowComparisons/'
#files = ['Calibration_set']
path = '/home/sarah/Documents/Spirometry/data/'
files = ['Normal_1.csv']

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

    offset = 26

    # measured
    pressure = dataset.pressure[0:-offset]
    flow = dataset.flow[offset:]
    fpressure = hamming(pressure, 40, 125, 4, plot=True)

    plt.plot(pressure, 'b',
             flow, 'g',
             fpressure, 'r-o',
             )

    plt.grid()
    plt.legend(["Pressure",
                "Flow"
                ])
    plt.show()

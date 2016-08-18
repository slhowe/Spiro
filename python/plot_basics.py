#!/bin/bash

# Import my extensions
import sys
sys.path.insert(0, '/home/sarah/Documents/Spirometry/python/extensions')

import csv
import matplotlib.pyplot as plt
from data_struct import DataStore as Data

path = '/home/sarah/Documents/Spirometry/python/calibration/calibration_files/'
files = [
        '100_0_Lmin-1',
        '21_5_Lmin-1',
        '43_0_Lmin-1',
        '61_0_Lmin-1',
        '79_8_Lmin-1',
        '11_85_Lmin-1',
        '30_6_Lmin-1',
        '43_75_Lmin-1',
        '67_45_Lmin-1',
        '8_45_Lmin-1',
        '18_55_Lmin-1',
        '30_8_Lmin-1',
        '51_0_Lmin-1',
        '70_18_Lmin-1',
        '88_1_Lmin-1',
        '2_11_Lmin-1',
        '39_8_Lmin-1',
        '55_85_Lmin-1',
        '78_15_Lmin-1',
        '90_5_Lmin-1'
         ]

# Create data classes
reset = Data('reset')
reset.clear_dataset_registry()

for i in range(len(files)):
    # Store pressure and flow data
    dataset = Data('dataset_1')
    filename = path + files[i]
    with open(filename, 'rb') as csvfile:
        reader = csv.reader(csvfile, delimiter = ',')

        # skip header line
        header = reader.next()

        for row in reader:
            pressure = float(row[1])
            flow = float(row[0])

            dataset.pressure.append(pressure)
            dataset.flow.append(flow)

    Fs = 90
    time = [t/float(Fs) for t in range(len(dataset.pressure))]

    plt.plot(time, dataset.pressure, 'b',
             time, dataset.flow, 'r',
             )

    plt.grid()
    plt.legend(["Pressure",
                "Flow",
                ])
    plt.show()

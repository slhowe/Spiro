#!/bin/bash

# Import my extensions
import sys
sys.path.insert(0, '/home/sarah/Documents/Spirometry/python/extensions')

import csv
import matplotlib.pyplot as plt
from data_struct import DataStore as Data
from calculus import integral

path = './'
files = [
        'data_recording.csv'
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
            pressure = float(row[0])
            flow = float(row[1])
            time = float(row[3])

            dataset.pressure.append(pressure)
            dataset.flow.append(flow)
            dataset.time.append(time)


#    plt.plot(dataset.time, dataset.pressure, 'b',
#             dataset.time, dataset.flow, 'r',
#             )
    volume = integral(dataset.flow, 60)
    plt.plot(dataset.pressure, 'b',
             dataset.flow, 'r',
             volume, 'm',
             )

    plt.grid()
    plt.legend(["Pressure",
                "Flow",
                "volume"
                ])
    plt.show()

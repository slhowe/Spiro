#!/bin/bash

# Import my extensions
import sys
sys.path.insert(0, '/home/sarah/Documents/Spirometry/python/extensions')

import csv
import matplotlib.pyplot as plt
from data_struct import DataStore as Data

path = '/home/sarah/Documents/Spirometry/python/'
files = [
        'spir_record.csv'
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

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
        'sensor_drift_tidal_breathing.csv',
        'test_syringe_good_out.csv',
        'test_syringe_good_out_reversed_sensor.csv',

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
    volume = integral(dataset.flow, 300)
    volume2 = integral(dataset.pressure, 300)
    plt.plot(dataset.time, dataset.pressure, 'b',
             dataset.time, dataset.flow, 'r',
             dataset.time, volume, 'm',
             dataset.time, volume2, 'orange',
             )

    plt.grid()
    plt.legend(["Pressure",
                "Flow",
                "integral of flow",
                "integral of pressure",
                ])
    plt.show()
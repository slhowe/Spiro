#!/bin/bash

# Import my extensions
import sys
sys.path.insert(0, '/home/sarah/Documents/Spirometry/python/extensions')

import numpy as np
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

mask_pressure_array = []
spir_pressure_array = []
relat_array = []

Fs = 300 #Hz

for i in range(len(files)):
    # Store pressure and flow data
    dataset = Data('dataset_1')
    filename = path + files[i]
    with open(filename, 'rb') as csvfile:
        reader = csv.reader(csvfile, delimiter = ',')

        # skip header line
        header = reader.next()

        for row in reader:
            spir_pressure = float(row[0])
            mask_pressure = float(row[1])
            time = float(row[2])

            mask_pressure_array.append(mask_pressure)
            spir_pressure_array.append(spir_pressure)
            dataset.time.append(time/1000)
            if(spir_pressure != 0):
                relat_array.append(mask_pressure/spir_pressure)
            else:
                relat_array.append(np.nan)


#    plt.plot(dataset.time, dataset.pressure, 'b',
#             dataset.time, dataset.flow, 'r',
#             )
    plt.plot(dataset.time, spir_pressure_array, 'b',
             dataset.time, mask_pressure_array, 'g',
             dataset.time, relat_array, 'r',
             )

    plt.grid()
    plt.legend(["spir_pressure",
                "mask_pressure",
                "mask/spir"
                ])
    plt.show()

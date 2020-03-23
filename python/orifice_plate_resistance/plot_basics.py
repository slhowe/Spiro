#!/bin/bash

# Import my extensions
import sys
sys.path.insert(0, '/home/sarah/Documents/Spirometry/python/extensions')

import csv
import numpy as np
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
high_pressure = []

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
            pres = float(row[2])
            time = float(row[3])

            dataset.pressure.append(pressure)
            dataset.flow.append(flow)
            high_pressure.append(pres)
            dataset.time.append(time)


    resistance = [0]*len(dataset.pressure)
    for j in range(len(resistance)):
        if(dataset.flow[j] != 0):
            resistance[j] = high_pressure[j]/dataset.flow[j]

    volume = integral(dataset.flow, 300)
    plt.plot(dataset.time, dataset.pressure, 'b',
             dataset.time, dataset.flow, 'r',
             dataset.time, high_pressure, 'k',
             dataset.time, volume, 'm',
             dataset.time, resistance, 'c',
             )

    plt.grid()
    plt.legend(["Pressure across spir",
                "Flow",
                "Pressure across orifice",
                "integral of flow",
                "resistance",
                ])
    plt.show()

    # Looking at breathing out (neg pressure)
    # Because for breathing in it is measuring drop across resistances going back in
    #(assuming lung at atmosphere, which it isn't)
    new_res = []
    new_pres = []
    new_time = []
    for j in range(len(resistance)):
        if high_pressure[j] > 0:
            new_res.append(resistance[j])
            new_pres.append(high_pressure[j])
            new_time.append(dataset.time[j])
    #    else:
    #        new_res.append(np.nan)
    #        new_pres.append(np.nan)
    #        new_time.append(dataset.time[j])

    plt.plot(new_time, new_pres, 'b',
             new_time, new_res, 'c',
             )

    plt.grid()
    plt.legend(["Pressure across orifice",
                "resistance",
                ])
    plt.show()

    mode = 'valid'
    filter_length = [250, 500, 1000, 2000, 3000]
    for m in filter_length:
        plt.plot(np.convolve(new_res, np.ones((m,))/m, mode=mode));

    plt.legend(filter_length, loc='lower center')
    plt.show()

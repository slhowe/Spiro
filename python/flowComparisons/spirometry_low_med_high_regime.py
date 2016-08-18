#!/bin/bash

# Import my extensions
import sys
sys.path.insert(0, '/home/sarah/Documents/Spirometry/python/extensions')

import csv
import matplotlib.pyplot as plt
from data_struct import DataStore as Data
from numpy import array, sign
from numpy.linalg import lstsq

path = '/home/sarah/Documents/Spirometry/data/'
plot_path = '/home/sarah/Documents/Spirometry/lungPlots/'
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

data = dataset_1
offset = 26
start = dataset_2_start_indices[0]
end = dataset_2_end_indices[0]
early_start = 1500
late_end = 1500
high_length = 14
low_length = 250

inh_start = -180
inh_end = -100

# measured
pressure = data.pressure[start-early_start : end+late_end]
flow = data.flow[start-early_start+offset : end+offset+late_end]

# inhalation
med_pressure = data.pressure[start+inh_start : start+inh_end]
med_flow = data.flow[start+inh_start+offset : start+inh_end+offset]

# turbulent
high_pressure = data.pressure[start : start+high_length]
high_flow = data.flow[start+offset : start+offset+high_length]

# laminar
low_pressure = data.pressure[end-low_length : end]
low_flow = data.flow[end+offset-low_length : end+offset]

# find Res
low_res = lstsq(array([low_flow]).T, array([low_pressure]).T)
LR = low_res[0][0][0]
print("LR: {}, residuals: {}".format(LR, low_res[1]))
high_res = lstsq(array([high_flow]).T, array([high_pressure]).T)
HR = high_res[0][0][0]
print("HR: {}, residuals: {}".format(HR, high_res[1]))
med_res = lstsq(array([med_flow]).T, array([med_pressure]).T)
MR = med_res[0][0][0]
print("MR: {}, residuals: {}".format(MR, med_res[1]))

# remake pressure
full_range = range(len(data.pressure)-offset)
mod_high_pressure = [HR*f for f in data.flow[offset::]]
mod_low_pressure = [LR*f for f in data.flow[offset::]]
mod_med_pressure = [MR*f for f in data.flow[offset::]]

plt.plot(full_range, data.pressure[0:-offset], 'b',
         full_range, mod_high_pressure, 'g',
         full_range, mod_low_pressure, 'r',
         full_range, mod_med_pressure, 'm',
         )

plt.grid()
plt.legend(["measured",
            "P = HR.Q",
            "P = LR.Q",
            "P = MR.Q",
            ])
plt.show()

#!/bin/bash

# Import my extensions
import sys
sys.path.insert(0, '/home/sarah/Documents/Spirometry/python/extensions')

import csv
import matplotlib.pyplot as plt
from data_struct import Data
from numpy import array
from numpy.linalg import lstsq

path = '/home/sarah/Documents/Spirometry/data/'
plot_path = '/home/sarah/Documents/Spirometry/lungPlots/'
files = 'DIFFUSTIK.csv'

# Create data classes
reset = Data('reset')
reset.clear_dataset_registry()
data = Data('dataset_1')

# Store pressure and flow data
file_index = 0
filename = path + files
with open(filename, 'rb') as csvfile:
    reader = csv.reader(csvfile, delimiter = ',')

    # skip header line
    header = reader.next()

    for row in reader:
        pressure = float(row[0])
        flow = float(row[1])

        data.pressure.append(pressure)
        data.flow.append(-1.0*flow)
file_index += 1

# pressure looks a bit wrong.
# shift down so 0 at 0 flow
data.pressure = [p-0.018 for p in data.pressure]

offset = +24
start = 0
end = len(data.pressure) - offset

med_pres_start = 0
low_pres_start = 4600
hf_low_start = 16520
hf_high_start = 17480
high_pres_start = 18300

# measured
pressure = data.pressure[start : end]
flow = data.flow[start+offset : end+offset]
sqr_flow = [f**2 for f in flow]

med_pres= data.pressure[med_pres_start: low_pres_start]
med_flow = data.flow[med_pres_start+offset : low_pres_start+offset]
med_sqr_flow = [f**2 for f in med_flow]

low_pres = data.pressure[low_pres_start : hf_low_start]
low_flow = data.flow[low_pres_start+offset : hf_low_start+offset]
low_sqr_flow = [f**2 for f in low_flow]

hf_low_pres = data.pressure[hf_low_start : hf_high_start]
hf_low_flow = data.flow[hf_low_start+offset : hf_high_start+offset]
hf_low_sqr_flow = [f**2 for f in hf_low_flow]

hf_high_pres = data.pressure[hf_high_start : high_pres_start]
hf_high_flow = data.flow[hf_high_start+offset : high_pres_start+offset]
hf_high_sqr_flow = [f**2 for f in hf_high_flow]

high_pres = data.pressure[high_pres_start : end]
high_flow = data.flow[high_pres_start+offset : end+offset]
high_sqr_flow = [f**2 for f in high_flow]


# find Res
med_res = lstsq(array([med_flow]).T, array([med_pres]).T)
MP = med_res[0][0][0]
print("MP: {}, residuals: {}".format(MP, med_res[1]))
low_res = lstsq(array([low_flow]).T, array([low_pres]).T)
LP = low_res[0][0][0]
print("LP: {}, residuals: {}".format(LP, low_res[1]))
hf_high_res = lstsq(array([hf_high_flow]).T, array([hf_high_pres]).T)
HF_HP = hf_high_res[0][0][0]
print("HF_HP: {}, residuals: {}".format(HF_HP, hf_high_res[1]))
hf_low_res = lstsq(array([hf_low_flow]).T, array([hf_low_pres]).T)
HF_LP = hf_low_res[0][0][0]
print("HF_LP: {}, residuals: {}".format(HF_LP, hf_low_res[1]))
high_res = lstsq(array([high_flow]).T, array([high_pres]).T)
HP = high_res[0][0][0]
print("HP: {}, residuals: {}".format(HP, high_res[1]))

# remake pressure
mod_med_pressure = [MP*f for f in flow]
mod_low_pressure = [LP*f for f in flow]
mod_high_pressure = [HP*f for f in flow]
mod_hf_low_pressure = [HF_LP*f for f in flow]
mod_hf_high_pressure = [HF_HP*f for f in flow]

full_range = range(len(pressure))
plt.plot(full_range, pressure, 'b',
         full_range, mod_med_pressure, 'g',
         full_range, mod_low_pressure, 'r',
         full_range, mod_hf_low_pressure, 'k',
         full_range, mod_hf_high_pressure, 'm',
         full_range, mod_high_pressure, 'c',
         )

plt.grid()
plt.legend(["measured",
            "P = MP.Q",
            "P = LP.Q",
            "P = HF_LP.Q",
            "P = HF_HP.Q",
            "P = HP.Q",
            ])
plt.show()

#!/bin/bash

# Import my extensions
import sys
sys.path.insert(0, '/home/sarah/Documents/Spirometry/python/extensions')

import csv
import matplotlib.pyplot as plt
from data_struct import DataStore as Data
from calculus import integral

path = '/home/sarah/Documents/Spirometry/python/'
files = [
        'data_recording.csv'
         ]

pressure_lowRes = []
pressure_highRes = []
time = []


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Extract pressure two data arrays (assumed pressure and flow,
# but could be pressure & pressure etc) and time stamp from csv
def extract_data_arrays_from_file(filename, pressure, flow, time):
    pressure = []
    flow = []
    time = []

    # Store pressure and flow data
    with open(filename, 'rb') as csvfile:
        reader = csv.reader(csvfile, delimiter = ',')

        # skip header line
        header = reader.next()

        for row in reader:
            pressure_point = float(row[0])
            flow_point = float(row[1])
            time_point = float(row[2])

            pressure.append(pressure_point)
            flow.append(flow_point)
            time.append(time_point)

    return(1)
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Calibration for flow sensor is in
# ~/Documents/Spirometry/python/calibration
# This may need to be re-checked
# NOTE Only valid up to 250Pa across spirometer
def convert_spirometer_pressure_to_flow(pressure):
    flow = [(-0.000475*p**2) + (0.427*p) for p in pressure]
    return flow
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Calibration for this conversion is in
# ~/Documents/Spirometry/python/maskLPF
def convert_spirometer_pressure_to_mouth_pressure(spir_pressure):
    pass
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


for i in range(len(files)):
    filename = path + files[i]
    extract_data_arrays_from_file(filename, pressure_lowRes, pressure_highRes, time)

#!/bin/bash

# Import my extensions
import sys
sys.path.insert(0, '/home/sarah/Documents/Spirometry/python/extensions')

import csv
import matplotlib.pyplot as plt
from data_struct import DataStore as Data
from calculus import integral

path = './trial_data/'
subject = [
          'F01',
          ]
files = [
        '_mp_noRx.csv',
        '_mp_Rx.csv',
         ]

# Coefficients for Ps -> Pm
# Pm = Ps*c[0] + c[1]
# [mp_RX, mp_woRX]
conversion = [
        [2.9, 0.05],
        [2.5, 0.03],
        ]



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Extract pressure two data arrays (assumed pressure and flow,
# but could be pressure & pressure etc) and time stamp from csv
# NOTE This function appends. Make sure empty arrays are passed in
def extract_data_arrays_from_file(filename, pressure, flow, time):
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
# ~/Documents/Spirometry/python/maskLPF/Ps_Pm_calibration
def convert_spirometer_pressure_to_mouth_pressure(spir_pressure, coeffs):
    mouth_pressure = [coeffs[0]*p + coeffs[1] for p in spir_pressure]
    return mouth_pressure
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


for j in range(len(subject)):
    for i in range(len(files)):
        pressure_lowRes = []
        pressure_highRes = []
        time = []

        filename = path + subject[j] + '/' + subject[j] + files[i]
        extract_data_arrays_from_file(filename, pressure_lowRes, pressure_highRes, time)

        flow = convert_spirometer_pressure_to_flow(pressure_highRes)
        mouth_pressure = convert_spirometer_pressure_to_mouth_pressure(pressure_highRes, conversion[i])

        plt.plot(time, pressure_lowRes, 'b', linewidth=2)
        plt.plot(time, pressure_highRes, 'k', linewidth=2)
        plt.xlabel('Time(s)')
        plt.ylabel('Pressure(Pa)')
        plt.legend(['low_res', 'high_res'])
        plt.show()

        plt.plot(time, flow, 'b', linewidth=2)
        plt.plot(time, pressure_highRes, 'm', linewidth=2)
        plt.plot(time, mouth_pressure, 'k', linewidth=2)
        plt.xlabel('Time(s)')
        plt.ylabel('Pressure(Pa)')
        plt.legend(['flow', 'high_res', 'mouth_pressure'])
        plt.show()

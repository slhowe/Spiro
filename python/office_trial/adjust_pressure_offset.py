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
          #'F01',
          #'F02',
          #'F03',
          #'F04',
          #'F05',
          #'F09',
          #'M01',
          'M02',
          'M03',
          'M04',
          'M05',
          'M06',
          'M07',
          'M08',
          'M09',
          'M10',
          ]
files = [
        #'_mask_noRx.csv',
        #'_mask_Rx.csv',
        #'_mask_peak.csv',
        #'_mask_peakRx.csv',
        #'_mp_noRx.csv',
        #'_mp_Rx.csv',
        '_mp_peak.csv',
        '_mp_peakRx.csv',
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
# Write two data arrays (assumed pressure and flow,
# but could be pressure & pressure etc) and time stamp from csv
def write_data_arrays_to_file(filename, pressure, flow, time):
    # Store pressure and flow data
    with open(filename, 'w') as csvfile:
        writer = csv.writer(csvfile, delimiter = ',')

        for data_point in range(len(pressure)):
            writer.writerow([pressure[data_point], flow[data_point], time[data_point]])
    return(1)
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~



for j in range(len(subject)):
    for i in range(len(files)):
        acceptable = False
        pressure_lowRes = []
        pressure_highRes = []
        time = []

        filename = path + subject[j] + '/' + subject[j] + files[i]
        try:
            extract_data_arrays_from_file(filename, pressure_lowRes, pressure_highRes, time)

            while not acceptable:
                print("Close plot and enter pressure offset")
                plt.plot(time, pressure_lowRes, 'b', linewidth=2)
                plt.xlabel('Time(s)')
                plt.ylabel('Pressure(Pa)')
                plt.grid()
                plt.show()
                lowRes_offset = input("Offset: ")

                print("Close plot and enter pressure offset")
                plt.plot(time, pressure_highRes, 'k', linewidth=2)
                plt.xlabel('Time(s)')
                plt.ylabel('Pressure(Pa)')
                plt.grid()
                plt.show()
                highRes_offset = input("Offset: ")

                new_pressure_lowRes = [p - lowRes_offset for p in pressure_lowRes]
                new_pressure_highRes = [p - highRes_offset for p in pressure_highRes]

                print("Close plot and decide if acceptable")
                plt.plot(time, new_pressure_lowRes, 'b', linewidth=2)
                plt.plot(time, new_pressure_highRes, 'k', linewidth=2)
                plt.xlabel('Time(s)')
                plt.ylabel('Pressure(Pa)')
                plt.legend(['lowRes', 'highRes'])
                plt.grid()
                plt.show()

                acceptable = input("Offset correct? (1 = yes, 0 = no) ")

            new_filename = path + subject[j] + '/' + 'normalised_' + subject[j] + files[i]
            write_data_arrays_to_file(new_filename, new_pressure_lowRes, new_pressure_highRes, time)

        except IOError:
            print("File {} does not exist".format(filename))

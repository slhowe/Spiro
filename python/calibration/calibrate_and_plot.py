#!/bin/bash

# Import my extensions
import sys
sys.path.insert(0, '/home/sarah/Documents/Spirometry/python/extensions')

from filters import semi_gauss_lp_filter
import csv
import matplotlib.pyplot as plt
from data_struct import DataStore as Data
from numpy import isnan, array
from numpy.linalg import lstsq

path = '/home/sarah/Documents/Spirometry/python/calibration/calibrationFiles/'
files = [
        '21_5_Lmin-1',
        '43_0_Lmin-1',
        '61_0_Lmin-1',
        '11_85_Lmin-1',
        '30_6_Lmin-1',
        '43_75_Lmin-1',
        '67_45_Lmin-1',
        '8_45_Lmin-1',
        '18_55_Lmin-1',
        '30_8_Lmin-1',
        '51_0_Lmin-1',
        '70_18_Lmin-1',
        '2_11_Lmin-1',
        '39_8_Lmin-1',
        '55_85_Lmin-1',
        '79_8_Lmin-1',
        '78_15_Lmin-1',
         ]

other_files = [
        '100_0_Lmin-1',
        '90_5_Lmin-1',
        '88_1_Lmin-1',
        #'79_8_Lmin-1',
        #'78_15_Lmin-1',
        ]

# Create data classes
reset = Data('reset')
reset.clear_dataset_registry()

def mean(array):
    sum_array = sum(array)
    len_array = len(array)
    return sum_array/len_array


flows = []
other_flows = []
pressures = []
other_pressures = []

for i in range(len(files)):
    # Store pressure and flow data
    dataset = Data('dataset_1')
    filename = path + files[i]

    # Get flow from filename
    split_filename = filename.split('_')
    split_split = split_filename[0].split('/')

    integer = split_split[-1]
    decimal = split_filename[1]

    flow = integer + '.' + decimal

    # read flow sensor data
    with open(filename, 'rb') as csvfile:
        reader = csv.reader(csvfile, delimiter = ',')

        for row in reader:
            data = float(row[0])
            dataset.flow.append(data*1000)

    # filter flow
    avg = mean(dataset.flow)

    flows.append(flow)
    pressures.append(avg)

    plotting = False
    if(plotting):
        Fs = 90
        time = [t/float(Fs) for t in range(len(dataset.flow))]

        plt.plot(
                time, dataset.flow, 'r',
                )
        plt.grid()
        plt.legend([
                    '{}'.format(flow),
                    ])
        plt.show()

saving = False
if(saving):
    with open('calibration_table.csv', 'w') as csvfile:
        csvwriter = csv.writer(csvfile)
        for i in range(len(pressures)):
            line = [pressures[i], flows[i]]
            csvwriter.writerow(line)

gp = []
gf = []
gf.append(100)
gf.append(90.5)
gf.append(88)
gp.append(340)
gp.append(280)
gp.append(320)

# Find an exponential fit to data
pressures_sqrd = [float(p)**2 for p in pressures]

dependent = array([flows])
independent = array([pressures_sqrd, pressures])
res = lstsq(independent.T, dependent.T)[0]

independent = array([pressures])
res2 = lstsq(independent.T, dependent.T)[0]

independent = array([pressures_sqrd])
res3 = lstsq(independent.T, dependent.T)[0]

for i in range(len(other_files)):
    # Store pressure and flow data
    dataset = Data('dataset_1')
    filename = path + other_files[i]

    # Get flow from filename
    split_filename = filename.split('_')
    split_split = split_filename[0].split('/')

    integer = split_split[-1]
    decimal = split_filename[1]

    flow = integer + '.' + decimal

    # read flow sensor data
    with open(filename, 'rb') as csvfile:
        reader = csv.reader(csvfile, delimiter = ',')

        for row in reader:
            data = float(row[0])
            dataset.flow.append(data*1000)

    # filter flow
    avg = mean(dataset.flow)

    other_flows.append(flow)
    other_pressures.append(avg)

    plotting = True
    if(plotting):
        Fs = 90
        time = [t/float(Fs) for t in range(len(dataset.flow))]

        plt.plot(
                time, dataset.flow, 'r',
                )
        plt.grid()
        plt.legend([
                    '{}'.format(flow),
                    ])
        plt.show()


# make line
print(res)
new_pressures = range(300)
new_pressures = [p for p in new_pressures]
line = [(p**2)*res[0][0] + p*res[1][0] for p in new_pressures]
line2 = [(p)*res2[0][0] for p in new_pressures]
line3 = [(p**2)*res3[0][0] for p in new_pressures]

plt.plot(pressures, flows, 'bo')
plt.plot(other_pressures, other_flows, 'ro')
plt.plot(gp, gf, 'yo')
plt.plot(new_pressures, line, 'g')
#plt.plot(new_pressures, line2, 'r')
#plt.plot(new_pressures, line3, 'k')
plt.legend([
     'data points used',
     'data points not used',
     'flow guess (top data cutoff)',
     'Q = aP^2 + bP',
#     'Q = cP',
#     'Q = dP^2'
    ])
plt.grid()
plt.show()

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

gp = [112, 32, 24, 9, 22, 26, 50, 54, 174, 203, 163, 105, 26, 116, 184, 48, 60]
gf = [45.9, 21.6, 7.8, 2.8, 11, 17.2, 28.2, 29.6, 59.5, 67.6, 54.9, 42.71, 18.6, 44.9, 63, 30.3, 32]

gp2 = [11, 12, 21, 42, 58, 59, 82, 13, 22, 50, 75, 83, 55, 36, 12, 16, 21, 29, 41, 61, 49, 72, 90, 87, 80, 45, 88, 32, 78, 80]
gf2 = [11.45, 24.2, 29.6, 39.2, 45.8, 51, 56, 16, 27.4, 39.5, 49.9, 56.2, 41.1, 34.2, 15.3, 22.7, 26.6, 31, 36.5, 43, 38.8, 47.7, 67.6, 63.7, 51.5, 37.3, 67, 35.5, 59.5, 78]

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
            data = float(row[1])
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
new_pressures = range(7000)
new_pressures = [p for p in new_pressures]
line = [(p**2)*res[0][0] + p*res[1][0] for p in new_pressures]
line2 = [(p)*res2[0][0] for p in new_pressures]
line3 = [(p**2)*res3[0][0] for p in new_pressures]

res_pts = [pressures[i]/float(flows[i])*0.6 for i in range(len(pressures))]

plt.plot(pressures, flows, 'bo')
plt.plot(other_pressures, other_flows, 'ro')
plt.plot(gp, gf, 'yo')
plt.plot(gp2, gf2, 'co')
plt.plot(new_pressures, line, 'g')
plt.plot(pressures, res_pts, 'rx')
#plt.plot(new_pressures, line2, 'r')
#plt.plot(new_pressures, line3, 'k')
plt.legend([
     'Flow sensor only calibration',
     'Data points not used',
     'Flow and Pressure sensor calibration points',
     'Flow and Pressure sensor calibration points 2',
     'Q = aP^2 + bP',
#     'Q = cP',
#     'Q = dP^2'
     'resistance points from data',
    ])
plt.grid()
plt.xlabel('Pressure Pa')
plt.ylabel('Flow L/min')
plt.show()

data = []
pressure = []
flow = []
for m in range(len(pressures)):
    # Pressures linked in tuples for later sorting
    data.append((pressures[m], flows[m]))
# Sorted in ascending order by mask pressure
# Because mask pressure will hit pressure cutoff first
data.sort(key=lambda tup: tup[1])
for i in range(len(data)):
    pressure.append(data[i][0])
    flow.append(data[i][1])

plt.plot(pressure, flow)
plt.show()

res = [pressure[m]/float(flow[m]) for m in range(len(pressures))]
plt.plot(flow, 'rx')
plt.plot(res)
plt.show()

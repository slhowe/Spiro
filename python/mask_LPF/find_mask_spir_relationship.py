#!/bin/bash

import csv
import matplotlib.pyplot as plt
from sympy import *
import numpy as np
from numpy import array
from numpy.linalg import lstsq

#############################################
# Map spirometer pressure and mask pressure #
# Find relationship between the two         #
#############################################
SENSOR_CUTOFF = 0.29

#filename = './data_recording.csv'
filename = './breath_example_Rx2.csv'

pressure = []

with open(filename, 'rb') as csvfile:
    reader = csv.reader(csvfile, delimiter = ',')

    # skip header line
    header = reader.next()

    for row in reader:
        # Values read as string. Need converting
        spir = float(row[0])
        mask = float(row[1])

        # Pressures linked in tuples for later sorting
        pressure.append((spir, mask))

# Sorted in ascending order by mask pressure
# Because mask pressure will hit pressure cutoff first
pressure.sort(key=lambda tup: tup[1])

spir_pressure = []
mask_pressure = []

# Sensor cuts off
for i in range(len(pressure)):
    spir_pressure.append(pressure[i][0])
    mask_pressure.append(pressure[i][1])

# relationship in expiration
# Find first index less than 0
# find first index less than -0.29
# relationship for all values in between
start = -1
end = -1
i = 0
num_datapoints = len(mask_pressure)
while i < num_datapoints:
    if(start < 0):
        if(mask_pressure[i] > -0.29):
            start = i
    else:
        if(mask_pressure[i] >= 0):
            end = i - 1
            i = num_datapoints
    i+=1

# Quadratic relationship?
one_array = [1]*(end-start)
spir_array = spir_pressure[start:end]
mask_array = mask_pressure[start:end]
spir_array_sqd = [x**2 for x in spir_array]

independents = array([spir_array_sqd, spir_array, one_array])
dependents = array(mask_array)

result = lstsq(independents.T, dependents.T)
print(result)

curve = range(0, -150, -1)
curve = [c/1000.0 for c in curve]
data = [result[0][0]*(c**2) + result[0][1]*c + result[0][2] for c in curve]
data2 = [-6.4*(c**2) + c for c in curve]

#linear relationship
start = -1
end = -1
i = 0
num_datapoints = len(mask_pressure)
while i < num_datapoints:
    if(start < 0):
        if(mask_pressure[i] > -0.20):
            start = i
    else:
        if(mask_pressure[i] >= -0.045):
            end = i - 1
            i = num_datapoints
    i+=1

one_array_l = [1]*(end-start)
spir_array_l = spir_pressure[start:end]
mask_array_l = mask_pressure[start:end]

independents = array([spir_array_l, one_array_l])
dependents = array(mask_array_l)

result = lstsq(independents.T, dependents.T)
print(result)

curve_l = range(-20, -150, -1)
curve_l = [c/1000.0 for c in curve_l]
data_l = [result[0][0]*c + result[0][1] for c in curve_l]
data_l2 = [1.8*c + 0.03 for c in curve_l]

# 1-1 line
curve_1t1 = range(-150, 260, 10)
curve_1t1 = [c/1000.0 for c in curve_1t1]

plt.plot(mask_pressure, spir_pressure, 'x')
plt.plot(mask_array, spir_array, 'mx')
plt.plot(data, curve, 'r')
plt.plot(data2, curve, 'g')
plt.plot(data_l, curve_l, 'k', linewidth=3)
plt.plot(data_l2, curve_l, ':k', linewidth=3)
plt.plot(curve_1t1, curve_1t1, 'k', linewidth=3)
plt.xlabel('mask')
plt.ylabel('spir')
plt.show()


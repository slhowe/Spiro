#!/bin/bash

# Import my extensions
import sys
sys.path.insert(0, '/home/sarah/Documents/Spirometry/python/extensions')

import numpy as np
from numpy import nan
import csv
import matplotlib.pyplot as plt
from data_struct import DataStore as Data
from calculus import integral
from numpy.linalg import lstsq
from numpy import array

path = './'
files = [
        'data_recording.csv'
         ]

# Create data classes
reset = Data('reset')
reset.clear_dataset_registry()

mask_pressure_array = []
spir_pressure_array = []
flow = []
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
            spir_pressure = -float(row[0])
            mask_pressure = -float(row[1])
            temp_flow = (spir_pressure*1000)**2*-0.000475 + (spir_pressure*1000)*0.427
            time = float(row[2])

            mask_pressure_array.append(mask_pressure)
            spir_pressure_array.append(spir_pressure)
            flow.append(temp_flow)
            dataset.time.append(time/1000)
            if(spir_pressure != 0):
                relat_array.append(mask_pressure/spir_pressure)
            else:
                relat_array.append(np.nan)

    volume = integral(flow, Fs)

#    plt.plot(dataset.time, dataset.pressure, 'b',
#             dataset.time, dataset.flow, 'r',
#             )
    plt.plot(dataset.time, spir_pressure_array, 'b',
             dataset.time, mask_pressure_array, 'g',
          #   dataset.time, flow, 'm',
             dataset.time, volume, 'c',
             #dataset.time, relat_array, 'r',
             )
    plt.show()

#    plt.plot(volume, mask_pressure_array)
#    plt.grid()
#    plt.legend(["spir_pressure",
#                "mask_pressure",
#                'flow',
#                'volume',
#                ])
#    plt.show()
#
#    def lazy_filt(x, y, step=1):
#        new_len = int(len(x)/step)
#        new_x = [0]*new_len
#        new_y = [0]*new_len
#        for i in range(new_len):
#            mean_value = lambda z, i, step: np.mean(z[i*step:i*step+step])
#            new_x[i] = mean_value(x, i, step)
#            new_y[i] = mean_value(y, i, step)
#        return(new_x, new_y)
#
#    def subsample(x, y, step=1):
#        new_len = int(len(x)/step)
#        new_x = [0]*new_len
#        new_y = [0]*new_len
#        for i in range(new_len):
#            new_x[i] = x[i*step+step/2]
#            new_y[i] = y[i*step+step/2]
#        return(new_x, new_y)
#
#
##    gaussian_func = lambda x, sigma: 1/np.sqrt(2*np.pi*sigma**2) * np.exp(-(x**2)/(2*sigma**2))
#    fig, ax = plt.subplots()
#    # Compute moving averages using different window sizes
#    x = dataset.time
#    y = mask_pressure_array
#    N = len(y)
#    sigma_lst = [1, 100, 200, 300]
# #   y_gau = np.zeros((len(sigma_lst), N))
#    for i, sigma in enumerate(sigma_lst):
#  #      gau_x = np.linspace(-2.7*sigma, 2.7*sigma, 6*sigma)
#  #      gau_mask = gaussian_func(gau_x, sigma)
#  #      y_gau[i, :] = np.convolve(y, gau_mask, 'same')
#        new_x, new_y = lazy_filt(x, y, sigma)
#        ax.plot(new_x, new_y, label=r"$steps = {}$".format(sigma))
#        new_x, new_y = subsample(x, y, sigma)
#        ax.plot(new_x, new_y, label=r"$subsample steps = {}$".format(sigma))
#  #      ax.plot(x, y_gau[i, :] + (i+1)*50, label=r"$\sigma = {}$, $points = {}$".format(sigma, len(gau_x)))
#    # Add legend to plot
#    ax.legend(loc='upper left')
#    plt.show()
#
#    new_time, new_pressure = lazy_filt(dataset.time, mask_pressure_array, 3)
#    new_time2, new_volume = lazy_filt(dataset.time, volume, 3)
#    delta_P = np.diff(new_pressure)
#    delta_V = np.diff(new_volume)
#
#    dependents = array([delta_P])
#    independents = array([delta_V])
#    result = lstsq(independents.T, dependents.T)
#    print(result)
#
#    x_axis = range(-8, 9, 1)
#    x_axis = [x/100.0 for x in x_axis]
#    line = [result[0][0][0]*V for V in x_axis]
#
#    plt.plot(new_pressure)
#    plt.plot(delta_P)
#    plt.plot(new_volume)
#    plt.plot(delta_V)
#    plt.show()
#
#    plt.plot(delta_V, delta_P, 'x')
#    plt.plot(x_axis, line, 'r')
#    plt.show()

    def find_peak_and_end(data, thresh=1.0):
        peaks = []
        ends = []
        count = 0
        start = nan
        for i in range(len(data)):
            if np.isnan(start):
                if(data[i] > thresh):
                    start = i
            else:
                if(data[i] < thresh):
                    peak_start = data[start:i-1].index(np.max(data[start:i-1]))
                    start += peak_start
                    peaks.append(start)
                    ends.append(start + ((i-1)-(start))*3/4)
                    start = nan
                    count += 1
        return(peaks, ends)
    p, e = find_peak_and_end(mask_pressure_array)
    print(p)
    print(e)

    linear_pressure = []
    linear_volume = []
    for n in range(len(p)):
        linear_pressure += (mask_pressure_array[p[n]:e[n]])
        linear_volume += (volume[p[n]:e[n]])

    P_lines = []
    V_lines = []
    for n in range(len(p)):
        delta_P = mask_pressure_array[p[n]:e[n]]
        delta_V = (volume[p[n]:e[n]])

        ones = [1]*len(delta_P)
        steps = range(len(delta_P))
        dependents = array([delta_P])
        independents = array([steps, ones])
        Presult = lstsq(independents.T, dependents.T)
        print(Presult)

        line = [Presult[0][0][0]*m + Presult[0][1][0] for m in range(len(delta_P))]
        P_lines+= line

        dependents = array([delta_V])
        Vresult = lstsq(independents.T, dependents.T)
        print(Vresult)

        line = [Vresult[0][0][0]*m + Vresult[0][1][0] for m in range(len(delta_P))]
        V_lines+= line

        print("~~~")
        print(Presult[0][0][0]/Vresult[0][0][0])


    plt.plot(mask_pressure_array)
    ps = [mask_pressure_array[x] for x in p]
    es = [mask_pressure_array[x] for x in e]
    plt.plot(p, ps, 'd')
    plt.plot(e, es, 's')
    plt.grid()
    plt.show()

    plt.plot(linear_pressure)
    plt.plot(linear_volume)
    plt.plot(P_lines)
    plt.plot(V_lines)
    plt.show()



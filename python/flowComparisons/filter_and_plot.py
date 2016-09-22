#!/bin/bash

# Import my extensions
import sys
sys.path.insert(0, '/home/sarah/Documents/Spirometry/python/extensions')

import csv
import matplotlib.pyplot as plt
from filters import hamming

path = '/home/sarah/Documents/Spirometry/data/'
plot_path = '/home/sarah/Documents/Spirometry/images/exponential_fit'
files = ['Loops_1.csv', 'Loops_3.csv']

class Data:
    def __init__(self):
        self.pressure = []
        self.flow = []

# Store pressure and flow data
dataset = Data()
filename = path + files[0]
with open(filename, 'rb') as csvfile:
    reader = csv.reader(csvfile, delimiter = ',')

    # skip header line
    header = reader.next()

    for row in reader:
        pressure = float(row[0])
        flow = float(row[1])

        dataset.pressure.append(pressure)
        dataset.flow.append(-1.0*flow)

Fs = 190
Fc = 10
bw = 11
length = len(dataset.pressure)
filt = hamming(dataset.pressure, Fc, Fs, bw, plot=True)

#~~~~~~ plot example ~~~~~~
plotting = False
if(plotting):
    times = [x/Fs for x in range(length)]

    plt.plot(dataset.pressure, 'b')
    plt.plot(filt, 'r')
    plt.plot(dataset.flow, 'g')

    plt.legend(['Input pressure',
                'Output pressure',
                'Modelled input',
                'Modelled input no RC'])
    plt.grid()

    saving = False
    if saving:
        print(fig_name)
        plt.savefig(fig_name)

    plt.show()
    plt.close()

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
import scipy.signal as sig

# Create data classes
P_pneu = []
P_tot = []
flow = []

plotting = True

# Store pressure and flow data
filename = './data_recording.csv'
with open(filename, 'rb') as csvfile:
    reader = csv.reader(csvfile, delimiter = ',')

    # skip header line
    header = reader.next()

    for row in reader:
        pressure_tot = float(row[0])
        flow_point = float(row[1])
        pressure_pneu = float(row[2])

        P_pneu.append(pressure_pneu)
        P_tot.append(pressure_tot)
        flow.append(flow_point)

datalength = len(flow)

#Pneumotach resistance across all data
R_pneu = [0 for k in range(datalength)]
for j in range(datalength):
    if(flow[j] != 0):
        R_pneu[j] = P_pneu[j]/flow[j]

#Total resistance across all data
R_tot = [0 for k in range(datalength)]
for j in range(datalength):
    if(flow[j] != 0):
        R_tot[j] = P_tot[j]/flow[j]

#Added resistance across all data
R_add= [0 for k in range(datalength)]
for j in range(datalength):
    if(flow[j] != 0):
        P_add = P_tot[j] - P_pneu[j]
        R_add[j] = P_add/flow[j]

# Medain filter the resistances
R_pneu_filt = sig.medfilt(R_pneu, 221)
R_tot_filt = sig.medfilt(R_tot, 221)
R_add_filt = sig.medfilt(R_add, 221)

#Plot
if(plotting):
    Fs = 90
    time = [t/float(Fs) for t in range(datalength)]

    plt.plot(
            time, R_pneu, 'xr',
            time, R_tot, 'ob',
            time, R_add, '+k',
            time, R_pneu_filt, 'r',
            time, R_tot_filt, 'b',
            time, R_add_filt, 'k',
            )
    plt.grid()
    plt.legend(['pneumotach', 'total', 'added'])
    plt.show()


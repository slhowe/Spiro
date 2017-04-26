#!/bin/bash

'''
Analysis of data from simulation.
Simulated output is for cpap.
'''

# Add my extensions to path
import sys
sys.path.insert(0, '/home/sarah/Documents/Spirometry/python/extensions')

# Import my extensions
from breath_analysis import BreathData, split_breaths, calc_flow_delay
from calculus import integral, derivative

#Import built-ins
import matplotlib.pyplot as plt
from math import sin, cos, pi, exp, log
from numpy import nan
import numpy as np
from numpy.linalg import lstsq
import csv

if(__name__ == '__main__'):
    filename = [
            'dataset_1.csv'
            ]
    file_location = '/home/sarah/Documents/Spirometry/python/qucs_cpap_analysis/'
    file_to_extract = file_location + filename[0]
    fs = 100

    time_data = []
    data = []

    # extract time and simulated data from csv
    with open(file_to_extract, 'r') as csvfile:
        reader = csv.reader(csvfile, delimiter=';')
        header = reader.next()
        for row in reader:
            time_data.append(float(row[0]))
            data.append(float(row[1]))

#    filename = [
#            'dataset_4_no_pump.csv'
#            ]
#    file_to_extract = file_location + filename[0]
#
#    data_no_pump = []
#
#    # extract time and simulated data from csv
#    with open(file_to_extract, 'r') as csvfile:
#        reader = csv.reader(csvfile, delimiter=';')
#        header = reader.next()
#        for row in reader:
#            data_no_pump.append(float(row[1]))
#

    # Things we know (see README)
    R_spir = 1
    R_pump = 9
    P_pump = 10

    P_pm = P_pump*(R_spir/float(R_spir+R_pump))
    Ip = P_pm/R_spir
    Q_pump = Ip

    # Split breaths
    breaths = [[data[0]]]
    times = [[time_data[0]]]
    for j in range(1, len(data)):
        if(data[j] >= Ip*R_spir and data[j-1] < Ip*R_spir ):
            breaths.append([data[j]])
            times.append([time_data[j]])
        else:
            breaths[-1].append(data[j])
            times[-1].append(time_data[j])

    # Every breath is the same so use the second one
    breath = breaths[2]
    time = times[2]

    # Find the point where breathing in stops
    half_breath = 0
    for j in range(len(breath)):
        if(breath[j] > Ip*R_spir and breath[j+1] < Ip*R_spir):
            half_breath = j
            j = len(breath)

    # Flow for breath
    Q_combined = [b/float(R_spir) for b in breath]
    Q_lung_only = [q-Q_pump for q in Q_combined]
    P_lung_only = [q*R_spir for q in Q_lung_only]
    Q_pump = [Q_pump]*len(breath)

    P_pm = [P_pump]*len(breath)

    # Volume
    V_combined = integral(Q_combined, fs)
    V_combined[0] = 0.001
    V_lung_only = integral(Q_lung_only, fs)
    V_lung_only[0] = 0.001
    V_pump = integral(Q_pump, fs)
    V_pump[0] = 0.001

    # E/R
    EoR_combined = [Q_combined[j]/V_combined[j] for j in range(len(Q_combined))]
    EoR_lung = [Q_lung_only[j]/V_lung_only[j] for j in range(len(Q_combined))]
    EoR_lung = [Q_lung_only[j]/V_combined[j] for j in range(len(Q_combined))]
    EoR_pump= [Q_pump[j]/(j+1) for j in range(len(Q_combined))]
    EoR_diff = [EoR_combined[j] - EoR_lung[j] for j in range(len(EoR_combined))]
    #EoR_diff = [EoR_combined[j] - EoR_pump[j] for j in range(len(EoR_combined))]

    thing = [EoR_combined[j]/EoR_diff[j] for j in range(len(EoR_diff))]

    EoR_pump_guess = [0.52*g for g in EoR_combined]

    iEoR_combined = integral(EoR_combined, 100)
    iEoR_lung = integral(EoR_lung, 100)


    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    plt.plot(EoR_combined[1:half_breath], 'k')
    plt.plot(EoR_lung[1:half_breath], 'rx')
    plt.plot(EoR_pump[1:half_breath], 'ro')
    plt.plot(EoR_pump_guess[1:half_breath], 'mo')
    plt.plot(EoR_diff[1:half_breath], 'b')
    plt.plot(thing[1:], 'g')


    plt.plot(Q_combined[1:], 'k')
    #plt.plot(V_combined[1:half_breath], 'k')
    plt.plot(Q_lung_only[1:], 'r')
    #plt.plot(V_lung_only[1:half_breath], 'r')
    plt.plot(Q_pump[1:], 'c')
    #plt.plot(V_pump[1:half_breath], 'm')
    plt.grid()
    plt.show()

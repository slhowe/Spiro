#!/bin/bash

# Add my extensions to path
import sys
sys.path.insert(0, '/home/sarah/Documents/Spirometry/python/extensions')

# Import my extensions
from data_struct import DataStore
from breath_analysis import BreathData, split_breaths, calc_flow_delay
from calculus import integral, derivative
from filters import hamming
from resistance_plots import ResistanceTable

#Import built-ins
import numpy as np
import matplotlib.pyplot as plt
from math import sin, cos, pi, exp, log
from numpy import isnan, array, real, nan
from numpy.linalg import lstsq

# Data path
#path = '/home/sarah/Documents/Spirometry/data/'
path = './'

#files = ['Normal_1.csv']
#files = ['example_breaths']
#files = ['example_breaths_size']
#files = ['example_breaths_tidal']
#files = ['Normal_1.csv',
#        'Banding_1.csv',
#        #'Normal_3.csv',
#        #'Banding_3.csv',
#        'Inflated_1.csv',
#        ]
files = ['mixed_results.csv',
        ]
#files = ['Inflated_3.csv']
#files = ['Loops_3.csv']
#files = [
#        'tidal_noR_1',
#        'tidal_lowR_1',
#        'tidal_noR_2',
#        'tidal_lowR_2',
#        'tidal_highR',
#        ]

files = [path + name for name in files]

# Place to save plots
plot_path = '/home/sarah/Documents/Spirometry/lungPlots/'

#f, (ax1, ax2, ax3) = plt.subplots(3, sharex=True)
f, (ax2) = plt.subplots(1, sharex=True)
for filename in files:
    print('\n...New File...')
    ''' Stepping to estimate RC

        (pg 15 of my little red book)

        Plots are saved if plot_path given
        Data in kPa and L and s
    '''

    # Definitions
    INSP = 0
    EXP = 1
    sampling_frequency = 85
    #sampling_frequency = 125

    # Remove any old data sets hanging around
    # (Needed when using ipython)
    reset = DataStore()
    reset.clear_dataset_registry()

    # Create data storage
    # Store pressure and flow data from files
    dataset = DataStore(Fs=sampling_frequency)
    dataset.store_data(filename)

    # Filter data
    # Use hamming filter
    # data, Fc, Fs, bw
    flow = dataset.flow
    pressure = dataset.pressure
    #pressure = [10*p for p in dataset.pressure]

    # Split flow into breaths
    # Find start, middle and end of each breath
    flow_splits = split_breaths(dataset.flow, plot=False)
    flow_starts = flow_splits[0]
    flow_middles = flow_splits[1]
    flow_stops = flow_splits[2]

    # Set up storage for breath data
    # Corrects for pressure delay
    # Not using pressure here so no delay needed
    data = BreathData(0)

    ##################################################
    print('number of full breaths: {}'.format(len(flow_starts)))
    number_of_breaths = len(flow_starts)-1
    number_of_breaths = 25
    ##################################################

    RCi = []
    RCi_norm = []
    RCe = []
    RCe_norm = []
    resistance = []

    dependent = np.array([pressure])
    independent = np.array([flow])
    spir_resistance = lstsq(independent.T, dependent.T)
    resistance.append(-spir_resistance[0][0][0])

    # Iterate through breaths
    for breath in range(15, number_of_breaths):
        # Extract breath data
        data.get_data(pressure,
                      flow,
                      flow_starts[breath],
                      flow_middles[breath],
                      flow_stops[breath],
                      dataset.sampling_frequency
                      )

        dependent = np.array([data.insp_pressure + data.exp_pressure])
        independent = np.array([data.insp_flow + data.exp_flow])
        spir_resistance = lstsq(independent.T, dependent.T)
        resistance.append(-spir_resistance[0][0][0])

        #plt.plot(data.insp_pressure)
        #plt.plot([f * spir_resistance[0][0][0] for f in data.insp_flow])
        #plt.show()

        # Go through each breath one half at a time
        # First: inspiration
        # Second: expiration
        for section in range(2):
            # Set data for insp/exp here
            if(section == INSP):
#                print('\nInsp parameters:')

                # Start at first real positive point
                start = 0
                i = data.insp_length/2
                while i >= 0:
                    if data.insp_flow[i] <= 0:
                        start = i + 1
                        i = 0
                    i -= 1

                end = min(data.insp_length-1, sampling_frequency*10)
                breath_length = end - start

                pres = [p for p in data.insp_pressure[start:end]]
                flw = data.insp_flow[start:end]
                vol = integral(flw, sampling_frequency)
                #time = data.time[start:end]

            elif(section == EXP):
#                print('\nExp parameters:')

                # Start at first real negative point
                start = 0
                i = 0
                while i < data.exp_length:
                    if data.exp_flow[i] < 0:
                        start = i
                        i = data.exp_length
                    i += 1

                end = min(data.exp_length, sampling_frequency*10)
                breath_length = end - start

                pres = [p for p in data.exp_pressure[start:end]]
                flw = [f for f in data.exp_flow[start:end]]
                vol = integral(flw, sampling_frequency)
                #time = data.time[-breath_length + start:end]

            for v in range(breath_length):
                if vol[v] == 0:
                    vol[v] = 1e-2

            curve = [flw[i]/vol[i] for i in range(breath_length)]

            for c in range(breath_length):
                if curve[c] == 0:
                    curve[c] = 1e-2

            int_curve = integral(curve, sampling_frequency)
            if(int_curve[-1] < 0):
                int_curve[-1] = 0
            t = len(int_curve)/float(sampling_frequency)

#            print('integral: {}'.format(int_curve[-1]))
#            print('time: {}'.format(len(int_curve)/float(sampling_frequency)))
#            print('integral norm: {}'.format(int_curve[-1]/t))
#            print('RC[-1] {}'.format(int_curve[-1]))

            if section == INSP:
                RCi.append(int_curve[-1])
                RCi_norm.append(int_curve[-1]/t)
            else:
                RCe.append(int_curve[-1])
                RCe_norm.append(int_curve[-1]/t)

            # frequency of fundamental sinusoidal part
            w = np.pi*sampling_frequency/len(flw)

            if(0):
                plt.plot(pres)
                plt.plot(flw)
                plt.plot(vol)
                plt.grid()
                plt.show()

            if(0):
                ''' Plot after changing data offsets'''
#                if(section == INSP):
#                    # plot V, Q, RC without any adjustment
#                    f, (ax1, ax2, ax3, ax4) = plt.subplots(4, sharex=True)
#                    RC = [0]*(breath_length)
#                    i = 0
#                    while i < (breath_length):
#                        v = vol[i]
#                        q = flw[i]
#                        if v == 0:
#                            v = 1e-2
#                        RC[i] = q/v
#                        i += 1
#
#                    time = [t/float(sampling_frequency) for t in range(breath_length)]
#                    wt = [t*w for t in time]
#                    ideal_flow = [np.sin(x) for x in wt]
#                    mxflow = max(flw)
#                    ideal_flow_scaled = [mxflow*np.sin(x) for x in wt]
#                    ideal_volume = [-np.cos(x) + 1 for x in wt]
#                    ideal_volume_scaled = [vol[-1]/2.0*(-np.cos(x) + 1) for x in wt]
#                    ideal_QV = [ideal_flow[j]/ideal_volume[j] for j in range(len(ideal_flow))]
#                    leftover_QV = [RC[j] / ideal_QV[j] for j in range(len(RC))]
#
#                    ax1.plot(time, flw, linewidth=3)
#                    ax1.plot(time, ideal_flow, linewidth=3)
#                    ax1.plot(time, ideal_flow_scaled, ':', linewidth=3)
#                    ax2.plot(time, vol, linewidth=3)
#                    ax2.plot(time, ideal_volume, linewidth=3)
#                    ax2.plot(time, ideal_volume_scaled, ':', linewidth=3)
#                    ax3.plot(time, RC, linewidth=3)
#                    ax3.plot(time, ideal_QV, linewidth=3)
#                    ax4.plot(time, leftover_QV, 'r', linewidth=3)
#
#                    ax1.set_title('Edrs in inspiration (Spirometry)', fontsize=18)
#                    ax1.set_ylabel('Flow (L/s)', fontsize=18)
#                    ax2.set_ylabel('Volume (L)', fontsize=18)
#                    ax3.set_ylim([-20, 160])
#                    ax3.set_ylabel('E(t)/R', fontsize=18)
#                    ax4.set_ylabel('E(t)/R', fontsize=18)
#                    ax4.set_xlabel('time (s)', fontsize=18)
#
#                    ax3.legend([
#                        "Result",
#                        "Ideal for E/R = 1",
#                        ])
#
#                    ax1.grid()
#                    ax2.grid()
#                    ax3.grid()
#                    ax4.grid()
#                    plt.show()

                if(section == EXP):
                    f, (ax1, ax2, ax3) = plt.subplots(3, sharex=True)
                    RC = [0]*(breath_length)
                    i = 0
                    while i < (breath_length):
                        v = vol[i]
                        q = flw[i]
                        if v == 0:
                            v = 1e-2
                        RC[i] = q/v
                        i += 1

                    time = [t/float(sampling_frequency) for t in range(breath_length)]
                    ax1.plot(time, flw, linewidth=3)
                    ax2.plot(time, vol, linewidth=3)
                    ax3.plot(time, RC, linewidth=3)
                    ax1.set_title('Edrs in exhalation (Spirometry)', fontsize=18)
                    ax1.set_ylabel('Flow (L/s)', fontsize=18)
                    ax2.set_ylabel('Volume (L)', fontsize=18)
                    ax3.set_ylim([-20, 160])
                    ax3.set_ylabel('1/RC(t)', fontsize=18)
                    ax3.set_xlabel('time (s)', fontsize=18)

                    ax1.grid()
                    ax2.grid()
                    ax3.grid()
                    plt.show()

                if(0):
                    # Say R = 2
                    f, (ax4) = plt.subplots(1, sharex=True)
                    R = 2
                    E = [rc*R for rc in RC]
                    P = [E[i]*vol[i] for i in range(len(vol))]
                    ax4.plot(time, P)
                    ax4.grid()
                    plt.show()


    if(1):
        #f, (ax1, ax2, ax3) = plt.subplots(3, sharex=True)

        print('min E/Ri: {}'.format(min(RCi)))
        print('min E/Re: {}'.format(min(RCe)))
        print('max E/Ri: {}'.format(max(RCi)))
        print('max E/Re: {}'.format(max(RCe)))
        print('E/Ri: {}'.format(np.percentile(RCi, [10, 25, 50, 75, 90])))
        print('E/Re: {}'.format(np.percentile(RCe, [10, 25, 50, 75, 90])))
        print('Resistance: {}'.format(np.percentile(resistance, [10, 25, 50, 75, 90])))
        full_res = resistance[0]
        resistance = resistance[1:]

        #ax1.plot(flow[:flow_stops[number_of_breaths-1]])
        #ax1.legend(['Flow'])
        #ax1.grid()

        inh_align = [flow_starts[i]+(flow_middles[i]-flow_starts[i])/2.0 for i in range(number_of_breaths)]
        exh_align = [flow_middles[i]+(flow_stops[i]-flow_middles[i])/2.0 for i in range(number_of_breaths)]
        mean_RCi = np.median(RCi)

        #ax2.plot(inh_align, RCi, 'o')
        ax2.plot(RCi, 'o')
        ax2.plot([mean_RCi]*len(RCi))
        #ax2.plot(exh_align, RCe, 'd')
        #ax2.plot(RCe, 'd')
        ax2.legend([
                'Normal Breathing',
                'Normal Breathing mean',
                'Banded breathing',
                'Banded breathing mean',
                 'integral(E/R) insp1',
                 'integral(E/R) exp1',
                 'integral(E/R) insp2',
                 'integral(E/R) exp2',
                 'integral(E/R) insp3',
                 'integral(E/R) exp3',
                 'integral(E/R) insp4',
                 'integral(E/R) exp4',
                 'integral(E/R) insp5',
                 'integral(E/R) exp5',
                 ])
        ax2.set_ylabel('Measure of effort (E/R)')
        ax2.set_xlabel('Breath')

        #ax3.plot(flow_middles[:number_of_breaths], resistance, 'd')
        #ax3.plot(0, full_res, 'd')
        #ax3.set_ylabel('Resistance of spirometer')
        #ax3.grid()
        #plt.show()


    # Final clean up
    dataset.clear_dataset_registry()
#ax1.grid()
ax2.grid()
#ax3.grid()
plt.show()



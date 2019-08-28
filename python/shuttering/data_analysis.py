#!/bin/bash

# Import my extensions
import sys
sys.path.insert(0, '/home/sarah/Documents/Spirometry/python/extensions')

import csv
import matplotlib.pyplot as plt
from calculus import integral
from breath_analysis import split_breaths
from scipy import stats
import numpy as np
from filters import hamming
from numpy.linalg import lstsq
from numpy import nan
from math import sin, cos, pi, exp, log

path = './'
files = [
        'human.csv',
        'human_normal_nopuff.csv',
        'human_puff.csv',
         ]


FS = 300

def quick_breath_split(pressure, cutoff, RISING_EDGE=False, timeout=0):
    starts = []
    for index in range(len(pressure)-1):
        if RISING_EDGE:
            if (pressure[index] < cutoff and pressure[index+1] >= cutoff):
                if (len(starts) != 0 and starts[-1] > index-timeout):
                    pass
                else:
                    starts.append(index)
        else:
            if (pressure[index] > cutoff and pressure[index+1] <= cutoff):
                if (len(starts) != 0 and starts[-1] > index-timeout):
                    pass
                else:
                    starts.append(index)
    return(starts)

gradients= [[] for i in range(len(files))]
decays = [[] for i in range(len(files))]

for i in range(len(files)):
    pressure_abs= []
    flow = []
    time = []
    pressure_diff = []
    # Create data classes
    #####################
    # GET DATA FROM FILE#
    #####################
    # Store pressure and flow data
    filename = path + files[i]
    print(' ')
    print(filename)
    with open(filename, 'rb') as csvfile:
        reader = csv.reader(csvfile, delimiter = ',')

        # skip header line
        header = reader.next()

        for row in reader:
            dp_pressure_diff = float(row[0]) # High res, differential pressure (for solving flow)
            #dp_flow = float(row[1])
            dp_pressure_abs= float(row[2]) # low res (abs pressure measured)
            dp_time = float(row[3]) # in ms

            pressure_diff.append(dp_pressure_diff)
            pressure_abs.append(dp_pressure_abs)
            time.append(dp_time/float(1000))

    # Filter pressure
    # Filter high resolution pressure data
    fc = 100#filter cutoff
    bw = 30#filter bandwidth (taps?)
    pressure = hamming(pressure_diff, fc, 300, bw, plot=False)
    pressure = np.real(pressure).tolist()


    # Flow calculation incorrect on arduino end
    # Fixed now but sets 1-4 are incorrect so
    # calculating new flow here instead :)
    for p in pressure:
        if p > 0:
            flow.append(0.0071*p)
        else:
            flow.append(0.0076*p)

    volume = integral(flow, FS)

    #######################
    # Split breaths apart #
    #######################
    flow_splits = quick_breath_split(flow, 0.05, timeout=400, RISING_EDGE=True)
    flow_splits.append(len(flow)-1)

    maxFlow = max(flow[:flow_splits[2]]) # max of first 3 breaths which should be tidal
    print(maxFlow)
    shutterTrigger = maxFlow*1.2

    shutter_splits = quick_breath_split(flow, shutterTrigger, timeout=600, RISING_EDGE=True)
    shutter_splits.append(len(flow)-1)

    if(1):
        plt.rc('legend',**{'fontsize':12})
        f, (ax1, ax2, ax3, ax4) = plt.subplots(4, sharex=True)

        ax1.plot(pressure_abs)
        ax1.plot(flow_splits, [pressure_abs[m] for m in flow_splits], 'or')
        ax1.legend(['Abs pressure'])

        ax2.plot(pressure_diff)
        ax2.plot(flow_splits, [pressure_diff[m] for m in flow_splits], 'or')
        ax2.legend(['Diff pressure'])

        ax3.plot(volume)
        ax3.plot(flow_splits, [volume[m] for m in flow_splits], 'or')
        ax3.legend(['Volume'])

        ax4.plot(flow)
        ax4.plot(flow_splits, [flow[m] for m in flow_splits], 'or')
        ax4.plot(shutter_splits, [flow[m] for m in shutter_splits], 'om')
        ax4.legend(['Flow'])
        plt.show()

        plt.plot(flow, linewidth=2)
        plt.ylabel('Flow (L/s)')
        plt.show()


    # Setup for analysis
    averageFlow = [0]*(FS*3)
    numberAveraged = 0
    nextShutter = 0
    ANALYSING = False

    for breath in range(len(flow_splits)-1):
        print('\nbreath: {}'.format(breath))

        #######################
        # Split breaths apart #
        #######################
        flw = flow[flow_splits[breath]:flow_splits[breath+1]]
        pres_abs = pressure_abs[flow_splits[breath]:flow_splits[breath+1]]
        pres_diff = pressure_diff[flow_splits[breath]:flow_splits[breath+1]]
        pres_filt = pressure[flow_splits[breath]:flow_splits[breath+1]]
        vol = integral(flw, FS)
        t = time[flow_splits[breath]:flow_splits[breath+1]]


        if(shutter_splits[nextShutter] > flow_splits[breath] and shutter_splits[nextShutter] < flow_splits[breath+1]):
            ANALYSING = True
            print('This breath is shuttered')
            nextShutter += 1

        if not ANALYSING:
            if(max(vol)>0.3): # get rid of inspiration-only breaths
                maxFlw = max(flw)
                maxFlwIndex = flw.index(maxFlw)
                newX = [vol[j]-vol[maxFlwIndex] for j in range(len(flw))]
                plt.plot(vol, flw)

                # Calculate the average breath profile
                numberAveraged += 1
                for k in range(min(len(flw), len(averageFlow))):
                    averageFlow[k] += flw[k]

        #####################
        # Down-stroke phase #
        #####################
        if(ANALYSING):
            # Setup
            ANALYSING=False
            maxFlw = max(flw)
            maxFlwIndex = flw.index(maxFlw)
            new_start = maxFlwIndex + 10
            decayLen = min(55, len(flw)-new_start)
            new_end = new_start + decayLen

            # Calc average flow
            shutterFlow = [nan]*len(averageFlow)
            shutterVol = [nan]*len(averageFlow)
            if numberAveraged:
                averageFlow = [f/numberAveraged for f in averageFlow]
                print('num averaged: {}'.format(numberAveraged))
                averageVol = integral(averageFlow, FS)

                # Estimate flow and volume only caused by shuttering (from average)
                shutterFlow = [flw[k]-averageFlow[k] for k in range(min(len(averageVol), len(flw)))]
                shutterVol = integral(shutterFlow, FS)


            # If there is some data past peak flow
            if (new_end):
                # Get flow and volume after shutter
                volume_curve = shutterVol[new_start:new_end]
                flow_curve = shutterFlow[new_start:new_end]

                # Directly measure gradient
                gradient_passive, intercept, r_value, p_value, std_err = stats.linregress(volume_curve, flow_curve)
                gradients[i].append([gradient_passive])
                print('direct gradient: {}'.format(gradient_passive))
                print('r value: {}'.format(r_value))

                if(abs(r_value) < 0.85):
                    print('R value too low. IGNORING DATA...')
                    gradients[i][-1] = [np.nan]

                # remake lines
                line = [intercept + gradient_passive*m for m in [shutterVol[new_start], shutterVol[new_end]]]
                lineLong = [intercept + gradient_passive*m for m in [min(shutterVol), max(shutterVol)]]


                ###########################
                # Plots and print results #
                ###########################
                if(1):
                    plt.plot(vol, flw, 'orange', linewidth=3)
                    plt.plot(averageVol, averageFlow, 'red', linewidth=3)
                    plt.plot(shutterVol, shutterFlow, 'pink', linewidth=3)
                    plt.plot(vol[new_start], flw[new_start], 'oy')
                    plt.plot(vol[new_end], flw[new_end], 'om')
                    plt.plot([vol[new_start], vol[new_end]], line, 'g')
                    plt.xlabel('Vol')
                    plt.ylabel('Flow')
                    plt.grid()
                    plt.show()

                    plt.rc('legend',**{'fontsize':12})
                    f, (ax1, ax2) = plt.subplots(2, sharex=False)

                    ax1.plot(shutterVol, shutterFlow)
                    ax1.plot([shutterVol[new_start], shutterVol[new_end]], line, 'g')
                    ax1.plot([min(shutterVol), max(shutterVol)], lineLong, 'pink')
                    ax1.plot([shutterVol[new_start], shutterVol[new_end]], line, 'g')
                    ax1.plot(shutterVol[new_start:new_end], shutterFlow[new_start:new_end], 'k')
                    ax1.plot(shutterVol[new_start], shutterFlow[new_start], 'ro')
                    ax1.plot(shutterVol[new_end], shutterFlow[new_end], 'mo')
                    ax1.legend(['Shuttered data (Measured-average)',
                                'Linear fit to shuttered region'])

                    ax2.plot(flw, 'blue')
                    ax2.plot(averageFlow, 'g')
                    ax2.plot(shutterFlow, 'red')
                    ax2.plot(new_start, flw[new_start], 'oy')
                    ax2.plot(new_end, flw[new_end], 'oy')

                    ax1.grid()
                    ax2.grid()
                    ax2.legend(['Measured flow',
                                'Flow due to body',
                                'Flow due to shutter',
                                ])
                    ax1.set_ylabel('Flow [L/s]', fontsize=30)
                    ax1.set_xlabel('Volume [L]', fontsize=30)
                    ax2.set_ylabel('Flow [L/s]', fontsize=30)
                    ax2.set_xlabel('Time (shutter=150ms)', fontsize=30)

                    plt.show()


                # Calc exp decay of flow
                ln_flw = [np.log(abs(f)) for f in flw[new_start:new_end]]
                times = t[new_start:new_end]
                ones = [1]*len(times)
                dependent = np.array([ln_flw])
                independent = np.array([ones, times])
                try:
                    res = lstsq(independent.T, dependent.T)
                    decay = (res[0][1][0])
                    offset = np.exp(res[0][0][0])
                except(ValueError):
                    print('ValueError: Data has nan?')
                    decay = nan
                    offset = nan

                print('     Decay rate stuff')
                print(decay, offset,)


            else:
                print("no occluded section end found")
                gradients.append([np.nan])

            averageFlow = [0]*(len(averageFlow))
            numberAveraged = 0



    print('\n Gradients')
    for g in gradients[i]:
        print(g)
    print('\n Mean gradient [passive, driven]')
    print(np.nanmean(gradients[i], axis=0))
    print(np.nanstd(gradients[i], axis=0))

print('\nFinal results:')
for i in range(len(files)):
    print('\n Gradients')
    for g in gradients[i]:
        print(g)
    print('\n Mean gradient [passive, driven]')
    print(np.nanmean(gradients[i], axis=0))
    print(np.nanstd(gradients[i], axis=0))
    print('\n')


#    print('\n decays')
#    for g in decays[i]:
#        print(g)
#    print('\n Mean  [passive, driven]')
#    print(np.nanmean(decays[i], axis=0))
#    print(np.nanstd(decays[i], axis=0))

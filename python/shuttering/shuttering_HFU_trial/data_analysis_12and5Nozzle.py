#!/bin/bash

# Import my extensions
import sys
sys.path.insert(0, '/home/sarah/Documents/Spirometry/python/extensions')

import csv
import matplotlib.pyplot as plt
from calculus import integral, derivative
from breath_analysis import split_breaths
from scipy import stats
import numpy as np
from filters import hamming
from numpy.linalg import lstsq
from numpy import nan
from math import sin, cos, pi, exp, log
import seaborn as sns
import pandas as pd

#NOTE
#
# 1 - no Rx
# 2- Rx
# A - nom compliance
# B - banded compliance
#
#flow = 0
#volume = 1
#time = 2
#phase = 3
#shutter = 4
#pressure = 5

DATASET = 0

if(DATASET == 0):
    # Banded data set, has FEV in other sets too
    path = '/home/sarah/Documents/Spirometry/data/HFU_Respiratory_Trial/P1-JGC-18032019/'
    files = [
            'RespiratoryOcclusionMeasurement_19-03-18_02-46-59_Meas2_Shutter-1A.csv',
            'RespiratoryOcclusionMeasurement_19-03-18_03-18-00_Meas8_Shutter-1B.csv',
            'RespiratoryOcclusionMeasurement_19-03-18_03-11-37_Meas6_Shutter-2A.csv',
            'RespiratoryOcclusionMeasurement_19-03-18_02-58-17_Meas4_Shutter-2B.csv',
             ]
        # For data
    xLabels = ['Normal','Chest banded', 'Rx Normal', 'Rx banded']
    dataLoc = [
        [1, 4, 7, 8, 13, 2],
        [1, 4, 7, 8, 13, 2],
        [1, 7, 10, 11, 13, 2],
        [1, 4, 7, 8, 13, 2],
    ]

elif(DATASET == 1):
    # Trial to separate mechanics. Uses 12.5mm nozzle resistance
    # Has weird drops in flow after shuttering
    path = '/home/sarah/Documents/Spirometry/data/HFU_Respiratory_Trial/12andHalfNozzle/'
    files = [
            'Corrected_data_RespiratoryOcclusionMeasurement_19-04-17_02-31-43_Trial-1.csv',
            'Corrected_data_RespiratoryOcclusionMeasurement_19-04-17_02-33-25_Trial-1.csv',
            'Corrected_data_RespiratoryOcclusionMeasurement_19-04-17_02-35-38_Trial-1.csv',
            'Corrected_data_RespiratoryOcclusionMeasurement_19-04-17_02-37-26_Trial-1.csv',
             ]
        # For data
    xLabels = ['Test 1 - Rx=0','Test 2 - Rx=0','Test 3 - Rx=0.4', 'Test 4 - Rx=0.4']
    dataLoc = [
        [1, 4, 7, 8, 13, 2],
        [1, 4, 7, 8, 13, 2],
        [1, 4, 7, 8, 13, 2],
        [1, 4, 7, 8, 13, 2],
    ]
elif(DATASET == 2):
    # Trial testing multiple nozzle sizes. Only low works because of taper requirement (orifice plate sucks)
    path = '/home/sarah/Documents/Spirometry/data/HFU_Respiratory_Trial/P1-JGC-21032019 - Shuttering/'
    files = [
            #'Corrected_data_RespiratoryOcclusionMeasurement_19-03-21_02-20-42-High_Trial-1.csv',   # V. noisy
            #'Corrected_data_RespiratoryOcclusionMeasurement_19-03-21_02-25-52-Medium_Trial-1.csv', # V. noisy
            'Corrected_data_RespiratoryOcclusionMeasurement_19-03-21_02-32-59-Low_Trial-1.csv',
            'Corrected_data_RespiratoryOcclusionMeasurement_19-03-21_02-41-51-Low_Trial-2.csv',
            'Corrected_data_RespiratoryOcclusionMeasurement_19-03-21_02-36-54-No_Trial-1.csv',
             ]
        # For data
    xLabels = ['Test 1 - Rx=0.08','Test 2 - Rx=0.08','Test 3 - Rx=0']
    dataLoc = [
        #[1, 4, 7, 8, 13],
        #[1, 4, 7, 8, 13],
        [1, 4, 7, 8, 13, 2],
        [1, 7, 10, 11, 13, 2],
        [1, 7, 10, 11, 13, 2],
    ]
elif(DATASET == 3):
    # Test whether strange drops in flow after shutter is consistent
    path = '/home/sarah/Documents/Spirometry/data/HFU_Respiratory_Trial/bumpTest/09052019/'
    files = [
            #'Corrected_data_RespiratoryOcclusionMeasurement_19-05-09_02-27-01_Trial-1.csv',
            #'Corrected_data_RespiratoryOcclusionMeasurement_19-05-09_02-28-48_Trial-1.csv',
            'Corrected_data_RespiratoryOcclusionMeasurement_19-05-09_02-30-48_Trial-1.csv',
            'Corrected_data_RespiratoryOcclusionMeasurement_19-05-09_02-33-46_Trial-1.csv',
             ]
        # For data
    xLabels = ['Test 1 - Rx=0','Test 2 - Rx=0','Test 3 - Rx=0.4', 'Test 4 - Rx=0.4']
    dataLoc = [
        [1, 4, 7, 8, 13, 2],
        [1, 7, 10, 11, 13, 2],
        [1, 4, 7, 8, 13, 2],
        [1, 7, 10, 11, 13, 2],
    ]

elif(DATASET == 4):
    path = '/home/sarah/Documents/Spirometry/data/HFU_Respiratory_Trial/resistanceTrial/subject_01/'
    files = [
            'Corrected_data_RespiratoryOcclusionMeasurement_19-06-19_11-15-23_Trial-1_NO.csv',
            'Corrected_data_RespiratoryOcclusionMeasurement_19-06-19_11-19-05_Trial-1_12.5mm.csv',
            'Corrected_data_RespiratoryOcclusionMeasurement_19-06-19_11-22-47_Trial-1_10.5mm.csv',
            'Corrected_data_RespiratoryOcclusionMeasurement_19-06-19_11-26-04_Trial-1_9.5mm.csv',
            #'Corrected_data_SlowSpirometryMeasurement_19-06-19_11-11-24_Trial-1.csv',
             ]
        # For data
    xLabels = ['Normal','12.5mm','10.5mm', '9.5mm']
    dataLoc = [
        [1, 7, 10, 11, 13, 2],
        [1, 4, 7, 8, 13, 2],
        [1, 7, 10, 11, 13, 2],
        [1, 7, 10, 11, 13, 2],
        [1, 6, 9, 10, 0, 0] # Forced expiration data
    ]

#############
# Functions #
#############

def find_first_index_above_value(data, startSearchIndex, target):
    # first index grater than target in range data[startIndex:]
    for searchIndex in range(startSearchIndex, len(data)):
        if(data[searchIndex] > target):
            return searchIndex
    return len(data)


def quick_breath_split(breath, cutoff, RISINGEDGE=True):
    # Finds an edge where the edge has gone above/below value cutoff
    # Rising/falling edge selection chosen as input. Default to rising edge
    starts = []
    for index in range(1, len(breath)):
        if RISINGEDGE:
            if (breath[index] >= cutoff and breath[index-1] < cutoff):
                starts.append(index-1)
        else:
            if (breath[index] <= cutoff and breath[index-1] > cutoff):
                starts.append(index-1)
    return(starts)


def time_dependent_exp_fit(flow, time):
    # Exponential decay fit to dataSets given
    # NOTE: takes abs of flow, so any neg flow will be flipped
    ln_flw = [np.log(abs(f)) for f in flow]
    ones = [1]*len(time)

    dependent = np.array([ln_flw])
    independent = np.array([ones, times])
    try:
        res = lstsq(independent.T, dependent.T)
        decay = (res[0][1][0])
        offset = np.exp(res[0][0][0])
    except(ValueError):
        print('ValueError in exp decay fit: Data has nan?\nSetting result to NAN')
        decay = nan
        offset = nan
    return([decay, offset])


numFiles= len(files)

colours = ['blue', 'green', 'red', 'black', 'magenta', 'pink', 'cyan', 'brown', 'orange', 'yellow']
colourIndex = -1

# For results
gradients = [[] for m in range(numFiles)]
gradientsShort = [[] for m in range(numFiles)]
decays = [[] for m in range(numFiles)]
prev_max_pres = [[] for m in range(numFiles)] # store maximum pressure in shuttered breaths
shutterVol = [[] for m in range(numFiles)] # store maximum pressure in shuttered breaths
allAverageVol = [[] for m in range(numFiles)] # store maximum pressure in shuttered breaths
shutterStarts = [[] for m in range(numFiles)] # store maximum pressure in shuttered breaths
Rguess = [[] for m in range(numFiles)] # store maximum pressure in shuttered breaths

#=============================================================
#for plotting
STACK = 0
FULLSET = 1
SHUTTERSET = 0
STACKPRESSURE = 1

if(STACK):
    # plot stacks of flow up to and including shutter
    plt.rc('legend',**{'fontsize':12})
    a, (ax1) = plt.subplots(1, sharex=False)

if(SHUTTERSET):
    # plot full dataset at start of program
    plt.rc('legend',**{'fontsize':12})
    k, (kx1, kx2) = plt.subplots(2, sharex=False)

if(STACKPRESSURE):
    # plot flow separated into components
    plt.rc('legend',**{'fontsize':12})
    g, (gax1, gax2) = plt.subplots(2, sharex=True)
    u, (uax1) = plt.subplots(1)
    #h, (hax1) = plt.subplots(1)

#=============================================================


for i in range(numFiles):
    colourIndex += 1

    # For data
    flow = [0]
    volume = [0]
    time = [0]
    pressure = [0]
    breathPhase = [0]
    shutterState = [0]


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

            try:
                dp_flow = -float(row[dataLoc[i][0]]) # in L/s
            except(ValueError):
                dp_flow = np.nan
            try:
                dp_pressure = float(row[dataLoc[i][5]]) # in L/s
            except(ValueError):
                dp_pressure= np.nan
            try:
                dp_volume = -float(row[dataLoc[i][1]]) # in L
            except(ValueError):
                dp_volume = np.nan
            try:
                dp_time = float(row[dataLoc[i][2]]) # in s
            except(ValueError):
                dp_time= np.nan
            try:
                dp_phase = float(row[dataLoc[i][3]]) # 1=in, 2=out
            except(ValueError):
                dp_phase = np.nan
            try:
                dp_shutter = float(row[dataLoc[i][4]]) # 50 or 1-0
            except(ValueError):
                dp_shutter = np.nan

            flow.append(dp_flow)
            volume.append(dp_volume)
            time.append(dp_time)
            breathPhase.append(dp_phase)
            pressure.append(dp_pressure)
            # Doesn't exist every second breath, so interpolate ;)
            if np.isnan(dp_shutter):
                shutterState.append(shutterState[-1])
            else:
                shutterState.append(dp_shutter)

#    # Filter pressure
#    # Filter high resolution pressure data
    #fc = 60#filter cutoff
    #bw = 20#filter bandwidth (taps?)
    #flow = hamming(flow, fc, 300, bw, plot=False)
    #flow = np.real(flow).tolist()
    #volume = hamming(volume, fc, 300, bw, plot=False)
    #volume = np.real(volume).tolist()


    #######################
    # Split breaths apart #
    #######################
    flow_splits = quick_breath_split(breathPhase, 1.5)
    flow_splits.append(len(flow)-1)
    flow_splits_end = quick_breath_split(breathPhase, 1.5, RISINGEDGE=False)
    flow_splits_end.append(len(flow)-1)

    shutter_splits = quick_breath_split(shutterState, 0.5, RISINGEDGE=False)
    shutter_splits.append(len(flow)-1)
    shutter_splits_start = quick_breath_split(shutterState, 0.5)
    shutter_splits_start.append(len(flow)-1)


    if(FULLSET):
        # Plot flow, volume, points of interest
        plt.rc('legend',**{'fontsize':12})
        q, (qax1, qax2, qax3) = plt.subplots(3, sharex=True)

        qax1.plot(flow)
        qax1.plot(flow_splits, [flow[m] for m in flow_splits], 'or')
        qax1.plot(flow_splits_end, [flow[m] for m in flow_splits_end], 'ok')
        qax1.legend(['flow', 'expiratory starts', 'ends'])

        qax2.plot(volume)
        qax2.plot(flow_splits, [volume[m] for m in flow_splits], 'or')
        qax2.legend(['Volume and expiratory starts'])

        qax3.plot(pressure)
        qax3.plot(shutterState)
        qax3.plot(shutter_splits, [flow[m] for m in shutter_splits], 'om')
        qax3.plot(shutter_splits_start, [flow[m] for m in shutter_splits_start], 'oy')
        qax3.legend(['Flow and shutter ends', 'shutter'])
        q.show()


    ############################
    # Correct for volume drift #
    ############################
    # set final expiratory volume to 0 :)
    # NOTE: ignoring, because QV fit only looks at shape and trend
    #       actual value relative to 0 doesn't matter


    ##################
    # Start analysis #
    ##################
    FS = 1/(time[3]-time[2])
    averageFlow = [0]*int(4*FS) #average first 4 seconds
    numberAveraged = 0

    nextShutter = 0
    ANALYSING = False
    for breath in range(len(flow_splits)-1):
        print('\nbreath: {}'.format(breath))

        if(shutter_splits[nextShutter] < flow_splits_end[breath]):
            ANALYSING = True
            print('This breath is shuttered')
            nextShutter += 1

        #######################
        # Split breaths apart #
        #######################
        flw = flow[flow_splits[breath]:flow_splits_end[breath]]
        vol = volume[flow_splits[breath]:flow_splits_end[breath]]
        vol = [v-vol[0] for v in vol] # shift volume curve down
        t = time[flow_splits[breath]:flow_splits_end[breath]]
        pres = pressure[flow_splits[breath]:flow_splits_end[breath]]

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~``
        ################################################################
        # Add flow data to averaging array when not a shuttered breath #
        # We are averaging the QV loop for the subject                 #
        ################################################################
        if not ANALYSING:
            for index in range(min(len(averageFlow), len(flw))):
                averageFlow[index] += flw[index]
            numberAveraged += 1

            if(SHUTTERSET):
                newt = [tim - t[0] for tim in t]
                #ax1.plot(vol, flw, 'k')
                #ax2.plot(newt, flw, 'k')

            if(STACK):
                ax1.plot(flw, color=colours[colourIndex])

            if(STACKPRESSURE):
                gax1.plot(flw, color='blue')
                #gax2.plot(pres, color='blue')


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~``
        ###########################################################
        # When analysing data:                                    #
        # 1) Get averaged data                                    #
        # 2) Calc QV loop for breath                              #
        # 3) remove averaged data to leave only shuttering effect #
        # 4) calculate gradient of QV loop                        #
        # 5) Success or cry!                                      #
        ###########################################################
        elif ANALYSING:
            # Step 1:
            averageFlow = [q/float(numberAveraged) for q in averageFlow]
            averageVol = integral(averageFlow, FS)
            allAverageVol[i].append(averageVol[-1])

            # Where is the shutter in the breath?
            shutterStartPoint = shutter_splits[nextShutter-1] - flow_splits[breath]
            shutterStarts[i].append(shutterStartPoint)

            if(STACK):
                ax1.plot(averageFlow, 'purple')
                ax1.plot(flw, 'r', linewidth=3)
                ax1.grid()
                ax1.set_ylabel('Flow', fontsize=30)
                a.show()
                plt.rc('legend',**{'fontsize':12})
                a, (ax1) = plt.subplots(1, sharex=False)

            ########################
            # Estimate flow due to shutter
            AminusBTrue = [flw[k]-averageFlow[k] for k in range(min(len(flw),len(averageFlow)))]
            AminusB = [q - AminusBTrue[shutterStartPoint] for q in AminusBTrue]

            fc = 50#filter cutoff
            bw = 10#filter bandwidth (taps?)
            AminusBFilt = hamming(AminusBTrue, fc, 300, bw, plot=False)
            AminusBFilt = np.real(AminusBFilt).tolist()
            AminusBVolFilt = integral(AminusBFilt, FS)

            AminusBVol = integral(AminusBTrue, FS)

            ########################
            # find max shutter flow, and the end point in volume
            trueMaxFlowIndex = AminusBFilt.index(max(AminusBFilt[min(shutterStartPoint,len(AminusBFilt)-2):min(shutterStartPoint+100,len(AminusBFilt)-1)]))
            maxFlowIndex = trueMaxFlowIndex + 1


            if(0):
                ####################
                # Check for peak after max.
                # Start looking for new peak from 10ml after maxFlowIndex
                # Shift maxFlowIndex to new peak, if there is one
                targetVolume = min(AminusBVol[trueMaxFlowIndex] + 0.01, AminusBVol[-2])
                startVol = find_first_index_above_value(AminusBVol, maxFlowIndex, targetVolume)

                nextPeak = AminusB.index(max(AminusB[startVol:startVol+50]))
                if(nextPeak > startVol):
                    # new peak found
                    maxFlowIndex = nextPeak


            ########################
            # find 30ml after shutter start
            # 2.5ml ~ 10% (assume)
            volShuttered = averageVol[shutterStartPoint] - averageVol[shutter_splits_start[nextShutter-1]-flow_splits[breath]]
            print(volShuttered)
            targetVolume = min(AminusBVol[trueMaxFlowIndex] + volShuttered/4, AminusBVol[-2])
            thirtyMlAfterMax = find_first_index_above_value(AminusBVol, maxFlowIndex, targetVolume)


            #######################
            # find end as function of (Q increase, V stagnation)
            flowInflection = len(AminusBVol)-1
            # Has the flow increased again in the future?
            # If so, inflection in flow (ideal: it would always decrease)
#            searchRange = 4
#            for index in range(maxFlowIndex, len(AminusBVolFilt)-searchRange):
#                if (AminusBFilt[index+searchRange] > AminusBFilt[index]):
#                    flowInflection = index
#                    break

            searchRange = 8
            volumeStagnation = len(AminusBVol)-1
            # Has the volume stagnated or started decreasing
            # expect large-ish increase every dp
            # call stagnation when volume increased less than a minimum amount
            for index in range(maxFlowIndex, len(AminusBVolFilt)-searchRange):
                if (AminusBVolFilt[index+searchRange] < AminusBVolFilt[index]+0.005): # 0.005ml change at 0.15L/s
                    volumeStagnation = index
                    break

            endFlowIndex = min(flowInflection, volumeStagnation, thirtyMlAfterMax)

            print(maxFlowIndex)
            print(endFlowIndex)
            print(shutterStartPoint)


            if(0):
                ###################
                # Calc grad from decay rate flow
                decayLen = min(25, len(flw)-maxFlowIndex-1)
                decayFlow = AminusB[maxFlowIndex:maxFlowIndex+decayLen]
                times = t[maxFlowIndex:maxFlowIndex+decayLen]

                decay, offset = time_dependent_exp_fit(decayFlow, times)

                flwPulseCurve = [offset*exp(decay*y) for y in times]
                print('\n     Decay rate stuff:')
                print('Qt decay rate: {}'.format(decay, offset,))



            #######################
            # Calc gradient of QV line
            if(endFlowIndex-maxFlowIndex > 3 and shutterStartPoint > 0): #min 3dp needed for least-squares

                #####################
                # fit to full identified linear length
                volume_curve = AminusBVol[maxFlowIndex:endFlowIndex]
                flow_curve = AminusB[maxFlowIndex:endFlowIndex]
                gradient_passive, intercept, r_value, p_value, std_err = stats.linregress(volume_curve, flow_curve)

                print('QV gradient full range: {}'.format(gradient_passive))
                print('r value full range: {}'.format(r_value))

                gradients[i].append(gradient_passive)
                line = [intercept + gradient_passive*m for m in [AminusBVol[endFlowIndex], AminusBVol[maxFlowIndex]]]
                lineFull = [intercept + gradient_passive*m for m in [min(AminusBVol), max(AminusBVol)]]

                #####################
                # fit to first 30ml identified linear length
                volume_curve = AminusBVol[maxFlowIndex:thirtyMlAfterMax]
                flow_curve = AminusB[maxFlowIndex:thirtyMlAfterMax]
                gradient_passive, intercept, r_value, p_value, std_err = stats.linregress(volume_curve, flow_curve)

                print('QV gradient 30ml: {}'.format(gradient_passive))
                print('r value 30ml: {}'.format(r_value))

                gradientsShort[i].append(gradient_passive)
                lineShort = [intercept + gradient_passive*m for m in [AminusBVol[thirtyMlAfterMax], AminusBVol[maxFlowIndex]]]
                lineFullShort = [intercept + gradient_passive*m for m in [min(AminusBVol), max(AminusBVol)]]


                ########################
                # EXCEPTIONS: SET NAN  #
                ########################
                if(maxFlowIndex > (shutterStartPoint + len(AminusB[shutterStartPoint:])/3)):
                    print('**\nSecond peak identified more than 1/3 data length after trueMaxIndex. Setting result to NAN\n')
                    if(gradients[i]):
                        gradients[i][-1] = np.nan
                    else:
                        gradients[i].append(np.nan)

#                if(shutterStartPoint > 75):
#                    print('**\nTempish test to adjust for inconsistent shuttering time. Shutter opened after 0.37s, result -> NAN\n')
#                    if(gradients[i]):
#                        gradients[i][-1] = np.nan
#                    else:
#                        gradients[i].append(np.nan)


                if(abs(r_value) < 0.85):
                    print('\tR value <0.85')

                if(gradient_passive < -30):
                    print('\tGradient >30. TOO BIG?')
            else:
                    print('**\nmaxFlowIndex to endFlowIndex less than 3. Cannot calculate results, recording NAN\n')
                    gradients[i].append(np.nan)

                    if(STACKPRESSURE):
                        uax1.plot(AminusBVol, AminusB, linewidth = 3)
                        uax1.plot(AminusBVol, AminusBTrue, ':', linewidth = 3)
                        uax1.plot(AminusBVolFilt, AminusBFilt, '.-', linewidth = 3)
                        uax1.plot(AminusBVol[maxFlowIndex], AminusB[maxFlowIndex], 'o', color='yellow')
                        uax1.plot(AminusBVol[endFlowIndex], AminusB[endFlowIndex], 'o', color='yellow')
                        uax1.plot(AminusBVol[trueMaxFlowIndex], AminusB[trueMaxFlowIndex], 'o', color='magenta')
                        uax1.plot(AminusBVol[thirtyMlAfterMax], AminusB[thirtyMlAfterMax], 'o', color='green')

                        uax1.set_ylabel('Flow', fontsize=30)
                        uax1.set_xlabel('Volume', fontsize=30)


            print('\nmax pressure measured in breath = {}'.format(max(pres)))
            print(('Prev pressures: {}'.format(prev_max_pres[i])))
            print(('NOTE: Increased max pressure indicates increased breathing effort'))
            prev_max_pres[i].append(max(pres))
            shutterVol[i].append(vol[-1] - vol[shutterStartPoint])
            #if(not np.isnan(gradients[i][-1])):
            #    shutterVol[i].append(vol[-1] - vol[shutterStartPoint])
            #else:
            #    shutterVol[i].append(np.nan)

#            print(('bundy airway resistance'))
#            print('pres just after shutter (like 2 dp later) : {}'.format(pres[maxFlowIndex]))
#            print('flow at same point : {}'.format(flw[maxFlowIndex]))
#            print('pressure change: {}'.format(prev_max_pres[i][-1]-pres[maxFlowIndex]))
#            print('so R guess is dp/dQ = {}'.format((prev_max_pres[i][-1]-pres[maxFlowIndex])/flw[maxFlowIndex]))
#            Rguess[i].append((prev_max_pres[i][-1]-pres[maxFlowIndex])/flw[maxFlowIndex])



            if(STACKPRESSURE):
                # flow and pressure plots
                gax1.plot(averageFlow[:len(flw)], color='purple', linewidth=1.5)
                gax1.plot(flw, color='red', linewidth=3)
                gax1.plot(AminusB, color='black', linewidth=3)
                gax1.plot(AminusBTrue,':', color='black', linewidth=3)
                #gax1.plot(range(maxFlowIndex, maxFlowIndex+decayLen), flwPulseCurve, color='orange', linewidth=2)
                gax1.plot(maxFlowIndex, flw[maxFlowIndex], 'o', color='yellow')
                gax1.plot(endFlowIndex, flw[endFlowIndex], 'o', color='yellow')
                gax1.grid()
                gax1.set_ylabel('Flow', fontsize=30)

                gax2.plot(pres, color='red', linewidth=3)
                gax2.grid()
                gax2.set_xlabel('DP', fontsize=31)
                gax2.set_ylabel('Pressure', fontsize=30)
                # flow vs volume plot
                if(not np.isnan(gradients[i][-1])):
                    uax1.plot(AminusBVol, AminusB, linewidth = 3)
                    uax1.plot(AminusBVol[maxFlowIndex:endFlowIndex], AminusB[maxFlowIndex:endFlowIndex], color='orange', linewidth=2)
                    uax1.plot(AminusBVolFilt, AminusBFilt, '.-', linewidth = 3)
                    uax1.plot([min(AminusBVol), max(AminusBVol)], lineFull, color='pink', linewidth=2)
                    uax1.plot([AminusBVol[endFlowIndex], AminusBVol[maxFlowIndex]], line, color='purple', linewidth=2)
                    uax1.plot(AminusBVol[maxFlowIndex], AminusB[maxFlowIndex], 'o', color='yellow')
                    uax1.plot(AminusBVol[endFlowIndex], AminusB[endFlowIndex], 'o', color='yellow')
                    uax1.plot(AminusBVol[trueMaxFlowIndex], AminusB[trueMaxFlowIndex], 'o', color='magenta')
                    uax1.plot(AminusBVol[thirtyMlAfterMax], AminusB[thirtyMlAfterMax], 'o', color='green')
                    uax1.grid()
                    uax1.set_ylabel('Flow', fontsize=30)
                    uax1.set_xlabel('Volume', fontsize=30)

                    #hax1.plot([0, (AminusBVol[thirtyMlAfterMax]-AminusBVol[trueMaxFlowIndex])],[0, (AminusB[thirtyMlAfterMax]-AminusB[trueMaxFlowIndex])], '-o', color='cyan')
                plt.show()
                plt.rc('legend',**{'fontsize':12})
                g, (gax1, gax2) = plt.subplots(2, sharex=True)
                u, (uax1) = plt.subplots(1)

                #plt.show()
                #g, (gax1, gax2) = plt.subplots(2, sharex=True)
                #u, (uax1) = plt.subplots(1)



            # look at missing volume
            #integrate flow between shutter_splits and shutter_splits_start
            #shutterVol = integral(oldFlw[shutterStartPoint:startPoint], FS)
            #peakyVol = integral(oldFlw[startPoint:endFlowIndex], FS)
            #peakyVol2 = integral(flw[startPoint:endFlowIndex], FS)
            #print('volMissing = {}'.format(oldVol[-1] - vol[-1]))
            #print('shutterVol: {}'.format(shutterVol[-1]))
            #print('PeakyVol = {}'.format(peakyVol2[-1] - peakyVol[-1]))
            #print('shutter-peaky:= {}'.format((shutterVol[-1])-(peakyVol2[-1] - peakyVol[-1])))


            if(SHUTTERSET):
                # Plot flow, volume, points of interest...
                ax1.set_ylabel('Flow', fontsize=30)
                ax1.set_xlabel('volume', fontsize=30)
                ax1.plot(AminusBVol, AminusB, 'b')
                #ax1.plot([AminusBVol[shutterStartPoint+decayLen], AminusBVol[shutterStartPoint]], line, 'g')
                ax1.plot(AminusBVol[shutterStartPoint], AminusB[shutterStartPoint], 'ro')
                #ax1.plot(AminusBVol[shutterStartPoint+decayLen], AminusB[shutterStartPoint+decayLen], 'mo')
                ax1.legend(['Shuttered data (Measured-average)',
                            'Linear fit to shuttered region'])

                newt = [tim - t[0] for tim in t]
                ax2.plot(newt, flw, 'b')
                ax2.plot(newt, averageFlow[:len(newt)], 'g')
                ax2.plot(newt[:len(AminusB)], AminusB, color='red')
                #ax2.plot(newt[shutterStartPoint:shutterStartPoint+decayLen], flwPulseCurve, color='purple')
                ax2.plot(newt[shutterStartPoint], AminusB[shutterStartPoint], 'or')
                #ax2.plot(newt[shutterStartPoint+decayLen], AminusB[shutterStartPoint+decayLen], 'mo')
                ax2.set_ylabel('Flow', fontsize=30)
                ax2.set_xlabel('Time', fontsize=30)
                ax2.legend(['Measured flow',
                            'Flow due to body',
                            'Flow due to shutter',
                            'Decay fits'
                            ])
                ax1.grid()
                ax2.grid()
                plt.show()

            averageFlow = [0]*len(averageFlow)
            numberAveraged = 0
            ANALYSING = False
            if(SHUTTERSET):
                plt.rc('legend',**{'fontsize':12})
                f, (ax1, ax2) = plt.subplots(2, sharex=False)

print(Rguess)
maxEdRecorded = [[] for m in range(numFiles)] # store maximum pressure in shuttered breaths
maxEdRecordedBoxy = [[] for m in range(numFiles)] # store maximum pressure in shuttered breaths
gradientsBoxy = [[] for m in range(numFiles)] # store maximum pressure in shuttered breaths
maxPressuresBoxy = [[] for m in range(numFiles)] # store maximum pressure in shuttered breaths
shutterVolBoxy = [[] for m in range(numFiles)] # store maximum pressure in shuttered breaths
for i in range(numFiles):
    maxEdRecorded[i] = [prev_max_pres[i][j]/shutterVol[i][j] for j in range(len(shutterVol[i]))]
    maxEdRecordedBoxy[i] = [v for v in maxEdRecorded[i] if not np.isnan(v)]
    maxEdRecordedBoxy[i] = [v for v in maxEdRecorded[i] if (v<1200)]
    gradientsBoxy[i] = [v for v in gradients[i] if not np.isnan(v)]
    maxPressuresBoxy[i] = [v for v in prev_max_pres[i] if not np.isnan(v)]
    shutterVolBoxy[i] = [v for v in shutterVol[i] if not np.isnan(v)]

plt.figure()
for i in range(numFiles):
    plt.plot(gradients[i], prev_max_pres[i], 'o', color=colours[i])
    plt.plot(gradients[i], maxEdRecorded[i], 'x', color=colours[i])
    plt.plot(shutterVol[i], prev_max_pres[i], 's', color=colours[i])
    #plt.plot(Rguess[i], prev_max_pres[i], 'd', color=colours[i])
plt.xlabel('grad')
plt.ylabel('max Pres')
plt.grid()

plt.figure()
plt.boxplot(allAverageVol, labels=xLabels)
for i in range(numFiles):
    y = allAverageVol[i]
    x = np.random.normal(1+i, 0.04, size=len(y))
    plt.plot(x, y, 'r.')
plt.ylabel('averageVol')
plt.grid()

plt.figure()
plt.boxplot(shutterStarts, labels=xLabels)
for i in range(numFiles):
    y = shutterStarts[i]
    x = np.random.normal(1+i, 0.04, size=len(y))
    plt.plot(x, y, 'r.')
plt.ylabel('shutterStarts')
plt.grid()

plt.figure()
plt.boxplot(shutterVolBoxy, labels=xLabels)
for i in range(numFiles):
    y = shutterVolBoxy[i]
    x = np.random.normal(1+i, 0.04, size=len(y))
    plt.plot(x, y, 'r.')
plt.ylabel('shutterVol')
plt.grid()

plt.figure()
plt.boxplot(maxPressuresBoxy, labels=xLabels)
for i in range(numFiles):
    y = maxPressuresBoxy[i]
    x = np.random.normal(1+i, 0.04, size=len(y))
    plt.plot(x, y, 'r.')
plt.ylabel('maxPressure')
plt.grid()

plt.figure()
plt.boxplot(maxEdRecordedBoxy, labels=xLabels)#, 'Labels',())
for i in range(numFiles):
    y = maxEdRecordedBoxy[i]
    x = np.random.normal(1+i, 0.04, size=len(y))
    plt.plot(x, y, 'r.')
plt.ylabel('Pseudo Elastance')
plt.grid()

plt.figure()
plt.boxplot(gradientsBoxy, labels=xLabels)#, 'Labels',())
for i in range(numFiles):
    y = gradientsBoxy[i]
    x = np.random.normal(1+i, 0.04, size=len(y))
    plt.plot(x, y, 'r.')
plt.ylabel('Decayt rates')
plt.grid()


dataBoxy = []
for i in range(numFiles):
    if i == 1:
        state = '1'
    else:
        state =  '0'
    dataBoxy+=([[v, state] for v in gradients[i] if not np.isnan(v)])
df = pd.DataFrame(dataBoxy, columns=['Decay rate', 'Banded'])
sns.set(style="whitegrid")
plt.figure()
sns.boxplot(x='Banded', y='Decay rate', data=df)#, 'Labels',())
plt.grid()

plt.show()


if(STACK):
    ax1.grid()
    ax1.set_xlabel('Volume', fontsize=31)
    ax1.set_ylabel('Flow', fontsize=30)
    plt.legend(['Blue = unbanded (1A)', 'Green = banded (1B)'])
    plt.show()


print('\n')
print('res:')
mean_grad = []
for i in range(numFiles):
    mean_grad.append(np.nanmean(gradients[i]))
    print(mean_grad[-1])
    print(np.nanstd(gradients[i], axis=0))

Rx = [#2.3,
        0.8,
        0.01,
        0.01,
        0
        ]
Es = [1]*numFiles

dependents = [[mean_grad[i]*Rx[i] for i in range(numFiles)]]
independents = [[-mean_grad[i] for i in range(numFiles)], Es]

print(gradients)
for i in range(numFiles):
    plt.plot(gradients[i], '-o', color=colours[i])
    plt.plot([mean_grad[i]]*10, '--', color=colours[i])
plt.show()

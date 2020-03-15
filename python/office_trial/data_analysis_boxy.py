#!/bin/bash

# Import my extensions
import sys
sys.path.insert(0, '/home/sarah/Documents/Spirometry/python/extensions')

import csv
import matplotlib.pyplot as plt
import numpy as np
from data_struct import DataStore as Data
from calculus import integral
from breath_analysis import split_breaths
from numpy import array
from numpy.linalg import lstsq
from scipy import stats
import locale
import seaborn as sns
import pandas as pd

# data has comma instead of decimal point

dataset = 0
if(dataset == 0):
    # Banded data set, has FEV in other sets too
    path = './trial_data/'
    subjects = [
            #  'Sarah',
              'F01',
              'F02',
              'F03',
              'F04',
              'F05',
              'F09',
              'M01',
              'M02',
              'M03',
              'M04',
              'M05',
              'M06',
              'M07',
              'M08',
              'M09',
              'M10',
              ]
    test_types = [
            #'_mp_peak.csv',
            #'_mp_peakRx.csv',
            '_mp_noRx.csv',
            '_mp_Rx.csv',
            ]
    tests = [
            ''
            ]


    scaling = [
            [[1,1,1]],
            ]
    FILEDELIMITER = ','
    indexing = [0,1]
    locale._override_localeconv["decimal_point"] = "."
    FS = 300
    PEAKHEIGHT = 0.0001
    MINFLOW = 0.2
    #MINFLOW = 0.06
    MINSHUTTERPRESSURE = 0.0006
    MINVOL = 0.5



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Extract pressure two data arrays (assumed pressure and flow,
# but could be pressure & pressure etc) and time stamp from csv
# NOTE This function appends. Make sure empty arrays are passed in
def extract_data_arrays_from_file(filename, indexing):
    # Store pressure and flow data
    # data:
    # flow, temperature, Patm(hPa), humidity(%), Pcabin, Pmouth, CO%, CH4%
    pressure = []
    flow = []
    #time = []

    with open(filename, 'rb') as csvfile:
        #reader = csv.reader(csvfile, delimiter = '\t')
        reader = csv.reader(csvfile, delimiter = FILEDELIMITER)

        # skip header line
        header = reader.next()

        for row in reader:

            try:
                dp_flow = locale.atof(row[indexing[0]]) # in L/s
            except(ValueError):
                dp_flow = np.nan
            try:
                dp_pressure = locale.atof(row[indexing[1]]) # in L/s
            except(ValueError):
                dp_pressure = np.nan
            #try:
            #    dp_time = float(row[]) # in s
            #except(ValueError):
            #    dp_time = np.nan

            flow.append(-dp_flow)
            #time.append(dp_time)
            pressure.append(-dp_pressure)
    return(pressure, flow)


def root_mean_squared_error(data, fit, peak):
    mse = 0
    N = len(data)
    for m in range(N):
        mse += (data[m]/peak - fit[m]/peak)**2
    mse /= N
    rmse = np.sqrt(mse)
    return rmse

def time_dependent_exp_fit(flow, times):
    # Exponential decay fit to dataSets given
    # NOTE: takes abs of flow, so any neg flow will be flipped
    ln_flw = [np.log(abs(f)) for f in flow]
    ones = [1]*len(times)

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

def check_StartMiddleStop(starts, middles, stops, minDistance):
    finStarts = [starts[0]]
    finMiddles= []
    finStops = []

    i = 0
    j=0
    k=0
    l=1
    #for i in range(min(len(starts)-1, len(middles), len(stops))):
    haveData = True
    while haveData:
        # find middle
        if(haveData):
            try:
                finMiddles.append(middles[[pt>minDistance for pt in [m-finStarts[-1] for m in middles]].index(True)])
                finMiddles
            except(ValueError):
                haveData = False
        # find end
        if(haveData):
            try:
                finStops.append(stops[[pt>minDistance for pt in [s-finMiddles[-1] for s in stops]].index(True)])
                finMiddles
            except(ValueError):
                haveData = False
        #find start
        if(haveData):
            try:
                finStarts.append(starts[[pt>0 for pt in [s-finStops[-1] for s in starts]].index(True)])
                finMiddles
            except(ValueError):
                haveData = False
    finalLength = min(len(finStarts), len(finMiddles), len(finStops))
    return(finStarts[:finalLength], finMiddles[:finalLength], finStops[:finalLength])


def shearTransformToFindMax(data):
    if(len(data) > 3):
        m = (data[-1] - data[0])/len(data)
        line = [i*m + data[0] for i in range(len(data))]

        shearData = [data[i]-line[i] for i in range(len(data))]
        shearMax = shearData.index(max(shearData))

        if shearMax > len(data)/2:
            shearMax = 0

        if(0):
            plt.figure()
            plt.plot(data)
            plt.plot(line, 'r')
            plt.plot(shearData, 'pink')
            plt.plot(shearMax, shearData[shearMax], 'o')
            plt.show()


    else:
        shearMax = 0
    return shearMax

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

numSubjects = len(subjects)
numTestFiles = len(test_types)
numTests = len(tests)

decaysFull = [[[[] for t in range(numTests)] for tf in range(numTestFiles)] for s in range(numSubjects)]
decays = [[[[] for t in range(numTests)] for tf in range(numTestFiles)] for s in range(numSubjects)]
prev_max_pres = [[[[] for t in range(numTests)] for tf in range(numTestFiles)] for s in range(numSubjects)]
shutterVol = [[[[] for t in range(numTests)] for tf in range(numTestFiles)] for s in range(numSubjects)]

dataBoxy = []

PLOTTING_BREATH = 0
PLOTTING_DATASET = 0


for subject in range(numSubjects):
    print('\nNEW subject\n--------\n')
    print('Subject: {}'.format(subjects[subject]))

    for testFile in range(numTestFiles):
        print('Test file: {}'.format(test_types[testFile]))

        for test in range(numTests):
            print('Test: {}\n'.format(tests[test]))

            if not(subjects == ['']):
                if not(tests == ['']):
                    filename = path + subjects[subject] +'/'+ test_types[testFile] +'/'+ tests[test]
                else:
                    filename = path + subjects[subject] +'/'+ subjects[subject]+test_types[testFile]
            else:
                filename = path + tests[test]
            print(filename)
            print('\n')

            try:
                # Extract data to empty arrays.
                # This will fail if file does not exist
                pressure = []
                flow = []
                if (dataset==3):
                    pressure, flow = extract_data_arrays_from_file(filename, indexing[test])
                else:
                    pressure, flow = extract_data_arrays_from_file(filename, indexing)

                if dataset==0:
                    pressure = [p*1000 for p in pressure]
                    flow = [p*0.0092 for p in pressure]


                # scale, some are in mPa, mL/s, some are Pa, L/s
              #  scale = scaling[subject][testFile][test]
                scale = -1
                if scale != 1:
                    pressure = [p/float(scale) for p in pressure]
                    flow = [f/float(scale) for f in flow]

                # Don't have time data, so guessing FS
                time = range(len(flow))
                for i in range(len(time)):
                    time[i] /= float(FS)


                # Find the start and end of expiration
                # make peak low. get as many splits as possible for shutter blocked flow
                # later code will filter it
#                flow_splits = split_breaths(flow, peak_height=PEAKHEIGHT, Fs=FS, plot=False)
#                flow_starts = flow_splits[0]
#                flow_middles = flow_splits[1]
#                flow_stops = flow_splits[2]

                flow_starts = quick_breath_split(flow, 0.1)
                flow_middles = quick_breath_split(flow, -0.1, RISINGEDGE=False)
                flow_stops = [s-1 for s in flow_starts[1:]]
                flow_stops.append(len(flow)-1)
                flow_starts, flow_middles, flow_stops = check_StartMiddleStop(flow_starts, flow_middles, flow_stops, minDistance=FS/4)




                if(PLOTTING_DATASET):
                    plt.plot(pressure, 'b')
                    plt.plot(flow, 'g')
                    plt.plot(flow_starts, [flow[s] for s in flow_starts], 'go')
                    plt.plot(flow_middles, [flow[s] for s in flow_middles], 'yo')
                    plt.plot(flow_stops, [flow[s] for s in flow_stops], 'ro')
                    plt.grid()
                    plt.show()


                for breath in range(len(flow_middles)):
                    # return all decay rates of expiration
                    flwExp = flow[flow_middles[breath]:flow_stops[breath]]
                    flwInsp = flow[flow_starts[breath]:flow_middles[breath]]
                    # Make all data +ve, all -ve data after flipping -> 0
                    for f in range(len(flwExp)):
                        if flwExp[f] >=0:
                            flwExp[f] = 0
                        else:
                            flwExp[f] *= -1
                    for f in range(len(flwInsp)):
                        if flwInsp[f] <=0:
                            flwInsp[f] = 0
                        else:
                            flwInsp[f] *= 1


                    # average flow last 1/4 of insp flow curve
                    inspMax = np.mean(flwInsp[int(len(flwInsp)-FS):])
                    if(inspMax>0.05): # if middle picked before pause
                        inspMax = np.mean(flwExp[:FS])

                    print('~~~')
                    print(inspMax)


                    startPt = min(flwExp.index(max(flwExp)),len(flwExp)-1)
                    if(len(flwExp) > 40 and startPt >= 0):
                        vol = integral(flwExp, FS)
                        timeExp = time[flow_middles[breath]:flow_stops[breath]]
                        timeExp = [t-timeExp[0] for t in timeExp]

                        # start at top. Full end at 1 second
                        startPt += shearTransformToFindMax(flwExp[startPt:min(len(flwExp)-2,startPt+FS)])

                        # Find index of first 20% volume drop
                        volMax = vol[startPt]
                        def find_volEnd(volEnd, startPt, vol):
                            i = startPt
                            endIndex = startPt
                            while (i <= len(vol)):
                                if vol[i] >= volEnd:
                                    endIndex = i
                                    i = len(vol)+ 1
                                i += 1
                            return endIndex


                        volEnd  = volMax + (vol[-1] - vol[startPt])*0.20
                        oneLitrePt = find_volEnd(volEnd, startPt, vol)
                        volEnd  = volMax + (vol[-1] - vol[startPt])*0.80
                        endPt = find_volEnd(volEnd, startPt, vol)
                        volEnd  = (vol[-1])*0.4
                        maxStartPt = find_volEnd(volEnd, startPt, vol)


                        if(PLOTTING_BREATH):
                            plt.plot(flow, 'g')
                            plt.plot(range(flow_middles[breath],flow_stops[breath]), flow[flow_middles[breath]:flow_stops[breath]], 'r')
                            plt.plot(range(flow_middles[breath]+startPt,flow_middles[breath]+endPt), flow[flow_middles[breath]+startPt:flow_middles[breath]+endPt], 'm')
                            plt.plot(range(flow_middles[breath]+startPt,flow_middles[breath]+oneLitrePt), flow[flow_middles[breath]+startPt:flow_middles[breath]+oneLitrePt], 'c')
                            plt.grid()
                            plt.figure()
                            plt.plot(vol[startPt:endPt], flwExp[startPt:endPt])
                            plt.plot(vol[startPt:oneLitrePt], flwExp[startPt:oneLitrePt], 'c')
                            plt.grid()
                            plt.show()

                        analysed = False
                        if(endPt-startPt)>3:

                            gradient_passive, intercept, r_value, p_value, std_err = stats.linregress(vol[startPt:endPt], flwExp[startPt:endPt])
                            #gradient_passive, intercept, r_value, p_value, std_err = stats.linregress(vol[startPt:oneLitrePt], flwExp[startPt:oneLitrePt])
                            if(max(np.abs(flwExp)) > MINFLOW and inspMax < 0.05 and startPt < maxStartPt):
                                # Look at 1s data first

                                if(0):
                                    plt.plot(flow, 'g')
                                    plt.plot(range(max(0,flow_middles[breath]-FS),flow_middles[breath]), flow[max(0,flow_middles[breath]-FS):flow_middles[breath]], 'y')
                                    plt.plot(range(flow_middles[breath],flow_stops[breath]), flow[flow_middles[breath]:flow_stops[breath]], 'r')
                                    plt.grid()
                                print('Ooh, looks like we\'ve got a good one!!\n')

                                analysed = True
                                print('Full range is a go!!')

                                print('QV gradient: {}'.format(gradient_passive))
                                print('r value: {}'.format(r_value))

                                decaysFull[subject][testFile][test].append(gradient_passive)
                                lineFull = [intercept + gradient_passive*m for m in vol[startPt:endPt]]

                                results = time_dependent_exp_fit(flwExp[startPt:endPt], timeExp[startPt:endPt])
                                print('decay from exp fit: {}'.format(results[0]))
                                print('offset from exp fit: {}'.format(results[1]))
                                lineDecayFull = [results[1]*np.exp(results[0]*t) for t in timeExp[startPt:endPt]]


                                # Look at 1 Litre data
                                if(oneLitrePt-startPt)>3:
                                    print('\nLooking at the first 20% volume too!!')
                                    gradient_passive, intercept, r_value, p_value, std_err = stats.linregress(vol[startPt:oneLitrePt], flwExp[startPt:oneLitrePt])

                                    print('QV gradient: {}'.format(gradient_passive))
                                    print('r value: {}'.format(r_value))

                                    decays[subject][testFile][test].append(gradient_passive)
                                    line = [intercept + gradient_passive*m for m in vol[startPt:oneLitrePt]]

                                    results = time_dependent_exp_fit(flwExp[startPt:oneLitrePt], timeExp[startPt:oneLitrePt])
                                    print('decay from exp fit: {}'.format(results[0]))
                                    print('offset from exp fit: {}'.format(results[1]))
                                    lineDecay = [results[1]*np.exp(results[0]*t) for t in timeExp[startPt:oneLitrePt]]


                                #dataBoxy.append([gradient_passive, subject, testFile, test])
                                if(testFile == 0):
                                    RxState = 'No external R'
                                else:
                                    RxState = 'External R = 0.8cmH2O/L/s'
                                if(subjects[subject]=='F09'):
                                    dataBoxy.append([decaysFull[subject][testFile][test][-1], 'F06', RxState, test])
                                else:
                                    dataBoxy.append([decaysFull[subject][testFile][test][-1], subjects[subject], RxState, test])

#                                if(max(np.abs(flow[max(0,flow_middles[breath]-FS):flow_middles[breath]])) > MINSHUTTERPRESSURE):
#                                    print('**hard breath')
#                                    dataBoxy[-1].append("Big effort")
#                                else:
#                                    dataBoxy[-1].append("Small effort")
#                                    print('**soft breath')
#
#                                if(vol[-1] > MINVOL):
#                                    print('**big breath')
#                                    dataBoxy[-1].append("Big volume")
#                                else:
#                                    dataBoxy[-1].append("Small volume")
#                                    print('**small breath')



                        else:
                            print('Yeah, kind of low flow. Probs not a FE breath')

                        if(analysed and PLOTTING_BREATH):
                            plt.figure()
                            plt.plot(flwInsp, 'orange')
                            plt.plot(flwExp, 'red')
                            if(analysed):
                                plt.plot(range(startPt,endPt), lineDecayFull, 'orange', linewidth=2)
                                if(oneLitrePt-startPt)>3:
                                    plt.plot(range(startPt,oneLitrePt), lineDecay, 'orange', linewidth=2)
                            plt.legend(['Insp', 'Exp'], fontsize=38)
                            plt.ylabel('Flow', fontsize=38)
                            plt.xlabel('dp', fontsize=38)
                            plt.grid()

                            plt.figure()
                            plt.plot(vol, flwExp, 'b', linewidth = 2)
                            plt.plot(vol[maxStartPt], flwExp[maxStartPt], 'ok')
                            if(analysed):
                                plt.plot(vol[startPt:endPt], lineFull, 'magenta', linewidth=2)
                                if(oneLitrePt-startPt)>3:
                                    plt.plot(vol[startPt:oneLitrePt], line, 'purple', linewidth=2)
                            plt.ylabel('Flow', fontsize=38)
                            plt.xlabel('Volume', fontsize=38)
                            plt.grid()

                            plt.show()

                    else:
                        print('Expiration data less than 40 long, os 0.4s')

                #Basic table print of all results (SF highlighted with *)
                print('\n')
                print('____ results for ____')
                print(filename)
                print('________')
                for n in range(len(decaysFull[subject][testFile][test])):
                    print(decaysFull[subject][testFile][test][n])


            except IOError:
                print("File {} does not exist".format(filename))
#
##
###
####
#####
######
#####
###
##
#

# Variation analysis



# Mean analysis
for s in range(len(subjects)):
    print('~~~~~ subject {} ~~~~~'.format(s))
    for tst in range(len(test_types)):
        arr = []
        for n in range(len(decaysFull[s][tst][0])):
            arr.append(decaysFull[s][tst][0][n])
        print('mean: {}'.format(np.mean(arr)))
        print('sd: {}'.format(np.std(arr)))




if(dataBoxy):
    plt.figure()
    sns.set(style="whitegrid", font_scale=2.5)
    df = pd.DataFrame(dataBoxy, columns=['Decay rate','Subject', 'External R added', 'Test'])#,'Effort','Volume'])
    ax=sns.boxplot(x='Subject', y="Decay rate", hue='External R added',  data=df)
    sns.stripplot(x="Subject", y="Decay rate", hue="External R added", data=df, size=7, color="black", split=True)
    plt.xlabel("Subject")

    # Get the handles and labels. For this example it'll be 2 tuples
    # of length 4 each.
    handles, labels = ax.get_legend_handles_labels()

    # When creating the legend, only use the first two elements
    # to effectively remove the last two.
    l = plt.legend(handles[0:2], labels[0:2])

    plt.show()

#!/bin/bash

# Import my extensions
import sys
sys.path.insert(0, '/home/sarah/Documents/Spirometry/python/extensions')

import csv
import copy
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
from dataStore import Datasets
import matplotlib.patches as mpatches
from matplotlib.lines import Line2D


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


#############
# Functions #
#############

def find_first_index_above_value(data, startSearchIndex, target):
    # first index grater than target in range data[startIndex:]
    for searchIndex in range(startSearchIndex, len(data)):
        if(data[searchIndex] > target):
            return searchIndex
    return len(data) -1


def find_first_index_below_value(data, startSearchIndex, target):
    # first index grater than target in range data[startIndex:]
    for searchIndex in range(startSearchIndex, len(data)):
        if(data[searchIndex] < target):
            return searchIndex
    return len(data) -1



def time_dependent_exp_fit(flow, time):
    # Exponential decay fit to dataSets given
    # NOTE: takes abs of flow, so any neg flow will be flipped
    ln_flw = [np.log(abs(f)) for f in flow]
    ones = [1]*len(time)

    dependent = np.array([ln_flw])
    independent = np.array([ones, time])
    try:
        res = lstsq(independent.T, dependent.T)
        decay = (res[0][1][0])
        offset = np.exp(res[0][0][0])
    except(ValueError):
        print('ValueError in exp decay fit: Data has nan?\nSetting result to NAN')
        decay = nan
        offset = nan
    return([decay, offset])



def identifyOutliers(x, outlierConstant):
    # Normally, an outlier is outside 1.5 * the IQR
    arr = np.array(x)
    upper_quartile = np.nanpercentile(arr, 75)
    lower_quartile = np.nanpercentile(arr, 25)
    arr.tolist()
    IQR = (upper_quartile - lower_quartile) * outlierConstant
    quartileSet = (lower_quartile - IQR, upper_quartile + IQR)

    resultList = []
    for y in range(len(arr)):
        if arr[y] <= quartileSet[0] or arr[y] >= quartileSet[1]:
            resultList.append(y)
    return resultList


#=============================================================


class Data:

    def __init__(self, filename):
        self.filename = filename
        self.dataLoc = []

        # For data
        self.flow = [0]
        self.volume = [0]
        self.time = [0]
        self.pressure = [0]
        self.breathPhase = [0]
        self.shutterState = [0]


    def __calcDataLoc(self, locations):
        dataLoc = []
        targets = ['Flow', 'Volume', 'TotalTime', 'BreathPhase', 'ShutterState', 'MouthPressureHighGain']

        for index in range(len(targets)):
            target = targets[index]

            locIndex = 0
            while locIndex < len(locations):
                if(locations[locIndex] == target):
                    dataLoc.append(locIndex)
                    locIndex = len(locations)
                locIndex += 1
        return dataLoc


    def extractData(self):

        with open(self.filename, 'rb') as csvfile:
           # reader = csv.reader(csvfile, delimiter = ',')
            reader = csv.reader(line.replace('\0','') for line in csvfile)

            # skip header line
            header = reader.next()
            self.dataLoc = self.__calcDataLoc(header)


            for row in reader:

                try:
                    dp_flow = -float(row[self.dataLoc[0]]) # in L/s
                except(ValueError):
                    dp_flow = np.nan
                try:
                    dp_pressure = float(row[self.dataLoc[5]]) # in L/s
                except(ValueError):
                    dp_pressure= np.nan
                try:
                    dp_volume = -float(row[self.dataLoc[1]]) # in L
                except(ValueError):
                    dp_volume = np.nan
                try:
                    dp_time = float(row[self.dataLoc[2]]) # in s
                except(ValueError):
                    dp_time= np.nan
                try:
                    dp_phase = float(row[self.dataLoc[3]]) # 1=in, 2=out
                except(ValueError):
                    dp_phase = np.nan
                try:
                    dp_shutter = float(row[self.dataLoc[4]]) # 50 or 1-0
                except(ValueError):
                    dp_shutter = np.nan

                self.flow.append(dp_flow)
                self.volume.append(dp_volume)
                self.time.append(dp_time)
                self.breathPhase.append(dp_phase)
                self.pressure.append(dp_pressure)
                # Doesn't exist every second breath, so interpolate ;)
                if np.isnan(dp_shutter):
                    self.shutterState.append(self.shutterState[-1])
                else:
                    self.shutterState.append(dp_shutter)


    def samplingFrequency(self):
        return(1/(self.time[3]-self.time[2]))


    def filterFlow(self):
        # Filter pressure
        # Filter high resolution pressure data
        fc = 220#filter cutoff
        bw = 40#filter bandwidth (taps?)
        flow = hamming(self.flow, fc, 300, bw, plot=False)
        self.flow = np.real(flow).tolist()
        volume = hamming(self.volume, fc, 300, bw, plot=False)
        self.volume = np.real(volume).tolist()




class Breaths():
    def __init__(self, fileData):
        self.fileData = fileData

        self.flow_splits = []
        self.flow_splits_end = []
        self.shutter_splits = []
        self.shutter_splits_start = []
        self.shutterIndex = 0

        self.breathNumber = 0
        self.flw = []
        self.vol = []
        self.t = []
        self.pres = []



    def __quick_breath_split(self, breath, cutoff, RISINGEDGE=True):
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



    def splitBreaths(self):
        # Split breaths apart
        # breath phase goes high at start of inspiration
        self.flow_splits = self.__quick_breath_split(self.fileData.breathPhase, 1.5)
        self.flow_splits.append(len(self.fileData.flow)-1)
        self.flow_splits_end = self.__quick_breath_split(self.fileData.breathPhase, 1.5, RISINGEDGE=False)
        self.flow_splits_end.append(len(self.fileData.flow)-1)

        #shutterState goes low at end of shuttering
        self.shutter_splits = self.__quick_breath_split(self.fileData.shutterState, 0.5, RISINGEDGE=False)
        self.shutter_splits.append(len(self.fileData.flow)-1)
        self.shutter_splits_start = self.__quick_breath_split(self.fileData.shutterState, 0.5)
        self.shutter_splits_start.append(len(self.fileData.flow)-1)


    def breathCount(self):
        return(len(self.flow_splits)-1)


    def extractSingleBreath(self, breathNumber):
        self.breathNumber = breathNumber

        self.flw = self.fileData.flow[self.flow_splits[breathNumber]:self.flow_splits_end[breathNumber]]

        self.vol = self.fileData.volume[self.flow_splits[breathNumber]:self.flow_splits_end[breathNumber]]
        self.vol = [v-self.vol[0] for v in self.vol] # shift volume curve down

        self.t = self.fileData.time[self.flow_splits[breathNumber]:self.flow_splits_end[breathNumber]]

        self.pres = self.fileData.pressure[self.flow_splits[breathNumber]:self.flow_splits_end[breathNumber]]


    def breathIsShuttered(self, breathNumber):
        #breath starts at start of inspiration
        breathStart =  self.flow_splits[breathNumber]
        breathMiddle = self.flow_splits_end[breathNumber]

        #breath ends when next breath starts, or data runs out (last pt of flow_split is final dp)
        if(breathNumber+1 <= self.breathCount()):
            breathEnd =  self.flow_splits[breathNumber+1]
        else:
            breathEnd =  self.flow_splits[-1]

        for shutterIndex in self.shutter_splits:
            # is a shutter started between breathing in and next breath in?
            if(shutterIndex > breathStart and shutterIndex < breathEnd):
#                print("This breath is shuttered")
                if(shutterIndex > (breathMiddle - (breathMiddle-breathStart)/4)):
#                    print("But shuttering started later than last quarter of expiration")
                    return False
                else:
                    self.shutterIndex = shutterIndex-breathStart # Save shutter start relative to breath start
                    return True
            # if shutter started beyond end-exp of this breath
            elif(shutterIndex > breathEnd):
                return False
        return False



class averageBreaths():
    def __init__(self):
        self.data = []

    def averageData(self):
        numberToAverage = len(self.data)
        shortestdata = min([len(array) for array in self.data])
        averagedata = [0]*shortestdata

        for index in range(shortestdata):
            averagedata[index] = (np.sum([array[index] for array in self.data]))/numberToAverage
        return averagedata


    def clean(self):
        self.data = []


class QVLoopAnalyser():
    def __init__(self, flow, volume, shutterStartPoint):
        self.flw = flow
        self.vol = volume
        self.shutterStartPoint = shutterStartPoint

        self.flwFilt = []
        self.volFilt = []

        self.maxFlwIndex = 0
        self.maxFlwFiltIndex = 0

        self.negFlowIndex = 0

        self.endFlowIndex = 0
        self.expEndIndex = 0

        self.gradient = 0
        self.intercept = 0
        self.rvalue = 0


    def __flowInflectionPoint(self):
        searchRange = 3
        # Has the flow increased again in the future?
        # If so, inflection in flow (ideal: it would always decrease)
        searchEnd = find_first_index_below_value(self.flwFilt, self.maxFlwFiltIndex, 0)
        for index in range(self.maxFlwFiltIndex, searchEnd-searchRange):
            if (self.flwFilt[index+searchRange] >= self.flwFilt[index]):
                return(index)
            elif (self.flwFilt[index] <= 0):
                return(index)
        return(searchEnd)


    def __volumeStagnationPoint(self):
        searchRange = 8
        # Has the volume stagnated or started decreasing
        # expect large-ish increase every dp
        # call stagnation when volume increased less than a minimum amount
        for index in range(self.maxFlwIndex, len(self.volFilt)-searchRange):
            if (self.volFilt[index+searchRange] < self.volFilt[index]+0.005): # 0.005ml change at 0.15L/s
                return(index)
        return(len(self.vol)-1)


    def calcLinearRange(self):
        # Filter flow to find inflection points easier
        fc = 50#filter cutoff
        bw = 20#filter bandwidth (taps?)
        AminusBFilt = hamming(self.flw, fc, 300, bw, plot=False)
        self.flwFilt = np.real(AminusBFilt).tolist()
        self.volFilt = integral(self.flwFilt, FS)


        # find max shutter flow, and the end point in volume
        self.maxFlwIndex = self.flw.index(max(self.flw[min(self.shutterStartPoint,len(self.flw)-2):min(self.shutterStartPoint+100,len(self.flw)-1)]))
        self.maxFlwFiltIndex = self.flwFilt.index(max(self.flwFilt[min(self.shutterStartPoint,len(self.flwFilt)-2):min(self.shutterStartPoint+100,len(self.flwFilt)-1)]))


        # find end as function of (Q increase, V stagnation)
        flowInflection = self.__flowInflectionPoint()
        volumeStagnation = self.__volumeStagnationPoint()

        # Find any point flow goes negative
        self.negFlowIndex = find_first_index_below_value(self.flw, self.maxFlwFiltIndex, 0)

#        print('Shutter range end at first of flow inflection and volume stagnation: {} and {}'.format(flowInflection, volumeStagnation))
        #self.endFlowIndex = min(flowInflection, volumeStagnation)
        self.endFlowIndex = flowInflection

        self.expEndIndex = self.endFlowIndex + (len(self.flw)-self.endFlowIndex)*1/3
        if(self.flw[self.expEndIndex] <= 0):
            self.expEndIndex = find_first_index_below_value(self.flw, self.endFlowIndex, 0)
#        print('Full expiration end point defined as 75% of data after end shutter range: {}'.format(self.expEndIndex))



    def calcDecayRate(self, start, end):
#        print(start)
#        print(end)
#        print(self.shutterStartPoint)
        # Calc gradient of QV line
        if(end-start > 3 and self.shutterStartPoint > 0): #min 3dp needed for least-squares
            # fit to full identified linear length
            volume_curve = self.vol[start:end]
            flow_curve = self.flw[start:end]
            gradient_passive, intercept, r_value, p_value, std_err = stats.linregress(volume_curve, flow_curve)

#            print('New linear fit made.')
#            print('Local gradient, intercept and r-value overwritten.')
            self.gradient = gradient_passive
            self.intercept = intercept
            self.rvalue = r_value
#            print('Gradient and intercept: {} and {}'.format(gradient_passive, intercept))

        else:
#            print('Unable to make linear fit. Not enough data points.')
#            print('Local gradient, intercept and r-value overwrittenas nan.')
            self.gradient = np.nan
            self.intercept = np.nan
            self.rvalue = np.nan


#            lineEnd = [intercept + gradient_passive*m for m in [self.vol[self.expEnd], self.vol[self.maxFlowIndex]]]
#            lineFullEnd = [intercept + gradient_passive*m for m in [min(AminusBVol), max(AminusBVol)]]
#            print('QV gradient End: {}'.format(gradient_passive))
#            print('r value End: {}'.format(r_value))
#
#
#                ########################
#                # EXCEPTIONS: SET NAN  #
#                ########################
#                if(maxFlowIndex > (shutterStartPoint + len(AminusB[shutterStartPoint:])/3)):
#                    print('**\nSecond peak identified more than 1/3 data length after trueMaxIndex. Setting result to NAN\n')
#                    if(gradients[i]):
#                        gradients[i][-1] = np.nan
#                    else:
#                        gradients[i].append(np.nan)

    def generateLines(self):
        line = [self.intercept + self.gradient*m for m in [self.vol[self.endFlowIndex], self.vol[self.maxFlwIndex]]]
        lineFull = [self.intercept + self.gradient*m for m in [self.vol[self.expEndIndex], self.vol[self.maxFlwIndex]]]
        return([line, lineFull])


class mechanicsAnalyser():
    def __init__(self):
        self.pressures = []
        self.flows = []
        self.Rxs = []
        self.timestep = 0


    def setTimestep(self, time):
        self.timestep = int((time[3]-time[2])*1000 + 0.5) #in ms (swedish rounding)



    def __calcShutterRange(self, pressure):
        # The exact rising and falling edge of the shutter
        start = find_first_index_above_value(pressure, 0, 100)
        end = find_first_index_below_value(pressure, start, 100)
        end = pressure.index(max(pressure[start:end]))
        return([start, end])


    def calcROCC(self, flow, pressure):
        # pressure 30-75ms after shutter back extrapolated to 15ms

        ms15 = int(15/self.timestep)
        ms30 = int(30/self.timestep)
        ms75 = int(75/self.timestep)

        shutterClose = self.__calcShutterRange(pressure)[0]

        grad = (pressure[shutterClose+ms75] - pressure[shutterClose+ms30])/(ms75-ms30)
        offs = pressure[shutterClose+ms30]

        p15 = grad*(ms15-ms30) + offs

        # flow 15ms before shutter
        # filter flow for smoother fit in resistance added data
        #fc = 120#filter cutoff
        #bw = 30#filter bandwidth (taps?)
        #filtFlow = hamming(flow, fc, 300, bw, plot=False)
        #filtFlow = np.real(filtFlow).tolist()

        if(0):
            plt.figure()
            plt.plot(flow)
            plt.plot(filtFlow)
            plt.show()

        qminus15 = flow[shutterClose-ms15]
        pminus15 = pressure[shutterClose-ms15]

        # Rocc as change in pressure due to shutter/flow pre shutter
        Rocc = (p15 - pminus15)/qminus15
       # print('Rocc: {}'.format(Rocc))

        if(1):
            print()
            plt.rc('legend',**{'fontsize':12})
            z, (zx1, zx2) = plt.subplots(2, sharex=False)

            zx1.plot(flow)
            zx1.plot(shutterClose-ms15, qminus15, 'bo')

            zx2.plot(pressure)
            zx2.plot(shutterClose+ms30, pressure[shutterClose+ms30], 'rx')
            zx2.plot(shutterClose+ms75, pressure[shutterClose+ms75], 'rx')
            zx2.plot(shutterClose-ms15, p15 - pressure[shutterClose-ms15], 'bo')
            plt.show()

        return(Rocc/100.0)



    def calcRelaxationGradient(self, pressure):
        # Rough estimation of gradient only using start and end points
        ms30 = int(30/self.timestep)
        start, end = self.__calcShutterRange(pressure)

        # Get in away from edges
        start += ms30
        end -= ms30

        grad = (pressure[end]-pressure[start])/(end-start)
        return grad


    def calcMechanicsFromDecay(self, decay, endShutterPressure, shutterVolume):

        divTau = decay/endShutterPressure
        Rrs = 1/(shutterVolume*divTau)
        Ed = decay*Rrs
        Ed = endShutterPressure/shutterVolume
        return(Rrs, Ed)


    def calcStartShutterElastance(self, startShutterPressure, shutterVolume):

        Ed = startShutterPressure/shutterVolume
        return(Ed)



#=============================================================
#=============================================================

#for plotting
colours = ['blue','cyan', 'green','red','black','pink','purple','magenta', 'orange', 'brown', 'yellow','gray','olive','lavender','gold','deeppink','orchid','cadetblue','palegreen','blue','cyan', 'green','red','black','pink','purple','magenta', 'orange', 'brown', 'yellow','gray','olive','lavender','gold','deeppink','orchid','cadetblue','palegreen']
symbols = ['o','o','o','o','o','o','o','o','o','o','o','o','o','o','o','o','o','o','o','d','d','d','d','d','d','d','d','d','d','d','d','d','d','d','d','d','d','d',]
#colours = ['blue', 'green', 'red', 'black', 'purple', 'cyan', 'magenta', 'brown', 'orange', 'pink', 'yellow', 'blue', 'green', 'red', 'black', 'purple', 'cyan', 'magenta']
colourIndex = 0

megaRelaxGrad = []
megaTau = []
megaRocc = []
megaE = []
megaPpk = []
megaPstart = []
megaVol = []

#==================
# Debug etc plots
STACK = 0
FULLSET = 0
STACKSHUTTER = 0
STACKPRESSURE = 1
BOXPLOTS = 0

if(STACK):
    # plot stacks of flow up to and including shutter
    plt.rc('legend',**{'fontsize':12})
    a, (ax1) = plt.subplots(1, sharex=False)

if(STACKPRESSURE):
    # plot flow separated into components
    plt.rc('legend',**{'fontsize':12})
    g, (gax1, gax2) = plt.subplots(2, sharex=True)
    u, (uax1) = plt.subplots(1)
    #h, (hax1) = plt.subplots(1)

if(STACKSHUTTER):
    allShutterFlw= []
    allShutterPressure = []

    plt.rc('legend',**{'fontsize':12})
    q, (qx1, qx2) = plt.subplots(2, sharex=False)
    w, (wx1) = plt.subplots(1)

#==================
#Results tables. Check wrappers of plots are enabled too
GenRoccTable = 0
if(GenRoccTable):
    RoccMeans = []
    RoccStds = []
    RoccExpected = []
    RoccError = []
    tableIDs = []

#=============================================================
dataNums = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17] #150]
#dataNums = [ 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 14, 15, 16, 17] #150]
#dataNums = [11, 12, 13] #150]
#dataNums = [9]
for DATASET in dataNums:
    dataInfo = Datasets(DATASET)
    dataID = 'V'+ str(DATASET)
    numFiles= len(dataInfo.files)

    #=============================================================
    # For results
    gradients = [[] for m in range(numFiles)]
    gradientsEnd = [[] for m in range(numFiles)]

    decays = [[] for m in range(numFiles)]
    prev_max_pres = [[] for m in range(numFiles)] # store maximum pressure in shuttered breaths
    shutterVol = [[] for m in range(numFiles)] # store maximum pressure in shuttered breaths
    shutterStarts = [[] for m in range(numFiles)] # store maximum pressure in shuttered breaths
    Rguess = [[] for m in range(numFiles)] # store maximum pressure in shuttered breaths


    # Mechanics analyser stuff
    allAverageFlow = []
    allGradients = []
    allGradientsEnd = []
    allRocc = []
    allRx = []
    allE = []
    allRrs = []
    allShutterVolume= []
    allRelaxGrad = []
    allPeakAvgFlow = []
    allPeakShutterFlow = []
    allPeakAminusB = []
    allStartShutterPressure = []
    allEndShutterPressure = []
    allStartShutterIndex = []


    #=============================================================
    measuredMechanics = mechanicsAnalyser()

    #numFiles = 2 # uncomment for only looking at no_Rx data
    for i in range(numFiles):
        #colourIndex += 1

        #####################
        # GET DATA FROM FILE#
        #####################

        filename = dataInfo.path + dataInfo.files[i]
        print("/n" + filename)

        fileData = Data(filename)
        fileData.extractData()

        RxLabel = dataInfo.xLabels[i]
        Rx = np.nan
        if RxLabel == 'Normal':
            Rx = 0
        elif RxLabel == '12.5mm':
            Rx = 0.4
        elif RxLabel =='10.5mm':
            Rx = 0.8
        elif RxLabel =='9.5mm':
            Rx =  1.2
        else:
            Rx = np.nan

        FS = fileData.samplingFrequency()

        fileBreaths = Breaths(fileData)
        fileBreaths.splitBreaths()


        averages = averageBreaths()
        averagePressures = averageBreaths()



        if(FULLSET):
            # Plot flow, volume, points of interest
            plt.rc('legend',**{'fontsize':12})
            q, (qax1, qax2, qax3) = plt.subplots(3, sharex=True)

            qax1.plot(fileData.flow)
            qax1.plot(fileBreaths.flow_splits, [fileData.flow[m] for m in fileBreaths.flow_splits], 'or')
            qax1.plot(fileBreaths.flow_splits_end, [fileData.flow[m] for m in fileBreaths.flow_splits_end], 'ok')
            qax1.legend(['flow', 'expiratory starts', 'ends'])

            qax2.plot(fileData.volume)
            qax2.plot(fileBreaths.flow_splits, [fileData.volume[m] for m in fileBreaths.flow_splits], 'or')
            qax2.legend(['Volume and expiratory starts'])

            qax3.plot(fileData.pressure)
            qax3.plot(fileData.shutterState)
            qax3.plot(fileBreaths.shutter_splits, [fileData.flow[m] for m in fileBreaths.shutter_splits], 'om')
            qax3.plot(fileBreaths.shutter_splits_start, [fileData.flow[m] for m in fileBreaths.shutter_splits_start], 'oy')
            qax3.legend(['Flow and shutter ends', 'shutter'])
            q.show()




        ANALYSING = False

        for breath in range(fileBreaths.breathCount()-1):
#            print('\nbreath: {}'.format(breath))
            fileBreaths.extractSingleBreath(breath)
            measuredMechanics.setTimestep(fileBreaths.t)

            if(fileBreaths.breathIsShuttered(breath)):
#                print('This breath is shuttered')
                ANALYSING = True



    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~``
            ################################################################
            # Add flow data to averaging array when not a shuttered breath #
            # We are averaging the QV loop for the subject                 #
            ################################################################
            if not ANALYSING:
                averages.data.append(fileBreaths.flw)


                if(STACK):
                    ax1.plot(fileBreaths.flw, color=colours[colourIndex])
                    colourIndex += 1

                if(STACKPRESSURE):
                    gax1.plot(fileBreaths.flw, color='blue')
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
                # Work out the average flow from prev breaths
                averageFlow = averages.averageData()
                averageVol = integral(averageFlow, FS)


                # Where is the shutter in the breath?
                shutterStartPoint = fileBreaths.shutterIndex

                # Find rising edge of shutter
                ms30 = int(30/measuredMechanics.timestep)
                startIndex = find_first_index_above_value(fileBreaths.pres, 0, 100) + ms30

                # find falling edge of shutter
                endIndex = find_first_index_below_value(fileBreaths.pres[startIndex:], 0, 100) + startIndex - int(1.5*ms30)
                peakPressureIndex = fileBreaths.pres.index(max(fileBreaths.pres))
                endIndex = max(peakPressureIndex, endIndex)

                if(0):
                    plt.figure()
                    plt.plot(fileBreaths.pres)
                    plt.plot(startIndex, fileBreaths.pres[startIndex], 'or')
                    plt.plot(endIndex, fileBreaths.pres[endIndex], 'og')
                    plt.show()


                # Estimate flow due to shutter
                AminusBTrue = [fileBreaths.flw[k]-averageFlow[k] for k in range(min(len(fileBreaths.flw),len(averageFlow)))]
                AminusBTrue = fileBreaths.flw

                # Volumes moved due to shutter
                AminusBVolTrue = integral(AminusBTrue, FS)



                # Determine linear range of post shutter flow
                measuredQVAnalysis = QVLoopAnalyser(AminusBTrue, AminusBVolTrue, shutterStartPoint)
                measuredQVAnalysis.calcLinearRange()


                # Gradient of shuttered region
                measuredQVAnalysis.calcDecayRate(measuredQVAnalysis.maxFlwIndex, measuredQVAnalysis.endFlowIndex)
                line, lineFull = measuredQVAnalysis.generateLines()
                gradients[i].append(measuredQVAnalysis.gradient)
                allGradients.append(measuredQVAnalysis.gradient)

                # Gradient of post-shuttered region
                # post shuter region either goes below
                if(abs(measuredQVAnalysis.expEndIndex - measuredQVAnalysis.negFlowIndex) > 10):
                    measuredQVAnalysis.calcDecayRate(measuredQVAnalysis.endFlowIndex, measuredQVAnalysis.expEndIndex)
                else:
                    measuredQVAnalysis.calcDecayRate(min(measuredQVAnalysis.negFlowIndex+15, len(fileBreaths.flw))-18, len(fileBreaths.flw) - 15)

                lineEnd, lineFullEnd = measuredQVAnalysis.generateLines()
                gradientsEnd[i].append(measuredQVAnalysis.gradient)
                allGradientsEnd.append(measuredQVAnalysis.gradient)

                # Save every shuttered breath from shutter start
                averagePressures.data.append(fileBreaths.pres)
                #allAverageFlow.append(averageFlow)
                allAverageFlow.append(fileBreaths.flw)

                # Calculate the relaxation gradient of pressure after initial rise
                relaxGrad = measuredMechanics.calcRelaxationGradient(fileBreaths.pres)
                peakAvgFlow = max(averageFlow)
                peakShutterFlow = max(fileBreaths.flw)
                peakAminusB = max(AminusBTrue)


                # Calc Rocc for this breath
                Rocc = measuredMechanics.calcROCC(fileBreaths.flw, fileBreaths.pres)

                # Calc mechanics from decay rate
                endShutterPressure = fileBreaths.pres[endIndex]/100
                startShutterPressure = fileBreaths.pres[startIndex]/100
                shutterVolume = fileBreaths.vol[-1] - fileBreaths.vol[startIndex]
                Rrs, Ed = measuredMechanics.calcMechanicsFromDecay(gradients[i][-1], endShutterPressure, shutterVolume)
                Estart = measuredMechanics.calcStartShutterElastance(startShutterPressure, shutterVolume)



                # Save stuff
                allRelaxGrad.append(relaxGrad)
                allPeakAvgFlow.append(peakAvgFlow)
                allPeakShutterFlow.append(peakShutterFlow)
                allPeakAminusB.append(peakAminusB)
                allStartShutterIndex.append(startIndex)
                allStartShutterPressure.append(fileBreaths.pres[startIndex]/100)
                allEndShutterPressure.append(endShutterPressure)
                allRocc.append(Rocc)
                allRx.append(Rx)
                allShutterVolume.append(shutterVolume)
                allE.append(Ed)
                allRrs.append(Rrs)


                # Plot stuff
                ########################
                if(STACK):
                    ax1.plot(averageFlow, 'purple')
                    ax1.plot(fileBreaths.flw, 'r', linewidth=3)
                    ax1.grid()
                    ax1.set_ylabel('Flow', fontsize=30)
                    a.show()
                    plt.rc('legend',**{'fontsize':12})
                    a, (ax1) = plt.subplots(1, sharex=False)


                if(STACKSHUTTER):
                    pass


                if(STACKPRESSURE):
                    # QV loops
                    minVol = min(AminusBVolTrue)
                    newAminusBVol = [v-minVol for v in AminusBVolTrue]
                    uax1.plot(newAminusBVol, AminusBTrue, linewidth = 3)
                    uax1.plot(newAminusBVol, AminusBTrue, ':', linewidth = 3)
                    uax1.plot(measuredQVAnalysis.volFilt, measuredQVAnalysis.flwFilt, '.-', linewidth = 3)
                    uax1.plot(newAminusBVol[measuredQVAnalysis.maxFlwIndex], AminusBTrue[measuredQVAnalysis.maxFlwIndex], 'o', color='yellow')
                    uax1.plot(newAminusBVol[measuredQVAnalysis.endFlowIndex], AminusBTrue[measuredQVAnalysis.endFlowIndex], 'o', color='yellow')
                    uax1.plot(newAminusBVol[measuredQVAnalysis.expEndIndex], AminusBTrue[measuredQVAnalysis.expEndIndex], 'o', color='yellow')
                    uax1.plot(newAminusBVol[measuredQVAnalysis.negFlowIndex], AminusBTrue[measuredQVAnalysis.negFlowIndex], 'o', color='cyan')

                    plt.figure()
                    plt.plot(averageVol, averageFlow, '-.', color='Red', linewidth = 3)
                    plt.plot([averageVol[measuredQVAnalysis.endFlowIndex], averageVol[measuredQVAnalysis.expEndIndex]],[averageFlow[measuredQVAnalysis.endFlowIndex], averageFlow[measuredQVAnalysis.expEndIndex]] )
                    plt.plot(averageVol[measuredQVAnalysis.maxFlwIndex], averageFlow[measuredQVAnalysis.maxFlwIndex], 'o', color='yellow')
                    plt.plot(averageVol[measuredQVAnalysis.endFlowIndex], averageFlow[measuredQVAnalysis.endFlowIndex], 'o', color='yellow')
                    plt.plot(averageVol[measuredQVAnalysis.expEndIndex], averageFlow[measuredQVAnalysis.expEndIndex], 'o', color='yellow')
                    plt.ylabel('Flow (L/s)', fontsize=32)
                    plt.xlabel('Volume (L)', fontsize=32)
                    plt.title('average QV loop')

                    if(not np.isnan(gradients[i][-1])):
                        uax1.plot([newAminusBVol[measuredQVAnalysis.expEndIndex], newAminusBVol[measuredQVAnalysis.maxFlwIndex]], lineFull, color='purple', linewidth=2)
                        uax1.plot([newAminusBVol[measuredQVAnalysis.expEndIndex], newAminusBVol[measuredQVAnalysis.maxFlwIndex]], lineFullEnd, color='purple', linewidth=2)

                    uax1.set_ylabel('Flow', fontsize=30)
                    uax1.set_xlabel('Volume', fontsize=30)
                    uax1.grid()


                    # flow and pressure plots
                    gax1.plot(averageFlow[:len(fileBreaths.flw)], color='purple', linewidth=1.5, label='Avgerage tidal')
                    gax1.plot(fileBreaths.flw, color='red', linewidth=3, label='During shuttering')
                    gax1.plot(AminusBTrue, '--', color='black', linewidth=3, label='Induced by shutter')

                    gax1.plot(measuredQVAnalysis.maxFlwIndex, fileBreaths.flw[measuredQVAnalysis.maxFlwIndex], 'o', color='yellow')
                    gax1.plot(measuredQVAnalysis.endFlowIndex, fileBreaths.flw[measuredQVAnalysis.endFlowIndex], 'o', color='yellow')
                    gax1.plot(measuredQVAnalysis.expEndIndex, fileBreaths.flw[measuredQVAnalysis.expEndIndex], 'o', color='yellow')

                    gax1.grid()
                    gax1.set_ylabel('Flow', fontsize=30)
                    gax1.legend()


                    gax2.plot(fileBreaths.pres, color='red', linewidth=3)
                    gax2.grid()
                    gax2.set_xlabel('Datapoint', fontsize=30)
                    gax2.set_ylabel('Pressure', fontsize=30)


                    plt.show()
                    g.show()
                    u.show()
                    plt.rc('legend',**{'fontsize':12})
                    g, (gax1, gax2) = plt.subplots(2, sharex=True)
                    u, (uax1) = plt.subplots(1)
                ########################



                averages.clean()
                ANALYSING = False

                #Print results table
#                print(gradients)

        measuredMechanics.pressures.append(list(averagePressures.data))
        measuredMechanics.flows.append(list(allAverageFlow))


#        if(STACKSHUTTER):
#            # for aceh test want to look at all the shutters
#            # have average flow and the shutter pressure for each shutter
#            average_pressure = averagePressures.averageData()
#
#            plt.figure(101)
#            plt.plot(average_pressure)
#
#            # Quick least squares between average pressure and average flow
#            # No good because of viscoelastic effect
#            for e in range(len(allShutterPressure)):
#               Pav = allShutterPressure[e]
#
#               flatStart = find_first_index_above_value(Pav, 0, 0.8*max(Pav))
#               flatEnd = Pav.index(max(Pav))
#
#               Pav = allShutterPressure[e][flatStart:flatEnd]
#               Fav = allShutterFlw[e][flatStart:flatEnd]
#
#              # plt.plot(allShutterPressure[e])
#              # plt.plot(range(flatStart,flatEnd), Pav)
#
#               dependent = np.array([Pav])
#               independent = np.array([Fav])
#               res = lstsq(independent.T, dependent.T)
#               gradient = res[0][0][0]
#               print(gradient)


        # Tidy up :D
        averagePressures.clean()
        allAverageFlow = []




    #=============================================================
    #=============================================================
    #=============================================================
    # Final analysis
    measuredMechanics.Rxs = dataInfo.xLabels
    #measuredMechanics.calcREFF()


    def extractRxArrays(allArray, allRx):

        Rxarray = [0, 0.4, 0.8, 1.2]
        arrayNone = []
        array12_5 = []
        array10_5 = []
        array9_5 = []
        for index in range(len(allArray)):
            if allRx[index] == Rxarray[0]:
                arrayNone.append(allArray[index])
            elif allRx[index] == Rxarray[1]:
                array12_5.append(allArray[index])
            elif allRx[index] == Rxarray[2]:
                array10_5.append(allArray[index])
            elif allRx[index] == Rxarray[3]:
                array9_5.append(allArray[index])
        return ([arrayNone, array12_5, array10_5, array9_5])



    gradsFullNone, gradsFull12_5, gradsFull10_5, gradsFull9_5 = extractRxArrays(allGradients, allRx)
    RoccArray = extractRxArrays(allRocc, allRx)
    volArray = extractRxArrays(allShutterVolume, allRx)
    pressureArray = extractRxArrays(allEndShutterPressure, allRx)

    gradsNone = []
    grads12_5 = []
    grads10_5 = []
    grads9_5 = []
    gradArray = [gradsNone, grads12_5, grads10_5, grads9_5]


    j = 0
    start = 0
    for grad in [gradsFullNone, gradsFull12_5, gradsFull10_5, gradsFull9_5]:
        end = len(grad)+start
        outlierIndices = identifyOutliers(grad, 1.5)
        print(outlierIndices)

        #calc mean wo outliers
        for thing in np.delete(grad, outlierIndices):
            gradArray[j].append(thing)

        volArray[j] = np.delete(volArray[j], outlierIndices).tolist()
        pressureArray[j] = np.delete(pressureArray[j], outlierIndices).tolist()

        j += 1



    RNone = [0 for x in range(len(gradsNone))]
    R12_5 = [0 for x in range(len(grads12_5))]
    R10_5 = [0 for x in range(len(grads10_5))]
    R9_5 = [0 for x in range(len(grads9_5))]
    ENone = [0 for x in range(len(gradsNone))]
    E12_5 = [0 for x in range(len(grads12_5))]
    E10_5 = [0 for x in range(len(grads10_5))]
    E9_5 = [0 for x in range(len(grads9_5))]
    RArray = [RNone, R12_5, R10_5, R9_5]
    EArray = [ENone, E12_5, E10_5, E9_5]

    j = 0
    for grad in gradArray:
        divTau = [0 for x in range(len(grad))]
        divVol = [0 for x in range(len(grad))]
        for k in range(len(grad)):
            divTau[k] = grad[k]/(pressureArray[j][k])
            divVol[k] = 1/volArray[j][k]

        # E and R estimates
        Rd = [0 for x in range(len(grad))]
        Ed = [0 for x in range(len(grad))]
        for k in range(len(grad)):
            Rd[k] = 1/(volArray[j][k]*divTau[k])
            Ed[k] = grad[k]*(Rd[k])
        RArray[j] = Rd
        EArray[j] = Ed
        j += 1

#    # for plotting all E vs all P0 for all subjects
#    for j in range(len(EArray)):
#        for v in (EArray[j][:]):
#            megaE.append(v)





    # Final table row printouts
    def MeansAndSTDevs(data):
        res = [[0]*2 for j in range(4)]
        j = 0
        for grad in data:
            mean = np.nanmean(grad)
            std = np.nanstd(grad)
            res[j][0] = mean
            res[j][1] = std

            j+=1
        return res


    # Tau: for all Rx, tau mean +/- std
    # where mean and std ignore outliers (values > 3*std)
    if(1):
        res = MeansAndSTDevs([gradsNone, grads12_5, grads10_5, grads9_5])
        print('Tau data')
        print(' & {:.2f} [{:.2f}] & {:.2f} [{:.2f}] & {:.2f} [{:.2f}] & {:.2f} [{:.2f}] \ \ '.format(res[0][0], res[0][1], res[1][0], res[1][1], res[2][0], res[2][1], res[3][0], res[3][1] ))

    # for mechanics
    if(1):
        res = MeansAndSTDevs(RArray)
        print('R data')
        print(' & {:.2f} [{:.2f}] & {:.2f} [{:.2f}] & {:.2f} [{:.2f}] & {:.2f} [{:.2f}] \ \ '.format(res[0][0], res[0][1], res[1][0], res[1][1], res[2][0], res[2][1], res[3][0], res[3][1] ))

        res = MeansAndSTDevs(EArray)
        print('E data')
        print(' & {:.2f} [{:.2f}] & {:.2f} [{:.2f}] & {:.2f} [{:.2f}] & {:.2f} [{:.2f}] \ \ '.format(res[0][0], res[0][1], res[1][0], res[1][1], res[2][0], res[2][1], res[3][0], res[3][1] ))

        res = MeansAndSTDevs(RoccArray)
        print('Rocc data')
        print(' & {:.2f} [{:.2f}] & {:.2f} [{:.2f}] & {:.2f} [{:.2f}] & {:.2f} [{:.2f}] \ \ '.format(res[0][0], res[0][1], res[1][0], res[1][1], res[2][0], res[2][1], res[3][0], res[3][1] ))


    #==========
    #==========
    # PLOTS PLOTS PLOTS
    #==========
    #==========
    # Here looking at time-wise pressure 30ms after shutter
    TimewisePressure = 0
    if(TimewisePressure):
        start=0
        plt.figure()
        plt.grid()

        # shutter pressure vs time
        for index in range(len(gradients)):

            plt.xlabel('Index (5ms div)')
            plt.ylabel('Shutter pressure')

            plt.plot(allStartShutterIndex[start:start+len(gradients[index])], allStartShutterPressure[start:start+len(gradients[index])], 'o')
            start+=len(gradients[index])-1


        plt.show()

    #==========
    #RELAXATION GRADIENT m
    # Looking at how relaxation gradient of shutter pressure is affected by effort
    # (Here using peak flow as surrogate for driving effort)
    RelaxationGradVsFlow = 0
    if(RelaxationGradVsFlow):

        rvalues = [0 for v in range(6)]
        index = 0
        plt.figure(301)
        plt.ylabel('Relaxation gradient', fontsize=20)
        plt.xlabel('Lung volume at shutter close (L)', fontsize=20)
        plt.plot(allShutterVolume, allRelaxGrad, ls='',marker=symbols[colourIndex], color=colours[colourIndex])

        gradient, intercept, r_value, p_value, std_err = stats.linregress(allShutterVolume, allRelaxGrad)
        line = [0*gradient + intercept, 2.5*gradient + intercept]
        plt.plot([0, 2.5], line, '--', color='lightgray')
        rvalues[index] = (r_value**2)
        index += 1


        plt.figure(303)
        plt.ylabel('Relaxation gradient', fontsize=20)
        plt.xlabel('Pressure 30ms after shutter start (cmH2O)', fontsize=20)
        plt.plot(allStartShutterPressure, allRelaxGrad, ls='',marker=symbols[colourIndex], color=colours[colourIndex])

        gradient, intercept, r_value, p_value, std_err = stats.linregress(allStartShutterPressure, allRelaxGrad)
        line = [0*gradient + intercept, 10*gradient + intercept]
        plt.plot([0, 10], line, '--', color='lightgray')
        rvalues[index] = (r_value**2)
        index += 1


        plt.figure(304)
        plt.ylabel('Relaxation gradient', fontsize=20)
        plt.xlabel('Peak pressure (cmH2O)', fontsize=20)
        plt.plot(allEndShutterPressure, allRelaxGrad, ls='',marker=symbols[colourIndex], color=colours[colourIndex])

        gradient, intercept, r_value, p_value, std_err = stats.linregress(allEndShutterPressure, allRelaxGrad)
        line = [0*gradient + intercept, 20*gradient + intercept]
        plt.plot([0, 20], line, '--', color='lightgray')
        rvalues[index] = (r_value**2)
        index += 1


#        plt.figure(305)
#        plt.ylabel('Relaxation gradient', fontsize=20)
#        plt.xlabel('Peak flow during shuttering (absolute, L/s)', fontsize=20)
#        plt.plot(allPeakShutterFlow, allRelaxGrad, ls='',marker=symbols[colourIndex], color=colours[colourIndex])
#
#        gradient, intercept, r_value, p_value, std_err = stats.linregress(allPeakShutterFlow, allRelaxGrad)
#        line = [0*gradient + intercept, 4.5*gradient + intercept]
#        plt.plot([0, 4.5], line, color=colours[colourIndex])
#        print('vol: {}'.format(r_value))
#        plt.grid()
#

#        plt.figure(306)
#        plt.ylabel('Relaxation gradient', fontsize=20)
#        plt.xlabel('Peak flow caused by shutter only (shutter-averageTidal)', fontsize=20)
#        plt.plot(allPeakAminusB, allRelaxGrad, ls='',marker=symbols[colourIndex], color=colours[colourIndex])
#
#        gradient, intercept, r_value, p_value, std_err = stats.linregress(allPeakAminusB, allRelaxGrad)
#        line = [0*gradient + intercept, 2.5*gradient + intercept]
#        plt.plot([0, 2.5], line, color=colours[colourIndex])
#        print('peakQ: {}'.format(r_value))
#        plt.grid()
#

        plt.figure(307)
        plt.ylabel('Relaxation gradient', fontsize=20)
        plt.xlabel('Flow decay rate', fontsize=20)
        plt.plot(allGradients, allRelaxGrad, ls='',marker=symbols[colourIndex], color=colours[colourIndex])

        gradient, intercept, r_value, p_value, std_err = stats.linregress(allGradients, allRelaxGrad)
        line = [0*gradient + intercept, -60*gradient + intercept]
        plt.plot([0, -60], line, '--', color='lightgray')
        rvalues[index] = (r_value**2)
        index += 1


        plt.figure(308)
        plt.ylabel('Relaxation gradient', fontsize=20)
        plt.xlabel('Elastance (cmH2O/L)', fontsize=20)
        plt.plot(allE, allRelaxGrad, ls='',marker=symbols[colourIndex], color=colours[colourIndex])

        gradient, intercept, r_value, p_value, std_err = stats.linregress(allE, allRelaxGrad)
        line = [0*gradient + intercept, 25*gradient + intercept]
        plt.plot([0, 25], line, '--', color='lightgray')
        rvalues[index] = (r_value**2)
        index += 1


        plt.figure(302)
        plt.ylabel('Relaxation gradient', fontsize=20)
        plt.xlabel('Rocc (cmH2Os/L)', fontsize=20)
        plt.plot(allRocc, allRelaxGrad, ls='',marker=symbols[colourIndex], color=colours[colourIndex])

        gradient, intercept, r_value, p_value, std_err = stats.linregress(allRocc, allRelaxGrad)
        line = [0*gradient + intercept, 10*gradient + intercept]
        plt.plot([0, 10], line, '--', color='lightgray')
        rvalues[index] = (r_value**2)
        index += 1

        colourIndex += 1

        for value in rvalues:
            print("& {:.2f}".format(value)),
        print(" \ \ \n")

    for j in range(len(allRelaxGrad)):
        megaRelaxGrad.append(allRelaxGrad[j])
        megaTau.append(allGradients[j])
        megaRocc.append(allRocc[j])
        megaE.append(allE[j])
        megaPpk.append(allEndShutterPressure[j])
        megaPstart.append(allStartShutterPressure[j])
        megaVol.append(allShutterVolume[j])



    #==========
    # Results lines for E, Rocc, Rrs
    customLines = [Line2D([0], [0], color='k', linestyle='-', marker='o'),
                   Line2D([0], [0], color='r', linestyle='-', marker='') ]


    mechLines = 1
    if(mechLines):
        xarray = ['0', '0.4', '0.8', '1.2']

        #tau
        plt.figure(0)
        res = MeansAndSTDevs(gradArray)
        line = [0 for v in range(4)]
        for index in range(4):
            line[index] = res[index][0] # grab mean value
        normaliser = np.copy(line[0]) # make all lines start at 0
        for index in range(4):
            line[index] -= normaliser
        plt.plot(line, '-', color='grey')
        plt.plot(line, 'o', color='black')
        plt.xticks(range(4), xarray)
        plt.xlabel('Added resistance (cmH2Os/L)', fontsize=30)
        plt.ylabel('Decay rate (1/s)', fontsize=30)
        plt.legend(customLines, ['Data'])

        #E
        plt.figure(1)
        res = MeansAndSTDevs(EArray)
        line = [0 for v in range(4)]
        for index in range(4):
            line[index] = res[index][0]
        normaliser = np.copy(line[0])
        for index in range(4):
            line[index] -= normaliser
        plt.plot(line, '-', color='grey')
        plt.plot(line, 'o', color='black')
        plt.xticks(range(4), xarray)
        plt.xlabel('Added resistance (cmH2Os/L)', fontsize=30)
        plt.ylabel('Elastance (cmH2O/L)', fontsize=30)
        plt.legend(customLines, ['Data'])

        #Rocc
        resistLine = [0, 0.4, 0.8, 1.2]

        plt.figure(2)
        res = MeansAndSTDevs(RoccArray)
        line = [0 for v in range(4)]
        for index in range(4):
            line[index] = res[index][0]
        normaliser = np.copy(line[0])
        for index in range(4):
            line[index] -= normaliser
        plt.plot(line, '-', color='grey', label='Data')
        plt.plot(line, 'o', color='black')
        plt.plot(resistLine, '-', color='red', label='Expected trend')
        plt.xticks(range(4), xarray)
        plt.xlabel('Added resistance (cmH2Os/L)', fontsize=30)
        plt.ylabel('Rocc (cmH2Os/L)', fontsize=30)
        plt.legend(customLines, ['Data', 'Expected trend'])

        #Rrs
        plt.figure(3)
        res = MeansAndSTDevs(RArray)
        line = [0 for v in range(4)]
        for index in range(4):
            line[index] = -res[index][0]
        normaliser = np.copy(line[0])
        for index in range(4):
            line[index] -= normaliser
        plt.plot(line, '-', color='grey', label='Data')
        plt.plot(line, 'o', color='black')
        plt.plot(resistLine, '-', color='red', label='Expected trend')
        plt.xticks(range(4), xarray)
        plt.xlabel('Added resistance (cmH2Os/L)', fontsize=30)
        plt.ylabel('Rrs (cmH2Os/L)', fontsize=30)
        plt.legend(customLines, ['Data', 'Expected trend'])






    #==========
    # Tau

    tauPlots = 0
    if(tauPlots):

#        # tau vs 1/Rocc grad gives E
#        plt.figure()
#        plt.ylabel('Tau')
#        plt.xlabel('1/Rocc')
#        j = 0
#        for grads in gradients:
#            invGrads = [1/v for v in allRocc[j:j+len(grads)]]
#            plt.plot(invGrads, grads,  'o')
#            j += len(grads)
#        plt.grid()


        # Tau vs time (sort of, tau all plotted in order of occurrance)
        plt.figure()
        plt.ylabel('Decay rate')
        plt.xlabel('Shutter repetition')
        plt.plot(allGradients, '-o')
        plt.grid()

#        plt.figure()
#        plt.ylabel('peakPressure')
#        plt.xlabel('peakFlow')
#        plt.plot(allPeakShutterFlow, allStartShutterPressure, 'o')
#        plt.grid()

#        plt.figure()
#        plt.ylabel('peakPressure/peakFlow')
#        plt.xlabel('Rx')
#        for j in range(len(allRx)):
#            Rvalue = allStartShutterPressure[j]/allPeakShutterFlow[j]
#            plt.plot(allRx[j], Rvalue, 'o')
#        plt.grid()


        # boxplots
        print('Gradients, end gradients:')
        print(gradients)
        print(gradientsEnd)

        # boxplot tau vs test number (1-16)
        plt.figure()
        dataBoxy = []
        for j in range(len(gradients)):
            array = ([[gradients[j][v], j] for v in range(len(gradients[j]))])
            for item in array:
                dataBoxy.append(item)
        df = pd.DataFrame(dataBoxy, columns=['Tau', 'Test'])
        sns.boxplot(x='Test', y='Tau', data=df)#, 'Labels',())
        sns.swarmplot(x='Test', y='Tau', data=df, color="grey")
        plt.grid()


        # tau vs Rx
        plt.figure()
        dataBoxy = [[allGradients[v], allRx[v]] for v in range(len(allRx))]
        df = pd.DataFrame(dataBoxy, columns=['Decay rate', 'Rx'])
        sns.boxplot(x='Rx', y='Decay rate', data=df)#, 'Labels',())
        sns.swarmplot(x='Rx', y='Decay rate', data=df, color="grey")
        plt.grid()


        # end of shutter pressure vs Rx
        plt.figure()
        dataBoxy = [[allEndShutterPressure[v], allRx[v]] for v in range(len(allRx))]
        df = pd.DataFrame(dataBoxy, columns=['Pressure (Pa)', 'Rx'])
        sns.boxplot(x='Rx', y='Pressure (Pa)', data=df)#, 'Labels',())
        sns.swarmplot(x='Rx', y='Pressure (Pa)', data=df, color="grey")
        plt.grid()


        # tau/P0 vs Rx --> should match 1/volume closely
        plt.figure()
        dataBoxy = []
        for j in range(len(allGradients)):
            item = allGradients[j]/allEndShutterPressure[j]
            dataBoxy.append([item, allRx[j]])
        df = pd.DataFrame(dataBoxy, columns=['Tau/P0', 'Test'])
        sns.boxplot(x='Test', y='Tau/P0', data=df)#, 'Labels',())
        sns.swarmplot(x='Test', y='Tau/P0', data=df, color="grey")
        plt.grid()

        # 1/volume
        plt.figure()
        plt.ylabel('(Tau/P0), (1/volume)')
        plt.xlabel('Shutter repetition')
        divTau = [0]*len(allGradients)
        divVol = [0]*len(allGradients)
        for j in range(len(allGradients)):
            divTau[j] = allGradients[j]/(allEndShutterPressure[j])
            divVol[j] = 1/allShutterVolume[j]
        plt.plot(divTau, '-x')
        plt.plot(divVol, '-d')



        # E and R estimates
        Rd = [0]*len(allGradients)
        Ed = [0]*len(allGradients)
        for j in range(len(allGradients)):
            Rd[j] = 1/(allShutterVolume[j]*divTau[j])
            Ed[j] = allGradients[j]*(Rd[j])

        plt.rc('legend',**{'fontsize':12})
        c, (cax1, cax2) = plt.subplots(2, sharex=True)
#        cax1.plot(Rd, '-og')
#        cax2.plot(Ed, '-dr')
#
#        cax2.set_xlabel('Shutter repetition', fontsize=31)
#        cax1.set_ylabel('Estimates from Tau: (R, green)', fontsize=30)
#        cax2.set_ylabel('Estimates from Tau: (Ed, red)', fontsize=30)
        cax1.grid()
        cax2.grid()


        dataBoxy = [[Rd[v], allRx[v]] for v in range(len(allRx))]
        df = pd.DataFrame(dataBoxy, columns=['R', 'Rx'])
        sns.boxplot(x='Rx', y='R', data=df, ax=cax1)
        sns.swarmplot(x='Rx', y='R', data=df, color="grey", ax=cax1)

        dataBoxy = [[Ed[v], allRx[v]] for v in range(len(allRx))]
        df = pd.DataFrame(dataBoxy, columns=['E', 'Rx'])
        sns.boxplot(x='Rx', y='E', data=df, ax=cax2)
        sns.swarmplot(x='Rx', y='E', data=df, color="grey", ax=cax2)


        # end of shutter pressure vs Rx
        plt.figure()
        dataBoxy = [[allRocc[v], allRx[v]] for v in range(len(allRx))]
        df = pd.DataFrame(dataBoxy, columns=['Rocc', 'Rx'])
        sns.boxplot(x='Rx', y='Rocc', data=df)#, 'Labels',())
        sns.swarmplot(x='Rx', y='Rocc', data=df, color="grey")
        plt.grid()



        plt.show()

        #==========
        # Edrs
    EdrsPlots = 0
    if(EdrsPlots):



        # Edrs vs shutter test


        j = 0
        for grads in gradients:
            plt.figure()
            plt.ylabel('Edrs (o b=decay r=end), Rocc (x)')
            plt.xlabel('Shutter repetition')
            Edrs = [grads[v]*allRocc[j:j+len(grads)][v] for v in range(len(grads))]
            EdrsEnd = [allGradientsEnd[j:j+len(grads)][v]*allRocc[j:j+len(grads)][v] for v in range(len(grads))]
            Roccs = [allRocc[j:j+len(grads)][v] for v in range(len(grads))]
            plt.plot(Edrs, '-ob')
            plt.plot(EdrsEnd, '-or')
            plt.plot(Roccs, '-x')
            j += len(grads)

            plt.grid()





#=============================================================
# Final printouts for all datasets

# Table of Rocc vs Rx
# With mean and std
if(GenRoccTable):
    for line in range(len(tableIDs)):
        info = str(tableIDs[line])
        info += ', '
        for index in range(4):
            info += '{:.2f} [{:.2f}] (expect {:.2f} error {:.1f}%), '.format(RoccMeans[4*line+index], RoccStds[4*line+index], RoccExpected[4*line+index], RoccError[4*line+index])
        print(info)



# linear regression for relaxation gradients
if(RelaxationGradVsFlow):
    rvalues = [0 for v in range(6)]
    mvalues = [0 for v in range(6)]
    index = 0
    #S13_patch = mpatches.Patch(color='purple', marker='d', label='Subject 13')

    S13_patch = [Line2D([0], [0], color='purple', linestyle='', marker='d')]

    plt.figure(301)
    gradient, intercept, r_value, p_value, std_err = stats.linregress(megaVol, megaRelaxGrad)
    line = [0*gradient + intercept, 2.5*gradient + intercept]
    plt.plot([0, 2.5], line, color='red', linewidth=3)
    rvalues[index] = (r_value**2)
    mvalues[index] = (gradient)
    index += 1
    plt.legend(S13_patch, ["Subject 13"])
    plt.grid()


    plt.figure(303)
    gradient, intercept, r_value, p_value, std_err = stats.linregress(megaPstart, megaRelaxGrad)
    line = [0*gradient + intercept, 10*gradient + intercept]
    plt.plot([0, 10], line, color='red', linewidth=3)
    rvalues[index] = (r_value**2)
    mvalues[index] = (gradient)
    index += 1
    plt.legend(S13_patch, ["Subject 13"])
    plt.grid()


    plt.figure(304)
    gradient, intercept, r_value, p_value, std_err = stats.linregress(megaPpk, megaRelaxGrad)
    line = [0*gradient + intercept, 20*gradient + intercept]
    plt.plot([0, 20], line, color='red', linewidth=3)
    rvalues[index] = (r_value**2)
    mvalues[index] = (gradient)
    index += 1
    plt.legend(S13_patch, ["Subject 13"])
    plt.grid()


    plt.figure(307)

    mask = ~np.isnan(megaTau)
    megaTauMask = np.array(megaTau)[mask]
    megaRelaxGradMask = np.array(megaRelaxGrad)[mask]

    gradient, intercept, r_value, p_value, std_err = stats.linregress(megaTauMask, megaRelaxGradMask)
    line = [0*gradient + intercept, -60*gradient + intercept]
    plt.plot([0, -60], line, color='red', linewidth=3)
    rvalues[index] = (r_value**2)
    mvalues[index] = (gradient)
    index += 1
    plt.legend(S13_patch, ["Subject 13"])
    plt.grid()
    #plt.show()


    plt.figure(308)

    mask = ~np.isnan(megaE)
    megaEMask = np.array(megaE)[mask]
    megaRelaxGradMask = np.array(megaRelaxGrad)[mask]

    gradient, intercept, r_value, p_value, std_err = stats.linregress(megaEMask, megaRelaxGradMask)
    line = [0*gradient + intercept, 25*gradient + intercept]
    plt.plot([0, 25], line, color='red', linewidth=3)
    rvalues[index] = (r_value**2)
    mvalues[index] = (gradient)
    index += 1
    plt.legend(S13_patch, ["Subject 13"])
    plt.grid()


    plt.figure(302)
    gradient, intercept, r_value, p_value, std_err = stats.linregress(megaRocc, megaRelaxGrad)
    line = [0*gradient + intercept, 10*gradient + intercept]
    plt.plot([0, 10], line, color='red', linewidth=3)
    rvalues[index] = (r_value**2)
    mvalues[index] = (gradient)
    index += 1
    plt.legend(S13_patch, ["Subject 13"])
    plt.grid()

    print("overall regression R2value")
    for value in rvalues:
        print("& {:.2f}".format(value)),
    print(" \ \ \n")

    print("overall regression gradient")
    for value in mvalues:
        print("& {:.2f}".format(value)),
    print(" \ \ \n")


plt.show()

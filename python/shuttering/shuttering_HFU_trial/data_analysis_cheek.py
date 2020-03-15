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



def identifyOutliers(x, outlierConstant=1.5):
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
        self.flow = []


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

        self.expEndIndex = self.endFlowIndex + (len(self.flw)-self.endFlowIndex)*2/3
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
        fc = 120#filter cutoff
        bw = 30#filter bandwidth (taps?)
        filtFlow = hamming(flow, fc, 300, bw, plot=False)
        filtFlow = np.real(filtFlow).tolist()

        if(0):
            plt.figure()
            plt.plot(flow)
            plt.plot(filtFlow)
            plt.show()

        qminus15 = filtFlow[shutterClose-ms15]
        pminus15 = pressure[shutterClose-ms15]

        # Rocc as change in pressure due to shutter/flow pre shutter
        Rocc = (p15 - pminus15)/qminus15
       # print('Rocc: {}'.format(Rocc))

        if(0):
            plt.figure(123)
            plt.plot(pressure)
            plt.plot(shutterClose, pressure[shutterClose], 'o')
            plt.plot(shutterClose+ms30, pressure[shutterClose+ms30], 'x')
            plt.plot(shutterClose+ms75, pressure[shutterClose+ms75], 'x')
            plt.plot(shutterClose+ms15, p15, 'd')
            plt.plot(shutterClose-ms15, pressure[shutterClose-ms15], 'd')
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

def fourierPlot(data, samplePeriod, plotID):
    import scipy.fftpack

    # Number of samplepoints
    N = len(data)
    # sample spacing
    T = samplePeriod
    yf = scipy.fftpack.fft(data)
    xf = np.linspace(0.0, 1.0/(2.0*T), N/2)

    plt.figure(007)
    if('a' in plotID):
        plt.plot(xf, 2.0/N * np.abs(yf[:N//2]), 'r')
    else:
        plt.plot(xf, 2.0/N * np.abs(yf[:N//2]), 'b')
    #plt.show()

#=============================================================
#=============================================================

#for plotting
colours = ['blue','cyan', 'green','red','black','pink','purple','magenta', 'orange', 'brown', 'yellow','gray','olive','lavender','gold','deeppink','orchid','cadetblue','palegreen','blue','cyan', 'green','red','black','pink','purple','magenta', 'orange', 'brown', 'yellow','gray','olive','lavender','gold','deeppink','orchid','cadetblue','palegreen']
symbols = ['o','o','o','o','o','o','o','o','o','o','o','o','o','o','o','o','o','o','o','d','d','d','d','d','d','d','d','d','d','d','d','d','d','d','d','d','d','d',]
#colours = ['blue', 'green', 'red', 'black', 'purple', 'cyan', 'magenta', 'brown', 'orange', 'pink', 'yellow', 'blue', 'green', 'red', 'black', 'purple', 'cyan', 'magenta']
colourIndex = -1

megaID = []
megaEstart = []
megaRocc = []
megaEend = []
megaRaws = []
megaRelaxGrad= []
megaVolume = []
megaPMax = []

#==================
# Debug etc plots
STACK = 0
FULLSET = 0
STACKSHUTTER = 0
STACKPRESSURE = 0
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
dataNums = [15, 150, 16, 160, 17, 170, 18, 180, 112, 1120, 19, 190, 20, 200, 21, 210, 22, 220]
printNames = ['1', '1a', '2', '2a', '3', '3a', '4', '4a', '5', '5a', '6', '6a', '7', '7a', '8', '8a', '9', '9a']
#dataNums = [15, 150]
#printNames = ['1', '1a']
printnameindex = 0
for DATASET in dataNums:
    dataInfo = Datasets(DATASET)
    dataID = 'V'+ str(DATASET)
    printID = printNames[printnameindex]
    printnameindex += 1
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
    allEstart = []
    allRrs = []
    allShutterVolume= []
    allRelaxGrad = []
    allPeakAvgFlow = []
    allPeakShutterFlow = []
    allPeakAminusB = []
    allPeakPressure = []
    allStartShutterPressure = []
    allEndShutterPressure = []
    allStartShutterIndex = []
    allFourier = []


    #=============================================================
    measuredMechanics = mechanicsAnalyser()

    numFiles = min(numFiles, 2) # uncomment for only looking at no_Rx data
    for i in range(numFiles):
        colourIndex += 1

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
                allPeakPressure.append(fileBreaths.pres[peakPressureIndex])
                allPeakShutterFlow.append(peakShutterFlow)
                allPeakAminusB.append(peakAminusB)
                allStartShutterIndex.append(startIndex)
                allStartShutterPressure.append(fileBreaths.pres[startIndex]/100)
                allEndShutterPressure.append(endShutterPressure)
                allRocc.append(Rocc)
                allRx.append(Rx)
                allShutterVolume.append(shutterVolume)
                allE.append(Ed)
                allEstart.append(Estart)
                allRrs.append(Rrs)

                for amb in AminusBTrue:
                    allFourier.append(amb)


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
                    uax1.plot(AminusBVolTrue, AminusBTrue, linewidth = 3)
                    uax1.plot(AminusBVolTrue, AminusBTrue, ':', linewidth = 3)
                    uax1.plot(measuredQVAnalysis.volFilt, measuredQVAnalysis.flwFilt, '.-', linewidth = 3)
                    uax1.plot(AminusBVolTrue[measuredQVAnalysis.maxFlwIndex], AminusBTrue[measuredQVAnalysis.maxFlwIndex], 'o', color='yellow')
                    uax1.plot(AminusBVolTrue[measuredQVAnalysis.endFlowIndex], AminusBTrue[measuredQVAnalysis.endFlowIndex], 'o', color='yellow')
                    uax1.plot(AminusBVolTrue[measuredQVAnalysis.expEndIndex], AminusBTrue[measuredQVAnalysis.expEndIndex], 'o', color='yellow')
                    uax1.plot(AminusBVolTrue[measuredQVAnalysis.negFlowIndex], AminusBTrue[measuredQVAnalysis.negFlowIndex], 'o', color='cyan')

                    uax1.plot(averageVol, averageFlow, '-.', color='pink', linewidth = 3)
                    uax1.plot(averageVol[measuredQVAnalysis.maxFlwIndex], averageFlow[measuredQVAnalysis.maxFlwIndex], 'o', color='yellow')
                    uax1.plot(averageVol[measuredQVAnalysis.endFlowIndex], averageFlow[measuredQVAnalysis.endFlowIndex], 'o', color='yellow')
                    uax1.plot(averageVol[measuredQVAnalysis.expEndIndex], averageFlow[measuredQVAnalysis.expEndIndex], 'o', color='yellow')

                    if(not np.isnan(gradients[i][-1])):
                        uax1.plot([AminusBVolTrue[measuredQVAnalysis.expEndIndex], AminusBVolTrue[measuredQVAnalysis.maxFlwIndex]], lineFull, color='purple', linewidth=2)
                        uax1.plot([AminusBVolTrue[measuredQVAnalysis.expEndIndex], AminusBVolTrue[measuredQVAnalysis.maxFlwIndex]], lineFullEnd, color='purple', linewidth=2)

                    uax1.set_ylabel('Flow', fontsize=30)
                    uax1.set_xlabel('Volume', fontsize=30)
                    uax1.grid()


                    # flow and pressure plots
                    gax1.plot(averageFlow[:len(fileBreaths.flw)], color='purple', linewidth=1.5)
                    gax1.plot(fileBreaths.flw, color='red', linewidth=3)
                    gax1.plot(AminusBTrue, color='black', linewidth=3)

                    gax1.plot(measuredQVAnalysis.maxFlwIndex, fileBreaths.flw[measuredQVAnalysis.maxFlwIndex], 'o', color='yellow')
                    gax1.plot(measuredQVAnalysis.endFlowIndex, fileBreaths.flw[measuredQVAnalysis.endFlowIndex], 'o', color='yellow')
                    gax1.plot(measuredQVAnalysis.expEndIndex, fileBreaths.flw[measuredQVAnalysis.expEndIndex], 'o', color='yellow')

                    gax1.grid()
                    gax1.set_ylabel('Flow', fontsize=30)


                    gax2.plot(fileBreaths.pres, color='red', linewidth=3)
                    gax2.grid()
                    gax2.set_xlabel('DP', fontsize=31)
                    gax2.set_ylabel('Pressure', fontsize=30)


                    plt.show()
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


        # Tidy up :D
        averagePressures.clean()
        allAverageFlow = []



    if(1):
        fourierPlot(allFourier, 1/FS, printID)


    #=============================================================
    #=============================================================
    #=============================================================
    # Final analysis
    measuredMechanics.Rxs = dataInfo.xLabels



    #==========
    #==========
    # PLOTS PLOTS PLOTS
    #==========
    #==========

    for j in range(len(allRelaxGrad)):
        megaID.append(printID)
        megaEend.append(allE[j])
        megaRocc.append(allRocc[j])
        megaEstart.append(allEstart[j])
        megaRaws.append(allRrs[j])
        megaRelaxGrad.append(allRelaxGrad[j])

    outliers = identifyOutliers(allPeakPressure)
    for index in sorted(outliers, reverse=True):
        del allPeakPressure[index]
    # Mean of remining points
    megaPMax.append(np.mean(allPeakPressure))
    allPeakPressure = []

    outliers = identifyOutliers(allShutterVolume)
    for index in sorted(outliers, reverse=True):
        del allShutterVolume[index]
    # Mean of remining points
    megaVolume.append(np.mean(allShutterVolume))
    allShutterVOlume= []




#=============================================================
# Final printouts for all datasets


rvalues = [0 for thing in range(5)]
mvalues = [0 for thing in range(5)]
index = 0
sns.set(font_scale=2.5)

# Rocc
plt.figure(301)
dataBoxy = [[megaRocc[v], megaID[v]] for v in range(len(megaID))]
df = pd.DataFrame(dataBoxy, columns=['Rocc', 'ID'])
sns.boxplot(x='ID', y='Rocc', data=df)#, 'Labels',())
sns.swarmplot(x='ID', y='Rocc', data=df, color="grey")
plt.grid()


# E shutter close
plt.figure(302)
dataBoxy = [[megaEstart[v], megaID[v]] for v in range(len(megaID))]
df = pd.DataFrame(dataBoxy, columns=['E(shutter-close)', 'ID'])
sns.boxplot(x='ID', y='E(shutter-close)', data=df)#, 'Labels',())
sns.swarmplot(x='ID', y='E(shutter-close)', data=df, color="grey")
plt.grid()


# E shutter open
plt.figure(303)
dataBoxy = [[megaEend[v], megaID[v]] for v in range(len(megaID))]
df = pd.DataFrame(dataBoxy, columns=['Ed', 'ID'])
sns.boxplot(x='ID', y='Ed', data=df)#, 'Labels',())
sns.swarmplot(x='ID', y='Ed', data=df, color="grey")
plt.grid()


# Rrs from decay
plt.figure(304)
dataBoxy = [[megaRaws[v], megaID[v]] for v in range(len(megaID))]
df = pd.DataFrame(dataBoxy, columns=['Rrs', 'ID'])
sns.boxplot(x='ID', y='Rrs', data=df)#, 'Labels',())
sns.swarmplot(x='ID', y='Rrs', data=df, color="grey")
plt.grid()


# relaxation gradient
plt.figure(305)
dataBoxy = [[megaRelaxGrad[v], megaID[v]] for v in range(len(megaID))]
df = pd.DataFrame(dataBoxy, columns=['Relaxation gradient', 'ID'])
sns.boxplot(x='ID', y='Relaxation gradient', data=df)#, 'Labels',())
sns.swarmplot(x='ID', y='Relaxation gradient', data=df, color="grey")
plt.grid()


# change in elastance vs change in resistsnce plot
# To assess whether two parameters are linked
def extractMeanMechanic(ID, megaID, mixedArray):
    Hold = []
    for item in range(len(megaID)):
        if megaID[item] == ID:
            Hold.append(mixedArray[item])
    # Find outliers and remove
    outliers = identifyOutliers(Hold)
    for index in sorted(outliers, reverse=True):
        del Hold[index]
    # Mean of remining points
    Holdmean = np.mean(Hold)
    return Hold


increasesInE = [0 for hh in range(len(printNames)/2)]
increasesInR = [0 for hh in range(len(printNames)/2)]
increasesInM = [0 for hh in range(len(printNames)/2)]
indexer = 0
while indexer < len(printNames):
    # always nohold first
    name = printNames[indexer]
    # Mean of remining points
    noholdE= extractMeanMechanic(name, megaID, megaEstart)
    noholdR= extractMeanMechanic(name, megaID, megaRocc)
    noholdM= extractMeanMechanic(name, megaID, megaRelaxGrad)

    noholdEmean = np.mean(noholdE)
    noholdRmean = np.mean(noholdR)
    noholdMmean = np.mean(noholdM)


    # now with holding
    indexer += 1
    name = printNames[indexer]

    holdE = extractMeanMechanic(name, megaID, megaEstart)
    holdR = extractMeanMechanic(name, megaID, megaRocc)
    holdM = extractMeanMechanic(name, megaID, megaRelaxGrad)

    holdEmean = np.mean(holdE)
    holdRmean = np.mean(holdR)
    holdMmean = np.mean(holdM)

    # Compare amount increase from non-ch to ch
    EincreasePercent = ((holdEmean-noholdEmean)/noholdEmean)*100
    RincreasePercent = ((holdRmean-noholdRmean)/noholdRmean)*100
    MincreasePercent = ((holdMmean-noholdMmean)/noholdMmean)*100


    #store results
    increasesInE[(indexer-1)/2] = EincreasePercent
    increasesInR[(indexer-1)/2] = RincreasePercent
    increasesInM[(indexer-1)/2] = MincreasePercent
    indexer += 1

# print and plot
print(increasesInE)
print(increasesInR)
print(increasesInM)
print(megaPMax)
print(megaVolume)

d = {'Elastance':increasesInE, 'Resistance':increasesInR, 'Viscoelastic':increasesInM}
df = pd.DataFrame(data=d)
corr = df.corr()
print(corr)

if (len(dataNums)>=6):
    g = sns.PairGrid(df)
    g.map_upper(plt.scatter)
    g.map_lower(sns.kdeplot)
    g.map_diag(sns.kdeplot, lw=3, legend=False)

plt.show()


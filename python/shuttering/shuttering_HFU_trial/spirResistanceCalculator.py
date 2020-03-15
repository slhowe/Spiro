
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
from dataStore import Datasets

#=============================================================
#=============================================================

#for plotting
colours = ['blue','cyan', 'green','green', 'red', 'red','black','black','blue','blue', 'green','green', 'red', 'red','black','black','pink','pink',]
#colours = ['blue', 'green', 'red', 'black', 'purple', 'cyan', 'magenta', 'brown', 'orange', 'pink', 'yellow', 'blue', 'green', 'red', 'black', 'purple', 'cyan', 'magenta']
colourIndex = -1


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
                print("This breath is shuttered")
                if(shutterIndex > (breathMiddle - (breathMiddle-breathStart)/4)):
                    print("But shuttering started later than last quarter of expiration")
                    return False
                else:
                    self.shutterIndex = shutterIndex-breathStart # Save shutter start relative to breath start
                    return True
            # if shutter started beyond end-exp of this breath
            elif(shutterIndex > breathEnd):
                return False
        return False




#=============================================================
dataNums = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11]
dataNums = [9]
for DATASET in dataNums:
    dataInfo = Datasets(DATASET)
    numFiles= len(dataInfo.files)


    #=============================================================

    for i in range(1):

        #####################
        # GET DATA FROM FILE#
        #####################

        filename = dataInfo.path + dataInfo.files[i]
        print("/n" + filename)

        fileData = Data(filename)
        fileData.extractData()

        fileBreaths = Breaths(fileData)
        fileBreaths.splitBreaths()


        CALCSPIROMETERRESISTANCE = 1
        if(CALCSPIROMETERRESISTANCE):
            fileBreaths.extractSingleBreath(0)
            if(fileBreaths.breathIsShuttered(0)):
                print('Ah, breath is shuttered.. unable to calc Rspr')
            else:
                Rspr = [0 for v in range(len(fileData.flow))]
                for index in range(len(Rspr)):
                    if(fileData.flow[index] == 0):
                        Rspr[index] = np.nan
                    else:
                        Rspr[index] = fileData.pressure[index]/fileData.flow[index]

                plt.figure()
                plt.plot(Rspr, 'ob')
                plt.grid()
                plt.show()

                plt.figure()
                plt.plot(fileData.pressure, 'b')
                plt.plot(fileData.flow, 'r')

                fc = 30#filter cutoff
                bw = 20#filter bandwidth (taps?)
                pres = hamming(fileData.pressure, fc, 300, bw, plot=False)
                pres = np.real(pres).tolist()
                plt.plot(pres, 'c')
                plt.grid()

                Rsprf = [0 for v in range(len(fileData.flow))]
                for index in range(len(Rspr)):
                    if(fileData.flow[index] == 0):
                        Rsprf[index] = np.nan
                    else:
                        Rsprf[index] = pres[index]/fileData.flow[index]

                plt.figure()
                plt.plot(Rsprf, 'or')
                plt.plot(Rspr, 'xk')
                plt.grid()
                plt.show()



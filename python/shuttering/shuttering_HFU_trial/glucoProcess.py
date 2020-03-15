#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jun 21 13:59:41 2019

@author: lui holder-pearson

    This file is part of Pulse Glucometery.

    Pulse Glucometery is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    Pulse Glucometery is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with Pulse Glucometery.  If not, see <https://www.gnu.org/licenses/>.
    

"""

#from heart_rate.py import HeartRate

import numpy as np
import matplotlib.pyplot as plt
from scipy import signal
import scipy as sp


#New main function will look kinda like this:
    
def main():
    #plt.close('all')
    
    file = 'test1/dataContPro.csv'
    fs = 66 # Hz
    plotlvl = 2
    '''
    Plot level:
        0 -> nothing
        1 -> meaningful results used for analysis
        2 -> intermediate data results (plus above
        3 -> details of individual steps (eg filter responses) plus above
        4 -> intricate details of individual steps (eg filter taps) plus above 
    '''
    
    data1 = glucoData(file, fs, plotlvl)
    
    data1.applyFilters()
    
    ignore1 = [1]
    
    ints = np.asarray([[502.1, 506.5],
                       [380.15, 383.71],
                       [821.097, 823.75],
                       [584.24, 587.921],
                       [916.1, 919.45],
                       [996.38, 1000.07],
                       [1039.15, 1042],
                       [1184.65, 1189.06],
                       [1593.2, 1597.8],
                       [1601.67, 1605.25],
                       [1676.11, 1686.2],
                       [1753.71, 1757.72],
                       [1817.35, 1820],
                       [1821.44, 1828.45],
                       [1855.06, 1857.94],
                       [1859.59, 1862.5],
                       [1933.94, 1937.27],
                       [1953, 1956.56],
                       [1960.55, 1963.07],
                       [1999.84, 2002.3]
                       ])
    
    numSects = np.shape(ints)[0]
    
    
    defaultLambda = ints[:,0] * 0
    '''
    This sets the default wavelength which we get the number of zerocrossings from
    '''
    defaultLambda[13] = -1
    defaultLambda[16] = 3
    
    defaultLambda = defaultLambda.astype(int)
        
    
    pkList = [peakData(data1, ints[i,:], ignore1) for i in range(numSects)]
    
    for i in range(numSects - 1):
        print("Computing peak set {}.".format(i))
        pkList[i].compute(defaultLambda[i])
        
    try:
        print("Computing peak set {}.".format(numSects - 1))
        pkList[-1].compute(defaultLambda[-1])
    except BaseException as e:
        print("Failed:")
        pkList[-2].plotAbs(ints[-1,:])
        raise Exception(e)
    else:
        pkList[-1].plotAbs()
        #pkList[-1].centPlot()

    
    

class glucoData ():
    
    def __init__ (self, filepath, fs, plotlevel = 0):
        
        # Settings
        self.plotting = plotlevel
        
        # Initial data read and subsequent things
        self.dataRaw = np.genfromtxt(filepath, delimiter=',')
        self.numLambda = np.shape(self.dataRaw)[1]
        self.fs = fs
        self.nyq = self.fs/2
        self.x = np.linspace(0, np.shape(self.dataRaw)[0]/self.fs, np.shape(self.dataRaw)[0])
        
        
        
    def applyFilters(self, ripple = 50, width = 1, cutoffLP = 5, order=7, cutoffBaseline = .8):
        self.lowPass(ripple, width, cutoffLP)
        self.highPass(order, cutoffBaseline)
        
        dataFiltSub = self.dataFilt[self.delay_discrete:,:]
        dataBaseLineSub = self.dataBaseLine[:-self.delay_discrete,:]
        
        self.dataOut = np.subtract(dataFiltSub, dataBaseLineSub)
        
        if self.plotting >= 2:
            plt.figure()
            plt.plot(self.x , self.dataRaw)
            plt.plot(self.x-self.delay, self.dataFilt, linewidth=3)
            plt.plot(self.x, self.dataBaseLine, linewidth=2)
            plt.legend(3*['660', '850', '940', '1450', '1550', '1650'])
            #plt.plot(x-self.delay, dataFilt)
            #plt.legend(['660', '940', '1450', '1550', '1650','850'])
            
            plt.grid(True)
        
        
    def lowPass(self, ripple, width, cutoff):
        self.numtaps, beta = signal.kaiserord(ripple, width/self.nyq)
        self.taps = signal.firwin(self.numtaps, cutoff, window = ('kaiser', beta), scale=False, nyq = .5*self.fs)
        
        self.dataFilt = signal.lfilter(self.taps, 1.0, self.dataRaw, axis=0)
        
        self.delay_discrete = int(0.5 *(self.numtaps-1))
        self.delay = self.delay_discrete / self.fs

        if self.plotting >= 3:
            plt.figure()
            w, h = signal.freqz(self.taps, worN=8000)
            plt.plot((w/np.pi)*self.fs*.5, np.abs(h), linewidth=2)
            plt.xlabel('Frequency (Hz)')
            plt.ylabel('Gain')
            plt.title('Frequency Response')
            plt.ylim(-0.05, 1.05)
            plt.grid(True)
        
            # Upper inset plot.
            plt.axes([0.42, 0.6, .45, .25])
            plt.plot((w/np.pi)*self.fs*.5, np.abs(h), linewidth=2)
            plt.xlim(0,8.0)
            plt.ylim(0.95, 1.05)
            plt.grid(True)
            
            # Lower inset plot
            plt.axes([0.42, 0.25, .45, .25])
            plt.plot((w/np.pi)*self.fs*.5, np.abs(h), linewidth=2)
            plt.xlim(6, 12)
            plt.ylim(0.0, 0.0105)
            plt.grid(True)
            
        if self.plotting >= 4:
            plt.figure()
            plt.plot(self.taps, 'bo-', linewidth=1)
            plt.title('Filter Coefficients (%d taps)' % self.numtaps)
            plt.grid(True)
    
     
    def highPass(self, order, cutoff):
        
        self.b, self.a = signal.butter(order, cutoff/self.nyq, 'lowpass')
        self.dataBaseLine = signal.filtfilt(self.b,self.a, self.dataRaw, axis = 0)
        
        if self.plotting >= 3:
            
            plt.figure()
            
            w, h = signal.freqz(self.b, self.a, worN = 8000)
            plt.plot(w*self.nyq, 20*np.log10(np.abs(h)))
            #plt.xscale('log')
            plt.title('Baseline filter characteristics')
            plt.xlabel('Frequency [radians / second] NOT ENTIRELY SURE OF UNITS!!')
            plt.ylabel('Amplitude [dB]')
            plt.margins(0, 0.1)
            plt.grid(which='both', axis='both')
            plt.axvline(cutoff, color='green') # cutoff frequency
            plt.show()
    
    
    def plotFilt(self):
        self.figAbs, self.axsAbs = plt.subplots(self.numLambda,1, sharex = True, num = 102)
        
        for row in range(self.numLambda):
            self.axsAbs[row].plot(self.x[:-self.delay_discrete], self.dataOut[:,row])
            
    
    
    def setPlotlevel(self, level):
        self.plotting = level    
        
class peakData():
    
    def __init__(self, gdata, tRange, ignore=np.empty([0,0])):
        
        self.tRange = tRange
        self.usedIndices = [x for x in range(gdata.numLambda) if x not in ignore]
        self.numLambda = np.size(self.usedIndices)
        
        self.int_dis = [int(np.floor(self.tRange[0]*gdata.fs) - 1), int(np.ceil(self.tRange[1]*gdata.fs) + 1)]
        self.gdata = gdata
        
    def compute(self, dLambda = 0, minTimeDisp = 0.4):
        self.mySection = self.gdata.dataOut[self.int_dis[0]:self.int_dis[1], self.usedIndices]
        
        self.getCrossings(dLambda, minTimeDisp)
        self.getPeakH()
        self.getPeakA()
        
        
        
    def getCrossings(self, dLambda, minTimeDisparity):
        self.numCrossings = np.zeros([np.shape(self.mySection)[1]])
        
        zeroCrossSingle = np.diff(np.sign(self.mySection[:,dLambda])).nonzero()[0]
        self.numCross = np.size(zeroCrossSingle)
        zeroCrossings_raw = np.empty([self.numCross, 0])
        
        for row in range(self.numLambda):
            zeroCrossSingle = np.diff(np.sign(self.mySection[:,row])).nonzero()[0]
            #if np.size(zeroCrossSingle) > self.numCross:
                # Then we have an issue!!
                #getRidSmallest(zeroCrossSingle)
                
            zeroCrossSingle = zeroCrossSingle[:self.numCross].reshape(self.numCross, 1)
            '''
            print(row)
            print(zeroCrossSingle)
            print(zeroCrossings_raw)
            print('')
            '''
            zeroCrossings_raw = np.append(zeroCrossings_raw, zeroCrossSingle, axis = 1)
            self.numCrossings[row] = np.size(zeroCrossSingle)
            
        assert np.all(self.numCrossings == self.numCrossings[0]), "Number of zero crossings isn't constant"
        self.numCross = int(self.numCrossings[0])
        
        assert np.all(self.mySection[0,:] < 0), "Not all staring before a trough"
        
        
        #if np.any(np.abs(np.diff(np.sum(zeroCrossings_raw, axis = 0)))>minTimeDisparity*self.gdata.fs):
        while np.any(np.abs(np.diff(zeroCrossings_raw[-1,:]))>minTimeDisparity*self.gdata.fs):
            # This means that there is quite a big time disparity, most likely due to a dichrotic notch
            zeroCrossings_raw, aRow, zC= self.getCrossings_removeDichrotic(zeroCrossings_raw)
            
        
        self.zeroCrossings = np.append(np.zeros([1,self.numLambda]), zeroCrossings_raw, axis = 0)
        
        
    def getCrossings_removeDichrotic(self, zC_raw):
        aRow = np.argmin(np.sum(zC_raw, axis = 0))
        
        
        zeroCrossSingle = np.diff(np.sign(self.mySection[:,aRow])).nonzero()[0]
        numToRem = int(np.floor((np.size(zeroCrossSingle) - np.shape(zC_raw)[0]) / 2))
        
        for i in range(numToRem):
            aCr = np.argmin(np.diff(zeroCrossSingle, axis = 0))
            
            ind = [w for w in range(np.size(zeroCrossSingle)) if w not in [aCr, aCr+1]]
            zeroCrossSingle = zeroCrossSingle[ind]
            print("Dichrotic notch detected and fixed on (used) LED {0}, ZC {1}".format(aRow, aCr))
        #zeroCrossSingle = zeroCrossSingle[:self.numCross].reshape(self.numCross, 1)
        
        
            
        
        zC_raw[:,aRow] = zeroCrossSingle
        
        return zC_raw, aRow, aCr
        
    def getPeakH(self):
        self.absolutes = np.zeros([self.numCross, self.numLambda])
        self.absTimes = self.absolutes + 0
        
        
        
        for row in range(self.numLambda):
            curSign = -1
            skipped = 0
            for ind in range(self.numCross):
                ind = ind - skipped
                self.absolutes[ind,row] = curSign * np.amax(np.abs(self.mySection[int(self.zeroCrossings[ind, row]):int(self.zeroCrossings[ind+1, row]),row]))
                
                self.absTimes[ind,row] = np.argmax(np.abs(self.mySection[int(self.zeroCrossings[ind, row]):int(self.zeroCrossings[ind+1, row]),row]))
                
                
                curSign = curSign * -1
        
        self.absTimes_abs_disc = self.absTimes + self.zeroCrossings[:-1,:] + self.int_dis[0]
        self.absTimes_abs_disc = self.absTimes_abs_disc.astype(int)
        
        self.absTimes_abs_disc_plt = np.empty([self.numCross, self.gdata.numLambda])
        self.absTimes_abs_disc_plt[:] = np.nan
        
        self.absolutes_plt = self.absTimes_abs_disc_plt + 0
        
        self.absTimes_abs_disc_plt[:,self.usedIndices] = self.absTimes_abs_disc
        self.absTimes_abs_disc_plt = self.absTimes_abs_disc_plt.astype(int)
        self.absolutes_plt[:,self.usedIndices] = self.absolutes
        
        self.heights = np.diff(self.absolutes, axis = 0)[0::2]
        
        if not hasattr(self, 'numPeaks'):
            self.numPeaks = int(np.floor(self.numCross/2))
            
        self.heights_plt = np.ones([self.numPeaks, self.gdata.numLambda])
        self.heights_plt[:] = np.nan
        
        self.heights_plt[:,self.usedIndices] = self.heights
        self.hAvg = np.mean(self.heights_plt, axis = 0)
        
        
    def getPeakA(self):
        
        if not hasattr(self, 'numPeaks'):
            self.numPeaks = int(np.floor(self.numCross/2))
            
        self.areas = np.ones([self.numPeaks, self.gdata.numLambda])
        self.areas[:] = np.nan
        
        for row in self.usedIndices:
            for peak in range(self.numPeaks):
                init = self.gdata.dataOut[self.absTimes_abs_disc_plt[peak*2, row], row]
                fin = self.gdata.dataOut[self.absTimes_abs_disc_plt[(peak+1)*2, row], row]
                baseL = np.linspace(init,fin, self.absTimes_abs_disc_plt[(peak+1)*2, row] - self.absTimes_abs_disc_plt[(peak)*2, row])
                y_local = self.gdata.dataOut[self.absTimes_abs_disc_plt[peak*2, row]:self.absTimes_abs_disc_plt[(peak+1)*2, row], row] - baseL
                x_local = self.gdata.x[self.absTimes_abs_disc_plt[peak*2, row]:self.absTimes_abs_disc_plt[(peak+1)*2, row]]
                
                self.areas[peak, row] = sp.integrate.cumtrapz(y_local,x_local)[-1]
        
        self.aAvg = np.mean(self.areas, axis = 0)
        
    def plotAbs(self,leg = ['660', '850', '940', '1450', '1550', '1650'],  extra = None):
        
        assert hasattr(self, 'absolutes_plt'), "No peaks to plot!"
        self.figAbs, self.axsAbs = plt.subplots(self.gdata.numLambda,1, sharex = True, num = 101)
        
        for row in range(self.gdata.numLambda):
            self.axsAbs[row].plot(self.gdata.x[:-self.gdata.delay_discrete], self.gdata.dataOut[:,row])
            self.axsAbs[row].plot(self.absTimes_abs_disc_plt[:,row]/self.gdata.fs, self.absolutes_plt[:,row], 'r+')
            self.axsAbs[row].axvline(x=self.tRange[0],color='k', linestyle='--')
            self.axsAbs[row].axvline(x=self.tRange[1],color='k', linestyle='--')
            
            if np.any(extra):
                self.axsAbs[row].axvline(x=extra[0],color='k', linestyle='--')
                self.axsAbs[row].axvline(x=extra[1],color='k', linestyle='--')        
            
            self.axsAbs[row].grid(True)
            self.axsAbs[row].set_ylabel(leg[row] + "nm")
            
    def centPlot (self):
        #print(self.figAbs)
        #plt.figure(self.figAbs)
        
        newx = np.add(self.tRange, [-2, 2])
        newymin = np.amin(self.absolutes_plt, axis = 0) * 1.4
        newymax = np.amax(self.absolutes_plt, axis = 0) * 1.4
        
        for row in self.usedIndices:
            self.axsAbs[row].set_xlim(newx[0], newx[1])
            self.axsAbs[row].set_ylim(newymin[row], newymax[row])
        

if __name__ == "__main__":
    main()
    
    
    
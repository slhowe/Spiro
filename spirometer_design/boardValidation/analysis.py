#!/bin/bash

# Import my extensions
import sys
sys.path.insert(0, '/home/sarah/Documents/Spirometry/python/extensions')
from scipy.fftpack import fft
from scipy.signal import butter
import scipy.signal as signal
import scipy.stats as stats
import csv
import matplotlib.pyplot as plt
import numpy as np
from calculus import integral

path = '/home/sarah/Documents/Spirometry/spirometer_design/boardValidation/'
files = [
        'NewFile1.csv', # Shutter
        'NewFile2.csv'  # NO shutter
         ]



# Create data classes
class Data:
    def __init__(self):
        self.solIn = []
        self.solOut = []

    def addSolIn(self, dp):
        self.solIn.append(dp)

    def addSolOut(self, dp):
        self.solOut.append(dp)



dataClasses = []

for i in range(len(files)):
    data = Data()

    # Store pressure and flow data
    filename = path + files[i]
    with open(filename, 'rb') as csvfile:
        reader = csv.reader(csvfile, delimiter = ',')

        # skip 2 header line
        header = reader.next()
        header = reader.next()

        for row in reader:
            CH1 = float(row[1])
            CH2 = float(row[2])

            data.addSolIn(CH1)
            data.addSolOut(CH2)
    dataClasses.append(data)


# finding shutter draw back time.
# manually found distance between peaks
y = dataClasses[0].solIn
N = len(y)
T = 1.0/5000.0

fs = 5000.0
fc = 600
nyq = 0.5 * fs
fc = fc/nyq
order = 3

b, a = butter(order, fc, 'low')
sigWithShutter =  signal.filtfilt(b, a, dataClasses[0].solIn)
sigNoShutter =  signal.filtfilt(b, a, dataClasses[1].solIn)
sigDif = [sigWithShutter[o]-sigNoShutter[o] for o in range(800)]

# find peaks
sigPeak = []
for i in range(len(sigDif)):
    if sigDif[i] >= 0.5:
        sigPeak.append(sigDif[i])

# pop values less than max within 10ms
Peaks = []
minWidth = int(fs*0.01)
peakIndex = 0
for value in range(len(sigPeak)-1): # are in time order
    #check next index is within peakWidth
    if(sigDif.index(sigPeak[peakIndex+1]) > (sigDif.index(sigPeak[peakIndex])+minWidth)):
        peakIndex += 1
    # if it is, pop lower of two values
    else:
        localPeak = sigPeak[peakIndex]
        if(localPeak < sigPeak[peakIndex]+1):
            sigPeak.pop(peakIndex)
        else:
            sigPeak.pop(peakIndex+1)

peaks = [sigDif.index(o) for o in sigPeak]
print('Time to close: {}ms'.format((peaks[1]-peaks[0])*1000/fs))

# Spring constant
F = [-10, -5, -2, 0, 2, 10]     # mass added
x = [27, 22, 19, 18, 15.5, 8]   # distance moved

F = [o*9.81/1000.0 for o in F]
x = [o/1000.0 for o in x]

grad, offs, r, p, e, = stats.linregress(x, F)
line = [offs + grad*i for i in x]
print('Spring const, k, :{}'.format(grad))



# simulate time for shutter re-opening
x0 = 13.5/1000.0 #m
xf = 10.5/1000.0 #m
M_plug = 9.6/1000.0 #kg
k = -grad

# x(t) = x0.cos(sqrt(k/m)t)
dt = 0.0001 #s = 0.1ms

simX = []
t = 0
xt = x0
while xt >= 0:
    xt = x0*np.cos(np.sqrt(k/M_plug)*t)
    simX.append(xt)
    t+=dt

finalPt = 0
for i in range(len(simX)):
    if simX[i] >= xf:
        finalPt = i
        i = len(simX) + 1
tf = finalPt*dt
print('Time to re-open: {}ms'.format(tf*1000))



plt.figure()
plt.plot(sigWithShutter, 'r', linewidth=2,)
plt.plot(sigNoShutter, 'c', linewidth=2,)
plt.plot(sigDif, color='purple', linewidth=2)

plt.grid()
plt.ylabel('Voltage (V)', fontsize=30)
plt.xlabel('Data point', fontsize=30)
plt.legend([
            "Voltage with shutter",
            "Voltage with NO shutter",
           ])


plt.figure()
plt.plot(x, F, 'o')
plt.plot(x, line, 'r', linewidth=2)

plt.grid()
plt.legend([
    'Measurements',
    'Regression'
    ])
plt.xlabel('Spring length (mm)', fontsize=30)
plt.ylabel('Mass (g)', fontsize=30)



plt.figure()
t = np.linspace(0, dt*len(simX), len(simX))
plt.plot(t, simX, linewidth=2)
plt.plot(tf, simX[finalPt], 'ro')

plt.grid()
plt.xlabel('Time (s)', fontsize=30)
plt.ylabel('Extension (m)', fontsize=30)
plt.legend(['Simulated spring position',
            'Open position'
            ])

plt.show()



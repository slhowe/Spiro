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

#
# Theory is:
# Combined data = sine_data + pulse_data
# This is for flow
#

path = './'
files = [
        'comb2.csv',
        'sine2.csv',
         ]


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

flow = [[],[]]
time = [[],[]]
volume = [[],[]]
pressure = []

for i in range(len(files)):

    #####################
    # GET DATA FROM FILE#
    #####################
    # Store pressure and flow data
    filename = path + files[i]
    print(' ')
    print(filename)
    if(i!=3):
        with open(filename, 'rb') as csvfile:
            reader = csv.reader(csvfile, delimiter = ';')

            # skip header line
            header = reader.next()

            for row in reader:
                dp_flow = float(row[1])
                dp_time = float(row[0]) # in ms

                flow[i].append(dp_flow)
                time[i].append(dp_time)
        fs = 1/(time[0][1] - time[0][0]) #filter cutoff
        volume[i] = integral(flow[i], fs)
    else:
        with open(filename, 'rb') as csvfile:
            reader = csv.reader(csvfile, delimiter = ';')

            # skip header line
            header = reader.next()

            for row in reader:
                dp_pres = float(row[1])
                dp_time = float(row[0]) # in ms

                pressure.append(dp_pres)

#######################
# Split breaths apart #
#######################
flow_splits = quick_breath_split(flow[0], 2, timeout=400, RISING_EDGE=True)
flow_splits = [f+2 for f in flow_splits]
flow_splits.append(len(flow[0])-1)


## Doubling flow to see what happens for higher efforts ##
flow[0] = [5.3*q for q in flow[0]]

##############################
# Calc sum of sine and pulse #
##############################
sumFlow = [flow[0][m] - flow[1][m] for m in range(len(flow[0]))]
sumVolume = integral(sumFlow, fs)

###############
# Decay rates #
###############
LEN = 20
flwPulse = sumFlow[flow_splits[0]:flow_splits[0]+LEN]
flwComb= flow[0][flow_splits[0]:flow_splits[0]+LEN]

volPulse = sumVolume[flow_splits[0]:flow_splits[0]+LEN]
volComb= volume[0][flow_splits[0]:flow_splits[0]+LEN]

t = time[0][flow_splits[0]:flow_splits[0]+LEN]

gradientPassivePulse, intercept, r_value, p_value, std_err = stats.linregress(volPulse, flwPulse)
print(r_value)
line = [intercept + gradientPassivePulse*m for m in [sumVolume[flow_splits[0]],sumVolume[flow_splits[0]+LEN]]]
gradientPassiveComb, intercept, r_value, p_value, std_err = stats.linregress(volComb, flwComb)
print('     Mechanics from linear')
print('Pulse: {}'.format(gradientPassivePulse))
print('Comb: {}'.format(gradientPassiveComb))


# Calc exp decay of flow
ln_flw = [np.log(abs(f)) for f in flwPulse]
times = t
ones = [1]*len(times)
dependent = np.array([ln_flw])
independent = np.array([ones, times])
try:
    res = lstsq(independent.T, dependent.T)
    decayPulse = (res[0][1][0])
    offset = np.exp(res[0][0][0])
except(ValueError):
    print('ValueError: Data has nan?')
    decayPulse = nan
    offset = nan

print('     Decay rate stuff')
print('Pulse: {}'.format(decayPulse, offset,))

flwPulseCurve = [offset*exp(decayPulse*y) for y in times]


# Calc exp decay of flow
ln_flw = [np.log(abs(f)) for f in volPulse]
times = t
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

print('PulseVol: {}'.format(decay, offset,))


ln_flw = [np.log(abs(f)) for f in flwComb]
times = t
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

print('Comb: {}'.format(decay, offset,))

flwCombCurve = [offset*exp(decay*y) for y in times]

###########################
# effective added voltage #
###########################
R = 2
C = 0.08
print('expected tau: {}'.format(1/(R*C)))
Veff = [sumFlow[i]*R + sumVolume[i]/C for i in range(len(sumFlow))]
Qeff = [(Veff[i]/R)*exp(-time[0][i]/R/C) for i in range(len(sumFlow))]
decay= [exp(-time[0][i]/R/C) for i in range(len(sumFlow))]

plt.plot(volume[0], 'b')
plt.plot(volume[1], 'g')
plt.plot(sumVolume, 'r')
plt.plot([volume[0][i] - volume[1][i] for i in range(len(volume[0]))], 'orange')
plt.title('Checking that vol(comb-sine) = vol(comb)-vol(sine)')
plt.show()

#####################
# Plots plots plots #
#####################
plt.rc('legend',**{'fontsize':12})
f, (ax1, ax2) = plt.subplots(2, sharex=False)

#ax1.plot(pressure)
ax1.plot(sumVolume, sumFlow)
ax1.plot([sumVolume[flow_splits[0]],sumVolume[flow_splits[0]+LEN]], line, 'g')
ax1.plot(sumVolume[flow_splits[0]], sumFlow[flow_splits[0]], 'oy')
ax1.plot(sumVolume[flow_splits[0]+LEN], sumFlow[flow_splits[0]+LEN], 'om')
#ax1.plot(flow_splits, [pressure[m] for m in flow_splits], 'or')

for i in range(len(files)):
    ax2.plot(flow[i])
ax2.plot(range(flow_splits[0], flow_splits[0]+LEN), flwPulseCurve, 'm', linewidth=2)
ax2.plot(sumFlow, ':', color='black')
#ax2.plot(range(flow_splits[0], flow_splits[0]+LEN), flwComb, color='purple', linewidth=2)
ax2.plot(flow_splits, [flow[i][m] for m in flow_splits], 'or')
ax1.legend(['Shuttered data', 'linear fit (m={})'.format(gradientPassivePulse)])
ax1.set_ylabel('Flow L/s', fontsize=30)
ax1.set_xlabel('Volume L', fontsize=30)
ax2.set_ylabel('Flow L/s', fontsize=30)
ax2.set_xlabel('Time (shutter=150ms)', fontsize=30)
ax2.legend(['Measured flow',
            'Flow due to body',
            'Decay fit (tau={})'.format(decayPulse),
            'Flow due to shutter',
            ])
plt.show()



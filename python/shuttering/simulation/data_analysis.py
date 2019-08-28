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
        'comb.csv',
        'sine.csv',
        'pulse.csv',
        'pressure.csv',
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

flow = [[],[],[]]
time = [[],[],[]]
volume = [[],[],[]]
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


##############################
# Calc sum of sine and pulse #
##############################
sumFlow = [flow[1][m] + flow[2][m] for m in range(len(flow[1]))]

difference = [sumFlow[m] - flow[0][m] for m in range(len(flow[0]))]
print('     Difference between comb and sine+pulse:')
print('mean: {}'.format(np.mean(difference)))
print('min: {}'.format(min(difference)))
print('max: {}'.format(max(difference)))


###############
# Decay rates #
###############
LEN = 100
flwPulse = flow[2][flow_splits[0]:flow_splits[0]+LEN]
flwComb= flow[0][flow_splits[0]:flow_splits[0]+LEN]

volPulse = volume[2][flow_splits[0]:flow_splits[0]+LEN]
volComb= volume[0][flow_splits[0]:flow_splits[0]+LEN]

t = time[0][flow_splits[0]:flow_splits[0]+LEN]

gradientPassivePulse, intercept, r_value, p_value, std_err = stats.linregress(volPulse, flwPulse)
gradientPassiveComb, intercept, r_value, p_value, std_err = stats.linregress(volComb, flwComb)
print('     Mechanics from linear')
print('Pulse: {}'.format(gradientPassivePulse))
print('Comb: {}'.format(gradientPassiveComb))
#line = [intercept + gradientPassivePulse*m for m in [vol[maxFlwIndex], vol[new_end]]]


# Calc exp decay of flow
ln_flw = [np.log(abs(f)) for f in flwPulse]
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

print('     Decay rate stuff')
print('Pulse: {}'.format(decay, offset,))

flwPulseCurve = [offset*exp(decay*y) for y in times]

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


plt.rc('legend',**{'fontsize':12})
f, (ax1, ax2) = plt.subplots(2, sharex=False)

#ax1.plot(pressure)
ax1.plot(volume[2], flow[2])
ax1.plot(volume[2][flow_splits[0]], flow[2][flow_splits[0]], 'oy')
ax1.plot(volume[2][flow_splits[0]+LEN], flow[2][flow_splits[0]+LEN], 'om')
#ax1.plot(flow_splits, [pressure[m] for m in flow_splits], 'or')
for i in range(len(files)-1):
    ax2.plot(flow[i])
ax2.plot(range(flow_splits[0], flow_splits[0]+LEN), flwPulseCurve, 'm', linewidth=2)
#ax2.plot(range(flow_splits[0], flow_splits[0]+LEN), flwComb, color='purple', linewidth=2)
ax2.plot(sumFlow, ':', color='black')
ax2.plot(flow_splits, [flow[i][m] for m in flow_splits], 'or')
ax1.legend(['Shuttered data'])
ax1.set_ylabel('Flow L/s', fontsize=30)
ax1.set_xlabel('Volume L', fontsize=30)
ax2.set_ylabel('Flow L/s', fontsize=30)
ax2.set_xlabel('Time (shutter=150ms)', fontsize=30)
ax2.legend(['Measured flow',
            'Flow due to body',
            'Flow due to shutter',
            'Decay fits'
            ])
plt.show()


 # Calc gradient directly
 # Find the closest volume index to the target one input
def vol_search(vol, value):
     for v in range(len(vol)-1):
         if vol[v] > value and vol[v+1]<= value:
             return v
     print('No volume index found, returning last volume')
     return len(vol)-1


# remake curves

volume_curve_passive = vol[new_start:new_end]
line_passive = [gradient_passive*V + intercept for V in volume_curve_passive]


if (decay > 0 or gradient_passive > 0):
    gradient_passive = nan
    decay = nan

print('     Passive (Q vs V), gradient is E/R estimate')
print(gradient_passive, intercept, r_value)


gradients.append([gradient_passive])
decays.append([decay])

print('\n Gradients')
for g in gradients:
    print(g)
print('\n Mean gradient [passive, driven]')
print(np.nanmean(gradients, axis=0))
print(np.nanstd(gradients, axis=0))

print('\n decays')
for g in decays:
    print(g)
print('\n Mean  [passive, driven]')
print(np.nanmean(decays, axis=0))
print(np.nanstd(decays, axis=0))

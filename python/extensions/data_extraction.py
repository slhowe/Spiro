#!/bin/bash
'''
Functions in here pull pressure and flow from
matlab files.
'''

import sys
sys.path.insert(0, '/home/sarah/Documents/Spirometry/python/extensions')

from breath_analysis import split_breaths

from numpy import array
from scipy import io
import matplotlib.pyplot as plt

def load_mat_file(filename):
    # Load matlab file
    mat = io.loadmat(filename)
    return mat


def ManualDetection_data(mat, breath):
    # Extract pressure and flow data
    data = mat['ManualDetection']
    pressure = data['Pressure'][0][breath]
    flow = data['Flow'][0][breath]

    # Pressure and flow are in a weird format
    # of array of lists of single values.
    # Convert array to list and pull values out
    # into main list.
    pressure = pressure.tolist()
    pressure = [p[0] for p in pressure] # every item is own list
    flow = flow.tolist()
    flow = [f[0] for f in flow]

    return([pressure, flow])

def dated_data(mat, breath):
    # Extract pressure and flow data
    data = mat['BreathData']
    pressure = data['Pressure'][0][breath]
    flow = data['Flow'][0][breath]

    # Pressure and flow are in a weird format
    # of array of lists of single values.
    # Convert array to list and pull values out
    # into main list.
    pressure = pressure.tolist()
    pressure = [p[0] for p in pressure] # every item is own list
    flow = flow.tolist()
    flow = [f[0] for f in flow]

    return([pressure, flow])


def PS_vs_NAVA_invasive_data(mat):
    # Extract pressure and flow data
    pressure = mat['press']
    flow = mat['flow']

    # Pressure and flow are in a weird format
    # of array of lists of single values.
    # Convert array to list and pull values out
    # into main list.
    pressure = pressure.tolist()
    pressure = [p[0]/100.0 for p in pressure] # every item is own list
    flow = flow.tolist()
    flow = [f[0]/10000.0 for f in flow]

    # All data is in one long list
    # Split into separate breaths
    splits = split_breaths(flow, peak_height=0.1, Fs=50, filt=False, plot=False)
    breath_starts = splits[0]

    # Now make list of lists
    # Where each sublist is a breath
    split_pressure = [0]*(len(breath_starts) - 1)
    split_flow = [0]*(len(breath_starts) - 1)

    i = 0
    while i < (len(breath_starts) - 1):
        start = breath_starts[i]
        end = breath_starts[i+1]
        split_pressure[i] = pressure[start:end]
        split_flow[i] = flow[start:end]
        i+=1

    return([split_pressure, split_flow])

def PS_vs_NAVA_noninvasive_data(mat):
    # Extract pressure and flow data
    data = mat[mat.keys()[0]].tolist() # Really just hoping the first one is data set
    pressure = [data[i][0] for i in range(len(data))]
    flow = [data[i][1] for i in range(len(data))]

    # Pressure and flow are in a weird format
    # of array of lists of single values.
    # Convert array to list and pull values out
    # into main list.
    pressure = [p/10000.0 for p in pressure] # every item is own list
    flow = [f/10000.0 for f in flow]

    # All data is in one long list
    # Split into separate breaths
    splits = split_breaths(flow, peak_height=0.3, Fs=50, filt=True, plot=False)
    breath_starts = splits[0]

    # Now make list of lists
    # Where each sublist is a breath
    split_pressure = [0]*(len(breath_starts) - 1)
    split_flow = [0]*(len(breath_starts) - 1)

    i = 0
    while i < (len(breath_starts) - 1):
        start = breath_starts[i]
        end = breath_starts[i+1]
        split_pressure[i] = pressure[start:end]
        split_flow[i] = flow[start:end]
        i+=1

    return([split_pressure, split_flow])

def sedation_changes_data(mat, key, breath):

    parent_data = mat[key]
    data = parent_data['Breath'][breath][0]

    pressure = [0]*len(data)
    flow = [0]*len(data)
    volume = [0]*len(data)

    for i in range(len(data)):
        pressure[i] = data[i][1]
        flow[i] = data[i][0]
        volume[i] = data[i][2]

    return([pressure, flow, volume])






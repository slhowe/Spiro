#!/bin/bash
'''
Functions in here pull pressure and flow from
matlab files.
'''

from numpy import array
from scipy import io

def load_mat_file(filename):
    # Load matlab file
    mat = io.loadmat(filename)
    data = mat['ManualDetection']
    return data


def ManualDetection_data(data, breath):
    # Extract pressure and flow data
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

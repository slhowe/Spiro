#!/bin/bash

from numpy import nan, matrix, array, log, exp
from numpy.linalg import lstsq
'''
Define data class with a meta class to make
data classes iterable.
'''
class IterRegistry(type):
    def __iter__(cls):
        return iter(cls._registry)

class Data:
    __metaclass__ = IterRegistry
    _registry = []

    SAMPLING_FREQUENCY = 125

    def __init__(self, name):
        self._registry.append(self)
        self.name = name
        self.test_count = 3
        self.finding_flow = False
        self.pressure = []
        self.flow = []
        self.indices = {0:[nan, nan], 1:[nan, nan], 2:[nan, nan]}
        self.decay  = {0:[nan, nan], 1:[nan, nan], 2:[nan, nan]}
        self.offset = {0:[nan, nan], 1:[nan, nan], 2:[nan, nan]}
        self.crossing = [nan, nan, nan]

    def init_start_indices(self, values):
        if len(values) != self.test_count:
            error("invalid start index values")
        else:
            for index in range(self.test_count):
                self.indices[index][0] = values[index]

    def init_end_indices(self, values):
        if len(values) != self.test_count:
            error("invalid end index values")
        else:
            for index in range(self.test_count):
                self.indices[index][1] = values[index]

    def clear_dataset_registry(self):
        del self._registry[:]

def find_drop_index(percentage, dataset, test):
    '''
    Finds the index where curve drops below
    percentage from start. Index relative to
    start of all curve data
    '''
    data_start = dataset.indices[test][0]
    data_stop = dataset.indices[test][1]
    if dataset.finding_flow:
        curve = dataset.flow[data_start:data_stop]
    else:
       curve = dataset.pressure[data_start:data_stop]

    drop = curve[0] - (percentage/100.0)*(curve[0] - curve[-1])

    i = data_start
    still_looking = True
    for value in curve:
        if(still_looking):
            if(value < drop):
                index = i
                still_looking = False
        i += 1
    if(index == data_start):
        error("Percentage drop specified not found in range")

    return index

def exponential_fit(start, stop, data):
    '''
    Exponential fit to data in specified range
    '''

    # Decide if fitting to flow or pressure
    if data.finding_flow:
        curve = data.flow[start:stop]
    else:
        curve = data.pressure[start:stop]

    # Set up array of known values
    Fs = float(data.SAMPLING_FREQUENCY)
    measurements = array([log(curve)])

    # set up array of data multiplied by constants
    one_array = [1 for val in range(len(measurements.T))]
    times = [i/Fs for i in range(len(measurements.T))]
    independents = array([one_array,times])

    # Least squares fit to data
    result = lstsq(independents.T, measurements.T)

    # Squared Euclidean 2-norm for each col in (b - a*x)
    residual = result[1]
    print('       residual on line fit is {}'.format(residual[0]))

    # Parameters - offset and decay constant
    constants = result[0]
    constants[0] = exp(constants[0])
    print('       offset is {}, decay is {}'.format(constants[0], constants[1]))

    return constants

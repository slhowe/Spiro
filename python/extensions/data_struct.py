#!/bin/bash

from numpy import nan
import csv
'''
Define data class with a meta class to make
data classes iterable.
'''
class IterRegistry(type):
    def __iter__(cls):
        return iter(cls._registry)

class DataStore:
    __metaclass__ = IterRegistry
    _registry = []

    def __init__(self, Fs=125):
        ''' Init empty data containers'''
        self._registry.append(self)
        self.pressure = []
        self.flow = []
        self.sampling_frequency = Fs

    def store_data(self, filename):
        ''' Store flow and pressure data from
            csv file in data struct'''
        with open(filename, 'rb') as csvfile:
           reader = csv.reader(csvfile, delimiter = ',')
           # skip header line
           header = reader.next()
           for row in reader:
               pressure = float(row[0])
               flow = float(row[1])
               self.pressure.append(pressure)
               self.flow.append(flow)

    def clear_dataset_registry(self):
        ''' Register must be cleared for use in
            programs with memory after execution
            (eg ipython)'''
        del self._registry[:]

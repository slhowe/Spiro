#!/bin/bash

# Import my extensions
import sys
sys.path.insert(0, '/home/sarah/Documents/Spirometry/python/extensions')

from breath_analysis import calc_flow_delay, split_breaths, BreathData
from filters import hamming
from calculus import integral

import csv
import matplotlib.pyplot as plt
import numpy as np
from numpy import isnan
from numpy import array
from numpy.linalg import lstsq

'''
Have two test sets.
One of just spirometer resistance
One of two spirometers in series
(Resistance doubles)
Testing if RC can be split into sensible result
NOTE: 2 spirometers in series is all kinds of hard to
      breathe through. Need to use lower series resistances
'''

# Locations etc
plot_path = '/home/sarah/Documents/Spirometry/images/exponential_fit'
files = ['test_single_1.csv', 'test_double_1.csv', 'test_double_2.csv']

# Create data classes
class Data:
    def __init__(self):
        self.pressure = []
        self.flow = []
        self.sampling_frequency = 125
        self.offset = []
        self.decay = []

    def init_start_indices(self, values):
        self.start_indices = values

    def init_end_indices(self, values):
        self.end_indices = values

# Some functions
def generate_curve(dataset, start, end, test):
    # Time
    Fs = float(dataset.sampling_frequency)
    times = [x/Fs for x in range((end - start))]
    # start point and decay
    start_point = dataset.offset[test][0]
    decay = dataset.decay[test][0]
    # make line
    curve = [start_point*np.exp(decay*times[x]) for x in range((end - start))]
    return curve

def find_drop_index(percentage, curve):
    '''
    Finds the index where curve drops below
    percentage from start. Index relative to
    start of all curve data
    '''
    curve_max = max(curve)
    curve_max_index = curve.index(curve_max)
    drop = curve_max - (percentage)*(curve_max)

    i = curve_max_index
    index = 0
    while(i < len(curve)):
        value = curve[i]
        if(value < drop):
            index = i
            i = len(curve)
        i += 1
    if(index == 0):
       raise ValueError("Percentage drop specified not found in range")

    return index

def exponential_fit(start, stop, curve):
    '''
    Exponential fit to data in specified range
    '''

    # Chop the set to the fitting range
    curve = curve[start:stop]

    # Set up array of known values
    Fs = float(125)
    measurements = array([np.log(curve)])

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
    constants[0] = np.exp(constants[0])
    print('       offset is {}, decay is {}'.format(constants[0][0], constants[1][0]))

    return constants
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~``
# Breaths to look at in each set
# found them earlier, make it easier
breaths_of_interest = [
                       [1, 4, 7, 10],
                       [1, 3, 6],
                       [1],
                      ]

# Go through each file
for file_index in range(len(files)):
    print('\n~~~ New File ~~~\n')
    # Create data classes
    dataset = Data()

    # Store pressure and flow data
    filename = files[file_index]
    with open(filename, 'rb') as csvfile:
        reader = csv.reader(csvfile, delimiter = ',')

        # skip header line
        header = reader.next()

        for row in reader:
            pressure = float(row[0])
            flow = float(row[1])

            dataset.pressure.append(pressure)
            dataset.flow.append(-1.0*flow)

    pressure = hamming(dataset.pressure, 15, dataset.sampling_frequency, 20, plot=False)
    pressure = np.real(pressure).tolist()
    flow = hamming(dataset.flow, 15, dataset.sampling_frequency, 20, plot=False)
    flow = np.real(flow).tolist()

    # Split all the breaths from the data
    data = BreathData(0)

    # Find starts, middles and ends of flow data
    flow_splits = split_breaths(dataset.pressure, peak_height=0.03, plot=True)
    flow_starts = flow_splits[0]
    flow_middles = flow_splits[1]
    flow_stops = flow_splits[2]

#    dots = [dataset.flow[i] for i in flow_starts]
#    plt.plot(dataset.flow)
#    plt.plot(flow_starts, dots, 'ro')
#    plt.show()

    # calculate decays and offsets for data ranges
    breath = 0
    for breath_num in breaths_of_interest[file_index]:
        # Get data for the breath
        # insp/exp pressure,flow,volume,length
        # time
        data.get_data(pressure,
                      flow,
                      flow_starts[breath_num],
                      flow_middles[breath_num],
                      flow_stops[breath_num],
                      dataset.sampling_frequency
                      )

        # Find Rsp
        Rsp = lstsq(np.array([data.exp_flow]).T, np.array([data.exp_pressure]).T)
        print(Rsp)

        plt.plot(data.insp_pressure)
        plt.plot(data.exp_pressure)
        plt.plot(data.exp_flow, '*')
        plt.plot([Rsp[0][0][0]*f for f in data.exp_flow])
        plt.show()

        print('Fitting to bottom half of data')
        # Crappy data at start, dunno why...
        # plastic flapping or sensor issues..?
        #P_max = max(data.insp_pressure)

        percentage_to_drop = 0.5
        line_fit_start = find_drop_index(percentage_to_drop, data.insp_pressure)

        # Values might go negative too close to end, breaks log(x)
        percentage_to_drop = 0.95
        line_fit_end = find_drop_index(percentage_to_drop, data.insp_pressure)

        plt.plot(data.insp_pressure)
        plt.plot(line_fit_start, data.insp_pressure[line_fit_start], 'go')
        plt.plot(line_fit_end, data.insp_pressure[line_fit_end], 'ro')
        plt.show()

        parameters = exponential_fit(line_fit_start, line_fit_end, data.insp_pressure)
        dataset.offset.append(parameters[0])
        dataset.decay.append(parameters[1])
        print(dataset.offset[-1])
        print(dataset.decay[-1])


        #~~~~~~ Remake curves ~~~~~~
        line_fit_end = data.insp_length
        H_curve = generate_curve(dataset, line_fit_start, line_fit_end, breath)
#        L_curve = generate_curve(dataset, line_fit_start, line_fit_end, breath)
#
#        # Find point where lines cross (relative to line_fit_start)
#        Fs = float(dataset.sampling_frequency)
#        for i in range(0, (line_fit_end - line_fit_start)):
#            if(H_curve[i] < L_curve[i]):
#                dataset.crossing[test] = i
#                break
#        if (isnan(dataset.crossing[test])):
#            print('exponential fit lines never intersect')

        #~~~~~~ plot example ~~~~~~
        plotting = 1
        if(plotting):
            times = [x/float(dataset.sampling_frequency) for x in range(data.insp_length)]

            f, (ax1) = plt.subplots(1, sharex=True)

            decay_info = ('Decay Rates\nHigh: {0:.2f}'.format(dataset.decay[breath][0])
                       + ' \nLow: {0:.2f}'.format(dataset.decay[breath][0]) )
            ax1.text(3.2, 0.49, decay_info, fontsize=12)
            ax1.set_ylabel("Pressure (kPa)")

            ax1.plot(
                    times, data.insp_pressure, 'b',
                    times[line_fit_start:line_fit_end], H_curve, 'r',
                    )
            ax1.legend(['Measured Data','Exponential fit (High)','Exponential fit (Low)'])
            ax1.set_xlabel('Time (s)')
            ax1.grid(True)

            saving = False
            if saving:
                print(fig_name)
                plt.savefig(fig_name)

            plt.show()
            plt.close()

            breath += 1

        print(' ')

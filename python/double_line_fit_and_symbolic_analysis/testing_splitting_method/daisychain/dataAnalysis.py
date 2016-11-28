#!/bin/bash

# Import my extensions
import sys
sys.path.insert(0, '/home/sarah/Documents/Spirometry/python/extensions')

from breath_analysis import calc_flow_delay, split_breaths, BreathData
from calculus import integral

import csv
import matplotlib.pyplot as plt
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


plot_path = '/home/sarah/Documents/Spirometry/images/exponential_fit'
files = ['test_single.csv', 'test_double.csv']

# Create data classes
reset = pe.Data('reset')
reset.clear_dataset_registry()
dataset_1 = pe.Data('dataset_1')
dataset_2 = pe.Data('dataset_2')

# Points I prepared earlier =D
dataset_1.init_start_indices([1190, 3824, 7128])
dataset_2.init_start_indices([2475, 5737, 8939])
dataset_1.init_end_indices([2000, 4634, 7938])
dataset_2.init_end_indices([3000, 6150, 9340])

finding_flow = False
dataset_1.finding_flow = finding_flow
dataset_2.finding_flow = finding_flow
if finding_flow:
    delay = 26;
    dataset_1.init_start_indices([1190+delay, 3824+delay, 7128+delay])
    dataset_2.init_start_indices([2475+delay, 5737+delay, 8939+delay])
    dataset_1.init_end_indices([2000+delay, 4634+delay, 7938+delay])
    dataset_2.init_end_indices([3000+delay, 6150+delay, 9340+delay])

# Store pressure and flow data
file_index = 0
for dataset in pe.Data:
    filename = path + files[file_index]
    with open(filename, 'rb') as csvfile:
        reader = csv.reader(csvfile, delimiter = ',')

        # skip header line
        header = reader.next()

        for row in reader:
            pressure = float(row[0])
            flow = float(row[1])

            dataset.pressure.append(pressure)
            dataset.flow.append(-1.0*flow)
    file_index += 1

def generate_curve(dataset, start, end, offset, test, index):
    Fs = float(dataset.SAMPLING_FREQUENCY)
    times = [x/Fs for x in range((start - offset),(end - offset))]
    start_point = dataset.offset[test][index]
    decay = dataset.decay[test][index]
    curve = [start_point*exp(times[x]*decay) for x in range(0, (end - start))]
    return curve

def remake_input_pressure():
    pass

# calculate decays and offsets for data ranges
test = 0
for dataset in pe.Data:
    print('~~~New Data Set~~~')
    for test in range(dataset.test_count):
        print('Test {}'.format(test+1))

        curve_start = dataset.indices[test][0]
        curve_end = dataset.indices[test][1]

        # First the top curve
        print('Fitting to top half of decay')
        percentage_to_drop = 50
        line_fit_end = pe.find_drop_index(percentage_to_drop, dataset, test)
        parameters = pe.exponential_fit(curve_start, line_fit_end, dataset)
        dataset.offset[test][0] = parameters[0]
        dataset.decay[test][0] = parameters[1]

        print('Fitting to bottom half of data')
        # Values go negative too close to end, breaks log(x)
        percentage_to_drop = 90
        line_fit_end = pe.find_drop_index(percentage_to_drop, dataset, test)
        line_fit_start = ((line_fit_end - curve_start)/2) + curve_start
        parameters = pe.exponential_fit(line_fit_start, line_fit_end, dataset)
        dataset.offset[test][1] = parameters[0]
        dataset.decay[test][1] = parameters[1]
        print(line_fit_start)


        #~~~~~~ Remake curves ~~~~~~
        H_curve = generate_curve(dataset, curve_start, curve_end, curve_start, test, 0)
        L_curve = generate_curve(dataset, curve_start, curve_end, line_fit_start, test, 1)

        # Find point where lines cross (relative to curve_start)
        Fs = float(dataset.SAMPLING_FREQUENCY)
        for i in range(0, (curve_end - curve_start)):
            if(H_curve[i] < L_curve[i]):
                dataset.crossing[test] = i
                break
        if (isnan(dataset.crossing[test])):
            print('exponential fit lines never intersect')


        #~~~~~~ plot example ~~~~~~
        plotting = True
        if(plotting):
            times = [x/Fs for x in range(0,(curve_end - curve_start))]

            f, (ax1) = plt.subplots(1, sharex=True)

            decay_info = 'Decay Rates\nHigh: {0:.2f}'.format(dataset.decay[test][0][0]) + ' \nLow: {0:.2f}'.format(dataset.decay[test][1][0])
            if dataset.finding_flow:
                curve = dataset.flow[curve_start:curve_end]
                ax1.text(3.2, 3.9, decay_info, fontsize=12)
                ax1.set_ylabel("Flow (L/s)")
                fig_name = '{}Flow_{}_{}.png'.format(plot_path, dataset.name, test+1)
            else:
                curve = dataset.pressure[curve_start:curve_end]
                ax1.text(3.2, 0.49, decay_info, fontsize=12)
                ax1.set_ylabel("Pressure (kPa)")
                fig_name = '{}Pressure_{}_{}.png'.format(plot_path, dataset.name, test+1)

            ax1.plot(times, curve, 'b', times, H_curve, 'r', times, L_curve, 'g')
            ax1.legend(['Measured Data','Exponential fit (High)','Exponential fit (Low)'])
            ax1.set_xlabel('Time (s)')
            ax1.grid(True)

            saving = False
            if saving:
                print(fig_name)
                plt.savefig(fig_name)

            plt.show()
            plt.close()

        print(' ')

    test += 1


# Clear all data
dataset_1.clear_dataset_registry()


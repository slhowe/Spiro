#!/bin/bash

import csv
import dill
import matplotlib.pyplot as plt
import pressureEquation as pe
import symbolicPressureEquation as spe
from sympy import *
from sympy.abc import s,t
from numpy import isnan

path = '/home/sarah/Documents/Spirometry/data/'
plot_path = '/home/sarah/Documents/Spirometry/lungPlots/'
files = ['Loops_1.csv', 'Loops_3.csv']

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

## Generate symbolic equations
#y = Function("y")(t)
#u = Function("u")(t)
#
#A1 = Symbol("A1")
#A2 = Symbol("A2")
#P1 = Symbol("P1")
#P2 = Symbol("P2")
#
## x is the time that the curves cross
#x = Symbol("x")
#
## Substitute into EoM_y here and get values
#y_of_t_minus_2= Symbol("y2t")
#y_of_t_minus_1 = Symbol("y1t")
#y_of_t = Symbol("yt")
#y_of_t_plus_1= Symbol("yt1")
#u_of_t = Symbol("ut")
#dt = Symbol("dt")
#
#print("Loading pre-calculated equations")
#[y_equation, u_equation] = dill.load(open('symbolicEquations.pk'))
#
#print(y_equation)
#print(u_equation)

def generate_curve(dataset, start, end, offset, test, index):
    Fs = float(dataset.SAMPLING_FREQUENCY)
    times = [x/Fs for x in range((start - offset),(end - offset))]
    start_point = dataset.offset[test][index]
    decay = dataset.decay[test][index]
    curve = [start_point*exp(times[x]*decay) for x in range(0, (end - start))]
    return curve

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

        # Remake curves
        H_curve = generate_curve(dataset, curve_start, curve_end, curve_start, test, 0)
        L_curve = generate_curve(dataset, curve_start, curve_end, line_fit_start, test, 1)

        # Find point where lines cross (relative to curve_start)
        Fs = float(dataset.SAMPLING_FREQUENCY)
        for i in range(0, (curve_end - curve_start)):
            if(H_curve[i] < L_curve[i]):
                dataset.crossing[test] = i
                break
        if (isnan(dataset.crossing[test])):
            error('exponential fit lines never intersect')

        # plot example
        plotting = True
        if(plotting):
            times = [x/Fs for x in range(0,(curve_end - curve_start))]

            decay_info = 'Decay Rates\nHigh: {0:.2f}'.format(dataset.decay[test][0][0]) + ' \nLow: {0:.2f}'.format(dataset.decay[test][1][0])
            if dataset.finding_flow:
                curve = dataset.flow[curve_start:curve_end]
                plt.axis([-0.01, 4, -0.01, 6])
                plt.text(3.2, 3.9, decay_info, fontsize=12)
                #crossing_info = 'Exp crossing: {0:.1f} % drop'.format((1-curve[dataset.crossing[test]]/curve[0])*100)
                #plt.text(2.2, 2.9, crossing_info, fontsize=12)
                plt.ylabel("Flow (L/s)")
                fig_name = '{}Flow_{}_{}.png'.format(plot_path, dataset.name, test+1)
            else:
                curve = dataset.pressure[curve_start:curve_end]
                plt.axis([-0.01, 4, -0.01, 0.75])
                plt.text(3.2, 0.49, decay_info, fontsize=12)
                #crossing_info = 'Exp crossing: {0:.1f} % drop'.format((1-curve[dataset.crossing[test]]/curve[0])*100)
                #plt.text(2.2, 0.39, crossing_info, fontsize=12)
                plt.ylabel("Pressure (cmH20)")
                fig_name = '{}Pressure_{}_{}.png'.format(plot_path, dataset.name, test+1)

            plt.plot(times, curve, 'b', times, H_curve, 'r', times, L_curve, 'g')
            plt.legend(['Measured Data','Exponential fit (High)','Exponential fit (Low)'])
            plt.xlabel('Time (s)')
            plt.grid(True)

            saving = True
            if saving:
                print(fig_name)
                plt.savefig(fig_name)

            plt.show()
            plt.close()

        print(' ')

    test += 1


# Clear all data
dataset_1.clear_dataset_registry()


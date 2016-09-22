#!/bin/bash

# Import my extensions
import sys
sys.path.insert(0, '/home/sarah/Documents/Spirometry/python/extensions')

import csv
import dill
import matplotlib.pyplot as plt
import pressureEquation as pe
from numpy import isnan
from math import atan, sin, cos, pi
from calculus import integral
from filters import semi_gauss_lp_filter
from breath_analysis import calc_flow_delay, split_breaths, BreathData
from numpy import array
from numpy import exp
from numpy.linalg import lstsq
from numpy.fft import ifft, fft

INSP = 0
EXP = 1
Fs = 125

path = '/home/sarah/Documents/Spirometry/data/'
plot_path = '/home/sarah/Documents/Spirometry/images/exponential_fit'
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

    for breath in range(2):#breaths):
        # Filter data and flip pressure
        pressure = semi_gauss_lp_filter(dataset.pressure,
                                        dataset.SAMPLING_FREQUENCY,
                                        3,
                                        plot=False)
        flip_pressure = [-p for p in pressure]

        # Split breaths and setup automatic alignment in storage
        delay = calc_flow_delay(flip_pressure,
                                dataset.flow,
                                dataset.SAMPLING_FREQUENCY,
                                plot=False
                                )
        print('Delay: {}'.format(delay))
        flow_splits = split_breaths(dataset.flow)
        breaths = len(flow_splits[0])
        data = BreathData(delay)

        data.get_data(pressure,
                      dataset.flow,
                      flow_splits[0][breath],
                      flow_splits[1][breath],
                      flow_splits[2][breath],
                      dataset.SAMPLING_FREQUENCY
                      )

        breath_pressure = data.insp_pressure + data.exp_pressure
        breath_flow = data.insp_flow + data.exp_flow
        breath_volume = data.insp_volume + data.exp_volume

        # Change the volume offset
        # so int(I) = vol
        # with vol0 at I(T/4) = insp/2
        offset_vol = data.insp_volume[data.insp_length/2]
        breath_volume = [v - offset_vol for v in breath_volume]

        # Find spirometer resistance
        dependent = array(breath_pressure)
        independent = array([breath_flow])
        res = lstsq(independent.T, dependent.T)
        Rsp = res[0][0]
        print('Rsp: {}'.format(Rsp))

        # What is the frequency?
        insp_f = 1/(2*data.insp_length/125.0)
        exp_f = 1/(2*data.exp_length/125.0)
        mean_f = (insp_f + exp_f)/2

        # Guess use the lower RC value for normal breathing
        # That way it matches with the flow decay too
        mean_RC = (dataset.decay[0][1][0]+dataset.decay[1][1][0]+dataset.decay[2][1][0])/3

        # Crazy way to look at w and theta
        wi = 2*pi*insp_f
        we = 2*pi*exp_f
        mean_w = 2*pi*mean_f
        thetai = atan(wi*mean_RC)
        thetae = atan(we*mean_RC)
        theta = atan(mean_w*mean_RC)

        # A ridiculous way to shift the signal by a few data points
        #shift = int(mean_f*125.0) #<-- Should just move everything by this
        shift = int(pi/2.0*len(breath_pressure)/2/pi)
        the_input = breath_pressure[-shift:] + breath_pressure[:-shift]
        #fft_p = fft(breath_pressure)
        #N = len(fft_p)
        #k = [r-len(fft_p)/2 for r in range(len(fft_p))]
        #shifter = [-1j*2*pi*shift*v/N for v in k]
        #freq_data = fft_p*exp(shifter)
        #the_input = ifft(freq_data)

        scaler = max(breath_pressure)

        for iteration in range(1):
            # Increase input to amplitude of 1
            input_est = [v/scaler for v in the_input]

            # Least squares
            # Pin = EV + RQ
            # Get E, R for input guess
            # know value of RQ
            # so scale Pin until RQ matches calculation
            dependent = array([input_est])
            independent = array([breath_volume, breath_flow])
            res = lstsq(independent.T, dependent.T)
            print(res)
            E = res[0][0][0]
            R = res[0][1][0]
            print('E: {}'.format(E))
            print('R: {}'.format(R))
            EoR = E/R
            print('E/R: {}'.format(EoR))

            # Simulate the total airway and breath breath pressure
            paw_sim = [R*breath_flow[i]
                            for i in range(len(breath_flow))]

            # Get spir pressure sim from total airway pressure
            factor = Rsp/R
            breath_pressure_sim = [factor*p for p in paw_sim]

            #~~~ plots ~~~
            line_in = [E*breath_volume[i] + R*breath_flow[i]
                            for i in range(len(breath_flow))]

            plt.plot(the_input, 'r')
            plt.plot(line_in, 'm')
            plt.plot(breath_pressure, 'g')
            plt.plot(breath_pressure_sim, 'xy')
            plt.plot(breath_volume)
            plt.legend(['shifted output',
                        'input sim',
                        'output',
                        'output sim',
                        'volume'
                        ])
            plt.grid()
            plt.show()

        #~~~~~~ plot example ~~~~~~
        plotting = False
        if(plotting):
            times = [x/Fs for x in range(0,(curve_end - curve_start))]

            f, (ax1, ax2) = plt.subplots(2, sharex=False)

            decay_info = ('Decay Rates\nHigh: {0:.2f}, {0:.2f}, {0:.2f}'.format(
                                                    dataset.decay[0][0][0],
                                                    dataset.decay[1][0][0],
                                                    dataset.decay[2][0][0]
                                                    )
                       + ' \nLow: {0:.2f}, {0:.2f}, {0:.2f}'.format(
                                                    dataset.decay[0][1][0],
                                                    dataset.decay[1][1][0],
                                                    dataset.decay[2][1][0]
                                                    ))

            curve = dataset.pressure[curve_start:curve_end]
            ax1.text(3.2, 0.49, decay_info, fontsize=12)
            ax1.set_ylabel("Pressure (kPa)")
            fig_name = '{}Pressure_{}_{}.png'.format(plot_path, dataset.name, test+1)

            ax1.plot(times, curve, 'b', times, H_curve, 'r', times, L_curve, 'g')
            ax1.legend(['Measured Data','Exponential fit (High)','Exponential fit (Low)'])
            ax1.set_xlabel('Time (s)')
            ax1.grid(True)

            times = [x/Fs for x in range(len(breath_pressure))]
            ax2.plot(times, the_input, 'r')
            ax2.plot(times, input_est, 'm')
            ax2.plot(times, breath_pressure, 'y')
            ax2.legend(['Input pressure',
                        'Input scaled',
                        'Output pressure',
                        ])
            ax2.grid()

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


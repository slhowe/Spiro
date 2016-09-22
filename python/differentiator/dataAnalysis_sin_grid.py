#!/bin/bash

# Import my extensions
import sys
sys.path.insert(0, '/home/sarah/Documents/Spirometry/python/extensions')

import csv
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from math import atan, sin, cos, pi
from calculus import integral
from filters import semi_gauss_lp_filter
from breath_analysis import calc_flow_delay, split_breaths, BreathData
from numpy import array, exp, isnan
from numpy.linalg import lstsq
from numpy.fft import ifft, fft

INSP = 0
EXP = 1
Fs = 190

path = '/home/sarah/Documents/Spirometry/data/'
plot_path = '/home/sarah/Documents/Spirometry/images/exponential_fit'
files = ['slow.csv','generated_breath.csv', 'Loops_1.csv', 'Loops_3.csv']

# Create data classes
class Data:
    def __init__(self, Fs):
        self.pressure = []
        self.flow = []
        self.SAMPLING_FREQUENCY = Fs

# Store pressure and flow data
dataset = Data(Fs)
#filename = path + files[0]
filename = files[0]

with open(filename, 'rb') as csvfile:
    reader = csv.reader(csvfile, delimiter = ',')

    # skip header line
    header = reader.next()

    for row in reader:
        pressure = float(row[0])
        flow = float(row[1])

        dataset.pressure.append(pressure)
        dataset.flow.append(1.0*flow)

# Filter data and flip pressure
pressure = semi_gauss_lp_filter(dataset.pressure,
                                dataset.SAMPLING_FREQUENCY,
                                3,
                                plot=False)

# Split breaths and setup automatic alignment in storage
delay = calc_flow_delay(pressure,
                        dataset.flow,
                        dataset.SAMPLING_FREQUENCY,
                        plot=False
                        )
print('Delay: {}'.format(delay))

pressure = [-p for p in pressure]
flow_splits = split_breaths(dataset.flow)
breaths = len(flow_splits[0])

print(flow_splits)

for breath in range(1):#breaths):
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

    plt.plot(breath_flow)
    plt.plot(breath_pressure)
    plt.show()

    # Change the volume offset
    # So even about zero
    offset_vol = data.insp_volume[data.insp_length/2]
    data.insp_volume = [v - offset_vol for v in data.insp_volume]

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
    w = 2*pi*insp_f
    wt = [w*t/float(Fs) for t in range(data.insp_length)]

    # Crazy way to look at w and theta
    thetas = [76, 78, 80, 82, 84, 86, 88]
    mags = [0.2, 0.9, 1.4, 2]
    errors = []
    x = []
    y = []

    for theta in thetas:
        # A ridiculous way to shift the signal by a few data points
        shift = theta*pi/180.0
        the_input = [sin(WT - shift) for WT in wt]
        print('~~~ {} ~~~'.format(shift))

        #normaliser = max(abs(min(the_input)), abs(max(the_input)))
        #the_input = [v/normaliser for v in the_input]

        for scaler in mags:
            x.append(theta)
            y.append(scaler)

            # Increase input to amplitude of 1
            input_est = [v*scaler for v in the_input]

            # Least squares
            # Pin = EV + RQ
            # Get E, R for input guess
            # know value of RQ
            # so scale Pin until RQ matches calculation
            dependent = array([input_est])
            independent = array([data.insp_volume, data.insp_flow])
            res = lstsq(independent.T, dependent.T)
            print(res)
            E = res[0][0][0]
            R = res[0][1][0]
            print('E: {}'.format(E))
            print('R: {}'.format(R))
            EoR = E/R
            print('E/R: {}'.format(EoR))

            # Simulate the total airway and breath breath pressure
            paw_sim = [R*data.insp_flow[i]
                            for i in range(len(data.insp_flow))]

            # Get spir pressure sim from total airway pressure
            factor = Rsp/R
            breath_pressure_sim = [factor*p for p in paw_sim]

            #~~~ plots ~~~
            line_in = [E*data.insp_volume[i] + R*data.insp_flow[i]
                            for i in range(data.insp_length)]

            line_R = [R*data.insp_flow[i]
                            for i in range(data.insp_length)]

            line_Rsp = [Rsp/R*v for v in line_R]

            # Look at error in Rsp
            error = [(line_Rsp[i] - data.insp_pressure[i])**2
                            for i in range(data.insp_length)]
            mean_squared_error = 0
            for e in error:
                mean_squared_error += e/len(error)

            errors.append(mean_squared_error)

            #time = [t/float(Fs) for t in range(data.insp_length)]
            #plt.plot(time, input_est, 'r')
            #plt.plot(time, line_in, 'm')
            #plt.plot(time, data.insp_pressure, 'g')
            #plt.plot(time, data.insp_volume, 'xb')
            #plt.plot(time, line_R, 'xk')
            #plt.plot(time, line_Rsp, 'gx')
            #plt.legend(['input guess',
                        #'input back simulated',
                        #'output measured (PRsp)',
                        #'volume',
                        #'PR est',
                        #'PRsp est'
                        #])
            #plt.grid()
            #plt.show()

    print(errors)
    fig = plt.figure()
    ax = Axes3D(fig)
    ax.scatter(x, y, errors)
    ax.set_zlabel('residuals in output')
    ax.set_xlabel('theta')
    ax.set_ylabel('mag')
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


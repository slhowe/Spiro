#!/bin/bash

# Import my extensions
import sys
sys.path.insert(0, '/home/sarah/Documents/Spirometry/python/extensions')

import csv
import dill
import matplotlib.pyplot as plt
from numpy import isnan
from calculus import integral
from data_struct import DataStore
from filters import semi_gauss_lp_filter
from breath_analysis import BreathData, split_breaths, calc_flow_delay
from numpy import array
from numpy import exp
from numpy.linalg import lstsq

path = '/home/sarah/Documents/Spirometry/data/'
plot_path = '/home/sarah/Documents/Spirometry/images/exponential_fit'
files = ['generated_breath.csv', 'Normal_1.csv', 'Banding_1.csv', 'Loops_3.csv']

# Create data classes
reset = DataStore()
reset.clear_dataset_registry()

dataset = DataStore(Fs=125)
dataset.store_data(files[0])

# Filter data and flip pressure
pressure = semi_gauss_lp_filter(dataset.pressure,
                                dataset.sampling_frequency,
                                100,
                                plot=True)
###pressure = [-p for p in pressure]

# Split breaths and setup automatic alignment in storage
delay = calc_flow_delay(pressure,
                        dataset.flow,
                        dataset.sampling_frequency,
                        plot=False
                        )
print('Delay: {}'.format(delay))
flow_splits = split_breaths(dataset.flow)
breaths = len(flow_splits[0])
data = BreathData(delay)

#~~~ Look at data ~~~
breaths = 1
for breath in range(breaths):
    print('\n\n~~~ New Breath ({}/{}) ~~~'.format(breath+1, breaths))

    data.get_data(pressure,
                dataset.flow,
                flow_splits[0][breath],
                flow_splits[1][breath],
                flow_splits[2][breath],
                dataset.sampling_frequency
                )

    breath_pressure = data.insp_pressure + data.exp_pressure
    breath_flow = data.insp_flow + data.exp_flow
    breath_volume = data.insp_volume + data.exp_volume

    # Find spirometer resistance
    dependent = array(breath_pressure)
    independent = array([breath_flow])
    res = lstsq(independent.T, dependent.T)
    Rsp = res[0][0]
    print('Rsp: {}'.format(Rsp))

    # Change the volume offset
    # so int(I) = vol
    # with vol0 at I(T/4) = insp/2
    offset_vol = data.insp_volume[data.insp_length/2]
    breath_volume = [v - offset_vol for v in breath_volume]

    # ~~~ Work in halves ~~~
    # initial guess
    # E/R = 1/RC
    RC = 1

    #~~~ PLOTS ~~~
    p_in = integral(breath_pressure,
                    dataset.sampling_frequency
                    )
    p_in = [p/RC for p in p_in]

    # Shift so P_in(0) at P_out(insp/2)
    offset_pressure = data.insp_pressure[data.insp_length/2]
    p_in = [p - offset_pressure for p in p_in]

#    plt.plot(breath_flow, 'b')
#    plt.plot(breath_volume, '.r')
#    plt.plot(p_in, '.k')
#    plt.plot(breath_pressure, '.g')

    p_act = integral(breath_pressure,
                    dataset.sampling_frequency
                    )
    p_act = [1.441/0.061*p/2.76 for p in p_act]

    # Shift so P_in(0) at P_out(insp/2)
    offset_pressure = p_act[data.insp_length/2]
    p_act = [p - offset_pressure for p in p_act]

    dependent = array(p_act)
    independent = array([breath_volume, breath_flow])
    res = lstsq(independent.T, dependent.T)
    print(res)

    for breath_half in range(2):
        if(breath_half == 0):
            print('\nInspiration:')
            half_pressure = p_in[0:data.insp_length]
            half_flow = breath_flow[0:data.insp_length]
            half_volume = breath_volume[0:data.insp_length]
            half_paw = breath_pressure[0:data.insp_length]
        else:
            print('\nExpiration:')
            half_pressure = p_in[data.insp_length:]
            half_flow = breath_flow[data.insp_length:]
            half_volume = breath_volume[data.insp_length:]
            half_paw = breath_pressure[data.insp_length:]

        #~~~ Start Iteration ~~~
        for iteration in range(3):

            # Calc Cr, Rsp
            dependent = array(half_pressure)
            print(dependent.shape)
            independent = array([half_volume, half_flow])
            print(independent.shape)
            res = lstsq(independent.T, dependent.T)
            Er = res[0][0]
            Rr = res[0][1]

            print('Er, Rr: {}'.format(res))
            RC = Rr/Er
            print('E/R: {}'.format(1/RC))

            #~~~ Plots ~~~
            plt.plot(p_act, 'k')
            plt.plot(half_pressure, 'b')
            plt.plot(half_flow, 'r')
            plt.plot(half_volume, 'm')
            plt.plot(half_paw, 'y')

            p_lung = [Er*v for v in half_volume]
            p_aw = [Rr*q for q in half_flow]
            p_in_guess = [p_lung[i] + p_aw[i] for i in range(len(p_aw))]

            plt.plot(p_lung, '.g')
            plt.plot(p_aw, '.y')
            plt.plot(p_in_guess, '.b')

            if Rr/Rsp > 1.1 or Rr/Rsp < 0.9:
                Rr = Rsp
                Er = 1/RC*Rr
                print('E/R updated: {}'.format(1/RC))
                p_aw = [Rsp*q for q in half_flow]
                p_lung = [Er*v for v in half_volume]
                half_pressure = [p_lung[i] + p_aw[i] for i in range(len(half_flow))]
                plt.plot(half_pressure, 'xb')
            else:
                # p_in = (E/R)* integral(p_out)
                p_in = integral(breath_pressure,
                                dataset.sampling_frequency
                                )
                p_in = [p/RC for p in p_in]

                # Shift so P_in(0) at P_out(insp/2)
                offset_pressure = p_in[data.insp_length/2]
                p_in = [p - offset_pressure for p in p_in]
                if(breath_half==0):
                    half_pressure = p_in[0:data.insp_length]
                else:
                    half_pressure = p_in[data.insp_length:]

            plt.show()

    #~~~ PLOTS ~~~
    plt.plot(breath_flow, 'b')
    plt.plot(breath_volume, '.r')
    plt.plot(breath_pressure, '.g')
    plt.plot(p_in, '.k')
    plt.show()

    saving = False
    if saving:
        print(fig_name)
        plt.savefig(fig_name)

    plt.show()
    plt.close()

# clear all data
dataset.clear_dataset_registry()


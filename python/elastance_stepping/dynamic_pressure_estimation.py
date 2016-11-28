#!/bin/bash

# Add my extensions to path
import sys
sys.path.insert(0, '/home/sarah/Documents/Spirometry/python/extensions')

# Import my extensions
from data_struct import DataStore
from breath_analysis import BreathData, split_breaths, calc_flow_delay
from calculus import integral, derivative
from filters import hamming
from resistance_plots import ResistanceTable

#Import built-ins
import matplotlib.pyplot as plt
from math import sin, cos, pi, exp, log
from numpy import isnan, array, real
from numpy.linalg import lstsq

def analysis(filename, plot=False, plot_path='\0'):
    ''' Using data files:
        find E = (resistance*flow)/(volume)

        Plots are saved if plot_path given
        Data in kPa and L and s'''

    # Definitions
    INSP = 0
    EXP = 1
    sampling_frequency = 125

    # Remove any old data sets hanging around
    # (Needed when using ipython)
    reset = DataStore()
    reset.clear_dataset_registry()

    # Create data storage
    # Store pressure and flow data from files
    dataset = DataStore(Fs=sampling_frequency)
    dataset.store_data(filename)

    # Filter data
    # Use hamming filter
    # data, Fc, Fs, bw
    flow = hamming(dataset.flow, 50, 125, 20, plot=False)
    flow = real(flow).tolist()
    pressure = dataset.pressure

    # Split flow into breaths
    # Find start, middle and end of each breath
    flow_splits = split_breaths(dataset.flow)
    flow_starts = flow_splits[0]
    flow_middles = flow_splits[1]
    flow_stops = flow_splits[2]

    # Set up storage for breath data
    # Corrects for pressure delay
    # Not using pressure here so no delay needed
    data = BreathData(0)

    print(len(flow_starts))
    print(flow_starts)
    RC_all = [0]*(2*len(flow_starts))
    vol_all = [0]*len(flow)

    # Iterate through breaths
    print('Remember:\nRC should be NEGATIVE\n*** V = PC - RCQ ***\nand we find RC assuming positive')
    for breath in range(1):#len(flow_starts)):
        # Extract breath data
        data.get_data(pressure,
                      flow,
                      flow_starts[breath],
                      flow_middles[breath],
                      flow_stops[breath],
                      dataset.sampling_frequency
                      )

        vol_all[flow_starts[breath]:flow_middles[breath]] = data.insp_volume
        vol_all[flow_middles[breath]:flow_stops[breath]] = data.exp_volume

        # Go through each breath one half at a time
        # First: inspiration
        # Second: expiration
        for section in range(2):
            # Set data for insp/exp here
            if(section == INSP):
                print('\nInsp parameters:')
                breath_length = data.insp_length
                vol = data.insp_volume
                flw = data.insp_flow
                time = data.time[:breath_length]
            elif(section == EXP):
                print('\nExp parameters:')
                breath_length = data.exp_length
                vol = data.exp_volume
                flw = data.exp_flow
                time = data.time[-breath_length:]

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#            y = range(450)
#            y = [i*0.01 for i in y]
#
#            flw = [3*exp(-i) for i in (y)]
#            vol = integral(flw, 125)
#            breath_length = len(y)
#            time = y
#
#            plt.plot(time, flw)
#            plt.show()
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

            # Storage for Parameter ID
            PC = [0]*breath_length
            RC = [0]*breath_length
            res = [0]*breath_length

            # Model pressure input with
            #  1) ln curve
            #  2) exponential decay
            #  3) constant
            #  4) linear section after 70%
            #  5) exponential growth

            # ln screws up for 0
            # Set to small number instead
            if(time[0] == 0):
                time[0] = 1e-3
            ln_curve = [log(t) for t in time]

            # Use same rate for both exponentials
            rate = .8
            decay = [exp(-rate*t) for t in time]
            growth = [exp(rate*t) for t in time]

            # Constant the whole time
            constant = [1]*breath_length

            # zero until a cutoff is met the linear
            cutoff = int(0.7*breath_length)
            linear = [0]*cutoff + range(breath_length-cutoff)

            # Least squares
            dependent = array([vol])
            independent = array([flw, constant])#, decay])#ln_curve, decay, constant, linear, growth])
            res = lstsq(independent.T, dependent.T)
            RC = res[0][0][0]
            a1 = res[0][1][0]
            #a2 = res[0][2][0]
            #a3 = res[0][3][0]
            #a4 = res[0][4][0]
            #a5 = res[0][5][0]

            print('RC: {}'.format(RC))
            print('a1: {}'.format(a1))
            #print('a2: {}'.format(a2))
            #print('a3: {}'.format(a3))
            #print('a4: {}'.format(a4))
            #print('a5: {}'.format(a5))

            RC_all[2*breath+section] = RC

            # Remake pressure/elastance
            PC = [#a1*ln_curve[i]
                   #+a2*decay[i]
                   +a1*constant[i]
                   #+a4*linear[i]
                   #+a5*growth[i]
                   for i in range(breath_length)]

            if(0):
                #c1 = [a1*ln_curve[i] for i in range(breath_length)]
                #c2 = [a2*decay[i] for i in range(breath_length)]
                c3 = [a1*constant[i] for i in range(breath_length)]
                #c4 = [a4*linear[i] for i in range(breath_length)]
                #c5 = [a5*growth[i] for i in range(breath_length)]

                plt.plot(PC)
                #plt.plot(c1, '--')
                #plt.plot(c2, '--')
                plt.plot(c3, '--')
                #plt.plot(c4, '--')
                #plt.plot(c5, '--')
                plt.legend(['PC', 'constant', 'decay'])#'ln','decay','const','linear','growth'])
                plt.grid()
                plt.show()

            if(plot):
                f, (ax2, ax3, ax4, ax5) = plt.subplots(4, sharex=True)
                ax2.set_title('Half Breath Data')

                ax2.plot(flw)
                ax2.set_ylabel('Flow (L/s)')
                ax2.grid()

                fwd_vol = [PC[i] + RC*flw[i] for i in range(breath_length)]
                ax3.plot(vol)
                ax3.plot(fwd_vol)
                ax3.set_ylabel('Volume (L)')
                ax3.legend(['Orig', 'Sim'])
                ax3.grid()

                PC_guess = [vol[i] + 0.3*flw[i] for i in range(len(vol))]
                ax4.plot(PC_guess)
                ax4.plot(PC)
                ax4.legend(['Direct guess (RC = 0.3)', 'Sim'])
                ax4.set_ylabel(['Pressure/Elastance (kPa/E)'])
                ax4.grid()

                ax5.plot([RC]*breath_length)
                ax5.set_ylabel(['Resistance/Elastance (R/E)'])
                ax5.grid()

                plt.show()

                # Plot saved if plot path given
                if(plot_path != '\0'):
                    fig_name = filename + '_iteration_{}'.format(iteration)
                    print(fig_name)
                    plt.savefig(fig_name)

                plt.close()

    t = [i*0.1 for i in range(100)]
    plt.plot([40 - 40*exp(-2*i) for i in t])

    f, (ax1, ax2, ax3) = plt.subplots(3, sharex=True)
    ax1.plot(flow_middles, RC_all[1::2], 'oc')
    ax1.plot(flow_middles, RC_all[0::2], 'xr')
    ax1.set_title('RC values')
    ax1.legend(['exh', 'inh'])
    ax1.grid()
    ax2.plot(flow)
    ax2.plot(flow_starts, [0]*len(flow_starts), 'og')
    ax2.plot(flow_stops, [0]*len(flow_stops), 'or')
    ax2.legend(['Data', 'Start point', 'End point'])
    ax2.grid()
    ax3.plot(vol_all)
    ax3.plot([0]*len(vol_all), 'r')
    ax3.set_ylabel('volume')
    ax3.grid()
    plt.show()

    # Final clean up
    dataset.clear_dataset_registry()

    return(PC, RC)

if __name__ == "__main__":
    # Data path
    path = '/home/sarah/Documents/Spirometry/data/'
    #files = ['Normal_1.csv']
    files = ['Loops_1.csv']
    files = [path + name for name in files]

    # Place to save plots
    plot_path = '/home/sarah/Documents/Spirometry/lungPlots/'

    filename = files[0]
    res = analysis(filename, plot=True)

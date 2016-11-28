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
from math import sin, cos, pi, exp
from numpy import isnan, array, real
from numpy.linalg import lstsq

class Chunk:
    def __init__(self, count, size, sampling_frequency):
        self.size = size
        self.volume = [0]*self.size
        self.flow = [0]*self.size
        self.resistance = [0]*self.size
        self.PC = 0
        self.RC = 0
        self.res = 0
        self.start = count*self.size
        self.Fs = sampling_frequency

    def get_PC_and_RC(self, vol, flw):
        # find the volume and flow in the chunk
        self.volume = vol[self.start : (self.start+self.size)]
        self.flow = flw[self.start : (self.start+self.size)]
        ones = [1]*self.size

        # Calc PC and RC for chunk
        # V = PC - RCQ
        dependent = array([self.volume])
        independent = array([ones, self.flow])
        res = lstsq(independent.T, dependent.T)

        # Results
        self.PC = res[0][0][0]
        self.RC = res[0][1][0]
        self.res = res[1]

    def try_another_thing(self, vol):
        time = range(self.start, self.start+self.size)
        time = [t/125.0 for t in time]

        ones = [1]*self.size
        self.volume = vol[self.start : (self.start+self.size)]

        dependent = array([self.volume])
        independent = array([ones, time])
        res = lstsq(independent.T, dependent.T)

        self.PC = exp(res[0][0][0])
        print(self.PC)
        self.RC = 1/(res[0][1][0]/res[0][0][0])
        print(self.RC)
        print(res)

def analysis(filename, plot=False, plot_path='\0'):
    ''' Using data files:
        find E = (resistance*flow)/(volume)

        Plots are saved if plot_path given
        Data in kPa and L and s'''

    # multiply kPa to get cmH20
    kPa_to_cmH20 = 10.197
    sampling_frequency = 125
    INSP = 0
    EXP = 1
    vol_offset = 0
    pres_offset = 0

    # Remove any old data sets hanging around
    # (Needed when using ipython)
    reset = DataStore()
    reset.clear_dataset_registry()

    # Create data storage
    # Store pressure and flow data from files
    dataset = DataStore(Fs=sampling_frequency)
    dataset.store_data(filename)

    print(len(dataset.flow))

    # Filter data
    #pressure = semi_gauss_lp_filter(dataset.pressure, dataset.sampling_frequency, 4, plot=False)
    #flow = semi_gauss_lp_filter(dataset.flow, dataset.sampling_frequency, 4, plot=False)
    flow = hamming(dataset.flow, 20, 125, 10, plot=True)
    flow = real(flow).tolist()
    pressure = dataset.pressure

    # Get the difference between flow and pressure
    # delay = calc_flow_delay(pressure, flow, dataset.sampling_frequency)

    # Find starts, middles and ends of flow data
    flow_splits = split_breaths(dataset.flow)
    flow_starts = flow_splits[0]
    flow_middles = flow_splits[1]
    flow_stops = flow_splits[2]

    # Set up storage for breath data to correct for pressure delay
    data = BreathData(0)

    # Iterate through each breath in dataset
    for breath in range(2):#len(flow_starts)-1):
        # Put data into breath
        data.get_data(pressure,
                      flow,
                      flow_starts[breath],
                      flow_middles[breath],
                      flow_stops[breath],
                      dataset.sampling_frequency
                      )

        for section in range(2):

            # Set data for insp/exp here
            if(section == INSP):
                breath_length = data.insp_length
                vol = data.insp_volume
                flw = data.insp_flow
            elif(section == EXP):
                breath_length = data.exp_length
                vol = data.exp_volume
                flw = data.exp_flow

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            t = range(450)
            t = [i*0.01 for i in t]

            flw = [3 for i in t]
            vol = integral(flw, 125)
           # vt = [3-exp(-i) for i in t]

            breath_length = len(t)

           # plt.plot(vol)
           # plt.plot(vt)
           # plt.plot(flw)
           # plt.show()
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

            # Storage for PID
            PC = [0]*breath_length
            RC = [0]*breath_length
            res = [0]*breath_length

            # Split breath into little chunks
            chunk_size = 3
            chunk_count = breath_length/chunk_size
            print('Percentage of breath per chunk: {}%'.format(100.0*chunk_size/breath_length))

            # note the data points leftover
            remainder = breath_length - chunk_count*chunk_size

            # Go through all chunks in a breath
            for count in range(chunk_count):
                chunk = Chunk(count, chunk_size, dataset.sampling_frequency)

                # Get results
                chunk.get_PC_and_RC(vol, flw)
                #chunk.try_another_thing(vol)

                # Store results
                for i in range(chunk.size):
                    PC[chunk.start + i] = chunk.PC
                    RC[chunk.start + i] = chunk.RC
                    res[chunk.start + i] = chunk.res

            # Use prev results for remainder
            for i in range(remainder):
                chunk_start = chunk.start + chunk.size
                PC[chunk_start + i] = chunk.PC
                RC[chunk_start + i] = chunk.RC
                res[chunk_start + i] = chunk.res


            if(plot):
                f, (ax2, ax3, ax4, ax5) = plt.subplots(4, sharex=True)
                ax2.set_title('Half Breath Data')

                ax2.plot(flw)
                ax2.set_ylabel('Flow (L/s)')
                ax2.grid()

                fwd_vol = [PC[i] + RC[i]*flw[i] for i in range(breath_length)]
                ax3.plot(vol)
                ax3.plot(fwd_vol)
                ax3.set_ylabel('Volume (L)')
                ax3.legend(['Orig', 'Sim'])
                ax3.grid()

                ax4.plot(PC)
                ax4.set_ylabel(['Pressure/Elastance (kPa/E)'])
                ax4.grid()

                ax5.plot(RC)
                ax5.set_ylabel(['Resistance/Elastance (R/E)'])
                ax5.grid()

                plt.show()

                # Plot saved if plot path given
                if(plot_path != '\0'):
                    fig_name = filename + '_iteration_{}'.format(iteration)
                    print(fig_name)
                    plt.savefig(fig_name)

                plt.close()

    # Final clean up
    dataset.clear_dataset_registry()

    return(PC, RC)

if __name__ == "__main__":
    # Data path
    path = '/home/sarah/Documents/Spirometry/data/'
    files = ['Normal_1.csv', 'loops_1.csv']
    files = [path + name for name in files]

    # Place to save plots
    plot_path = '/home/sarah/Documents/Spirometry/lungPlots/'

    filename = files[0]
    res = analysis(filename, plot=True)

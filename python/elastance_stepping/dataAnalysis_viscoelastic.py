#!/bin/bash

# Add my extensions to path
import sys
sys.path.insert(0, '/home/sarah/Documents/Spirometry/python/extensions')

# Import my extensions
from data_struct import DataStore
from breath_analysis import BreathData, split_breaths, calc_flow_delay
from calculus import integral, derivative
from filters import semi_gauss_lp_filter
from resistance_plots import ResistanceTable

#Import built-ins
import matplotlib.pyplot as plt
from math import sin, cos, pi, exp
from numpy import isnan, array
from numpy.linalg import lstsq

class Chunk:
    def __init__(self,count, size, sampling_frequency):
        self.size = size
        self.volume = [0]*self.size
        self.flow = [0]*self.size
        self.resistance = [0]*self.size
        self.elastance_RQ = 0
        self.elastance_P = 0
        self.start = count*self.size
        self.Fs = sampling_frequency

        self.MIN_V = 0.2 #L

    def get_elastance_and_resistance(self, table, vol, flw, pres):
        # find the volume and flow in the chunk
        self.volume = vol[self.start : (self.start+self.size)]
        self.flow = flw[self.start : (self.start+self.size)]
        self.pressure = [1]*self.size
        if(min(self.volume) > self.MIN_V):

            # Use look-up table to find resistance in the chunk
            for i in range(self.size):
                resistance = table.look_up_resistance(self.flow[i])
                self.resistance[i] = resistance

            # Time varying elastance = integral(resistance*flow)/integral(volume)
            int_V = integral(self.volume, self.Fs)
            int_P =  integral(self.pressure, self.Fs)
            RQ = [(self.resistance[i]*self.flow[i]) for i in range(self.size)]
            int_RQ = integral(RQ, self.Fs)

            # Calculate constant chunk elastance
            self.elastance_RQ = int_RQ[-1]/int_V[-1]
            self.elastance_P = int_P[-1]/int_V[-1]
        else:
            self.elastance_RQ = 0
            self.elastance_P = 0


    def store_results_of_analysis(self, elastance_RQ, res, elastance_P):
        for point in range(self.size):
            elastance_RQ.append(self.elastance_RQ)
            res.append(self.resistance[point])
            elastance_P.append(self.elastance_P)


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
    vol_offset = 2
    pres_offset = 0

    # Remove any old data sets hanging around
    # (Needed when using ipython)
    reset = DataStore()
    reset.clear_dataset_registry()

    # Create data storage
    # Store pressure and flow data from files
    dataset = DataStore(Fs=sampling_frequency)
    dataset.store_data(filename)

    # Filter data
    pressure = semi_gauss_lp_filter(dataset.pressure, dataset.sampling_frequency, 4, plot=False)
    flow = semi_gauss_lp_filter(dataset.flow, dataset.sampling_frequency, 4, plot=False)
    if(plot):
        plt.plot(dataset.pressure, 'g')
        plt.plot(pressure, 'r')
        plt.show()

    # Spirometer resistance
    insp_file = '/home/sarah/Documents/Spirometry/python/extensions/resistance_table_insp.csv'
    exp_file = '/home/sarah/Documents/Spirometry/python/extensions/resistance_table_exp.csv'
    insp_res_table = ResistanceTable(insp_file)
    exp_res_table = ResistanceTable(exp_file)

    # Get the difference between flow and pressure
    delay = calc_flow_delay(pressure, flow, dataset.sampling_frequency)

    # Find starts, middles and ends of flow data
    flow_splits = split_breaths(flow)
    flow_starts = flow_splits[0]
    flow_middles = flow_splits[1]
    flow_stops = flow_splits[2]

    # Set up storage for breath data to correct for pressure delay
    data = BreathData(delay)

    # Iterate through each breath in dataset
    for breath in range(5):#len(flow_starts)-1):
        # Put data into breath
        data.get_data(pressure,
                      flow,
                      flow_starts[breath],
                      flow_middles[breath],
                      flow_stops[breath],
                      dataset.sampling_frequency
                      )

        for section in range(2):
            # Elastance
            elastance_RQ = []
            elastance_P = []
            res = []

            # Set data for insp/exp here
            if(section == INSP):
                breath_length = data.insp_length
                vol = data.insp_volume
                flw = data.insp_flow
                pres = data.insp_pressure
                table = insp_res_table
            elif(section == EXP):
                vol = data.exp_volume
                flw = data.exp_flow
                pres = data.exp_pressure
                breath_length = data.exp_length
                table = exp_res_table

            # Offset volume by estimate tidal breathing residualcapacity
            vol = [vol_offset + v for v in vol]

            # Add atmospheric pressure
            pres = [pres_offset + p for p in pres]

            # Split breath into little chunks
            chunk_size = 5
            chunk_count = breath_length/chunk_size

            # note the data points leftover
            remainder = breath_length - chunk_count*chunk_size

            # Go through all chunks in a breath
            for count in range(chunk_count):
                # Find the start of the chunk
                chunk = Chunk(count, chunk_size, dataset.sampling_frequency)
                chunk.get_elastance_and_resistance(table, vol, flw, pres)
                chunk.store_results_of_analysis(elastance_RQ, res, elastance_P)

            # Add on elastance for remainder
            chunk.size = remainder
            chunk.store_results_of_analysis(elastance_RQ, res, elastance_P)

            if(plot):
                f, (ax3, ax4, ax5, ax6) = plt.subplots(4)
                ax3.set_title('Half Breath Data')

                ax3.plot(pres)
                ax3.set_ylabel('Pressure (kPa)')
                ax3.grid()

                ax4.plot(vol)
                ax4.set_ylabel('Volume (L)')
                ax4.grid()

                ax5.plot(elastance_P)
                ax5.grid()
                ax5.legend(['Elastance (kPa/L)'])
                ax5.set_ylabel('Elastance kPa/L')
                ax5.set_xlabel('Datapoint')

                ax6.plot(flow)
                if(section == INSP):
                    data_range = range(flow_starts[breath], flow_middles[breath])
                elif(section == EXP):
                    data_range = range(flow_middles[breath], flow_stops[breath])
                ax6.plot(data_range, flw, 'r')
                ax6.legend(['Flow (current in red)'])

                plt.show()

                # Plot saved if plot path given
                if(plot_path != '\0'):
                    fig_name = filename + '_iteration_{}'.format(iteration)
                    print(fig_name)
                    plt.savefig(fig_name)

                plt.close()

    # Final clean up
    dataset.clear_dataset_registry()

if __name__ == "__main__":
    # Data path
    path = '/home/sarah/Documents/Spirometry/data/'
    files = ['Loops_1.csv', 'Normal_1.csv', 'Normal_3.csv']
    files = [path + name for name in files]

    # Place to save plots
    plot_path = '/home/sarah/Documents/Spirometry/lungPlots/'

    filename = files[0]
    analysis(filename, plot=True)

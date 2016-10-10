#! /bin/bash

from calculus import integral
from filters import semi_gauss_lp_filter, hamming
import matplotlib.pyplot as plt
from math import pi
from collections import Counter
from numpy import real

class BreathData:
    ''' All the things you could need
        to know about a breath'''
    def __init__(self, flow_delay=28):
        self.FLOW_DELAY = flow_delay
        self.insp_pressure = []
        self.exp_pressure = []
        self.insp_flow = []
        self.exp_flow = []
        self.insp_volume = []
        self.exp_volume = []
        self.insp_length = 0
        self.exp_length = 0
        self.time = []

    def get_data(self, pressure, flow, flow_start, flow_middle, flow_stop, Fs):
        ''' Extract breath data given start and end indices of
            breath and the flow and pressure data'''
        pressure_start = flow_start - self.FLOW_DELAY
        pressure_middle = flow_middle - self.FLOW_DELAY
        pressure_stop = flow_stop - self.FLOW_DELAY

        self.insp_length = pressure_middle - pressure_start
        self.exp_length = pressure_stop - pressure_middle

        self.insp_pressure = pressure[pressure_start:pressure_middle]
        self.exp_pressure = pressure[pressure_middle:pressure_stop]

        self.insp_flow = flow[flow_start:flow_middle]
        self.exp_flow = flow[flow_middle:flow_stop]

        full_flow = self.insp_flow + self.exp_flow
        temp_volume = integral(full_flow, Fs)
        self.insp_volume = temp_volume[:self.insp_length]
        self.exp_volume = temp_volume[self.insp_length:(self.insp_length + self.exp_length)]

        self.time = [t/float(Fs) for t in range(self.insp_length + self.exp_length)]

def split_breaths(data, peak_height=0.1, plot=False):
    ''' Returns array of start and end indices for
        whole breaths in the data set.
        Function doesn't handle data that is noisy
        around zero crossing well. Filter data in
        this case.
        Use data data as input when no peak_height specified'''

    def peak_in_range(data, start, end, peak_height):
        ''' Look for a peak greater than expected noise '''
        found_peak = False

        # Look for a peak
        for i in range(start,end):
            if((abs(data[i]) > peak_height)):
                found_peak = True
        return found_peak

    def check_crossings_valid(data, startpoints, midpoints, endpoints, peak_height, Fs=125):
        ''' Check the breaths are sine-y and are long enough '''
        MIN_RANGE = Fs

        good_start = []
        good_middle = []
        good_end = []
        num_breaths = len(endpoints)

        # Look for peaks in inhalation and exhalation
        # also check the breath is at least 1 second long
        for breath in range(num_breaths):
            if((peak_in_range(data, startpoints[breath], midpoints[breath], peak_height))
            and(peak_in_range(data, midpoints[breath], endpoints[breath], peak_height))
            and(endpoints[breath] - startpoints[breath] > MIN_RANGE)):
                good_start.append(startpoints[breath])
                good_middle.append(midpoints[breath])
                good_end.append(endpoints[breath])
        good_breaths = [good_start, good_middle, good_end]
        return good_breaths

    def find_crossings_and_midpoints(data, startpoints, midpoints, endpoints):
        ''' Has to be for flow data, or else volume will be wrong later'''
        # Set the three areas of the breath to find
        sections = {'start':0, 'middle':1, 'end':2}
        target = sections['start']

        # go through every datapoint
        for current_index in range(1, len(data)-1):
            # Looking for start of breath
            if(target == sections['start']):
                # Find point where zero crossed from below to above
                crossing_found = (data[current_index-1] <= 0 and data[current_index] > 0)
                if(crossing_found):
                    startpoints.append(current_index)
                    target = sections['middle']

            # Looking for midpoint of breath
            elif(target == sections['middle']):
                # find point where zero crossed from above to below
                crossing_found = data[current_index-1] >= 0 and data[current_index] < 0
                if(crossing_found):
                    midpoints.append(current_index)
                    target = sections['end']

            # Looking for end of breath
            elif(target == sections['end']):
                # Find point where zero will cross from below to above
                crossing_found = (data[current_index] <= 0 and data[current_index+1] > 0)
                if(crossing_found):
                    endpoints.append(current_index)
                    target = sections['start']

    # Storage arrays
    startpoints = []
    midpoints = []
    endpoints = []

    # Filter the shit out of the signal
    data = hamming(data, 5, 125, 10, plot)
    data = real(data).tolist()
    # Find the crossing points and check they are good
    find_crossings_and_midpoints(data, startpoints, midpoints, endpoints)
    final_breaths = check_crossings_valid(data, startpoints, midpoints, endpoints, peak_height)

    if(plot):
        fss = [data[o] for o in startpoints]
        print(startpoints)
        fms = [data[o] for o in midpoints]
        fes = [data[o] for o in endpoints]
        fts = [t for t in startpoints]
        ftm = [t for t in midpoints]
        fte = [t for t in endpoints]
        plt.plot(range(len(data)), data, 'b',
                fts, fss, 'og',
                ftm, fms, 'oy',
                fte, fes, 'or',
                )
        plt.legend(['data', 'start', 'middle', 'end'])
        plt.grid()
        plt.show()

    return final_breaths

def calc_flow_delay(pressure, flow, Fs=125, plot=False):
    # Filter the data heaps and flip pressure
    flow = semi_gauss_lp_filter(flow, 125, 0.5)
    flow[0] = 0
    pressure = semi_gauss_lp_filter(pressure, 125, 0.5)
    pressure = [-p for p in pressure]
    pressure[0] = 0

    # Get the crossing points
    splits = split_breaths(flow)
    p_splits = split_breaths(pressure, 0.01)

    # Compare position of centre crossings
    # If result > sampling frequency,
    # a breath has been skipped so don't record
    length = min(len(p_splits[1]), len(splits[1]))
    diff  = []
    p_offset = 0
    f_offset = 0
    i = 0
    while (i + max(p_offset, f_offset)) < length:
        difference = (splits[1][i + f_offset] - p_splits[1][i + p_offset])
        if(plot):
            print(difference)
        if(abs(difference) < Fs):
            diff.append(difference)
        elif(difference > Fs):
            # Skip pressure forward to catch up
            if(plot):
                print('missed flow')
            p_offset += 1
        else:
            # Skip flow forward to catch up
            if(plot):
                print('missed pressure')
            f_offset += 1
        i += 1

    def mean(item_list):
        num_items = (len(item_list))
        mean = 0
        for i in range(num_items):
            mean += item_list[i]/float(num_items)
        return int(round(mean))

    # Take mean as the shift
    #counter = Counter(diff)
    #max_count = max(counter.values())
    #mode = [k for k,v in counter.items() if v == max_count]
    #shift = mean(mode)
    shift = mean(diff)

    if plot:
        print('{}'.format(diff))
        print(shift)
        fss = [flow[o] for o in splits[0]]
        fms = [flow[o] for o in splits[1]]
        fes = [flow[o] for o in splits[2]]
        fts = [t for t in splits[0]]
        ftm = [t for t in splits[1]]
        fte = [t for t in splits[2]]
        ss = [pressure[o] for o in p_splits[0]]
        ms = [pressure[o] for o in p_splits[1]]
        es = [pressure[o] for o in p_splits[2]]
        ts = [t for t in p_splits[0]]
        tm = [t for t in p_splits[1]]
        te = [t for t in p_splits[2]]
        plt.plot(range(len(pressure)), pressure, 'b',
                range(len(flow)), flow, 'r',
                ts, ss, 'sg',
                tm, ms, 'sy',
                fts, fss, 'og',
                ftm, fms, 'oy',
                fte, fes, 'or',
                te, es, 'sr')
        plt.grid()
        plt.show()

    return shift

#!/bin/bash

from scipy.integrate import cumtrapz

def derivative(data, Fs):
    '''Calc derivative of list using difference'''
    dir_list = [0]*(len(data)-1)
    for i in range(1, len(data)):
        dir_list[i-1] = (data[i] - data[i-1])*Fs
    return dir_list

def integral(data, Fs):
    '''Calc integral using cumtrapz'''
    int_data = cumtrapz(data).tolist()
    int_data = [0.0] + [int_data[i]/float(Fs) for i in range(len(int_data))]
    return int_data

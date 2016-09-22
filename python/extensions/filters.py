#!/bin/bash

import numpy
import numpy as np
import time
from numpy.fft import ifft, fft, fftshift
from math import exp, sqrt, pi, cos, sin
import matplotlib.pyplot as plt

def gauss_lp_filter(data, Fs, sigma, plot=False):
    ''' Uses gauss function as low pass filter
        Inputs: data, sampling frequency, sigma value for gauss function'''
    # Fourier transform of data
    fft_data = fft(data)
    datalength = len(data)

    # Remove mirrored half
    freq_data = fft_data[0:(datalength+1)/2]

    # Normalise data
    freq_data = [freq_data[0]] + [f*2 for f in freq_data[1:]]

    # Gauss function centred at x=0
    data_max = max(freq_data)
    print(abs(data_max))
    CENTRE = freq_data.index(data_max)/float(Fs)
    print(CENTRE)

    # Get x values of frequency
    x = [v/float(Fs) for v in range(datalength)]

    # Make gauss function
    factor = (1/(sigma*sqrt(2*pi)))
    power = [((-(x[i]-CENTRE)**2)/(2*sigma**2)) for i in range(datalength)]
    fx = [factor*exp(power[i]) for i in range(datalength)]

    # Maximum value in gauss function is set to 1
    scaling = max(fx)
    fx = [v/scaling for v in fx]

    if(plot):
        fd = [abs(d) for d in freq_data]
        plt.plot(fd)
        gauss_line = [f*abs(data_max) for f in fx]
        plt.plot(gauss_line)
        plt.legend(['Freq domain', 'gauss filter'])
        plt.grid()
        plt.show()

    # apply filter to data
    filtering_data = freq_data + [0]*len(freq_data)
    filtered_data = [filtering_data[i]*fx[i] for i in range(datalength)]
    filtered_data = ifft(filtered_data)

    # Only grab the real part if imaginary components left over
    try:
        filtered_data = filtered_data.real.tolist()
    except AttributeError:
        pass

    return filtered_data

def semi_gauss_lp_filter(data, Fs, sigma, plot=False):
    ''' Uses gauss function as low pass filter
        Inputs: data, sampling frequency, sigma value for gauss function'''
    # Fourier transform of data
    fft_data = fft(data)
    datalength = len(data)

    # Remove mirrored half
    freq_data = fft_data[0:(datalength+1)/2]

    # Normalise data
    freq_data = [freq_data[0]] + [f*2 for f in freq_data[1:]]

    # Gauss function centred at x=0
    data_max = max(freq_data)
    max_location = freq_data.index(data_max)
    CENTRE = max_location/float(Fs)

    # Get x values of frequency
    x = [v/float(Fs) for v in range(datalength)]

    # Make gauss function
    factor = (1/(sigma*sqrt(2*pi)))
    power = [((-(x[i]-CENTRE)**2)/(2*sigma**2)) for i in range(datalength)]
    fx = [factor*exp(power[i]) for i in range(datalength)]

    # Maximum value in gauss function is set to 1
    scaling = max(fx)
    fx = [v/scaling for v in fx]

    # Set all up to the centre point to 1
    for i in range(max_location):
        fx[i] = 1

    if(plot):
        fd = [abs(d) for d in freq_data]
        plt.plot(fd)
        gauss_line = [f*abs(data_max) for f in fx]
        plt.plot(gauss_line)
        plt.legend(['Freq domain', 'gauss filter'])
        plt.grid()
        plt.show()

    # apply filter to data
    filtering_data = freq_data + [0]*len(freq_data)
    filtered_data = [filtering_data[i]*fx[i] for i in range(datalength)]
    filtered_data = ifft(filtered_data)

    # Only grab the real part if imaginary components left over
    try:
        filtered_data = filtered_data.real.tolist()
    except AttributeError:
        pass

def hamming(data, fc, fs, bw, plot=False):
    # Length of the filter.
    data_length = len(data)

    fc = fc/float(fs)/2.0
    tw = bw/float(fs)
    N = int(np.ceil((4/tw)))

    if not N%2: N+=1
    #n = np.arrange(N)

    # Calculates points in time domain hamming window
    def hamming_value(n, N):
        value = 0.54 - 0.46*cos(2*pi*n / float(N-1))
        return(value)

    hamming = [0]*N
    trunc_sinc =[0]*N
    for n in range(N):
        hamming[n] = hamming_value(n, N)
        if(n - (N-1)/2.0 == 0):
            trunc_sinc[n] = 2*np.pi*fc
        else:
            i = (n - (N-1)/2.0)
            trunc_sinc[n] = sin(2 * np.pi * fc * i)/(i)

    window = [trunc_sinc[i] * hamming[i] for i in range(N)]
    window = window/np.sum(window)


    fft_window = fft(window, data_length)
    fft_signal = fft(data)
    fft_filtered = [fft_window[i] * fft_signal[i] for i in range(data_length)]
    filtered_signal= ifft(fft_filtered)

    samples_shift = int((N-1)/2)
    filtered_signal = filtered_signal[samples_shift:len(filtered_signal)]

    if(plot):
        plt.plot(window, 'r')
        plt.show()

        plt.plot(np.abs(fft(window)))
        plt.show()

        plt.plot(freq, mag, 'm')
        plt.show()

        plt.plot(freq, response, 'r')
        plt.show()

    return(filtered_signal)

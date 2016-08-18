#!/bin/bash

import numpy
from numpy.fft import ifft, fft
from math import exp, sqrt, pi
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

    return filtered_data

def hamming_low_pass_filter(fc, transition_width, signal, dt, time_shift_removed='n'):
    """
    Created on Wed Jun 22 12:40:00 2016

    @author: jcb137

    ---____ Creating a Hamming Filter_________________________________

    INPUTS:
    * fc is the cutoff frequency
    * transition_width is the width of the transition band in Hz.
        NB it is then normalized by dividing it by the sampling frequency
    * signal is the time domain signal you want to filter once the filter is created.
    * dt is the step size between successive samples in the time domain for the signal.
    * time_shift_removed is a yes 'y', no 'n' input which the user enter to indicate whether
        they want the output filtered signal to have been corrected so it has zero time shift

    OUTPUTS:
    filtered_signal is the filtered signal
        NB! Filtering results in Nf unusable points at the start of the signal in the time domain.
        This is becauase we need Nf previous points to filter
        (see Nf output below) any given data point, hence the first Nf points in a signal can not
        actually be filtered. So unusable points are set to zero and the filtered signal output
        is kept the same length as the input signal. The exception to this when the output has
        its time_shift_removed,(see input above) this leads an output which is
        (Nf-1)/2 = samples_shift shorter.

    * h is the filter impulse responce in the time domain
    * Nf is the filter length ie the number of filter coefficients
    * freq_spec_filt is an array of the discrete frequencies for which filter coefficients
        exist, making up the filter frequency spectrum (useful for plotting)
    * t_array_filt is an array of the discrete times for which filter impulse responce is known.
        (useful for plotting)
    * t_shift is the time shift or delay in seconds caused by the filter on the output when
        compared to the input
    * samples_shift is the shift or delay in SAMPLES caused by the filter on the output when
        compared to the input
    * w is the hamming window function in the time domain
    * trunc_sinc is our truncated sinc function aka our truncated ideal filter responce in the
        time domamin
    """

    # Generating the low pass filter impulse responce_____________
    # sampling frequency of the input signal
    fs = 1/dt

    # number of points in signal
    N = len(signal)

    # normalized cutoff frequency
    fc_norm = fc/fs

    # normalized transition width
    tw_norm = transition_width/fs

    # working out the number of filter coefficients (filter length)
    # required for the hamming window.
    Nf = int(numpy.ceil(3.3/tw_norm))

    # if Nf is worked out to be even, we add one to make it odd.
    # This isnt always necessary, but even and odd filter impulse responces are
    # slightly different and the odd type is more versatile see Ifeachor Emmanuel,
    # Jervis Barrie - 1993 - Digital Signal Processing Chpt on FIR filter design.
    # (pg 284 in 1st Edition)
    if Nf%2 == 0:
        Nf = Nf+1

    # Working out the Delay caused by the filter
     # time shift, see Ifeachor Emmanuel, Jervis Barrie - 1993 - Digital Signal
     # Processing Chpt on FIR filter design.
    t_shift = ((Nf-1)/2)*dt

    # working out the shift as the equivalent of sample rather than seconds.
    samples_shift = int((Nf-1)/2)

    # Initilizing arrays
    # y_filtered will be the output signal from the filter
    filtered_signal = numpy.zeros(N)

    # h will be the filter impulse responce
    h = numpy.zeros(Nf)

    # w will be the filter window
    w = numpy.zeros(Nf)

    # trunc_sinc is our ideal low pass filter impulse responce truncated to Nf points
    trunc_sinc = numpy.zeros(Nf)

    # working out frequency spectrum and time arrays
    df_filt = 1/(Nf*dt)
    t_array_filt = numpy.arange(0, Nf*dt, dt)[:Nf]
    freq_spec_filt = numpy.arange(0, Nf*df_filt, df_filt)[:Nf]

    for n in range(0,Nf):
        w[n] = (25/46) - (21/46)*numpy.cos(2*numpy.pi*n/(Nf-1))
        if n-(Nf-1)/2 == 0:
            trunc_sinc[n] = 2*numpy.pi*fc_norm
            h[n] = w[n]*trunc_sinc[n]
        else:
            trunc_sinc[n] = numpy.sin(2*numpy.pi*fc_norm*(n-(Nf-1)/2))/(n-(Nf-1)/2)
            h[n] = w[n]*trunc_sinc[n]

    # before I normalize h (necessary for using filter) I want to keep h as it was calculated,
    # in case the person wants to actually look at the calculated h=w*trunc_sinc,
    # could be useful for debugging
    h_unnorm = numpy.copy(h)

    # normalize the low pass filter impulse responce
    h = h/numpy.sum(h)

    # Performing the Convolution

    # NB with a stopping condition of N in the range function,
    # the final computed itereation will be j= N-1
    for j in range(Nf-1, N):
        filtered_signal[j] = 0
     # NB with a stopping condition of Nf in the range function,
     # the final computed iteration will be n= Nf-1
        for n in range(0,Nf):
            filtered_signal[j] = filtered_signal[j] + (signal[j-n]*h[n])

    # Accounting for the time shift caused by filtering______________

    # The user can choose to remove the time shift caused by filtering. If they define
    # time_shift_removed = 'y' the output filtered_signal
    # will be have zero phase delay/time shift. The output signal is shifted to have zero
    # time shift and as a result is 'samples_shift' shorter.

    def yes_no_input_checker(user_input):
        ''' used to check the user input is either 'y' for yes or 'n' for no. '''
        if user_input != 'y' and user_input != 'n':
            user_input = input('Invalid input, please re-enter y for yes or n for no: ')
            return(yes_no_input_checker(user_input))
        else:
            return(user_input)

    time_shift_removed = yes_no_input_checker(time_shift_removed)

    if time_shift_removed == 'y':
        # shifting the filtered_signal 'samples_shift' left so it now has zero time delay.
        filtered_signal = filtered_signal[samples_shift:len(filtered_signal)]

    return(filtered_signal,
           h,
           Nf,
           freq_spec_filt,
           t_array_filt, t_shift,
           samples_shift,
           w,
           trunc_sinc)

# -*- coding: utf-8 -*-
"""
Created on Tue Jun 28 11:06:02 2016

@author: jcb137
"""

# ---_______________________________________ TEST FILTER 2 __________________________________________________

# Creating filter based on 'Digital Signal Processing: a practical guide for engineers and scientists' By Steven 2003 pg 294

# Low pass Hamming Filter

# NB NB NB NB NB NB
# One of the most important things I recognised from this example is:
# 1) You MUST use the NORMALISED cut off frequency in calculation your filter impulse responce.
# 2) fc_norm = fc/fs, so the higher your sampling frequency the smaller fc_norm will be. This can 
#    lead to issues. For example:
#
#           Say you want to generate y = sin(2*pi*0.5*t) + sin(2*pi*5*t) so your highest frequency is 5Hz
#           
#           to make the signal smooth you use fs = 100*highest_freq_interest = 100*5 = 500Hz.
#           now say you want to filter out the 5Hz component from your generated signal. Let
#           fc = 3.5Hz. In order to not affect our 0.5Hz component and fully remove the 5Hz using fc=3.5Hz
#           we would need a transition width of tw <= 3Hz, meaning the f_stop <= 5Hz (where f_stop = 3.5Hz + 3Hz/2) 
#           while f_pass <= 2Hz (whre f_pass = 3.5Hz - 3Hz/2).
#
#           The next step is to take our filter characteristics and normalize them to find the number of filter
#           coefficients we need.  fc_norm = 3.5Hz/500Hz = 0.007 (eq from above). Now this is where we start to 
#           see the issue with using very high fs, our normalized 0.5Hz component is 0.001 (=0.5Hz/500Hz), which 
#           now appears VERY close to our fc_norm of 0.007 when inspecting these in the normalized frequency
#           spectrum. This means normalized transition band tw_norm = 0.006 (= 3Hz/500Hz) is VERY VERY NARROW!
#           It is narrow because the roll off of the filter reponce needs to be between 0.001 (0.5Hz) and 0.01
#           (5Hz). In order to get such incredibly fast roll off we will need a large number of filter coefficients
#           for example the hamming filter would need around Nf = 3.3/tw_norm = 550. Now this isnt huge by absolute
#           standards, BUT when you consider 10sec of our signal is 5000 points, we require 11% of our signal. 
#
#
#           What if we wanted to fix the filter length:
#           Now lets say we wanted to only use a filter length of Nf = 101, our tw_norm would need to be approx
#           0.033 --> tw = tw_norm*fs = 16.5Hz! So our transition band is waaaaaaayyyyyy to big! Both the frequency
#           we want to filter (5Hz) and the one we want to preserve (0.5Hz) will now be in the transition band!
#
#           To fix this we could change our fs to 25Hz, still 5 times higher than the highest frequency in our
#           original signal y. Then our new fc_norm will be fc_norm = 3.5Hz/25Hz = 0.14, new tw_norm = 3Hz/25Hz
#           = 0.12, 20 times wider than when tw_norm = 0.006! So now we can get away with far less filter coefficients
#           and our filtered signal frequency of 0.5Hz is still 50 times smaller than fs so our filtered signal
#           will still look good as well!

# Importing modules
import numpy
import matplotlib.pyplot as plt

# Initializing Arrays
N = 5000    # defining signal length
Nf = 101    # defining Filter length
y = numpy.zeros(N)   # y will be the input signal to the filter
y_filtered = numpy.zeros(N)   # y_filtered will be the output signal from the filter
h = numpy.zeros(Nf)    # h will be the filter impulse responce
w = numpy.zeros(Nf)    # w will be the filter window
trunc_sinc = numpy.zeros(Nf)    # trunc_sinc is our ideal low pass filter impulse responce truncated to Nf points

# Filter characteristics
fc_norm = 0.14         # normalized cutoff frequency
tw_norm = 3.3/Nf       # normalized tranisiton width tw_norm = 5.5/101 based on Ifeachor Signal Processing book

# generating input signal _________________________
fa = 0.5    # low frequency component of signal input y
fb = 5      # high frequency component of signal input y

# Working out sampling frequency, frequency spectrum and time arrays
fs = 5*(max(fa,fb))     # fs = 25, see 'Effects of Sampling Frequency on Signal Processing' in Onedrive python folder

dt = 1/fs       # 0.02sec per step
df = 1/(N*dt)   # 0.01Hz per step
t = numpy.arange(0, N*dt, dt)  # time array starting at t[0] = 0sec to t[4999] = 9.998sec
freq_spec = numpy.arange(0, N*df, df)

for i in range(0,N):
    y[i] = numpy.sin(2*numpy.pi*fa*t[i]) + numpy.sin(2*numpy.pi*fb*t[i])
    
# Generating the low pass filter impulse responce_____________
    
    # working out frequency spectrum and time arrays
    df_filt = 1/(Nf*dt)
    t_filter = numpy.arange(0, Nf*dt, dt)
    freq_spec_filt = numpy.arange(0, Nf*df_filt, df_filt)
    
for n in range(0,Nf):
    w[n] = (25/46) - (21/46)*numpy.cos(2*numpy.pi*n/(Nf-1))
    if n-(Nf-1)/2 == 0:
        trunc_sinc[n] = 2*numpy.pi*fc_norm
        h[n] = w[n]*trunc_sinc[n]
    else:
        trunc_sinc[n] = numpy.sin(2*numpy.pi*fc_norm*(n-(Nf-1)/2))/(n-(Nf-1)/2)
        h[n] = w[n]*trunc_sinc[n]

h_before_normalisation = numpy.copy(h)  # before I normalize h (necessary for using filter) I want to keep h as it was calculated,
                                        # this is because later I plot h, along with trunc_sinc and w.
        
# normalize the low pass filter impulse responce
h = h/numpy.sum(h)

# Performing the Convolution
for j in range(Nf-1, N):      # NB with a stopping condition of N in the range function, the final computed itereation will be j= N-1
    y_filtered[j] = 0
    for n in range(0,Nf):     # NB with a stopping condition of Nf in the range function, the final computed itereation will be n= Nf-1
        y_filtered[j] = y_filtered[j] + (y[j-n]*h[n])
        
#%% Plotting Results
Y = numpy.fft.fft(y)     
H = numpy.fft.fft(h)
Y_filtered = numpy.fft.fft(y_filtered)

fig, axs = plt.subplots(nrows= 1, ncols = 2)
fig.suptitle('Filtering the 5Hz component out of: \n $y = sin(2 \pi t 0.5Hz) + sin(2 \pi t 5Hz)$', fontsize = 20)

# Plotting Time Domain Signals
axs[0].set_title('Signal in the Time Domain')
axs[0].set_xlabel('t (sec)')
axs[0].set_ylabel('y')
axs[0].plot(t, y, color = 'b', label = 'Original Signal')
axs[0].plot(t, y_filtered, color = 'r', label = 'Filtered Signal')
axs[0].legend(fontsize = 'small')

# Plotting Frequency Responce
axs[1].set_title('Normalised Signals in Frequency Domain')
axs[1].set_xlabel('f (Hz)')
axs[1].set_ylabel('Magnitude')
axs[1].stem(freq_spec,abs(Y)/max(abs(Y)), linefmt = 'b', basefmt = 'b', markerfmt = ' ', label = 'Orignal signal')
axs[1].stem(freq_spec,abs(Y_filtered)/max(abs(Y_filtered)),linefmt = 'r', basefmt = 'r', markerfmt = ' ', label = 'Filtered Signal' )
axs[1].plot(freq_spec_filt, abs(H)/max(abs(H)), color = 'g', label = 'Filter Freq Responce' )
axs[1].legend(fontsize = 'small')

# Plotting Hamming Window and Truncated Ideal Filter
plt.figure()
plt.title('Filter Impulse Responce and its Components')
plt.xlabel('t (sec)')
plt.ylabel('y')
plt.plot(numpy.arange(0, Nf*dt, dt) ,trunc_sinc, color = 'orange', label = 'Truncated Ideal Filter Impulse Responce')
plt.plot(numpy.arange(0, Nf*dt, dt), w, color = 'purple', label = 'Hamming Window')
plt.plot(numpy.arange(0, Nf*dt, dt), h_before_normalisation, color = 'g', label = 'Filter Impulse Responce (h = ideal x window)')
plt.legend(fontsize = 'small')

#%% Working out the Bode Plot for the filter

# Working out the magnitude plot
fig2, axs2 = plt.subplots(nrows = 2, ncols = 1, sharex = True)
axs2[0].set_title('Magnitude Responce of the Filter')
axs2[1].set_title('Phase Responce of the Filter')
axs2[1].set_xlabel('Frequency $f$ (Hz)')
axs2[0].set_ylabel('Mag |H|')
axs2[1].set_ylabel(r'Phase $\theta$ (deg)')

axs2[0].plot(freq_spec_filt[0:int(numpy.floor(len(freq_spec_filt)/2))], abs(H)[0:int(numpy.floor(len(freq_spec_filt)/2))])     # Plotting the magnitude responce of the filter

Phase_H = numpy.arctan(numpy.imag(H)/numpy.real(H))*(180/numpy.pi) # working out the phase responce of the filter
axs2[1].plot(freq_spec_filt[0:int(numpy.floor(len(freq_spec_filt)/2))], Phase_H[0:int(numpy.floor(len(freq_spec_filt)/2))], label = r'Phase $\theta = 7.2f\rightarrow$ time delay of $t_d=\frac{\theta}{f\times(360\degree/cycle)}=0.02sec$')
axs2[1].plot([0,0.5,0.5],[3.6,3.6,0],linewidth = 2, color = 'r', label = r'$\theta(f=0.5Hz)=3.6\degree$')
axs2[1].legend()          
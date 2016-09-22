from numpy import sqrt, exp, pi
from numpy import fft
from math import atan, cos, sin
import csv
import matplotlib.pyplot as plt

T = 4 #s
C = 2
R = 1.38
Rsp = 0.061
Fs = 125

w = 2*pi/T
print('w: {}'.format(w))
time = [x/float(Fs) for x in range(0,T*2*Fs)]
wRC = w*R*C

P_in = [2*cos(w*t) for t in time]
plt.plot(P_in, 'k')

Gc = 1/(sqrt(1+(wRC)**2))
theta = atan(-wRC)
print('C theta: {}'.format(theta))
amp = Gc * max(P_in)
wave = [cos(w*t + theta) for t in time]
Vc = [w*amp for w in wave]
plt.plot(Vc, 'r')

Vr = [P_in[i] - Vc[i] for i in range(len(P_in))]
theta_R = atan(1/wRC)
print('R theta: {}'.format(theta_R))
plt.plot(Vr, 'g')

Vsp = [(Rsp/R)*v for v in Vr]
plt.plot(Vsp, 'm')

I = [v/R for v in Vr]
plt.plot(I, 'y')

plt.legend(['P_in', 'P_C', 'P_R', 'P_Rsp', 'I'])
plt.grid()
plt.show()

saving = False
if(saving):
    with open('generated_breath.csv', 'w') as csvfile:
        csvwriter = csv.writer(csvfile)
        for i in range(len(Vsp)):
            line = [Vsp[i], I[i]]
            csvwriter.writerow(line)



#!/bin/bash

# Import my extensions
import sys
sys.path.insert(0, '/home/sarah/Documents/Spirometry/python/extensions')

import csv
import matplotlib.pyplot as plt
from data_struct import DataStore as Data
from calculus import integral
import numpy as np

path = '/home/sarah/Documents/Spirometry/python/office_trial/trial_data/'
subject = [
          'F01',
          'F02',
          'F03',
          'F04',
          'F06',
          'F09',
          'M01',
          'M02',
          'M03',
          'M04',
          'M05',
          'M06',
          'M07',
          'M08',
          'M09',
          'M10',

          ]
test_type = [
        '_mp_noRx.csv',
        '_mp_Rx.csv',
        ]

# Coefficients for Ps -> Pm
# Pm = Ps*c[X][0] + c[X][1]
# X depends on pressure range [low, high]
# c = [crossing_point, [mp_RX, mp_woRX]]
conversion = [
        [-6, 1.5],
        [-12, 1.3],
        #[0.027, [1.1, 0], [2.9, 0.05]],
        #[0.023, [1.1, 0], [2.5, 0.03]],
        ]



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Extract pressure two data arrays (assumed pressure and flow,
# but could be pressure & pressure etc) and time stamp from csv
# NOTE This function appends. Make sure empty arrays are passed in
def extract_data_arrays_from_file(filename, pressure, flow, time):
    # Store pressure and flow data
    with open(filename, 'rb') as csvfile:
        reader = csv.reader(csvfile, delimiter = ',')

        # skip header line
        header = reader.next()

        for row in reader:
            # Flipping data because looking at inspiration and calibration
            # also flipped the data.
            pressure_point = -float(row[0])
            flow_point = -float(row[1])
            time_point = float(row[2])

            pressure.append(pressure_point)
            flow.append(flow_point)
            time.append(time_point)

    return(1)
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Calibration for flow sensor is in
# ~/Documents/Spirometry/python/calibration
# This may need to be re-checked
# NOTE Only valid up to 250Pa across spirometer
def convert_spirometer_pressure_to_flow(pressure):
    flow = [(-0.000475*p**2) + (0.427*p) for p in pressure]
    return flow
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Calibration for this conversion is in
# ~/Documents/Spirometry/python/maskLPF/Ps_Pm_calibration
# This conversion rate only applies to expiration
def convert_spirometer_pressure_to_mouth_pressure(spir_pressure, coeffs):
    mouth_pressure = [coeffs[0]*p**2 + coeffs[1]*p for p in spir_pressure]
    return mouth_pressure
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Resistance from mouth onwards is calculated
def calculate_external_resistance(mouth_pressure, flow):
    resistance = []*len(flow)
    for m in range(len(mouth_pressure)):
        if mouth_pressure[m] < -0.1 and mouth_pressure[m] > -0.15:
            try:
                R = mouth_pressure[m]/flow[m]
                if R > 0:
                    resistance.append( R )
            except ZeroDivisionError:
                R = np.nan
        else:
            R = np.nan
        #resistance[m] = R

    return resistance
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


for j in range(len(subject)):
    resistances = [[], []]
    flows = [[], []]
    pressures = [[], []]

    for k in range(2):
        filename = path + subject[j] + '/normalised_' + subject[j] + test_type[k]

        try:
            # Extract data to empty arrays.
            # This will fail if file does not exist
            pressure_lowRes = []
            pressure_highRes = []
            time = []
            extract_data_arrays_from_file(filename, pressure_lowRes, pressure_highRes, time)

            flow = convert_spirometer_pressure_to_flow(pressure_highRes)
            mouth_pressure = convert_spirometer_pressure_to_mouth_pressure(pressure_highRes, conversion[k])

            flows[k] = flow
            resistances[k] = calculate_external_resistance(mouth_pressure, flow)
            pressures[k] = mouth_pressure

            plt.plot(pressure_highRes)
            plt.plot(mouth_pressure)
            plt.plot(flow)
            plt.plot(resistances[k])
            plt.legend(['P_s', 'P_m', 'Q', 'R'])
            plt.show()

        except IOError:
            pass


    print('res 1')
    print(max(resistances[0]))
    print(min(resistances[0]))
    print(np.mean(resistances[0]))
    print(np.median(resistances[0]))

    print('res 2')
    print(max(resistances[1]))
    print(min(resistances[1]))
    print(np.mean(resistances[1]))
    print(np.median(resistances[1]))

    print('difference')
    print(np.mean(resistances[1]) - np.mean(resistances[0]))
    print(np.median(resistances[1]) - np.median(resistances[0]))

    plt.plot(pressures[0], flows[0], 'rx')
    plt.plot(pressures[1], flows[1], 'bo')
    plt.xlabel('pressure')
    plt.ylabel('flow')
    plt.legend(['woRx', 'Rx'])
    plt.show()

    resistances[0].sort()
    resistances[1].sort()
    plt.plot(resistances[0], 'x')
    plt.plot(resistances[1], 'o')
    plt.show()


    Ps = range(-50, -70, -1)
    Ps = [p/1000.0 for p in Ps]

    Q = convert_spirometer_pressure_to_flow(Ps)

    Pm_Rx = [2.7*p + 0.04 for p in Ps]
    Pm_woRx = [2.4*p + 0.03 for p in Ps]
    PmQ_Rx = [-12*p**2 + 1.3*p for p in Ps]
    PmQ_woRx = [-6.5*p**2 + 1.5*p for p in Ps]

    delta_Pm = [0.40*p + 0.00 for p in Ps]
    delta_rm = [Pm_Rx[i] - Pm_woRx[i] for i in range(len(Q))]
    delta_rmQ = [PmQ_Rx[i] - PmQ_woRx[i] for i in range(len(Q))]

    Rx = [0]*len(Q)
    Rm_Rx = [0]*len(Q)
    RmQ_Rx = [0]*len(Q)
    for index in range(len(Q)):
        Rx[index] = delta_Pm[index]/Q[index]
        Rm_Rx[index] = delta_rm[index]/Q[index]
        RmQ_Rx[index] = delta_rmQ[index]/Q[index]


    mean_Rx = np.mean(Rx)
    print(mean_Rx)
    print(mean_Rx/6.0)
   # mean_Rx = [mean_Rx]*len(Q)

    plt.plot(Ps, 'b')
    plt.plot(Pm_Rx, 'r')
    plt.plot(Pm_woRx, 'g')
    plt.plot(PmQ_Rx, 'y')
    plt.plot(PmQ_woRx, 'm')
    plt.plot(delta_Pm, 'k')
    plt.plot(delta_rm, 'k:')
    plt.plot(delta_rmQ, 'k.')

    plt.plot(Q, 'orange')

    plt.plot(Rx, 'c.-')
    plt.plot(Rm_Rx, 'r.-')
    plt.plot(RmQ_Rx, 'y.-')
    plt.legend(['Ps', 'pm_Rx', 'pm_woRx', 'pmQ_Rx', 'pmQ_woRx', 'delta_pm',
                'Q',
                'Rx', 'rm_Rx', 'rmQ_Rx',
                'delta rm', 'delta rmQ', 'mean_Rx'])
    plt.grid()
    plt.show()


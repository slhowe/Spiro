#!/bin/bash

# Import my extensions
import sys
sys.path.insert(0, '/home/sarah/Documents/Spirometry/python/extensions')

import csv
import matplotlib.pyplot as plt
from data_struct import DataStore as Data
from calculus import integral
import numpy as np

path = './resistance_data/'
subject = [
          'normalised_resistance_data_short_tube_flow_reference.csv',
          'normalised_resistance_data_long_tube_flow_reference.csv',
          'normalised_resistance_data_peak_long_tube_flow_reference.csv',
          ]

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# NOTE
# Three data sets per file
# First is pressure after resistance
# Second is pressure before resistance
# Third is time
#
# All pressures are relative to end of second (after res) pneumotach.
# Flow calculated form second pneumotach pressure
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Extract pressure two data arrays (assumed pressure and flow,
# but could be pressure & pressure etc) stamp from csv and sort
# NOTE This function appends. Make sure empty arrays are passed in
def extract_data_arrays_from_file(filename, pressure, flow):
    data = []
    # Store pressure and flow data
    with open(filename, 'rb') as csvfile:
        reader = csv.reader(csvfile, delimiter = ',')

        # skip header line
        header = reader.next()

        for row in reader:
            # convert pressure to Pa
            pressure_point = float(row[0])*1000
            flow_point = float(row[1])*1000
            # Pressures linked in tuples for later sorting
            data.append((pressure_point, flow_point))

        # Sorted in ascending order by mask pressure
        # Because mask pressure will hit pressure cutoff first
        data.sort(key=lambda tup: tup[1])
        for i in range(len(data)):
            pressure.append(data[i][0])
            flow.append(data[i][1])


    return(1)
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Calibration for flow sensor is in
# ~/Documents/Spirometry/python/calibration
# This may need to be re-checked
# Returns flow in L/s
# NOTE Only valid up to 250Pa across spirometer
def convert_spirometer_pressure_to_flow(pressure):
    flow = [np.nan]*len(pressure)
    for m in range(len(pressure)):
        if pressure[m] > 0:
            flow[m] = (-0.000475*(pressure[m])**2 + 0.427*(pressure[m]))
        else:
            p = abs(pressure[m])
            flow[m] = -((-0.000475*(p)**2) + (0.427*(p)))

    flow = [f/60.0 for f in flow]
    return flow
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Resistance from mouth onwards is calculated
def calculate_external_resistance(pressure, flow):
    resistance = [np.nan]*len(flow)
    for m in range(len(pressure)):
        #if(abs(flow[m])>0.83 and abs(flow[m])<1.17):
        if(abs(flow[m])>0.3):
            try:
                R = pressure[m]/flow[m]
            except ZeroDivisionError:
                R = np.nan
            resistance[m] = R
    return resistance
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
def calculate_constant_resistance_sections(resistance):
    constant_resistance = [np.nan]*len(resistance)
    g = (i for i, e in enumerate(delta_pressure) if e >= 0.2)
    start = next(g)
    step = 1000
    stop = len(resistance) - step + 1
    current = start

    while current < stop:
        res = resistance[current:current+step]
        constant_resistance[current:current+step] = [np.mean(res)]*step
        current += step

    res = resistance[current:len(resistance)]
    constant_resistance[current:len(resistance)] = [np.mean(res)]*(len(resistance)-current)

    return constant_resistance
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Discharge coefficient for orifice plate
# mass_flowrate = ( (A_2*C)*(sqrt(2*roh*delta_P)) )/( sqrt(1-(A_2/A_1)^2) )
# for A1 = 15mm dia
#     A2 = 10 mm dia
#     roh - 1.225 kg/m3
#     deltaP in Pa
#
#     den = (sqrt(1-(A_2/A_1)^2))
#     num = A_2*(sqrt(2*roh*delta_P))
def calculate_discharge_coefficient(flow, delta_pressure):
    coeff = [np.nan]*len(flow)
    for m in range(len(coeff)):
        if(flow[m]>0.0 and flow[m]<0.8):
            den = 0.896
            num = 7.85e-5 * np.sqrt(2.45*(delta_pressure[m])*100)
            fact = num/den
            coeff[m] = (flow[m]*1e-3)/fact
    return coeff
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
def calculate_friction_coefficient():
    beta = (0.50*7.85e-5)/(np.sqrt(1 - (0.785/1.77)**2))
    alpha = 2.45
    den = beta**2 * alpha
    coeff = 1.0/den
    return coeff
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
def estimate_flow(pressure, K):
    flow = [np.nan]*len(pressure)
    for m in range(len(pressure)):
        if(pressure[m] > 0):
            # convert cmH2O to Pa
            f = np.sqrt(100*pressure[m]/K)
            flow[m] = f*1e3
    return flow
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

for j in range(len(subject)):
    filename = path + subject[j]

    try:
        # Extract data to empty arrays.
        # This will fail if file does not exist
        pressure_after= []
        pressure_before = []
        extract_data_arrays_from_file(filename, pressure_after, pressure_before)

        plt.plot(pressure_before)
        plt.plot(pressure_after)
        plt.show()


        flow = convert_spirometer_pressure_to_flow(pressure_after)

        delta_pressure = [pressure_before[i] - pressure_after[i] for i in range(len(pressure_after))]
        # convert from Pa to cmH20
        delta_pressure = [p/100.0 for p in delta_pressure]

        resistance = calculate_external_resistance(delta_pressure, flow)
        const_resistance = calculate_constant_resistance_sections(resistance)

        C_coeff = calculate_discharge_coefficient(flow, delta_pressure)
        print(np.nanmean(C_coeff))

        K_coeff = calculate_friction_coefficient()
        print(np.nanmean(K_coeff))

        flow_estimate = estimate_flow(delta_pressure, K_coeff)


        #plt.plot(pressure_after, 'b:')
        #plt.plot(pressure_before, 'r:')
        plt.plot(delta_pressure, 'b')
        plt.plot(flow, 'g')
        plt.plot(resistance, 'r')
        plt.plot(const_resistance, 'm')
        plt.plot(C_coeff, 'c')
        plt.plot(flow_estimate, 'k', linewidth=0.5)
        plt.ylabel('L/s, cmH2O, cmH2Os/L')
        plt.legend(['P', 'Q', 'R'])
        plt.show()

    except IOError:
        pass



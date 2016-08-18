#!/bin/bash

# Import my extensions
import sys
sys.path.insert(0, '/home/sarah/Documents/Spirometry/python/extensions')

from calculus import integral
from filters import gauss_lp_filter
from data_struct import DataStore as Data
from breath_analysis import split_breaths, calc_flow_delay, BreathData

import csv
import matplotlib.pyplot as plt
import numpy
from numpy import array, sign
from numpy.linalg import lstsq

def _crossing_flow(equations):
    ''' Work out where the two lines cross '''
    # Set initial
    crossing_flow = -1

    # Pull info from equations
    eq_1 = equations[0]
    eq_2 = equations[1]
    low_gradient = eq_1[0][0]
    low_offset = eq_1[1][0]
    high_gradient = eq_2[0][0]
    high_offset = eq_2[1][0]

    # Find point where they cross in flow
    crossing_flow = (high_offset - low_offset)/(low_gradient - high_gradient)
    return crossing_flow

def _calc_resistance(flow, grad, offset):
    ''' Get resistance from line equation '''
    resistance = grad*flow + offset
    return resistance

def _csv_write(filename, equations, section):
    ''' Write equations to csv file
        Equations expanded so csv used as look-up table'''
    # Get flow_range
    # Steps of 0.01
    if(section == 'insp'):
        flow_range = range(0, 10*100)
        flow_range = [f/100.0 for f in flow_range]
    if(section == 'exp'):
        flow_range = range(-10*100, 1)
        flow_range = [f/100.0 for f in flow_range]

    # Get flow crossing point
    xing = _crossing_flow(equations)

    with open(filename, 'w') as csvfile:
        csvwriter = csv.writer(csvfile, delimiter=',')

        # Header line
        csvwriter.writerow(['Flow', 'Resistance'])

        # Find resistance from lines
        for flow in flow_range:
            if(abs(flow) < abs(xing)):
                grad = equations[0][0][0]
                offset = equations[0][1][0]
                resistance = _calc_resistance(flow, grad, offset)
            else:
                grad = equations[1][0][0]
                offset = equations[1][1][0]
                resistance = _calc_resistance(flow, grad, offset)

            # Save flow to 2dp
            flow = '{:.2f}'.format(flow)
            # Write data to file
            csvwriter.writerow([flow, resistance])

def _append_data_to_lists(flow, pressure, Q, P):
    ''' make new list from data, removing low values '''
    min_flow = 0.1

    for i in range(len(flow)):
        if(abs(flow[i]) >= min_flow):
            Q.append(flow[i])
            P.append(pressure[i])

def _make_line_from_data(Q, R, start, stop):
    ''' Find gradient and offset of data '''
    the_R = []
    the_Q = []
    for i in range(len(Q)):
        if(abs(Q[i]) < stop
        and abs(Q[i]) > start):
            the_R.append(R[i])
            the_Q.append(Q[i])

    if(len(the_Q) > 1):
        length = len(the_R)
        ones = [1]*length
        dependent = array([the_R])
        independent = array([the_Q, ones])

        # least squares
        res = lstsq(independent.T, dependent.T)
        parameters = res[0].tolist()
    else:
        parameters = [[numpy.nan],[numpy.nan]]
    return parameters

def _remake_line(line_params, start=0, stop=10):
    ''' Returns a line for plotting '''
    line = range(int(start), int(stop))
    grad = line_params[0][0]
    off = line_params[1][0]
    for i in range(len(line)):
        line[i] = line[i]*grad + off
    return line

def resistance_equation(pressure, flow, tablename, section,
                        offset=28, plot=False):
    ''' Find resistance of input dataset.
        Save to look-up table in csv. Filename is set through 'tablename'
        Offset between pressure and flow corrected by 'offset' (flow leads)
        'plot' set to true to display plots
        'section' set to 'insp' for inspiration resistance, 'exp' for expiration resistance
        '''
    if plot:
        # Three subplots created
        f, (ax1, ax2) = plt.subplots(2)

    # filtered
    pressure = gauss_lp_filter(pressure, 125, 2)
    flow = gauss_lp_filter(flow, 125, 4)

    # Split breaths apart
    splits = split_breaths(flow)

    # Initialise storage
    grads = []
    offs = []
    Q = []
    P = []
    full_Q = []
    full_P = []

    # For each breath
    data = BreathData(flow_delay=offset)
    for breath in range(len(splits[2])):
        # Get breath info
        data.get_data(pressure, flow,
                    splits[0][breath],
                    splits[1][breath],
                    splits[2][breath],
                    125)

        # Resistance funny at start and end of section
        # Only use middle half to get linear relationship
        if(section == 'insp'):
            start = 3
            factor = 3
            insp_flow = data.insp_flow[data.insp_length/start:-data.insp_length/factor]
            insp_pressure = data.insp_pressure[data.insp_length/start:-data.insp_length/factor]

            # Record flow, pressure, volume
            _append_data_to_lists(insp_flow, insp_pressure, Q, P)
            _append_data_to_lists(data.insp_flow, data.insp_pressure, full_Q, full_P)

        elif(section == 'exp'):
            start = 3
            factor = 3
            exp_flow = data.exp_flow[data.exp_length/start:-data.exp_length/factor]
            exp_pressure = data.exp_pressure[data.exp_length/start:-data.exp_length/factor]

            # Record flow, pressure, volume
            _append_data_to_lists(exp_flow, exp_pressure, Q, P)
            _append_data_to_lists(data.exp_flow, data.exp_pressure, full_Q, full_P)

        else:
            raise ValueError('Invalid section value. Must be either \'insp\' or \'exp\'')

    # R = P/Q
    R = [P[i]/Q[i] for i in range(len(P))]
    plt.plot(full_P)
    plt.plot(full_Q, '.')
    full_R = [full_P[i]/full_Q[i] for i in range(len(full_P))]

    # Try to use very start and end data.
    # If the flow is too low, use mid-range data
    params_low = _make_line_from_data(Q, R, 0.115, 0.9)
    params_high = _make_line_from_data(Q, R, 4, 8)
    if(numpy.isnan(params_high[0][0])):
        params_high = _make_line_from_data(Q, R, 1, 4)

    # Collect equation information
    equations = [params_low, params_high]

    # Save result to csv as look-up table
    _csv_write(tablename, equations, section)

    if(plot):
        # output line values
        print('{}'.format(params_low))
        print('{}'.format(params_high))

        # align entire sets
        pressure = pressure[:len(pressure)-delay]
        flow = flow[delay:]

        # make the lines to draw
        if(section == 'insp'):
            low_line = _remake_line(params_low)
            high_line = _remake_line(params_high)
        elif(section == 'exp'):
            low_line = _remake_line(params_low, start=-10, stop = 1)
            high_line = _remake_line(params_high, start=-10, stop = 1)

        ax1.plot(full_Q, full_R, 'bd', markersize=4)
        ax1.plot(Q, R, 'yd', markersize=4)
        if(section=='insp'):
            ax1.plot(range(len(low_line)), low_line, '-m')
            ax1.plot(range(len(high_line)), high_line, '-r')
        elif(section=='exp'):
            ax1.plot(range(-1*(len(low_line)-1), 1), low_line, '-m')
            ax1.plot(range(-1*(len(high_line)-1), 1), high_line, '-r')
        ax1.grid()
        ax1.set_title('Resistance(Y) vs flow(X)')
        low_info = 'low flow (magenta): grad = {:.2}, offset = {:.2}'.format(params_low[0][0], params_low[1][0])
        high_info = 'high flow (red): grad = {:.2}, offset = {:.2}'.format(params_high[0][0], params_high[1][0])
        info_string = low_info + '\n' + high_info
        ax1.text(3, -0.06, info_string, fontsize=15)
        ax2.plot(pressure, 'b')
        ax2.grid()
        ax2.plot(flow, 'r')
        ax2.set_title('Breath pressure (b) and flow (r)')
        plt.show()

    return equations

class ResistanceTable():
    def __init__(self, filename):
        self.table = {}
        with open(filename) as csvfile:
            csvreader = csv.reader(csvfile, delimiter=',')
            for row in csvreader:
                flow = row[0]
                resistance = row[1]
                self.table[flow] = resistance

    def look_up_resistance(self, Q):
        key = '%.2f' % Q
        # if is a very small negative number
        # will break because flow table starts at 0
        if(key == '-0.00'):
            key = '0.00'
        resistance = self.table[key]
        return float(resistance)

if __name__ == "__main__":
    path = '/home/sarah/Documents/Spirometry/data/'
    plot_path = '/home/sarah/Documents/Spirometry/lungPlots/'
    files = ['DIFFUSTIK.csv']

    # Clear out registry
    reset = Data('reset')
    reset.clear_dataset_registry()

    # Create data classes
    dataset = Data('dataset_1')

    # Store pressure and flow data
    filename = path + files[0]
    dataset.store_data(filename)

    flow = dataset.flow
    pressure = dataset.pressure

    # pressure not centred at zero so shift down
    pressure = [p-0.018 for p in pressure]

    # Find the shift between the two data sets in time
    delay = calc_flow_delay(pressure, flow, plot=True)

    # Get resistance
    resistance_equation(pressure,
            flow,
            offset=delay,
            plot=True,
            tablename='resistance_table_insp.csv',
            section='insp')

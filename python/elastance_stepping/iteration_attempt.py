#!/bin/bash

# Add my extensions to path
import sys
sys.path.insert(0, '/home/sarah/Documents/Spirometry/python/extensions')

# Import my extensions
from data_struct import DataStore
from breath_analysis import BreathData, split_breaths
from calculus import integral, derivative
from filters import gauss_lp_filter

#Import built-ins
import matplotlib.pyplot as plt
from math import sin, cos, pi, exp
from numpy import isnan, array
from numpy.linalg import lstsq

def simulate_time_varying_elastance(pressure, flow, volume, R_aw):
    length = len(pressure)
    P_change = [(pressure[i] + R_aw*flow[i]) for i in range(length)]
    E_t = [flow[i]/volume[i] for i in range(length)]
    return E_t

def simulate_driving_pressure(flow, volume, E_t, R_aw, R_sp):
    length = len(flow)
    P_driving = [E_t[i]*volume[i] + (R_aw + R_sp)*flow[i]
                 for i in range(length)]
    return P_driving

def identify_R_aw_and_static_E(P_driving, flow, volume, R_sp):
    ''' Returns parameters [E, R_aw]'''
    ## Set up arrays for least squares
    length = len(P_driving)
    P_dependent = [P_driving[i] - R_sp*flow[i] for i in range(length)]
    dependent = array([P_dependent])
    independent = array([volume, flow])

    # least squares
    res = lstsq(independent.T, dependent.T)
    parameters = res[0].tolist()
    residuals = res[1].tolist()
    print('Parameters: {}'.format(parameters))
    print('Residual: {}'.format(residuals))
    return parameters

def identify_spirometer_resistance(pressure, flow):
    dependent = array([pressure])
    independent = array([flow])
    res = lstsq(independent.T, dependent.T)
    R_sp = res[0].tolist()
    R_sp = R_sp[0][0]
    return R_sp

def set_zeros_to_0_001(data_list):
    length = len(data_list)
    for i in range(length):
        if data_list[i] == 0:
            data_list[i] = 1e-3

def iterative_analysis(files, plot_path='\0'):
    ''' Using data files:
        Try estimate R
        -> calc E(t)
        -> calc Pd
        -> reestimate R and E
        -> recalc E(t)
        -> repeat until breaks or converges

        Plots are saved if plot_path given
        Data in kPa and L and s'''

    # Some constants
    # atmospheric pressure
    P_ATM = 101.3 #kPa
    # multiply kPa to get cmH2O
    kPa_to_cmH2O = 10.197
    # Guess for R_aw
    # from D = 25mm and r = 12.5mm
    # and assuming flow all nice and lovely
    INIT_R_AW = -3

    # Remove any old data sets hanging around
    # (Needed when using ipython)
    reset = DataStore()
    reset.clear_dataset_registry()

    num_files = len(files)
    for file_index in range(num_files):
        print(' ')
        print('~~~New Data Set~~~')

        # Create data storage
        dataset = DataStore()

        # Store pressure and flow data from files
        filename = files[file_index]
        dataset.store_data(filename)
        dataset.pressure = [p*kPa_to_cmH2O for p in dataset.pressure]

        # Filter data
        pressure = gauss_lp_filter(dataset.pressure, dataset.sampling_frequency, 6)
        flow = gauss_lp_filter(dataset.flow, dataset.sampling_frequency, 3)

        # Set up data arrays
        flow_splits = split_breaths(flow)
        flow_starts = flow_splits[0]
        flow_middles = flow_splits[1]
        flow_stops = flow_splits[2]

        # Iterate through each breath in dataset
        data = BreathData()
        for breath in range(1):#len(pressure_starts)-1):
            # Put data into breath
            data.get_data(pressure,
                          flow,
                          flow_starts[breath],
                          flow_middles[breath],
                          flow_stops[breath],
                          dataset.sampling_frequency
                          )

            # Set zero volume to a small number
            set_zeros_to_0_001(data.insp_volume)
            set_zeros_to_0_001(data.exp_volume)

            # work out spirometer resistance
            R_sp = identify_spirometer_resistance(data.insp_pressure, data.insp_flow)

            # initial guess of R_aw
            R_aw_insp = INIT_R_AW
            R_aw_exp = INIT_R_AW

            # initial guess of E_t
            E_t_insp = [1 for t in data.time]
            E_t_exp = [-1*e for e in reversed(E_t_insp)]


            max_iterations = 2
            for iteration in range(max_iterations):
                # Step 2
                # Simulate driving pressure
                P_driving_insp = simulate_driving_pressure(
                                data.insp_flow,
                                data.insp_volume,
                                E_t_insp,
                                R_aw_insp,
                                R_sp
                                )

                P_driving_exp = simulate_driving_pressure(
                                data.exp_flow,
                                data.exp_volume,
                                E_t_exp,
                                R_aw_exp,
                                R_sp
                                )

                # Step 3
                # Identify new R_aw and static E
                parameters_insp = identify_R_aw_and_static_E(
                                P_driving_insp,
                                data.insp_flow,
                                data.insp_volume,
                                R_sp
                                )

                parameters_exp = identify_R_aw_and_static_E(
                                P_driving_exp,
                                data.exp_flow,
                                data.exp_volume,
                                R_sp
                                )

                R_aw_insp = parameters_insp[1][0]
                E_insp = parameters_insp[0][0]
                R_aw_exp = parameters_exp[1][0]
                E_exp = parameters_exp[0][0]

                # Step 1
                # Forward simulate time varying elastance
                E_t_insp = simulate_time_varying_elastance(
                        data.insp_pressure,
                        data.insp_flow,
                        data.insp_volume,
                        R_aw_insp
                        )
                E_t_exp = simulate_time_varying_elastance(
                        data.exp_pressure,
                        data.exp_flow,
                        data.exp_volume,
                        R_aw_exp
                        )

                # Plot all the things
                plotting = True
                if(plotting):
                    # Combine half breath data into full breath data
                    full_pressure = data.insp_pressure + data.exp_pressure
                    full_P_driving = P_driving_insp + P_driving_exp
                    full_E = [E_insp]*data.insp_length + [E_exp]*data.exp_length
                    full_E_t = E_t_insp + E_t_exp
                    full_R_aw = [R_aw_insp]*data.insp_length + [R_aw_insp]*data.exp_length
                    full_R_sp = [R_sp]*data.insp_length + [R_sp]*data.exp_length

                    # 3 subplots sharing x-axis
                    f, (ax1, ax2, ax3) = plt.subplots(3, sharex=True)

                    ax1.plot(data.time, full_pressure, 'b',
                             data.time, full_P_driving, 'k')
                    ax1.grid()
                    ax1.legend(['Measured pressure',
                                'Simulated driving pressure'])

                    ax2.plot(data.time, full_E, 'b',
                             data.time, full_E_t, 'r')
                    ax2.grid()
                    ax2.legend(['Static elastance = {:.2},{:.2}'.format(E_insp, E_exp),
                                'Time varying elastance'])

                    ax3.plot(data.time, full_R_aw, 'b',
                             data.time, full_R_sp, 'g')
                    ax3.grid()
                    ax3.legend(['Airway resistance = {:.2}, {:.2}'.format(R_aw_insp, R_aw_exp),
                                'Spirometer resistance = {:.2}'.format(R_sp)])

                    #plt.ylabel("Unit")
                    #plt.xlabel("Time")
                    #plt.grid()
                    #plt.ylim([-1, 1])
                    #plt.text(3.2, 3.9, decay_info, fontsize=12)
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
    files = [path + 'Normal_1.csv', path + 'Normal_3.csv']

    # Place to save plots
    plot_path = '/home/sarah/Documents/Spirometry/lungPlots/'

    iterative_analysis(files)

#!/usr/bin/python
import time
import serial
import signal
import csv
import subprocess

def read_serial(ser, data_array):
    i = 0
    while(i < len(data_array)):
        char = ser.readline().encode('utf-8')
        data_array[i] = char
        i += 1
    return data_array

def convert_diff_voltage_to_pressure(voltage):
    kPa_per_V = 0.125
    offset_voltage = 2.244
    voltage -= offset_voltage
    pressure = voltage*kPa_per_V
    return pressure

def convert_mouth_voltage_to_pressure(voltage):
    kPa_per_V = 1.724
    offset_voltage = 0.502
    voltage -= offset_voltage
    pressure = voltage*kPa_per_V
    return pressure

def write_to_csv(csvwriter, data_array):
        for item in data_array:
            clean_item = str(item).rstrip('\r\n')
            clean_item = clean_item.split(',')

            diff_voltage = float(clean_item[0])/1000
            diff_pressure = convert_diff_voltage_to_pressure(diff_voltage)
            mouth_voltage = float(clean_item[1])/1000
            mouth_pressure = convert_mouth_voltage_to_pressure(mouth_voltage)

            csvwriter.writerow([diff_pressure, mouth_pressure])

def create_csv(filename, data_array):
    with open(filename, 'w') as csvfile:
        csvwriter = csv.writer(csvfile)
        write_to_csv(csvwriter, data_array)

def append_csv(filename, data_array):
    with open(filename, 'a') as csvfile:
        csvwriter = csv.writer(csvfile)
        write_to_csv(csvwriter, data_array)

def main():
    ser = serial.Serial('/dev/ttyUSB0', 9600)

    Fs = 90 # Hz
    record_length = 1 # seconds
    array_length = Fs*record_length
    data_array = [0]*array_length

    filename = 'Calibration_set'
    print('Reading serial')
    data_array = read_serial(ser, data_array)
    print('Writing to {}'.format(filename))
    create_csv(filename, data_array)

    for i in range(4):
        print('Reading serial')
        data_array = read_serial(ser, data_array)
        print('Writing to {}'.format(filename))
        append_csv(filename, data_array)

    data_array = [0]*array_length

if __name__ == '__main__':
    main()

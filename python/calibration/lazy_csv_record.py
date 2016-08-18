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
            clean_item = str(item).split(',')

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

    array_length = 10#samples
    data_array = [0]*array_length

    filename = 'Calibration_set'
    print('Attempting to read serial')
    reading = False
    while not reading:
        try:
            data_array = read_serial(ser, data_array)
            print('Reading serial')
            reading = True
        except UnicodeDecodeError:
            pass

    create_csv(filename, data_array)

    for i in range(10000):
        data_array = read_serial(ser, data_array)
        append_csv(filename, data_array)

    print('Max file length reached.\n Stopped reading from serial')

if __name__ == '__main__':
    try:
        main()
    except KeyboardInterrupt:
        print('Stopping serial read')

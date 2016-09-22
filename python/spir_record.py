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

def write_to_csv(csvwriter, data_array):
        for item in data_array:
            clean_item = str(item).split(',')

            mouth_pressure = float(clean_item[0])/1000.0
            flow = float(clean_item[1])/60.0

            csvwriter.writerow([mouth_pressure, flow])

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

    filename = 'spir_record.csv'
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

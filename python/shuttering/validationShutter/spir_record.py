#!/usr/bin/python
#import time
import serial
#import signal
import csv
#import subprocess

def read_serial(ser, data_array):
    i = 0
    while(i < len(data_array)):
        char = ser.readline().encode('utf-8')
        data_array[i] = char
        i += 1
    return data_array

def get_sampling_frequency(ser, max_samples):
    i = 0
    # Flush serial
    while(i < 10):
        ser.readline()
        i += 1

    i = 0
    period = 0
    count = 0
    time_last = ser.readline().decode().strip()
    print(time_last)
    time_last = time_last.split(",")
    while(i < max_samples):
        time_now = ser.readline().decode().strip()
        time_now = time_now.split(",")
        period += float(time_now[3]) - float(time_last[3])
        time_last = time_now
        count += 1
        i += 1
    freq = 1.0/(period/float(count)/1000.0)
    return freq

def write_to_csv(csvwriter, data_array):
        for item in data_array:
            clean_item = str(item).split(',')

            mouth_pressure = float(clean_item[0])/1.0
            flow = float(clean_item[1])
            pneu_pressure = float(clean_item[2])/1.0
            time = float(clean_item[3])/1000

            csvwriter.writerow([mouth_pressure, flow, pneu_pressure, time])

def create_csv(filename, data_array):
    with open(filename, 'w') as csvfile:
        csvwriter = csv.writer(csvfile)
        write_to_csv(csvwriter, data_array)

def append_csv(filename, data_array):
    with open(filename, 'a') as csvfile:
        csvwriter = csv.writer(csvfile)
        write_to_csv(csvwriter, data_array)

def main():
    #ser = serial.Serial('/dev/ttyUSB0', 115200)
    ser = serial.Serial('/dev/ttyACM0', 115200)

    array_length = 10#samples
    data_array = [0]*array_length

    filename = 'data_recording.csv'
    print('Attempting to read serial...')
    reading = False
    while not reading:
        try:
            sampling_frequency = get_sampling_frequency(ser, 80)
            print('Connection established')
            print('Sampling frequency: {}'.format(sampling_frequency))
            print('Reading serial')
            reading = True
        except UnicodeDecodeError:
            pass

    data_array = read_serial(ser, data_array)
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

"""
A serial port monitor that plots live data from
the spirometer.

The monitor receives data in the form of pascals
and mL/s
"""

import Queue
import numpy as np
import serial
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from collections import deque

from queue_manager import get_queue_item
from serial_monitor import SerialMonitorThread
from livedatafeed import LiveDataFeed

'''
TODO:
    Get Serial thread running
    Get plotter reading from serial queue
    Get plotter plotting
'''

class AnalogPlot():
    def __init__(self, max_len):
        self.ax = deque([0.0]*max_len)
        self.ay = deque([0.0]*max_len)
        self.max_len = max_len

    def addToBuf(self, buf, val):
        # Adding to right and removing from left
        # so image scrolls nicely
        if len(buf) < self.maxLen:
            buf.append(val)
        else:
            buf.popleft()
            buf.append(val)

    def add(self, data):
        # add data to 2 buffers
        assert(len(data) == 2)
        self.addToBuf(self.ax, data[0])
        self.addToBuf(self.ay, data[1])

    def update(self, frameNum, data_q, a0, a1):
        serial_item = get_queue_item(data_q)
        if serial_item is not None:
            print(serial_item)

class PlottingDataMonitor():
    def __init__(self):
        self.serial_monitor = None
        self.serial_data_q = None
        self.serial_error_q = None
        self.portname = 'tty/USB0'

    def start_up(self):
        '''
        Start serial coms
        TODO:
            want this on start/stop button?
        '''
        # Check there isn't a monitor running alreay
        # Also check there is a port specified
        if self.serial_monitor is not None or self.portname == '':
            return

        # Data and error queues
        self.serial_data_q = Queue.Queue()
        self.serial_error_q = Queue.Queue()

        # Start the serial monitor
        self.serial_monitor = SerialMonitorThread(
            self.serial_data_q,
            self.serial_error_q,
            self.portname,
            9600)
        self.serial_monitor.run()

        # If there was an error starting, shut down again
        serial_error = get_queue_item(self.serial_error_q)
        if serial_error is not None:
            print('SerialMonitorThread error')
            self.serial_monitor = None

        print('Monitor running')

    def run_the_thing(self):
        plotter = AnalogPlot(1000)

        print('Plotting data')

        # Set up the animation
        fig = plt.figure()
        ax = plt.axes(xlim=(0, 1000), ylim=(-580, 580))
        a0, = ax.plot([], [])
        a1, = ax.plot([], [])
        anim = animation.FuncAnimation(fig, plotter.update,
            fargs=(self.serial_data_q, a0, a1),
            interval=1,
            blit=True)

        # show plot
        plt.show()

        # Close serial therad
        self.serial_monitor.join()
        print('exiting.')



def main():
    monitor = PlottingDataMonitor()
    monitor.start_up()
    monitor.run_the_thing()

if __name__ == "__main__":
    main()

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
import cProfile, pstats, StringIO
import time

from queue_manager import get_queue_item
from serial_monitor import SerialMonitorThread
from livedatafeed import LiveDataFeed
from analogue_plotter import AnalogPlot

'''
TODO:
    Get Serial thread running
    Get plotter reading from serial queue
    Get plotter plotting
'''

class PlottingDataMonitor():
    def __init__(self):
        self.serial_monitor = None
        self.serial_data_q = None
        self.serial_error_q = None
        self.plotter = None
        self.analyser = None
        self.portname = '/dev/ttyUSB0'

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
        self.serial_data_q = Queue.LifoQueue()
        self.serial_error_q = Queue.Queue()

        # Start the serial monitor
        self.serial_monitor = SerialMonitorThread(
            self.serial_data_q,
            self.serial_error_q,
            self.portname,
            9600)
        self.serial_monitor.start()

        # If there was an error starting, shut down again
        serial_error = get_queue_item(self.serial_error_q)
        if serial_error is not None:
            print('SerialMonitorThread error')
            self.serial_monitor = None

        print('Monitor running')

    def run_the_thing(self):
        self.plotter = AnalogPlot(1000, 25)

        print('Plotting data')

        # Set up the animation
        #plt.ion()
        fig = plt.figure()
        ax = plt.axes()
        ax.set_xlim(0, 1000)
        ax.set_ylim(-580, 580)
        ax.grid()
        a0, = ax.plot([], [])
        a1, = ax.plot([], [])
        print('Begin animation')

        try:
#            while(1):
#                time.sleep(1)
#                print('awake')
            anim = animation.FuncAnimation(fig, self.plotter.update_serial,
                fargs=(self.serial_data_q, a0, a1),
                interval=195,
                blit=True)

            # show plot
            plt.show()
#
#            # show plot
#            plt.show()
#            while(1):
#                # This takes about 10 ms to run
#                # Updating at 25 HZ gives 40 ms period
#                # --> Sleep 30 ms then wake and plot
#                self.plotter.update_serial(fig, self.serial_data_q, a0, a1)
#                fig.canvas.draw()

        except(KeyboardInterrupt):
            # Close serial therad
            print('exiting.')

        while(self.serial_monitor.is_alive()):
            self.serial_monitor.join()



def main():
    monitor = PlottingDataMonitor()
    monitor.start_up()
    monitor.run_the_thing()

if __name__ == "__main__":
    main()

#!/usr/bin/python

import Queue
import numpy as np
import serial
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import cProfile, pstats, StringIO
import time
#import Tkinter
try:
    import wx
except ImportError:
    raise ImportError, "The wxPython module is required to run this program. \
If \"ImportError: /usr/lib/libpangoft2-1.0.so.0: undefined symbol: FcWeightToOpenType\" \
comes up after installation with conda, run: conda install -c asmeurer pango"

from queue_manager import get_queue_item
from serial_monitor import SerialMonitorThread
from livedatafeed import LiveDataFeed
from analogue_plotter import AnalogPlot
from analyser import ExprDecayAnalyser
from breath_splitter import BreathSplitter

'''
TODO:
    Get Serial thread running
    Get plotter reading from serial queue
    Get plotter plotting
'''

"""
Plotting data monitor:
    Plots data stream from serial to screen
    Data is analysed in the background and results displayed

    Gui managed by Tkinter with wxPython
"""
class PlottingDataMonitor(wx.Frame):
    def __init__(self):

        self.serial_monitor = None
        self.serial_data_q = None
        self.serial_error_q = None
        self.plotter = None
        self.analyser = None
        self.analysis_results_q = None
        self.analyser_data_q = None
        self.breath_splitter = None
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

        # Breath splitting queue
        self.breath_splitting_q = Queue.Queue()

        # Analyser input/output queues
        self.analyser_data_q = Queue.Queue()
        self.analyser_results_q = Queue.Queue()

        # Start the serial monitor
        self.serial_monitor = SerialMonitorThread(
            self.breath_splitting_q,
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

#        # Start splitting up breaths
#        self.breath_splitter = BreathSplitter(self.breath_splitting_q,
#                self.analyser_data_q)
#        self.breath_splitter.start()
#
#        # Start data analysis
#        self.analyser = ExprDecayAnalyser(self.analyser_data_q,
#                self.analyser_results_q
#                )
#        self.analyser.start()

    def run_the_thing(self):
        self.plotter = AnalogPlot(1000)

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
            anim = animation.FuncAnimation(fig, self.plotter.animate,
                fargs=(self.serial_data_q, a0, a1),
                interval=95,
                blit=True)

            # show plot
            plt.show()

        except(KeyboardInterrupt):
            # Close serial therad
            print('exiting.')

        while(self.serial_monitor.is_alive()):
            self.serial_monitor.join()
            self.breath_splitter.join()
            self.analyser.join()



def main():
    frame = PlottingDataMonitor()
    frame.start_up()
    frame.run_the_thing()

if __name__ == "__main__":
    main()

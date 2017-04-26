#!/usr/bin/python
import Queue
import threading
import time

class BreathSplitter(threading.Thread):
    """
    A thread for analysing expiration. Looks for
    exponential decays in expiration data.

    data_q:
            Contains data stream from serial

    breath_q:
            Contains lists of data for each breath segment.
            Currently only has expiration for each breath.
    """
    def __init__(self,
                 breath_q,
                 output_q
                 ):
        threading.Thread.__init__(self)

        self.breath_q = breath_q
        self.output_q = output_q

        # Internal flag for inter-thread communication
        self.alive = threading.Event()
        self.alive.set()

        print('Analyser setup')

    def run(self):
        # Any setup required at runtime goes here

        while self.alive.isSet():
            print('I\'m alive!!!!')
            time.sleep(100)

    def join(self, timeout=None):
        self.alive.clear()
        threading.Thread.join(self, timeout)

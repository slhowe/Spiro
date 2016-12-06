import Queue
from collections import deque
import time

from queue_manager import get_queue_item, clear_queue

"""
Functions to plot serial data using animation.

Takes data queue with array of two values and timestamp.
Values are plotted as two concurrent lines
"""

class AnalogPlot():
    def __init__(self, max_len, sampling_freq):
        self.pressure = deque([0.0]*max_len)
        self.flow = deque([0.0]*max_len)
        self.time = deque([0.0]*max_len)
        self.max_len = max_len
        self.update_period = 1.0/sampling_freq
        self.last_update = 0

    def addToBuf(self, buf, val):
        # Adding to right and removing from left
        # so image scrolls nicely
        if len(buf) < self.max_len:
            buf.append(val)
        else:
            buf.popleft()
            buf.append(val)

    def add(self, data):
        # add data to 2 buffers
        assert(len(data) == 2)
        self.addToBuf(self.pressure, data[0])
        self.addToBuf(self.flow, data[1])

    def update_serial(self, frameNum, data_q, a0, a1):
        #if(time.time() - self.last_update >= self.update_period):
        serial_item = get_queue_item(data_q)
        if serial_item is not None:
            self.last_update = time.time()
            print(serial_item)
            # Get time
            #self.addToBuf(self.time, serial_item[1])

            # Extract data for lines
            self.add(serial_item[0])

            #time = list(self.time)
            a0.set_data(range(self.max_len), self.pressure)
            a1.set_data(range(self.max_len), self.flow)

        #    clear_queue(data_q)
        return a0, a1

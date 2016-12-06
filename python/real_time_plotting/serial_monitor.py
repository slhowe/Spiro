import Queue
import threading
import time

import serial


class SerialMonitorThread(threading.Thread):
    """ A thread for monitoring a COM port. The COM port is
        opened when the thread is started.

        data_q:
            Queue for received data. Items in the queue are
            (data, timestamp) pairs, where data is a list
            of 2 floats received from serial, and timestamp
            is the time elapsed from the thread's start (in
            seconds).

        error_q:
            Queue for error messages. In particular, if the
            serial port fails to open for some reason, an error
            is placed into this queue.

        port:
            The COM port to open. Must be recognized by the
            system.

        port_baud/stopbits/parity:
            Serial communication parameters

        port_timeout:
            The timeout used for reading the COM port. If this
            value is low, the thread will return data in finer
            grained chunks, with more accurate timestamps, but
            it will also consume more CPU.
    """
    def __init__(   self,
                    data_q, error_q,
                    port_num,
                    port_baud,
                    port_stopbits=serial.STOPBITS_ONE,
                    port_parity=serial.PARITY_NONE,
                    port_timeout=0.05# Waits for \r\n at eol, timeout if not in 50 ms
                    ):
        threading.Thread.__init__(self)

        '''
        TODO: How to check port actually exists
        '''
        self.serial_port = None
        self.serial_arg = dict( port=port_num,
                                baudrate=port_baud,
                                stopbits=port_stopbits,
                                parity=port_parity,
                                timeout=port_timeout
                              )

        self.data_q = data_q
        self.error_q = error_q

        self.alive = threading.Event()
        self.alive.set()

        self.sampling_period = 1.0/30
        self.last_time = 0
        self.start_time = time.time()

        print('Thread Setup')

    def run(self):
        # Try connect to serial
        try:
            if self.serial_port:
                self.serial_port.close()
            self.serial_port = serial.Serial(**self.serial_arg)
        except serial.SerialException, e:
            self.error_q.put(e.message)
            return

        timestamp = 0

        while self.alive.isSet():
            '''
            TODO: Check SerialException
                  Triggered if port already used or suddenly disconnected
            '''
            try:
                line = self.serial_port.readline()
                current_time = time.time()
                if(current_time - self.last_time >= self.sampling_period):
                    data = [float(val) for val in line.split(",")]
                    if len(data) == 2:
                        timestamp = current_time - self.start_time
                        print(1/(current_time - self.last_time))
                        self.last_time = current_time
                        self.data_q.put((data, timestamp))
            except(ValueError):
                print('Could not convert serial data to float')


        # clean up
        if self.serial_port:
            self.serial_port.close()

    def join(self, timeout=None):
        self.alive.clear()
        threading.Thread.join(self, timeout)

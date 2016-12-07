import Queue
import threading

class ExprDecayAnalyser(threading.Thread):
    """
    A thread for analysing expiration. Looks for
    exponential decays in expiration data.
    """
    def __init__(self,
                 data_q,
                 results_q
                 ):
        threading.Thread.__init__(self)

        self.data_q = data_q

        # Internal flag for inter-thread communication
        self.alive = threading.Event()
        self.alive.set()

        print('Analyser setup')

    def run(self):
        # Any setup required at runtime goes here

        while self.alive.isSet():
            # Read queue to list

            # Find peak

            # Find low point

            # Fit exp decay

            # Results
            print('I\'m alive!!!!')

    def join(self, timeout=None):
        self.alive.clear()
        threading.Thread.join(self, timeout)

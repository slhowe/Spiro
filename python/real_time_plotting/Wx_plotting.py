
import matplotlib
matplotlib.use("TkAgg")
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2TkAgg
from matplotlib.figure import Figure
import matplotlib.animation as animation
from matplotlib import style

import Tkinter as tk
import ttk

import Queue

from queue_manager import get_queue_item
from serial_monitor import SerialMonitorThread
from livedatafeed import LiveDataFeed
from analogue_plotter import AnalogPlot
from analyser import ExprDecayAnalyser
from breath_splitter import BreathSplitter


class SeaofB(tk.Tk):

    def __init__(self, *args, **kwargs):

        tk.Tk.__init__(self, *args, **kwargs)
        tk.Tk.wm_title(self, "Sea of BTC client")

        self.serial_monitor = None
        self.serial_data_q = None
        self.serial_error_q = None
        self.plotter = None
        self.analyser = None
        self.analysis_results_q = None
        self.analyser_data_q = None
        self.breath_splitter = None
        self.portname = '/dev/ttyUSB0'

        container = tk.Frame(self)
        container.pack(side="top", fill="both", expand = True)
        #container.grid(row = 1, sticky='SE')
        container.grid_rowconfigure(0, weight=1)
        container.grid_columnconfigure(0, weight=1)

        self.frames = {}
        for F in (PlotPage, ResultsPage):
            # frame parent is container, controller is SeaofB
            frame = F(container, self)
            self.frames[F] = frame
            frame.grid(row=0, column=0, sticky="nsew")

        self.show_frame(PlotPage)

    def start_serial_monitor(self):
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

    def shutdown_serial_monitor(self):
        self.serial_monitor.join()

    def initialise_plot(self, length):
        # input plot window length
        self.plotter = AnalogPlot(length)

  # Select frame to show and raise it
    def show_frame(self, cont):
        frame = self.frames[cont]
        frame.tkraise()


class ResultsPage(tk.Frame):

    def __init__(self, parent, controller):
        tk.Frame.__init__(self,parent)
        label = tk.Label(self, text="Results Page", font=LARGE_FONT)
        label.grid(column=0, row=0, sticky='EW')

        button = ttk.Button(self, text="Graph page",
                            command=lambda: controller.show_frame(PlotPage))
        button.grid(column=1, row=1, sticky='EW')

        self.grid_columnconfigure(0, weight=1)
        self.grid_rowconfigure(0, weight=1)


class PlotPage(tk.Frame):

    def __init__(self, parent, controller):
        tk.Frame.__init__(self, parent)
        label = tk.Label(self, text="Graph Page! XD", font=LARGE_FONT)
        label.pack(pady=10,padx=10)

        button1 = ttk.Button(self, text="Results",
                            command=lambda: controller.show_frame(ResultsPage))
        button1.pack()


        canvas = FigureCanvasTkAgg(fig, self)
        canvas.show()
        canvas.get_tk_widget().pack(side=tk.BOTTOM, fill=tk.BOTH, expand=True)

        toolbar = NavigationToolbar2TkAgg(canvas, self)
        toolbar.update()
        canvas._tkcanvas.pack(side=tk.TOP, fill=tk.BOTH, expand=True)


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
LARGE_FONT= ("Verdana", 12)
style.use("ggplot")

fig = Figure(figsize=(5,5), dpi=100)
ax = fig.add_subplot(111)
plot_length = 100
ax.set_xlim(0, plot_length)
ax.set_ylim(-580, 580)
a0, = ax.plot([], [])
a1, = ax.plot([], [])

app = SeaofB()
app.start_serial_monitor()
app.initialise_plot(plot_length)
ani = animation.FuncAnimation(fig, app.plotter.animate,
                            fargs=(app.serial_data_q, a0, a1),
                            interval=100,
                            blit=True
                            )
app.mainloop()
app.shutdown_serial_monitor()

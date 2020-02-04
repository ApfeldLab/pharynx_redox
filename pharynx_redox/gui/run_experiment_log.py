import sys
from pathlib import Path
import logging

import numpy as np
import xarray as xr
from matplotlib import cm
from skimage.morphology import disk
from PyQt5 import QtCore, QtGui, QtWidgets
from PyQt5.QtWidgets import QFileDialog
import pyqtgraph as pg

from .. import pharynx_io as pio
from .. import image_processing as ip
from .. import experiment
from .. import utils
from .qt_py_files.experiment_run_log import Ui_Form


class QTextEditLogger(logging.Handler, QtCore.QObject):

    log_signal = QtCore.pyqtSignal(str)

    def __init__(self):
        logging.Handler.__init__(self)
        QtCore.QObject.__init__(self)

    def emit(self, record):
        msg = self.format(record)
        self.log_signal.emit(msg)


class ExperimentRunWidget(QtWidgets.QWidget):
    def __init__(self, experiment_dir):
        super(ExperimentRunWidget, self).__init__()

        # Build UI
        self.ui = Ui_Form()
        self.ui.setupUi(self)

        # Configure logging
        logTextBox = QTextEditLogger()
        logging.getLogger().addHandler(logTextBox)
        logging.getLogger().setLevel(logging.DEBUG)

        logTextBox.log_signal.connect(self.on_new_log_msg)

        self.exp_thread = RunExperimentThread(experiment_dir)
        self.exp_thread.start()

    def on_new_log_msg(self, msg):
        self.ui.logTextBox.appendPlainText(msg)


class RunExperimentThread(QtCore.QThread):
    def __init__(self, experiment_dir):
        super(RunExperimentThread, self).__init__()
        self.exp = experiment.Experiment(Path(experiment_dir))

    def run(self):
        self.exp.full_pipeline()


class Worker(QtCore.QRunnable):
    def __init__(self, fn, *args, **kwargs):
        super(Worker, self).__init__()
        self.fn = fn
        self.args = args
        self.kwargs = kwargs

    @QtCore.pyqtSlot()
    def run(self):
        self.fn(args, kwargs)


import sys
from PyQt5 import QtWidgets
import logging

# Uncomment below for terminal log messages
# logging.basicConfig(level=logging.DEBUG, format=' %(asctime)s - %(name)s - %(levelname)s - %(message)s')


app = QtWidgets.QApplication(sys.argv)
dlg = ExperimentRunWidget(
    "/Users/sean/code/pharynx_redox/data/paired_ratio/2017_02_22-HD233_SAY47"
)
dlg.show()
dlg.raise_()
sys.exit(app.exec_())

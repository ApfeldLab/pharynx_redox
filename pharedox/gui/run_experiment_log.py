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

from pharedox import io as pio
from pharedox import image_processing as ip
from pharedox import experiment
from pharedox import utils
from pharedox.qt_py_files.experiment_run_log import Ui_Form


class QTextEditLogger(logging.Handler, QtCore.QObject):

    log_signal = QtCore.pyqtSignal(str)

    def __init__(self):
        logging.Handler.__init__(self)
        QtCore.QObject.__init__(self)

    def emit(self, record):
        msg = self.format(record)
        self.log_signal.emit(msg)


class ExperimentRunWidget(QtWidgets.QWidget):
    def __init__(self):
        super(ExperimentRunWidget, self).__init__()

        # Build UI
        self.ui = Ui_Form()
        self.ui.setupUi(self)

        # Configure logging
        logTextBox = QTextEditLogger()
        logging.getLogger().addHandler(logTextBox)
        logging.getLogger().setLevel(logging.DEBUG)

        logTextBox.log_signal.connect(self.on_new_log_msg)

    def set_exp_dir(self, exp_dir):
        self.exp_dir = exp_dir
        self.exp_thread = RunExperimentThread(exp_dir)

    def on_new_log_msg(self, msg):
        self.ui.logTextBox.appendPlainText(msg)


class RunExperimentThread(QtCore.QThread):
    def __init__(self, experiment_dir):
        super(RunExperimentThread, self).__init__()
        self.exp = experiment.Experiment(Path(experiment_dir))

    def run(self):
        self.exp.full_pipeline()


if __name__ == "__main__":
    app = QtWidgets.QApplication(sys.argv)
    dlg = ExperimentRunWidget(
        "/Users/sean/code/pharedox/data/paired_ratio/2017_02_22-HD233_SAY47"
    )
    dlg.show()
    dlg.raise_()
    sys.exit(app.exec_())

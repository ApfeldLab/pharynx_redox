"""
Run the analysis with a GUI
"""

import sys
from pathlib import Path

print("SYS.PATH: ", sys.path)
import numpy as np
import xarray as xr
from matplotlib import cm
from skimage.morphology import disk
from PyQt5 import QtCore, QtGui, QtWidgets
from PyQt5.QtWidgets import QFileDialog
import pyqtgraph as pg

from pharynx_redox import pharynx_io as pio
from pharynx_redox import image_processing as ip
from pharynx_redox import utils
from pharynx_redox.qt_py_files import load_dialog
from pharynx_redox.run_experiment_log import ExperimentRunWidget


class RunExperimentDialogWidget(QtWidgets.QWidget):
    def __init__(self):
        super(RunExperimentDialogWidget, self).__init__()

        self.selected_directory = ""

        # Build UI
        self.ui = load_dialog.Ui_Dialog()
        self.ui.setupUi(self)

        self.exp_run_widget = ExperimentRunWidget()

        # Connect Signals
        self.ui.dirSelectToolButton.pressed.connect(self.handle_dir_select_pressed)
        self.ui.cancelButton.pressed.connect(self.cancel)
        self.ui.runButton.pressed.connect(self.run)
        self.ui.experimentDirectoryLineEdit.textChanged.connect(
            self.handle_exp_dir_text_changed
        )

    def handle_exp_dir_text_changed(self):
        self.selected_directory = self.ui.experimentDirectoryLineEdit.text()

    def handle_dir_select_pressed(self):
        fname = QFileDialog.getExistingDirectory(self, "Select Experiment Directory")
        if fname:
            self.selected_directory = fname
            self.ui.experimentDirectoryLineEdit.setText(fname)

    def run(self):
        self.exp_run_widget.set_exp_dir(self.selected_directory)
        self.exp_run_widget.show()
        self.exp_run_widget.exp_thread.start()

    def cancel(self):
        self.close()


if __name__ == "__main__":
    qapp = QtWidgets.QApplication(sys.argv)

    main_dialog = RunExperimentDialogWidget()
    main_dialog.show()

    qapp.exec_()

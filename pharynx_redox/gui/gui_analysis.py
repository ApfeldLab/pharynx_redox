"""
Run the analysis with a GUI
"""

import sys
from pathlib import Path

import numpy as np
import xarray as xr
from matplotlib import cm
from skimage.morphology import disk
from PyQt5 import QtCore, QtGui, QtWidgets
from PyQt5.QtWidgets import QFileDialog
import pyqtgraph as pg

from .. import pharynx_io as pio
from .. import image_processing as ip
from .. import utils
from .qt_py_files.load_dialog import Ui_Dialog


class RunExperimentDialogWidget(QtWidgets.QWidget):
    def __init__(self):
        super(RunExperimentDialogWidget, self).__init__()

        self.selected_directory = ""

        # Build UI
        self.ui = Ui_Dialog()
        self.ui.setupUi(self)

        # Connect Signals
        self.ui.dirSelectToolButton.pressed.connect(self.handle_dir_select_pressed)

    def handle_dir_select_pressed(self):
        if self.selected_directory:
            dir_ = self.selected_directory
        else:
            dir_ = "/"
        fname = QFileDialog.getExistingDirectory(
            self, "Select Experiment Directory", dir_
        )
        if fname:
            self.selected_directory = fname
            self.ui.experimentDirectoryLineEdit.setText(fname)

    def accept(self):
        pass

    def reject(self):
        self.close()


if __name__ == "__main__":
    qapp = QtWidgets.QApplication(sys.argv)

    main_dialog = RunExperimentDialogWidget()
    main_dialog.show()

    qapp.exec_()

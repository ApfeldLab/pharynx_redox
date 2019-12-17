from pathlib import Path

import numpy as np
import pyqtgraph as pg
from PyQt5 import QtWidgets, QtGui, QtCore
from PyQt5.QtWidgets import QFileDialog

from pharynx_redox import pharynx_io as pio

imgs = pio.load_images('/Users/sean/Desktop/2019-12-13 calibration 1-6 baseline 10 fixed.tif', strain_map=np.repeat('HD233', 6))

qapp = QtWidgets.QApplication([])


# qapp.exec_()

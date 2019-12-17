import sys
from pathlib import Path

import numpy as np
import pyqtgraph as pg
import xarray as xr
from PyQt5 import QtCore, QtGui, QtWidgets
from PyQt5.QtWidgets import QFileDialog
from .qt_py_files.image_stack_widget import Ui_XArrayDisplayWidget

from .. import pharynx_io as pio
from .. import utils


class ImageWidget(QtWidgets.QWidget):
    def __init__(self, imgs):
        super(ImageWidget, self).__init__()
        
        self.imgs = imgs
        self.current_animal = 0

        # Build UI
        self.ui = Ui_XArrayDisplayWidget()
        self.ui.setupUi(self)

        self.set_image_stack(imgs.isel(pair=0, wavelength=0), update_scale=True, autoLevels=True) 
        self.ui.pairSlider.setMinimum(0)
        self.ui.pairSlider.setMaximum(self.imgs.pair.size - 1)
        self.ui.wvlBox.addItems(imgs.wavelength.values)

        # connect signals
        self.ui.pairSlider.valueChanged.connect(self.handle_pairSlider_changed)
        self.ui.wvlBox.currentIndexChanged.connect(self.handle_wvlBox_change)
    
    def set_current_animal(self, animal_no):
        self.current_animal = animal_no
        self.ui.ImageViewBox.setCurrentIndex(animal_no)

    def handle_animal_slider_changed(self):
        self.current_animal = self.ui.ImageViewBox.currentIndex

    def handle_pairSlider_changed(self):
        self.set_image_stack(self.imgs.sel(pair=self.ui.pairSlider.value(), wavelength=self.get_current_wvl()))

    def handle_wvlBox_change(self):
        self.set_image_stack(self.imgs.sel(pair=self.get_current_pair(), wavelength=self.get_current_wvl()))

    def get_current_wvl(self):
        return self.ui.wvlBox.currentText()

    def set_current_pair(self, pair):
        self.ui.pairSlider.setValue(pair)
        self.set_image_stack(self.imgs.sel(wavelength=self.get_current_wvl(), pair=pair))
    
    def get_current_pair(self):
        return self.ui.pairSlider.value()

    def set_image_stack(self, img_stack: xr.DataArray, update_scale=False, autoLevels=False, transpose=True):
        # Save scale
        _view = self.ui.ImageViewBox.getView()
        _state = _view.getState()
        _animal_idx = self.current_animal

        if autoLevels:
            levels = None
        else:
            levels = self.ui.ImageViewBox.getHistogramWidget().getLevels()

        if transpose:
            imgdata = img_stack.transpose('animal', 'x', 'y').values
        else:
            imgdata = img_stack.values

        self.ui.ImageViewBox.setImage(imgdata, autoLevels=autoLevels, levels=levels)

        if not update_scale:
            # Restore scale
            _view.setState(_state)
        
        self.set_current_animal(_animal_idx)



# imgs = pio.load_images(
#     "/Users/sean/code/pharynx_redox/data/timeseries/2019-12-13 calibration 1-6 diamide 30/2019-12-13 calibration 1-6 diamide 30.tif",
#     strain_map=np.repeat("HD233", 6),
# )
# imgs[0:5].to_netcdf('~/Desktop/test_imgs.nc')

imgs = xr.load_dataarray('~/Desktop/test_seg.nc')

qapp = QtWidgets.QApplication(sys.argv)

xarray_widget = ImageWidget(imgs)
xarray_widget.show()

qapp.exec_()

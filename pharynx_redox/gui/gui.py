"""
A Graphical User Interface for the pipeline
"""

import sys
from pathlib import Path

import numpy as np
import xarray as xr
from PyQt5 import QtCore, QtGui, QtWidgets
from PyQt5.QtWidgets import QFileDialog

from .. import pharynx_io as pio
from .. import utils
from .qt_py_files.image_stack_widget import Ui_XArrayDisplayWidget


class ImageStackWidget(QtWidgets.QWidget):
    def __init__(self, imgs):
        super(ImageStackWidget, self).__init__()

        self.imgs = imgs
        self.current_animal = 0

        # Build UI
        self.ui = Ui_XArrayDisplayWidget()
        self.ui.setupUi(self)

        self.ui.pairSlider.setMinimum(0)
        self.ui.pairSlider.setMaximum(self.imgs.pair.size - 1)
        self.ui.wvlBox.addItems(imgs.wavelength.values)

        self.set_image_stack(
            imgs.isel(pair=0, wavelength=0), update_scale=True, autoLevels=True
        )

        # connect signals
        self.ui.pairSlider.valueChanged.connect(self.handle_pairSlider_changed)
        self.ui.wvlBox.currentIndexChanged.connect(self.handle_wvlBox_change)

    def get_current_img_frame(self):
        return self.imgs.sel(
            wavelength=self.get_current_wvl(),
            pair=int(self.get_current_pair()),
            animal=self.current_animal,
        )

    def update_img_properties_tree(self):
        img = self.get_current_img_frame()
        frame_data = {coord: img.coords[coord].values[()] for coord in img.coords}
        self.ui.propertiesDataTreeWidget.setData(frame_data)

    def set_current_animal(self, animal_no):
        self.current_animal = animal_no
        self.ui.ImageViewBox.setCurrentIndex(animal_no)
        self.update_img_properties_tree()

    def handle_animal_slider_changed(self):
        self.current_animal = self.ui.ImageViewBox.currentIndex
        self.update_img_properties_tree()

    def handle_pairSlider_changed(self):
        self.set_image_stack(
            self.imgs.sel(
                pair=self.ui.pairSlider.value(), wavelength=self.get_current_wvl()
            )
        )

    def handle_wvlBox_change(self):
        self.set_image_stack(
            self.imgs.sel(
                pair=self.get_current_pair(), wavelength=self.get_current_wvl()
            )
        )

    def get_current_wvl(self) -> str:
        return self.ui.wvlBox.currentText()

    def set_current_pair(self, pair):
        self.ui.pairSlider.setValue(pair)
        self.set_image_stack(
            self.imgs.sel(wavelength=self.get_current_wvl(), pair=pair)
        )

    def get_current_pair(self) -> int:
        return int(self.ui.pairSlider.value())

    def set_image_stack(
        self,
        img_stack: xr.DataArray,
        update_scale=False,
        autoLevels=False,
        transpose=True,
    ):
        # Save scale
        _view = self.ui.ImageViewBox.getView()
        _state = _view.getState()
        _animal_idx = self.current_animal

        if autoLevels:
            levels = None
        else:
            levels = self.ui.ImageViewBox.getHistogramWidget().getLevels()

        if transpose:
            imgdata = img_stack.transpose("animal", "x", "y").values
        else:
            imgdata = img_stack.values

        self.ui.ImageViewBox.setImage(imgdata, autoLevels=autoLevels, levels=levels)

        if not update_scale:
            # Restore scale
            _view.setState(_state)

        self.set_current_animal(_animal_idx)


if __name__ == "__main__":
    imgs = pio.load_images(
        "/Users/sean/code/pharynx_redox/data/timeseries/2019-12-13 calibration 1-6 baseline 10/2019-12-13 calibration 1-6 baseline 10.tif",
        strain_map=np.repeat("HD233", 6),
    )
    # imgs[0:5].to_netcdf('~/Desktop/test_imgs.nc')

    # imgs = xr.load_dataarray('~/Desktop/test_seg.nc')

    qapp = QtWidgets.QApplication(sys.argv)

    xarray_widget = ImageStackWidget(imgs)
    xarray_widget.show()

    qapp.exec_()

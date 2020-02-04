"""
A Graphical User Interface for the pipeline
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
from .qt_py_files.image_stack_widget import Ui_XArrayDisplayWidget


class ImageStackWidget(QtWidgets.QWidget):
    def __init__(self, imgs, masks):
        super(ImageStackWidget, self).__init__()

        self.imgs = imgs.transpose(
            "animal", "wavelength", "pair", "x", "y", transpose_coords=True
        )
        self.masks = masks.transpose(
            "animal", "wavelength", "pair", "x", "y", transpose_coords=True
        )
        self.current_animal = 0
        self.display_mask = True

        pos = np.array([0, 1])
        color = np.array([[0, 0, 0, 0], [255, 0, 0, 128]], dtype=np.ubyte)
        cmap = pg.ColorMap(pos, color)
        self.lut = cmap.getLookupTable(0, 1, 2)

        # EDITING STATE
        self.is_editing_mask = False
        self.is_erasing = False
        self.mask_draw_radius = 5
        self.kernel = disk(self.mask_draw_radius)

        # Build UI
        self.ui = Ui_XArrayDisplayWidget()
        self.ui.setupUi(self)

        self.ui.pairSlider.setMinimum(0)
        self.ui.pairSlider.setMaximum(self.imgs.pair.size - 1)
        self.ui.wvlBox.addItems(imgs.wavelength.values)
        self.ui.drawToolButton.setCheckable(True)

        self.mask_item = pg.ImageItem(
            image=self.masks.sel(
                wavelength=self.get_current_wvl(),
                pair=self.get_current_pair(),
                animal=0,
            ).values
        )
        self.mask_item.setLookupTable(self.lut)

        self.set_image_stack(
            imgs.isel(pair=0, wavelength=0), update_scale=True, autoLevels=True
        )

        # connect signals
        self.ui.pairSlider.valueChanged.connect(self.handle_pairSlider_changed)
        self.ui.wvlBox.currentIndexChanged.connect(self.handle_wvlBox_change)
        self.ui.displayMaskCheckbox.stateChanged.connect(
            self.handle_display_mask_pressed
        )
        self.ui.editMaskCheckBox.stateChanged.connect(self.handle_edit_mask_check)
        self.ui.ImageViewBox.sigTimeChanged.connect(self.handle_animal_slider_changed)
        self.ui.drawToolButton.pressed.connect(self.handle_pressed_draw_mask)

    def handle_edit_mask_check(self):
        self.is_editing_mask = self.ui.editMaskCheckBox.isChecked()
        if self.is_editing_mask:
            if self.is_erasing:
                kernel = 1 - self.kernel
            else:
                kernel = self.kernel
            self.mask_item.setDrawKernel(
                kernel,
                mask=self.kernel,
                center=(self.mask_draw_radius, self.mask_draw_radius),
                mode="set",
            )
        else:
            self.mask_item.setDrawKernel(None)

    def keyPressEvent(self, ev):
        print(ev.key())

        if ev.key() == QtCore.Qt.Key_M:
            self.is_editing_mask = not self.is_editing_mask
            self.ui.drawToolButton.setChecked(self.is_editing_mask)
            self.ui.drawToolButton.pressed.emit()
        if ev.key() == QtCore.Qt.Key_Minus:
            self.is_erasing = not self.is_erasing
            self.ui.drawToolButton.pressed.emit()
            print(f"is erasing? {self.is_erasing}")
        else:
            super(ImageStackWidget, self).keyPressEvent(ev)

    def handle_pressed_draw_mask(self):
        print(f"drawing? {self.is_editing_mask}")
        print(f"erasing? {self.is_erasing}")

    def handle_display_mask_pressed(self):
        self.display_mask = self.ui.displayMaskCheckbox.isChecked()
        print(f"display_mask={self.display_mask}")
        if self.display_mask:
            self.ui.ImageViewBox.addItem(self.mask_item)
            self.ui.drawToolButton.setEnabled(True)
        else:
            self.ui.ImageViewBox.removeItem(self.mask_item)
            self.ui.drawToolButton.setChecked(False)
            self.ui.drawToolButton.setEnabled(False)
            self.is_editing_mask = False

    def get_current_img_frame(self):
        return self.imgs.sel(
            wavelength=self.get_current_wvl(),
            pair=int(self.get_current_pair()),
            animal=self.current_animal,
        )

    def get_current_mask_frame(self):
        return self.masks.sel(
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
        self.mask_item.setImage(self.get_current_mask_frame().values)

    def handle_animal_slider_changed(self):
        self.current_animal = self.ui.ImageViewBox.currentIndex
        print(f"current animal: {self.current_animal}")
        self.mask_item.setImage(self.get_current_mask_frame().values)

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
        self, img_stack: xr.DataArray, update_scale=False, autoLevels=False
    ):
        # Save scale
        _view = self.ui.ImageViewBox.getView()
        _state = _view.getState()
        _animal_idx = self.current_animal

        if autoLevels:
            levels = None
        else:
            levels = self.ui.ImageViewBox.getHistogramWidget().getLevels()

        imgdata = img_stack.values

        self.ui.ImageViewBox.setImage(
            imgdata, autoLevels=autoLevels, levels=levels,
        )
        self.ui.ImageViewBox.getImageItem().setBorder({"color": "FF0", "width": 2})

        # if self.display_mask:

        if not update_scale:
            # Restore scale
            _view.setState(_state)

        self.set_current_animal(_animal_idx)


if __name__ == "__main__":
    imgs = pio.load_images(
        "/Users/sean/code/pharynx_redox/data/paired_ratio/2017_02_22-HD233_SAY47/2017_02_22-HD233_SAY47.tif",
        indexer_path="/Users/sean/code/pharynx_redox/data/paired_ratio/2017_02_22-HD233_SAY47/2017_02_22-HD233_SAY47-indexer.csv",
    )
    masks = ip.segment_pharynxes(imgs, threshold=2000)

    qapp = QtWidgets.QApplication(sys.argv)

    xarray_widget = ImageStackWidget(imgs, masks)
    xarray_widget.show()

    qapp.exec_()

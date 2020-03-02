import sys
import numpy as np
import xarray as xr
from matplotlib import cm
from skimage.morphology import disk
from PyQt5 import QtCore, QtGui, QtWidgets
from PyQt5.QtWidgets import QFileDialog
import pyqtgraph as pg
from dataclasses import dataclass

from pharynx_redox import pharynx_io as pio

from .qt_py_files.image_stack_widget import Ui_XArrayDisplayWidget

import logging

logging.basicConfig(
    stream=sys.stdout,
    level=logging.INFO,
    format="%(levelname)s:%(name)s:%(funcName)s:%(message)s",
)


@dataclass
class ImgStackState:
    wvl: str
    frame: int = 0
    pair: int = 0


class NDImageStackWidget(QtWidgets.QWidget):
    def __init__(self, imgs, midlines=None):
        """
        A widget to display n-dimensional image stacks
        
        Parameters
        ----------
        imgs : dict
            a dictionary with the following structure:
                {
                    'raw images': {
                        'data': xr.DataArray,
                        'masks': xr.DataArray,
                    },
                    'rotated images': {
                        'data': xr.DataArray,
                        'masks': xr.DataArray,
                    }
                }
            
            The `data` and `masks` DataArrays must have the same shape.
            There may be an arbitrary number of entries.
        midlines : xr.DataArray
            The midlines associated with the images. Must have the same shape as the
            image data.
        """
        super(NDImageStackWidget, self).__init__()

        self.imgs = imgs
        self.midlines = midlines

        self.img_stack_state = {}
        for stack_name in imgs.keys():
            self.img_stack_state[stack_name] = ImgStackState(
                wvl=self.imgs[stack_name]["data"].wavelength.values[0]
            )

        # Build UI
        self.ui = Ui_XArrayDisplayWidget()
        self.ui.setupUi(self)

        self.set_tabs()
        self.update_wvl_box()

        # Threshold line
        threshold_line = pg.InfiniteLine(angle=0, movable=True, pen="g")
        self.ui.tabWidget.currentWidget().getHistogramWidget().vb.addItem(
            threshold_line
        )
        threshold_line.setValue(1000)
        threshold_line.setZValue(1000)

        # Editing State Machine
        self.state_machine = QtCore.QStateMachine()

        # initialize states
        base_state = QtCore.QState()
        mask_editing_state = QtCore.QState()
        midline_editing_state = QtCore.QState()

        # connect state transitions
        # base_state.addTransition(base_state, self.ui.maskCheckBox.)

        # Connect signals
        self.ui.tabWidget.currentChanged.connect(self.handle_tab_switch)
        self.ui.wvlBox.currentIndexChanged.connect(self.set_wavelength)
        self.ui.flipButton.clicked.connect(self.handle_flip_clicked)
        self.ui.excludeButton.clicked.connect(self.handle_exclude_clicked)
        self.ui.paintButton.clicked.connect(self.handle_paint_clicked)
        self.ui.maskCheckBox.stateChanged.connect(
            self.handle_mask_checkbox_state_changed
        )
        self.ui.midlineCheckBox.stateChanged.connect(
            self.handle_midline_checkbox_state_changed
        )
        # for i in range(self.ui.tabWidget.count()):
        # img_view = self.ui.tabWidget.widget(i)
        # img_view.sigTimeChanged.connect(self.handle_frame_slider_changed)

    def enter_base_state(self):
        print("entering base state")

    def enter_mask_editing_state(self):
        print("entering mask editing state")

    def enter_midline_editing_state(self):
        print("entering midline editing state")

    def handle_midline_checkbox_state_changed(self):
        print(f"midline checkbox is {self.ui.midlineCheckBox.isChecked()}")

    def handle_mask_checkbox_state_changed(self):
        print(f"mask checkbox is {self.ui.maskCheckBox.isChecked()}")

    def handle_exclude_clicked(self):
        print("exclude clicked")

    def handle_paint_clicked(self):
        print("paint clicked")

    def handle_flip_clicked(self):
        print("flip clicked")

    def set_wavelength(self, new_wvl_idx: int, stack_name: str = None):
        """
        Set the wavelength of the image stack with the given stack name.
        
        Parameters
        ----------
        new_wvl_idx : int
            the index of the new wavelength with respect to the image stack that is
            being updated
        stack_name : str, optional
            the name of the image stack to be updated (in the `self.imgs` dict), maps
            onto `self.ui.tabWidget`'s tab names. If `None`, uses the current tab's
            label.
        """
        if stack_name is None:
            stack_name = self.get_current_stack_name()
        img_data = self.imgs[stack_name]["data"]

        self.img_stack_state[stack_name].wvl = img_data.wavelength.values[new_wvl_idx]
        self.update_img_stacks()

    def handle_frame_slider_changed(self, ind, time):
        """
        This function is called whenever the frame slider is changed
        """
        pass

    def handle_tab_switch(self, idx):
        print(f"switching to tab: {self.ui.tabWidget.tabText(idx)}")
        self.update_wvl_box()

    def update_img_stacks(self, link_views=False):
        """
        Update the wavelength/pair of the hyper stacks according to their states as
        kept track of by `self.img_stack_state`
        
        Parameters
        ----------
        link_views : bool, optional
            if X/Y axis of the views should be linked together, by default False. 

            If True, will override the user-input of linking/unlinking that can be done
            by right-clicking the image in the GUI.

            If False, the user's preferences are maintained.
        """
        for i in range(self.ui.tabWidget.count()):
            img_view = self.ui.tabWidget.widget(i)
            stack_name = self.ui.tabWidget.tabText(i)
            stack_state = self.img_stack_state[stack_name]
            selector = {
                "wavelength": stack_state.wvl,
                "pair": stack_state.pair,
            }
            img_view.setImage(self.imgs[stack_name]["data"].sel(**selector).values)
            img_view.getImageItem().setBorder(dict(color="FF0", width=1))
            img_view.setCurrentIndex(stack_state.frame)

    def update_wvl_box(self):
        """
        Get the wavelengths from the current tab's image data and set the wavelength
        selection box accordingly
        """
        stack_name = self.get_current_stack_name()
        print(f"updating wvl box to reflect {stack_name}")

        img_data = self.get_current_img_hyperstack()
        wvls = img_data.wavelength.values

        self.set_wvl_box_list(wvls)

        wvl = self.img_stack_state[stack_name].wvl
        idx = self.ui.wvlBox.findText(wvl)
        if idx != -1:
            self.ui.wvlBox.setCurrentIndex(idx)

    def get_current_img_hyperstack(self) -> xr.DataArray:
        """
        Returns the "hyperstack" that the user is currently looking at.         

        """
        stack_name = self.ui.tabWidget.tabText(self.ui.tabWidget.currentIndex())
        return self.imgs[stack_name]["data"]

    def set_wvl_box_list(self, wvls):
        self.ui.wvlBox.clear()
        self.ui.wvlBox.addItems(wvls)

    def set_tabs(self):
        # First, clear all tabs
        self.ui.tabWidget.clear()

        for stack_name, stack_data in self.imgs.items():
            iv = pg.ImageView(name=stack_name)
            self.ui.tabWidget.addTab(iv, stack_name)

        self.update_img_stacks(link_views=True)

    def get_current_stack_name(self):
        """
        Returns the `stack_name` from the label of the current tab. Lets you index into
        `self.imgs` and `self.img_stack_state`.
        """
        return self.ui.tabWidget.tabText(self.ui.tabWidget.currentIndex())


if __name__ == "__main__":

    pg.setConfigOptions(imageAxisOrder="row-major")

    imgs = pio.load_images(
        "/Users/sean/code/pharynx_redox/data/paired_ratio/2017_02_22-HD233_SAY47/2017_02_22-HD233_SAY47.tif",
        indexer_path="/Users/sean/code/pharynx_redox/data/paired_ratio/2017_02_22-HD233_SAY47/2017_02_22-HD233_SAY47-indexer.csv",
    )
    imgs[:10].to_netcdf("~/Desktop/testing_img.nc")
    # imgs = xr.load_dataarray("~/Desktop/testing_img.nc")

    qapp = QtWidgets.QApplication(sys.argv)

    widget = NDImageStackWidget(
        imgs={"raw": {"data": imgs}, "rotated": {"data": imgs},}
    )
    widget.show()

    qapp.exec_()

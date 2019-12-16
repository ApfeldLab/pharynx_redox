from PyQt5 import QtGui
from PyQt5 import QtWidgets
import pyqtgraph as pg

# from pharynx_redox import pharynx_io as pio
from pathlib import Path
import xarray as xr
import sys


class PolySplineROI(pg.PolyLineROI):
    pass


class ImageStackWidget(PyQt5.QWidget):
    def __init__(self, image_data: xr.DataArray, *args, **kwargs):
        super(ImageStackWidget, self).__init__(*args, **kwargs)

        self.layout = QtWidgets.QLayout()

        self.imv = pg.ImageView()

        self._img_data = image_data
        self.pair = image_data.pair.values[0]
        self.wavelength = image_data.wavelength.values[0]
        self.selected_stack_data = self.image_data.sel(
            wavelength=self.wavelength, pair=self.pair
        ).values

        self.imv.setImage(self.selected_stack_data)


img_data = xr.load_dataarray(
    "/Users/sean/code/pharynx_redox/data/paired_ratio/all_rot_fl.nc"
)

## Always start by initializing Qt (only once per application)
app = QtGui.QApplication([])

## Define a top-level widget to hold everything
w = QtGui.QWidget()

## Create some widgets to be placed inside
wvl_selector = pg.ComboBox(items=list(img_data.wavelength.values))
pair_selector = pg.ComboBox(items=list(img_data.pair.values.astype(str)))

# ROI
midline_roi = PolySplineROI([(10, 10), (20, 20), (30, 30)])
midline_roi

imv = pg.ImageView()
animal_label = pg.LabelItem(text="Animal 0 | Pair 0")
midline_plot = pg.PlotDataItem(x=[50, 120], y=[60, 60])

imv.getView().addItem(animal_label)
imv.getView().addItem(midline_plot)
imv.getView().addItem(midline_roi)

I = (
    (img_data.sel(wavelength="410") / img_data.sel(wavelength="470"))
    .sel(pair=0)
    .values
    # .transpose("animal", "x", "y")
)


def show_error_msg(message):
    msg = QtWidgets.QMessageBox()
    msg.setIcon(QtWidgets.QMessageBox.Critical)
    msg.setText(message)
    msg.setInformativeText("More information")
    msg.setWindowTitle("Error")
    msg.exec_()


def update_img_data():
    idx = get_current_idx()
    wvl = get_current_wvl()
    pair = get_current_pair()
    try:
        imv.setImage(
            img_data.sel(wavelength=wvl, pair=pair).transpose("animal", "x", "y").values
        )
        imv.setCurrentIndex(idx)
        set_img_label()
    except:
        print("caught an error")


def get_current_wvl():
    return wvl_selector.value()


def get_current_pair():
    return int(pair_selector.value())


def get_current_idx():
    return imv.currentIndex


def set_img_label():
    animal_label.setText(
        f"Animal {get_current_idx()} | {get_current_wvl()} | Pair {get_current_pair()}"
    )


# SIGNALS
imv.timeLine.sigPositionChanged.connect(set_img_label)
wvl_selector.currentIndexChanged.connect(update_img_data)
pair_selector.currentIndexChanged.connect(update_img_data)

## Create a grid layout to manage the widgets size and position
layout = QtGui.QGridLayout()
w.setLayout(layout)

## Add widgets to the layout in their proper positions
layout.addWidget(wvl_selector, 0, 0)
layout.addWidget(pair_selector, 1, 0)
layout.addWidget(imv, 0, 1, 3, 1)

# Initialize
update_img_data()

## Display the widget as a new window
w.show()

## Start the Qt event loop
app.exec_()

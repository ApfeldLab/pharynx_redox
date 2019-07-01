from PyQt5 import QtWidgets, QtCore
import pyqtgraph as pg
import numpy as np
from PyQt5.QtCore import Qt
from PyQt5.QtWidgets import QGroupBox, QRadioButton, QSlider, QVBoxLayout


def create_initialization_input_form(default_imaging_scheme_str):
    # Create the input form
    input_form_layout = QtWidgets.QFormLayout()
    imaging_scheme_input = QtWidgets.QLineEdit(default_imaging_scheme_str)

    # select file button
    btn = QtWidgets.QPushButton('Select Raw Image Path')
    raw_img_path_input = QtWidgets.QLineEdit('')

    input_form_layout.addRow('Imaging Scheme', imaging_scheme_input)
    input_form_layout.addRow(btn, raw_img_path_input)

    return input_form_layout


def create_strain_map_table():
    # TODO: get adding rows working
    tbl = QtWidgets.QTableWidget()
    tbl.setColumnCount(3)
    tbl.setHorizontalHeaderLabels(['Strain', 'Start Animal', 'End Animal'])
    tbl.insertRow(0)
    return tbl


def create_input_group_box():
    gb = QtWidgets.QGroupBox("Input && Initialization")

    vbox_layout = QtWidgets.QVBoxLayout()

    vbox_layout.addLayout(create_initialization_input_form())
    vbox_layout.addWidget(create_strain_map_table())

    gb.setLayout(vbox_layout)
    return gb


def create_left_column():
    left_column_layout = QtWidgets.QVBoxLayout()
    left_column_layout.addWidget(create_input_group_box(), alignment=QtCore.Qt.AlignLeft)

    return left_column_layout


def create_linked_image_view(raw_imgs):
    linked_image_view_layout = QtWidgets.QVBoxLayout()

    # Images
    multi_image_layout = pg.GraphicsLayoutWidget()
    n_cols = 2
    vboxes = []
    for i, wvl in enumerate(np.unique(raw_imgs.wavelength)):
        img_ = np.swapaxes(raw_imgs.sel(wavelength=wvl).data[0], 0, 1)
        imi = pg.ImageItem(image=img_)
        vb = pg.ViewBox(lockAspect=True)
        vb.addItem(imi)
        vboxes.append(vb)
        vb.scaleBy(1)
        if i > 0:
            vb.linkView(pg.ViewBox.YAxis, vboxes[0])
            vb.linkView(pg.ViewBox.XAxis, vboxes[0])
        multi_image_layout.addItem(vb, col=i % n_cols, row=np.floor(i / n_cols))

    # Slider
    slider = QtWidgets.QSlider(QtCore.Qt.Horizontal)

    # build layout
    linked_image_view_layout.addWidget(multi_image_layout)
    linked_image_view_layout.addWidget(slider)

    return linked_image_view_layout


class LinkedImageWidget(QtWidgets.QWidget):

    def __init__(self, image_data, wavelengths, parent=None, flags=None):
        QtWidgets.QWidget.__init__(self, parent, flags)
        self.image_data = image_data
        self.wavelengths = wavelengths
        self.frame = 0
        self.layout = QtWidgets.QGraphicsGridLayout()

        # View Initialization
        multi_image_layout = pg.GraphicsLayoutWidget()
        n_cols = 2
        v_boxes = []
        for i, wvl in enumerate(np.unique(image_data.wavelength)):
            img_ = np.swapaxes(image_data.sel(wavelength=wvl).data[0], 0, 1)
            imi = pg.ImageItem(image=img_)
            vb = pg.ViewBox(lockAspect=True)
            vb.addItem(imi)
            v_boxes.append(vb)
            vb.scaleBy(1)
            if i > 0:
                vb.linkView(pg.ViewBox.YAxis, v_boxes[0])
                vb.linkView(pg.ViewBox.XAxis, v_boxes[0])
            multi_image_layout.addItem(vb, col=i % n_cols, row=np.floor(i / n_cols))

        # Slider
        self.slider = QtWidgets.QSlider(QtCore.Qt.Horizontal)

    def update_frame(self, frame):
        self.frame = frame

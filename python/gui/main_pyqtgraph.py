import pickle

import numpy as np
import pyqtgraph as pg
from PyQt5 import QtWidgets, QtCore
from pyqtgraph import GraphicsLayoutWidget, ImageItem, PlotWidget

from gui.gui_pyqtgraph import Ui_MainWindow


class ImageGridWidget(GraphicsLayoutWidget):

    def __init__(self, image_stack_set, **kwargs):
        super(ImageGridWidget, self).__init__(**kwargs)
        self.image_stack_set = image_stack_set
        self.wavelengths = image_stack_set.wavelength.data
        self.frame = 0

        self.image_items = {
            wvl: ImageItem(image=image_stack_set.sel(wavelength=wvl)[self.frame].data, border='w')
            for wvl in self.wavelengths
        }

        self.image_viewboxes = {}
        for i, (wvl, image_item) in enumerate(self.image_items.items()):
            if (i > 0) and (i % 2 == 0):
                self.nextRow()
            self.image_viewboxes[wvl] = self.addViewBox(name=wvl)
            self.image_viewboxes[wvl].addItem(self.image_items[wvl])
            self.image_viewboxes[wvl].setAspectLocked(True)

    def set_frame(self, frame):
        self.frame = frame
        for wvl in self.wavelengths:
            self.image_items[wvl].setImage(self.image_stack_set.sel(wavelength=wvl)[self.frame].data)


class ProfilePlotGridWidget(GraphicsLayoutWidget):
    def __init__(self, profile_data, **kwargs):
        super(ProfilePlotGridWidget, self).__init__(**kwargs)
        self.frame = 0
        self.profile_data = profile_data
        self.wavelengths = self.profile_data.wavelength.data
        self.xs = np.linspace(1, 100, self.profile_data.shape[2])
        self.plots = {}
        for i, wvl in enumerate(self.wavelengths):
            if (i>0) and (i % 2 == 0):
                self.nextRow()
            self.plots[wvl] = self.addPlot(x=self.xs, y=self.profile_data.sel(wavelength=wvl)[self.frame].data)

    def set_frame(self, frame):
        self.frame = frame
        for wvl in self.wavelengths:
            self.plots[wvl].plot(x=self.xs, y=self.profile_data.sel(wavelength=wvl)[self.frame].data, clear=True)


    def _mean_data(self):
        pass


class MainWindow(QtWidgets.QMainWindow):

    def __init__(self):
        super(MainWindow, self).__init__()

        self.ui = Ui_MainWindow()
        self.ui.setupUi(self)

        self.frame = 0

        # self.experiment = experiment.PairExperiment(img_path, "TL/470_1/410_1/470_2/410_2", strains)
        self.experiment = pickle.load(open('/Users/sean/code/wormAnalysis/data/experiment.pickle', 'rb'))

        self.ui.horizontalSlider.setMinimum(0)
        self.ui.horizontalSlider.setMaximum(self.experiment.raw_image_data.shape[0] - 1)

        self.ui.horizontalSlider.valueChanged.connect(self.handle_slider_changed)

        # Set up images
        self.rot_image_grid = ImageGridWidget(self.experiment.rot_fl)
        self.raw_image_grid = ImageGridWidget(self.experiment.fl_images)

        self.ui.rotatedImagesBox.layout().addWidget(self.rot_image_grid)
        self.ui.rawImagesBox.layout().addWidget(self.raw_image_grid)

        # Set up plots
        # intensity_plot_widget = PlotWidget(background='w')
        redox_plot_widget = PlotWidget(background='w')

        self.intensity_plot_widget = ProfilePlotGridWidget(self.experiment.raw_intensity_data)
        self.ui.intensityPlotBox.layout().addWidget(self.intensity_plot_widget)
        self.ui.redoxPlotBox.layout().addWidget(redox_plot_widget)

        # Event Handlers
        self.ui.horizontalSlider.valueChanged.connect(self.handle_slider_changed)

    def handle_slider_changed(self):
        self.frame = int(self.ui.horizontalSlider.value())
        self.ui.label.setText(str(self.frame))
        self.rot_image_grid.set_frame(self.frame)
        self.raw_image_grid.set_frame(self.frame)
        self.intensity_plot_widget.set_frame(self.frame)


if __name__ == '__main__':
    pg.setConfigOption('imageAxisOrder', 'row-major')
    pg.setConfigOptions(antialias=True)
    qapp = QtWidgets.QApplication([])
    qapp.setAttribute(QtCore.Qt.AA_EnableHighDpiScaling)
    window = MainWindow()
    window.show()
    qapp.exec_()

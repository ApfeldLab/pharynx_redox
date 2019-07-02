import pickle

import numpy as np
import pyqtgraph as pg
from PyQt5 import QtWidgets
from pyqtgraph import GraphicsLayoutWidget, ImageItem, PlotWidget

import pharynx_analysis.pharynx_io as pio
from gui.gui_pyqtgraph import Ui_MainWindow
from pharynx_analysis import experiment


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
            self.image_items[wvl].setImage(self.image_stack_set.sel(wavelength=wvl)[self.frame].data,
                                           autoDownsample=True)


class ProfilePlotGridWidget(GraphicsLayoutWidget):
    def __init__(self, profile_data, **kwargs):
        super(ProfilePlotGridWidget, self).__init__(**kwargs)
        self.frame = 0
        self.profile_data = profile_data
        self.wavelengths = self.profile_data.wavelength.data
        self.xs = np.linspace(1, 100, self.profile_data.shape[2])
        self.plots = {}
        self.means = {}
        self.current_strain = self.profile_data.strain.data[self.frame]

        for wvl in self.wavelengths:
            self.means[wvl] = {}
            for strain in self.profile_data.strain.data:
                self.means[wvl][strain] = self.profile_data.sel(wavelength=wvl, strain=strain).mean(dim='strain')
        self.legends = {}

        self.idx_plot = {}
        self.mean_plot = {}

        for i, wvl in enumerate(self.wavelengths):
            if (i > 0) and (i % 2 == 0):
                self.nextRow()
            self.plots[wvl] = self.addPlot(title=wvl)
            self.plots[wvl].setYRange(0, 1.55e4)
            self.plots[wvl].setXRange(0, 100)
            self.plots[wvl].disableAutoRange()
            self.idx_plot[wvl] = self.plots[wvl].plot(x=self.xs,
                                                      y=self.profile_data.sel(wavelength=wvl)[self.frame].data)
            self.mean_plot[wvl] = self.plots[wvl].plot(x=self.xs, y=self.means[wvl][self.current_strain].data,
                                                       pen={'color': 'b'})

        self.set_frame(self.frame)

    def set_frame(self, frame):
        self.frame = frame
        for wvl in self.wavelengths:
            strain = self.profile_data.strain.data[self.frame]
            if strain is not self.current_strain:
                self.mean_plot[wvl].setData(x=self.xs, y=self.means[wvl][strain].data, pen={'color': 'b'})
                self.current_strain = strain
            self.idx_plot[wvl].setData(x=self.xs, y=self.profile_data.sel(wavelength=wvl)[self.frame].data)

    def _mean_wvl_by_strain(self, wvl):
        return self.profile_data.sel(wavelength=wvl).groupby('strain', restore_coord_dims=False).mean(dim='strain')


class MainWindow(QtWidgets.QMainWindow):

    def __init__(self, reload=False):
        super(MainWindow, self).__init__()

        self.ui = Ui_MainWindow()
        self.ui.setupUi(self)

        self.frame = 0

        if reload:
            img_path = "/Users/sean/code/wormAnalysis/data/paired_ratio_movement_data_sean/2017_02_22-HD233_SAY47/2017_02_22-HD233_SAY47.tif"
            strain_map_path = "/Users/sean/code/wormAnalysis/data/paired_ratio_movement_data_sean/2017_02_22-HD233_SAY47/indexer.csv"

            strains = pio.load_strain_map(strain_map_path)
            self.experiment = experiment.PairExperiment(img_path, "TL/470_1/410_1/470_2/410_2", strains)
            pickle.dump(self.experiment, open('/Users/sean/code/wormAnalysis/data/experiment.pickle', 'wb'))
        else:
            self.experiment = pickle.load(open('/Users/sean/code/wormAnalysis/data/experiment.pickle', 'rb'))

        self.ui.horizontalSlider.setMinimum(0)
        self.ui.horizontalSlider.setMaximum(self.experiment.raw_image_data.shape[0] - 1)

        self.ui.horizontalSlider.valueChanged.connect(self.handle_slider_changed)

        # Set up images
        self.rot_image_grid = ImageGridWidget(self.experiment.rot_fl)
        # self.raw_image_grid = ImageGridWidget(self.experiment.fl_images)

        self.ui.rotatedImagesBox.layout().addWidget(self.rot_image_grid)
        # self.ui.rawImagesBox.layout().addWidget(self.raw_image_grid)

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
        # self.raw_image_grid.set_frame(self.frame)
        self.intensity_plot_widget.set_frame(self.frame)


if __name__ == '__main__':
    pg.setConfigOption('imageAxisOrder', 'row-major')
    qapp = QtWidgets.QApplication([])
    window = MainWindow(reload=True)
    window.show()
    qapp.exec_()

import pickle

import numpy as np
import pyqtgraph as pg
from PyQt5 import QtWidgets, QtGui, QtCore
from pyqtgraph import GraphicsLayoutWidget, ImageItem, PlotWidget, ViewBox

import pharynx_analysis.pharynx_io as pio
from gui.gui_pyqtgraph import Ui_MainWindow
from pharynx_analysis import experiment


class ImageGridWidget(GraphicsLayoutWidget):
    # TODO Docs

    def __init__(self, fl_images, ratio_images, image_display_order, midlines=None, **kwargs):
        super(ImageGridWidget, self).__init__(**kwargs)
        self.image_stack_set = fl_images
        self.ratio_images = ratio_images
        self.wavelengths = fl_images.wavelength.data
        self.image_display_order = image_display_order
        self.frame = 0
        self.midlines = midlines
        self.midline_plots = {}
        self.n_animals = self.image_stack_set.strain.size

        self.image_items = {
            wvl: ImageItem(image=fl_images.sel(wavelength=wvl)[self.frame].data, border='w')
            for wvl in self.wavelengths
        }

        for i in [0, 1]:
            self.image_items[f'r{i + 1}'] = ImageItem(image=ratio_images.isel(pair=i)[self.frame].data, border='w')

        self.image_viewboxes = {}
        for i, wvl in enumerate(image_display_order):
            if (i > 0) and (i % 3 == 0):
                self.nextRow()
            self.image_viewboxes[wvl] = self.addViewBox(name=wvl)
            self.image_viewboxes[wvl].addItem(self.image_items[wvl])
            self.image_viewboxes[wvl].setAspectLocked(True)
            if i > 0:
                self.image_viewboxes[wvl].linkView(ViewBox.XAxis, self.image_viewboxes[self.wavelengths[0]])
                self.image_viewboxes[wvl].linkView(ViewBox.YAxis, self.image_viewboxes[self.wavelengths[0]])

        # Draw Midlines
        if self.midlines:
            self.midline_xs = {
                wvl: np.tile(np.linspace(50, 120), (self.n_animals, 1))
                for wvl in self.wavelengths
            }
            self.midline_ys = {
                wvl: np.array(
                    [self.midlines[wvl][frame](self.midline_xs[wvl][frame]) for frame in range(self.n_animals)])
                for wvl in self.wavelengths
            }
            for wvl in self.wavelengths:
                self.midline_plots[wvl] = pg.PlotDataItem(pen={'color': 'r', 'width': 2})
                self.image_viewboxes[wvl].addItem(self.midline_plots[wvl])

    def update_midlines(self):
        for wvl in self.wavelengths:
            self.midline_plots[wvl].setData(x=self.midline_xs[wvl][self.frame], y=self.midline_ys[wvl][self.frame])

    def update_images(self):
        for wvl in self.wavelengths:
            self.image_items[wvl].setImage(self.image_stack_set.sel(wavelength=wvl)[self.frame].data)
        for i in [0, 1]:
            wvl = f'r{i + 1}'
            self.image_items[wvl].setImage(self.ratio_images.isel(pair=i)[self.frame].data)

    def set_frame(self, frame):
        self.frame = frame
        self.update_images()
        if self.midlines:
            self.update_midlines()


class ProfilePlotGridWidget(GraphicsLayoutWidget):
    def __init__(self, profile_data, regions, **kwargs):
        super(ProfilePlotGridWidget, self).__init__(**kwargs)
        self.frame = 0
        self.profile_data = profile_data
        self.regions = regions
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
        self.linear_regions = {}
        for i, wvl in enumerate(self.wavelengths):
            if (i > 0) and (i % 2 == 0):
                self.nextRow()
            self.plots[wvl] = self.addPlot(title=wvl, background='w')
            self.plots[wvl].setYRange(2e3, 1.55e4)
            self.plots[wvl].setXRange(0, 100)
            self.plots[wvl].disableAutoRange()
            self.idx_plot[wvl] = self.plots[wvl].plot(x=self.xs,
                                                      y=self.profile_data.sel(wavelength=wvl)[self.frame].data)
            self.mean_plot[wvl] = self.plots[wvl].plot(x=self.xs, y=self.means[wvl][self.current_strain].data,
                                                       pen={'color': 'r'})
            self.linear_regions[wvl] = pg.LinearRegionItem()

        # TODO address region plot boundaries
        # This code *works*, but looks messy on the screen... maybe just put the boundaries on one of the subplots?
        # Maybe have a flag in the UI to enable/disable them?

        # self.region_handles = {}
        # for wvl in self.wavelengths:
        #     self.region_handles[wvl] = {}
        #     for i, (region_name, region_boundaries) in enumerate(self.regions.items()):
        #         self.region_handles[wvl][region_name] = pg.LinearRegionItem(values=self.regions[region_name], brush=pg.intColor(i, alpha=50))
        #         self.plots[wvl].addItem(self.region_handles[wvl][region_name])

        self.set_frame(self.frame)

    def set_frame(self, frame):
        self.frame = frame
        for wvl in self.wavelengths:
            strain = self.profile_data.strain.data[self.frame]
            if strain is not self.current_strain:
                self.mean_plot[wvl].setData(x=self.xs, y=self.means[wvl][strain].data, pen={'color': 'r'})
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
        self.rot_image_grid = ImageGridWidget(self.experiment.rot_fl, midlines=self.experiment.midlines,
                                              ratio_images=self.experiment.rot_ratios,
                                              image_display_order=self.experiment.image_display_order)
        self.raw_image_grid = ImageGridWidget(self.experiment.fl_images,
                                              ratio_images=self.experiment.rot_ratios,
                                              image_display_order=self.experiment.image_display_order)

        self.ui.rotatedImagesBox.layout().addWidget(self.rot_image_grid)
        self.ui.rawImagesBox.layout().addWidget(self.raw_image_grid)

        # Set up plots
        redox_plot_widget = PlotWidget(background='w')

        self.intensity_plot_widget = ProfilePlotGridWidget(
            self.experiment.trimmed_intensity_data, self.experiment.get_scaled_region_boundaries()
        )
        self.ui.intensityPlotBox.layout().addWidget(self.intensity_plot_widget)
        self.ui.redoxPlotBox.layout().addWidget(redox_plot_widget)

        # Event Handlers
        self.ui.horizontalSlider.valueChanged.connect(self.handle_slider_changed)

    def keyPressEvent(self, a0: QtGui.QKeyEvent) -> None:
        k = a0.key()
        if k == QtCore.Qt.Key_Right:
            self.set_frame(self.frame + 1)

        if k == QtCore.Qt.Key_Left:
            self.set_frame(self.frame - 1)

    def set_frame(self, new_frame):
        self.frame = max(0, min(self.experiment.raw_image_data.strain.size - 1, new_frame))
        self.ui.label.setText(str(self.frame))
        self.rot_image_grid.set_frame(self.frame)
        self.raw_image_grid.set_frame(self.frame)
        self.intensity_plot_widget.set_frame(self.frame)
        self.ui.horizontalSlider.setValue(self.frame)

    def handle_slider_changed(self):
        self.set_frame(int(self.ui.horizontalSlider.value()))


if __name__ == '__main__':
    pg.setConfigOptions(imageAxisOrder='row-major', antialias=True)
    qapp = QtWidgets.QApplication([])
    window = MainWindow(reload=True)
    window.show()
    qapp.exec_()

import os
import pickle
import sys

import matplotlib.projections as proj
import matplotlib.pyplot as plt
import numpy as np
from PyQt5 import QtWidgets, QtCore
from bidict import bidict
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.backends.backend_qt5agg import NavigationToolbar2QT as NavigationToolbar
from matplotlib.figure import Figure
from matplotlib.lines import Line2D
from scipy.interpolate import UnivariateSpline

import pharynx_analysis.pharynx_io as pio
from gui import plots
from gui.gui_ui import Ui_MainWindow

config_dict = {
    'default_imaging_scheme': 'TL/470_1/410_1/470_2/410_2'
}

img_path = "/Users/sean/code/wormAnalysis/data/paired_ratio_movement_data_sean/2017_02_22-HD233_SAY47/2017_02_22-HD233_SAY47.tif"
strain_map_path = "/Users/sean/code/wormAnalysis/data/paired_ratio_movement_data_sean/2017_02_22-HD233_SAY47/indexer.csv"

strains = pio.load_strain_map(strain_map_path)

dir_path = os.path.dirname(os.path.realpath(__file__))


class MatplotlibWidget(QtWidgets.QWidget):
    """
    Implements a Matplotlib figure inside a QWidget.
    Use getFigure() and redraw() to interact with matplotlib.

    Example::

        mw = MatplotlibWidget()
        subplot = mw.getFigure().add_subplot(111)
        subplot.plot(x,y)
        mw.draw()
    """

    def __init__(self, size=(5.0, 4.0), dpi=100):
        QtWidgets.QWidget.__init__(self)
        self.fig = Figure(size, dpi=dpi)
        self.canvas = FigureCanvas(self.fig)
        self.canvas.setParent(self)
        self.toolbar = NavigationToolbar(self.canvas, self)

        self.vbox = QtWidgets.QVBoxLayout()
        self.vbox.addWidget(self.toolbar, alignment=QtCore.Qt.AlignTop)
        self.vbox.addWidget(self.canvas, alignment=QtCore.Qt.AlignCenter)

        self.setLayout(self.vbox)

    def get_figure(self):
        return self.fig

    def draw(self):
        self.canvas.draw()


class ImageGroupLayout(QtWidgets.QGridLayout):

    def __init__(self, images, midlines=None, n_spline_control_points=5, crop=False):
        super(ImageGroupLayout, self).__init__()
        proj.register_projection(plots.ThinAxes)

        self.frame = 0
        self.images = images
        self.n_animals = self.images.shape[1]

        self.active_vertex_idx = None
        self.epsilon = 3

        self.main_widget = MatplotlibWidget()
        self.addWidget(self.main_widget)
        self.fig = self.main_widget.get_figure()
        self.canvas = self.fig.canvas

        self.wavelengths = self.images.wavelength.data
        self.backgrounds = {}

        self.midlines = midlines
        if self.midlines:
            self.spline_resolution = 100
            self.spline_xs = {
                wvl: np.tile(np.linspace(50, 120, self.spline_resolution), (self.n_animals, 1))
                for wvl in self.wavelengths
            }
            self.spline_ys = {
                wvl: np.array(
                    [self.midlines[wvl][frame](self.spline_xs[wvl][frame]) for frame in range(self.n_animals)])
                for wvl in self.wavelengths
            }
            self.n_spline_control_points = n_spline_control_points

            self.spline_control_xs = {
                wvl: np.tile(
                    np.linspace(np.min(self.spline_xs[wvl]), np.max(self.spline_xs[wvl]), n_spline_control_points),
                    (self.n_animals, 1)
                ) for wvl in self.wavelengths
            }
            self.spline_control_ys = {
                wvl: [self.midlines[wvl][frame](self.spline_control_xs[wvl][frame]) for frame in range(self.n_animals)]
                for wvl in self.wavelengths
            }
            self.pts = {
                wvl: Line2D(
                    self.spline_control_xs[wvl][self.frame],
                    self.midlines[wvl][self.frame](self.spline_control_xs[wvl][self.frame]),
                    color='white', ls='', marker='o', markersize=self.epsilon, animated=True, zorder=100)
                for wvl in self.wavelengths
            }
            self.lines = {
                wvl: Line2D(
                    self.spline_xs[wvl][self.frame],
                    self.spline_ys[wvl][self.frame],
                    color='red', animated=True, alpha=0.3, linewidth=2)
                for wvl in self.wavelengths
            }

        # Set up axes
        xlim = (45, 125)
        ylim = (75, 50)

        if crop:
            kwargs = {'xlim': xlim, 'ylim': ylim}
        else:
            kwargs = {}
        self.axes = bidict({
            '410_1': self.main_widget.get_figure().add_subplot(321, title="$I_{410_1}$", **kwargs),
            '410_2': self.main_widget.get_figure().add_subplot(322, title="$I_{410_2}$", **kwargs),
            '470_1': self.main_widget.get_figure().add_subplot(323, title="$I_{470_1}$", **kwargs),
            '470_2': self.main_widget.get_figure().add_subplot(324, title="$I_{470_2}$", **kwargs),
            'R1': self.main_widget.get_figure().add_subplot(325, title="$I_{R_1}$", **kwargs),
            'R2': self.main_widget.get_figure().add_subplot(326, title="$I_{R_2}$", **kwargs)
        })

        for i, ax in enumerate(self.axes.values()):
            ax.get_shared_x_axes().join(ax, self.axes['410_1'])
            ax.get_shared_y_axes().join(ax, self.axes['410_1'])

        if midlines:
            for wvl in self.wavelengths:
                self.axes[wvl].add_artist(self.pts[wvl])
                self.axes[wvl].add_artist(self.lines[wvl])

        self.ims = {
            wvl: self.axes[wvl].imshow(self.images.sel(wavelength=wvl)[self.frame]) for wvl in self.wavelengths
        }

        self.ims['R1'] = self.axes['R1'].imshow(
            (self.images.sel(wavelength='410_1') / self.images.sel(wavelength='470_1'))[self.frame]
        )
        self.ims['R2'] = self.axes['R2'].imshow(
            (self.images.sel(wavelength='410_2') / self.images.sel(wavelength='470_2'))[self.frame]
        )

        for wvl in self.wavelengths:
            self.ims[wvl].set_clim(0, 1.5e4)
        for ax in self.axes.values():
            ax.set_axis_off()

        # self.fig.subplots_adjust(right=0.8)
        # self.cbar_ax = self.fig.add_axes([0.85, .10, .05, .7])
        # self.fig.colorbar(self.ims['410_1'], cax=self.cbar_ax)

        self.background = None

        self.canvas.mpl_connect('draw_event', self.draw_callback)
        self.canvas.mpl_connect('button_press_event', self.button_press_callback)
        self.canvas.mpl_connect('button_release_event', self.button_release_callback)
        self.canvas.mpl_connect('motion_notify_event', self.motion_notify_callback)

        plt.tight_layout()

    def button_press_callback(self, event):
        """Called when a mouse button is pressed"""
        if self.midlines is None:
            return
        if event.inaxes is None:
            return
        if event.button != 1:  # check that it is left mouse
            return
        self.active_vertex_idx = self.get_ind_under_point(event)

    def button_release_callback(self, event):
        """Called when a mouse button is released"""
        if event.button != 1:
            return
        self.active_vertex_idx = None

    def motion_notify_callback(self, event):
        if self.active_vertex_idx is None:
            return
        if self.midlines is None:
            return
        ax = event.inaxes
        if ax is None:
            return

        wvl = self.axes.inverse[ax]

        try:
            midline = self.midlines[wvl][self.frame]
        except KeyError:
            return

        mouse_x, mouse_y = event.xdata, event.ydata

        old_ctrl_xs = self.spline_control_xs[wvl][self.frame]
        old_ctrl_ys = self.spline_control_ys[wvl][self.frame]

        self.spline_control_xs[wvl][self.frame][self.active_vertex_idx] = mouse_x
        self.spline_control_ys[wvl][self.frame][self.active_vertex_idx] = mouse_y

        ctrl_xs = self.spline_control_xs[wvl][self.frame]
        ctrl_ys = self.spline_control_ys[wvl][self.frame]

        self.midlines[wvl][self.frame] = self.make_spline(ctrl_xs, ctrl_ys)

        self.spline_xs[wvl][self.frame] = np.linspace(np.min(ctrl_xs), np.max(ctrl_xs), self.spline_resolution)
        self.spline_ys[wvl][self.frame] = midline(self.spline_xs[wvl][self.frame])

        self.lines[wvl].set_data(self.spline_xs[wvl][self.frame], self.spline_ys[wvl][self.frame])
        self.pts[wvl].set_data(ctrl_xs, ctrl_ys)

        self.canvas.restore_region(self.backgrounds[wvl])
        ax.draw_artist(self.lines[wvl])
        ax.draw_artist(self.pts[wvl])
        self.canvas.draw()

    def make_spline(self, xs, ys):
        return UnivariateSpline(xs, ys, ext=0, s=0, k=3)

    def get_ind_under_point(self, event):
        ax = event.inaxes
        wvl = self.axes.inverse[ax]

        try:
            midline = self.midlines[wvl][self.frame]
        except KeyError:
            return None
        ctrl_xs = self.spline_control_xs[wvl][self.frame]

        d = np.hypot(ctrl_xs - event.xdata, midline(ctrl_xs) - event.ydata)
        index_sequence, = np.nonzero(d == d.min())
        ind = index_sequence[0]
        if d[ind] >= self.epsilon:
            ind = None
        return ind

    def set_frame(self, frame):
        self.frame = frame

        for wvl in self.wavelengths:
            self.ims[wvl].set_data(self.images.sel(wavelength=wvl)[self.frame])

        self.ims['R1'].set_data((self.images.sel(wavelength='410_1') / self.images.sel(wavelength='470_1'))[self.frame])
        self.ims['R2'].set_data((self.images.sel(wavelength='410_2') / self.images.sel(wavelength='470_2'))[self.frame])

        if self.midlines:
            for wvl in self.wavelengths:
                midline = self.midlines[wvl][self.frame]
                ctrl_xs = self.spline_control_xs[wvl][self.frame]
                ctrl_ys = midline(ctrl_xs)
                self.pts[wvl].set_data(ctrl_xs, ctrl_ys)
                self.lines[wvl].set_data(self.spline_xs[wvl][self.frame], self.spline_ys[wvl][self.frame])

        self.canvas.draw()

    def draw_callback(self, e):
        self.backgrounds = {
            wvl: self.canvas.copy_from_bbox(self.axes[wvl].bbox) for wvl in self.wavelengths
        }
        for wvl, ax in self.axes.items():
            self.axes[wvl].draw_artist(self.ims[wvl])
        if self.midlines:
            for wvl in self.wavelengths:
                self.axes[wvl].draw_artist(self.lines[wvl])
                self.axes[wvl].draw_artist(self.pts[wvl])


class MPLLayout(QtWidgets.QVBoxLayout):

    def __init__(self, data):
        super(MPLLayout, self).__init__()
        self.data = data

        self.frame = 0
        self.strains = data.strain.data

        self.main_widget = MatplotlibWidget()
        self.addWidget(self.main_widget)
        self.fig = self.main_widget.get_figure()
        self.wavelengths = data.wavelength.data

        if len(self.wavelengths) > 2:
            n_cols = 2
        else:
            n_cols = 1

        n_rows = len(self.wavelengths) // n_cols
        n_rows += len(self.wavelengths) % n_cols

        pos = range(1, len(self.wavelengths) + 1)

        self.axes = bidict({
            wvl: self.main_widget.get_figure().add_subplot(n_rows, n_cols, pos[i])
            for i, wvl in enumerate(self.wavelengths)
        })

        self.xs = np.linspace(0, 500, data.shape[2])
        self.individual_measurements = {
            wvl: Line2D(self.xs, self.data.sel(wavelength=wvl)[self.frame], linestyle='--', color='k')
            for wvl in self.wavelengths
        }

        for i, wvl in enumerate(self.data.wavelength.data):
            ax = self.axes[wvl]
            d = self.data.sel(wavelength=wvl).groupby('strain').mean(dim='strain').T
            ax.plot(d)
            ax.legend(d.strain.data)
            ax.set_ylim((0, 1.6e4))
            ax.set_title(f'${wvl}$')
            ax.add_artist(self.individual_measurements[wvl])

        plt.tight_layout()

    def set_frame(self, frame):
        self.frame = frame

        backgrounds = {
            wvl: self.fig.canvas.copy_from_bbox(self.axes[wvl].bbox) for wvl in self.wavelengths
        }

        for wvl in self.wavelengths:
            ax = self.axes[wvl]
            self.individual_measurements[wvl].set_data(self.xs, self.data.sel(wavelength=wvl)[self.frame])
            self.fig.canvas.restore_region(backgrounds[wvl])
            ax.draw_artist(self.individual_measurements[wvl])

        self.fig.canvas.draw()


class MainWindow(QtWidgets.QMainWindow):

    def __init__(self):
        super(MainWindow, self).__init__()

        self.ui = Ui_MainWindow()
        self.ui.setupUi(self)

        self.frame = 0

        # self.experiment = experiment.PairExperiment(img_path, "TL/470_1/410_1/470_2/410_2", strains)
        self.experiment = pickle.load(open('/Users/sean/code/wormAnalysis/data/experiment.pickle', 'rb'))

        self.ui.rotatedImagesTab.setLayout(
            ImageGroupLayout(self.experiment.rot_fl, midlines=self.experiment.midlines, crop=True))
        self.ui.rawImagesTab.setLayout(ImageGroupLayout(self.experiment.fl_images))

        self.ui.intensityPlotTab.setLayout(MPLLayout(self.experiment.raw_intensity_data))
        self.ui.redoxPlotTab.setLayout(MPLLayout(self.experiment.raw_intensity_data))

        self.ui.horizontalSlider.setMaximum(self.experiment.raw_image_data.shape[0] - 1)

        self.ui.flipLRpushButton.clicked.connect(self.handle_flip_lr)
        self.ui.horizontalSlider.valueChanged.connect(self.handle_slider_changed)
        self.ui.imagesTabWidget.currentChanged.connect(self.handle_images_tab_change)

    def handle_flip_lr(self):
        print(f'flip lr was pressed. frame{self.frame}')
        self.experiment.flip_at(self.frame)

    def handle_slider_changed(self):
        self.frame = int(self.ui.horizontalSlider.value())

        self.ui.label.setText(str(self.frame))

        self.ui.imagesTabWidget.currentWidget().layout().set_frame(self.frame)
        self.ui.intensityPlotTab.layout().set_frame(self.frame)

    def handle_images_tab_change(self):
        self.ui.imagesTabWidget.currentWidget().layout().set_frame(self.frame)


if __name__ == '__main__':
    qapp = QtWidgets.QApplication([])
    qapp.setAttribute(QtCore.Qt.AA_EnableHighDpiScaling)
    window = MainWindow()
    window.show()

    sys.exit(qapp.exec_())

    # import pickle
    #
    # pe = experiment.PairExperiment(img_path, "TL/470_1/410_1/470_2/410_2", strains)
    # pickle.dump(pe, open('/Users/sean/Desktop/experiment.pickle', 'wb'))

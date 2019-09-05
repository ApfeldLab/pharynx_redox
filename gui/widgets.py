import numpy as np
import pyqtgraph as pg
from PyQt5 import QtCore
from pyqtgraph import GraphicsLayoutWidget, ImageItem, ViewBox


class ImageGridWidget(GraphicsLayoutWidget):
    # TODO Docs

    def __init__(
        self, fl_images, ratio_images, image_display_order, midlines=None, **kwargs
    ):
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
            wvl: ImageItem(
                image=fl_images.sel(wavelength=wvl)[self.frame].data, border="w"
            )
            for wvl in self.wavelengths
        }

        for i in [0, 1]:
            self.image_items[f"r{i + 1}"] = ImageItem(
                image=ratio_images.isel(pair=i)[self.frame].data, border="w"
            )

        self.image_viewboxes = {}
        for i, wvl in enumerate(image_display_order):
            if (i > 0) and (i % 3 == 0):
                self.nextRow()
            self.image_viewboxes[wvl] = self.addViewBox(name=wvl)
            self.image_viewboxes[wvl].addItem(self.image_items[wvl])
            self.image_viewboxes[wvl].setAspectLocked(True)

            # add label
            # TODO: make image legends nicer
            wvl_label = pg.TextItem(text=f"{wvl}")
            self.image_viewboxes[wvl].addItem(wvl_label)

            if i > 0:
                self.image_viewboxes[wvl].linkView(
                    ViewBox.XAxis, self.image_viewboxes[self.wavelengths[0]]
                )
                self.image_viewboxes[wvl].linkView(
                    ViewBox.YAxis, self.image_viewboxes[self.wavelengths[0]]
                )

        # Draw Midlines
        if self.midlines:
            self.midline_xs = {
                wvl: np.tile(np.linspace(50, 120), (self.n_animals, 1))
                for wvl in self.wavelengths
            }
            self.midline_ys = {
                wvl: np.array(
                    [
                        self.midlines[wvl][frame](self.midline_xs[wvl][frame])
                        for frame in range(self.n_animals)
                    ]
                )
                for wvl in self.wavelengths
            }
            for wvl in self.wavelengths:
                self.midline_plots[wvl] = pg.PlotDataItem(
                    pen={"color": "r", "width": 2}
                )
                self.image_viewboxes[wvl].addItem(self.midline_plots[wvl])

            self.ratio_midline_xs = {
                i: {
                    wvl: self.midline_xs[wvl].copy()
                    for wvl in [f"410_{i + 1}", f"470_{i + 1}"]
                }
                for i in [0, 1]
            }
            self.ratio_midline_ys = {
                i: {
                    wvl: self.midline_ys[wvl].copy()
                    for wvl in [f"410_{i + 1}", f"470_{i + 1}"]
                }
                for i in [0, 1]
            }

            for i in [0, 1]:
                for j, wvl in enumerate(self.ratio_midline_xs[i].keys()):
                    self.midline_plots[f"r{i + 1}_{wvl}"] = pg.PlotDataItem(
                        pen=pg.intColor(j)
                    )
                    self.image_viewboxes[f"r{i + 1}"].addItem(
                        self.midline_plots[f"r{i + 1}_{wvl}"]
                    )

        self.set_frame(self.frame)

    def update_midlines(self):
        for wvl in self.wavelengths:
            self.midline_plots[wvl].setData(
                x=self.midline_xs[wvl][self.frame], y=self.midline_ys[wvl][self.frame]
            )

        for i in [0, 1]:
            for j, wvl in enumerate(self.ratio_midline_xs[i].keys()):
                self.midline_plots[f"r{i + 1}_{wvl}"].setData(
                    x=self.ratio_midline_xs[i][wvl][self.frame],
                    y=self.ratio_midline_ys[i][wvl][self.frame],
                )

    def update_images(self):
        for wvl in self.wavelengths:
            self.image_items[wvl].setImage(
                self.image_stack_set.sel(wavelength=wvl)[self.frame].data
            )
        for i in [0, 1]:
            wvl = f"r{i + 1}"
            self.image_items[wvl].setImage(
                self.ratio_images.isel(pair=i)[self.frame].data
            )

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
                self.means[wvl][strain] = self.profile_data.sel(
                    wavelength=wvl, strain=strain
                ).mean(dim="strain")
        self.legends = {}

        self.idx_plot = {}
        self.mean_plot = {}
        self.linear_regions = {}
        for i, wvl in enumerate(self.wavelengths):
            if (i > 0) and (i % 2 == 0):
                self.nextRow()
            self.plots[wvl] = self.addPlot(title=wvl, background="w")
            self.plots[wvl].setYRange(2e3, 1.55e4)
            self.plots[wvl].setXRange(0, 100)
            self.plots[wvl].disableAutoRange()
            self.idx_plot[wvl] = self.plots[wvl].plot(
                x=self.xs, y=self.profile_data.sel(wavelength=wvl)[self.frame].data
            )
            self.mean_plot[wvl] = self.plots[wvl].plot(
                x=self.xs,
                y=self.means[wvl][self.current_strain].data,
                pen={"color": "r"},
            )
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
                self.mean_plot[wvl].setData(
                    x=self.xs, y=self.means[wvl][strain].data, pen={"color": "r"}
                )
                self.current_strain = strain
            self.idx_plot[wvl].setData(
                x=self.xs, y=self.profile_data.sel(wavelength=wvl)[self.frame].data
            )

    def _mean_wvl_by_strain(self, wvl):
        return (
            self.profile_data.sel(wavelength=wvl)
            .groupby("strain", restore_coord_dims=False)
            .mean(dim="strain")
        )


class PandasModel(QtCore.QAbstractTableModel):
    """
    Class to populate a table view with a pandas dataframe
    """

    def __init__(self, data, parent=None):
        QtCore.QAbstractTableModel.__init__(self, parent)
        self._data = data

    def rowCount(self, parent=None):
        return len(self._data.values)

    def columnCount(self, parent=None):
        return self._data.columns.size

    def data(self, index, role=QtCore.Qt.DisplayRole):
        if index.isValid():
            if role == QtCore.Qt.DisplayRole:
                return str(self._data.values[index.row()][index.column()])
        return None

    def headerData(self, col, orientation, role):
        if orientation == QtCore.Qt.Horizontal and role == QtCore.Qt.DisplayRole:
            return self._data.columns[col]
        return None

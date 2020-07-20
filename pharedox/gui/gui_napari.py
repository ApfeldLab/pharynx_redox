"""
A Basic GUI based on napari
"""
import logging
import sys
from pathlib import Path
from typing import Optional

import napari
import numpy as np
import xarray as xr
from PyQt5.QtCore import pyqtSignal
from qtpy.QtWidgets import QMessageBox, QWidget
from skimage import morphology
from skimage.measure import label

from pharedox import experiment, plots, utils
from pharedox.gui.qt_py_files.pipeline_buttons import Ui_Form


def segment_pharynxes(imgs, t, wvl):
    segs = imgs.sel(wavelength=wvl) > t

    segs = xr.apply_ufunc(
        lambda x: label(x),
        segs,
        vectorize=True,
        input_core_dims=[["y", "x"]],
        output_core_dims=[["y", "x"]],
    )
    return segs


def remove_small_objects(label_data, min_obj_size):
    return xr.apply_ufunc(
        lambda x: morphology.remove_small_objects(x, min_obj_size),
        label_data,
        vectorize=True,
        input_core_dims=[["y", "x"]],
        output_core_dims=[["y", "x"]],
    )


class PipelineButtonsWidget(QWidget):

    t_slider_changed = pyqtSignal()

    def __init__(self, *args, **kwargs):
        super(PipelineButtonsWidget, self).__init__(*args, **kwargs)

        self.ui = Ui_Form()
        self.ui.setupUi(self)

        self.ui.thresholdSlider.setMinimum(np.iinfo(np.uint16).min)
        self.ui.thresholdSlider.setMaximum(np.iinfo(np.uint16).max)

        self.ui.thresholdSlider.valueChanged.connect(self.handle_t_slider_changed)
        self.ui.thresholdSpinBox.valueChanged.connect(
            self.handle_threshold_spin_box_changed
        )

    def handle_t_slider_changed(self):
        self.ui.thresholdSpinBox.setValue(self.ui.thresholdSlider.value())
        self.t_slider_changed.emit()

    def handle_threshold_spin_box_changed(self):
        self.ui.thresholdSlider.setValue(self.ui.thresholdSpinBox.value())


class App:

    viewer = None
    buttons = None

    def __init__(self, exp_):
        self.experiment = exp_

    def set_up_viewer(self):
        self.viewer = napari.Viewer()
        self.buttons = PipelineButtonsWidget()
        self.viewer.window.add_dock_widget(self.buttons, name="pipeline", area="left")

        # connect signals/slots
        self.buttons.t_slider_changed.connect(self.handle_t_slider_changed)
        self.buttons.ui.removeObjectsButton.pressed.connect(
            self.handle_remove_objects_pressed
        )
        self.buttons.ui.runNeuronsButton.pressed.connect(self.run_neuron_analysis)
        self.buttons.ui.runPharynxButton.pressed.connect(self.run_pharynx_analysis)

    def run_pharynx_analysis(self):
        if self.experiment.seg_images is None:
            self.show_simple_dialog("no masks")
        else:
            self.experiment.full_pipeline()
            self.show_dialog("Analysis finished!")

    def run_neuron_analysis(self):
        if self.experiment.seg_images is not None:
            self.experiment.run_neuron_pipeline()
            self.show_dialog("Analysis finished!")

    def show_dialog(self, message, title=""):
        msg_box = QMessageBox()
        msg_box.setIcon(QMessageBox.Information)
        msg_box.setText(message)
        msg_box.setWindowTitle(title)
        msg_box.setStandardButtons(QMessageBox.Open | QMessageBox.Ok)

        return_value = msg_box.exec()
        if return_value == QMessageBox.Open:
            utils.open_folder(self.experiment.analysis_dir)

    @staticmethod
    def show_simple_dialog(message, title=""):
        msg_box = QMessageBox()
        msg_box.setIcon(QMessageBox.Information)
        msg_box.setText(message)
        msg_box.setWindowTitle(title)
        msg_box.setStandardButtons(QMessageBox.Ok)
        msg_box.exec()

    def get_layer(self, name):
        for layer in self.viewer.layers:
            if layer.name == name:
                return layer
        return None

    def handle_remove_objects_pressed(self):
        if self.experiment.seg_images is None:
            return

        layer = self.get_layer("masks")

        min_obj_size = self.buttons.ui.smallObjectSizeSpinBox.value()
        self.experiment.seg_images = remove_small_objects(
            self.experiment.seg_images, min_obj_size
        )

        layer.data = self.experiment.seg_images.values
        layer.refresh()

    def get_current_wvl(self) -> Optional[str]:
        """
        Get the wavelength of the active layer in Napari, or `None` if the active layer
        name does not correspond to a wavelength in the experiment's images.
        """
        wvl_candidate = self.viewer.active_layer.name
        true_wvls = self.experiment.images.wavelength.values
        if wvl_candidate in true_wvls:
            return wvl_candidate
        else:
            return None

    def segment_pharynxes(self, t) -> xr.DataArray:
        wvl = self.get_current_wvl()
        if wvl is None:
            self.show_simple_dialog(
                message="The active layer does not correspond to a wavelength in the data set.",
                title="Invalid Wavelength Selected",
            )
            return

        masks = self.experiment.images.sel(wavelength=wvl) > t

        masks = xr.apply_ufunc(
            lambda x: label(x),
            masks,
            vectorize=True,
            input_core_dims=[["y", "x"]],
            output_core_dims=[["y", "x"]],
        )

        return masks

    def handle_t_slider_changed(self):
        t = self.buttons.ui.thresholdSlider.value()
        self.update_threshold(t)

    def update_threshold(self, t):
        masks = self.segment_pharynxes(t)
        if masks is None:
            return
        self.experiment.seg_images = masks
        self.get_layer("masks").data = masks
        self.get_layer("masks").refresh()

    def run(self):
        with napari.gui_qt():
            self.set_up_viewer()

            if self.experiment.images is not None:
                self.viewer.add_image(
                    plots.imgs_to_rgb(
                        self.experiment.images, r_min=0.9, r_max=1.9, i_max=1500
                    )
                )
                for wvl in self.experiment.images.wavelength.values:
                    self.viewer.add_image(
                        self.experiment.images.sel(wavelength=wvl), name=wvl
                    )

            if self.experiment.seg_images is not None:
                self.viewer.add_labels(
                    self.experiment.seg_images, name="masks",
                )


if __name__ == "__main__":

    logging.basicConfig(
        format="%(asctime)s %(levelname)s:%(message)s",
        level=logging.INFO,
        datefmt="%I:%M:%S",
    )

    exp_dir = sys.argv[1]
    exp = experiment.Experiment(Path(exp_dir))

    app = App(exp_=exp)
    app.run()

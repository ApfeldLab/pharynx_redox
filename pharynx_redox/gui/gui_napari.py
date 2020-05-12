"""
A Basic GUI based on napari
"""
import sys
import numpy as np
import napari
from skimage.measure import label
from skimage import morphology
import xarray as xr

from pathlib import Path
from pharynx_redox import image_processing as ip
from qtpy.QtWidgets import QWidget, QMessageBox
from PyQt5.QtCore import pyqtSignal
from qtpy.QtWidgets import QApplication, QSplashScreen
from pharynx_redox.gui.qt_py_files.pipeline_buttons import Ui_Form
from pharynx_redox import experiment, utils

import logging


def segment_pharynxes(imgs, t, skip_wvl=["TL"], ref_wvl="410"):
    segs = imgs > t
    for wvl in skip_wvl:
        try:
            segs.loc[dict(wavelength=wvl)] = False
        except:
            continue

    segs = xr.apply_ufunc(
        lambda x: label(x),
        segs,
        vectorize=True,
        input_core_dims=[["y", "x"]],
        output_core_dims=[["y", "x"]],
    )
    return segs


def remove_small_objects(label_data, min_obj_size=5):
    return xr.apply_ufunc(
        lambda x: morphology.remove_small_objects(x, min_obj_size),
        label_data,
        vectorize=True,
        input_core_dims=[["y", "x"]],
        output_core_dims=[["y", "x"]],
    )


# napari_hook_implementation = HookimplMarker("napari")
# readable_extensions = tuple(set(x for f in formats for x in f.extensions))


# @napari_hook_implementation
# def napari_get_reader(path):
#     """A basic implementation of the napari_get_reader hook specification."""
#     # if we know we cannot read the file, we immediately return None.
#     if not path.endswith(readable_extensions):
#         return None
#     # otherwise we return the *function* that can read ``path``.
#     return reader_function


# def reader_function(path):
#     """Take a path and returns a list of LayerData tuples."""
#     data = imread(path)
#     # Readers are expected to return data as a list of tuples, where each tuple
#     # is (data, [meta_dict, [layer_type]])
#     return [(data,)]


class PipelineButtonsWidget(QWidget):

    segment_sig = pyqtSignal()
    t_slider_changed = pyqtSignal()

    def __init__(self, *args, **kwargs):
        super(PipelineButtonsWidget, self).__init__(*args, **kwargs)

        self.ui = Ui_Form()
        self.ui.setupUi(self)

        self.ui.thresholdSlider.setMinimum(np.iinfo(np.uint16).min)
        self.ui.thresholdSlider.setMaximum(np.iinfo(np.uint16).max)

        self.ui.segmentButton.pressed.connect(self.segment_sig.emit)

        self.ui.thresholdSlider.valueChanged.connect(self.handle_t_slider_changed)
        self.ui.thresholdSpinBox.valueChanged.connect(
            self.handle_thresholdSpinBox_changed
        )

    def handle_t_slider_changed(self):
        self.ui.thresholdSpinBox.setValue(self.ui.thresholdSlider.value())
        self.t_slider_changed.emit()

    def handle_thresholdSpinBox_changed(self):
        self.ui.thresholdSlider.setValue(self.ui.thresholdSpinBox.value())


class App:

    viewer = None
    buttons = None

    def __init__(self, experiment):
        self.experiment = experiment

    def set_up_viewer(self):
        self.viewer = napari.Viewer()
        self.buttons = PipelineButtonsWidget()
        self.viewer.window.add_dock_widget(self.buttons, name="pipeline", area="left")

        # connect signals/slots
        self.buttons.segment_sig.connect(self.handle_segment_pressed)
        self.buttons.t_slider_changed.connect(self.handle_t_slider_changed)
        self.buttons.ui.removeObjectsButton.pressed.connect(
            self.handle_remove_objects_pressed
        )
        # self.viewer.layers.events.changed.connect(self.on_layers_change)
        self.buttons.ui.runNeuronsButton.pressed.connect(self.run_neuron_analysis)
        self.buttons.ui.runPharynxButton.pressed.connect(self.run_pharynx_analysis)

    def run_pharynx_analysis(self):
        if self.experiment.seg_images is not None:
            self.experiment.full_pipeline()
            self.showDialog("Analysis finished!")

    def run_neuron_analysis(self):
        if self.experiment.seg_images is not None:
            self.experiment.run_neuron_pipeline()
            self.showDialog("Analysis finished!")

    def showDialog(self, message, title=""):
        msg_box = QMessageBox()
        msg_box.setIcon(QMessageBox.Information)
        msg_box.setText(message)
        msg_box.setWindowTitle(title)
        msg_box.setStandardButtons(QMessageBox.Open | QMessageBox.Ok)

        return_value = msg_box.exec()
        if return_value == QMessageBox.Open:
            utils.open_folder(self.experiment.analysis_dir)

    # def on_layers_change(self, event):
    #     if self.get_layer("masks") is None:
    #         self.experiment.seg_images = None

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

        layer.data = self.experiment.seg_images.sel(wavelength="410").values
        layer.refresh()

    def handle_segment_pressed(self):
        t = self.buttons.ui.thresholdSpinBox.value()
        masks = segment_pharynxes(self.experiment.images, t)

        if self.experiment.seg_images is None:
            self.experiment.seg_images = masks
            self.viewer.add_labels(
                self.experiment.seg_images.sel(wavelength="410"), name="masks"
            )
        else:
            self.update_threshold(t)

    def handle_t_slider_changed(self):
        t = self.buttons.ui.thresholdSlider.value()
        self.update_threshold(t)

    def update_threshold(self, t):
        masks = segment_pharynxes(self.experiment.images, t)
        if self.experiment.seg_images is None:
            return
        else:
            self.experiment.seg_images = masks
            self.get_layer("masks").data = masks.sel(wavelength="410")
            self.get_layer("masks").refresh()

    def run(self):
        with napari.gui_qt():
            self.set_up_viewer()

            if self.experiment.images is not None:
                for wvl in self.experiment.images.wavelength.values:
                    self.viewer.add_image(
                        self.experiment.images.sel(wavelength=wvl), name=wvl
                    )
            if self.experiment.seg_images is not None:
                masks = self.experiment.seg_images
                self.viewer.add_labels(
                    self.experiment.seg_images.sel(wavelength="410"), name="masks"
                )


if __name__ == "__main__":

    logging.basicConfig(
        format="%(asctime)s %(levelname)s:%(message)s",
        level=logging.INFO,
        datefmt="%I:%M:%S",
    )

    exp_dir = sys.argv[1]
    exp = experiment.Experiment(Path(exp_dir))

    app = App(experiment=exp)
    app.run()

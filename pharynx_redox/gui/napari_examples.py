import numpy as np
import napari
from skimage import data
from qtpy.QtWidgets import QWidget
from pharynx_redox.gui.qt_py_files.pipeline_buttons import Ui_Form


blobs_raw = data.binary_blobs(length=64, n_dim=2, volume_fraction=0.1)


class PipelineButtonsWidget(QWidget):
    def __init__(self, *args, **kwargs):
        super(PipelineButtonsWidget, self).__init__(*args, **kwargs)

        self.ui = Ui_Form()
        self.ui.setupUi(self)


with napari.gui_qt():
    viewer = napari.view_image(blobs_raw)
    viewer.window.add_dock_widget(PipelineButtonsWidget(), name="pipeline", area="left")

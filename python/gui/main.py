import pickle
from pathlib import Path

import pandas as pd
import pyqtgraph as pg
from PyQt5 import QtWidgets, QtGui, QtCore
from PyQt5.QtWidgets import QFileDialog

import pharynx_analysis.pharynx_io as pio
from gui.qt_py_files.gui_pyqtgraph import Ui_MainWindow
from gui.qt_py_files.load_raw_image import Ui_LoadRawImageFileDialog
from gui.qt_py_files.meta_loader import Ui_MetaLoader
from gui.widgets import ImageGridWidget, ProfilePlotGridWidget
from pharynx_analysis import experiment


class MainWindow(QtWidgets.QMainWindow):

    def __init__(self, reload=False):
        super(MainWindow, self).__init__()
        self.reload = reload

        self.ui = Ui_MainWindow()
        self.ui.setupUi(self)

        self.frame = 0

        # Load Experiment
        if self.reload:
            img_path = "/Users/sean/code/wormAnalysis/data/paired_ratio_movement_data_sean/2017_02_22-HD233_SAY47/2017_02_22-HD233_SAY47.tif"
            strain_map_path = "/Users/sean/code/wormAnalysis/data/paired_ratio_movement_data_sean/2017_02_22-HD233_SAY47/indexer.csv"

            strains = pio.load_strain_map_from_disk(strain_map_path)
            self.experiment = experiment.PairExperiment(img_path, "TL/470_1/410_1/470_2/410_2", strains)
            pickle.dump(self.experiment, open('/Users/sean/code/wormAnalysis/data/experiment.pickle', 'wb'))
        else:
            self.experiment = pickle.load(open('/Users/sean/code/wormAnalysis/data/experiment.pickle', 'rb'))

        self.initialize_slider()
        self.initialize_table()

        # Set up images
        self.rot_image_grid = ImageGridWidget(self.experiment.rot_fl,
                                              midlines=self.experiment.midlines,
                                              ratio_images=self.experiment.rot_ratios,
                                              image_display_order=self.experiment.image_display_order)
        self.raw_image_grid = ImageGridWidget(self.experiment.fl_images,
                                              ratio_images=self.experiment.rot_ratios,
                                              image_display_order=self.experiment.image_display_order)

        self.ui.rotatedImagesBox.layout().addWidget(self.rot_image_grid)
        self.ui.rawImagesBox.layout().addWidget(self.raw_image_grid)

        # Set up plots
        self.intensity_plot_widget = ProfilePlotGridWidget(
            self.experiment.trimmed_intensity_data, self.experiment.get_scaled_region_boundaries()
        )
        self.ui.intensityPlotBox.layout().addWidget(self.intensity_plot_widget)

        # Event Handlers
        self.ui.horizontalSlider.valueChanged.connect(self.handle_slider_changed)

        self.set_frame(self.frame)

    def keyPressEvent(self, a0: QtGui.QKeyEvent) -> None:
        k = a0.key()
        if k == QtCore.Qt.Key_Right:
            self.set_frame(self.frame + 1)

        if k == QtCore.Qt.Key_Left:
            self.set_frame(self.frame - 1)

    def set_frame(self, new_frame):
        self.frame = max(0, min(self.experiment.raw_image_data.strain.size - 1, new_frame))
        self.ui.horizontalSlider.setValue(self.frame)
        self.ui.label.setText(str(self.frame))
        self.rot_image_grid.set_frame(self.frame)
        self.raw_image_grid.set_frame(self.frame)
        self.intensity_plot_widget.set_frame(self.frame)
        self.ui.tableWidget.selectRow(self.frame)

    def handle_slider_changed(self):
        self.set_frame(int(self.ui.horizontalSlider.value()))

    def initialize_slider(self):
        self.ui.horizontalSlider.setMinimum(0)
        self.ui.horizontalSlider.setMaximum(self.experiment.raw_image_data.shape[0] - 1)

        self.ui.horizontalSlider.valueChanged.connect(self.handle_slider_changed)

    def initialize_table(self):
        table_data = [
            {'Strain': self.experiment.strain_map[i]}
            for i in range(len(self.experiment.strain_map) - 1)
        ]
        self.ui.tableWidget.setData(table_data)
        self.ui.tableWidget.setEditable()
        self.ui.tableWidget.setSortingEnabled(False)


class MetaLoaderWindow(QtWidgets.QMainWindow):

    def __init__(self):
        super(MetaLoaderWindow, self).__init__()

        self.ui = Ui_MetaLoader()
        self.ui.setupUi(self)

        self.ui.loadRawImageButton.clicked.connect(self.show_raw_image_loader)
        self.ui.loadExperimentButton.clicked.connect(self.get_directory)

    def show_raw_image_loader(self):
        raw_image_loader = LoadRawImageWindow()
        raw_image_loader.exec_()
        raw_image_loader.show()
        self.hide()

    def get_directory(self):
        fname = QFileDialog.getExistingDirectory(self, 'Select Experiment Directory')


class LoadRawImageWindow(QtWidgets.QDialog):

    def __init__(self):
        super(LoadRawImageWindow, self).__init__()

        self.ui = Ui_LoadRawImageFileDialog()
        self.ui.setupUi(self)

        self.ui.buttonBox.accepted.connect(self.accept)
        # self.ui.buttonBox.rejected.connect(self.reject)

        self.ui.selectImageFilePushButton.clicked.connect(self.get_raw_image_file_name)
        self.ui.addRowPushButton.clicked.connect(self.handle_add_new_row)
        self.ui.deleteRowPushButton.clicked.connect(self.handle_delete_row)

        self.ui.strainTable.setSortingEnabled(False)
        self.ui.strainTable.setEditable(True)

        self.ui.strainTable.itemChanged.connect(self.handle_table_changed)

        self.indexer_df = pd.DataFrame()

    def get_raw_image_file_name(self):
        f_name = QFileDialog.getOpenFileName(self, filter='Image Files (*.tif *.tiff)')
        file_path = Path(f_name[0])
        for fp in file_path.parent.iterdir():
            if 'indexer' in fp.name:
                self.indexer_df = pd.read_csv(fp.absolute())
                self.update_strain_table()
        self.ui.imageFileLineEdit.setText(str(file_path.absolute()))

    def handle_table_changed(self, item):
        header = self.ui.strainTable.horizontalHeaderItem(item.column()).text()
        self.indexer_df.at[item.row(), header] = item.text()
        return

    def update_strain_table(self):
        self.ui.strainTable.setData(self.indexer_df.to_dict(orient='records'))

    def handle_add_new_row(self):
        self.indexer_df = self.indexer_df.append({'Strain': '', 'Start Animal': 0, 'End Animal': 0}, ignore_index=True)
        self.update_strain_table()
        return

    def handle_delete_row(self):
        self.indexer_df.drop(self.indexer_df.tail(1).index, inplace=True)
        self.update_strain_table()

    def verify_inputs(self):
        # TODO: Do this here, or when creating an experiment I can just throw an error and quit?
        # Check image stack divisible by len(strategy)

        # Check animals start @ 1

        # Check
        return

    def accept(self):
        self.verify_inputs()


if __name__ == '__main__':
    # TODO: resize recursion bug
    pg.setConfigOptions(imageAxisOrder='row-major', antialias=True)
    qapp = QtWidgets.QApplication([])
    # window = MainWindow(reload=False)
    window = MetaLoaderWindow()
    window.show()
    qapp.exec_()

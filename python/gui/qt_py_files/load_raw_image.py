# -*- coding: utf-8 -*-

# Form implementation generated from reading ui file '/Users/sean/code/wormAnalysis/python/gui/qt_ui_files/load_raw_image.ui'
#
# Created by: PyQt5 UI code generator 5.12.1
#
# WARNING! All changes made in this file will be lost!

from PyQt5 import QtCore, QtGui, QtWidgets


class Ui_LoadExperimentDialog(object):
    def setupUi(self, LoadExperimentDialog):
        LoadExperimentDialog.setObjectName("LoadExperimentDialog")
        LoadExperimentDialog.resize(583, 331)
        self.gridLayout_2 = QtWidgets.QGridLayout(LoadExperimentDialog)
        self.gridLayout_2.setObjectName("gridLayout_2")
        self.buttonBox = QtWidgets.QDialogButtonBox(LoadExperimentDialog)
        self.buttonBox.setOrientation(QtCore.Qt.Horizontal)
        self.buttonBox.setStandardButtons(QtWidgets.QDialogButtonBox.Cancel|QtWidgets.QDialogButtonBox.Ok)
        self.buttonBox.setObjectName("buttonBox")
        self.gridLayout_2.addWidget(self.buttonBox, 2, 0, 1, 1)
        self.gridLayout = QtWidgets.QGridLayout()
        self.gridLayout.setObjectName("gridLayout")
        self.imagingStrategyLineEdit = QtWidgets.QLineEdit(LoadExperimentDialog)
        self.imagingStrategyLineEdit.setObjectName("imagingStrategyLineEdit")
        self.gridLayout.addWidget(self.imagingStrategyLineEdit, 3, 1, 1, 1)
        self.label_2 = QtWidgets.QLabel(LoadExperimentDialog)
        self.label_2.setObjectName("label_2")
        self.gridLayout.addWidget(self.label_2, 3, 0, 1, 1)
        self.selectImageFilePushButton = QtWidgets.QPushButton(LoadExperimentDialog)
        self.selectImageFilePushButton.setObjectName("selectImageFilePushButton")
        self.gridLayout.addWidget(self.selectImageFilePushButton, 0, 0, 1, 1)
        self.verticalLayout = QtWidgets.QVBoxLayout()
        self.verticalLayout.setObjectName("verticalLayout")
        spacerItem = QtWidgets.QSpacerItem(20, 40, QtWidgets.QSizePolicy.Minimum, QtWidgets.QSizePolicy.Expanding)
        self.verticalLayout.addItem(spacerItem)
        self.addRowPushButton = QtWidgets.QPushButton(LoadExperimentDialog)
        self.addRowPushButton.setObjectName("addRowPushButton")
        self.verticalLayout.addWidget(self.addRowPushButton)
        self.deleteRowPushButton = QtWidgets.QPushButton(LoadExperimentDialog)
        self.deleteRowPushButton.setObjectName("deleteRowPushButton")
        self.verticalLayout.addWidget(self.deleteRowPushButton)
        spacerItem1 = QtWidgets.QSpacerItem(20, 40, QtWidgets.QSizePolicy.Minimum, QtWidgets.QSizePolicy.Expanding)
        self.verticalLayout.addItem(spacerItem1)
        self.gridLayout.addLayout(self.verticalLayout, 2, 0, 1, 1)
        self.imageFileLineEdit = QtWidgets.QLineEdit(LoadExperimentDialog)
        self.imageFileLineEdit.setObjectName("imageFileLineEdit")
        self.gridLayout.addWidget(self.imageFileLineEdit, 0, 1, 1, 1)
        self.strainTable = TableWidget(LoadExperimentDialog)
        self.strainTable.setObjectName("strainTable")
        self.strainTable.setColumnCount(3)
        self.strainTable.setRowCount(0)
        item = QtWidgets.QTableWidgetItem()
        self.strainTable.setHorizontalHeaderItem(0, item)
        item = QtWidgets.QTableWidgetItem()
        self.strainTable.setHorizontalHeaderItem(1, item)
        item = QtWidgets.QTableWidgetItem()
        self.strainTable.setHorizontalHeaderItem(2, item)
        self.strainTable.horizontalHeader().setStretchLastSection(True)
        self.gridLayout.addWidget(self.strainTable, 2, 1, 1, 1)
        self.gridLayout_2.addLayout(self.gridLayout, 1, 0, 1, 1)
        self.label = QtWidgets.QLabel(LoadExperimentDialog)
        self.label.setObjectName("label")
        self.gridLayout_2.addWidget(self.label, 0, 0, 1, 1)

        self.retranslateUi(LoadExperimentDialog)
        self.buttonBox.accepted.connect(LoadExperimentDialog.accept)
        self.buttonBox.rejected.connect(LoadExperimentDialog.reject)
        QtCore.QMetaObject.connectSlotsByName(LoadExperimentDialog)

    def retranslateUi(self, LoadExperimentDialog):
        _translate = QtCore.QCoreApplication.translate
        LoadExperimentDialog.setWindowTitle(_translate("LoadExperimentDialog", "Load Experiment"))
        self.imagingStrategyLineEdit.setText(_translate("LoadExperimentDialog", "TL/410/470/410/470"))
        self.label_2.setText(_translate("LoadExperimentDialog", "Imaging Strategy"))
        self.selectImageFilePushButton.setText(_translate("LoadExperimentDialog", "Select Image File"))
        self.addRowPushButton.setText(_translate("LoadExperimentDialog", "Add Row"))
        self.deleteRowPushButton.setText(_translate("LoadExperimentDialog", "Delete Row"))
        item = self.strainTable.horizontalHeaderItem(0)
        item.setText(_translate("LoadExperimentDialog", "Strain"))
        item = self.strainTable.horizontalHeaderItem(1)
        item.setText(_translate("LoadExperimentDialog", "Start Animal"))
        item = self.strainTable.horizontalHeaderItem(2)
        item.setText(_translate("LoadExperimentDialog", "End Animal"))
        self.label.setText(_translate("LoadExperimentDialog", "Load Raw Image Stack"))


from pyqtgraph import TableWidget

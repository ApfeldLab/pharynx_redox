# -*- coding: utf-8 -*-

# Form implementation generated from reading ui file '/Users/sean/code/pharynx_redox/python/gui/qt_ui_files/load_raw_image.ui'
#
# Created by: PyQt5 UI code generator 5.12.1
#
# WARNING! All changes made in this file will be lost!

from PyQt5 import QtCore, QtGui, QtWidgets


class Ui_LoadRawImageFileDialog(object):
    def setupUi(self, LoadRawImageFileDialog):
        LoadRawImageFileDialog.setObjectName("LoadRawImageFileDialog")
        LoadRawImageFileDialog.resize(583, 344)
        self.gridLayout_2 = QtWidgets.QGridLayout(LoadRawImageFileDialog)
        self.gridLayout_2.setObjectName("gridLayout_2")
        self.buttonBox = QtWidgets.QDialogButtonBox(LoadRawImageFileDialog)
        self.buttonBox.setOrientation(QtCore.Qt.Horizontal)
        self.buttonBox.setStandardButtons(
            QtWidgets.QDialogButtonBox.Cancel | QtWidgets.QDialogButtonBox.Ok
        )
        self.buttonBox.setObjectName("buttonBox")
        self.gridLayout_2.addWidget(self.buttonBox, 2, 0, 1, 1)
        self.gridLayout = QtWidgets.QGridLayout()
        self.gridLayout.setObjectName("gridLayout")
        self.imagingStrategyLineEdit = QtWidgets.QLineEdit(LoadRawImageFileDialog)
        self.imagingStrategyLineEdit.setObjectName("imagingStrategyLineEdit")
        self.gridLayout.addWidget(self.imagingStrategyLineEdit, 4, 1, 1, 1)
        self.label_2 = QtWidgets.QLabel(LoadRawImageFileDialog)
        self.label_2.setObjectName("label_2")
        self.gridLayout.addWidget(self.label_2, 4, 0, 1, 1)
        self.selectImageFilePushButton = QtWidgets.QPushButton(LoadRawImageFileDialog)
        self.selectImageFilePushButton.setObjectName("selectImageFilePushButton")
        self.gridLayout.addWidget(self.selectImageFilePushButton, 1, 0, 1, 1)
        self.verticalLayout = QtWidgets.QVBoxLayout()
        self.verticalLayout.setObjectName("verticalLayout")
        spacerItem = QtWidgets.QSpacerItem(
            20, 40, QtWidgets.QSizePolicy.Minimum, QtWidgets.QSizePolicy.Expanding
        )
        self.verticalLayout.addItem(spacerItem)
        self.addRowPushButton = QtWidgets.QPushButton(LoadRawImageFileDialog)
        self.addRowPushButton.setObjectName("addRowPushButton")
        self.verticalLayout.addWidget(self.addRowPushButton)
        self.deleteRowPushButton = QtWidgets.QPushButton(LoadRawImageFileDialog)
        self.deleteRowPushButton.setObjectName("deleteRowPushButton")
        self.verticalLayout.addWidget(self.deleteRowPushButton)
        spacerItem1 = QtWidgets.QSpacerItem(
            20, 40, QtWidgets.QSizePolicy.Minimum, QtWidgets.QSizePolicy.Expanding
        )
        self.verticalLayout.addItem(spacerItem1)
        self.gridLayout.addLayout(self.verticalLayout, 3, 0, 1, 1)
        self.imageFileLineEdit = QtWidgets.QLineEdit(LoadRawImageFileDialog)
        self.imageFileLineEdit.setObjectName("imageFileLineEdit")
        self.gridLayout.addWidget(self.imageFileLineEdit, 1, 1, 1, 1)
        self.horizontalLayout_2 = QtWidgets.QHBoxLayout()
        self.horizontalLayout_2.setObjectName("horizontalLayout_2")
        self.strainTable = TableWidget(LoadRawImageFileDialog)
        sizePolicy = QtWidgets.QSizePolicy(
            QtWidgets.QSizePolicy.Fixed, QtWidgets.QSizePolicy.Expanding
        )
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.strainTable.sizePolicy().hasHeightForWidth())
        self.strainTable.setSizePolicy(sizePolicy)
        self.strainTable.setObjectName("strainTable")
        self.strainTable.setColumnCount(3)
        self.strainTable.setRowCount(0)
        item = QtWidgets.QTableWidgetItem()
        self.strainTable.setHorizontalHeaderItem(0, item)
        item = QtWidgets.QTableWidgetItem()
        self.strainTable.setHorizontalHeaderItem(1, item)
        item = QtWidgets.QTableWidgetItem()
        self.strainTable.setHorizontalHeaderItem(2, item)
        self.strainTable.horizontalHeader().setCascadingSectionResizes(False)
        self.strainTable.horizontalHeader().setDefaultSectionSize(80)
        self.strainTable.horizontalHeader().setSortIndicatorShown(False)
        self.strainTable.horizontalHeader().setStretchLastSection(True)
        self.strainTable.verticalHeader().setCascadingSectionResizes(False)
        self.strainTable.verticalHeader().setSortIndicatorShown(False)
        self.horizontalLayout_2.addWidget(self.strainTable)
        spacerItem2 = QtWidgets.QSpacerItem(
            40, 20, QtWidgets.QSizePolicy.Expanding, QtWidgets.QSizePolicy.Minimum
        )
        self.horizontalLayout_2.addItem(spacerItem2)
        self.gridLayout.addLayout(self.horizontalLayout_2, 3, 1, 1, 1)
        self.gridLayout_2.addLayout(self.gridLayout, 1, 0, 1, 1)
        self.label = QtWidgets.QLabel(LoadRawImageFileDialog)
        self.label.setObjectName("label")
        self.gridLayout_2.addWidget(self.label, 0, 0, 1, 1)

        self.retranslateUi(LoadRawImageFileDialog)
        self.buttonBox.accepted.connect(LoadRawImageFileDialog.accept)
        self.buttonBox.rejected.connect(LoadRawImageFileDialog.reject)
        QtCore.QMetaObject.connectSlotsByName(LoadRawImageFileDialog)

    def retranslateUi(self, LoadRawImageFileDialog):
        _translate = QtCore.QCoreApplication.translate
        LoadRawImageFileDialog.setWindowTitle(
            _translate("LoadRawImageFileDialog", "Load Raw Image File")
        )
        self.imagingStrategyLineEdit.setText(
            _translate("LoadRawImageFileDialog", "TL/410/470/410/470")
        )
        self.label_2.setText(_translate("LoadRawImageFileDialog", "Imaging Strategy"))
        self.selectImageFilePushButton.setText(
            _translate("LoadRawImageFileDialog", "Select Image File")
        )
        self.addRowPushButton.setText(_translate("LoadRawImageFileDialog", "Add Row"))
        self.deleteRowPushButton.setText(
            _translate("LoadRawImageFileDialog", "Delete Row")
        )
        self.strainTable.setSortingEnabled(False)
        item = self.strainTable.horizontalHeaderItem(0)
        item.setText(_translate("LoadRawImageFileDialog", "Strain"))
        item = self.strainTable.horizontalHeaderItem(1)
        item.setText(_translate("LoadRawImageFileDialog", "Start Animal"))
        item = self.strainTable.horizontalHeaderItem(2)
        item.setText(_translate("LoadRawImageFileDialog", "End Animal"))
        self.label.setText(_translate("LoadRawImageFileDialog", "Load Raw Image Stack"))


from pyqtgraph import TableWidget

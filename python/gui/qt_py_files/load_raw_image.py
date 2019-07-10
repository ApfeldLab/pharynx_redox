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
        LoadExperimentDialog.resize(494, 311)
        self.gridLayout_2 = QtWidgets.QGridLayout(LoadExperimentDialog)
        self.gridLayout_2.setObjectName("gridLayout_2")
        self.buttonBox = QtWidgets.QDialogButtonBox(LoadExperimentDialog)
        self.buttonBox.setOrientation(QtCore.Qt.Horizontal)
        self.buttonBox.setStandardButtons(QtWidgets.QDialogButtonBox.Cancel|QtWidgets.QDialogButtonBox.Ok)
        self.buttonBox.setObjectName("buttonBox")
        self.gridLayout_2.addWidget(self.buttonBox, 2, 0, 1, 1)
        self.gridLayout = QtWidgets.QGridLayout()
        self.gridLayout.setObjectName("gridLayout")
        self.tableWidget = QtWidgets.QTableWidget(LoadExperimentDialog)
        self.tableWidget.setObjectName("tableWidget")
        self.tableWidget.setColumnCount(3)
        self.tableWidget.setRowCount(0)
        item = QtWidgets.QTableWidgetItem()
        self.tableWidget.setHorizontalHeaderItem(0, item)
        item = QtWidgets.QTableWidgetItem()
        self.tableWidget.setHorizontalHeaderItem(1, item)
        item = QtWidgets.QTableWidgetItem()
        self.tableWidget.setHorizontalHeaderItem(2, item)
        self.tableWidget.horizontalHeader().setStretchLastSection(True)
        self.gridLayout.addWidget(self.tableWidget, 2, 1, 1, 1)
        self.lineEdit = QtWidgets.QLineEdit(LoadExperimentDialog)
        self.lineEdit.setObjectName("lineEdit")
        self.gridLayout.addWidget(self.lineEdit, 0, 1, 1, 1)
        self.pushButton = QtWidgets.QPushButton(LoadExperimentDialog)
        self.pushButton.setObjectName("pushButton")
        self.gridLayout.addWidget(self.pushButton, 0, 0, 1, 1)
        self.label_2 = QtWidgets.QLabel(LoadExperimentDialog)
        self.label_2.setObjectName("label_2")
        self.gridLayout.addWidget(self.label_2, 3, 0, 1, 1)
        self.comboBox = QtWidgets.QComboBox(LoadExperimentDialog)
        self.comboBox.setObjectName("comboBox")
        self.comboBox.addItem("")
        self.comboBox.addItem("")
        self.comboBox.addItem("")
        self.comboBox.addItem("")
        self.comboBox.addItem("")
        self.comboBox.addItem("")
        self.gridLayout.addWidget(self.comboBox, 3, 1, 1, 1)
        self.verticalLayout = QtWidgets.QVBoxLayout()
        self.verticalLayout.setObjectName("verticalLayout")
        spacerItem = QtWidgets.QSpacerItem(20, 40, QtWidgets.QSizePolicy.Minimum, QtWidgets.QSizePolicy.Expanding)
        self.verticalLayout.addItem(spacerItem)
        self.pushButton_2 = QtWidgets.QPushButton(LoadExperimentDialog)
        self.pushButton_2.setObjectName("pushButton_2")
        self.verticalLayout.addWidget(self.pushButton_2)
        self.pushButton_3 = QtWidgets.QPushButton(LoadExperimentDialog)
        self.pushButton_3.setObjectName("pushButton_3")
        self.verticalLayout.addWidget(self.pushButton_3)
        spacerItem1 = QtWidgets.QSpacerItem(20, 40, QtWidgets.QSizePolicy.Minimum, QtWidgets.QSizePolicy.Expanding)
        self.verticalLayout.addItem(spacerItem1)
        self.gridLayout.addLayout(self.verticalLayout, 2, 0, 1, 1)
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
        item = self.tableWidget.horizontalHeaderItem(0)
        item.setText(_translate("LoadExperimentDialog", "Strain"))
        item = self.tableWidget.horizontalHeaderItem(1)
        item.setText(_translate("LoadExperimentDialog", "Start Animal"))
        item = self.tableWidget.horizontalHeaderItem(2)
        item.setText(_translate("LoadExperimentDialog", "End Animal"))
        self.pushButton.setText(_translate("LoadExperimentDialog", "Select Image File"))
        self.label_2.setText(_translate("LoadExperimentDialog", "Imaging Strategy"))
        self.comboBox.setItemText(0, _translate("LoadExperimentDialog", "TL/410/470/410/470"))
        self.comboBox.setItemText(1, _translate("LoadExperimentDialog", "TL/470/410/470/410"))
        self.comboBox.setItemText(2, _translate("LoadExperimentDialog", "410/470/410/470"))
        self.comboBox.setItemText(3, _translate("LoadExperimentDialog", "470/410/470/410"))
        self.comboBox.setItemText(4, _translate("LoadExperimentDialog", "410/470"))
        self.comboBox.setItemText(5, _translate("LoadExperimentDialog", "470/410"))
        self.pushButton_2.setText(_translate("LoadExperimentDialog", "Add Row"))
        self.pushButton_3.setText(_translate("LoadExperimentDialog", "Delete Row"))
        self.label.setText(_translate("LoadExperimentDialog", "Load Raw Image Stack"))



# -*- coding: utf-8 -*-

# Form implementation generated from reading ui file '/Users/sean/code/wormAnalysis/python/gui/qt_ui_files/meta_loader.ui'
#
# Created by: PyQt5 UI code generator 5.12.1
#
# WARNING! All changes made in this file will be lost!

from PyQt5 import QtCore, QtGui, QtWidgets


class Ui_MetaLoader(object):
    def setupUi(self, MetaLoader):
        MetaLoader.setObjectName("MetaLoader")
        MetaLoader.resize(238, 84)
        self.loadRawImageButton = QtWidgets.QPushButton(MetaLoader)
        self.loadRawImageButton.setGeometry(QtCore.QRect(6, 42, 167, 32))
        self.loadRawImageButton.setObjectName("loadRawImageButton")
        self.loadExperimentButton = QtWidgets.QPushButton(MetaLoader)
        self.loadExperimentButton.setGeometry(QtCore.QRect(6, 8, 204, 32))
        self.loadExperimentButton.setObjectName("loadExperimentButton")

        self.retranslateUi(MetaLoader)
        QtCore.QMetaObject.connectSlotsByName(MetaLoader)

    def retranslateUi(self, MetaLoader):
        _translate = QtCore.QCoreApplication.translate
        MetaLoader.setWindowTitle(_translate("MetaLoader", "Load"))
        self.loadRawImageButton.setText(_translate("MetaLoader", "Load Raw Image File"))
        self.loadExperimentButton.setText(
            _translate("MetaLoader", "Load Experiment Directory")
        )

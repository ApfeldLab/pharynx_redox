# -*- coding: utf-8 -*-

# Form implementation generated from reading ui file '/Users/sean/code/pharedox/pharedox/gui/qt_ui_files/pipeline_buttons.ui'
#
# Created by: PyQt5 UI code generator 5.12.3
#
# WARNING! All changes made in this file will be lost!


from PyQt5 import QtCore, QtGui, QtWidgets


class Ui_Form(object):
    def setupUi(self, Form):
        Form.setObjectName("Form")
        Form.resize(287, 425)
        Form.setLayoutDirection(QtCore.Qt.LeftToRight)
        self.verticalLayout = QtWidgets.QVBoxLayout(Form)
        self.verticalLayout.setObjectName("verticalLayout")
        self.groupBox = QtWidgets.QGroupBox(Form)
        self.groupBox.setObjectName("groupBox")
        self.verticalLayout_2 = QtWidgets.QVBoxLayout(self.groupBox)
        self.verticalLayout_2.setObjectName("verticalLayout_2")
        self.segmentButton = QtWidgets.QPushButton(self.groupBox)
        self.segmentButton.setObjectName("segmentButton")
        self.verticalLayout_2.addWidget(self.segmentButton)
        self.thresholdSlider = QtWidgets.QSlider(self.groupBox)
        self.thresholdSlider.setOrientation(QtCore.Qt.Horizontal)
        self.thresholdSlider.setObjectName("thresholdSlider")
        self.verticalLayout_2.addWidget(self.thresholdSlider)
        self.thresholdSpinBox = QtWidgets.QSpinBox(self.groupBox)
        sizePolicy = QtWidgets.QSizePolicy(
            QtWidgets.QSizePolicy.Expanding, QtWidgets.QSizePolicy.Fixed
        )
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(
            self.thresholdSpinBox.sizePolicy().hasHeightForWidth()
        )
        self.thresholdSpinBox.setSizePolicy(sizePolicy)
        self.thresholdSpinBox.setMinimumSize(QtCore.QSize(150, 0))
        self.thresholdSpinBox.setMaximumSize(QtCore.QSize(16777215, 16777215))
        self.thresholdSpinBox.setAlignment(
            QtCore.Qt.AlignRight | QtCore.Qt.AlignTrailing | QtCore.Qt.AlignVCenter
        )
        self.thresholdSpinBox.setMaximum(65535)
        self.thresholdSpinBox.setSingleStep(100)
        self.thresholdSpinBox.setObjectName("thresholdSpinBox")
        self.verticalLayout_2.addWidget(self.thresholdSpinBox)
        self.horizontalLayout = QtWidgets.QHBoxLayout()
        self.horizontalLayout.setObjectName("horizontalLayout")
        self.removeObjectsButton = QtWidgets.QPushButton(self.groupBox)
        self.removeObjectsButton.setObjectName("removeObjectsButton")
        self.horizontalLayout.addWidget(self.removeObjectsButton)
        self.smallObjectSizeSpinBox = QtWidgets.QSpinBox(self.groupBox)
        self.smallObjectSizeSpinBox.setMinimum(1)
        self.smallObjectSizeSpinBox.setMaximum(10000)
        self.smallObjectSizeSpinBox.setProperty("value", 5)
        self.smallObjectSizeSpinBox.setObjectName("smallObjectSizeSpinBox")
        self.horizontalLayout.addWidget(self.smallObjectSizeSpinBox)
        self.verticalLayout_2.addLayout(self.horizontalLayout)
        self.verticalLayout.addWidget(self.groupBox)
        spacerItem = QtWidgets.QSpacerItem(
            20, 40, QtWidgets.QSizePolicy.Minimum, QtWidgets.QSizePolicy.Expanding
        )
        self.verticalLayout.addItem(spacerItem)
        self.runPharynxButton = QtWidgets.QPushButton(Form)
        self.runPharynxButton.setObjectName("runPharynxButton")
        self.verticalLayout.addWidget(self.runPharynxButton)
        self.runNeuronsButton = QtWidgets.QPushButton(Form)
        self.runNeuronsButton.setObjectName("runNeuronsButton")
        self.verticalLayout.addWidget(self.runNeuronsButton)

        self.retranslateUi(Form)
        QtCore.QMetaObject.connectSlotsByName(Form)

    def retranslateUi(self, Form):
        _translate = QtCore.QCoreApplication.translate
        Form.setWindowTitle(_translate("Form", "Form"))
        self.groupBox.setTitle(_translate("Form", "Segmentation"))
        self.segmentButton.setText(_translate("Form", "Segment"))
        self.removeObjectsButton.setText(_translate("Form", "Remove Objects <"))
        self.runPharynxButton.setText(_translate("Form", "Analyze Pharynxes"))
        self.runNeuronsButton.setText(_translate("Form", "Analyze Neurons"))

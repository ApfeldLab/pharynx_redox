# -*- coding: utf-8 -*-

# Form implementation generated from reading ui file '/Users/sean/code/pharynx_redox/pharynx_redox/gui/qt_ui_files/image_stack_widget.ui'
#
# Created by: PyQt5 UI code generator 5.9.2
#
# WARNING! All changes made in this file will be lost!

from PyQt5 import QtCore, QtGui, QtWidgets


class Ui_XArrayDisplayWidget(object):
    def setupUi(self, XArrayDisplayWidget):
        XArrayDisplayWidget.setObjectName("XArrayDisplayWidget")
        XArrayDisplayWidget.resize(1680, 1005)
        self.gridLayout_2 = QtWidgets.QGridLayout(XArrayDisplayWidget)
        self.gridLayout_2.setObjectName("gridLayout_2")
        self.widget = QtWidgets.QWidget(XArrayDisplayWidget)
        self.widget.setMinimumSize(QtCore.QSize(300, 0))
        self.widget.setMaximumSize(QtCore.QSize(200, 16777215))
        self.widget.setObjectName("widget")
        self.verticalLayout_2 = QtWidgets.QVBoxLayout(self.widget)
        self.verticalLayout_2.setObjectName("verticalLayout_2")
        self.groupBox_2 = QtWidgets.QGroupBox(self.widget)
        self.groupBox_2.setObjectName("groupBox_2")
        self.verticalLayout_3 = QtWidgets.QVBoxLayout(self.groupBox_2)
        self.verticalLayout_3.setObjectName("verticalLayout_3")
        self.label = QtWidgets.QLabel(self.groupBox_2)
        self.label.setMaximumSize(QtCore.QSize(16777215, 16777215))
        self.label.setObjectName("label")
        self.verticalLayout_3.addWidget(self.label)
        self.wvlBox = QtWidgets.QComboBox(self.groupBox_2)
        self.wvlBox.setMaximumSize(QtCore.QSize(16777215, 16777215))
        self.wvlBox.setObjectName("wvlBox")
        self.verticalLayout_3.addWidget(self.wvlBox)
        self.label_2 = QtWidgets.QLabel(self.groupBox_2)
        self.label_2.setMaximumSize(QtCore.QSize(16777215, 16777215))
        self.label_2.setObjectName("label_2")
        self.verticalLayout_3.addWidget(self.label_2)
        self.pairSlider = QtWidgets.QSlider(self.groupBox_2)
        self.pairSlider.setMaximumSize(QtCore.QSize(16777215, 16777215))
        self.pairSlider.setOrientation(QtCore.Qt.Horizontal)
        self.pairSlider.setObjectName("pairSlider")
        self.verticalLayout_3.addWidget(self.pairSlider)
        self.displayMaskCheckbox = QtWidgets.QCheckBox(self.groupBox_2)
        self.displayMaskCheckbox.setObjectName("displayMaskCheckbox")
        self.verticalLayout_3.addWidget(self.displayMaskCheckbox)
        self.horizontalLayout = QtWidgets.QHBoxLayout()
        self.horizontalLayout.setObjectName("horizontalLayout")
        self.editMaskCheckBox = QtWidgets.QCheckBox(self.groupBox_2)
        self.editMaskCheckBox.setObjectName("editMaskCheckBox")
        self.horizontalLayout.addWidget(self.editMaskCheckBox)
        self.drawToolButton = QtWidgets.QToolButton(self.groupBox_2)
        self.drawToolButton.setObjectName("drawToolButton")
        self.horizontalLayout.addWidget(self.drawToolButton)
        self.kernelRadiusSpinBox = QtWidgets.QSpinBox(self.groupBox_2)
        self.kernelRadiusSpinBox.setProperty("value", 5)
        self.kernelRadiusSpinBox.setObjectName("kernelRadiusSpinBox")
        self.horizontalLayout.addWidget(self.kernelRadiusSpinBox)
        self.verticalLayout_3.addLayout(self.horizontalLayout)
        self.verticalLayout_2.addWidget(self.groupBox_2)
        self.line_2 = QtWidgets.QFrame(self.widget)
        self.line_2.setFrameShape(QtWidgets.QFrame.HLine)
        self.line_2.setFrameShadow(QtWidgets.QFrame.Sunken)
        self.line_2.setObjectName("line_2")
        self.verticalLayout_2.addWidget(self.line_2)
        spacerItem = QtWidgets.QSpacerItem(
            20, 40, QtWidgets.QSizePolicy.Minimum, QtWidgets.QSizePolicy.Expanding
        )
        self.verticalLayout_2.addItem(spacerItem)
        self.line = QtWidgets.QFrame(self.widget)
        self.line.setFrameShape(QtWidgets.QFrame.HLine)
        self.line.setFrameShadow(QtWidgets.QFrame.Sunken)
        self.line.setObjectName("line")
        self.verticalLayout_2.addWidget(self.line)
        self.propertiesGroupBox = QtWidgets.QGroupBox(self.widget)
        self.propertiesGroupBox.setMaximumSize(QtCore.QSize(300, 16777215))
        self.propertiesGroupBox.setObjectName("propertiesGroupBox")
        self.gridLayout_3 = QtWidgets.QGridLayout(self.propertiesGroupBox)
        self.gridLayout_3.setObjectName("gridLayout_3")
        self.propertiesDataTreeWidget = DataTreeWidget(self.propertiesGroupBox)
        self.propertiesDataTreeWidget.setObjectName("propertiesDataTreeWidget")
        self.gridLayout_3.addWidget(self.propertiesDataTreeWidget, 0, 0, 1, 1)
        self.verticalLayout_2.addWidget(self.propertiesGroupBox)
        self.gridLayout_2.addWidget(self.widget, 0, 0, 1, 1)
        self.verticalLayout = QtWidgets.QVBoxLayout()
        self.verticalLayout.setSpacing(1)
        self.verticalLayout.setObjectName("verticalLayout")
        self.tabWidget = QtWidgets.QTabWidget(XArrayDisplayWidget)
        self.tabWidget.setObjectName("tabWidget")
        self.tab = QtWidgets.QWidget()
        self.tab.setObjectName("tab")
        self.gridLayout = QtWidgets.QGridLayout(self.tab)
        self.gridLayout.setObjectName("gridLayout")
        self.ImageViewBox = ImageView(self.tab)
        self.ImageViewBox.setMinimumSize(QtCore.QSize(400, 0))
        self.ImageViewBox.setObjectName("ImageViewBox")
        self.gridLayout.addWidget(self.ImageViewBox, 0, 0, 1, 1)
        self.tabWidget.addTab(self.tab, "")
        self.tab_2 = QtWidgets.QWidget()
        self.tab_2.setObjectName("tab_2")
        self.tabWidget.addTab(self.tab_2, "")
        self.verticalLayout.addWidget(self.tabWidget)
        self.gridLayout_2.addLayout(self.verticalLayout, 0, 1, 1, 1)

        self.retranslateUi(XArrayDisplayWidget)
        self.tabWidget.setCurrentIndex(0)
        QtCore.QMetaObject.connectSlotsByName(XArrayDisplayWidget)

    def retranslateUi(self, XArrayDisplayWidget):
        _translate = QtCore.QCoreApplication.translate
        XArrayDisplayWidget.setWindowTitle(_translate("XArrayDisplayWidget", "Form"))
        self.groupBox_2.setTitle(_translate("XArrayDisplayWidget", "Controls"))
        self.label.setText(_translate("XArrayDisplayWidget", "Wavelength"))
        self.label_2.setText(_translate("XArrayDisplayWidget", "Pair"))
        self.displayMaskCheckbox.setText(
            _translate("XArrayDisplayWidget", "Display Mask")
        )
        self.editMaskCheckBox.setText(
            _translate("XArrayDisplayWidget", "Edit Mask (m)")
        )
        self.drawToolButton.setText(_translate("XArrayDisplayWidget", "Draw (d)"))
        self.propertiesGroupBox.setTitle(
            _translate("XArrayDisplayWidget", "Properties")
        )
        self.tabWidget.setTabText(
            self.tabWidget.indexOf(self.tab), _translate("XArrayDisplayWidget", "Raw")
        )
        self.tabWidget.setTabText(
            self.tabWidget.indexOf(self.tab_2),
            _translate("XArrayDisplayWidget", "Segmented"),
        )


from pyqtgraph import DataTreeWidget, ImageView

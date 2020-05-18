# -*- coding: utf-8 -*-

# Form implementation generated from reading ui file '/Users/sean/code/pharedox/pharedox/gui/qt_ui_files/mask_painter_widget.ui'
#
# Created by: PyQt5 UI code generator 5.12.3
#
# WARNING! All changes made in this file will be lost!


from PyQt5 import QtCore, QtGui, QtWidgets


class Ui_Form(object):
    def setupUi(self, Form):
        Form.setObjectName("Form")
        Form.resize(559, 347)
        self.horizontalLayout = QtWidgets.QHBoxLayout(Form)
        self.horizontalLayout.setObjectName("horizontalLayout")
        self.frame = QtWidgets.QFrame(Form)
        self.frame.setMaximumSize(QtCore.QSize(200, 16777215))
        self.frame.setFrameShape(QtWidgets.QFrame.StyledPanel)
        self.frame.setFrameShadow(QtWidgets.QFrame.Raised)
        self.frame.setObjectName("frame")
        self.verticalLayout = QtWidgets.QVBoxLayout(self.frame)
        self.verticalLayout.setObjectName("verticalLayout")
        self.checkBox = QtWidgets.QCheckBox(self.frame)
        self.checkBox.setObjectName("checkBox")
        self.verticalLayout.addWidget(self.checkBox)
        self.horizontalSlider = QtWidgets.QSlider(self.frame)
        self.horizontalSlider.setOrientation(QtCore.Qt.Horizontal)
        self.horizontalSlider.setObjectName("horizontalSlider")
        self.verticalLayout.addWidget(self.horizontalSlider)
        self.toolButton = QtWidgets.QToolButton(self.frame)
        self.toolButton.setObjectName("toolButton")
        self.verticalLayout.addWidget(self.toolButton)
        self.horizontalLayout.addWidget(self.frame)
        self.widget = QtWidgets.QWidget(Form)
        self.widget.setObjectName("widget")
        self.horizontalLayout.addWidget(self.widget)

        self.retranslateUi(Form)
        QtCore.QMetaObject.connectSlotsByName(Form)

    def retranslateUi(self, Form):
        _translate = QtCore.QCoreApplication.translate
        Form.setWindowTitle(_translate("Form", "Form"))
        self.checkBox.setText(_translate("Form", "Show / Edit Mask"))
        self.toolButton.setText(_translate("Form", "paint"))

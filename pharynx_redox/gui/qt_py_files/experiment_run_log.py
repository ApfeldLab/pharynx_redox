# -*- coding: utf-8 -*-

# Form implementation generated from reading ui file '/Users/sean/code/pharynx_redox/pharynx_redox/gui/qt_ui_files/experiment_run_log.ui'
#
# Created by: PyQt5 UI code generator 5.12.3
#
# WARNING! All changes made in this file will be lost!


from PyQt5 import QtCore, QtGui, QtWidgets


class Ui_Form(object):
    def setupUi(self, Form):
        Form.setObjectName("Form")
        Form.resize(570, 423)
        self.verticalLayout = QtWidgets.QVBoxLayout(Form)
        self.verticalLayout.setObjectName("verticalLayout")
        self.label = QtWidgets.QLabel(Form)
        self.label.setObjectName("label")
        self.verticalLayout.addWidget(self.label)
        self.logTextBox = QtWidgets.QPlainTextEdit(Form)
        self.logTextBox.setReadOnly(True)
        self.logTextBox.setObjectName("logTextBox")
        self.verticalLayout.addWidget(self.logTextBox)

        self.retranslateUi(Form)
        QtCore.QMetaObject.connectSlotsByName(Form)

    def retranslateUi(self, Form):
        _translate = QtCore.QCoreApplication.translate
        Form.setWindowTitle(_translate("Form", "Form"))
        self.label.setText(_translate("Form", "Running Analysis"))

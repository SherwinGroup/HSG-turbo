# -*- coding: utf-8 -*-

# Form implementation generated from reading ui file 'C:\Users\FELLab\Documents\GitHub\Interactivepg-waffle\interactivePG\fixes\LegendSettings.ui'
#
# Created by: PyQt5 UI code generator 5.6
#
# WARNING! All changes made in this file will be lost!

from PyQt5 import QtCore, QtGui, QtWidgets

class Ui_LegendSettingsDialog(object):
    def setupUi(self, LegendSettingsDialog):
        LegendSettingsDialog.setObjectName("LegendSettingsDialog")
        LegendSettingsDialog.resize(207, 319)
        self.verticalLayout = QtWidgets.QVBoxLayout(LegendSettingsDialog)
        self.verticalLayout.setObjectName("verticalLayout")
        self.formLayout = QtWidgets.QFormLayout()
        self.formLayout.setObjectName("formLayout")
        self.bBGColor = ColorButton(LegendSettingsDialog)
        self.bBGColor.setText("")
        self.bBGColor.setObjectName("bBGColor")
        self.formLayout.setWidget(0, QtWidgets.QFormLayout.FieldRole, self.bBGColor)
        self.label = QtWidgets.QLabel(LegendSettingsDialog)
        self.label.setObjectName("label")
        self.formLayout.setWidget(0, QtWidgets.QFormLayout.LabelRole, self.label)
        self.label_2 = QtWidgets.QLabel(LegendSettingsDialog)
        self.label_2.setObjectName("label_2")
        self.formLayout.setWidget(1, QtWidgets.QFormLayout.LabelRole, self.label_2)
        self.bBorderColor = ColorButton(LegendSettingsDialog)
        self.bBorderColor.setText("")
        self.bBorderColor.setObjectName("bBorderColor")
        self.formLayout.setWidget(1, QtWidgets.QFormLayout.FieldRole, self.bBorderColor)
        self.label_3 = QtWidgets.QLabel(LegendSettingsDialog)
        self.label_3.setObjectName("label_3")
        self.formLayout.setWidget(2, QtWidgets.QFormLayout.LabelRole, self.label_3)
        self.sbFontSize = SpinBox(LegendSettingsDialog)
        self.sbFontSize.setObjectName("sbFontSize")
        self.formLayout.setWidget(2, QtWidgets.QFormLayout.FieldRole, self.sbFontSize)
        self.verticalLayout.addLayout(self.formLayout)
        self.groupBox = QtWidgets.QGroupBox(LegendSettingsDialog)
        self.groupBox.setFlat(True)
        self.groupBox.setObjectName("groupBox")
        self.verticalLayout_2 = QtWidgets.QVBoxLayout(self.groupBox)
        self.verticalLayout_2.setContentsMargins(0, 10, 0, 0)
        self.verticalLayout_2.setObjectName("verticalLayout_2")
        self.teDesc = QtWidgets.QTextEdit(self.groupBox)
        self.teDesc.setObjectName("teDesc")
        self.verticalLayout_2.addWidget(self.teDesc)
        self.widget = QtWidgets.QWidget(self.groupBox)
        self.widget.setObjectName("widget")
        self.verticalLayout_2.addWidget(self.widget)
        self.verticalLayout.addWidget(self.groupBox)
        self.buttonBox = QtWidgets.QDialogButtonBox(LegendSettingsDialog)
        self.buttonBox.setOrientation(QtCore.Qt.Horizontal)
        self.buttonBox.setStandardButtons(QtWidgets.QDialogButtonBox.Cancel|QtWidgets.QDialogButtonBox.Ok)
        self.buttonBox.setObjectName("buttonBox")
        self.verticalLayout.addWidget(self.buttonBox)

        self.retranslateUi(LegendSettingsDialog)
        self.buttonBox.accepted.connect(LegendSettingsDialog.accept)
        self.buttonBox.rejected.connect(LegendSettingsDialog.reject)
        QtCore.QMetaObject.connectSlotsByName(LegendSettingsDialog)

    def retranslateUi(self, LegendSettingsDialog):
        _translate = QtCore.QCoreApplication.translate
        LegendSettingsDialog.setWindowTitle(_translate("LegendSettingsDialog", "Legend Settings"))
        self.label.setText(_translate("LegendSettingsDialog", "Background Color:"))
        self.label_2.setText(_translate("LegendSettingsDialog", "Border Color:"))
        self.label_3.setText(_translate("LegendSettingsDialog", "Font Size"))
        self.groupBox.setTitle(_translate("LegendSettingsDialog", "Item Description"))

from pyqtgraph import ColorButton, SpinBox

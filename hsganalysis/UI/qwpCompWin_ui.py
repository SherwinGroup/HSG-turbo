# -*- coding: utf-8 -*-

# Form implementation generated from reading ui file 'C:\Users\dvalovcin\Documents\GitHub\HSG-turbo\hsganalysis\UI\qwpCompWin.ui'
#
# Created by: PyQt5 UI code generator 5.6
#
# WARNING! All changes made in this file will be lost!

from PyQt5 import QtCore, QtGui, QtWidgets

class Ui_PolarimeterWindow(object):
    def setupUi(self, PolarimeterWindow):
        PolarimeterWindow.setObjectName("PolarimeterWindow")
        PolarimeterWindow.resize(800, 600)
        self.centralwidget = QtWidgets.QWidget(PolarimeterWindow)
        self.centralwidget.setObjectName("centralwidget")
        self.horizontalLayout = QtWidgets.QHBoxLayout(self.centralwidget)
        self.horizontalLayout.setObjectName("horizontalLayout")
        self.tabWidget = QtWidgets.QTabWidget(self.centralwidget)
        self.tabWidget.setObjectName("tabWidget")
        self.tabAngles = QtWidgets.QWidget()
        self.tabAngles.setObjectName("tabAngles")
        self.horizontalLayout_2 = QtWidgets.QHBoxLayout(self.tabAngles)
        self.horizontalLayout_2.setContentsMargins(0, 0, 0, 0)
        self.horizontalLayout_2.setObjectName("horizontalLayout_2")
        self.splitter = QtWidgets.QSplitter(self.tabAngles)
        self.splitter.setOrientation(QtCore.Qt.Horizontal)
        self.splitter.setObjectName("splitter")
        self.ptFile = ParameterTree(self.splitter)
        self.ptFile.setObjectName("ptFile")
        self.gAngles = DraggablePlotWidget(self.splitter)
        self.gAngles.setObjectName("gAngles")
        self.horizontalLayout_2.addWidget(self.splitter)
        self.tabWidget.addTab(self.tabAngles, "")
        self.tabRaw = QtWidgets.QWidget()
        self.tabRaw.setObjectName("tabRaw")
        self.horizontalLayout_3 = QtWidgets.QHBoxLayout(self.tabRaw)
        self.horizontalLayout_3.setContentsMargins(0, 0, 0, 0)
        self.horizontalLayout_3.setObjectName("horizontalLayout_3")
        self.splitFits = QtWidgets.QSplitter(self.tabRaw)
        self.splitFits.setOrientation(QtCore.Qt.Vertical)
        self.splitFits.setObjectName("splitFits")
        self.gRaw = DraggablePlotWidget(self.splitFits)
        self.gRaw.setObjectName("gRaw")
        self.horizontalLayout_3.addWidget(self.splitFits)
        self.tabWidget.addTab(self.tabRaw, "")
        self.horizontalLayout.addWidget(self.tabWidget)
        PolarimeterWindow.setCentralWidget(self.centralwidget)
        self.menubar = QtWidgets.QMenuBar(PolarimeterWindow)
        self.menubar.setGeometry(QtCore.QRect(0, 0, 800, 21))
        self.menubar.setObjectName("menubar")
        PolarimeterWindow.setMenuBar(self.menubar)
        self.statusbar = QtWidgets.QStatusBar(PolarimeterWindow)
        self.statusbar.setObjectName("statusbar")
        PolarimeterWindow.setStatusBar(self.statusbar)

        self.retranslateUi(PolarimeterWindow)
        self.tabWidget.setCurrentIndex(0)
        QtCore.QMetaObject.connectSlotsByName(PolarimeterWindow)

    def retranslateUi(self, PolarimeterWindow):
        _translate = QtCore.QCoreApplication.translate
        PolarimeterWindow.setWindowTitle(_translate("PolarimeterWindow", "HSG Analysis Tool"))
        self.tabWidget.setTabText(self.tabWidget.indexOf(self.tabAngles), _translate("PolarimeterWindow", "Angles"))
        self.tabWidget.setTabText(self.tabWidget.indexOf(self.tabRaw), _translate("PolarimeterWindow", "Raw Curves"))

from .draggablePlotWidget import DraggablePlotWidget
from pyqtgraph.parametertree import ParameterTree

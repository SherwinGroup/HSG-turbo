# -*- coding: utf-8 -*-

# Form implementation generated from reading ui file 'C:\Users\dvalovcin\Documents\GitHub\HSG-turbo\hsganalysis\UI\mainWin.ui'
#
# Created: Thu Jun 11 11:25:36 2015
#      by: PyQt4 UI code generator 4.10.4
#
# WARNING! All changes made in this file will be lost!

from PyQt5 import QtCore, QtWidgets

try:
    _fromUtf8 = QtCore.QString.fromUtf8
except AttributeError:
    def _fromUtf8(s):
        return s

try:
    _encoding = QtWidgets.QApplication.UnicodeUTF8
    def _translate(context, text, disambig):
        return QtWidgets.QApplication.translate(context, text, disambig, _encoding)
except AttributeError:
    def _translate(context, text, disambig):
        return QtWidgets.QApplication.translate(context, text, disambig)

class Ui_MainWindow(object):
    def setupUi(self, MainWindow):
        MainWindow.setObjectName(_fromUtf8("MainWindow"))
        MainWindow.resize(800, 600)
        self.centralwidget = QtWidgets.QWidget(MainWindow)
        self.centralwidget.setObjectName(_fromUtf8("centralwidget"))
        self.horizontalLayout = QtWidgets.QHBoxLayout(self.centralwidget)
        self.horizontalLayout.setObjectName(_fromUtf8("horizontalLayout"))
        self.tabWidget = QtWidgets.QTabWidget(self.centralwidget)
        self.tabWidget.setObjectName(_fromUtf8("tabWidget"))
        self.tabSpectrum = QtWidgets.QWidget()
        self.tabSpectrum.setObjectName(_fromUtf8("tabSpectrum"))
        self.horizontalLayout_2 = QtWidgets.QHBoxLayout(self.tabSpectrum)
        self.horizontalLayout_2.setObjectName(_fromUtf8("horizontalLayout_2"))
        self.splitSpectrum = QtWidgets.QSplitter(self.tabSpectrum)
        self.splitSpectrum.setOrientation(QtCore.Qt.Horizontal)
        self.splitSpectrum.setObjectName(_fromUtf8("splitSpectrum"))
        self.ptFile = ParameterTree(self.splitSpectrum)
        self.ptFile.setObjectName(_fromUtf8("ptFile"))
        self.gSpectrum = DraggablePlotWidget(self.splitSpectrum)
        self.gSpectrum.setObjectName(_fromUtf8("gSpectrum"))
        self.horizontalLayout_2.addWidget(self.splitSpectrum)
        self.tabWidget.addTab(self.tabSpectrum, _fromUtf8(""))
        self.tabFits = QtWidgets.QWidget()
        self.tabFits.setObjectName(_fromUtf8("tabFits"))
        self.horizontalLayout_3 = QtWidgets.QHBoxLayout(self.tabFits)
        self.horizontalLayout_3.setObjectName(_fromUtf8("horizontalLayout_3"))
        self.splitFits = QtWidgets.QSplitter(self.tabFits)
        self.splitFits.setOrientation(QtCore.Qt.Vertical)
        self.splitFits.setObjectName(_fromUtf8("splitFits"))
        self.gFits = DraggablePlotWidget(self.splitFits)
        self.gFits.setObjectName(_fromUtf8("gFits"))
        self.ptFits = ParameterTree(self.splitFits)
        self.ptFits.setObjectName(_fromUtf8("ptFits"))
        self.horizontalLayout_3.addWidget(self.splitFits)
        self.tabWidget.addTab(self.tabFits, _fromUtf8(""))
        self.horizontalLayout.addWidget(self.tabWidget)
        MainWindow.setCentralWidget(self.centralwidget)
        self.menubar = QtWidgets.QMenuBar(MainWindow)
        self.menubar.setGeometry(QtCore.QRect(0, 0, 800, 21))
        self.menubar.setObjectName(_fromUtf8("menubar"))
        MainWindow.setMenuBar(self.menubar)
        self.statusbar = QtWidgets.QStatusBar(MainWindow)
        self.statusbar.setObjectName(_fromUtf8("statusbar"))
        MainWindow.setStatusBar(self.statusbar)

        self.retranslateUi(MainWindow)
        self.tabWidget.setCurrentIndex(0)
        QtCore.QMetaObject.connectSlotsByName(MainWindow)

    def retranslateUi(self, MainWindow):
        MainWindow.setWindowTitle(_translate("MainWindow", "HSG Analysis Tool", None))
        self.tabWidget.setTabText(self.tabWidget.indexOf(self.tabSpectrum), _translate("MainWindow", "Spectrum", None))
        self.tabWidget.setTabText(self.tabWidget.indexOf(self.tabFits), _translate("MainWindow", "Fits", None))

from .draggablePlotWidget import DraggablePlotWidget
from pyqtgraph.parametertree import ParameterTree

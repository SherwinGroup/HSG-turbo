# -*- coding: utf-8 -*-
"""
Created on Mon Jun 01 16:01:43 2015

@author: dvalovcin
"""
import numpy as np
from PyQt4 import QtCore, QtGui
import pyqtgraph as pg
import pyqtgraph.parametertree.parameterTypes as pTypes
from pyqtgraph.parametertree import Parameter, ParameterTree, ParameterItem, registerParameterType
import os
import hsganalysis as hsg # local file, not global, no hsg.hsg

from UI.mainWin_ui import Ui_MainWindow

fileList = dict()


class MainWindow(QtGui.QMainWindow):
    """
    A note on parameter trees:
    They're not the most convenient things to navigate. You can't reference by name,
    you have to cheat, call their childs parameter, and know which index you're looking at

    if you change the order of the parameters at initialization, then things move and they aren'
    where you say they are. That's really fucking terrible. I could maybe subclass it
    to add my own features, which could be fun, but  I do not have time right now
    """
    def __init__(self, inp = None):
        super(MainWindow, self).__init__()
        self.initUI()
        self.__class__.dragEnterEvent = self.dragEnterEvent
        self.__class__.dragMoveEvent = self.dragEnterEvent
        self.__class__.dropEvent = self.drop
        self.setAcceptDrops(True)

        if inp is not None:
            if type(inp) is str: # a filename to open
                pass
            elif type(inp) is type(hsg.CCD):
                pass
            else:
                print "unknown type, ", type(inp)

        self.hsgObj = None

        import pyqtgraph.console as pgc
        self.consoleWindow = pgc.ConsoleWidget(namespace={"self": self, "np": np})
        self.consoleWindow.show()
        
    def initUI(self):
        self.ui = Ui_MainWindow()
        self.ui.setupUi(self)
        self.spectrumPlot = self.ui.gSpectrum.plot()
        
        params = self.genParameters()
        p = Parameter.create(name='Open a File', type='group', children=params)
        self.ui.ptFile.setParameters(p, showTop=True)
        
        self.sbLine = pg.InfiniteLine(pos = 750, movable=True)
        self.ui.gSpectrum.addItem(self.sbLine)
        self.sbLine.sigPositionChanged.connect(self.updateSBCalc)
        self.ui.splitSpectrum.setStretchFactor(1, 10)

        self.fitsPlot = self.ui.gFits.plot()
        self.ui.splitFits.setStretchFactor(0, 100)
        self.ui.splitFits.setStretchFactor(1, 1)



        self.menuSpec = QtGui.QMenu("Plot")

        self.menuSpecX = QtGui.QMenu("x")
        act = self.menuSpecX.addAction("nm")
        act.setCheckable(True)
        act.setChecked(True)
        act.triggered.connect(self.parseSpecXChange)
        act = self.menuSpecX.addAction("eV")
        act.setCheckable(True)
        act.triggered.connect(self.parseSpecXChange)
        act = self.menuSpecX.addAction("wavenumber")
        act.setCheckable(True)
        act.triggered.connect(self.parseSpecXChange)

        self.menuSpecY = QtGui.QMenu("y")
        act = self.menuSpecY.addAction("Log")
        act.setCheckable(True)
        act.triggered.connect(self.parseSpecYChange)

        self.menuSpec.addMenu(self.menuSpecX)
        self.menuSpec.addMenu(self.menuSpecY)



        self.menuFit = QtGui.QMenu("Plot")
        self.menuFitX = QtGui.QMenu("x")
        act = self.menuFitX.addAction("nm")
        act.setCheckable(True)
        act.setChecked(True)
        act.triggered.connect(self.parseFitXChange)
        act = self.menuFitX.addAction("eV")
        act.setCheckable(True)
        act.triggered.connect(self.parseFitXChange)
        act = self.menuFitX.addAction("wavenumber")
        act.setCheckable(True)
        act.triggered.connect(self.parseFitXChange)

        self.menuFitY = QtGui.QMenu("y")
        act = self.menuFitY.addAction("Log")
        act.setCheckable(True)
        act.triggered.connect(self.parseFitYChange)
        self.menuFitY.addSeparator()

        act = self.menuFitY.addAction("Height")
        act.setCheckable(True)
        act.triggered.connect(self.parseFitYChange)

        act = self.menuFitY.addAction("Sigma")
        act.setCheckable(True)
        act.triggered.connect(self.parseFitYChange)

        act = self.menuFitY.addAction("Position")
        act.setCheckable(True)
        act.triggered.connect(self.parseFitYChange)

        self.menuFit.addMenu(self.menuFitX)
        self.menuFit.addMenu(self.menuFitY)

        self.ui.tabWidget.currentChanged.connect(self.changeToolbars)
        # self.menuSpec.setVisible(False)
        self.ui.menubar.addMenu(self.menuSpec)
        # self.menuFit.setVisible(False)
        self.ui.menubar.addMenu(self.menuFit)
        self.changeToolbars()
        
        
        self.show()

    def changeToolbars(self):
        isSpec = self.ui.tabWidget.currentIndex() == 0
        self.menuFit.menuAction().setVisible(not isSpec)
        self.menuSpec.menuAction().setVisible(isSpec)

        # if self.ui.tabWidget.currentIndex() == 0:
        #     self.menuFit.menuAction().setVisible
        #     # self.ui.menubar.addMenu(self.menuSpec)
        #     self.ui.menubar.removeAction(self.menuFit)
        #     self.ui.menubar.addMenu(self.menuSpec)
        # else:
        #     self.ui.menubar.removeAction(self.menuSpec)
        #     self.ui.menubar.addMenu(self.menuFit)

    @staticmethod
    def __OPENTHINGS(): pass
        
    @staticmethod
    def genParameters(**kwargs):
        params = [
        {'name': "General Settings", "type":"group", "children":[
            {"name": "Exposure (s)", "type": "float", "value": kwargs.get("exposure", 0)},
            {"name": "Gain", "type": "int", "value": kwargs.get("gain", 0)},
            {"name": "CCD Temp (C)", "type": "str", "value": kwargs.get("ccd_temperature", 0)},
            {"name": "Spec Grating", "type": "int", "value": kwargs.get("grating", 0)},
            {"name": "Series", "type": "str", "value": kwargs.get("series", 0)},
            ]},
        {"name":"Laser Settings", "type":"group", "children":[
            {"name":"NIR Power (W)", "type":"float", "value":kwargs.get("nir_power", 0)},
            {"name":"NIR Frequency (nm)", "type":"float", "value":kwargs.get("nir_lambda", 0)},
            {"name":"Center Lambda (nm)", "type":"float", "value":kwargs.get("center_lambda", 0)}
            ]},
        {"name":"FEL Settings", "type":"group", "children":[
            {"name":"SB Number", "type":"float", "value":0},
            {"name":"FEL Energy (mJ)", "type":"float", "value":kwargs.get("fel_power", 0)},
            {"name":"FEL Frequency (cm-1)", "type":"float", "value":kwargs.get("fel_lambda", 0)},
            {"name":"Pulses", "type":"int", "value":kwargs.get("fel_pulses", 0)},
            {"name":"Pulse RR (Hz)", "type":"float", "value":kwargs.get("fel_reprate", 0)}
            ]}
        ]
        if "fieldInt" in kwargs:
            ptreelist = []
            for i, inten in enumerate(kwargs["fieldInt"]):
                d = {"name":"Pulse {:}".format(i), "type":"float", "value":inten}
                ptreelist.append(d)
            d = {"name":"FEL Intensity {:.3f} W/cm2".format(np.mean(kwargs["fieldInt"])),
                 "expanded":False, "type":"group", "children":ptreelist}
            params.append(d)
        if "fieldStrength" in kwargs:
            ptreelist = []
            for i, inten in enumerate(kwargs["fieldStrength"]):
                d = {"name":"Pulse {:}".format(i), "type":"float", "value":inten}
                ptreelist.append(d)
            d = {"name":"FEL Strength {:.3f} V/cm".format(np.mean(kwargs["fieldStrength"])),
                 "expanded":False, "type":"group", "children":ptreelist}
            params.append(d)

        return params

    def genFitParams(self, sbList):
        p = []
        for sb in sbList:
            d = {"name":"Sideband {}".format(sb[0]), "type":"group", "children":[
                {"name":"Peak Pos (eV): {}".format(sb[1]), "type":"group", "children":[
                    {"name":"nm", "type":"float", "value":1239.84/sb[1]},
                    {"name":"wavenumber", "type":"float", "value":10000000*sb[1]/1239.84}
                ], "expanded":False },
                {"name":"Height", "type":"float","value":sb[3]},
                {"name":"Sigma (1/eV): {}".format(sb[5]), "type":"group", "children":[
                    {"name":"1/nm", "type":"float", "value":sb[5]/1239.84},
                    {"name":"1/wavenumber", "type":"float", "value":1239.84/(10000000*sb[5])}
                ], "expanded":False }
            ]}
            p.append(d)
        return p



    def openFile(self, filename):
        filename = str(filename)
        self.setWindowTitle(os.path.basename(filename))
        data = np.genfromtxt(filename, delimiter=',', comments='#')
        # remove nans because genfromtxt doesn't read comments properly
        data = data[~np.all(np.isnan(data), axis=1)]

        with open(filename) as f:
            header = f.readline()
        header = header[1:]
        import json
        header = json.loads(header)

        self.ui.ptFile.clear()
        params = self.genParameters(**header)
        self.specParams = Parameter.create(name=filename, type='group', children=params)
        self.ui.ptFile.setParameters(self.specParams, showTop=True)

        hsgObj = hsg.CCD(filename)
        hsgObj = hsg.sum_spectra([hsgObj])[0]
        hsgObj.guess_better(cutoff = 3.5)
        hsgObj.fit_sidebands(plot=False)
        self.hsgObj = hsgObj

        self.plotSBFits(hsgObj)
        self.plotSpectrum(hsgObj.hsg_data)

        params = self.genFitParams(hsgObj.sb_results)
        self.fitParams = Parameter.create(name="Fit Results. FEL Freq: {}".format(self.calcFELFreq(hsgObj.sb_results)), type="group", children=params)
        self.ui.ptFits.setParameters(self.fitParams)







    @staticmethod
    def __PLOT_SPECT_THINGS(): pass


    def plotSBFits(self, hsgObj):
        for fit in hsgObj.sb_results:
            x = np.linspace(fit[1]*0.9995, fit[1]*1.0005, num=300)
            args = list(fit[1::2])
            args.append(0)
            y = hsg.gauss(x, *args)
            xType = [str(i.text()) for i in self.menuSpecX.actions() if i.isChecked()][0]
            if xType == "nm":
                x = 1239.84 / x
            elif xType == "wavenumber":
                x = 10000000*x/1239.84

            self.ui.gSpectrum.plot(x, y, pen=pg.mkPen('g'))

    def plotSpectrum(self, data):
        # Assumes data[:,0] is in eV
        x = data[:,0]
        xType = [str(i.text()) for i in self.menuSpecX.actions() if i.isChecked()][0]
        if xType == "nm":
            x = 1239.84 / x
        elif xType == "wavenumber":
            x = 10000000*x/1239.84
        print x
        self.ui.gSpectrum.plot(x, data[:,1])




    @staticmethod
    def __CHANGING_PLOT_CONTROLS(): pass

    def parseSpecXChange(self, val=None):
        sent = self.sender()
        if not val: # ignore if you uncheck
            sent.setChecked(True)
            return
        [i.setChecked(False) for i in self.menuSpecX.actions() if i is not sent] #uncheck the other one

        if self.hsgObj is not None: #reupdate plots if we've already loaded some data
            self.ui.gSpectrum.getPlotItem().clear()
            self.plotSpectrum(self.hsgObj.hsg_data)
            self.plotSBFits(self.hsgObj)


    def parseSpecYChange(self, val=None):
        self.ui.gSpectrum.getPlotItem().setLogMode(x=False, y=val)

    def parseFitYChange(self, val=None):
        sent = self.sender()
        if str(sent.text()) == "Log":
            self.ui.gFits.getPlotItem().setLogMode(x=False, y=val)
            return
        if not val: #untoggled an option
            sent.setChecked(True)

    def parseFitXChange(self, val=None):
        pass


    @staticmethod
    def __DRAGGING_CONTROLS(): pass
        
    def dragEnterEvent(self, event):
        event.accept()
        
    def drop(self, event):
        global fileList
        if not event.mimeData().hasUrls():
            event.reject()
            return
        filename = event.mimeData().urls()[0].toLocalFile()
        if not fileList: # no files have been opened
            self.openFile(filename)
        else:
            pass



    @staticmethod
    def __PROCESSING_CONTROLS(): pass

    def updateSBCalc(self):
        params = self.specParams.getValues()
        NIRL = float(params["Laser Settings"][1]["NIR Frequency (nm)"][0])
        FELL = float(self.specParams.childs[2].childs[2].value())

        laserL = 10000000./NIRL
        wantedWN = 10000000./self.sbLine.value()
        sbn = (wantedWN-laserL)/FELL

        self.specParams.childs[2].childs[0].setValue(sbn)

    def calcFELFreq(self, sbList):
        params = self.specParams.getValues()
        NIRL = float(params["Laser Settings"][1]["NIR Frequency (nm)"][0])
        laserL = 10000000./NIRL

        spacings = []

        for sb in sbList:
            peakPos = 10000000*sb[1]/1239.84
            spacings.append(round(peakPos-laserL)/sb[0])

        return np.mean(spacings)


        
        
        
        
        




if __name__=="__main__":
    import sys
    ex = QtGui.QApplication(sys.argv)
    win = MainWindow()
    sys.exit(ex.exec_())
    
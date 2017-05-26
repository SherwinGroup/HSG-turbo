# -*- coding: utf-8 -*-
"""
Created on Mon Jun 01 16:01:43 2015

@author: dvalovcin
"""
import os

import numpy as np
import pyqtgraph as pg
import pyqtgraph.parametertree.parameterTypes as pTypes
from PyQt4 import QtCore, QtGui
from pyqtgraph.parametertree import Parameter

# import hsganalysis as hsg # local file, not global, no hsg.hsg
import newhsganalysis as hsg

from UI.mainWin_ui import Ui_MainWindow
from UI.draggablePlotWidget import DraggablePlotWidget
try:
    import interactivePG as ipg
except:
    raise
# fileList = dict()
fileList = []
combinedWindowList = []

saveLoc = r"Z:\~Darren\Analysis\2017"
class ComboParameter(pTypes.WidgetParameterItem):
    # sigTreeStateChanged = QtCore.pyqtSignal(object, object)
    def __init__(self, param, depth, itemList, name=None):
        self.itemList = itemList
        if self.name is None:
            name = ""
        self._name = name
        super(ComboParameter, self).__init__(param, depth)
    def makeWidget(self):
        w = QtGui.QComboBox()
        [w.addItem(str(i)) for i in self.itemList]
        w.sigChanged = w.currentIndexChanged
        w.value = lambda: w.currentText()
        w.setValue = lambda x: w.setCurrentIndex(w.findText(str(x)))
        return w
    def name(self):
        return self._name
    def parentChanged(self, parent):
        pass


class BaseWindow(QtGui.QMainWindow):
    """
    A note on parameter trees:
    They're not the most convenient things to navigate. You can't reference by name,
    you have to cheat, call their childs parameter, and know which index you're looking at

    if you change the order of the parameters at initialization, then things move and they aren'
    where you say they are. That's really fucking terrible. I could maybe subclass it
    to add my own features, which could be fun, but  I do not have time right now
    """

    # Class to be instantiated with data
    dataClass = hsg.CCD
    dataSubclasses = []

    sigClosed = QtCore.pyqtSignal(object)
    def __init__(self, inp = None, parentWin = None):
        super(BaseWindow, self).__init__()
        self.initUI()
        # I must have grabbed these from stackoverflow
        # Sorry, don't really know why they're needed
        self.__class__.dragEnterEvent = self.dragEnterEvent
        self.__class__.dragMoveEvent = self.dragEnterEvent
        self.__class__.dropEvent = self.drop
        self.setAcceptDrops(True)
        self.dataObj = None
        self.titlePath = None

        # Handle being instantiated with various input
        if inp is not None:
            if type(inp) is str: # a filename to open
                self.openFile(inp)
            elif isinstance(inp, self.dataClass):
                self.processSingleData(inp)
            elif type(inp) in self.dataSubclasses:
                # to handle FullHighSideband class
                self.processSingleData(inp)
            else:
                print "unknown type, ", type(inp), self.__class__.__name__

        self.sigClosed.connect(updateFileClose)

        if parentWin is not None:
            self.inheritParent(parentWin)


    def initUI(self):
        self.ui = Ui_MainWindow()
        self.ui.setupUi(self)
        self.spectrumPlot = self.ui.gSpectrum.plot()

        # params = self.genParameters()
        # p = Parameter.create(name='Open a File', type='group', children=params)
        p = Parameter.create(name='Open a File', type='group')
        params = self.genParameters(parent=p)
        # p.addChildren(params)
        self.ui.ptFile.setParameters(p, showTop=True)
        self.specParams = p

        self.sbLine = pg.InfiniteLine(pos = 750, movable=True)
        self.ui.gSpectrum.addItem(self.sbLine)
        self.sbLine.sigPositionChanged.connect(self.updateSBCalc)
        self.ui.splitSpectrum.setStretchFactor(1, 10)

        self.fitsPlot = self.ui.gFits.plot(fmt='o')
        self.ui.splitFits.setStretchFactor(0, 100)
        self.ui.splitFits.setStretchFactor(1, 1)


        self.ui.gFits.plotItem.vb.sigDropEvent.connect(self.handleFitDragEvent)
        self.ui.gSpectrum.plotItem.vb.sigDropEvent.connect(self.handleSpecDragEvent)

        self.curveSpectrum = self.ui.gSpectrum.plot(pen='k', name='Spectrum')



        self.menuSpec = QtGui.QMenu("Plot")


        #################################
        #
        # Spectrum menu, x-axis
        #
        #################################
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


        #################################
        #
        # Spectrum menu, y-axis
        #
        #################################

        self.menuSpecY = QtGui.QMenu("y")
        act = self.menuSpecY.addAction("Log")
        act.setCheckable(True)
        act.triggered.connect(self.parseSpecYChange)

        self.menuSpecScaleY = QtGui.QMenu("Divide By")

        self.uisbDivideBy = pg.SpinBox(value=1, step=1.0)
        self.uisbDivideByAction = QtGui.QWidgetAction(None)
        self.uisbDivideByAction.setDefaultWidget(self.uisbDivideBy)
        self.menuSpecScaleY.addAction(self.uisbDivideByAction)
        self.menuSpecY.addMenu(self.menuSpecScaleY)
        self.uisbDivideBy.sigValueChanged.connect(self.scaleData)

        self.menuSpec.addMenu(self.menuSpecX)
        self.menuSpec.addMenu(self.menuSpecY)



        self.menuFit = QtGui.QMenu("Plot")


        #################################
        #
        # Fit menu, x-axis
        #
        #################################
        self.menuFitX = QtGui.QMenu("x")
        act = self.menuFitX.addAction("SB Num")
        act.setCheckable(True)
        act.setChecked(True)
        act.triggered.connect(self.parseFitXChange)
        act = self.menuFitX.addAction("nm")
        act.setCheckable(True)
        act.triggered.connect(self.parseFitXChange)
        act = self.menuFitX.addAction("eV")
        act.setCheckable(True)
        act.triggered.connect(self.parseFitXChange)
        act = self.menuFitX.addAction("wavenumber")
        act.setCheckable(True)
        act.triggered.connect(self.parseFitXChange)





        #################################
        #
        # Fit menu, y-axis
        #
        #################################

        self.menuFitY = QtGui.QMenu("y")
        act = self.menuFitY.addAction("Log")
        act.setCheckable(True)
        act.triggered.connect(self.parseFitYChange)

        # self.menuFitScaleY = QtGui.QMenu("Divide By")
        # self.uisbDivideByF = pg.SpinBox()
        # self.uisbDivideByActionF = QtGui.QWidgetAction(None)
        # self.uisbDivideByActionF.setDefaultWidget(self.uisbDivideByF)
        # self.menuFitScaleY.addAction(self.uisbDivideByActionF)
        # self.menuFitY.addMenu(self.menuFitScaleY)

        self.menuFitY.addSeparator()

        act = self.menuFitY.addAction("Height")
        act.setCheckable(True)
        act.setChecked(True)
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
        self.menuBar().setNativeMenuBar(False) # otherwise the below breeaks on OS X
        # self.menuSpec.setVisible(False)
        self.ui.menubar.addMenu(self.menuSpec)

        # self.menuFit.setVisible(False)
        self.ui.menubar.addMenu(self.menuFit)
        self.changeToolbars()


        #################################
        #
        # setup new thing
        #
        #################################
        ret = self.makeTitleActionList()
        self.ui.menubar.addMenu(ret)
        self.ret = ret


        save = self.ui.menubar.addAction("Save Processed...")
        save.triggered.connect(self.saveProcessed)

        self.show()

    def changeToolbars(self):
        isSpec = self.ui.tabWidget.currentIndex() == 0
        self.menuFit.menuAction().setVisible(not isSpec)
        self.menuSpec.menuAction().setVisible(isSpec)

    def inheritParent(self, parent):

        # set the same window title
        if parent.titlePath is not None:
            self.updateTitle(path = parent.titlePath)

    @staticmethod
    def __OPENTHINGS(): pass

    def genParameters(self, **kwargs):
        params = [
        {'name': "General Settings", "type":"group", "children":[
            {"name": "Date", "type": "str", "value": kwargs.get("date", "01/01/79"), "readonly":True},
            {"name": "Series", "type": "str", "value": kwargs.get("series", 0), "readonly":True},
            {"name": "Comments", "type": "text", "value": kwargs.get("comments", "NA"), "readonly":True},
            {"name": "Exposure (s)", "type": "float", "value": kwargs.get("exposure", 0), "readonly":True},
            {"name": "Gain", "type": "int", "value": kwargs.get("gain", 0), "readonly":True},
            {"name": "CCD Image", "type": "str", "value": "[{}, {}, {}]".format(
                *np.array(kwargs.get("ccd_image_settings", [-1]*6))[[1, 4, 5]]), "readonly":True},
            {"name": "CCD Temp (C)", "type": "str", "value": kwargs.get("ccd_temperature", 0), "readonly":True},
            {"name": "Spec Grating", "type": "int", "value": kwargs.get("grating", 0), "readonly":True},
            {"name": "Spec Sweep", "type": "int", "value": kwargs.get("spec_step", 0), "readonly":True},
            {"name":"Center Lambda (nm)", "type":"float", "value":kwargs.get("center_lambda", 0), "readonly":True},
            {"name":"Offset (nm)", "type":"float", "value":kwargs.get("offset", 0)}
            ]},
        ]
        try:
            if self.__class__ == BaseWindow:
                kwargs["parent"].addChildren(params)
                kwargs["parent"].child("General Settings", "Offset (nm)").sigValueChanged.connect(
                    self.updateOffset
                )
        except KeyError:
            pass
        return params

    def openFile(self, filename):
        filename = str(filename)

        dataObj = self.dataClass(filename)


        if dataObj.parameters["series"]:
            self.setWindowTitle(dataObj.parameters["series"])
        else:
            self.setWindowTitle(os.path.basename(filename))

        self.ui.ptFile.clear()
        # params = self.genParameters(**dataObj.parameters)
        # self.specParams = Parameter.create(name=filename, type='group', children=params)
        # self.ui.ptFile.setParameters(self.specParams, showTop=True)

        self.processSingleData(dataObj)
        self.ui.ptFile.setToolTip(filename)

    def processSingleData(self, dataObj):
        self.dataObj = dataObj
        self.plotSpectrum()
        p = Parameter.create(
            name=dataObj.fname,
            type='group')
        params = self.genParameters(parent=p, **dataObj.parameters)
        self.ui.ptFile.setParameters(p, showTop=True)


    @staticmethod
    def getWindowClass(fname):
        if 'hsg' in fname:
            if "Images" in fname:
                return HSGImageWindow
            return HSGWindow
        elif 'abs' in fname:
            return AbsWindow
        elif 'pl' in fname:
            return PLWindow
        else:
            return BaseWindow

    @staticmethod
    def getDataClass(fname):
        if 'hsg' in fname:
            return hsg.HighSidebandCCD
        elif 'abs' in fname:
            return hsg.Absorbance
        elif 'pl' in fname:
            return hsg.Photoluminescence
        else:
            return hsg.CCD

    @staticmethod
    def __PLOT_SPECT_THINGS(): pass

    def plotSpectrum(self):
        self.curveSpectrum.setData(*self.convertDataForPlot(self.dataObj.proc_data))
    def convertDataForPlot(self, data):
        """
        Return the proper units on the x-values to be used for plotting.
        Takes the desired values from the GUI selection.
        :param data:
        :return:
        """
        # Assumes data[:,0] is in eV
        x = data[:,0].copy()
        xType = [str(i.text()) for i in self.menuSpecX.actions() if i.isChecked()][0]
        x = converter["eV"][xType](x)

        return [x, data[:,1]/self.uisbDivideBy.value()]

    @staticmethod
    def __PLOTTING_FIT_RESUTLS(): pass
    def plotFits(self):
        data = self.dataObj.sb_results
        wantX = [str(i.text()) for i in self.menuFitX.actions() if i.isChecked()][0]
        if wantX == "SB Num":
            x = data[:,0]
        else:
            x = data[:,1]
            x = converter["eV"][wantX](x)

        wantY = [str(i.text()) for i in self.menuFitY.actions() if
                 i.isChecked() and str(i.text()) != "Log"][0]
        y = data[:, {"Height":   3,
                     "Sigma":    5,
                     "Position": 1}[wantY]]


        self.fitsPlot.setData(x, y/self.uisbDivideBy.value())

    @staticmethod
    def __CHANGING_PLOT_CONTROLS(): pass

    def parseSpecXChange(self, val=None):
        sent = self.sender()
        if not val: # ignore if you uncheck
            sent.setChecked(True)
            return
        oldName = [str(i.text()) for i in self.menuSpecX.actions() if i is not sent and i.isChecked()][0]
        newName = str(sent.text())
        [i.setChecked(False) for i in self.menuSpecX.actions() if i is not sent] #uncheck the other one

        if self.dataObj is not None: #reupdate plots if we've already loaded some data
            sb = self.sbLine.value()
            # self.ui.gSpectrum.getPlotItem().clear()

            self.plotSpectrum()
            self.plotSBFits()


            sb = converter[oldName][newName](sb)
            self.sbLine.setValue(sb)

    def parseSpecYChange(self, val=None):
        self.ui.gSpectrum.getPlotItem().setLogMode(x=False, y=val)

    def parseFitYChange(self, val=None):
        sent = self.sender()
        if str(sent.text()) == "Log":
            self.ui.gFits.getPlotItem().setLogMode(x=False, y=val)
            return
        if not val: #untoggled an option
            sent.setChecked(True)

        # uncheck what was previously checked
        [i.setChecked(False) for i in self.menuFitY.actions() if i is not sent and str(i.text())!="Log"]
        if self.dataObj is not None:
            self.ui.gFits.plotItem.clear()
            self.plotFits()

    def parseFitXChange(self, val=None):
        sent = self.sender()
        if not val: # ignore if you uncheck
            sent.setChecked(True)
            return
        oldName = [str(i.text()) for i in self.menuFitX.actions() if i is not sent and i.isChecked()][0]
        newName = str(sent.text())
        [i.setChecked(False) for i in self.menuFitX.actions() if i is not sent]
        if self.dataObj is not None:
            self.plotFits()

    def scaleData(self, sb):
        if self.dataObj is not None:
            # self.ui.gSpectrum.getPlotItem().clear()
            # self.ui.gSpectrum.plotItem.removeItem(self.sbLine)
            self.plotSpectrum()
            self.plotSBFits()
            # self.ui.gSpectrum.addItem(self.sbLine)

    def updateOffset(self, paramObj, val):
        if self.dataObj is None: return
        self.dataObj.parameters["offset"] = val

        # convert back to nm, add the offset, reconvert
        self.dataObj.ccd_data[:, 0] = 1239.84/ (self.dataObj.raw_data[:,0] + val)

        if hasattr(self.dataObj, "proc_data"):
            self.dataObj.proc_data[:, 0] = 1239.84 / (self.dataObj.raw_data[:, 0] + val)

        self.processSingleData(self.dataObj)
        if self.titlePath is not None:
            self.updateTitle(path=self.titlePath)


    @staticmethod
    def __CHANGING_WINDOW_TITLE(): pass
    def makeTitleActionList(self, group=None, actions=None):
        if group is None:
            group = self.specParams
            actions = QtGui.QMenu("Set Window Title")
            a = actions.addAction("Filename")
            a.triggered.connect(self.updateTitle)
            a = actions.addAction("Other...")
            a.triggered.connect(self.updateTitle)


        for ch in group.children():
            if isinstance(ch, pTypes.GroupParameter):
                self.makeTitleActionList(group = ch, actions=actions.addMenu(ch.name()))
            else:
                act = actions.addAction(ch.name())
                act.triggered.connect(self.updateTitle)
        return actions

    def updateTitle(self, bool=True, path=None):
        """
        Update the window title, either called from signals
        from the menu (where bool=False is due to the signals
        being sent with a value)
        Can be called direclty
        :param bool:
        :param path: list of heirarchy to set title value to
        :return:
        """
        if path is None:
            child = self.sender()
            if child.parentWidget().title() == "Set Window Title":
                if child.text() == "Filename" and self.dataObj is not None:
                    self.setWindowTitle(os.path.basename(self.dataObj.fname))
                elif child.text() == "Other...":
                    text, ok = QtGui.QInputDialog.getText(self,
                                                      "Enter Window Name",
                                                      "Title:",
                                                      QtGui.QLineEdit.Normal,
                                                      self.windowTitle())
                    if ok: self.setWindowTitle(text)
                self.path=None
                return
            path = [str(child.parentWidget().title()), str(child.text())]

        self.titlePath = path
        params = self.specParams
        try:
            val = params.child(*path).value()
        except Exception as e:
            print "ERROR updating window title"
            print "\tbool=",bool
            print "path=",path
        name = path[1]
        # cut out units if they're present
        pref = name[:name.find("(")-1] if "(" in name else name

        self.setWindowTitle("{}: {}".format(pref, val))

    @staticmethod
    def __DRAGGING_CONTROLS(): pass

    def dragEnterEvent(self, event):
        """

        :param event:
        :type event: QtGui.QDragEnterEvent
        :return:
        """
        event.accept()

    def drop(self, event):
        """
        :param event:
        :type event: QtGui.QDropEvent
        :return:
        """
        global fileList
        # Dropped something that isn't a link
        if not event.mimeData().hasUrls():
            event.reject()
            return

        # Force Qt to "Copy" the file instead of
        # "Move", preventing it from removing the
        # file from the directory.
        event.setDropAction(QtCore.Qt.CopyAction)

        if event.keyboardModifiers() & QtCore.Qt.ShiftModifier:
            print "held shift"


        # Only one file was dropped
        if len(event.mimeData().urls()) == 1:
            filename = str(event.mimeData().urls()[0].toLocalFile())
            c = BaseWindow.getWindowClass(filename)
            try:
                a = c(filename, parentWin=self)
            except Exception as e:
                raise
            else:
                fileList.append(a)
        # Multiple files were dropped
        else:
            # Make a list of all of them, cuttingout the "seriesed" ones
            # Because we'll use the HSGanalysis methods to handle it
            filelist = [str(i.toLocalFile()) for i in event.mimeData().urls()]
            filelist = [i for i in filelist if "seriesed" not in i.lower()]
            # Make a CCD obj of each file
            objlist = [BaseWindow.getDataClass(i)(i) for i in filelist]
            if np.all(["hsg" in ii.fname for ii in objlist]):
                # Sum them all up if every file is hsg, otherwise
                # do nothing
                # series = hsg.hsg_sum_spectra(objlist)
                series = hsg.proc_n_plotCCD(filelist)
                # series = hsg.hsg_combine_spectra(series)
            else:
                series = objlist


            for obj in series:
                # If there's no series tag, it was a single
                # spectrum. Setting the window title
                # to the spectrometer wavelength is kinda random
                if obj.parameters["series"] == "":
                    c = BaseWindow.getWindowClass(obj.fname)
                    fileList[obj.fname] = c(obj, parentWin=self)
                    if self.titlePath is None:
                        # default to center for window title unless we've been messing around
                        fileList[obj.fname].setWindowTitle(str(obj.parameters["center_lambda"]))
                else:
                    c = BaseWindow.getWindowClass(obj.fname)
                    a = c(obj, parentWin=self)
                    fileList.append(a)
                    if self.titlePath is None:
                        # default to series for window title unless we've been messing around
                        a.setWindowTitle("Series: {}".format(obj.parameters["series"]))
        if self.dataObj is None: # I'm the first window. Get outta here!
            self.close()

    def handleFitDragEvent(self, obj, val):
        """
        called when the plot of fit results is dragged/dropped
        :param obj: The thing dragged
        :param val:  the pyqtgraph coordinate of the drop point
        :return:
        """
        # Get the xy data from the plot directly,
        # may break if pyqtgraph changes indexing of curves
        # I think the other indices include axes and things?
        # could chagne this by making the data a class
        # member, but this keeps it clear what units are used
		
        if self.dataObj is None:
            return
        # data = self.dataObj.sb_results
        # want = [str(i.text()) for i in self.menuFitY.actions() if i.isChecked() and str(i.text())!="Log"][0]
        # x = data[:,0]
        # y = data[:, {"Height":   3,
        #              "Sigma":    5,
        #              "Position": 1}[want]]
        #
        # #self.ui.gFits.plot(x, y/self.uisbDivideBy.value(), pen=pg.mkPen("w"), symbol="o")
        #
        # #d = [self.ui.gFits.plotItem.curves[2].xData,
        # #     self.ui.gFits.plotItem.curves[2].yData]
        # d = [x, y/self.uisbDivideBy.value()]
        d = self.fitsPlot.getData()
        self.createCompWindow(data = d, p = val)

    def handleSpecDragEvent(self, obj, val):
        """
        See comments on handleFitDragEvent
        :param obj:
        :param val:
        :return:
        """
        # d = [self.ui.gSpectrum.plotItem.curves[1].xData,
        #      self.ui.gSpectrum.plotItem.curves[1].yData]
        self.createCompWindow(data = self.convertDataForPlot(self.dataObj.proc_data), p = val)

    def createCompWindow(self, data, p, label=None):
        # If there's already open windows, see if the drop point
        # is over an already created window. If it is,
        # add it to that window
        if combinedWindowList:
            inners = [i for i in combinedWindowList if i.containsPoint(p.toQPoint())]
            if inners: #tests if it's empty
                a = inners[-1] # last focused window
                # Data may not be passed if you make a comparision window
                # Before loading any files, don't throw errors
                if not None in data:
                    a.addCurve(data, str(self.windowTitle()))
                return

        # If there's no window or not dragged on top of something,
        # create a new window
        a = ComparisonWindow()

        # Connect the signals for when it gets closed or focused
        a.sigClosed.connect(updateCompClose)
        a.sigGotFocus.connect(updateFocusList)

        # Data may not be passed if you make a comparision window
        # Before loading any files, don't throw errors
        if not None in data:
            a.addCurve(data, str(self.windowTitle()))
        # Move the window to the drag point
        a.move(p.toQPoint() - QtCore.QPoint(a.geometry().width()/2, a.geometry().height()/2))
        # Add it to the list of opened windows
        combinedWindowList.append(a)

    @staticmethod
    def __PROCESSING_CONTROLS(): pass

    def updateSBCalc(self):
        params = self.specParams
        # see notes in calcFELFreq for explaination
        # of following two lines
        NIRL = float(
            params.child("Laser Settings", "NIR Frequency").opts["values"][0].split()[0])
        FELL = float(
            params.child("FEL Settings", "FEL Frequency (cm-1)").opts["value"])
        rawVal = self.sbLine.value()
        units = [str(i.text()) for i in self.menuSpecX.actions() if i.isChecked()][0]
        linenm = converter[units]["nm"](rawVal)



        laserL = 10000000./NIRL
        wantedWN = 10000000./linenm
        sbn = (wantedWN-laserL)/FELL

        self.specParams.child("General Settings", "SB Number").setValue(sbn)

    def saveProcessed(self):
        global saveLoc
        print "Saved Clicked"
        if not self.dataObj:
            print "Load a file first"
        path = QtGui.QFileDialog.getSaveFileName(self, "Save File", saveLoc, "Text File (*.txt)")
        if not path:
            return
        path = str(path)
        saveLoc = path
        self.dataObj.save_processing(os.path.basename(path), os.path.dirname(path))

    def calcFELFreq(self, sbList):
        params = self.specParams
        #            v  - gets the child
        laserL = float(
            params.child("Laser Settings", "NIR Frequency").opts["values"][1].split()[0])
        #            ^                                       ^             ^    ^ split the unit off
        #            |                                       |             |
        #            |                                       |             L- grab the cm-1 value
        #            |                                       L- internal dict to grab the values of the list
        #            L- search the param tree for NIr Freq under the laser settings
        if laserL==10000000.:
            # default value whe windows is first made
            return -1

        spacings = []

        for sb in sbList:
            peakPos = 10000000*sb[1]/1239.84
            spacings.append(round(peakPos-laserL)/sb[0])

        return np.mean(spacings)

    @staticmethod
    def __CLEANING_UP(): pass

    def closeEvent(self, event):
        self.sigClosed.emit(self)
        super(BaseWindow, self).closeEvent(event)


class HSGWindow(BaseWindow):
    dataClass = hsg.HighSidebandCCD
    dataSubclasses = [hsg.FullHighSideband]
    def __init__(self, *args, **kwargs):
        self.curveFits = {}
        super(HSGWindow, self).__init__(*args, **kwargs)


    def inheritParent(self, parent):
        if isinstance(parent, HSGWindow):
            # set the same units on NIR frequency
            params = parent.specParams
            idx = params.child("Laser Settings", "NIR Frequency").opts["values"].index(
                params.child("Laser Settings", "NIR Frequency").opts["value"]
            )
            params = self.specParams
            params.child("Laser Settings", "NIR Frequency").setValue(
                params.child("Laser Settings", "NIR Frequency").opts["values"][idx]
            )
            self.ui.tabWidget.setCurrentIndex(parent.ui.tabWidget.currentIndex())
            [i.setChecked(False) for i in self.menuFitX.actions()]
            [i.setChecked(True) for (i,j) in
                zip(self.menuFitX.actions(), parent.menuFitX.actions()) if j.isChecked()]
            self.plotFits()

            # [i.setChecked(False) for i in self.menuSpecX.actions()]
            # [i.setChecked(True) for (i,j) in
            #     zip(self.menuSpecX.actions(), parent.menuSpecX.actions()) if j.isChecked()]
            # self.plotFits()
        super(HSGWindow, self).inheritParent(parent)

    def processSingleData(self, dataObj):
        super(HSGWindow, self).processSingleData(dataObj)

        self.dataObj.guess_sidebands()
        self.dataObj.fit_sidebands()
        self.dataObj.infer_frequencies()


        self.specParams = Parameter.create(
            name=dataObj.fname,
            type='group')
        params = self.genParameters(parent=self.specParams, **dataObj.parameters)
        self.ui.ptFile.setParameters(self.specParams, showTop=True)

        self.specParams.child("General Settings", "SB Number").sigValueChanging.connect(
            self.updateSBLine
        )

        self.curveFits = {ii: self.ui.gSpectrum.plot(pen='g', name=ii) for
                          ii in self.dataObj.full_dict}


        self.plotSBFits()

        params = self.genFitParams(dataObj.sb_results)
        self.fitParams = Parameter.create(
            name="Fit Results. FEL Freq: {}".format(self.calcFELFreq(dataObj.sb_results)),
            type="group", children=params)
        self.ui.ptFits.setParameters(self.fitParams)

        self.plotFits()

    def plotSBFits(self):
        xType = [str(i.text()) for i in self.menuSpecX.actions() if i.isChecked()][0]
        units = converter["eV"][xType]

        # for sb, curve in self.curveFits.items():
        for sb in self.curveFits:
            fit = self.dataObj.full_dict[sb]
            curve = self.curveFits[sb]

            x = np.linspace(fit[0]-5*fit[4], fit[0]+5*fit[4], num=300)
            args = list(fit[0::2])
            args.append(0)
            y = hsg.gauss(x, *args)
            curve.setData(units(x), y/self.uisbDivideBy.value())

    def updateSBLine(self, paramObj, val):
        # print "updating line", self.sender(), args, kwargs

        params = self.specParams
        # see notes in calcFELFreq for explaination
        # of following two lines
        NIRL = float(
            params.child("Laser Settings", "NIR Frequency").opts["values"][1].split()[0])
        FELL = float(
            params.child("FEL Settings", "FEL Frequency (cm-1)").opts["value"])

        units = [str(i.text()) for i in self.menuSpecX.actions() if i.isChecked()][0]

        sbWN = NIRL + FELL*val
        line = converter["wavenumber"][units](sbWN)

        paramObj.blockSignals(True)
        self.sbLine.setValue(line)
        paramObj.blockSignals(False)

    def genParameters(self, **kwargs):
        if kwargs.get("nir_lambda", None) is not None:
            kwargs["nir_lambda"] = float(kwargs["nir_lambda"])
        params = super(HSGWindow, self).genParameters(**kwargs)

        # print "22recasting NIR_LAMBDA"
        # kwargs["nir_lambda"] = float(kwargs["nir_lambda"])
        # print 22, kwargs["nir_lambda"], type(kwargs["nir_lambda"])
        params[-1]["children"].append(
            {"name":"SB Number", "type":"float", "value":0, "step":0.05}
        )
        params.append(
        {"name":"Sample Parameters", "type":"group", "children":[
            {"name":"Sample", "type":"str", "value":kwargs.get("sample_name", "None?"), "readonly":True},
            {"name":"Sample Temp", "type":"float", "value":kwargs.get("sample_Temp", -1), "readonly":True}
            ]})
        params.append(
        {"name":"Laser Settings", "type":"group", "children":[
            {"name":"NIR Power (mW)", "type":"float", "value":kwargs.get("nir_power", 0), "readonly":True},
            {"name":"NIR Frequency", "type":"list",
                "values":["{:.3f} nm".format(kwargs.get("nir_lambda", 0)),
                          "{:.2f} cm-1".format(1e7/kwargs.get("nir_lambda", 1))]},
            {"name":"NIR Polarization", "type":"str", "value":kwargs.get("nir_pol", "?"), "readonly":True},
            {"name": "Detector Pol", "type": "float",
                "value": kwargs.get("detectorHWP", "-1"), "readonly":True},
            {"name": "Fit NIR (cm-1)", "type": "float",
                "value": kwargs.get("calculated NIR freq (cm-1)", 0), "readonly":True},
            {"name": "Fit THz (cm-1)", "type": "float",
                "value": kwargs.get("calculated THz freq (cm-1)", 0), "readonly":True}
            ]})
        params.append(
        {"name":"FEL Settings", "type":"group", "children":[
            {"name":"Pulses", "type":"int", "value":np.mean(kwargs.get("fel_pulses", 0)), "readonly":True},
            {"name":"FEL Energy (mJ)", "type": "str",
                "value": "{:.1f} +/- {:.1f}".format(
                    kwargs.get("pulseEnergies", {"mean":-1})["mean"],
                    kwargs.get("pulseEnergies", {"std": -1})["std"]), "readonly":True},
            {"name":"Field Strength (kV/cm)", "type": "str",
                "value": "{:.2f} +/- {:.2f}".format(
                    kwargs.get("fieldStrength", {"mean":-1})["mean"],
                    kwargs.get("fieldStrength", {"std": -1})["std"]), "readonly":True},
            {"name": "Transmission", "type": "float", "value": kwargs.get("fel_transmission", 0), "readonly":True},
            {"name": "Pyro Voltage", "type": "str",
                "value": "{:.1f} +/- {:.1f} mV".format(
                    kwargs.get("pyroVoltage", {"mean":-1})["mean"]*1e3,
                    kwargs.get("pyroVoltage", {"std": -1})["std"]*1e3), "readonly":True},
            {"name": "CD Ratio", "type": "str",
                "value": "{:.1f} +/- {:.1f}".format(
                    kwargs.get("cdRatios", {"mean":-1})["mean"]*1e2,
                    kwargs.get("cdRatios", {"std": -1})["std"]*1e2), "readonly":True},
            {"name": "FP Time", "type": "str",
                "value": "{:.0f} +/- {:.0f} ns".format(
                    kwargs.get("fpTime", {"mean":-1})["mean"]*1e3,
                    kwargs.get("fpTime", {"std": -1})["std"]*1e3), "readonly":True},
            {"name":"Pulse RR (Hz)", "type":"float", "value":kwargs.get("fel_reprate", 0), "readonly":True},
            {"name":"FEL Frequency (cm-1)", "type":"float", "value":kwargs.get("fel_lambda", 0), "readonly":True}
            ]})
        try:
            kwargs["parent"].addChildren(params)
            kwargs["parent"].child("General Settings", "Offset (nm)").sigValueChanged.connect(
                self.updateOffset
            )
        except KeyError:
            raise
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

    def updateOffset(self, paramObj, val):
        [self.ui.gSpectrum.getPlotItem().removeItem(ii) for ii in self.curveFits.values()]
        super(HSGWindow, self).updateOffset(paramObj, val)
        self.updateSBCalc()


class AbsWindow(BaseWindow):
    dataClass = hsg.Absorbance

class PLWindow(BaseWindow):
    dataClass = hsg.Photoluminescence

penList = [pg.intColor(i, 20) for i in range(20)]

class ComparisonWindow(QtGui.QMainWindow):
    sigClosed = QtCore.pyqtSignal(object) #emit oneself when closed
    sigGotFocus = QtCore.pyqtSignal(object)
    def __init__(self, parent = None,  data = None):
        super(ComparisonWindow, self).__init__()
        self.initUI()
        self.curveList = {} # a way to keep track of what's been added

        self.selectedList = [] # List of which lines have been selected

        self.setFocusPolicy(QtCore.Qt.StrongFocus)
        self.gPlot.setFocusPolicy(QtCore.Qt.StrongFocus)

        self.gPlot.focusInEvent = self.focusInEvent
        self.gPlot.changeEvent = self.changeEvent

        # self.gPlot.plotItem.vb.sigClickedEvent.connect(self.handleMouseClick)

        #line = pg.LineSegmentROI()

    def initUI(self):
        # self.gPlot = DraggablePlotWidget()
        self.gPlot = ipg.PlotWidget()
        self.gPlot.addLegend()
        self.setCentralWidget(self.gPlot)
        # self.legend = pg.LegendItem()
        # self.legend.setParentItem(self.gPlot.plotItem)
        self.menuBar().setNativeMenuBar(False)

        # removeItems = QtGui.QAction("Remove Selected Items", self.menuBar())
        # removeItems.triggered.connect(self.removeSelectedLines)
        # self.menuBar().addAction(removeItems)

        # changeColor = QtGui.QAction("Change Line Color", self.menuBar())
        # changeColor.triggered.connect(self.changeLineColor)
        # self.menuBar().addAction(changeColor)

        self.show()

    def containsPoint(self, p):
        """
        Calculates whether a specified QPoint (in abs coords) is
        within the bounds of myself
        :param p: the q point
        :return: True if it is within my bounds, else false
        """
        return self.frameGeometry().contains(p)

    def addCurve(self, data, label=None):
        symbol = None
        if len(data[0])<100:
            symbol = 'o'
        # if label is not None:
            # self.legend.addItem(p, label)
        color = pg.intColor(len(self.curveList), hues=20)
        p = self.gPlot.plot(data[0], data[1],
                                     pen=color,
                                     symbol=symbol,
                                    symbolPen=color,
                                    symbolBrush=color,
                                    name=label)
        self.curveList[p] = label
        # p = self.gPlot.plotItem.plot(data[0], data[1],
        #                              pen=pg.mkPen(pg.intColor(len(self.curveList), hues=20)),
        #                              symbol=symbol)
        # p.curve.setClickable(True)
        # p.sigClicked.connect(self.handleMouseClick)

    def removeSelectedLines(self):
        for line in self.selectedList:
            # Want to remove the line from both the internal
            # curveList and the legend
            self.legend.removeItem(self.curveList.pop(line))
            self.gPlot.plotItem.removeItem(line)
            self.selectedList.remove(line)
        self.legend.updateSize()

    def closeEvent(self, event):
        self.sigClosed.emit(self)
        super(ComparisonWindow, self).closeEvent(event)

    def focusInEvent(self, *args, **kwargs):
        print "got focus", args, kwargs
        self.sigGotFocus.emit(self)

    def changeEvent(self, ev):
        if ev.type() == QtCore.QEvent.ActivationChange:
            ev.accept()
            if self.isActiveWindow():
                self.sigGotFocus.emit(self)
        else:
            ev.ignore()
            ipg.PlotWidget.changeEvent(self.gPlot, ev)

    def handleMouseClick(self, obj, pos = None):
        """
        Handle highlighting curves when they're clicked
        :param obj: a PlotCurveitem if a line is selected, else
                    the ViewBox item of the gPlot
        :param pos: If line selected, None
                    Else, the position of the click
        :return:
        """
        # If a line was clicked
        if pos is None:
            # Either select or deselect,
            # adding/removing from the list
            # and setting the width as appropriate
            if obj in self.selectedList:
                width = 1
                self.selectedList.remove(obj)
            else:
                width = 3
                self.selectedList.append(obj)
            # Get the pen so that everything is the same
            # And we just change the width, and set it back
            pen = obj.opts['pen']
            pen.setWidth(width)
            obj.setPen(pen)
        # Clicked not a curve, so un-bold all the plots
        else:
            for obj in self.curveList.keys():
                pen = obj.opts['pen']
                pen.setWidth(1)
                obj.setPen(pen)
            self.selectedList = []
        # Force the legend to update to redraw with
        # the new pen, too
        self.legend.updateSize()

    def changeLineColor(self):
        if not self.selectedList:
            return
        color = QtGui.QColorDialog.getColor()
        p = self.selectedList[-1]
        pen = p.opts["pen"]
        pen.setColor(color)
        p.setPen(pen)

class HSGImageWindow(HSGWindow):
    # dataClass = object
    def __init__(self, *args, **kwargs):
        super(HSGImageWindow, self).__init__(*args, **kwargs)

    def inheritParent(self, parent):
        if isinstance(parent, HSGImageWindow):
            self.ui.gSpectrum.setLevels(*parent.ui.gSpectrum.ui.histogram.getLevels())


    def initUI(self):
        self.ui = Ui_MainWindow()
        self.ui.setupUi(self)

        # Set arbitrary parent for garbage colection
        self.ui.gSpectrum.setParent(QtGui.QWidget())

        # add image plot to it.
        self.ui.gSpectrum = ipg.ImageView(self.ui.splitSpectrum)
        self.ui.gSpectrum.view.setAspectLocked(False)
        self.ui.splitSpectrum.setStretchFactor(1, 10)
        p = Parameter.create(name='Open a File', type='group')
        params = self.genParameters(parent=p)
        # p.addChildren(params)
        self.ui.ptFile.setParameters(p, showTop=True)
        self.specParams = p
        ret = self.makeTitleActionList()
        self.ui.menubar.addMenu(ret)
        self.ret = ret


        self.show()

    def processSingleData(self, dataObj):
        super(HSGWindow, self).processSingleData(dataObj)
        # Compensate for the baseclass (which I'm too lazy to
        # subclass and correct) converting nm to eV in a typical file
        self.dataObj.proc_data[:, 0] = 1239.84/self.dataObj.proc_data[:,0]
        self.plotSpectrum()

        # params = self.genParameters(**dataObj.parameters)
        # self.specParams = Parameter.create(name=dataObj.fname, type='group',
        #                                    children=params)
        # self.ui.ptFile.setParameters(self.specParams, showTop=True)


    def plotSpectrum(self):
        self.ui.gSpectrum.setImage(self.dataObj.proc_data)

def updateFileClose(obj):
    try:
        fileList.remove(obj)
    except Exception as e:
        print "Error removing file from fileList:", e, obj
    if not fileList and not combinedWindowList:
        ex.exit(0)

def updateCompClose(obj):
    try:
        combinedWindowList.remove(obj)
    except Exception as e:
        print "error removign from list,", e
    if not fileList and not combinedWindowList:
        ex.exit(0)

def updateFocusList(obj):
    if obj not in combinedWindowList:
        return # thsi gets called when window is closed
    try:
        combinedWindowList.remove(obj)
        combinedWindowList.append(obj)
    except Exception as e:
        print "error updating focus list", e


# converter[<inp>][<out>]
#                                           I N P U T
converterArr = [               #nm                        #eV                          #wn
                 [lambda x: x,           lambda x:1239.84/x,            lambda x: 10000000./x          ], #nm
                 [lambda x: 1239.84/x,   lambda x: x,                   lambda x: 10000000. * x/1239.84], #eV
                 [lambda x: 10000000./x, lambda x: 1239.84/x/10000000., lambda x: x                   ]  # wn
]

converter = {
    "nm":         {"nm": lambda x: x,           "eV": lambda x:1239.84/x,            "wavenumber": lambda x: 10000000./x},
    "eV":         {"nm": lambda x: 1239.84/x,   "eV": lambda x: x,                   "wavenumber":lambda x: 10000000. * x/1239.84},
    "wavenumber": {"nm": lambda x: 10000000./x, "eV": lambda x: 1239.84/x/10000000., "wavenumber": lambda x: x}
}



if __name__=="__main__":
    import sys
    print "argv:", sys.argv
    print len(sys.argv)
    ex = QtGui.QApplication(sys.argv)
    win = BaseWindow()
    # win.openFile(sys.argv[1])
    print "made window"

    import pyqtgraph.console as pgc
    consoleWindow = pgc.ConsoleWidget(namespace={"fl":fileList,"np": np, "cl":combinedWindowList, "pg":pg})
    consoleWindow.show()
    consoleWindow.lower()

    sys.exit(ex.exec_())
    # ex.exec_()

import pyqtgraph as pg
import numpy as np
import sys
from PyQt5 import QtCore, QtGui
from .clickablePlotSettings_ui import Ui_LineSettingsDialog
from .PlotDataErrorItem import *
from ..packageSettings import config_options
from .fittingHelpers import *
from ..helperFuncs import *
from scipy.optimize import curve_fit as spcf



class signalBlocker(object):
    def __init__(self, toBlock):
        # print("Init. To block", toBlock)
        self.toBlock = toBlock
    def __call__(self, f):
        # print("called, f", f)
        def wrappedf(*args, **kwargs):
            args[0].ui.__getattribute__(self.toBlock).blockSignals(True)
            ret = f(*args, **kwargs)
            args[0].ui.__getattribute__(self.toBlock).blockSignals(False)
            return ret

        return wrappedf

class ClickablePlotWidget(pg.PlotWidget):
    # emitted when the window is closed
    sigPlotClosed = QtCore.pyqtSignal(object)
    sigMouseMoved = QtCore.pyqtSignal(object)

    # Emitted when an arrow key is pressed. Emits a
    # 2x1 np.array, first index is left/right (-1/+1) and
    # second is down/up (-1/+1, resp)
    sigArrowKeyPressed = QtCore.pyqtSignal(object)

    sigFitSettingsChange = QtCore.pyqtSignal(object)


    sigCrosshairsMoved = QtCore.pyqtSignal(object)

    def __init__(self, *args, **kwargs):
        super(ClickablePlotWidget, self).__init__(*args, **kwargs)

        self.opts = config_options.copy()
        self.updateOpts()

        self.selectedCurves = {}

        self.fitSettings = {"linearRegion": pg.LinearRegionItem(),
                            "button": None,
                            "curveToFit": None,
                            "fitFunc": None,
                            "p0": None,
                            "menu": None,
                            "fitCurve": pg.PlotDataItem(pen=pg.mkPen('k', width=3)),
                            "pText": pg.TextItem("Select a function", color=(0, 0, 0))}

        self.crosshairSettings = {"crosshairsVisible": False,
                                  "freeCrosshairs": False,
                                  "dataCrosshairs": False,
                                  "dataCrosshairsIdx": 0,
                                  "dataCrosshairsSource": 0}

        # Map these up here so I don't need to constantly figure out
        # how to add on to the contet menu of the view
        self.addContextMenu = self.plotItem.vb.menu.addMenu
        self.addContextAction = self.plotItem.vb.menu.addAction

        # set up a timer to tell when you've double
        # clicked a
        self.doubleClickTimer = QtCore.QTimer()
        self.doubleClickTimer.setInterval(250)
        self.doubleClickTimer.setSingleShot(True)

        self.sigFitSettingsChange.connect(self.updateFitSettings)

        self.scene().sigMouseMoved.connect(self.mousemove)
        self.scene().sigMouseClicked.connect(self.mouseClickedEvent)
        self.sigArrowKeyPressed.connect(self.updateFreeCrosshairs)
        self.sigArrowKeyPressed.connect(self.updateDataCrosshairs)

        self.doubleClickCurve = None

        p = pg.mkPen(color='r', width=3)
        self.crosshairs = [pg.InfiniteLine(pen=p), pg.InfiniteLine(angle=0, pen=p)]
        self.crosshairWindow = QtGui.QMessageBox(self)
        self.crosshairWindow.setWindowFlags(QtCore.Qt.WindowStaysOnTopHint)
        self.crosshairWindow.setWindowFlags(QtCore.Qt.Tool)
        self.crosshairWindow.setModal(False)
        self.plotItem.addItem(self.crosshairs[0])
        self.plotItem.addItem(self.crosshairs[1])
        self.sigCrosshairsMoved.connect(self.updateCrosshairWindow)
        # Hide these things
        self.removeFreeCrosshairs()

        self.plotItem.addItem(self.fitSettings["linearRegion"])
        self.fitSettings["button"] = QtGui.QToolButton()
        self.fitSettings["button"].setArrowType(QtCore.Qt.RightArrow)
        self.fitSettings["button"].clicked.connect(self.popupFitMenu)
        self.plotItem.vb.scene().addWidget(self.fitSettings["button"])

        self.fitSettings["linearRegion"].sigRegionChanged.connect(self.updateFitButtonPosition)
        self.fitSettings["linearRegion"].sigRegionChanged.connect(self.attemptFitting)
        self.plotItem.sigRangeChanged.connect(self.updateFitButtonPosition)
        self.plotItem.sigRangeChanged.connect(self.updateFitTextPosition)

        self.plotItem.addItem(self.fitSettings["pText"])
        self.fitSettings["pText"].hide()
        self.fitSettings["pText"].setFont(QtGui.QFont("", 24))
        self.fitSettings["pText"].setZValue(1001)


        # make it sit on top
        self.fitSettings["fitCurve"].setZValue(1000)
        # add it directly to the viewBox, NOT the plotItem, to avoid
        # adding it to the list of plotItem.curves, etc, since a lot
        # of my stuff depends heavily on plotItem.curves and I don't want
        # the fit curve to be seen with it.
        # Also, ignore the bounds so it can be ignored when calculating the autorange
        self.plotItem.vb.addItem(self.fitSettings["fitCurve"], ignoreBounds=True)

        self.plotItem.ctrl.logXCheck.toggled.connect(
            lambda: self.fitSettings["fitCurve"].setLogMode(
                self.plotItem.ctrl.logXCheck.isChecked(),
                self.plotItem.ctrl.logYCheck.isChecked()
            )
        )
        self.plotItem.ctrl.logYCheck.toggled.connect(
            lambda: self.fitSettings["fitCurve"].setLogMode(
                self.plotItem.ctrl.logXCheck.isChecked(),
                self.plotItem.ctrl.logYCheck.isChecked()
            )
        )
        self.fitSettings["fitCurve"].hide()

        self.removeFitRegion()

    def updateOpts(self, **opts):
        self.opts.update(opts)
        if self.opts["boundedAxes"]:
            self.plotItem.showAxis("top")

            self.plotItem.showAxis("right")
        for pos in ["top", "left", "right", "bottom"]:
            axis = self.plotItem.getAxis(pos)
            axis.setTickFont(QtGui.QFont("", self.opts["axisFontSize"]))
            axis.setLabel(**{"font-size":"{}pt".format(self.opts["axisLabelFontSize"])})
            axis.setPen(color=self.opts["foreground"],width=self.opts["axisThickness"])
            axis.setStyle(tickLength=self.opts["axisTickLength"])
            if axis.orientation in ["left", "right"]:
                # axis.setWidth(100)
                pass
                # axis.setStyle(tickTextOffset=20, tickTextWidth=20, autoExpandTextSpace=True)
            elif axis.orientation in ["bottom", "top"]:
                pass
            # axis.setMinimumWidth(300)
            # axis.setMinimumHeight(300)
            pg.ImageView

    def plot(self, *args, **kwargs):
        aargs, kkwargs = getPlotPens(pw=self, *args, **kwargs)
        kwargs.update(kkwargs)
        args = aargs
        p = self.plotItem.plot(*args, **kwargs)
        self.setItemClickable(p)
        return p

    def brazilPlot(self,*args, **kwargs):
        """
        Pass an nx3 as first arg, or kwargs data:
            Plot d[:,0]vs d[:,1]+d[:,2] :: d[:,1]-d[:,2]
            as d[:,1] is a data value, and d[:,2] is errors


        Pass two Nx2

        :param args:
        :param kwargs:
        :return:
        """
        data = kwargs.pop("data", None)
        data1 = kwargs.pop("data1", None)
        data2 = kwargs.pop("data2", None)
        args = list(args)

        if (None is data1 or None is data2) and None is data:
            # Assume it's in the *args... Figure out how
            # many numpy arrays there are.
            numnp = [ii for ii in args if type(ii) is np.ndarray]
            if len(numnp) == 0:
                raise RuntimeError("I don't know what to plot")
            if len(numnp) == 1:
                # Follow it down to below parsing
                if args[0].shape[1] in [3, ]:
                    data = args[0]
                    args.pop(0)
            if len(numnp) == 2:
                if args[0].shape[1] == 2 and args[1].shape[1] == 2:
                    data1 = args[0]
                    data2 = args[1]
                    args.pop(1)
                    args.pop(0)
                elif args[0].shape[1] == 1 and args[1].shape[1] == 2:
                    # args[0] is x,
                    # args[1] is [y1, y2]
                    data1 = np.column_stack((args[0], args[1][:, 0]))
                    data2 = np.column_stack((args[0], args[1][:, 1]))
                    args.pop(1)
                    args.pop(0)
            if len(numnp) == 3:
                if args[0].ndim == 1 and args[1].ndim == 1 and args[2].ndim == 1:
                    data = np.column_stack((args[0], args[1], args[2]))
                    args.pop(2)
                    args.pop(1)
                    args.pop(0)




        if data is not None:
            # either Nx3 or Nx4?
            # Is Nx3 [x, y1, y2], or [x, y, yerr]? need a flag?
            # I'm gonna assume the later cause that's how I'd use it...
            if data.shape[1] == 3:
                data1 = np.column_stack((data[:, 0], data[:, 1] + data[:, 2]))
                data2 = np.column_stack((data[:, 0], data[:, 1] - data[:, 2]))
            else:
                raise RuntimeError("Not sure how to parse this brazil input: {}".format(data.shape))


        aargs, kkwargs = getPlotPens(self, *args, **kwargs)
        kwargs.update(kkwargs)
        args = aargs

        c1 = pg.PlotDataItem(data1, **kwargs)
        c2 = pg.PlotDataItem(data2, **kwargs)
        c1.hide()
        c2.hide()
        # c1 = self.plot(data1, *args, **kwargs)
        # c2 = self.plot(data2, *args, **kwargs)
        kwargs.pop("symbolPen") #not used in fillbetween
        kwargs.pop("symbolBrush")  # not used in fillbetween
        kwargs.pop("symbol")  # not used in fillbetween
        brush = kwargs.pop("brush", None)
        if brush is None:
            color = kwargs.get("pen").color()

            alpha = kwargs.pop("alpha", None)
            if alpha is None:
                alpha = 0.5
            color.setAlphaF(alpha)

            brush = pg.mkBrush(color)
            kwargs["brush"] = brush
            # kwargs["pen"] = pg.mkPen(color)
        kwargs["pen"] = None
        fill = pg.FillBetweenItem(c1, c2, **kwargs)
        # p = self.plotItem.plot(*args, **kwargs)
        # self.setItemClickable(p)
        self.plotItem.addItem(fill)
        # self.setItemClickable(c1)
        # self.setItemClickable(c2)
        # self.plotItem.addItem(c1)
        # self.plotItem.addItem(c2)

        return fill


    def errorbars(self, *args, **kwargs):
        # """
        # create and add an errordataitem
        #
        # errorbars(<xdata>, <ydata>, <error>, *kwargs)
        # errorbars(<ydata>, <errors>, **kwargs) to assume
        #     x as an index
        # errorbars(x=<xdata>, y=<ydata>, errorbars=<errors>, **kwargs)
        # errorbars(Nx(2, 3), **kwargs) to unpack as
        #     x, y, yerr
        #
        # kwargs:
        # capWidth: widths of bar caps
        # """
        #
        # aargs, kkwargs = getPlotPens(self, *args, **kwargs)
        # kwargs.update(kkwargs)
        # args = aargs
        #
        #
        # if len(args)==2:
        #     kwargs.update(y=args[0])
        #     kwargs.update(errorbars=args[1])
        # elif len(args)==3:
        #     kwargs.update(x=args[0])
        #     kwargs.update(y=args[1])
        #     kwargs.update(errorbars=args[2])
        # elif len(args) == 1 and isinstance(args[0], np.ndarray):
        #     if args[0].shape[1] == 2:
        #         kwargs.update(y=args[0][:,0])
        #         kwargs.update(errorbars=args[0][:,1])
        #     elif args[0].shape[1] == 3:
        #         kwargs.update(x=args[0][:,0])
        #         kwargs.update(y=args[0][:,1])
        #         kwargs.update(errorbars=args[0][:,2])
        #     else:
        #         raise RuntimeError("I do not know how to parse this np.ndarray")
        # elif len(args)==1:
        #     raise RuntimeError("I do not know how to parse this argument", args)
        #
        # assert 'y' in kwargs
        # assert 'errorbars' in kwargs
        #
        # if 'capWidth' not in kwargs:
        #     kwargs['capWidth'] = config_options["errorBarCapWidth"]
        # kwargs['beam'] = kwargs.pop('capWidth')
        #
        # args, kwargs = getPlotPens(self, *args, **kwargs)
        erroritem = PlotDataErrorItem(pw = self, *args, **kwargs)

        self.addItem(erroritem)
        self.setItemClickable(erroritem)
        return erroritem

    def setItemClickable(self, item):
        """
        item: a data item which should be modified
        to add behavior for when it can be clicked
        and modified by the UI
        """

        item.curve.setClickable(True)
        item.sigClicked.connect(self._updateCurveClicked)
        item.opts['isEnabled'] = True
        item.opts['isSelected'] = False
        item.opts["plotWidget"] = self
        item.opts["xOffset"] = 0
        item.opts["yOffset"] = 0

    def updateLegendNames(self, curveItem):
        if curveItem.name() is None:
            return
        for sample, label in self.plotItem.legend.items:
            if sample.item is curveItem:
                label.setText(curveItem.opts["name"])
                # The legend will resize to fit new text,
                # but will not shrink back down after changing
                # to a shorter name.
                # For some reason. LegendItem.updateSize() doesn't
                # work properly.
                # For some reason, setting the geometry to 0 causes
                # it to resize properly.
                # Not sure what the miscommunication is from
                # PyQt/graph is.
                self.plotItem.legend.setGeometry(0, 0, 0, 0)
                break
        else:
            self.plotItem.legend.addItem(curveItem, name=curveItem.opts["name"])

    def _updateCurveClicked(self, plotDataItem):
        pen = plotDataItem.opts['pen']
        pen = pg.mkPen(pen)
        width = pen.width()
        # print("Curve clicked. Current width: {}".format(pen.width()))
        if self.doubleClickTimer.isActive():
            self.doubleClickTimer.stop()
            CurveItemSettings.getNewParameters(self, self.plotItem.curves, plotDataItem)
        elif plotDataItem.opts["isSelected"]:
            pen.setWidth(width - config_options["selectionThickness"])
            plotDataItem.setPen(pen)
            self.doubleClickTimer.start()
            plotDataItem.opts["isSelected"] = False
        else:
            pen.setWidth(width + config_options["selectionThickness"])
            plotDataItem.setPen(pen)
            plotDataItem.opts["isSelected"] = True
            self.doubleClickTimer.start()

    def closeEvent(self, *args, **kwargs):
        self.crosshairWindow.close()
        self.sigPlotClosed.emit(self)

    def mousemove(self, *args, **kwargs):
        self.sigMouseMoved.emit(self.plotItem.vb.mapSceneToView(args[0]))

    def keyPressEvent(self, ev):
        arr = [0, 0]
        if ev.key() == QtCore.Qt.Key_Left:
            arr[0] = -1
        elif ev.key() == QtCore.Qt.Key_Right:
            arr[0] = 1
        elif ev.key() == QtCore.Qt.Key_Up:
            arr[1] = 1
        elif ev.key() == QtCore.Qt.Key_Down:
            arr[1] = -1
        self.sigArrowKeyPressed.emit(np.array(arr))

        super(ClickablePlotWidget, self).keyPressEvent(ev)

    def updateFreeCrosshairs(self, newPos):
        if not self.crosshairSettings["freeCrosshairs"]:
            return
        drange = self.plotItem.vb.viewRange()
        xsp = np.abs(np.diff(drange[0])/100.)[0]
        ysp = np.abs(np.diff(drange[1])/100.)[0]

        change = np.array(newPos) * np.array([xsp, ysp])

        for i in [0, 1]:
            self.crosshairs[i].setPos(
                self.crosshairs[i].value() + change[i]
            )
        self.sigCrosshairsMoved.emit([i.value() for i in self.crosshairs])

    def updateDataCrosshairs(self, direction=np.array([0,0])):
        if not self.crosshairSettings["dataCrosshairs"]:
            return

        self.crosshairSettings["dataCrosshairsSource"] += direction[1]
        self.crosshairSettings["dataCrosshairsSource"] %= len(self.plotItem.curves)

        dataCurve = self.plotItem.curves[self.crosshairSettings["dataCrosshairsSource"]]
        xData, yData = dataCurve.getData()
        if xData is None or yData is None:
            print("it's nonE?")

        self.crosshairSettings["dataCrosshairsIdx"] += direction[0]
        self.crosshairSettings["dataCrosshairsIdx"] %= len(xData)

        self.crosshairs[0].setPos(
            xData[self.crosshairSettings["dataCrosshairsIdx"]]
        )
        self.crosshairs[1].setPos(
            yData[self.crosshairSettings["dataCrosshairsIdx"]]
        )
        self.sigCrosshairsMoved.emit([i.value() for i in self.crosshairs])

    def updateCrosshairWindow(self):
        x = self.crosshairs[0].value()
        y = self.crosshairs[1].value()
        if self.crosshairSettings["dataCrosshairs"]:
            curve = self.plotItem.curves[
                self.crosshairSettings["dataCrosshairsSource"]]
            idx =self.crosshairSettings["dataCrosshairsIdx"]
            xD, yD = curve.getData()
            name = curve.name()
            if name is None: name="Data"
            txt = "{}:\nx[{}]={}, y[{}]={}".format(name, idx, xD[idx], idx, yD[idx])
        else:
            txt = "x={}, y={}".format(x, y)
        self.crosshairWindow.setText(txt)

    def addFreeCrosshairs(self):
        self.crosshairSettings["crosshairsVisible"] = True
        self.crosshairSettings["freeCrosshairs"] = True
        drange = self.plotItem.vb.viewRange()
        self.crosshairs[0].setValue(np.mean(drange[0]))
        self.crosshairs[1].setValue(np.mean(drange[1]))

        # self.addItem(self.crosshairs[0])
        # self.addItem(self.crosshairs[1])
        self.crosshairs[0].show()
        self.crosshairs[1].show()

        self.crosshairWindow.show()
        self.updateCrosshairWindow()

    def removeFreeCrosshairs(self):
        self.crosshairSettings["crosshairsVisible"] = False
        self.crosshairSettings["freeCrosshairs"] = False
        self.crosshairWindow.hide()
        self.crosshairs[0].hide()
        self.crosshairs[1].hide()

    def addDataCrosshairs(self):
        self.crosshairSettings["crosshairsVisible"] = True
        self.crosshairSettings["dataCrosshairs"] = True

        self.updateDataCrosshairs()

        self.crosshairs[0].show()
        self.crosshairs[1].show()

        self.crosshairWindow.show()
        self.updateCrosshairWindow()

    def removeDataCrosshairs(self):
        self.crosshairSettings["crosshairsVisible"] = False
        self.crosshairSettings["dataCrosshairs"] = False
        self.crosshairWindow.hide()
        self.crosshairs[0].hide()
        self.crosshairs[1].hide()

    def addFitRegion(self):
        xrange = self.plotItem.vb.viewRange()[0]
        d = xrange[1]-xrange[0]
        self.fitSettings["linearRegion"].setRegion((xrange[0] + d*0.2, xrange[1] - d*0.2))
        self.updateFitButtonPosition()

        if self.fitSettings["curveToFit"] is None:
            self.fitSettings["curveToFit"] = self.plotItem.curves[0]

        self.fitSettings["linearRegion"].show()
        self.fitSettings["button"].show()
        self.fitSettings["fitCurve"].show()
        self.fitSettings["pText"].show()
        self.updateFitTextPosition()

    def updateFitButtonPosition(self):
        if not self.fitSettings["button"].isVisible(): return
        yrange = self.plotItem.vb.viewRange()[1]
        x = self.fitSettings["linearRegion"].getRegion()[1]
        p = self.plotItem.vb.mapViewToDevice(QtCore.QPointF(x, yrange[1]))

        self.fitSettings["button"].setGeometry(
            p.x(), p.y(), self.fitSettings["button"].geometry().width(), self.fitSettings["button"].geometry().height()
        )

    def updateFitTextPosition(self):
        if not self.fitSettings["pText"].isVisible(): return
        vrange = self.plotItem.vb.viewRange()
        dx = vrange[0][1] - vrange[0][0]
        dy = vrange[1][1] - vrange[1][0]

        x = np.sum(vrange[0])/2.
        self.fitSettings["pText"].setPos(vrange[0][0] + dx*0.05, vrange[1][1] - dy*0.05)

    def removeFitRegion(self):
        # self.plotItem.removeItem(self.fitSettings["linearRegion"])
        self.fitSettings["linearRegion"].hide()
        self.fitSettings["button"].hide()
        self.fitSettings["fitCurve"].hide()
        self.fitSettings["pText"].hide()

    def popupFitMenu(self):
        menu = self.getFitMenu()
        # Get the cursor position to poppup
        # the menu
        p = QtGui.QCursor.pos()
        menu.popup(p)

    def getFitMenu(self):
        menu = QtGui.QMenu()
        mData = menu.addMenu("Select Data...")
        mFit = menu.addMenu("Set Function...")

        for c in self.plotItem.curves:
            name = str(c.name())
            if name is None:
                name = str(c)
            a = mData.addAction(name)
            a.setCheckable(True)
            if c is self.fitSettings["curveToFit"]:
                a.setChecked(True)
            # Need the c=c to overcome a scoping thing in python
            # when you loop through lambdas like this
            a.triggered.connect(lambda x, c=c: self.sigFitSettingsChange.emit(c))



        funcs = list(getFuncs.keys())
        funcs.append("Other...")
        for n in funcs:
            a = mFit.addAction(n)
            a.setCheckable(True)
            if n == self.fitSettings["fitFunc"]:
                a.setChecked(True)
            a.triggered.connect(lambda x, n=n: self.sigFitSettingsChange.emit(n))

        self.fitSettings["menu"] = menu

        return menu

    def updateFitSettings(self, *args, **kwargs):
        if isinstance(args[0], pg.PlotDataItem):
            # Been told to change the data source
            self.fitSettings["curveToFit"] = args[0]
            self.fitSettings["p0"] = None
        elif isinstance(args[0], str):
            # been told to change the fit function
            if args[0] == "Other...": return
            self.fitSettings["fitFunc"] = args[0]
            self.fitSettings["p0"] = None
        self.attemptFitting()

    def attemptFitting(self):
        # get the full data
        try:
            x = self.fitSettings["curveToFit"].xData
            y = self.fitSettings["curveToFit"].yData
        except AttributeError:
            # no dataset has been selected, abort
            return
        if self.fitSettings["fitFunc"] is None: return

        # now need to slice it by the region of itnerest
        rg = self.fitSettings["linearRegion"].getRegion()
        xidx0 = np.where(x>rg[0])[0]
        xidx1 = np.where(x<rg[1])[0]
        xidx = list(set(xidx0) & set(xidx1))

        x = x[xidx]
        y = y[xidx]

        d = np.column_stack((x, y))
        d = d[d[:,0].argsort()]
        x = d[:,0]
        y = d[:,1]
        f = getFuncs[self.fitSettings["fitFunc"]]

        if self.fitSettings["p0"] is None:
            self.fitSettings["p0"] = f.guessParameters(x, y)
            # print("return from guess:", self.fitSettings["p0"])

        try:
            p, _ = spcf(f.__call__, x, y, p0 = self.fitSettings["p0"])
        except Exception as e:
            print("Error fitting function,",e)
            return
        self.fitSettings["p0"] = p

        # I really wanted to put this into a UI setting, but I can't
        # remember how things are laid out for this
        # This one will plot the full view range
        xfit = self.plotItem.vb.viewRange()[0]
        xfit = np.linspace(xfit[0], xfit[1])
        # This one will plot only within the linearregion
        ### xfit = np.linspace(x[0], x[-1], 250)


        y = f(xfit, *self.fitSettings["p0"])
        self.fitSettings["fitCurve"].setData(xfit, y)
        st = str(f).format(*p)

        self.fitSettings["pText"].setText(st, color=(0,0,0))

    def addFunctionLine(self):
        pass

    def addLegend(self, *args, **kwargs):
        if self.plotItem.legend is None:
            self.plotItem.addLegend()
        return self.plotItem.legend

    def mouseClickedEvent(self, ev):
        # Intercept if the crosshairs
        # are visible so we can move their position
        if self.crosshairSettings["crosshairsVisible"] and\
                        ev.button() == QtCore.Qt.LeftButton:

            ev.accept()
            p = self.plotItem.vb.mapSceneToView(ev.scenePos())

            # if free, set the crosshairs to the clicked position
            if self.crosshairSettings["freeCrosshairs"]:
                self.crosshairs[0].setPos(p.x())
                self.crosshairs[1].setPos(p.y())
            # otherwise, find the nearest datapoint and go to it.
            elif self.crosshairSettings["dataCrosshairs"]:
                xD, yD = self.plotItem.curves[
                    self.crosshairSettings["dataCrosshairsSource"]].getData()
                xD = np.nan_to_num(xD)
                yD = np.nan_to_num(yD)
                dist = (p.x()-xD)**2 + (p.y()-yD)**2
                self.crosshairSettings["dataCrosshairsIdx"] = np.argmin(dist)
                self.updateDataCrosshairs()
            self.updateCrosshairWindow()


class CurveItemSettings(QtGui.QDialog):
    # QListWidgets allows multiple selection of listItems
    # However, if you toggle the checkbox of a set of selected items,
    # it will de-select all but the checked one, which is frustrating
    # if you want to toggle it (i.e. to quickly see differences).
    # I put this in place to prevent this. If you check multiple
    # items, start a quick single-shot timer, which another
    # function will check to see if it's running to reselect those
    # items.
    _multiSelectHelper = {"timer": QtCore.QTimer(),
                          "selectionList": []}
    def __init__(self, *args, **kwargs):
        curveList = kwargs.pop('curves', None)
        clickedCurve = kwargs.pop('clickedCurve', None)
        super(CurveItemSettings, self).__init__(*args, **kwargs)
        if curveList is None:
            curveList = [pg.PlotDataItem()]
        self.listViewCurveItems = {}
        self.initUI(curveList, clickedCurve)

        self._multiSelectHelper["timer"].setInterval(100)
        self._multiSelectHelper["timer"].setSingleShot(True)


    def initUI(self, curveList, firstCurve):
        """

        :param curveItem: curveitim
        :type curveItem: pg.PlotDataItem
        :return:
        """
        self.ui = Ui_LineSettingsDialog()
        self.ui.setupUi(self)
        self.ui.sbLineWidth.setMinimum(0)
        self.ui.sbLineWidth.setSingleStep(1)
        self.ui.sbLineWidth.setOpts(int=True)

        self.ui.sbMarkerSize.setMinimum(0)
        self.ui.sbMarkerSize.setSingleStep(1)
        self.ui.sbMarkerSize.setOpts(int=True)
        self.originalSettings = {}

        initialListItem = None

        for curveItem in curveList:
            listItem = QtGui.QListWidgetItem(self.ui.lwCurves)
            listItem.setFlags(QtCore.Qt.ItemIsUserCheckable | QtCore.Qt.ItemIsSelectable |
                              QtCore.Qt.ItemIsEnabled | QtCore.Qt.ItemIsEditable)
            listItem.setCheckState(QtCore.Qt.Unchecked)
            listItem.setSelected(curveItem.opts["isSelected"])
            if curveItem is firstCurve:
                initialListItem = listItem

            name = curveItem.name()
            if name is None:
                name = curveItem
            name = str(name)
            listItem.setText(name)

            cs = QtCore.Qt.Checked if curveItem.opts["isEnabled"] else QtCore.Qt.Unchecked
            listItem.setCheckState(cs)


            self.ui.lwCurves.addItem(listItem)
            self.listViewCurveItems[curveItem] = listItem
            self.originalSettings[curveItem] = curveItem.opts.copy()


        self.ui.lwCurves.currentItemChanged.connect(self.handlecurrentItemChanged)
        self.ui.lwCurves.itemChanged.connect(self.handleitemChanged)
        self.ui.lwCurves.itemSelectionChanged.connect(self.handleitemSelectionChanged)

        self.ui.colLine.sigColorChanging.connect(self.updateLineColors)
        self.ui.cbLineStyle.currentIndexChanged.connect(self.updateLineStyle)
        self.ui.sbLineWidth.sigValueChanging.connect(self.updateLineWidth)

        self.ui.colMarker.sigColorChanging.connect(self.updateMarkerBrushColor)
        self.ui.cbMarkerStyle.currentIndexChanged.connect(self.updateMarkerStyle)
        self.ui.sbMarkerSize.sigValueChanging.connect(self.updateMarkerSize)


        self.show()
        # toggle it to set the interface to have the settings
        # of the clicked line
        self.handlecurrentItemChanged(initialListItem)

    def getCurveFromItem(self, item):
        for k, v in list(self.listViewCurveItems.items()):
            if v is item: return k

    def handlecurrentItemChanged(self, *args, **kwargs):
        item = args[0]
        curve = self.getCurveFromItem(item)
        pen = pg.mkPen(curve.opts["pen"])
        opts = curve.opts

        self.ui.colLine.blockSignals(True)
        self.ui.colLine.setColor(pen.color())
        self.ui.colLine.blockSignals(False)

        self.ui.cbLineStyle.blockSignals(True)
        self.ui.cbLineStyle.setCurrentIndex(pen.style())
        self.ui.cbLineStyle.blockSignals(False)

        self.ui.sbLineWidth.blockSignals(True)
        self.ui.sbLineWidth.setValue(pen.width() -
                                     config_options["selectionThickness"]*int(curve.opts["isSelected"]))
        self.ui.sbLineWidth.blockSignals(False)

        self.ui.cbMarkerStyle.blockSignals(True)
        self.ui.cbMarkerStyle.setCurrentIndex(
            self.ui.cbMarkerStyle.findText(str(opts["symbol"]))
        )
        self.ui.cbMarkerStyle.blockSignals(False)

        self.ui.sbMarkerSize.blockSignals(True)
        self.ui.sbMarkerSize.setValue(opts["symbolSize"])
        self.ui.sbMarkerSize.blockSignals(False)


        self.ui.colMarker.blockSignals(True)
        self.ui.colMarker.setColor(pg.mkBrush(opts["symbolBrush"]).color())
        self.ui.colMarker.blockSignals(False)

    @signalBlocker("lwCurves")
    def handleitemChanged(self, *args, **kwargs):
        clickedItem = args[0]
        if clickedItem not in self.ui.lwCurves.selectedItems():
            self.ui.lwCurves.setCurrentItem(clickedItem)
        for listItem in self.ui.lwCurves.selectedItems():
            curve = self.getCurveFromItem(listItem)
            if listItem is clickedItem:
                pass
            else:
                # Want to toggle other things if they're not in this list
                listItem.setCheckState(2*(not listItem.checkState()))
                listItem.setSelected(True)
            # if curve.opts["isEnabled"] ^ (listItem.checkState() == QtCore.Qt.Checked):
            if True:
                curve.opts["isEnabled"] = (listItem.checkState() == QtCore.Qt.Checked)
                if not curve.opts["isEnabled"]:
                    curve.curve.hide()
                    curve.scatter.hide()
                    if hasattr(curve, 'errorbars'):
                        curve.errorbars.hide()
                else:
                    curve.curve.show()
                    curve.scatter.show()
                    if hasattr(curve, 'errorbars'):
                        curve.errorbars.show()
        item = args[0]
        curve = self.getCurveFromItem(item)

        if curve.name() != str(item.text()):
            curve.opts["name"] = str(item.text())
            if curve.opts["plotWidget"].plotItem.legend is not None:
                curve.opts["plotWidget"].updateLegendNames(curve)

        self._multiSelectHelper["selectionList"] = self.ui.lwCurves.selectedItems()
        self._multiSelectHelper["timer"].start()

    @signalBlocker("lwCurves")
    def handleitemSelectionChanged(self, *args, **kwargs):
        if self._multiSelectHelper["timer"].isActive():
            for item in self._multiSelectHelper["selectionList"]:
                item.setSelected(True)
        # print("Timer active?", self._multiSelectHelper["timer"].isActive())
        # print("Selectionlist:", self._multiSelectHelper["selectionList"])

    def updateLineColors(self, colorButton):
        for listItem in self.ui.lwCurves.selectedItems():
            curve = self.getCurveFromItem(listItem)
            p = pg.mkPen(curve.opts["pen"])
            p.setColor(colorButton.color())
            curve.setPen(p)
        self.ui.colMarker.setColor(colorButton.color())
        self.updateMarkerBrushColor(colorButton)

    def updateLineStyle(self, newIdx):
        for listItem in self.ui.lwCurves.selectedItems():
            curve  = self.getCurveFromItem(listItem)
            p = pg.mkPen(curve.opts["pen"])
            p.setStyle(newIdx)
            curve.setPen(p)

    def updateLineWidth(self, sb):
        for listItem in self.ui.lwCurves.selectedItems():
            curve  = self.getCurveFromItem(listItem)
            p = pg.mkPen(curve.opts["pen"])
            p.setWidth(sb.value() +
                       config_options["selectionThickness"]*int(curve.opts["isSelected"]))
            curve.setPen(p)

    def updateMarkerStyle(self, newIdx):
        for listItem in self.ui.lwCurves.selectedItems():
            curve  = self.getCurveFromItem(listItem)
            if str(self.ui.cbMarkerStyle.currentText()) == "None":
                curve.setSymbol(None)
            else:
                curve.setSymbol(
                    str(self.ui.cbMarkerStyle.currentText())
                )

    def updateMarkerSize(self, sb):
        for listItem in self.ui.lwCurves.selectedItems():
            curve  = self.getCurveFromItem(listItem)
            curve.setSymbolSize(sb.value())

    def updateMarkerBrushColor(self, colorButton):
        for listItem in self.ui.lwCurves.selectedItems():
            curve  = self.getCurveFromItem(listItem)
            brush = pg.mkBrush(curve.opts["symbolBrush"])
            pen = pg.mkPen(curve.opts["symbolPen"])
            brush.setColor(colorButton.color())
            pen.setColor(colorButton.color())
            curve.setSymbolBrush(brush)
            curve.setSymbolPen(pen)

    @staticmethod
    def getNewParameters(parent, curves, clickedCurve = None):
        dialog = CurveItemSettings(curves=curves, clickedCurve = clickedCurve)
        parent.sigPlotClosed.connect(dialog.close)
        ok = dialog.exec_()
        if not ok:
            for curve in curves:
                # curve.opts = dialog.originalSettings[curve].copy()
                curve.opts.update(dialog.originalSettings[curve])
                curve.updateItems()
                if not curve.opts["isEnabled"]:
                    curve.curve.hide()
                    curve.scatter.hide()
                    if hasattr(curve, 'errorbars'):
                        curve.errorbars.hide()
                if curve.opts["plotWidget"].plotItem.legend is not None:
                    curve.opts["plotWidget"].updateLegendNames(curve)
        else:
            for curve in curves:
                if not curve.opts["isEnabled"]:
                    curve.opts["isSelected"] = False


if __name__=='__main__':
    ex = QtGui.QApplication([])
    wid = ClickablePlotWidget()
    wid.plotItem.addLegend()
    wid.plot([0, 1, 2, 3, 4], pen='g', width=5, name='green')
    p = wid.plot([1, -2, 4, 8, 2], pen='r', name='r')
    wid.plot([1]*5, pen='y')
    x = np.arange(5)
    y = np.random.normal(1, 1, 5)
    err = np.ones(5)

    err = wid.errorbars(np.column_stack((x, y, err)), pen='c', name='errors')

    wid.show()
    # sys.exit()
    sys.exit(ex.exec_())

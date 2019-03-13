import numpy as np
import pyqtgraph as pg
from .images.ImageViewWithPlotItemContainer import ImageViewWithPlotItemContainer
from .curves.clickablePlotWidget import ClickablePlotWidget
from .widgets.LabviewSlider import LabviewSlider as LS

from PyQt5 import QtCore, QtGui, QtWidgets


class BaseIcon(QtGui.QIcon):
    def __init__(self, *args, **kwargs):
        img = QtGui.QImage(200, 200, QtGui.QImage.Format_ARGB32)
        img.fill(QtCore.Qt.transparent)
        painter = QtGui.QPainter()
        painter.begin(img)
        self.drawIconPainter(painter, img.height(), img.width())
        painter.end()
        pxmap = QtGui.QPixmap.fromImage(img)
        super(BaseIcon, self).__init__(pxmap)

    def drawIconPainter(self, p, h, w):
        """
        :type p: QtGui.QPainter
        """
        raise NotImplementedError

class CrosshairIcon(BaseIcon):
    def drawIconPainter(self, p, h, w):
        """
        :type p: QtGui.QPainter
        """

        p.setPen(pg.mkPen('r', width=8))


        path = QtGui.QPainterPath()
        path.moveTo(w/2, h*.2)
        path.lineTo(w/2, h*.8)
        path.moveTo(w*.2, h/2)
        path.lineTo(w*.8, h/2)

        p.drawPath(path)

        return p

class DataBrowseIcon(BaseIcon):
    def drawIconPainter(self, p, h, w):
        """
        :type p: QtGui.QPainter
        """

        p.setPen(pg.mkPen('r', width=8))


        path = QtGui.QPainterPath()
        path.moveTo(w/2, h*.2)
        path.lineTo(w/2, h*.8)
        path.moveTo(w*.2, h/2)
        path.lineTo(w*.8, h/2)
        path.addRect(w*.35, h*.35, w*.3, h*.3)

        p.drawPath(path)

        return p

class FitFunctionIcon(BaseIcon):
    def drawIconPainter(self, p, height, width):
        """
        :type p: QtGui.QPainter
        """
        pen = QtGui.QPen(QtCore.Qt.red)
        pen.setWidth(15)
        p.setPen(pen)

        x = np.linspace(0, 2.5*np.pi, num=100)
        yorig = -np.sin(x) + x/(1*2.5*np.pi)
        yerr = np.random.normal(0, 0.25, size=len(x))

        y = yorig + yerr
        x = x/np.max(x) * width
        yorig -= np.min(y)
        y -= np.min(y)
        yorig = yorig/np.max(y) * height
        y = y/np.max(y) * height

        yorig = height - yorig
        y = height-y

        # add some padding
        x = x*0.6 + 0.2*width
        y = y*0.6 + 0.2*height

        # path = QtGui.QPainterPath()
        pnts = [ QtCore.QPointF(xi, yi) for (xi, yi) in zip(x[::5], y[::5])]
        lns = [QtCore.QPointF(xi, yi) for (xi, yi) in zip(x, yorig)]
        p.drawPoints(*pnts)

        pen.setColor(QtCore.Qt.black)
        pen.setWidth(8)
        p.setPen(pen)
        p.drawPolyline(*lns)

class AddFunctionIcon(BaseIcon):
    def drawIconPainter(self, p, h, w):
        """
        :type p: QtGui.QPainter
        """
        font = p.font()
        font.setPointSize(60)
        p.setFont(font)
        rect = QtCore.QRectF(w*0.2, h*0.2, w*0.6, h*0.6)
        p.drawText(rect, QtCore.Qt.AlignCenter, "f(x)")


class PlotContainerWindow(QtGui.QMainWindow):
    def __init__(self, *args, **kwargs):
        super(PlotContainerWindow, self).__init__()
        self.plotWidget = kwargs.pop("plotWidget", None)
        if self.plotWidget is None:
            self.plotWidget = ClickablePlotWidget()
        self.setCentralWidget(self.plotWidget)
        self.plotWidget.sigMouseMoved.connect(self.updateMousePos)

        self.xPosStatus = QtGui.QLabel(self)
        self.xPosStatus.setText("x=")
        self.yPosStatus = QtGui.QLabel(self)
        self.yPosStatus.setText("y=")
        self.statusBar().addPermanentWidget(self.xPosStatus)
        self.statusBar().addPermanentWidget(self.yPosStatus)

        self.plotTools = QtGui.QToolBar()
        self.plotTools.setMovable(False)

        a = self.plotTools.addAction(DataBrowseIcon(), 'Data Crosshairs')
        a.setCheckable(True)
        a.toggled.connect(self.toggleDataCrosshairs)

        a = self.plotTools.addAction(CrosshairIcon(), 'Screen Crosshairs')
        a.setCheckable(True)
        a.toggled.connect(self.togglePlotCrosshairs)

        a = self.plotTools.addAction(FitFunctionIcon(), 'Fit...')
        a.setCheckable(True)
        a.toggled.connect(self.toggleFitRegion)

        a = self.plotTools.addAction(AddFunctionIcon(), 'Add function...')
        a.setCheckable(True)
        a.toggled.connect(self.toggleFunctionLine)




        self.addToolBar(QtCore.Qt.LeftToolBarArea, self.plotTools)

    def __getattr__(self, item):
        try:
            return getattr(self.plotWidget, item)
        except Exception as e:
            pass
            # print(("Does it not work like this?", item, e, self.plotWidget))
            # print((hasattr(self.plotWidget, item)))

    def _repr_png_(self):
        QtGui.QApplication.processEvents()
        size = self.viewRect().size().toSize()
        size = self.geometry().size()
        # size.setHeight(300)
        print("size", size)
        try:
            self.image = QtGui.QImage(size,
                                QtGui.QImage.Format_RGB32)
        except Exception as e:
            print("Exception is ", e)
        # except AttributeError:
        #     self._widget.updateGL()
        #     self.image = self._widget.grabFrameBuffer()

        painter = QtGui.QPainter(self.image)
        self.render(painter)

        byte_array = QtCore.QByteArray()
        buffer = QtCore.QBuffer(byte_array)
        buffer.open(QtCore.QIODevice.ReadWrite)
        self.image.save(buffer, 'PNG')
        buffer.close()

        return bytes(byte_array)


    def updateMousePos(self, pos):
        self.xPosStatus.setText("x={}".format(pos.x()))
        self.yPosStatus.setText("y={}".format(pos.y()))

    def toggleDataCrosshairs(self, enabled):
        if enabled:
            for act in self.plotTools.actions():
                if act is self.sender(): continue
                act.setChecked(False)
            self.plotWidget.addDataCrosshairs()
        if not enabled:
            self.plotWidget.removeDataCrosshairs()

    def togglePlotCrosshairs(self, enabled):
        if enabled:
            for act in self.plotTools.actions():
                if act is self.sender(): continue
                act.setChecked(False)
            self.plotWidget.addFreeCrosshairs()
        if not enabled:
            self.plotWidget.removeFreeCrosshairs()

    def toggleFitRegion(self, enabled):
        if enabled:
            for act in self.plotTools.actions():
                if act is self.sender(): continue
                act.setChecked(False)
            self.plotWidget.addFitRegion()
        else:
            self.plotWidget.removeFitRegion()

    def toggleFunctionLine(self, enabled):
        if enabled:
            for act in self.plotTools.actions():
                # if act is self.sender(): continue
                act.setChecked(False)
            self.plotWidget.addFunctionLine()

    def maximize(self):
        self.showMaximized()

class ManipulateWindow(QtGui.QMainWindow):
    def __init__(self, *args, **kwargs):
        super(ManipulateWindow, self).__init__()
        centralWid = QtWidgets.QWidget()

        layout = QtWidgets.QVBoxLayout(self)

        self._manipulatorLayout = QtWidgets.QGridLayout()
        layout.addLayout(self._manipulatorLayout)

        self._plotWindow = PlotContainerWindow()

        # keep a list f the current items so we have a reference to which curves should
        # be updated
        self._updateCurves = self._plotWindow.plotWidget.plotItem.curves.copy()


        layout.addWidget(self._plotWindow)

        centralWid.setLayout(layout)
        self.setCentralWidget(centralWid)

        self._callBack = None



    def setManipulators(self, manipulations):
        """
        Pass manipulators as

        [
            ("name1", lowerBound, upperBound, <startVal>, <step>),
            ("name2", lowerBound, upperBound, <startVal>, <step>),            ...
        ]
        can optionally pass the first argument as a callback function
        :param manipulations:
        :return:
        """

        ## TODO: clear the layout? Or assume this only gets
        ## set once per instance?
        if manipulations is None: return



        if callable(manipulations[0]):
            self.setCallable(manipulations.pop(0)) # send it off to be the callbac,


        for idx, (lbl, *bnds) in enumerate(manipulations):
            self._manipulatorLayout.addWidget(QtWidgets.QLabel(lbl), idx, 0)

            slider = LS(range = bnds[:2])
            if len(bnds)>2:
                slider.setValue(bnds[2])
                if len(bnds)>3:
                    slider.setStep(bnds[3])
            # slider.sigValueChanging.connect(self.recalculate)

            slider.sigValueChanging.connect(self.updateCurve)
            self._manipulatorLayout.addWidget(slider, idx, 1)

    def setCallable(self, callback):
        """
        callback function should return a list of values,
        lst[0] is the x values, lst[1:] are the y values, one for each curve to be
        updated

        ManipulateWindow.plot() must be called witih as many curve which get manipulated
        before a call to setcallable. This allows plotting other functions on top of
        the manipulate curve
        :param callback:
        :return:
        """

        self._updateCurves = self._plotWindow.plotWidget.plotItem.curves.copy()
        self._callBack = callback

    def recalculate(self):
        if self._callBack is None: return
        # rc = self._manipulatorLayout.rowCount()
        # callVals = []
        # for idx in range(rc):
        #     callVals.append(
        #         self._manipulatorLayout.itemAtPosition(idx, 1).widget().value()
        #     )

        callVals = [ii[1] for ii in self.getCurrentValues()]
        ret = self._callBack(*callVals)
        return ret

    def getCurrentValues(self):
        out = []
        rc = self._manipulatorLayout.rowCount()

        for idx in range(rc):
            out.append([
            self._manipulatorLayout.itemAtPosition(idx, 0).widget().text(),
            self._manipulatorLayout.itemAtPosition(idx, 1).widget().value(),
                ])

        return out


        # if len(self._updateCurves) != len(ret)-1:
        #     raise RuntimeError("Expected updates for {} curves, got {}".format(
        #         len(self._updateCurves),len(ret)-1))
        #
        # for idx, curve in enumerate(self._updateCurves):
        #     curve.setData(ret[0], ret[idx+1])

    def updateCurve(self):
        if self._callBack is None: return

        ret = self.recalculate()

        # if len(self._updateCurves) != len(ret)-1:
        #     raise RuntimeError("Expected updates for {} curves, got {}".format(
        #         len(self._updateCurves),len(ret)-1))
        #
        # for idx, curve in enumerate(self._updateCurves):
        #     curve.setData(ret[0], ret[idx+1])
        if len(self._updateCurves) == len(ret) - 1:
            # passed as a list where the first element is x, the rest is y
            for idx, curve in enumerate(self._updateCurves):
                curve.setData(ret[0], ret[idx+1])
        elif len(self._updateCurves) == len(ret) and isinstance(ret[0], np.ndarray):
            # passed where each is a set of x,y datasets
            for idx, curve in enumerate(self._updateCurves):
                curve.setData(ret[idx][:,0], ret[idx][:,1])
        else:
            raise RuntimeError("Expected updates for {} curves, got {}".format(
                    len(self._updateCurves), len(ret) - 1))







    def __getattr__(self, item):
        try:
            return getattr(self._plotWindow, item)
        except Exception as e:
            print(("Does it not work like this?", item, e, self._plotWindow))
            print((hasattr(self._plotWindow, item)))

































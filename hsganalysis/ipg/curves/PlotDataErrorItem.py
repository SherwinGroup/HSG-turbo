
import pyqtgraph as pg
import numpy as np
import sys
from PyQt5 import QtCore, QtGui
from .clickablePlotSettings_ui import Ui_LineSettingsDialog
from ..packageSettings import config_options
from ..helperFuncs import *


class PlotDataErrorItem(pg.PlotDataItem):
    """
    This class will wrap a pg.PlotDataItem
    (which, in turn, wraps a curveItem and scatterPlotItem)
    and an errorbarItem, so they can all be wrapped together.
    Really hate that this isn't a standard feature.

    Need to make this a subclass of PlotDataItem, because
    some stuff was hardcoded in LegendItems which would
    be a huge headache to work around if you don't subclass.
    """
    def __init__(self, *args, **kwargs):
        """
        create and add an errordataitem

        errorbars(<xdata>, <ydata>, <error>, *kwargs)
        errorbars(<ydata>, <errors>, **kwargs) to assume
            x as an index
        errorbars(x=<xdata>, y=<ydata>, errorbars=<errors>, **kwargs)
        errorbars(Nx(2, 3), **kwargs) to unpack as
            x, y, yerr

        kwargs:
        capWidth: widths of bar caps
        """

        aargs, kkwargs = getPlotPens(*args, **kwargs)
        kwargs.update(kkwargs)
        args = aargs

        if len(args) == 2:
            kwargs.update(y=args[0])
            kwargs.update(errorbars=args[1])
        elif len(args) == 3:
            kwargs.update(x=args[0])
            kwargs.update(y=args[1])
            kwargs.update(errorbars=args[2])
        elif len(args) == 1 and isinstance(args[0], np.ndarray):
            if args[0].shape[1] == 2:
                kwargs.update(y=args[0][:, 0])
                kwargs.update(errorbars=args[0][:, 1])
            elif args[0].shape[1] == 3:
                kwargs.update(x=args[0][:, 0])
                kwargs.update(y=args[0][:, 1])
                kwargs.update(errorbars=args[0][:, 2])
                args = list(args)
                args[0]=args[0][:, :2]
            else:
                raise RuntimeError("I do not know how to parse this np.ndarray")
        elif len(args) == 1:
            raise RuntimeError("I do not know how to parse this argument", args)

        assert 'y' in kwargs
        assert 'errorbars' in kwargs

        if 'capWidth' not in kwargs:
            kwargs['capWidth'] = config_options["errorBarCapWidth"]
        kwargs['beam'] = kwargs.pop('capWidth')

        if "label" in kwargs:
            kwargs["name"] = kwargs.pop("label")


        ## 3/21/18
        ## I changed this argument call. It was throwing errors on
        ## passing an Nx3 np array.
        # args, kwargs = getPlotPens(self, *args, **kwargs)
        args, kwargs = getPlotPens(*args, **kwargs)
        self.errorbars = pg.ErrorBarItem(**kwargs)
        super(PlotDataErrorItem, self).__init__(*args, **kwargs)
        self.errorbars.setParentItem(self)

    def setData(self, *args, **kwargs):
        # do I need to keep a reference to this?
        self.errorData = np.array(kwargs.pop("errorbars", None))
        super(PlotDataErrorItem, self).setData(*args, **kwargs)
        kwargs.update({'x': self.xData})
        kwargs.update({'y': self.yData})
        kwargs.update({'height': 2*self.errorData})

        self.errorbars.setData(**kwargs)

        # I don't want error bars to follow a dashed style.
        # This isn't a good way to do it, but I'm busy right now.
        pen = kwargs.pop("pen", None)
        if pen is not None:
            self.setPen(pen)

    def setPen(self, *args, **kwargs):
        super(PlotDataErrorItem, self).setPen(*args, **kwargs)

        # force solid lines. Is this always wanted?
        pen = pg.mkPen(self.opts["pen"])
        pen.setStyle(1)
        self.errorbars.setOpts(pen=pen)
        # self.errorbars.setOpts(pen=self.opts['pen'])

    def setLogMode(self, xMode, yMode):
        super(PlotDataErrorItem, self).setLogMode(xMode, yMode)

        # errrobaritem doesn't have an implementation
        # for doing log things. Try to do some hackish
        # stuff here so it at least won't break.
        kwargs = {}
        kwargs.update({'x': self.xDisp})
        kwargs.update({'y': self.yDisp})

        if yMode:
            # kwargs.update({'height': np.log10(2*self.errorData)})
            kwargs.update({'height': 0.434*self.errorData/self.yData})
        else:
            kwargs.update({'height': 2*self.errorData})
        self.errorbars.setData(**kwargs)

    def updateItems(self):
        try:
            super(PlotDataErrorItem, self).updateItems()
        except IndexError:
            # Causes an index error in pyqtgraph built in things for
            # creating shapes. Not sure how to fix this right now
            pass
            # raise
        except np.VisibleDeprecationWarning:
            pass
        self.errorbars.setData(**self.opts)


class PlotDataErrorItemOld(pg.GraphicsObject):
    """
    This class will wrap a pg.PlotDataItem
    (which, in turn, wraps a curveItem and scatterPlotItem)
    and an errorbarItem, so they can all be wrapped together.
    Really hate that this isn't a standard feature.
    """
    sigClicked = QtCore.pyqtSignal(object)

    def __init__(self, *args, **kwargs):
        super(PlotDataErrorItem, self).__init__()

        self.setFlag(self.ItemHasNoContents)

        self.plotDataItem = pg.PlotDataItem()
        self.errorbarsItem = pg.ErrorBarItem()
        self.plotDataItem.setParentItem(self)
        self.errorbarsItem.setParentItem(self)

        self.plotDataItem.sigClicked.connect(self.curveClicked)

        # Wrap a ton of things to the plotDataItem
        for m in ['dataBounds', 'pixelPadding', 'setDownsampling',
                  'setSymbolSize', 'setSymbolBrush', 'setSymbolPen',
                  'setSymbol', 'name', 'opts', 'curve', 'scatter',
                  'opts'
                  ]:
            setattr(self, m, getattr(self.plotDataItem, m))

        self.setData(*args, **kwargs)


    def boundingRect(self):
        return QtCore.QRectF()  ## let child items handle this

    def setData(self, *args, **kwargs):
        # do I need to keep a reference to this?
        self.errorbars = np.array(kwargs.pop("errorbars", None))
        self.plotDataItem.setData(*args, **kwargs)
        kwargs.update({'x': self.plotDataItem.xData})
        kwargs.update({'y': self.plotDataItem.yData})
        kwargs.update({'height': 2*self.errorbars})
        self.errorbarsItem.setData(**kwargs)

    def setPen(self, *args, **kwargs):
        p = pg.mkPen(*args, **kwargs)
        # print("Setting pen. Color: {}, width: {}".format(p.color().name(), p.width()))
        self.plotDataItem.setPen(p)
        self.errorbarsItem.setOpts(pen=p)
        # print("After setting pens. self.width: {}".format(self.opts["pen"].width()))
        # print("After setting pens. curve.width: {}".format(self.plotDataItem.opts["pen"].width()))
        # print("self.opts is curve.opts?", self.opts is self.plotDataItem.opts)

    def setLogMode(self, xMode, yMode):
        self.plotDataItem.setLogMode(xMode, yMode)

        # errrobaritem doesn't have an implementation
        # for doing log things. Try to do some hackish
        # stuff here so it at least won't break.
        kwargs = {}
        kwargs.update({'x': self.plotDataItem.xDisp})
        kwargs.update({'y': self.plotDataItem.yDisp})
        # kwargs.update({'top': np.log10(self.errorbars+self.plotDataItem.yData) -
        #                np.log10(self.plotDataItem.yData) })
        # kwargs.update({'bottom': np.log10(self.errorbars-self.plotDataItem.yData) -
        #                np.log10(self.plotDataItem.yData) })
        kwargs.update({'height': np.log10(2*self.errorbars)})
        self.errorbarsItem.setData(**kwargs)

    def curveClicked(self):
        # print("\ncurve is clicked")
        self.sigClicked.emit(self)

    def updateItems(self):
        p = self.opts["pen"]
        # print("UpdateItemsColor: {}, width: {}".format(p.color().name(), p.width()))
        self.plotDataItem.updateItems()
        self.errorbarsItem.setData(**self.opts)

    def implements(self, interface=None):
        # pyqtgraph internals which I'd like to mimic
        ints = ['plotData']
        if interface is None:
            return ints
        return interface in ints


if __name__=='__main__':
    ex = QtGui.QApplication([])
    wid = pg.PlotWidget()
    # wid.plotItem.addLegend()


    it = PlotDataErrorItem([0, 1, 2, 3, 4], [1e2, 2e3, 1e1, 5e4, 2e2], errorbars=[1]*5)
    it.setPen('r', width=3)
    print(it.opts["pen"])
    # it = pg.PlotDataItem([0, 1, 2, 3, 4], [1, -2, 4, 8, 2])
    wid.addItem(it)

    wid.show()
    sys.exit(ex.exec_())
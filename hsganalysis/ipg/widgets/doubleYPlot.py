# -*- coding: utf-8 -*-
"""
Created on Thu Sep 25 15:41:33 2014

@author: dvalovcin
"""

import numpy as np
from PyQt5 import QtWidgets, QtCore
from ..curves.clickablePlotWidget import ClickablePlotWidget
import pyqtgraph as pg


class DoubleYPlot(ClickablePlotWidget):
    """
    Often want to have a graph which has two independent
    y axes. it's a bit more work to always handle, so have
    a simple function with conveniences for easy calling
    and manipulating

    Add linear regions to plotItem1 or p2
    """
    def __init__(self, *args, **kwargs):
        super(DoubleYPlot, self).__init__(*args, **kwargs)

        self.plotOne = self.plot()
        self.plotItem1 = self.plotItem

        #http://bazaar.launchpad.net/~luke-campagnola/pyqtgraph/inp/view/head:/examples/MultiplePlotAxes.py
        #Need to do all this nonsense to make it plot on two different axes.
        #Also note the self.updatePhase plot which shows how to update the data.
        self.p1 = self.plotItem1.vb
        self.p2 = pg.ViewBox()
        self.plotItem1.showAxis('right')
        self.plotItem1.scene().addItem(self.p2)
        self.plotItem1.getAxis('right').linkToView(self.p2)
        self.p2.setXLink(self.plotItem1)
        self.plotTwo = pg.PlotDataItem()
        self.p2.addItem(self.plotTwo)

        #need to set it up so that when the window (and thus main plotItem) is
        #resized, it informs the ViewBox for the second plot that it must update itself.
        self.plotItem1.vb.sigResized.connect(lambda: self.p2.setGeometry(self.plotItem1.vb.sceneBoundingRect()))
        self.setY1Color('k')
        self.setY2Color('r')
        print(type(self.plotOne), type(self.plotTwo))

    def setXLabel(self, label="X", units=""):
        self.plotItem1.setLabel('bottom',text=label, units=units)

    def setY1Label(self, label="Y1", units=""):
        self.plotItem1.setLabel('left', text=label, units=units, color = self.y1Pen.color().name())

    def setY2Label(self, label="Y2", units=""):
        self.plotItem1.getAxis('right').setLabel(label, units=units, color=self.y2Pen.color().name())

    def setTitle(self, title="Title"):
        self.plotItem1.setTitle(title)

    def setY1Data(self, *data):
        if len(data)==1:
            self.plotOne.setData(data[0])
        elif len(data)==2:
            self.plotOne.setData(*data)
        else:
            raise ValueError("I don't know what you want me to plot {}".format(data))

    def setY2Data(self, *data):
        if len(data)==1:
            self.plotTwo.setData(data[0])
        elif len(data)==2:
            self.plotTwo.setData(*data)
        self.p2.setGeometry(self.plotItem1.vb.sceneBoundingRect())

    def setY2Pen(self, *args, **kwargs):
        self.y2Pen = pg.mkPen(*args, **kwargs)
        self.plotItem1.getAxis("right").setPen(self.y2Pen)
        self.plotTwo.setPen(self.y2Pen)

    def setY1Color(self, color='k'):
        self.y1Pen = pg.mkPen(color)
        self.plotItem1.getAxis("left").setPen(self.y1Pen)
        self.plotOne.setPen(self.y1Pen)

    def setY1Pen(self, *args, **kwargs):
        self.y1Pen = pg.mkPen(*args, **kwargs)
        self.plotItem1.getAxis("left").setPen(self.y1Pen)
        self.plotOne.setPen(self.y1Pen)

    def setY2Color(self, color='r'):
        self.y2Pen = pg.mkPen(color)
        self.plotItem1.getAxis("right").setPen(self.y2Pen)
        self.plotTwo.setPen(self.y2Pen)







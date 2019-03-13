from __future__ import division
import numpy as np
import hsganalysis.ipg as pg
# import pyqtgraph as pg
import matplotlib.pylab as plt
import glob
import os
import json
import hsganalysis as hsg

def makePolarPlot():
    plot = pg.figure()
    plot.setAspectLocked()
    # Add polar grid lines
    plot.addLine(x=0, pen=0.2)
    plot.addLine(y=0, pen=0.2)
    for r in range(4, 38, 2):
        circle = pg.QtGui.QGraphicsEllipseItem(-r, -r, r * 2, r * 2)
        circle.setPen(pg.pg.mkPen(0.2, width=1+1*(r==18)))
        plot.addItem(circle)
    return plot


def polarPlot(*args, **kwargs):

    plt = pg.gcf()

    crv = plt.plot(*args, **kwargs)

    old = crv.setData
    def newSetData(r, theta, *args, **kwargs):

        x = r * np.cos(theta)
        y = r * np.sin(theta)
        return old(x, y, *args, **kwargs)

    crv.setData = newSetData

    # crv.setData(r, theta)

    return crv
def polarPlotBad(*args, **kwargs):

    plt = pg.gcf()

    crv = plt.plot()

    old = crv.setData
    def newSetData(self, *args, **kwargs):
        old(self, *args, **kwargs)
        x = self.xData * np.cos(self.yData)
        y = self.xData * np.sin(self.yData)

    crv.setData = newSetData

    crv.setData(*args, **kwargs)

    return crv



#
# plot = pg.plot()
# plot.setAspectLocked()
#
# # Add polar grid lines
# plot.addLine(x=0, pen=0.2)
# plot.addLine(y=0, pen=0.2)
# for r in range(2, 20, 2):
#     circle = pg.QtGui.QGraphicsEllipseItem(-r, -r, r*2, r*2)
#     circle.setPen(pg.mkPen(0.2))
#     plot.addItem(circle)
#
# # make polar data
# import numpy as np
# theta = np.linspace(0, 2*np.pi, 100)
# radius = np.random.normal(loc=10, size=100)
#
# # Transform to cartesian and plot
# x = radius * np.cos(theta)
# y = radius * np.sin(theta)
# plot.plot(x, y)









































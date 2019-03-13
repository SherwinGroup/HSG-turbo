import numpy as np
import pyqtgraph
from PyQt5 import QtCore, QtGui
from pyqtgraph.graphicsItems.ScatterPlotItem import drawSymbol
from pyqtgraph.graphicsItems.LegendItem import ItemSample
from pyqtgraph.graphicsItems.ScatterPlotItem import ScatterPlotItem
from pyqtgraph.graphicsItems.PlotDataItem import PlotDataItem
import pyqtgraph.functions as fn

# I really dislike how pyqtgraph puts its
# symbol for the legend into the center because it
# blocks the shape of the line ( you can't distinguish
# line from dashed). I like origin's much more
# where it caps the sample with the symbol instead of
# centering it

def paint(self, p, *args):
    # p.setBrush(QtGui.QBrush())
    # p.setPen(fn.mkPen('r'))
    # p.drawRect(self.boundingRect())
    # p.setRenderHint(p.Antialiasing)  # only if the data is antialiased.
    opts = self.item.opts

    if opts.get('fillLevel', None) is not None and opts.get('fillBrush',
                                                            None) is not None:
        p.setBrush(fn.mkBrush(opts['fillBrush']))
        p.setPen(fn.mkPen(None))
        p.drawPolygon(QtGui.QPolygonF(
            [QtCore.QPointF(2, 18), QtCore.QPointF(18, 2), QtCore.QPointF(18, 18)]))

    if not isinstance(self.item, ScatterPlotItem):
        p.setPen(fn.mkPen(opts['pen']))
        p.drawLine(2, 18, 18, 2)

    symbol = opts.get('symbol', None)
    if symbol is not None:
        if isinstance(self.item, PlotDataItem):
            opts = self.item.scatter.opts

        pen = fn.mkPen(opts['pen'])
        brush = fn.mkBrush(opts['brush'])
        size = opts['size']

        size = int(3*size/4)
        # save and reset to undo the forced resize that
        # drawSymbol causes
        p.save()
        p.translate(2, 18)
        path = drawSymbol(p, symbol, size, pen, brush)
        p.restore()
        p.translate(18, 2)
        path = drawSymbol(p, symbol, size, pen, brush)

ItemSample.paint = paint














































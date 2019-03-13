import numpy as np
import weakref
import pyqtgraph
from PyQt5 import QtCore, QtGui
from pyqtgraph.exporters import CSVExporter
from io import StringIO
from pyqtgraph import PlotItem

# I want to be able to copy data in the plot to the clipboard.

def export(self, fileName=None, toBytes=False, copy=False):
    if not copy:
        oldExport(self, fileName)
        return

    # Base code uses a python builtin open(), which I can't seem to figure out
    # how to monkey patch in to a StringIO object. I could try working with temporary
    # disk files, but I was running into issues that those get default opened
    # and aren't compatible with open(...). I got tired of trying to figure out
    # how to do it, so I sadly had to copy/paste his code just to remove that line,
    # then handle dumping it to the clipboard.

    # NOTE:
    # I also hardcode tab delimited so it plays better with Origin, which is where
    # I dump it half the time
    self.params['separator'] = 'tab'
    fd = StringIO()


    if not isinstance(self.item, PlotItem):
        raise Exception("Must have a PlotItem selected for CSV export.")

    data = []
    header = []

    appendAllX = self.params['columnMode'] == '(x,y) per plot'

    for i, c in enumerate(self.item.curves):
        cd = c.getData()
        if cd[0] is None:
            continue
        data.append(cd)
        if hasattr(c, 'implements') and c.implements(
                'plotData') and c.name() is not None:
            name = c.name().replace('"', '""') + '_'
            xName, yName = '"' + name + 'x"', '"' + name + 'y"'
        else:
            xName = 'x%04d' % i
            yName = 'y%04d' % i
        if appendAllX or i == 0:
            header.extend([xName, yName])
        else:
            header.extend([yName])

    if self.params['separator'] == 'comma':
        sep = ','
    else:
        sep = '\t'

    fd.write(sep.join(header) + '\n')
    i = 0
    numFormat = '%%0.%dg' % self.params['precision']
    numRows = max([len(d[0]) for d in data])
    for i in range(numRows):
        for j, d in enumerate(data):
            # write x value if this is the first column, or if we want x
            # for all rows
            if appendAllX or j == 0:
                if d is not None and i < len(d[0]):
                    fd.write(numFormat % d[0][i] + sep)
                else:
                    fd.write(' %s' % sep)

            # write y value
            if d is not None and i < len(d[1]):
                fd.write(numFormat % d[1][i] + sep)
            else:
                fd.write(' %s' % sep)
        fd.write('\n')
    # fd.close()

    fd.seek(0)
    QtGui.QGuiApplication.clipboard().setText(fd.read())



CSVExporter.allowCopy = True
oldExport = CSVExporter.export
CSVExporter.export = export











































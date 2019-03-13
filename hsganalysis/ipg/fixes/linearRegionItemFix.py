import numpy as np
import pyqtgraph
from pyqtgraph import mkPen, mkBrush
from PyQt5 import QtCore, QtGui
from ..packageSettings import config_options
from .LegendSettings_ui import Ui_LegendSettingsDialog
"""

pyqtgraph's legend item has a hard-coded slightly transparent,
grey background. Modify the class methods to allow you to change
that.

"""

oldHandler = pyqtgraph.LinearRegionItem.mouseDoubleClickEvent
def mouseDoubleClickEvent(self, *args, **kwargs):
    oldHandler(self, *args, **kwargs)
    print("I got double clicked", args, kwargs)

pyqtgraph.LinearRegionItem.mouseDoubleClickEvent = mouseDoubleClickEvent

def setLogMode(self, xMode, yMode):
    try:
        curModes = self.logMode
    except AttributeError:
        self.logMode = curModes = [False, False]
    if curModes == [xMode, yMode]: return
    self.logMode = [xMode, yMode]
    print("setLogMode called", xMode, yMode, self.angle, self.pos(), end=' ')
    if self.angle == 90:
        if xMode:
            self.setPos(np.log10(self.pos()[0]))
        else:
            self.setPos(10**(self.pos()[0]))
    elif self.angle == 0:
        if yMode:
            self.setPos(np.log10(self.pos()[1]))
        else:
            self.setPos(10**(self.pos()[1]))
    print(self.pos())
# pyqtgraph.InfiniteLine.setLogMode = setLogMode






































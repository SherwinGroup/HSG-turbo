__author__ = 'dvalovcin'
from PyQt4 import QtGui, QtCore
import pyqtgraph as pg


class DraggablePlotWidget(pg.PlotWidget):
    def __init__(self, parent=None, background='default', **kargs):
        vb = DraggableViewBox()
        kargs["viewBox"] = vb
        super(DraggablePlotWidget, self).__init__(parent, background, **kargs)


class DraggableViewBox(pg.ViewBox):
    """
    Subclassing which allows me to have a viewbox
    where I can overwrite the dragging controls
    """
    # emits (<self>, <drop pos>)
    sigDropEvent = QtCore.pyqtSignal(object, object)
    def __init__(self, parent=None, border=None, lockAspect=False, enableMouse=True, invertY=False, enableMenu=True, name=None, invertX=False):
        super(DraggableViewBox, self).__init__(parent, border, lockAspect, enableMouse, invertY, enableMenu, name, invertX)

    def mouseDragEvent(self, ev, axis=None):
        # if QtGui.QApplication.queryKeyboardModifiers() & QtCore.Qt.ShiftModifier:
        if ev.modifiers() & QtCore.Qt.ShiftModifier:
            ev.accept()
            if not ev.isFinish(): return
            self.sigDropEvent.emit(self, ev.lastScreenPos())
        else:
            super(DraggableViewBox, self).mouseDragEvent(ev, axis)

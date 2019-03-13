__author__ = 'dvalovcin'
from PyQt5 import QtWidgets, QtCore
import pyqtgraph as pg
# from hsganalysis import ipg as cpg
import hsganalysis.ipg as cpg

# class DraggablePlotWidget(pg.PlotWidget):
class DraggablePlotWidget(cpg.PlotWidget):
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

    #emits (<self>, <click pos>)
    sigClickedEvent = QtCore.pyqtSignal(object, object)

    def __init__(self, parent=None, border=None, lockAspect=False, enableMouse=True, invertY=False, enableMenu=True, name=None, invertX=False):
        super(DraggableViewBox, self).__init__(parent, border, lockAspect, enableMouse, invertY, enableMenu, name, invertX)

    def mouseDragEvent(self, ev, axis=None):
        # if QtWidgets.QApplication.queryKeyboardModifiers() & QtCore.Qt.ShiftModifier:
        if ev.modifiers() & QtCore.Qt.ShiftModifier:
            ev.accept()
            if not ev.isFinish(): return
            self.sigDropEvent.emit(self, ev.lastScreenPos())
        else:
            super(DraggableViewBox, self).mouseDragEvent(ev, axis)

    def mouseClickEvent(self, ev):
        if ev.button() == QtCore.Qt.LeftButton:
            # ev.accept()
            self.sigClickedEvent.emit(self, ev.pos())
        else:
            super(DraggableViewBox, self).mouseClickEvent(ev)

from PyQt5 import QtCore
import numpy as np
import pyqtgraph as pg


class PolarizationEllipseItem(pg.ViewBox):
    """
    A view box which can be added into a pyqtgrpah plot or image (mostly into any other view box)

    It plots a single curve which is used to show a polarization ellipse. Update vlaues with
    PolarizationEllipseItem.setEllipseCurve(alpha, gamma).

    There's a lot to be improved:
      Current positioning is done via PolarizationEllipseItem.setGeometry(), would be better
         to somehow map it to the parent window. Maybe add mouse control for moving/resizing
      Add functions for controlling the style of the curve
      Add an arrow showing direction (and change style, etc)
         Add animation for it, if desired

    """
    def __init__(self, parent=None, border=None, lockAspect=False, enableMouse=True, invertY=False, enableMenu=True, name=None, invertX=False):
        super(PolarizationEllipseItem, self).__init__(parent=parent, border=border, lockAspect=True,
                                                      enableMouse=False, invertY=invertY, enableMenu=False, name=name, invertX=invertX)


        self.ellipseCurve = pg.PlotCurveItem([], [], pen=pg.mkPen("k", width=2))

        self.addItem(self.ellipseCurve)
        self.setEllipseCurve(45, 0)
        self.setZValue(1000)

        self.setGeometry(QtCore.QRectF(0, 0, 50, 50))

    @staticmethod
    def updateEllipseValues(phi, delta):
        mag = [np.cos(phi), np.sin(phi)]
        Ex, Ey = mag
        ang = [0, delta * np.pi / 180.]
        alpha = np.arctan2(2 * mag[0] * mag[1] * np.cos(ang[1] - ang[0]),
                           (mag[0] ** 2 - mag[1] ** 2)) / 2 * 180 / 3.14159
        gamma = np.arcsin(2 * mag[0] * mag[1] * np.sin(ang[1] - ang[0]) / (
                mag[0] ** 2 + mag[1] ** 2)) / 2 * 180 / 3.14159

        return alpha, gamma
    @staticmethod
    def updateEFieldValues(alpha, gamma):
        ell = gamma * np.pi / 180.
        sampleAngle = alpha * np.pi / 180.

        vec = np.array([np.cos(ell), 1j * np.sin(ell)])[:, None]
        rot = np.array([[np.cos(sampleAngle), -np.sin(sampleAngle)],
                        [np.sin(sampleAngle), np.cos(sampleAngle)]])
        vec = np.dot(rot, vec)

        phi = np.arctan2(np.abs(vec[1]), np.abs(vec[0])) * 180. / np.pi
        delta = np.angle(vec[1], deg=True) - np.angle(vec[0], deg=True)

        return phi, delta

    def setEllipseCurve(self, alpha, gamma):
        phi, delta = self.updateEFieldValues(alpha, gamma)

        t = np.linspace(0, 2 * np.pi, 100)
        amp1 = np.cos(phi * np.pi / 180.) * np.sin(t)
        amp2 = np.sin(phi * np.pi / 180.) * np.sin(t + delta * np.pi / 180.)

        self.ellipseCurve.setData(amp1, amp2)
        return amp1, amp2

    def mouseDragEvent(self, ev, axis=None):
        ev.ignore()

    def addArrow(self, **kwargs):
        defaults = {
            "index": 24,
            "pen": "k",
            "brush": "k",
            "headLen":16
        }
        defaults.update(kwargs)
        arr = pg.CurveArrow(self.ellipseCurve, **defaults)


        return arr



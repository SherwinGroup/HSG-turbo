from PyQt5 import QtGui, QtCore, QtWidgets
import pyqtgraph as pg
import numpy as np
import collections


class PathItem(pg.GraphicsObject):

    """
    I wanted to make a fancy animation for my fan diagrams. The way people suggested
    doing it was with a QPropertyAnimation. But to do it for a graphicsitem, you have
    to make it subclass QObject. But PyQt doesn't allow this. The best I could do was
    mimic what pyqtgraph does.

    I feel like i should have actually just used a graphicsitemanimation, but I'm not
    goign to backtrack it now.
    """
    def __init__(self, *args, **kwargs):
        pg.GraphicsObject.__init__(self)
        self.item = QtWidgets.QGraphicsPathItem(*args, **kwargs)
        self.item.setParentItem(self)

    def __getattr__(self, item):
        return self.item.__getattribute__(item)

    def paint(self, *args):
        pass

    def boundingRect(self):
        return QtCore.QRectF()

def cosd(x):
    return np.cos(x * np.pi/180)

def sind(x):
    return np.sin(x * np.pi/180)

class PolarImagePlot(pg.ImageView):
    """
    ImageView subclass that will default set the ImageItem to the polar plot version
    of below.

    It'd be nice to extend it to add radial/azimuthal axes and stuff.
    """
    def __init__(self, r=None, theta=None, imageData = None, name="PolarImageView", parent=None, view=None, *args):
        item = PolarImageItem(r = r, theta=theta)
        super(PolarImagePlot, self).__init__(parent=parent, view=view, name=name, imageItem = item, *args)
        if imageData is not None:
            self.setImage(imageData)




class PolarPointInfo(object):
    def __init__(self, ridx=0, tidx=0, r=0, t=0, val=0):
        self.ridx = ridx
        self.tidx = tidx
        self.r = r
        self.t = t
        self.val = val



class PolarImageItem(pg.ImageItem):
    """
    Reimpleted image item that draws on a polar coordinate systems

    Makes use of the histogram/coloring/LUT information from a standard
    ImageItem. Some weird bug with red/blue channels, you can see the
    comment in render() if you care.

    Rendering is done by drawing a bunch of QPainterPath's, each path corresponding to
    a single pixel. I'd tried doing a few ways of making these paths
    (remnents are commented out in the render() function), but the one currently
    used is the only one seeming to still work. The paths are currently cached,
    and some amount of work needs to be done to figure out the best place to invalidate
    that cache, such as if the image is updated. However, I couldn't immediately see
    a way how to do that only when the internal image is changed, and not when the
    LUT table is changed or anything.

    The paths are draw onto a QImage, which has a compication. In the default
    ImageItem, there's a 1-to-1 of QImage pixel to underlaying data points. However,
    we need higher QImage resolution than this to get appropriate, anti-aliased curves.
    This is done by setting the _scaleFactor, but extensive testing hasn't really be done

    Note: The violation of the 1-to-1 correspondence from above is why the width(),
    height(), and paint() functions have been overridden, to compensate for the
    default ImageItem depending on the internal self.image data shape
    """

    # emits the above PolarPointInfo Class
    sigPointClicked = QtCore.pyqtSignal(object)
    def __init__(self, r, theta, image=None, **kwargs):
        # Enforce r to be monotonically increasing to handle the
        # rendering
        if np.all(np.diff(r))<0:
            r = r[::-1]
            raise NotImplementedError("You need to figure out how to shift the image data"
                                      "to account for this ")
            # image[]
        if not np.all(np.diff(r)):
            raise RuntimeError("Radial coordinate must be monotonic")

        self.r = r
        self.theta = theta
        self._scaleFactor = 1
        self._paintingPath = None
        self._paintingPathItems = None

        super(PolarImageItem, self).__init__(image, **kwargs)


        # set it so the center of the imageItem corresponds to the view
        self.translate(-self.width()/2, -self.height()/2)

        # self.getViewBox().setRenderHint(QtGui.QPainter.Antialiasing)



        self.allowMouseClicks = True
        self._previousClickObject = None
        # self.scene().sigMouseClicked.connect(self.handleMouseClicks)

    def setImage(self, image=None, autoLevels=None, **kargs):
        super(PolarImageItem, self).setImage(image, autoLevels, **kargs)

    def genPaintingPaths(self, radii, ang):
        self._paintingPath = []
        self._paintingPathItems = []
        dr = np.diff(radii)[0]  # ASSUMES MONOTONICITY/EQUAL SPACING
        dang = np.diff(ang)[0]  ## ASSUMES MONOTONICITY/EQUAL SPACING
        rect = lambda rval: QtCore.QRectF(-rval, -rval, 2 * rval, 2 * rval)
        for ridx in range(len(radii)):
            innerPaths = []
            innerItems = []
            for tidx in range(len(ang)):
                path = QtGui.QPainterPath()
                # path.arcTo(rect(radii[ridx] + dr / 2), ang[tidx] + dang / 2, -dang)
                # path.lineTo(0, 0)
                # path.arcTo(rect(radii[ridx] - dr / 2), ang[tidx] - dang / 2, dang)

                # start outer radius, outter angle
                path.arcMoveTo(rect(radii[ridx] + dr / 2), ang[tidx] + dang / 2)

                # to outter radius, inner angle
                path.arcTo(rect(radii[ridx] + dr / 2), ang[tidx] + dang / 2, -dang)

                # to inner radius, inner angle
                path.arcTo(rect(radii[ridx] - dr / 2), ang[tidx] - dang / 2, 0)

                # to inner radius, outer angle
                path.arcTo(rect(radii[ridx] - dr / 2), ang[tidx] - dang / 2, dang)

                #back to outer radius, outter angle
                path.arcTo(rect(radii[ridx] + dr / 2), ang[tidx] + dang / 2, 0)

                # outterArc = QtGui.QPainterPath()
                # outterArc.arcTo(rect(radii[ridx] + dr / 2), ang[tidx] + dang / 2, -dang)
                # innerArc = QtGui.QPainterPath()
                # innerArc.arcTo(rect(radii[ridx] - dr / 2), ang[tidx] + dang / 2, -dang)
                #
                #
                # path = outterArc.subtracted(innerArc)
                # path.setFillRule(QtCore.Qt.WindingFill)

                innerPaths.append(path)
                # item = QtWidgets.QGraphicsPathItem(path)
                item = PathItem(path)
                item.setPen(QtGui.QPen(QtCore.Qt.NoPen))

                innerItems.append(item)
                #
                # I don't know why, but when a PolarImageItem is added to a pg.ViewBox(),
                # a lot of calls to ChildrenBounds were being made, which was making
                # this function take forever. But if I just set the bounds to be ignored,
                # it resolves this issue.
                #
                # It doesn't seem to happen when the vb is an ImageView,
                self.getViewBox().addItem(item, ignoreBounds=True)

            self._paintingPath.append(innerPaths)
            self._paintingPathItems.append(innerItems)

    def width(self):
        if self.qimage is None:
            return 1
        return self.qimage.width()/self._scaleFactor

    def height(self):
        if self.qimage is None:
            return 1
        return self.qimage.height()/self._scaleFactor

    def render(self):
        if self.image is None or self.image.size == 0:
            return
        if isinstance(self.lut, collections.Callable):
            lut = self.lut(self.image)
        else:
            lut = self.lut
        argb, alpha = pg.fn.makeARGB(self.image, lut=lut, levels=self.levels)


        # I honestly don't know what's going on. the makeARGB function has a comment
        # that the R/B channels are swapped. I think I agree when playing around with
        # a straight default ImageItem that his swap is valid for that.
        # But for some reason, they shouldn't be swapped here. I don't know what's
        # going on, just that I need to unswap to get the colors to match.
        swap = argb[..., 2].copy();argb[..., 2] = argb[..., 0];argb[..., 0] = swap

        # I want to set np.nan values to be fully transparent (i.e. not rendered)
        nanargs = np.where(np.isnan(self.image))
        argb[nanargs[0], nanargs[1], 3] = 0


        # print(argb)
        # print(argb.shape)
        # return
        # qimage = QtGui.QImage(self.image.shape[0], self.image.shape[1], QtGui.QImage.Format_RGB32)

        radii = self.r
        radii = radii * self._scaleFactor
        dr = np.diff(radii)[0] # ASSUMES MONOTONICITY/EQUAL SPACING
        # print("Dr is", dr)



        ang = self.theta # * 16 # qt angles are in 1/16th degree
        dang = np.diff(ang)
        dang = np.diff(ang)[0]  ## ASSUMES MONOTONICITY/EQUAL SPACING
        rect = lambda rval: QtCore.QRectF(-rval, -rval, 2*rval, 2*rval)


        if self._paintingPath is None:
            self.genPaintingPaths(radii, ang)

        dim = int(max(abs(radii)))*2 + dr
        # qimage = QtGui.QImage(dim, dim, QtGui.QImage.Format_RGB32)
        qimage = QtGui.QImage(dim, dim, QtGui.QImage.Format_ARGB32)

        ## Todo: set this to the base color?
        qimage.fill(QtGui.QColor.fromRgbF(1, 1, 1, 0))


        painter = QtGui.QPainter(qimage)
        painter.setRenderHint(QtGui.QPainter.Antialiasing)


        # Make sure to reference to the center of the image to make things easier to
        # reference
        painter.translate(qimage.width() / 2, qimage.height() / 2)
        painter.setPen(pg.mkPen("r", width=10))

        try:
            for ridx in range(len(radii)):
                for tidx in range(len(ang)):
                    # print("make color", ridx, tidx,)
                    # print("make color", ridx, tidx, argb[ridx, tidx, :])
                    color = QtGui.QColor(*argb[ridx, tidx, :].tolist())
                    # path = self._paintingPath[ridx][tidx]
                    # painter.fillPath(path, pg.mkBrush(color))
                    item = self._paintingPathItems[ridx][tidx]
                    item.setBrush(pg.mkBrush(color))
                    item.setPen(pg.mkPen(color))
                    # item.setPen(QtCore.Qt.NoPen)
                    # if color.alpha()!=0:
                    #     painter.strokePath(path, pg.mkPen("k", width=1))
                    #     item.setPen()


        finally:
            painter.end()
        self.qimage = qimage


        # I want to make it so the center of the polar plot is the center of the viewbox...
        self.setPos(-self.width()/2, -self.height()/2)

    def paint(self, p, *args):
        """
        Need to override paint method because default argument
        fits the image into a rectangle whose dimensions match
        the image size, which isn't valid here.
        :param p:
        :param args:
        :return:
        """
        if self.image is None:
            return
        if self.qimage is None:
            self.render()
            if self.qimage is None:
                return
        if self.paintMode is not None:
            p.setCompositionMode(self.paintMode)

        # shape = self.image.shape[
        #         :2] if self.axisOrder == 'col-major' else self.image.shape[:2][::-1]
        shape = [self.width(), self.height()]
        p.drawImage(QtCore.QRectF(0, 0, *shape), self.qimage)
        if self.border is not None:
            p.setPen(self.border)
            p.drawRect(self.boundingRect())

    def mouseClickEvent(self, ev):

        if self._paintingPath is None: return
        if ev.button() != QtCore.Qt.LeftButton:
            ev.ignore()
            return
        if not self.allowMouseClicks:
            ev.ignore()
            return
        try:
            for ridx in range(len(self._paintingPath)):
                for tidx in range(len(self._paintingPath[ridx])):
                    if self._paintingPath[ridx][tidx].contains(
                            (ev.pos()+self.pos())*self._scaleFactor
                    ):
                        raise StopIteration()
            if self._previousClickObject is not None:
                self.getViewBox().removeItem(self._previousClickObject)
                self._previousClickObject = None
        except StopIteration:
            # print("Found idx", ridx, tidx, self.r[ridx], self.theta[tidx])
            ev.accept()
            if self._previousClickObject is None:
                self._previousClickObject = QtWidgets.QGraphicsPathItem()
                self.getViewBox().addItem(self._previousClickObject)
                self._previousClickObject.setBrush(pg.mkBrush("y", width=5))
                self._previousClickObject.setPen(QtGui.QPen(QtCore.Qt.NoPen))
                self._previousClickObject.setScale(1./self._scaleFactor)
            self._previousClickObject.setPath(self._paintingPath[ridx][
                                                  tidx].simplified())

            obj = PolarPointInfo(ridx, tidx, self.r[ridx], self.theta[tidx],
                                 self.image[ridx, tidx])
            self.sigPointClicked.emit(obj)
            return
        ev.ignore()

        # print("Not found\n")






    def mapToPolar(self, p):
        r = np.sqrt(p.x()**2 + p.y()**2)
        t = np.arctan2(-p.y(), p.x()) * 180/np.pi
        return r, t

pg.GraphicsScene
__author__ = 'Home'

import pyqtgraph as pg
import numpy as np

from pyqtgraph.graphicsItems.GradientEditorItem import Gradients
Gradients["rgbSpec"] = {'ticks': [
                            (0.0, (112, 0, 118, 0)),
                            (0.0001, (112, 0, 118, 255)),
                            (0.094, (108, 0, 255, 255)),
                            (0.141, (17, 0, 255, 255)),
                            (0.145, (0, 8, 255, 255)),
                            (0.266, (0, 248, 255, 255)),
                            (0.267, (0, 255, 246, 255)),
                            (0.322, (6, 255, 9, 255)),
                            (0.494, (255, 255, 0, 255)),
                            (0.663, (255, 0, 0, 255)),
                            (0.8, (255, 0, 0, 255)),
                            (1, (99, 0, 0, 255))
            ], 'mode': 'rgb'}


Gradients["alphas"] = {
            "mode": "rgb",
            "ticks": [
                (0, (0, 0, 0, 255)),
                (.25, (128, 128, 0, 255)),
                (.5, (255, 255, 255, 255)),
                (.75, (0, 128, 0, 255)),
                (1, (0, 0, 0, 255))
            ]
        }


Gradients["gammas"] = {
            "mode": "rgb",
            "ticks": [
                (0, (255, 0, 0, 255)),
                (.5, (255, 255, 255, 255)),
                (1, (0, 0, 255, 255))
            ]
        }


class FakePlotWidget(object):
    # stupid hackey way to get this to play nice with the plot container
    pass

class ImageViewWithPlotItemContainer(pg.ImageView):
    """
    I make all my uis with Qt Designer.
    I want an image view class so I can use the time-trace
    portion to display multiple ccd images as a funciton of
    "time". But it's practically necessary to have the pixel
    number labeled to compare them. THis requires passing
    the view to the imageView as a plotItem at construction

    Qt Designer doesn't allow that functionality. One option
    is to change the pyqtgrpah source so that the default
    view is a plotitem, not a viewbox, but that's super
    cumbersome to carry through to every single computer

    The alternative is to subclass it so I can force
    the view kwarg to be a plotitem, as I want

    Can set colormap with
    ['thermal', 'flame', 'yellowy', 'bipolar', 'spectrum', 'cyclic', 'greyclip', 'grey']

    """
    def __init__(self, *args, **kwargs):
        if kwargs.get("view", None) is None:
            kwargs["view"] = pg.PlotItem()
        self.plotWidget = FakePlotWidget()
        self.plotWidget.plotItem = kwargs["view"]
        super(ImageViewWithPlotItemContainer, self).__init__(*args, **kwargs)
        self.timeLine.setPen(pg.mkPen('k'))

        auto = self.ui.histogram.item.vb.menu.addAction("Autoscale Histogram")
        auto.triggered.connect(lambda: self.ui.histogram.item.imageChanged(True))

        self.view.removeItem(self.roi)
        self.roi = pg.LineSegmentROI([(0, 0), (0, 10)])
        self.roi.setZValue(20)
        self.view.addItem(self.roi)
        self.roi.hide()
        self.roi.sigRegionChanged.connect(self.roiChanged)
        self.roiCurve.setPen('k')

    def setImage(self, *args, **kwargs):
        """
        New arg features:

        if a single string is passed, it's assumed to be a CCD image file name.
        Load the file and plot it

        can pass x, y, z as kwargs and it will plot a scaled image. Works for
        nonlinear spacing in either direction

        "cmap" can set the colormap with
        ['thermal', 'flame', 'yellowy', 'bipolar',
            'spectrum', 'cyclic', 'greyclip', 'grey']

        :param args:
        :param kwargs:
        :return:
        """
        if 'z' in kwargs:
            img = kwargs.pop('z')
        elif isinstance(args[0], str):
            img = np.genfromtxt(args[0])
            if np.all(np.isnan(img)):
                img = np.genfromtxt(args[0], delimiter=',')
        elif isinstance(args[0], np.ndarray):
            args = list(args)
            img = args.pop(0)
        else:
            raise TypeError(
                "I don't know how to interpret your img arg, {}, type: {}".format(
                    arg[0], type(arg[0])
                ))
        img = img.copy()
        x, y = kwargs.pop('x', None), kwargs.pop('y', None)
        if x is not None and y is not None:
            x, y = np.array(x), np.array(y)
            if kwargs.pop("unequalSpacings", False):
                # adjust the img vector to account forunequal spacings
                newX = np.arange(x.min(), x.max(), np.diff(x).min() / 3.)
                newY = np.arange(y.min(), y.max(), np.diff(y).min() / 3.)
                newZ = np.zeros((newX.size, newY.size))

                # Define the new bins to be as wide as the sum of half the upper
                # and lower difference between neighboring points.
                # if [a, b, c, d], then pixel b would span from
                # b-(b-a)/2 to b+ (c-b)/2
                # This mimics origin's methods
                # First pixel is set as wide as the difference between
                # the first and second pixel
                diffsx = np.diff(x)[:-1]/2+np.diff(x)[1:]/2
                diffsx = np.array([0, np.diff(x)[0]]+list(diffsx))
                diffsx = np.cumsum(diffsx)
                xbins = x[0]-np.diff(x)[0]/2 + diffsx

                diffsy = np.diff(y)[:-1] / 2 + np.diff(y)[1:] / 2
                diffsy = np.array([0, np.diff(y)[0]] + list(diffsy))
                diffsy = np.cumsum(diffsy)
                ybins = y[0] - np.diff(y)[0] / 2 + diffsy

                for ii, xval in enumerate(newX):
                    for jj, yval in enumerate(newY):
                        try:
                            oldii = np.where(xval < xbins)[0][0]
                        except IndexError:
                            oldii = -1
                        try:
                            oldjj = np.where(yval < ybins)[0][0]
                        except IndexError:
                            oldjj = -1
                        newZ[ii, jj] = img[oldii, oldjj]
                img = newZ
            else:
                newX, newY = x, y

            if kwargs.get("scale", None) is None:
                kwargs["scale"] = [(newX.max()-newX.min())/newX.size, (newY.max()-newY.min())/newY.size]
                # print "scaled to", kwargs["scale"]
            if kwargs.get("pos", None) is None:
                kwargs["pos"] = [newX.min(), newY.min()]
                # print "moved to", kwargs["pos"]

        else:
            # if you pass a 3D image to pg.ImageView
            # of img, with data.shape=(L, M, N)
            # with N<=4, he interprets the axes differently
            # than N>4. I don't want to do this
            axes = kwargs.get("axes", None)
            if axes is None:
                kwargs["axes"] = {'t': 0, 'x': 1, 'y': 2, 'c': None}


            # I was having a lot of trouble getting ImageView
            # to accept an axes argument. That is, it was just ignoring it.
            # I couldn't figure out how to handle it, so this convenience
            # kwarg will transpose the image from my typical
            # considerations (time, x, y) to his bizarre standard,
            # (time, y, x)

            resize = kwargs.pop("transpose", True)
            if resize:
                if img.ndim == 3:
                    img = np.transpose(img, (0, 2, 1))
                elif img.ndim == 2:
                    img = img.T
                    # need to set it to 3 axes to deal with
                    # setting the axes as above

                    img = img[None,:,:]
        cmap = kwargs.pop("cmap", "rgbSpec")
        # print("plot args", args, kwargs)
        super(ImageViewWithPlotItemContainer, self).setImage(img, *args, **kwargs)
        if cmap is not None:
            self.ui.histogram.item.gradient.loadPreset(cmap)
        # For some reason, it hasn't been updating the color levels.
        self.ui.histogram.item.imageChanged(True)

    def roiChanged(self):
        if self.image is None:
            return

        image = self.getProcessedImage()
        if image.ndim == 2:
            axes = (0, 1)
        elif image.ndim == 3:
            axes = (1, 2)
        else:
            return
        data = self.roi.getArrayRegion(image.view(np.ndarray), self.imageItem, axes)[0]

        self.roiCurve.setData(y=data)

    def xlim(self, xmin=None, xmax=None):
        self.view.setXRange(xmin, xmax, padding=0)

    def ylim(self, ymin=None, ymax=None):
        self.view.setYRange(ymin, ymax, padding=0)

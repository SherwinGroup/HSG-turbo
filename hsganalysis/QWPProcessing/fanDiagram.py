import numpy as np
from PyQt5 import QtGui, QtCore, QtWidgets
from hsganalysis.ipg import PolarImageItem, PolarAxis
from pyqtgraph import GraphicsView, HistogramLUTWidget, ViewBox, mkColor, TextItem, Point
from .extractMatrices import *
from .expFanCompiler import *


class FanDiagram(QtGui.QWidget):
    """
    Produces a fan diagram
    alphas/gammas should be the matrices saved from the FanCompiler, of form:

    arb     | niralpha1 | niralpha2 | niralpha3 | niralpha4 | ...
    1st sb  | 1sb alpha | 1sb alpha |     .
    2nd sb  | 2sb alpha | 2sb alpha |     .
    3rd sb  | 3sb alpha | 3sb alpha |     .
      .
      .
      .

    Assumes both alphas/gammas are the same shape

    Alternatively, pass nirAlpha and SBs as 1st/2nd args (as 1D arrays) to have it
    create the fan without any data

    Or, if alphaData and gammaData are strings, assumes they're paths to data files to plot
    :param alphas:
    :param gammas:
    :param kwargs:
    :return:
    """
    def __init__(self, alphaData, gammaData=None, view=None):
        super(FanDiagram, self).__init__()
        if gammaData is None and isinstance(alphaData, FanCompiler):
            alphaData, gammaData = alphaData.build(withErrors=False)[:2]

        self.layout = QtWidgets.QHBoxLayout()

        self.histAlpha = HistogramLUTWidget(self)
        self.centralView = GraphicsView()
        self.histGamma = HistogramLUTWidget(self)

        self.histAlpha.setMinimumWidth(150)
        self.histGamma.setMinimumWidth(150)

        self.layout.addWidget(self.histGamma)
        self.layout.addWidget(self.centralView)
        self.layout.addWidget(self.histAlpha)
        self.layout.setContentsMargins(0,0,0,0)
        self.layout.setSpacing(0)
        self.setLayout(self.layout)

        if view is None:
            self.view = ViewBox()
        else:
            self.view = view


        self.centralView.setCentralItem(self.view)
        self.view.setAspectLocked(True)
        self.view.invertY(True)

        if isinstance(alphaData, str) and isinstance(gammaData, str):
            alphaData = np.genfromtxt(alphaData, delimiter=',')
            gammaData = np.genfromtxt(gammaData, delimiter=',')

        if alphaData.ndim == gammaData.ndim == 1:
            # Assume you just want it to be created, and will later populate it
            nirAlphas = alphaData
            sbs = gammaData

            alphaData = np.ones((sbs.shape[0] + 1, nirAlphas.shape[0] + 1)) * -1
            alphaData[1:, 0] = sbs
            alphaData[0, 1:] = nirAlphas

            gammas = np.ones((sbs.shape[0] + 1, nirAlphas.shape[0] + 1)) * -1
            gammas[1:, 0] = sbs
            gammas[0, 1:] = nirAlphas

        sbs = alphaData[1:, 0]
        maxSB = sbs.max()
        nirAlphas = alphaData[0, 1:]


        self.alphaItem = PolarImageItem(r=sbs, theta=nirAlphas)
        self.alphaItem.setImage(alphaData[1:,1:])
        # nirAlphas+180 is what causes the gamma angles to appear on the left side of the
        # fan. This seemed easier than doing some sort of coordinate inversion/flipping
        # on the plot itself.
        self.gammaItem = PolarImageItem(sbs, nirAlphas+180, gammaData[1:,1:])

        self.view.addItem(self.alphaItem)
        self.view.addItem(self.gammaItem)

        self.histAlpha.setImageItem(self.alphaItem)
        self.histGamma.setImageItem(self.gammaItem)

        # manually set the default state to the black-gold-white-green-black. Not sure
        # if it's necessary to have this be a free parameter vs being hardcoded
        self.histAlpha.gradient.restoreState({
            "mode": "rgb",
            "ticks": [
                (0, (0, 0, 0, 255)),
                (.25, (128, 128, 0, 255)),
                (.5, (255, 255, 255, 255)),
                (.75, (0, 128, 0, 255)),
                (1, (0, 0, 0, 255))
            ]
        })

        # Set the default spacings for the alpha color axis. Again, not sure if it's
        # better to leave the 18pt font hard-coded or not, but I am
        self.histAlpha.axis.setTickFont(QtGui.QFont("Arial", 18))
        self.histAlpha.axis.setTickSpacing(30, 15)
        self.histAlpha.axis.setLabel("&alpha; (&deg;)", **{'font-family': 'Times',
                                                           "font-size":  "18pt"})

        # As with alpha, hard code the initial color space for gamma (blue-white-red)
        # and the font spacings and stuff
        self.histGamma.gradient.restoreState({
            "mode": "rgb",
            "ticks": [
                (0, (255, 0, 0, 255)),
                (.5, (255, 255, 255, 255)),
                (1, (0, 0, 255, 255))
            ]
        })
        self.histGamma.axis.setTickFont(QtGui.QFont("Arial", 18))
        self.histGamma.axis.setTickSpacing(15, 15)
        self.histGamma.axis.setLabel("&gamma; (&deg;)", **{'font-family': 'Times',
                                                           "font-size":  "18pt"})

        self.histAlpha.item.setLevels(-90, 90)
        self.histGamma.item.setLevels(-45, 45)

        self.histAlpha.autoHistogramRange()
        self.histGamma.autoHistogramRange()

        # Make it the right dimensions, making sure that the width is appropriate.
        # This makes it easier to automate plotting/saving fans and making sure
        # their dimensions are consistent.
        g = self.geometry()
        # I found these by eye, there's not very much important about them
        g.setWidth(773)
        g.setHeight(480)
        # Manually center it on the screen, since geometry isn't well defined at this point
        # before events are processed
        g.moveCenter(QtWidgets.QApplication.desktop().screenGeometry().center())
        self.setGeometry(g)

        # Add in the radial axes for it
        self.axes = {
            "radial": PolarAxis("radial"),
            "azimuthal": PolarAxis("azimuthal")
        }
        # Lighten the radial font to make it distinct from the other
        p = self.axes["radial"].pen()
        p.setColor(mkColor("#666666"))
        self.axes["radial"].setPen(p)

        for a in self.axes.values():
            # Make sure the axes sit on top of all other items
            a.setZValue(10000)
            # make sure that they scale appropriately, instead of just floating on top
            a.linkToView(self.view)
            # Ignore bounds prevents the window from resizing to try and fit in
            # the axes items
            self.addItem(a, ignoreBounds=True)

        # manually set the positions and string values for alpha angles. [-90, 90] work
        # well. The other half needs the +-180 to make sure the gamma angles have the
        # correctly labeled with respect to alpha_nir
        self.axes["azimuthal"].setTicks(
            [
                [(ii, str(ii)) for ii in np.arange(-90, 91, 30)] + # alpha side (Q1+Q4)
                [(ii, str(ii + 180)) for ii in np.arange(-180, -91, 30)] + #Q3
                [(ii, str(ii - 180)) for ii in np.arange(120, 151, 30)], #Q1
            ]
        )

        # add a title (without text)
        self.titleItem = TextItem()
        self.titleItem.setAnchor(Point(0.5, 1)) # anchor on bottom-center
        # Again, not sure if it's necessary to have the font color/size being
        # a free parameter
        self.titleItem.setColor("k")
        self.titleItem.setFont(QtGui.QFont("Arial", 15))
        # Ignore bounds so that the view won't try to account for it (which
        # causes a conflict because the title is placed with respect to the
        # view region)
        self.view.addItem(self.titleItem, ignoreBounds=True)

        self.show()
        # Arbitrary forcing updates to try and track down why some things don't
        # update correctly
        QtWidgets.QApplication.processEvents()
        self.view.updateViewRange(True, True)

    def setAlphaImage(self, img):
        self.alphaItem.setImage(img)

    def setGammaImage(self, img):
        self.gammaItem.setImage(img)

    def setImages(self, alpha, gamma):
        self.setAlphaImage(alpha)
        self.setGammaImage(gamma)

    def export(self, fname, hideHistograms=True,
               pngScale = 4):
        """
        Save fan diagrams to file, with the full image, and color bars on the alpha/gamma
        values
        :param fname: the fname to save as
        hideHistograms - (True) Prevent rendering the histograms, often ncier for
            figures/presentations

        If fname.endswith(".svg"), it outputs as an SVG. Howver, it's not the cleanest
            thing (the files are quite large/unoptimized, and I can't think of an
            easy way to correct that). Also, when the svg is converted to pdf via
            Inkscape, things get fucked up for some reasons (axes get thicker, fonts
            get borked, pixels get messed up). So, it kinda works, but there's
            stuff seriously wrong.

        One thing to make things cleaner is to use this site:
            https://jakearchibald.github.io/svgomg/
        which optimizies the svg and makes it a lot easier to work with
        :return:
        """

        # defaults = {
        #     "hideHistograms": False
        # }
        #
        # defaults.update(kwargs)

        doSvg = fname.endswith(".svg")

        if hideHistograms:
            # Hide the histogram data (and shrink the plot)
            # to avoid confusing people
            self.histAlpha.plot.hide()
            self.histAlpha.vb.setMaximumWidth(20)
            self.histGamma.plot.hide()
            self.histGamma.vb.setMaximumWidth(20)


        QtWidgets.QApplication.processEvents()

        self.histGamma.axis.setRange(-46.75, 46.75)
        self.histAlpha.axis.setRange(-94, 94)
        width, height = self.width(), self.height()


        if doSvg:
            from PyQt5 import QtSvg
            outputImage = QtSvg.QSvgGenerator()
            outputImage.setFileName(fname)
            outputImage.setSize(QtCore.QSize(int(width), int(height)))
            # I'm not sure why it has to be this, but the axis on the histogrm
            # were fuckingup without it
            outputImage.setResolution(96)
        else:
            outputImage = QtGui.QImage(width * pngScale, height * pngScale,
                                       QtGui.QImage.Format_ARGB32)
            # outputImage.setDotsPerMeterX(650 * 100 / 2.54)
            # outputImage.setDotsPerMeterY(650 * 100 / 2.54)
            # this gives a moderatly high quality image
            outputImage.setDevicePixelRatio(pngScale)
            outputImage.fill(QtGui.QColor("white"))

        outputPainter = QtGui.QPainter(outputImage)

        self.render(outputPainter)

        if not doSvg:
            ret = outputImage.save(fname)

        outputPainter.end()

    def addItem(self, item, ignoreBounds=False):
        self.view.addItem(item, ignoreBounds)

    def setViewRadius(self, r):
        # Set the view range of the fan diagram such that radius r is visible
        self.view.setRange(QtCore.QRect(-r, -r, 2*r, 2*r), padding=0)

    def hideHistogramAxes(self, hideTicks=True):
        # Hide the histogram region item and plots and all that for
        # less cluttered plots. Definitely useful if export is called with
        # hideHistograms=True, where the regions are useless.

        # Hide the linear regions
        self.histGamma.region.hide()
        self.histAlpha.region.hide()
        # Keep a reference to the old paint methods so you can reverse it if desired
        # This stops the painting of the bars which go from the linear region to the
        # gradient editor
        self.histGamma.item.oldPaint = self.histGamma.item.paint
        self.histAlpha.item.oldPaint = self.histAlpha.item.paint
        # Overwriting the functions to return None causes all the other rendering
        # things to abort
        self.histGamma.item.paint = lambda *x: None
        self.histAlpha.item.paint = lambda *x: None

        if hideTicks:
            # Hide the ticks which can be used for changing the stops/colors of
            # the gradients, which are rather ugly
            # Note: Since this only hides ticks which are present, I don't think
            [ii.hide() for ii in self.histAlpha.item.gradient.ticks.keys()]
            [ii.hide() for ii in self.histGamma.item.gradient.ticks.keys()]

        QtWidgets.QApplication.processEvents()
        # Hard coded numbers which make it look like the axes values line up with
        # the gradient item, which is more in-line with how color bars are interpreted
        self.histGamma.axis.setRange(-46.75, 46.75)
        self.histAlpha.axis.setRange(-94, 94)

    def showHistogramAxes(self, showTicks=True):
        try:
            self.histGamma.item.paint = self.histGamma.item.oldPaint
            self.histAlpha.item.paint = self.histAlpha.item.oldPaint
            del self.histAlpha.item.oldPaint
            del self.histGamma.item.oldPaint
        except AttributeError:
            # You didn't hide them first (or at least not here
            return
        self.histGamma.region.show()
        self.histAlpha.region.show()


        if showTicks:
            [ii.show() for ii in self.histAlpha.item.gradient.ticks.keys()]
            [ii.show() for ii in self.histGamma.item.gradient.ticks.keys()]

    @staticmethod
    def fromTMatrix(tMatrix, angle = 45, sbs=None):
        """
        Create a fan diagram from T matrices directly. The angle needs to be specified
        so the T matrix can be converted to a J matrix. The angle is relative to what's
        specified in the qwp.extractMatrices.makeU function.

        if you pass a string, it assumes it's a file name from a saved one. It'll load
        that and plot it. If you also pass values to sbs, it'll make sure only the passed
        values are plotted. Otherwise, it'll plot all the sbs in the file

        If you pass a tMatrix as returned from the fitting routines, you also need
        to pass the sbs directly in this case, since the tMatrices don't include them.
        :param tMatrix:
        :param angle:
        :param sbs:
        :return:
        """
        if isinstance(tMatrix, str):
            # a file is passed
            if sbs is not None:
                # Pass an array of sbs with a string, and this'll parse
                # out the sidebands which aren't included in  the passed array
                wantsbs = sbs
            else:
                wantsbs = None
            tMatrix, sbs = loadT(tMatrix)
            # Handle if only a select number of sidebands is specified
            if wantsbs is not None:
                try:
                    # Find the indices of the desired sidebands within the array of
                    # sidebands actually loaded
                    wantIdx = [sbs.tolist().index(ii) for ii in wantsbs]
                    # Cull it to only those specified
                    sbs = sbs[wantIdx]
                    # tMatrix is multidimensional (tMatrix.ndim>2), so ellipses cut
                    # out the other axes
                    tMatrix = tMatrix[..., wantIdx]
                    # Ensure that you got everything you want. Could happen if sidebands
                    # are requested (passed to the function) and not found
                    assert np.all(wantsbs == sbs)
                except ValueError as e:
                    raise IndexError("Invalid sideband requested ({} is not in loaded)".format(
                        e.args[0].split(' ')[0]
                    ))
                except AssertionError:
                    raise IndexError("Invalid sideband requested")


        jMatrix = makeJfromT(tMatrix, angle)
        if sbs is None:
            raise RuntimeWarning("Desired sidebands to plot should be specified as kwarg sbs")
            sbs = np.arange(8, 38, 2)

        alpha, gamma = jonesToFans(sbs = sbs, J=jMatrix)
        return FanDiagram(alpha, gamma)

    def setTitle(self, title="", adjustBounds=True):
        """
        Sets the title of the fan diagram, positioning the text right above the center
        of the fan
        :param title:
        :param adjustBounds:
        :return:
        """
        self.titleItem.setText(title)

        # Move the title so the bottom is at the top of the outer axis
        self.titleItem.setPos(0, self.axes["azimuthal"].fullBoundingRect.top())
        QtWidgets.QApplication.processEvents()
        # Double up because of some weird fucking issue with Qt not appropriately
        # updating things when requested
        self.titleItem.setPos(0, self.axes["azimuthal"].fullBoundingRect.top())
        QtWidgets.QApplication.processEvents()
        # print(self.titleItem.mapRectToView(self.titleItem.boundingRect()))

        if adjustBounds:
            # Readjust the viewbox to frame the fan better

            # Find the top, based on the coordinates of the top of the title
            top = self.titleItem.mapRectToView(self.titleItem.boundingRect()).top()
            # Bottom is defiend by the bottom of the axes (includes the text)
            # Note: this assumes the
            bottom = self.axes["azimuthal"].fullBoundingRect.bottom()
            # print("bottom", bottom)
            w = abs(top-bottom)
            # print("new rect", QtCore.QRectF(-w/2, top, w, w))
            self.view.setRange(QtCore.QRectF(-w/2, top, w, w), padding=0)
            self.view.setRange(QtCore.QRectF(-w/2, top, w, w), padding=0)
            # self.view.update()

    def setMaxRadius(self, radius=40):
        # Set the maximum value for both of the axes to the value specified.
        # The 1e-6 is to prevent it from producing an "r=0" label and stuff
        self.axes["azimuthal"]._bounds["radial"] = [1e-6, radius]
        self.axes["radial"]._bounds["radial"] = [1e-6, radius]

        # Need to invalidate the cache for the axes, forcing it to redraw and update
        # the bounding rect and stuff
        self.axes["azimuthal"].picture = None
        self.axes["radial"].picture = None



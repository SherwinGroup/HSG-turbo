import numpy as np

from PyQt5 import QtGui, QtCore, QtWidgets
from interactivePG import PolarImageItem, PolarAxis
from pyqtgraph import GraphicsView, HistogramLUTWidget, ViewBox, mkColor


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
    def __init__(self, alphaData, gammaData, view=None):
        super(FanDiagram, self).__init__()

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

        # self.alphaItem = PolarImageItem(sbs, nirAlphas, alphaData[1:,1:])
        self.alphaItem = PolarImageItem(r=sbs, theta=nirAlphas)
        self.alphaItem.setImage(alphaData[1:,1:])
        self.gammaItem = PolarImageItem(sbs, nirAlphas+180, gammaData[1:,1:])

        self.view.addItem(self.alphaItem)
        self.view.addItem(self.gammaItem)

        self.histAlpha.setImageItem(self.alphaItem)
        self.histGamma.setImageItem(self.gammaItem)

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

        self.histAlpha.axis.setTickFont(QtGui.QFont("Arial", 18))
        self.histAlpha.axis.setTickSpacing(30, 15)
        self.histAlpha.axis.setLabel("&alpha; (&deg;)", **{'font-family': 'Times',
                                                           "font-size":  "18pt"})

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

        # Make it the rightd imensions
        g = self.geometry()
        g.setWidth(773)
        g.setHeight(480)
        # Manually center it on the screen, since geometry isn't well defined at this point
        # before events are processed
        g.moveCenter(QtWidgets.QApplication.desktop().screenGeometry().center())
        self.setGeometry(g)


        self.axes = {
            "radial": PolarAxis("radial"),
            "azimuthal": PolarAxis("azimuthal")
        }
        # Lighten the radial font to make it distinct from the other
        p = self.axes["radial"].pen()
        p.setColor(mkColor("#666666"))
        self.axes["radial"].setPen(p)

        for a in self.axes.values():
            a.setZValue(10000)
            a.linkToView(self.view)
            self.addItem(a, ignoreBounds=True)

        self.axes["azimuthal"].setTicks(
            [
                [(ii, str(ii)) for ii in np.arange(-90, 91, 30)] + # alpha side (Q1+Q4)
                [(ii, str(ii + 180)) for ii in np.arange(-180, -91, 30)] + #Q3
                [(ii, str(ii - 180)) for ii in np.arange(120, 151, 30)], #Q1
            ]
        )


    def export(self, fname, hideHistograms=True):
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
            outputImage = QtGui.QImage(width * 4, height * 4,
                                       QtGui.QImage.Format_ARGB32)
            # outputImage.setDotsPerMeterX(650 * 100 / 2.54)
            # outputImage.setDotsPerMeterY(650 * 100 / 2.54)
            # this gives a moderatly high quality image
            outputImage.setDevicePixelRatio(4)
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
        self.view.setRange(QtCore.QRect(-r, -r, 2*r, 2*r))

    def hideHistogramAxes(self):
        # Hide the histogram region item and plots and all that for
        # less excessive plots
        self.histGamma.region.hide()
        self.histGamma.item.oldPaint = self.histGamma.item.paint
        self.histGamma.item.paint = lambda *x: None
        self.histAlpha.region.hide()
        self.histAlpha.item.oldPaint = self.histAlpha.item.paint
        self.histAlpha.item.paint = lambda *x: None

        # Hard coded numbers which make it look like the axes values line up with
        # the gradient item
        self.histGamma.axis.setRange(-46.75, 46.75)
        self.histAlpha.axis.setRange(-94, 94)

    def showHistogramAxes(self):
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


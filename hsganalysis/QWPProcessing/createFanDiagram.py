import numpy as np
from PyQt5 import QtCore, QtGui, QtWidgets
from hsganalysis.ipg import PolarImagePlot
import hsganalysis.ipg as pg


"""
12/17/18
This is an older version of making fan diagrams, it would create indpedent polarImagePlots
and then add them on top of each other, and things were a little gross.

the fanDiagram class is a little bit cleaner, keeping everything together and
is a bit tidier. At least, it makes much nicer plots for figures
"""


def createFan(alphas, gammas,
              plotEllipseHelpers=True,
              showGamma=True,
              showCenterEllipse=True,
              showInfoText=True,
              **kwargs):
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
    :param alphas:
    :param gammas:
    :param kwargs:
    :return:
    """

    defaults = {
        "plotEllipseHelpers": plotEllipseHelpers,
        "showGamma": showGamma,
        "showCenterEllipse": showCenterEllipse,
        "showInfoText": showInfoText,
    }
    defaults.update(kwargs)

    alphas = np.array(alphas)
    gammas = np.array(gammas)

    if alphas.ndim == gammas.ndim == 1:
        # Assume you just want it to be created, and will later populate it
        nirAlphas = alphas
        sbs = gammas

        alphas = np.ones((sbs.shape[0] + 1, nirAlphas.shape[0] + 1)) * -1
        alphas[1:, 0] = sbs
        alphas[0, 1:] = nirAlphas

        gammas = np.ones((sbs.shape[0] + 1, nirAlphas.shape[0] + 1)) * -1
        gammas[1:, 0] = sbs
        gammas[0, 1:] = nirAlphas


    sbs = alphas[1:, 0]
    maxSB = sbs.max()
    nirAlphas = alphas[0, 1:]

    # p(olarplot)alp(ha)
    palp = PolarImagePlot(r=sbs, theta=nirAlphas, imageData=alphas[1:, 1:])

    # Set the color on the la
    palp.ui.histogram.gradient.restoreState({
        "mode": "rgb",
        "ticks": [
            (0, (0, 0, 0, 255)),
            (.25, (128, 128, 0, 255)),
            (.5, (255, 255, 255, 255)),
            (.75, (0, 128, 0, 255)),
            (1, (0, 0, 0, 255))
        ]
    })
    palp.ui.histogram.axis.setTickFont(QtGui.QFont("Arial", 18))
    palp.view.setXRange(-maxSB*1.1, maxSB*1.1)
    palp.view.setYRange(-maxSB*1.1, maxSB*1.1)
    palp.setLevels(-90, 90)

    # p(olarplot)gam(ma)
    pgam = PolarImagePlot(r=sbs, theta=180 + nirAlphas, imageData=gammas[1:, 1:])
    pgam.ui.histogram.gradient.restoreState({
        "mode": "rgb",
        "ticks": [
            (0, (255, 0, 0, 255)),
            (.5, (255, 255, 255, 255)),
            (1, (0, 0, 255, 255))
        ]
    })
    pgam.ui.histogram.axis.setTickFont(QtGui.QFont("Arial", 18))
    pgam.setLevels(-45, 45)
    palp.addItem(pgam.imageItem)

    palp.ui.histogram.axis.setTickSpacing(30, 15)
    palp.ui.histogram.axis.setLabel("&alpha; (&deg;)", **{'font-family': 'Times',
                                                          "font-size":  "18pt"})

    pgam.ui.histogram.axis.setTickSpacing(15, 5)
    pgam.ui.histogram.axis.setLabel("&gamma; (&deg;)", **{'font-family': 'Times',
                                                          "font-size":  "18pt"})

    # For some reason, this is important. It doesn't seem to re-render it properly
    # showing the full gamma range.
    pgam.ui.histogram.autoHistogramRange()
    palp.ui.histogram.autoHistogramRange()


    #p(olarziation)e(llipse)
    pe = pg.PolarizationEllipseItem()
    pe.setGeometry(-7, -7, 14, 14)
    pe.setEllipseCurve(45, 45)
    textItem = pg.pg.TextItem("Start", color=(0, 0, 0))
    palp.textItem = textItem
    palp.pe = pe
    if not defaults["showCenterEllipse"]:
        pe.hide()
    if not defaults["showInfoText"]:
        textItem.hide()

    palp.addItem(textItem)
    textItem.setPos(-maxSB*1.1, -maxSB)
    palp.imageItem.sigPointClicked.connect(
        lambda x: textItem.setHtml(
            f"r={x.r:.0f}, &theta;={x.t:.0f},<br>f(r, &theta;)={x.val:.3f}"
        )
    )

    def updateCurve(info):
        # a = alphaData[info.ridx, info.tidx]
        # gamma = gammaData[info.ridx, info.tidx]
        a = palp.imageItem.image[info.ridx, info.tidx]
        gamma = pgam.imageItem.image[info.ridx, info.tidx]
        pe.setEllipseCurve(a, gamma)
        textItem.setHtml(
            f"r={info.r:.0f}, &theta;={info.t:.0f},<br>"
            f"  &alpha;(r, &theta;)={a:.3f}<br>  &gamma;(r, &theta;)={gamma:.3f}"
        )

    palp.imageItem.sigPointClicked.connect(updateCurve)
    pgam.imageItem.sigPointClicked.connect(updateCurve)
    palp.addItem(pe)

    if defaults["plotEllipseHelpers"]:
        # Add helper lines to indicate the polarization states

        pgam.ui.histogram.axis.setWidth(85)
        palp.ui.histogram.axis.setWidth(85)
        for a in [-90, -60, -30, 30, 60, 90]:
            e = pg.PolarizationEllipseItem()
            palp.ui.histogram.addItem(e)
            e.setEllipseCurve(a, 0)

            # Hardcoded guess-and-check
            e.setPos(0, 20 + 78*(90-a)/30)
            # e.setScale(0.7)

            # Stupid pos curve item has some weird race conditions, and it rarely
            # orients the arrow correctly. So, disable it and set the rotations
            # manually.
            arr = e.addArrow(rotate=False)
            arr.setIndex(24)
            arr.rotate(-2*a)

            arr = e.addArrow()
            arr._rotate = False
            arr.setIndex(74)
            arr.rotate(-a)

        for g in [45, 30, 15, -15, -30, -45]:
            e = pg.PolarizationEllipseItem()
            pgam.ui.histogram.addItem(e)
            e.setEllipseCurve(0, g)
            # Hardcoded guess-and-check
            e.setPos(-0, 0 + 170*(45-g)/30)
            # e.setScale(0.7)


            arr = e.addArrow(rotate=False)
            arr.setIndex(0)
            # arr.rotate(-2*a)

            arr = e.addArrow()
            arr._rotate = False
            arr.setIndex(50)



    palp.show()
    if defaults["showGamma"]:
        pgam.show()
    return palp, pgam

def saveAndRenderFan(p1, p2, fname, hideHistograms=False):
    """
    Save fan diagrams to file, with the full image, and color bars on the alpha/gamma
    values
    :param p1: palp coming from the createFan function above
    :param p2: pgam coming from the createFan function above
    :param fname: the fname to save as

    potential kwargs:
        hideHistograms - (True) Prevent rendering the histograms
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
        p1.ui.histogram.plot.hide()
        p1.ui.histogram.vb.setMaximumWidth(20)
        p2.ui.histogram.plot.hide()
        p2.ui.histogram.vb.setMaximumWidth(20)

    hist1 = p1.ui.histogram.scene().getViewWidget()
    hist2 = p2.ui.histogram.scene().getViewWidget()
    center = p1.scene.getViewWidget()

    r1 = QtCore.QRectF(hist1.viewRect())
    r2 = QtCore.QRectF(hist2.viewRect())
    rc = QtCore.QRectF(center.viewRect())

    QtWidgets.QApplication.processEvents()

    width = r1.width() + r2.width() + rc.width()

    height = rc.height()
    if doSvg:
        from PyQt5 import QtSvg
        outputImage = QtSvg.QSvgGenerator()
        outputImage.setFileName(fname)
        outputImage.setSize(QtCore.QSize(int(width), int(height)))

        outputImage.setResolution(96) # I'm not sure why it has to be this...
    else:
        outputImage = QtGui.QImage(width*4, height*4, QtGui.QImage.Format_ARGB32)
        # outputImage.setDotsPerMeterX(650 * 100 / 2.54)
        # outputImage.setDotsPerMeterY(650 * 100 / 2.54)
        outputImage.setDevicePixelRatio(4)
        outputImage.fill(QtGui.QColor("white"))

    outputPainter = QtGui.QPainter(outputImage)

    
    r2.setHeight(height)
    r1.setHeight(height)



    hist2.render(outputPainter, r2)
    rc.moveLeft(r2.width())
    center.render(outputPainter, rc)
    r1.moveLeft(rc.width()+r2.width())
    hist1.render(outputPainter, r1)

    if not doSvg:
        ret = outputImage.save(fname)

    outputPainter.end()
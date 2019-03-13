from PyQt5 import QtCore, QtWidgets
import time
import numpy as np
import pyqtgraph as pg


class DateAxis(pg.AxisItem):
    def __init__(self, *args, **kwargs):
        """
        Override the init so that, if you pass a plotWidget,
        or some larger plot container (plotWidget, plotItem, plotWindow, etc)
        it will put itself in the proper place. There's some hackey-ness to make
        this work, and you should _probably_ make the DateAxis and pass
        it as an argument to the plotItem constructor (which pyqtgraph
        supports), but then... I dunno, maybe I should actually
        change this code to handle that, but then how would I specify
        both top and bottom being DateAxes?

        Example usage:
        plot = pg.plot()

        pg.DateAxis(plot, orientation="bottom")

        :param args:
        :param kwargs:
        """

        ## Try to parse if the first argument is some plotItem-like
        ## thing to try and add myself to
        if len(args) == 0 or isinstance(args[0], str):
            pi = None
        elif isinstance(args[0], pg.PlotItem):
            pi = args[0]
            args = list(args)
            args = tuple(args[1:])
        else:
            try:
                pi = args[0].plotItem
                args = list(args)
                args = tuple(args[1:])
            except:
                raise RuntimeError("I don't know how to parse this argument,",
                                   type(args[0]), args[0])
        ## Call super only after we've removed a potentially conflicting
        ## arg
        super(DateAxis, self).__init__(*args, **kwargs)
        self.autoSIPrefix = False ## DateAxis shouldn't need SI units
        if pi is not None:
            ## Handle replacing this axis with
            ## what was sepcified there.
            ## Need to remove the initial item that plotItem
            ## originally created
            pi.layout.removeItem(pi.getAxis(self.orientation))
            self.setParentItem(pi)
            ## link to scale the axis properly
            self.linkToView(pi.vb)
            pi.axes[self.orientation]["item"] = self
            pi.layout.addItem(self, *pi.axes[self.orientation]["pos"])
            ## The following two lines are set in the plotItem
            ## init when added axes. I'm not 100% sure what
            ## it's for
            # self.setZValue(-1000)
            # self.setFlag(self.ItemNegativeZStacksBehindParent)

    def tickStrings(self, values, scale, spacing):
        """
        Overwrite the default tickStrings to return
        the values in clock-time
        :param values:
        :param scale:
        :param spacing:
        :return:
        """
        strns = []
        for x in values:
            try:
                strns.append(time.strftime("%X", time.localtime(x)))
            except (ValueError, OSError):  ## Windows can't handle dates before 1970
                strns.append('')
        return strns

    def drawPicture(self, p, axisSpec, tickSpecs, textSpecs):
        """
        Rotate the tick text labels to make them fit a little
        better.
        :param p:
        :type p: QtGui.QPainter
        :param axisSpec:
        :param tickSpecs:
        :param textSpecs:
        :return:
        """

        p.setRenderHint(p.Antialiasing, False)
        p.setRenderHint(p.TextAntialiasing, True)

        ## draw long line along axis
        pen, p1, p2 = axisSpec
        p.setPen(pen)
        p.drawLine(p1, p2)
        p.translate(0.5,0)  ## resolves some damn pixel ambiguity

        ## draw ticks
        for pen, p1, p2 in tickSpecs:
            p.setPen(pen)
            p.drawLine(p1, p2)

        ## Draw all text
        if self.tickFont is not None:
            p.setFont(self.tickFont)
        p.setPen(self.pen())
        ## TODO: Figure this out, I want this text rotate, ffs.
        for rect, flags, text in textSpecs:
            p.save()
            p.translate(rect.center())
            p.rotate(40)
            # height = rect.width()*0.8
            p.translate(-rect.center())
            p.drawText(rect, flags, text)
            # p.rotate(2)
            # p.drawRect(rect)
            p.restore()
        # self.setStyle(tickTextHeight=int(height)*5)

    def tickSpacing(self, minVal, maxVal, size):
        """
        I delcare that my desired possible spacings are every
        1s, 5s, 10s, 15s, 30s,
        1m, 5m, 10m, 15m, 30m,
        1h, 5h, 10h
        And hopefully no further than that?
        """
        superRet = super(DateAxis, self).tickSpacing(minVal, maxVal, size)
        # return ret
        # Todo set spacing to be reasonable
        dif = abs(maxVal - minVal)
        spacings = np.array([1, 5, 10, 15, 30,
                             1*60, 5*60, 10*60, 15*60, 30*60,
                             1*60*60, 5*60*60, 10*60*60])
        numTicks = (maxVal - minVal)/spacings

        # I really want to know where this comes from,
        # I just nabbed it from Luke's code
        optimalTickCount = max(2., np.log(size))

        bestValidx = np.abs(numTicks-optimalTickCount).argmin()
        desiredTicks = numTicks[bestValidx]
        if desiredTicks > 20:
            # Too many ticks to plot, would cause the render engine to break
            # Cutoff is arbitrary
            # todo: set it to a density (px/spacing) to handle it better
            return superRet
        bestVal = spacings[bestValidx]
        # todo: set better minor tick spacings
        ret =  [
            (bestVal, 0),
            (bestVal/5., 0),
            (bestVal/10., 0)
        ]
        return ret

        return super(DateAxis, self).tickSpacing(minVal, maxVal, size)


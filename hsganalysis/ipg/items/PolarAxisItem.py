from pyqtgraph import Point, AxisItem, mkBrush, mkPen
from PyQt5 import QtGui, QtCore, QtWidgets
import numpy as np
class PolarAxis(AxisItem):
    def __init__(self, orientation, pen=None, linkView=None, parent=None,
                 maxTickLength=-5, showValues=True,
                 radialBounds = None, azimuthalBounds = None):
        # force default bottom orientation since it can be "azimuthal" or "radial" here
        super().__init__("bottom", pen, linkView, parent, maxTickLength, showValues)
        if orientation not in ["radial", "azimuthal"]:
            raise RuntimeError("Polar axis must be radial or azimutal, not", orientation)

        self.pOrientation = orientation
        # self._maxValues = {
        #     "radial": maxRadial,
        #     "azimuthal": maxAngle
        # }

        if radialBounds is None:
            radialBounds = [1e-6, 40]
        if azimuthalBounds is None:
            azimuthalBounds = [-180, 180]
        self._bounds = {
            "radial": radialBounds,
            "azimuthal": azimuthalBounds
        }

        if self.pOrientation == "radial":
            self._tickSpacing = 10
        elif self.pOrientation == "azimuthal":
            self._tickSpacing = 30


        # I want to know what the full bounding rect would be, but
        # .boundingRect() is overridden to only give the view range
        self.fullBoundingRect = QtCore.QRectF(0,0,radialBounds[1]*2,radialBounds[1]*2)
        self.fullBoundingRect.moveCenter(Point(0,0))


        # At what azimuthal angle should I put the radial axis labels?
        self._radialLabelPosition = 75
        self.tickFont = QtGui.QFont("Arial", 2)
        # self.tickFont.setPointSizeF(1.2)
        # self.tickFont.setPixelSize(2)

    def paint(self, p, opt, widget):
        if self.picture is None:
            try:
                picture = QtGui.QPicture()
                painter = QtGui.QPainter(picture)
                specs = self.generateDrawSpecs(painter)
                if specs is not None:
                    self.drawPicture(painter, opt, widget, *specs)
            finally:
                painter.end()
            self.picture = picture
        #p.setRenderHint(p.Antialiasing, False)   ## Sometimes we get a segfault here ???
        #p.setRenderHint(p.TextAntialiasing, True)
        self.picture.play(p)

    def drawPicture(self, p, opt, widget, axisSpec, tickSpecs, textSpecs):
        p.setRenderHint(p.Antialiasing, True)
        p.setRenderHint(p.TextAntialiasing, True)
        orig = QtCore.QPointF(0, 0)
        if self.pOrientation is "radial":
            for pen, p1, p2 in tickSpecs:
                p.setPen(pen)
                p.drawEllipse(orig, p1.y(), p1.y())
        else:
            for pen, p1, p2 in tickSpecs:
                p.setPen(pen)
                p.drawLine(p1, p2)


        # # item = QtWidgets.QGraphicsSimpleTextItem()
        # item = QtWidgets.QGraphicsTextItem()
        # # form = QtGui.QTextCharFormat()
        # # form.setTextOutline()
        # # item.setFont(self.tickFont)
        # p.setPen(self.pen())
        # for rect, flags, text in textSpecs:
        #     # item.setText(text)
        #     item.setPlainText(text)
        #     # item.setPos(rect.topLeft())
        #     p.translate(rect.x(), rect.y())
        #     p.scale(0.2, 0.2)
        #     item.paint(p, opt, widget)
        #     p.scale(5,5)
        #     p.translate(-rect.x(), -rect.y())
        #     # p.drawRect(rect)



        if self.tickFont is not None:
            p.setFont(self.tickFont)
        # p.setOpacity(0.25)
        for rect, flags, text in textSpecs:
            p.setBrush(mkBrush("#FFFFFFAA"))
            p.setPen(mkPen(None))
            p.drawRect(rect)
            p.setPen(self.pen())
            p.drawText(rect, flags, text)
        p.setOpacity(1)

    def tickSpacing(self, minVal, maxVal, size):
        # if self.pOrientation is "radial":
        #     return [
        #         (10, 0),
        #         (10, 0)
        #     ]
        return [
            (self._tickSpacing, 0),
            (self._tickSpacing, 0)
        ]

    def tickValues(self, minVal, maxVal, size):
        tickLevels = self.tickSpacing(minVal, maxVal, size)

        maxVal = max(np.abs([minVal, maxVal]))
        maxVal = min(maxVal, self._bounds[self.pOrientation][1])

        ticks = []

        for spacing, offset in tickLevels:
            ticks.append([spacing, np.arange(minVal//spacing, maxVal//spacing + 1)*spacing])

        # make sure it's in there?
        # if maxVal != np.inf and maxVal not in ticks[0][1]:
        #     ticks[0][1].append(maxVal)
        return ticks

    def boundingRect(self):
        # return the full view
        # This works for the radial grid, but I dunno about azimuthal
        linkedView = self.linkedView()
        # not sure hwy he has the default with geometry
        # return self.mapRectFromParent(self.geometry()) | linkedView.mapRectToItem(self,
        #                                                                       linkedView.boundingRect())
        return self.mapRectFromParent(linkedView.mapRectToItem(self, linkedView.boundingRect()))

    def getAzimuthalLines(self, rng, ang):
        """
        Way too over-the-top rendering of lines to make them start and end at the bounds
        of the view rect
        :param rng:
        :param ang:
        :return:
        """
        c = np.cos(ang * np.pi / 180)
        # Don't forget to flip the axis
        s = np.sin(-ang * np.pi / 180)
        # get the center of the viewrect and the
        # _radial_ diagonal of the rectangle
        cx, cy = rng.center().x(), rng.center().y()
        diagonal = Point(
            (-1) ** (cx < 0) * rng.width()  / 2,
            (-1) ** (cy < 0) * rng.height() / 2)

        minX, minY = rng.center() - diagonal
        maxX, maxY = rng.center() + diagonal

        # Stop it at the origin (so the user can specify 0->180 instead?)
        if rng.contains(0,0):
            r1 = 0
        else:
            r1 = np.sqrt(minX**2 + minY**2)

        # Check the bounds of the rect, the minimum needs to be one of the edges
        # if your viewRect is orthogonal to nPi/2
        r1 = min((r1, abs(rng.left()), abs(rng.right()), abs(rng.top()), abs(rng.bottom())))
        r2 = np.sqrt(maxX**2 + maxY**2)
        r2 = min(r2, self._bounds["radial"][1])

        p1 = [r1 * c, r1 * s]
        p2 = [r2 * c, r2 * s]

        return p1, p2

    def generateDrawSpecs(self, p):
        """
        Calls tickValues() and tickStrings() to determine where and how ticks should
        be drawn, then generates from this a set of drawing commands to be
        interpreted by drawPicture().
        """

        # bounds = self.boundingRect()

        self.fullBoundingRect = QtCore.QRectF(0,0,
                                      self._bounds["radial"][1]*2,
                                      self._bounds["radial"][1]*2)
        self.fullBoundingRect.moveCenter(Point(0,0))


        bounds = self.mapRectFromParent(self.geometry())

        linkedView = self.linkedView()
        tickBounds = linkedView.mapRectToItem(self, linkedView.boundingRect())

        span = (bounds.bottomLeft(), bounds.topRight())
        tickStart = 0
        tickStop = max(bounds.width(), bounds.height())
        tickDir = 0
        axis = 0
        # print tickStart, tickStop, span

        ## determine size of this item in pixels
        points = list(map(self.mapToDevice, span))
        if None in points:
            return
        lengthInPixels = Point(points[1] - points[0]).length()
        if lengthInPixels == 0:
            return

        # Determine major / minor / subminor axis ticks
        if self._tickLevels is None:
            rng = self.linkedView().viewRange()
            if self.pOrientation is "radial":
                tickLevels = self.tickValues(self.range[0], np.sqrt(max(np.abs(rng[0]))**2 + max(np.abs(rng[1]))**2), lengthInPixels)
            else:
                # todo undo hardcoded angles
                tickLevels = self.tickValues(self._bounds["azimuthal"][0],
                                             self._bounds["azimuthal"][1],
                                             lengthInPixels)
            tickStrings = None
        else:
            ## parse self.tickLevels into the formats returned by tickLevels() and tickStrings()
            tickLevels = []
            tickStrings = []
            for level in self._tickLevels:
                values = []
                strings = []
                tickLevels.append((None, values))
                tickStrings.append(strings)
                for val, strn in level:
                    values.append(val)
                    strings.append(strn)

        # Not needed for radial plots
        xScale = 1
        offset = 0

        if self.pOrientation is "radial":
            rng = self.linkedView().viewRange()
            xRange = [0] + [abs(x * xScale - offset) for x in self.range]
            xMin = self._bounds["radial"][0]
            xMax = np.sqrt(max(np.abs(rng[0]))**2 + max(np.abs(rng[1]))**2)
            xMax = max(xMax, self._bounds["radial"][0])
        else:
            xRange = [x * xScale - offset for x in self.range]
            xMin = self._bounds["azimuthal"][0]
            xMax = self._bounds["azimuthal"][1]




        tickPositions = []  # remembers positions of previously drawn ticks

        ## compute coordinates to draw ticks
        ## draw three different intervals, long ticks first
        tickSpecs = []
        for i in range(len(tickLevels)):
            tickPositions.append([])
            ticks = tickLevels[i][1]

            ## length of tick
            tickLength = self.style['tickLength'] / ((i * 0.5) + 1.0)

            lineAlpha = 255 / (i + 3)
            if self.grid is not False:
                lineAlpha *= self.grid / 255. * np.clip(
                    (0.05 * lengthInPixels / (len(ticks) + 1)), 0., 1.)

            for v in ticks:
                ## determine actual position to draw this tick
                x = (v * xScale) - offset
                if x < xMin or x >= xMax:  ## last check to make sure no out-of-bounds ticks are drawn
                    tickPositions[i].append(None)
                    continue

                tickPositions[i].append(x)

                if self.pOrientation is "radial":
                    p1 = [x, x]
                    p2 = [x, x]
                    p1[axis] = tickStart
                    p2[axis] = tickStop
                    if self.grid is False:
                        p2[axis] += tickLength * tickDir
                else:
                    p1, p2 = self.getAzimuthalLines(self.boundingRect(), x)


                tickPen = self.pen()
                color = tickPen.color()
                color.setAlpha(lineAlpha)
                tickPen.setColor(color)
                tickSpecs.append((tickPen, Point(p1), Point(p2)))

        if self.style['stopAxisAtTick'][0] is True:
            stop = max(span[0].y(), min(map(min, tickPositions)))
            if axis == 0:
                span[0].setY(stop)
            else:
                span[0].setX(stop)
        if self.style['stopAxisAtTick'][1] is True:
            stop = min(span[1].y(), max(map(max, tickPositions)))
            if axis == 0:
                span[1].setY(stop)
            else:
                span[1].setX(stop)
        axisSpec = (self.pen(), span[0], span[1])

        textOffset = self.style['tickTextOffset'][
            axis]  ## spacing between axis and text
        # if self.style['autoExpandTextSpace'] is True:
        # textWidth = self.textWidth
        # textHeight = self.textHeight
        # else:
        # textWidth = self.style['tickTextWidth'] ## space allocated for horizontal text
        # textHeight = self.style['tickTextHeight'] ## space allocated for horizontal text

        textSize2 = 0
        textRects = []
        textSpecs = []  ## list of draw

        # If values are hidden, return early
        if not self.style['showValues']:
            return (axisSpec, tickSpecs, textSpecs)

        # For appropriate sizing of boxes
        if self.tickFont is None: return (axisSpec, tickSpecs, textSpecs)
        p.setFont(self.tickFont)
        for i in range(min(len(tickLevels), self.style['maxTextLevel'] + 1)):
            ## Get the list of strings to display for this level
            if tickStrings is None:
                spacing, values = tickLevels[i]
                strings = self.tickStrings(values, self.autoSIPrefixScale * self.scale,
                                           spacing)
            else:
                strings = tickStrings[i]

            if len(strings) == 0:
                continue

            ## ignore strings belonging to ticks that were previously ignored
            for j in range(len(strings)):
                if tickPositions[i][j] is None:
                    strings[j] = None

            ## Measure density of text; decide whether to draw this level
            rects = []
            for s in strings:
                if s is None:
                    rects.append(None)
                else:
                    br = p.boundingRect(QtCore.QRectF(0, 0, 100, 100),
                                        QtCore.Qt.AlignCenter, s)
                    ## boundingRect is usually just a bit too large
                    ## (but this probably depends on per-font metrics?)
                    br.setHeight(br.height() * 0.8)

                    rects.append(br)
                    textRects.append(rects[-1])

            if len(textRects) > 0:
                ## measure all text, make sure there's enough room
                if axis == 0:
                    textSize = np.sum([r.height() for r in textRects])
                    textSize2 = np.max([r.width() for r in textRects])
                else:
                    textSize = np.sum([r.width() for r in textRects])
                    textSize2 = np.max([r.height() for r in textRects])
            else:
                textSize = 0
                textSize2 = 0

            if i > 0:  ## always draw top level
                ## If the strings are too crowded, stop drawing text now.
                ## We use three different crowding limits based on the number
                ## of texts drawn so far.
                textFillRatio = float(textSize) / lengthInPixels
                finished = False
                for nTexts, limit in self.style['textFillLimits']:
                    if len(textSpecs) >= nTexts and textFillRatio >= limit:
                        finished = True
                        break
                if finished:
                    break

            # spacing, values = tickLevels[best]
            # strings = self.tickStrings(values, self.scale, spacing)
            # Determine exactly where tick text should be drawn
            for j in range(len(strings)):
                vstr = strings[j]
                if vstr is None:  ## this tick was ignored because it is out of bounds
                    continue
                x = tickPositions[i][j]
                # textRect = p.boundingRect(QtCore.QRectF(0, 0, 100, 100), QtCore.Qt.AlignCenter, vstr)
                textRect = rects[j]
                height = textRect.height()
                width = textRect.width()
                # self.textHeight = height
                offset = max(0, self.style['tickLength']) + textOffset
                textFlags = QtCore.Qt.TextDontClip | QtCore.Qt.AlignCenter | QtCore.Qt.AlignVCenter




                if self.pOrientation is "azimuthal":

                    # Originally wanted to put the text on the outside of
                    # the viewbox, but that was breaking things, especially
                    # because it wasn't handling intialization correctly, when the
                    # viewbox is first being setup.
                    #pp = tickSpecs[j][-1] # Take the outside coordinate
                    pp = self._bounds["radial"][1] if np.isfinite(self._bounds["radial"][1]) else tickSpecs[j][-1]
                    pp = Point(pp*np.cos(x*3.14159/180), -pp*np.sin(x*3.14159/180))

                    rect = QtCore.QRectF(0,0, width, height)
                    rect.moveCenter(Point(0,0))
                    xm = x%360
                    if xm == 0:
                        rect.moveLeft(pp.x())
                        # print("00pp", pp)
                    elif 0<xm<90:
                        rect.moveBottomLeft(pp)
                    elif xm == 90:
                        rect.moveBottom(pp.y())
                        # print("90pp", pp)
                    elif 90<xm<180:
                        rect.moveBottomRight(pp)
                    elif xm == 180:
                        rect.moveRight(pp.x())
                    elif 180<xm<270:
                        rect.moveTopRight(pp)
                    elif xm == 270:
                        rect.moveTop(pp.y())
                    elif 270<xm<360:
                        rect.moveTopLeft(pp)
                elif self.pOrientation is "radial":
                    rect = QtCore.QRectF(0, 0, width, height)
                    # rect.moveCenter(Point(0, 0))
                    rect.moveBottomLeft(Point(
                        x*np.cos(self._radialLabelPosition * np.pi / 180),
                    # Don't forget to flip the axis
                        x*np.sin(-self._radialLabelPosition * np.pi / 180)
                    ))
                self.fullBoundingRect |= rect

                # p.setPen(self.pen())
                # p.drawText(rect, textFlags, vstr)
                textSpecs.append((rect, textFlags, vstr))

        ## update max text size if needed.
        self._updateMaxTextSize(textSize2)

        return (axisSpec, tickSpecs, textSpecs)

    def mouseDragEvent(self, event):
        event.ignore()

    def mouseClickEvent(self, event):
        event.ignore()


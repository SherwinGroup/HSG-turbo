import numpy as np
import pyqtgraph
from pyqtgraph import mkPen, mkBrush
from PyQt5 import QtCore, QtGui
from ..packageSettings import config_options
from .LegendSettings_ui import Ui_LegendSettingsDialog
"""

pyqtgraph's legend item has a hard-coded slightly transparent,
grey background. Modify the class methods to allow you to change
that.

"""

def createLegendList(descriptor):
    """
    Given a string with item/name descriptions and
    create pairs of ItemSamples/name to be populated into the legend
    """
    pass

def createDescriptorFromLegend(legend):
    """
    Given a legend item, create the string descriptor
    which would create it.
    """
    pass


class LegendSettingsDialog(QtGui.QDialog):
    def __init__(self, *args, **kwargs):
        self.legendItem = kwargs.pop("legendItem", pyqtgraph.LegendItem())
        super(LegendSettingsDialog, self).__init__(*args, **kwargs)
        self.initialSettings = {
            "pen": self.legendItem.backgroundPen.color(),
            "brush": self.legendItem.backgroundBrush.color(),
            "font": self.legendItem.items[-1][1].opts.get("size", config_options["legendFontSize"]).split('p')[0]
        }

        self.ui = Ui_LegendSettingsDialog()
        self.ui.setupUi(self)
        self.initUI()

        self.ui.bBGColor.sigColorChanging.connect(self.updateLegend)
        self.ui.bBorderColor.sigColorChanging.connect(self.updateLegend)
        self.ui.sbFontSize.setOpts(int=True, step=1, bounds=(1,100))
        self.ui.sbFontSize.sigValueChanging.connect(self.updateLegend)

    def initUI(self):
        self.ui.bBGColor.setColor(self.initialSettings["brush"])
        self.ui.bBorderColor.setColor(self.initialSettings["pen"])
        self.ui.sbFontSize.setValue(int(self.initialSettings["font"]))

    def updateLegend(self):
        self.legendItem.setBackgroundPen(mkPen(self.ui.bBorderColor.color()))
        self.legendItem.setBackgroundBrush(mkBrush(self.ui.bBGColor.color()))
        for sample, label in self.legendItem.items:
            # label.setAttr("size", "{}pt".format(self.ui.sbFontSize.value()))
            label.setText(label.text, size="{}pt".format(self.ui.sbFontSize.value()))
            sample.setScale(self.ui.sbFontSize.value()/10.)
            label.setGeometry(label.item.boundingRect())
        self.legendItem.layout.setColumnMinimumWidth(0, 10*self.ui.sbFontSize.value()/10.)
        self.legendItem.update()
        self.legendItem.updateSize()

    @staticmethod
    def makeSettings(legendItem):
        dialog = LegendSettingsDialog(legendItem=legendItem)
        ok = dialog.exec_()
        if not ok:
            dialog.initUI()
            dialog.updateLegend()

oldinit = pyqtgraph.LegendItem.__init__
def __init__(self, size=None, offset=None):
    oldinit(self, size, offset)
    self.backgroundBrush = mkBrush(config_options["legendBackground"])
    self.backgroundPen = mkPen(config_options["legendBorder"])
    self.layout.setColumnStretchFactor(0, 5)
    self.layout.setColumnStretchFactor(1, 5)

oldparent = pyqtgraph.LegendItem.setParentItem
def setParentItem(self, p):
    # print "set parent item in legenditemfix", p
    # print "current parent", self.parent()
    ret = oldparent(self, p)
    self.setParent(p)
    # print "new parent", self.parent()
    self.scene().sigMouseClicked.connect(self.mouseClickedEvent)
    return ret

def paint(self, p, *args):
    p.setPen(self.backgroundPen)
    p.setBrush(self.backgroundBrush)
    p.drawRect(self.boundingRect())
    # for samp, label in self.items:
    #     p.translate(label.geometry().x(), label.geometry().y())
    #     p.drawRect(label.itemRect())
    #     p.translate(-label.geometry().x(), -label.geometry().y())

def setBackgroundPen(self, p):
    self.backgroundPen = mkPen(p)

def setBackgroundBrush(self, b):
    self.backgroundBrush = mkBrush(b)

def mouseClickedEvent(self, ev):
    if self.contains(self.mapFromScene(ev.scenePos())) and ev.double():
        self.openSettings()

def openSettings(self):
    import pyqtgraph.console as pgc
    self.a = pgc.ConsoleWidget(namespace={"self":self})
    self.a.show()
    LegendSettingsDialog.makeSettings(self)

def updateSize(self):
    if self.size is not None:
        return

    height = 0
    width = 0
    symWidth = 0
    labelWidth = 0
    #print("-------")
    for sample, label in self.items:
        height += max(sample.height(), label.height()) + 3
        width = max(width, sample.width()+label.width())
        symWidth = max(symWidth, sample.boundingRect().width()*sample.scale())
        labelWidth = max(labelWidth, label.boundingRect().width())
        #print(width, height)
    #print width, height
    self.setGeometry(0, 0, width+25, height)
    self.layout.setColumnMaximumWidth(0,symWidth)
    self.layout.setColumnMinimumWidth(0,symWidth)
    self.layout.setColumnMaximumWidth(1, labelWidth)
    self.layout.setColumnMinimumWidth(1, labelWidth)

def setFont(self, **kwargs):
    for sample, label in self.items:
        label.setText(label.text, **kwargs)
        # sample.setScale(self.ui.sbFontSize.value() / 10.)
        label.setGeometry(label.item.boundingRect())

# Drop-in support for matplotlib legends.
def draggable(self, bl=False):
    pass

pyqtgraph.ViewBox
pyqtgraph.GraphicsScene
pyqtgraph.LegendItem.__init__ = __init__
pyqtgraph.LegendItem.paint = paint
pyqtgraph.LegendItem.setParentItem = setParentItem
pyqtgraph.LegendItem.updateSize = updateSize
pyqtgraph.LegendItem.setBackgroundPen = setBackgroundPen
pyqtgraph.LegendItem.setBackgroundBrush = setBackgroundBrush

pyqtgraph.LegendItem.setFont = setFont

pyqtgraph.LegendItem.openSettings = openSettings
pyqtgraph.LegendItem.draggable = draggable
pyqtgraph.LegendItem.mouseClickedEvent = mouseClickedEvent

# pyqtgraph.LegendItem.__init__ = LegendItemFix.__init__
# pyqtgraph.LegendItem.paint = LegendItemFix.paint
# pyqtgraph.LegendItem.setParentItem = LegendItemFix.setParentItem
# pyqtgraph.LegendItem.updateSize = LegendItemFix.updateSize
# pyqtgraph.LegendItem.setBackgroundPen = LegendItemFix.setBackgroundPen
# pyqtgraph.LegendItem.setBackgroundBrush = LegendItemFix.setBackgroundBrush
#
# pyqtgraph.LegendItem.openSettings = LegendItemFix.openSettings
# pyqtgraph.LegendItem.mouseClickedEvent = LegendItemFix.mouseClickedEvent





# pyqtgraph.LegendItem = LegendItemFix

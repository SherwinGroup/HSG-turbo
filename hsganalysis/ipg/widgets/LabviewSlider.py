# -*- coding: utf-8 -*-
"""
Created on Thu Sep 25 15:41:33 2014

@author: dvalovcin
"""

import numpy as np
from PyQt5 import QtWidgets, QtCore
import pyqtgraph as pg


class LabviewSlider(QtWidgets.QWidget):
    sigValueChanging = QtCore.pyqtSignal(object)
    def __init__(self, parent=None, **kwargs):
        super(LabviewSlider, self).__init__(parent)
        self._spinbox = pg.SpinBox(self)

        self._lMin = QtWidgets.QLineEdit(self)
        self._lMax = QtWidgets.QLineEdit(self)
        for wid in [self._lMin, self._lMax]:
            wid.setStyleSheet(" QLineEdit { "
                                " background-color: rgba(0, 0, 0, 0); }")
            wid.setAttribute(QtCore.Qt.WA_MacShowFocusRect, 0)
            # wid.setStyleSheet(" QLineEdit { "
            #             " color:rgb(0, 0, 0); "
            #             " background-color: rgb(255, 255, 255); } "
            #
            #             " QLineEdit:focus { "
            #             " border:4px outset; "
            #             " border-radius: 8px; "
            #             " border-color: rgb(41, 237, 215); "
            #             " color:rgb(0, 0, 0); "
            #             " background-color: rgb(150, 150, 150); } "
            #        )
            wid.editingFinished.connect(self._parseBoundsChange)
            wid.setFrame(False)

        self._lMax.setAlignment(QtCore.Qt.AlignRight)
        self._lMin.setAlignment(QtCore.Qt.AlignLeft)

        sliderLayout = QtWidgets.QGridLayout()

        orientation=kwargs.get("orientation", 0)
        if orientation:
            self._slider = QtWidgets.QSlider(self, orientation=QtCore.Qt.Vertical)
            layout = QtWidgets.QVBoxLayout(self)
            layout.addWidget(self._spinbox)
            layout.addWidget(self._slider)
            layout.setStretch(1, 10)
            layout.setStretch(0, 1)
        else:
            self._slider = QtWidgets.QSlider(self, orientation=QtCore.Qt.Horizontal)
            sliderLayout.addWidget(self._slider, 0, 0, 1, 3)
            sliderLayout.addWidget(self._lMin, 1, 0, QtCore.Qt.AlignLeft)
            sliderLayout.addWidget(self._lMax, 1, 2, QtCore.Qt.AlignRight)
            # sliderLayout.setColumnStretch(0, 1)
            # sliderLayout.setColumnStretch(1, 10)
            # sliderLayout.setColumnStretch(2, 1)
            layout = QtWidgets.QHBoxLayout(self)
            layout.addLayout(sliderLayout)
            layout.addWidget(self._spinbox)
            layout.setStretch(0, 5)
            layout.setStretch(1, 1)

        self.setLayout(layout)

        self._slider.valueChanged.connect(self._updateSpinboxValue)
        self._spinbox.sigValueChanging.connect(self._updateSliderValue)
        self._slider.setTickPosition(QtWidgets.QSlider.TicksBothSides)

        self._slider.setRange(0, 1000)
        self._slider.setMinimumHeight(20) #eyeballed it, default was way too large

        self.opts = {"orientation": 0,
                     "range": (1, 2),
                     "step": 0.1,
                     "value": 0}

        self.opts.update(kwargs)
        self.updateSettings()




    def _updateSpinboxValue(self, val=None):
        self._spinbox.blockSignals(True)
        minVal, maxVal = self.range()
        step = self.step()
        newVal = minVal + val*step
        self._spinbox.setValue(newVal)
        self._spinbox.blockSignals(False)
        self.sigValueChanging.emit(self.value())

    def _updateSliderValue(self, sb=None):
        self._slider.blockSignals(True)

        val = self._spinbox.value()
        minVal, maxVal = self.range()
        step = self.step()
        newVal = np.round((val-minVal)/step)
        self._slider.setValue(newVal)
        self._slider.blockSignals(False)
        self.sigValueChanging.emit(self.value())

    def _updateSliderRange(self):
        """
        QSlider only allows ints. This function
        will map the min/max and step of the desired
        range onto ints from 0 to N, where N = (max-min)/step
        """
        minVal, maxVal = self.range()
        step = self.step()
        self._slider.setMinimum(0)
        self._slider.setMaximum(round((maxVal-minVal)/step))
        self._slider.setTickInterval(self._slider.maximum()/10)

        self._updateSliderValue()

    def _parseBoundsChange(self):
        mn = float(self._lMin.text())
        mx = float(self._lMax.text())
        self._lMin.clearFocus()
        self._lMax.clearFocus()
        if mn>=mx:
            print("failed update")
            self._lMin.blockSignals(True)
            self._lMax.blockSignals(True)
            self._lMin.setText(str(self.opts["range"][0]))
            self._lMax.setText(str(self.opts["range"][1]))
            self._lMin.blockSignals(False)
            self._lMax.blockSignals(False)
        else:
            self.setRange((mn, mx))


    def updateSettings(self, **kwargs):
        self.opts.update(kwargs)

        self.setRange(self.opts["range"])
        self.setValue(self.opts["value"])
        self.setStep(self.opts["step"])

    def setRange(self, bounds, step=None):
        self._spinbox.setRange(*bounds)
        self.opts["range"] = tuple(bounds)

        self._lMax.blockSignals(True)
        self._lMin.blockSignals(True)
        self._lMin.setText("{}".format(bounds[0]))
        self._lMax.setText("{}".format(bounds[1]))
        self._lMax.blockSignals(False)
        self._lMin.blockSignals(False)
        self._updateSliderRange()
        if step is None:
            self.setStep((bounds[1]-bounds[0])/100)
        else:
            self.setStep(step)

    def range(self):
        return self.opts["range"]

    def setStep(self, value):
        self._spinbox.setSingleStep(value)
        self.opts["step"] = value
        self._updateSliderRange()

    def step(self):
        return self.opts["step"]

    def setValue(self, value):
        self._spinbox.setValue(value)
        self.opts["value"] = self._spinbox.value()

    def value(self):
        return self._spinbox.value()



if __name__ == '__main__':
    import sys
    app = QtWidgets.QApplication([])
    wid = LabviewSlider()
    wid.show()
    sys.exit(app.exec_())








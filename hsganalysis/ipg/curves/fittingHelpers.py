
import pyqtgraph as pg
import numpy as np
import sys
from PyQt5 import QtCore, QtGui
from .clickablePlotSettings_ui import Ui_LineSettingsDialog
from .PlotDataErrorItem import *
from ..packageSettings import config_options

from scipy.fftpack import rfft, rfftfreq



def group_data(data, cutoff = 7):
    """
    Given a list of numbers, group them up based on where
    they're continuous
    """
    groupedData = []
    curData = []
    for d in list(data):
        if curData:
            if np.abs(d-curData[-1])>cutoff:
                groupedData.append(curData)
                curData = [d]
            else:
                curData.append(d)
        else:
            curData.append(d)
    groupedData.append(curData)
    return groupedData

def getFrequency(x, y):
    """
    try to estimate the frequency of a data set based
    on the extrema of the data.
    """


    ymx, ymn = np.max(y), np.min(y)
    ymaxIdx = np.where(y>(ymx - 0.05*(ymx-ymn)))[0]
    yminIdx = np.where(y<(ymn + 0.05*(ymx-ymn)))[0]

    ymaxIdxGroups = group_data(ymaxIdx)
    yminIdxGroups = group_data(yminIdx)


    yMaxLocs = []
    for lst in ymaxIdxGroups:
        yMaxLocs.append(np.mean(x[lst]))

    yMinLocs = []
    for lst in yminIdxGroups:
        yMinLocs.append(np.mean(x[lst]))

    yMaxLocs.extend(yMinLocs)
    yExtremaLoc = np.array(yMaxLocs)
    yExtremaLoc = yExtremaLoc[yExtremaLoc.argsort()]
    T = np.polyfit(np.arange(yExtremaLoc.shape[0]), yExtremaLoc, deg=1)[0]/np.pi
    return 1./T


class FuncHelper(object):
    def __call__(self, *args, **kwargs):
        pass
    def __str__(self):
        return self.__repr__()

    def __repr__(self):
        pass

    def guessParameters(self, x, y):
        x = x.copy()
        y = y.copy()


class sine(FuncHelper):
    def __call__(self, x, *p):
        A, mu, f, offset = p
        return A * np.sin(f*(x-mu)) + offset
    def __repr__(self):
        return "A={:.5g}, mu={:.5g},\n  f={:.5g}, y0={:.5g}"
    def guessParameters(self, x, y):
        super(sine, self).guessParameters(x, y)
        A = np.max(y) - np.min(y)
        offset = np.min(y) + 1
        f = getFrequency(x, y)
        mu = 0
        return [A, mu, f, offset]


class cosine(FuncHelper):
    def __call__(self, x, *p):
        A, mu, f, offset = p
        return A * np.cos(f*(x-mu)) + offset
    def __repr__(self):
        return "A={:.5g}, mu={:.5g},\n  f={:.5g}, y0={:.5g}"
    def guessParameters(self, x, y):
        super(cosine, self).guessParameters(x, y)
        A = np.max(y) - np.min(y)
        offset = np.min(y) + 1
        f = getFrequency(x, y)
        mu = 0
        return [A, mu, f, offset]

class exp(FuncHelper):
    def __call__(self, x, *p):
        A, l, offset = p
        return A * np.exp(l*x) + offset
    def __repr__(self):
        return "A={:.5g},\n  lamb={:.5g}, y0={:.5g}"
    def guessParameters(self, x, y):
        super(exp, self).guessParameters(x, y)
        y = y.copy()
        x = x.copy()

        offset = np.min(y)
        y -= offset
        x1 = x[0]
        x2 = x[-1]
        y1 = y[0]
        y2 = y[-1]

        yi = y[y.size // 2]
        xi = x[y.size // 2]

        d2 = (y2-yi)/(x2-xi) - (yi-y1)/(xi-x1)
        d2 = np.sign(d2)

        l = np.log(y1/y2)/(x1-x2)
        A = np.exp(np.log(y2)-l*x2)

        l *= d2
        A *= d2 * np.sign(y2-y1)


        return [A, l, offset]

class gauss(FuncHelper):
    def __call__(self, x, *p):
        mu, A, sigma, y0 = p
        return (A / sigma) * np.exp(-(x - mu)**2 / (2. * sigma**2)) + y0
    def __repr__(self):
        return "mu={:.5g}, A={:.5g},\n  sig={:.5g}, y0={:.5g}"
    def guessParameters(self, x, y):
        super(gauss, self).guessParameters(x, y)
        x = x.copy()
        y = y.copy()
        # Get height and offset
        bottom = np.min(y)

        height = np.max(y)
        y -= bottom

        # find FWHM
        hm = height/2.
        argHeight = np.argmax(y)

        argBelowHalfHeight = np.where(y<hm)[0] - argHeight
        # Find out where the sign changes
        argProd = np.sign(argBelowHalfHeight[1:]*argBelowHalfHeight[:-1])

        zeroCrossing = np.argmin(argProd)

        argBelowHalfHeight+= argHeight
        argLowHM = argBelowHalfHeight[zeroCrossing]
        argHighHM = argBelowHalfHeight[zeroCrossing+1]

        FWHM = x[argHighHM] - x[argLowHM]
        sigma = FWHM / 2.35 # ~ 2sqrt(2ln(2)) to convert between them

        mu = x[argHeight]
        return [mu, height, sigma, bottom]

class lorentzian(FuncHelper):
    def __call__(self, x, *p):
        mu, A, gamma, y0 = p
        return (A / np.pi) * (gamma / ((x - mu)**2 + gamma**2)) + y0
    def __repr__(self):
        return "mu={:.5g}, A={:.5g},\n  gam={:.5g}, y0={:.5g}"
    def guessParameters(self, x, y):# Get height and offset
        super(lorentzian, self).guessParameters(x, y)
        x = x.copy()
        y = y.copy()
        bottom = np.min(y)

        height = np.max(y)
        y -= bottom

        # find FWHM
        hm = height/2.
        argHeight = np.argmax(y)

        argBelowHalfHeight = np.where(y<hm)[0] - argHeight
        # Find out where the sign changes
        argProd = np.sign(argBelowHalfHeight[1:]*argBelowHalfHeight[:-1])

        zeroCrossing = np.argmin(argProd)

        argBelowHalfHeight+= argHeight
        argLowHM = argBelowHalfHeight[zeroCrossing]
        argHighHM = argBelowHalfHeight[zeroCrossing+1]

        FWHM = x[argHighHM] - x[argLowHM]

        mu = x[argHeight]
        return [mu, height, FWHM, bottom]


class polynomial(FuncHelper):
    def __init__(self, deg=1):
        self._deg = deg
        super(polynomial, self).__init__()

    def __call__(self, x, *p):
        return np.polyval(p, x)
    def __repr__(self):
        if self._deg==1:
            return "m={:.5g}, b={:.5g}"
        if self._deg==2:
            return "a={:.5g}, b={:.5g}, c={:.5g}"
        if self._deg==3:
            return "a={:.5g}, b={:.5g},\n  c={:.5g}, d={:.5g}"
    def guessParameters(self, x, y):
        super(polynomial, self).guessParameters(x, y)
        return [1] * (self._deg + 1)

"""
def sine(x, *p):
    A, mu, f, offset = p
    return A * np.sin(f*(x-mu)) + offset

def cosine(x, *p):
    A, mu, f, offset = p
    return A * np.cos(f*(x-mu)) + offset

def exp(x, *p):
    A, mu, l, offset = p
    return A * np.exp(l*(x-mu)) + offset

def gauss(x, *p):
    mu, A, sigma, y0 = p
    return (A / sigma) * np.exp(-(x - mu)**2 / (2. * sigma**2)) + y0

def lorentzian(x, *p):
    mu, A, gamma, y0 = p
    return (A / np.pi) * (gamma / ((x - mu)**2 + gamma**2)) + y0

def polynomial(x, *p):
    return np.polyval(p, x)
"""
def guessP0(func, x, y):
    if func=="Sine" or func=="Cosine":
        A = np.max(y) - np.min(y)
        offset = np.min(y) + 1
        f = getFrequency(x, y)
        mu = 0
        return [A, mu, f, offset]
    elif func == "Exp":
        return [1, 1, 1, 1]
    elif func == "Gaussian" or func == "Lorentzian":
        # Get height and offset
        bottom = np.min(y)

        height = np.max(y)
        y -= bottom

        # find FWHM
        hm = height/2.
        argHeight = np.argmax(y)

        argBelowHalfHeight = np.where(y<hm)[0] - argHeight
        # Find out where the sign changes
        argProd = np.sign(argBelowHalfHeight[1:]*argBelowHalfHeight[:-1])

        zeroCrossing = np.argmin(argProd)

        argBelowHalfHeight+= argHeight
        argLowHM = argBelowHalfHeight[zeroCrossing]
        argHighHM = argBelowHalfHeight[zeroCrossing+1]

        FWHM = x[argHighHM] - x[argLowHM]
        sigma = FWHM / 2.35 # ~ 2sqrt(2ln(2)) to convert between them

        mu = x[argHeight]
        if func == "Gaussian":
            return [mu, height, sigma, bottom]
        return [mu, height, FWHM, bottom]


    elif func == "Linear":
        return [1, 1]
    elif func == "Quadratic":
        return [1, 1, 1]
    elif func == "Cubic":
        return [1, 1, 1, 1]

def getPStringFormat(func, *p):
    """
    Will generate a pretty string for use by the
    plot widget to display results of fit
    """
    if func=="Sine" or func=="Cosine":
        return "A={:.5g}, mu={:.5g},\n  f={:.5g}, y0={:.5g}"
    if func=="Exp":
        return "A={:.5g}, mu={:.5g},\n  lamb={:.5g}, y0={:.5g}"
    if func=="Gaussian":
        return "mu={:.5g}, A={:.5g},\n  sig={:.5g}, y0={:.5g}"
    if func=="Lorentzian":
        return "mu={:.5g}, A={:.5g},\n  gam={:.5g}, y0={:.5g}"
    if func=="Linear":
        return "m={:.5g}, b={:.5g}"
    if func=="Quadratic":
        return "a={:.5g}, b={:.5g}, c={:.5g}"
    if func=="Cubic":
        return "a={:.5g}, b={:.5g},\n  c={:.5g}, d={:.5g}"


getFuncs = {"Sine": sine(),
            "Cosine": cosine(),
            "Exp": exp(),
            "Gaussian": gauss(),
            "Lorentzian": lorentzian(),
            "Linear": polynomial(1),
            "Quadratic": polynomial(2),
            "Cubic": polynomial(3)}
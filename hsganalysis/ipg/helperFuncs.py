import numpy as np
import pyqtgraph as pg
import re
from pyqtgraph.graphicsItems.ScatterPlotItem import Symbols
from .packageSettings import *

def parsePenFormatString(st):
    settingDict = {}
    colPat = 'b|g|r|c|k|m'
    match = re.search(colPat, st)
    if match:
        settingDict["color"] = match.group()

    styPat = '(?:_)|(?:\-\-)|(?:\-\.\.)|(?:\-\.)|(?:\-)|(?:\.)'
    match = re.search(styPat, st)
    if match:
        settingDict["style"] = match.group()
    else:
        settingDict["style"] = "-"

    symPat = 'o|s|t|d|\+|x'
    match = re.search(symPat, st)
    if match:
        settingDict["symbol"] = match.group()
    return settingDict


def plotArgsParser(*args):
    """
    Take the args passed to a pg.plot command
    to determine if things what things look like
    use regex to compare to possible

    Need to remove this string from the args list as well,
    because a lot of pg breaks if len(args) != 1,2
    """
    settingDict = {}
    args = list(args)
    for idx, arg in enumerate(args):
        if not isinstance(arg, str):
            continue
        # colPat = 'b|g|r|c|k|m'
        # match = re.search(colPat, arg)
        # if match:
        #     settingDict["color"] = match.group()
        #
        # styPat = '(?:_)|(?:\-\-)|(?:\-\.\.)|(?:\-\.)|(?:\-)|(?:\.)'
        # match = re.search(styPat, arg)
        # if match:
        #     settingDict["style"] = match.group()
        #
        # symPat = 'o|s|t|d|\+|x'
        # match = re.search(symPat, arg)
        # if match:
        #     settingDict["symbol"] = match.group()
        pens = parsePenFormatString(arg)
        settingDict.update(pens)
        # args.remove(arg)
        # remove it from the args list

        del args[idx]
        break
    return tuple(args), settingDict

def getPlotPens(*args, **kwargs):

    kwargs = kwargs.copy()
    pw = kwargs.pop("pw", None)
    if 'marker' in kwargs:
        kwargs["symbol"] = kwargs.pop('marker')
    pen = kwargs.get('pen', None)
    if pen is None:
        pen = pg.mkPen()

        if "fmt" in kwargs and kwargs["fmt"] is not None:
            kwargs.update(parsePenFormatString(kwargs.pop("fmt")))


        newArgs, newKwargs = plotArgsParser(*args)
        kwargs.update(newKwargs)
        args = newArgs

        if isinstance(config_options["standardColors"], int):
            numCols = config_options["standardColors"]
        else:
            numCols = len(config_options["standardColors"])

        color = kwargs.get("color", None)
        if color is None:
            if pw is None:
                color = 'b'
            else:
                if isinstance(config_options["standardColors"], int):
                    idx = len(pw.plotItem.curves) % numCols
                    color = pg.intColor(idx, config_options["standardColors"])
                else:
                    idx = len(pw.plotItem.curves) % numCols
                    color = config_options["standardColors"][idx]

        style = kwargs.pop('style', None)
        if style is None:
            if pw is None:
                style = config_options["linestyleChars"].index('-')
            else:
                idx = (len(pw.plotItem.curves) // numCols) % config_options["standardLineshapes"] + 1
                style = idx
        else:
            # pass int and use the qt pen value for it
            # otherwise, parse it.
            if not isinstance(style, int):
                style = config_options["linestyleChars"].index(style)

        if 'symbol' in kwargs and 'symbolPen' not in kwargs:
            kwargs['symbolPen'] = pg.mkPen(color=color)
        if 'symbol' in kwargs and 'symbolBrush' not in kwargs:
            kwargs['symbolBrush'] = pg.mkBrush(color=color)

        if "symbol" not in kwargs:
            # If you don't specify a symbol, default the colors to
            # the line values to match if you turn them on.
            kwargs['symbolPen'] = pg.mkPen(color=color)
            kwargs['symbolBrush'] = pg.mkBrush(color=color)
            # Not quite sure why this has to get put in. Base pyqtgraph
            # may put circle default if pens are specified?
            kwargs["symbol"] = None


        if "width" in kwargs:
            kwargs["linewidth"] = kwargs.pop("width")
        width = kwargs.get('linewidth', config_options["linewidth"])

        pen.setColor(pg.mkColor(color))
        pen.setWidth(width)
        try:
            pen.setStyle(style)
        except TypeError:
            print("styleerror", style)
            raise
        kwargs['pen'] = pen
    return args, kwargs























































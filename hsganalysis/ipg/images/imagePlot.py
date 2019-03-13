import numpy as np
import pyqtgraph as pg
from .ImageViewWithPlotItemContainer import ImageViewWithPlotItemContainer

def image(*args, **kwargs):
    img = ImageViewWithPlotItemContainer()
    img.view.setAspectLocked(False)
    img.view.invertY(False)
    img.setImage(*args, **kwargs)


    return img



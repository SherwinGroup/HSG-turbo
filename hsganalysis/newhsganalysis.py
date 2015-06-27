# -*- coding: utf-8 -*-
"""
Created on Fri Jun 19 13:16:25 2015

@author: hbanks

Brevity required, prurience preferred
"""

from __future__ import division
import os
import glob
import numpy as np
import matplotlib.pyplot as plt

####################
# Objects 
####################

class CCD(object):
    def __init__(self):
    
class Photoluminescence(CCD):
    def __init__(self):

    def save_processing(self):
        
class Absorbance(CCD):
    def __init__(self):
    
    def save_processing(self):

class HighSideband(CCD):
    def __init__(self):
    
    def __add__(self, other):
    
    def __sub__(self, other):
    
    def add_std_error(self):
    
    def calc_approx_sb_order(self):
    
    def image_normalize(self):
    
    def guess_sidebands(self):
    
    def fit_sidebands(self):
        
    def save_processing(self):

class PMT(object):
    def __init__(self):
    
class HighSideband(PMT):
    def __init__(self):
    
    def laser_line(self):
        
    def fit_sidebands(self):
    
    def save_processing(self):

class FullSpectrum(object):
    def __init__(self):
        
class FullAbsorbance(FullSpectrum):
    """
    I'm imagining this will sew up absorption spectra, but I'm not at all sure
    how to do that at the moment.
    """
    def __init__(self):

class FullHighSideband(FullSpectrum):
    """
    I'm imagining this class is created with a base CCD file, then gobbles up
    other spectra that belong with it, then grabs the PMT object to normalize
    everything, assuming that PMT object exists.
    """
    def __init__(self):
    
    def add_CCD(self):
        
    def add_PMT(self):
####################
# Fitting functions 
####################

####################
# Collection functions 
####################

def sum_spectra():
    """
    This function will add together individual images appropriately.  It will
    not handle stitching.
    """
def stitch_results():
    """
    This will make FullSpectrum objects out of a list of processed spectra.
    """
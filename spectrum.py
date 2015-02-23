# -*- coding: utf-8 -*-
"""
Created on Sun Feb 15 15:22:30 2015

@author: dreadnought

"Brevity required, prurience preferred"
"""

from __future__ import division
import os, errno
import copy
import json
import numpy as np
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt

class Base_spectrum(object):
    '''
    This class will be for loading hsg spectra.  Currently, this will only load
    data output from the EMCCD.  I think future verions should be able to 
    combine PMT and EMCCD data, as well as stitching scans.  It will contain 
    all the methods useful for loading this kind of data.
    '''
    
    def __init__(self, fname):
        '''
        This will read the appropriate file.  The header needs to be fixed to
        reflect the changes to the output header from the Andor file.  Because
        another helper file will do the cleaning and background subtraction,
        those are no longer part of this init.
        
        fname = file name of the hsg spectrum
        '''
        self.fname = fname
        
        f = open(fname,'rU')
            
        self.description = f.readline()
        self.description = self.description[1:-3]
        parameters_str = f.readline()
        self.parameters = json.loads(parameters_str[1:])

        f.close()
        
        self.spectrum = np.genfromtxt(fname, comments='#', delimiter=',')
        
        # These will keep track of what operations have been performed
        self.shot_normalized = False
        self.addenda = self.parameters['addenda']
        self.subtrahenda = self.parameters['subtrahenda']
    
    def __add__(self, other):
        '''
        Add together the image data from self.hsg_data, or add a constant to 
        that np.array.  
        '''
        ret = copy.deepcopy(self)
        
        # Add a constant offset to the data
        if type(other) in (int, float):
            ret.hsg_data[:,1] = self.hsg_data[:,1] + other # What shall the name be?
            ret.addenda[0] = ret.addenda[0] + other
        
        # or add the data of two hsg_spectra together
        else:
            if np.isclose(ret.hsg_data[0,0], other.hsg_data[0,0]):
                ret.hsg_data[:,1] = self.hsg_data[:,1] + other.hsg_data[:,1] # Again, need to choose a name
                ret.addenda[0] = ret.addenda[0] + other.addenda[0]
                ret.addenda.extend(other.addenda[1:])
                ret.subtrahenda.extend(other.subtrahenda)
            else:
                raise Exception('Source: Spectrum.__add__\nThese are not from the same grating settings')
        return ret
        
    def __sub__(self, other):
        '''
        This subtracts constants or other data sets between self.hsg_data.  I 
        think it even keeps track of what data sets are in the file and how 
        they got there
        '''
        ret = copy.deepcopy(self)
        if type(other) in (int, float):
            ret.hsg_data[:,1] = self.hsg_data[:,1] - other # Need to choose a name
            ret.addenda[0] = ret.addenda[0] - other
        else:
            if np.isclose(ret.hsg_data[0,0], other.hsg_data[0,0]):
                ret.hsg_data[:,1] = self.hsg_data[:,1] - other.hsg_data[:,1]
                ret.subtrahenda.extend(other.addenda[1:])
                ret.addenda.extend(other.subtrahenda)
            else:
                raise Exception('These are not from the same grating settings')
        return ret
    
#    def inspect_dark_regions(self):
#        '''
#        The idea behind this method is to look at a region where there is no
#        input light and measure the mean and the statistical details of the 
#        noise there.  The mean will then be subtracted off so that there are 
#        minimal errors 
#        '''
#        raise NotImplementedError
    
    def guess_sidebands(self, mycutoff=4):
        '''
        This just implements the peak_detect function defined below.
        '''
        self.sb_index, sb_loc, sb_amp = peak_detect(self.hsg_data, cut_off=mycutoff)
        self.sb_guess = np.array([np.asarray(sb_loc), np.asarray(sb_amp)]).T
    
    def fit_sidebands(self, plot=False):
        '''
        This takes self.sb_guess and fits to each maxima to get the details of
        each sideband.
        '''
        self.sb_fits = []
        
        for index in self.sb_index:
            data_temp = self.hsg_data[index - 25:index + 25, :]
            p0 = [data_temp[25, 1], data_temp[25, 0], 0.1, 0.0]
            coeff, var_list = curve_fit(gauss, data_temp[:, 0], data_temp[:, 1], p0=p0)
            self.sb_fits.append(np.hstack((coeff, np.sqrt(np.diag(var_list)))))
            if plot:
                plt.plot(data_temp[:, 0], gauss(data_temp[:, 0], *coeff))
        
        self.sb_fits = np.asarray(self.sb_fits)
        
    def stitch_spectra(self):
        '''
        I really hope this is important later!
        '''        
        raise NotImplementedError

def peak_detect(data, window=20, cut_off=4):
    '''
    This function finds local maxima by comparing a point to the maximum and
    average of the absolute value of a neighborhood the size of 2*window.  I'm 
    sure this will take some tweaking.  
    
    data[0] = wavelength axis for the peaked data
    data[1] = signal we're looking for peaks in
    window = half the size of the neighborhood we look in.  
    cut_off = multiplier of the average of the absolute value for peak 
              identification.  This should probably be some noise parameter 
              calculated from the data.
    '''
    x_axis = data[:,0]
    y_axis = data[:,1]
    
    pre = np.empty(window)
    post = np.empty(window)
    
    pre[:] = -np.Inf
    post[:] = -np.Inf
    y_axis_temp = np.concatenate((pre, y_axis, post))
    
    max_index = []
    max_x = []
    max_y = []
    index = 0
    while index < len(x_axis):
        test_value = index + window
        
        check_y = y_axis_temp[test_value - window:test_value + window]
        check_max = check_y.max()
        check_ave = np.mean(abs(check_y[np.isfinite(check_y)])) # The inf's will dominate the mean
        check_value = y_axis_temp[test_value]
        
        if check_value == check_max and check_value > cut_off * check_ave:
            max_index.append(index)
            max_x.append(x_axis[index])
            max_y.append(y_axis[index])
            if check_value > 2 * cut_off * check_ave:
                index += window
            else:
                index += 1
        else:
            index += 1
        
    return max_index, max_x, max_y


def gauss(x, *p):
    A, mu, sigma, y0 = p
    return A*np.exp(-(x-mu)**2/(2.*sigma**2)) + y0
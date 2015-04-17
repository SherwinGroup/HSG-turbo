# -*- coding: utf-8 -*-
"""
Created on Sun Feb 15 15:22:30 2015

@author: dreadnought

"Brevity required, prurience preferred"
"""

from __future__ import division
import os
import errno
import copy
import json
import numpy as np
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt

class Spectrum(object):
    """
    This class will be for loading hsg spectra.  Currently, this will only load
    data output from the EMCCD.  I think future versions should be able to
    combine PMT and EMCCD data, as well as stitching scans.  It will contain 
    all the methods useful for loading this kind of data.
    """
    
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

        parameters_str = f.readline()
        self.parameters = json.loads(parameters_str[1:])
        self.description = ''
        read_description = True
        while read_description:
            line = f.readline()
            if line[0] == '#':
                self.description += line[1:]
            else:
                read_description = False

        f.close()

        
        self.raw_data = np.genfromtxt(fname, comments='#', delimiter=',')
        self.hsg_data = None
        
        # These will keep track of what operations have been performed
        self.initial_processing = False
        self.addenda = self.parameters['addenda']
        self.subtrahenda = self.parameters['subtrahenda']
    
    def __add__(self, other):
        """
        Add together the image data from self.hsg_data, or add a constant to 
        that np.array.  
        """
        if self.initial_processing == True:
            raise Exception('Source: Spectrum.__add__:\nToo much processing already!')
        ret = copy.deepcopy(self)

        # Add a constant offset to the data
        if type(other) in (int, float):
            ret.raw_data[:, 1] = self.raw_data[:, 1] + other # What shall the name be?
            ret.addenda[0] = ret.addenda[0] + other
        
        # or add the data of two hsg_spectra together
        else:
            if np.isclose(ret.parameters['center_lambda'], other.parameters['center_lambda']):
                ret.raw_data[:, 1] = self.raw_data[:, 1] + other.raw_data[:, 1] # Again, need to choose a name
                ret.addenda[0] = ret.addenda[0] + other.addenda[0]
                ret.addenda.extend(other.addenda[1:])
                ret.subtrahenda.extend(other.subtrahenda)
                ret.parameters['FEL_pulses'] += other.parameters['FEL_pulses']
            else:
                raise Exception('Source: Spectrum.__add__:\nThese are not from the same grating settings')
        return ret
        
    def __sub__(self, other):
        '''
        This subtracts constants or other data sets between self.hsg_data.  I 
        think it even keeps track of what data sets are in the file and how 
        they got there.

        I have no idea if this is something that will be useful?
        '''
        if self.initial_processing == True:
            raise Exception('Source: Spectrum.__sub__:\nToo much processing already!')
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
                raise Exception('Source: Spectrum.__sub__:\nThese are not from the same grating settings')
        return ret

    def __str__(self):
        return self.description

    def initial_process(self):
        """
        This method will divide everything by the number of FEL shots that contributed, as well as change wavelength
        into energy, nm to eV
        :return:
        """
        self.hsg_data = np.array(self.raw_data)
        self.hsg_data[:, 0] = 1239.84 / self.raw_data[:, 0]
        if self.parameters['FEL_pulses'] > 0:
            self.hsg_data[:, 1] = self.raw_data[:, 1] / self.parameters['FEL_pulses']
        else:
            self.parameters['FELP'] = "0"
        self.initial_processing = True
        self.parameters['shot_normalized'] = True

    def guess_sidebands(self, mycutoff=4):
        '''
        This just implements the peak_detect function defined below.
        '''
        self.sb_index, sb_loc, sb_amp = peak_detect(self.hsg_data, cut_off=mycutoff)
        self.sb_guess = np.array([np.asarray(sb_loc), np.asarray(sb_amp)]).T
        print self.sb_index
        print len(self.sb_guess[:, 0])
    
    def fit_sidebands(self, sensitivity=2.5, plot=False):
        '''
        This takes self.sb_guess and fits to each maxima to get the details of
        each sideband.
        
        sensitivity - the multiplier of the noise stdev required to save a peak
        plot - if you want to plot the fits on the same plot as before
        '''
        self.sb_fits = []
        
        for elem in xrange(len(self.sb_index)):
            data_temp = self.hsg_data[self.sb_index[elem] - 50:self.sb_index[elem] + 50, :]
            p0 = [self.sb_guess[elem, 0], self.sb_guess[elem, 1], 0.0001, 1.0]
            #print "This is the p0:", p0
            #print "Let's fit this shit!"
            try:
                coeff, var_list = curve_fit(gauss, data_temp[:, 0], data_temp[:, 1], p0=p0)
                coeff[2] = abs(coeff[2]) # The linewidth shouldn't be negative
                #print coeff
                #coeff = coeff[:4]
                noise_stdev = (np.std(data_temp[:45]) + np.std(data_temp[55:])) / 2
                print coeff[1], " vs. ", noise_stdev
                if 1e-4 > coeff[2] > 20e-6 and coeff[1] > sensitivity * noise_stdev:
                    self.sb_fits.append(np.hstack((coeff, np.sqrt(np.diag(var_list)))))
                    if plot:
                        x_vals = np.linspace(data_temp[0, 0], data_temp[-1, 0], num=200)
                        plt.plot(x_vals, gauss(x_vals, *coeff))
            except:
                print "I couldn't fit that"        
        sb_fits_temp = np.asarray(self.sb_fits)
        reorder = [0, 4, 1, 5, 2, 6, 3, 7]
        #print "The temp fits list", sb_fits_temp
        try:
            self.sb_fits = sb_fits_temp[:, reorder]
        except:
            self.sb_fits = list(sb_fits_temp)
    
    def fit_sidebands_for_NIR_freq(self, sensitivity=2.5, plot=False):
        '''
        This takes self.sb_guess and fits to each maxima to get the details of
        each sideband.
        
        sensitivity - the multiplier of the noise stdev required to save a peak
        plot - if you want to plot the fits on the same plot as before
        '''
        self.sb_fits = []
        
        for elem in xrange(len(self.sb_index)):
            data_temp = self.hsg_data[self.sb_index[elem] - 50:self.sb_index[elem] + 50, :]
            p0 = [self.sb_guess[elem, 0], self.sb_guess[elem, 1], 0.0001, 1.0]
            #print "This is the p0:", p0
            #print "Let's fit this shit!"
            try:
                coeff, var_list = curve_fit(gauss, data_temp[:, 0], data_temp[:, 1], p0=p0)
                coeff[2] = abs(coeff[2]) # The linewidth shouldn't be negative
                #print coeff
                #coeff = coeff[:4]
                noise_stdev = (np.std(data_temp[:45]) + np.std(data_temp[55:])) / 2
                #print coeff[1], " vs. ", noise_stdev
                if 1e-4 > coeff[2] > 20e-6 and coeff[1] > sensitivity * noise_stdev:
                    #print "Going to redefine stuff", coeff[0], 1239.84 / float(self.parameters["NIR_lambda"])
                    coeff[0] = round((coeff[0] - (1239.84 / float(self.parameters["NIR_lambda"]))) / .002515, 1)
                    print "New coeff[0] ", coeff[0]
                    self.sb_fits.append(np.hstack((coeff, np.sqrt(np.diag(var_list)))))
                    if plot:
                        x_vals = np.linspace(data_temp[0, 0], data_temp[-1, 0], num=200)
                        plt.plot(x_vals, gauss(x_vals, *coeff))
            except:
                print "I couldn't fit that"        
        sb_fits_temp = np.asarray(self.sb_fits)
        reorder = [0, 4, 1, 5, 2, 6, 3, 7]
        #print "The temp fits list", sb_fits_temp
        try:
            self.sb_fits = sb_fits_temp[:, reorder]
        except:
            self.sb_fits = list(sb_fits_temp)

    
    def save_processing(self, file_name, folder_str, marker, index):
        """
        This will save all of the results from data processing
        :param file_prefix:
        :return:
        """
        try:
            os.mkdir(folder_str)
        except OSError, e:
            if e.errno == errno.EEXIST:
                pass
            else:
                raise

        spectra_fname = file_name + '_' + marker + '_' + str(index) + '.txt'
        fit_fname = file_name + '_' + marker + '_' + str(index) + '_fits.txt'
        self.parameters['addenda'] = self.addenda
        self.parameters['subtrahenda'] = self.subtrahenda
        try:
            parameter_str = json.dumps(self.parameters, sort_keys=True)
        except:
            print "Source: EMCCD_image.save_images\nJSON FAILED"
            print self.parameters
            return

        origin_import_spec = '\nWavelength,Signal\neV,arb. u.'
        spec_header = '#' + parameter_str + '\n' + '#' + self.description[:-2] + origin_import_spec
        #print "Spec header: ", spec_header
        origin_import_fits = '\nCenter energy,error,Amplitude,error,Linewidth,error,Constant offset,error\neV,,arb. u.,,eV,,arb. u.,\n,,' + marker
        fits_header = '#' + parameter_str + '\n' + '#' + self.description[:-2] + origin_import_fits
        #print "Fits header: ", fits_header
        np.savetxt(os.path.join(folder_str, spectra_fname), self.hsg_data, delimiter=',',
                   header=spec_header, comments='', fmt='%f')
        np.savetxt(os.path.join(folder_str, fit_fname), self.sb_fits, delimiter=',',
                   header=fits_header, comments='', fmt='%f')

        print "Save image.\nDirectory: {}".format(os.path.join(folder_str, spectra_fname))

    def stitch_spectra(self):
        '''
        I really hope this is important later!
        '''        
        raise NotImplementedError

def peak_detect(data, window=30, cut_off=5.5):
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
        test_value = index + window # The value at test_value will be checked to see if it is the max
        
        check_y = y_axis_temp[test_value - window:test_value + window] # This is the region to test over
        check_max = check_y.max()

        check_ave = np.mean(abs(check_y[np.isfinite(check_y)])) # The inf's will dominate the mean
        check_value = y_axis_temp[test_value]
        check_value_neighborhood = sum(y_axis_temp[test_value - 2:test_value + 2]) # Since peaks are broader than one pixel
        if check_value == check_max and check_value_neighborhood > cut_off * check_ave:
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
    mu, A, sigma, y0 = p
    return A * np.exp(-(x - mu)**2 / (2. * sigma**2)) + y0

def lingauss(x, *p):
    mu, A, sigma, y0, m = p
    return A * np.exp(-(x - mu)**2 / (2. * sigma**2)) + y0 + m*x

def lorentzian(x, *p):
    mu, A, gamma, y0 = p
    return (A * gamma**2) / ((x - mu)**2 + gamma**2) + y0

def sum_spectra(object_list):
    """
    This function will add all the things that should be added.  Obvs.

    object_list: A list of spectrum objects
    :param object_list:
    :return:
    """
    good_list = []
    for index in xrange(len(object_list)):
        try:
            temp = object_list.pop(0)
        except:
            break
#        print "temp has series: {}.\ttemp has cl: {}.\ttemp has fn: {}".format(temp.parameters['series'], temp.parameters['center_lambda'], temp.fname[-16:-13])
        for spec in list(object_list):
#            print "\tspec has series: {}.\tspec has cl: {}.\tspec has fn: {}".format(spec.parameters['series'], spec.parameters['center_lambda'], spec.fname[-16:-13])
#            print "I am trying to add", temp.parameters['FELP'], spec.parameters['FELP']
            if temp.parameters['series'] == spec.parameters['series']:
                if temp.parameters['center_lambda'] == spec.parameters['center_lambda']:
                    temp += spec
#                    print "\t\tadded"
                    #print "I ADDED", temp.parameters['FELP'], spec.parameters['FELP']
                    object_list.remove(spec)
        good_list.append(temp)
    return good_list




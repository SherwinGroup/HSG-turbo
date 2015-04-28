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
        """
        This will read the appropriate file.  The header needs to be fixed to
        reflect the changes to the output header from the Andor file.  Because
        another helper file will do the cleaning and background subtraction,
        those are no longer part of this init.  This also turns all wavelengths
        from nm (NIR ones) or cm-1 (THz ones) into eV.
        
        fname = file name of the hsg spectrum
        """
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

        self.parameters["NIR_freq"] = 1239.84 / float(self.parameters["NIR_lambda"])
        self.parameters["THz_freq"] = 0.000123984 * float(self.parameters["FEL_lambda"])
        self.raw_data = np.flipud(np.genfromtxt(fname, comments='#', delimiter=','))
        self.raw_data = np.array(self.raw_data[:1600,:])
        self.raw_data[:, 0] = 1239.84 / self.raw_data[:, 0]
        
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
        """
        This subtracts constants or other data sets between self.hsg_data.  I 
        think it even keeps track of what data sets are in the file and how 
        they got there.

        I have no idea if this is something that will be useful?
        """
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
        This method will divide everything by the number of FEL shots that contributed
        :return:
        """
        self.hsg_data = np.array(self.raw_data)

        if self.parameters['FEL_pulses'] > 0:
            self.hsg_data[:, 1] = self.raw_data[:, 1] / self.parameters['FEL_pulses']
        else:
            self.parameters['FELP'] = "0"
        self.initial_processing = True
        self.parameters['shot_normalized'] = True

    def guess_sidebands(self, cutoff=4):
        """
        This method finds all the sideband candidates in hsg_data.  It first
        finds the lowest order sideband that could be in the data by looking at
        the first wavelength (eV) in the x_axis and calls that sb_init. 
        
        NB: Currently the lowest order sideband that we can measure is the laser,
        so order will be initialized to zero and the first peak that we'd measure
        would be the laser peak. 
        
        Next, the algorithm will look for that lowest sideband in a 30% window
        around where the NIR frequency should be.  It will take the maximum of
        that window and add it with the two neighboring points.  If this is 
        larger than the cutoff, then it will save that sideband order, the index
        of the point, the frequency and the amplitude.  It will look for odd
        sidebands as well as even ones, and if it cannot find two sidebands in 
        a row (usually the second missed sideband will be an even one), then it
        will stop looking.
        
        inputs:
        cutoff: the size of check_max_area must be relative to check_ave to be
                called a sideband.
        
        important local variables:
        x_axis: a np.array consisting of the frequency axis of self.hsg_data.
        y_axis: a np.array consisting of the signal axis of self.hsg_data.
        sb_init: the smallest-order sideband that the frequency axis could see.
        check_max_index: index of maximum value in the 30% region
        check_max: amplitude of the sideband candidate, important for fitting
        check_max_area: approximate integral of the sideband candidate
        
        attribute outputs:
        self.sb_list: list of sideband orders detected by this method
        self.sb_index: list of the index of the sideband
        self.sb_guess: 2D array that contains the frequency and amplitude of
                       the found sidebands.
        """
        x_axis = np.array(self.hsg_data[:, 0])
        y_axis = np.array(self.hsg_data[:, 1])
        
        window = 20
        
        pre = np.empty(window)
        post = np.empty(window)
        
        pre[:] = -np.Inf
        post[:] = -np.Inf
        y_axis_temp = np.concatenate((pre, y_axis, post))
        
        NIR_freq = self.parameters["NIR_freq"]
        THz_freq = self.parameters["THz_freq"]
        sb_init = False
        for order in xrange(50):
            print "I'm checking for sb energy", NIR_freq + order * THz_freq
            print "That's order", order
            if x_axis[0] < (NIR_freq + order * THz_freq): #and (x_axis[0] > NIR_freq + (order - 1) * THz):
                print "Lowest x_axis value is", x_axis[0]
                sb_init = order
                break
            elif order == 49:
                raise Exception("Source: self.guess_sidebands.\nCouldn't find sb_init.")
        print "We think the lowest sideband order is", sb_init
        
        self.sb_list = []
        self.sb_index = []
        sb_freq_guess = []
        sb_amp_guess = []
        
        last_sb = NIR_freq + THz_freq * (sb_init - 1)
        index_guess = 0
        consecutive_null_sb = 0
        for order in xrange(sb_init, 50):
            lo_freq_bound = last_sb + THz_freq * (1 - 0.15)
            hi_freq_bound = last_sb + THz_freq * (1 + 0.15)
            start_index = False
            end_index = False
            print "\nSideband", order, "\n"            
            for i in xrange(index_guess, 1600):
                if start_index == False and x_axis[i] > lo_freq_bound:
                    print "start_index is", i
                    start_index = i
                elif end_index == False and x_axis[i] > hi_freq_bound:
                    end_index = i
                    print "end_index is", i
                    index_guess = i
                    break
                
            check_y = y_axis_temp[start_index + window:end_index + window]
            print "check_y is", check_y
            check_max = check_y.max()
            print "check_max is", check_max
            check_max_index = np.argmax(check_y) # This assumes that two floats won't be identical
            check_max_area = np.sum(check_y[check_max_index - 1:check_max_index + 1])
            check_ave = np.mean(abs(check_y[np.isfinite(check_y)]))
            print "check_ave is", check_ave
            if check_max_area > cutoff * check_ave:
                found_index = np.argmax(check_y) + start_index
                self.sb_index.append(found_index)
                last_sb = x_axis[found_index]
                print "I just found", last_sb
                
                sb_freq_guess.append(x_axis[found_index])
                sb_amp_guess.append(y_axis[found_index])
                self.sb_list.append(order)
                consecutive_null_sb = 0
            else:
                print "I could not find sideband with order", order
                last_sb = last_sb + THz_freq
                consecutive_null_sb += 1
            
            if consecutive_null_sb == 2:
                print "I can't find any more sidebands"
                break
        
        print "I found these sidebands:", self.sb_list
        self.sb_guess = np.array([np.asarray(sb_freq_guess), np.asarray(sb_amp_guess)]).T
    
    def fit_sidebands(self, plot=False):
        """
        This takes self.sb_guess and fits to each maxima to get the details of
        each sideband.
        """
        self.sb_fits = []
        
        for elem in xrange(len(self.sb_index)):
            data_temp = self.hsg_data[self.sb_index[elem] - 25:self.sb_index[elem] + 25, :]
            p0 = [self.sb_guess[elem, 0], self.sb_guess[elem, 1], 0.0001, 1.0]
            #print "This is the p0:", p0
            #print "Let's fit this shit!"
            try:
                coeff, var_list = curve_fit(gauss, data_temp[:, 0], data_temp[:, 1], p0=p0)
            except:
                coeff = np.array(p0)
                var_list = np.array([p0, p0, p0, p0])
                print "Had to make coeff equal to p0"
            coeff[2] = abs(coeff[2]) # The linewidth shouldn't be negative
            #print coeff
            #coeff = coeff[:4]
            if coeff[0] > 0.:
                self.sb_fits.append(np.hstack((coeff, np.sqrt(np.diag(var_list)))))
                if plot:
                    x_vals = np.linspace(data_temp[0, 0], data_temp[-1, 0], num=200)
                    plt.plot(x_vals, gauss(x_vals, *coeff))
        
        sb_fits_temp = np.asarray(self.sb_fits)
        reorder = [0, 4, 1, 5, 2, 6, 3, 7]
        #print "The temp fits list", sb_fits_temp
        try:
            self.sb_fits = sb_fits_temp[:, reorder]
        except:
            self.sb_fits = list(sb_fits_temp)
            print "\n!!!!!!!!!\nSome shit went wrong here!!!\n!!!!!!!!\n"
        
        # Going to label the appropriate row with the sideband
        sb_names = np.vstack(self.sb_list)
        self.sb_results = np.hstack((sb_names, self.sb_fits[:,:6]))

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

        spectra_fname = file_name + str(index) + '.txt'
        self.save_name = spectra_fname
        fit_fname = file_name + str(index) + '_fits.txt'
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
        """
        I really hope this is important later!
        """        
        raise NotImplementedError

####################
# Useful functions 
####################

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

def save_parameter_sweep(spectrum_list, file_name, folder_str, param_name, unit):
    """
    This function will take a fully processed list of spectrum objects and 
    slice Spectrum.sb_fits appropriately to get an output like:
    
    "Parameter" | SB1 freq | err | SB1 amp | error | SB1 linewidth | error | SB2...| SBn...|
    param1      |    .     |
    param2      |    .     |
      .
      .
      .
    
    Currently I'm thinking fuck the offset y0
    After constructing this large matrix, it will save it somewhere.
    """
    #included_spectra = dict()
    param_array = None
    sb_included = None    
    for spectrum in spectrum_list:
        temp_1d = np.hstack((spectrum.parameters[param_name], spectrum.sb_results.ravel()))
        #included_spectra[spectrum.save_name] = spectrum.parameters[param_name]
        if param_array is None:
            param_array = np.array(temp_1d)
            sb_included = list(spectrum.sb_list)
            print "param_array initialized:", param_array
        else:
            if sb_included == spectrum.sb_list:
                param_array = np.vstack((param_array, temp_1d))
                print "param_array added to the easy way:", param_array
            else:   
                spec_list = list(spectrum.sb_list)
                perm_list = list(sb_included)
                print "spec_list:", spec_list
                print "perm_list:", perm_list
                temp_order = 0
                while temp_order < 50:
                    print "\ntemp order is:", temp_order
                    print "spec_list:", spec_list
                    print "perm_list:", perm_list
                    print param_array
                    if spec_list == [] and perm_list == []:
                        break
                    elif temp_order == spec_list[0] and temp_order == perm_list[0]:
                        spec_list.pop(0)
                        perm_list.pop(0)
                        #temp_order += 1
                        #continue
                    elif temp_order == spec_list[0] and perm_list == []:
                        blank = np.zeros((len(param_array[:, 0]), 7))
                        blank[:, 0] = temp_order
                        param_array = np.hstack((param_array, blank))
                        sb_included = sorted([temp_order] + sb_included)
                    elif temp_order == perm_list[0] and spec_list == []:
                        blank = np.zeros(7)
                        blank[0] = temp_order
                        temp_1d = np.hstack((temp_1d, blank))
                        spec_list = [temp_order] + spec_list
                    elif temp_order == spec_list[0] and temp_order < perm_list[0]:
                        param_array_temp = np.hsplit(param_array, [np.nonzero(int(round(param_array[0, :])) == perm_list[0])[0]])
                        blank = np.zeros((len(param_array[:, 0]), 7))
                        blank[:, 0] = temp_order
                        param_array = np.hstack((param_array_temp[0], blank, param_array_temp[1]))
                        perm_list = [temp_order] + perm_list
                        sb_included = sorted([temp_order] + sb_included)
                        #continue
                    elif temp_order == perm_list[0] and temp_order < spec_list[0]:
                        temp = np.hsplit(temp_1d, [np.nonzero(int(round(temp_1d)) == spec_list[0])])
                        blank = np.zeros(7)
                        blank[0] = temp_order
                        temp_1d = np.hstack((temp[0], blank, temp[1]))
                        spec_list = [temp_order] + spec_list
                        #continue
                    temp_order += 1
                param_array = np.vstack((param_array, temp_1d))
                print "param_array added to the hard way:", param_array
    print param_array
    '''
    try:
        os.mkdir(folder_str)
    except OSError, e:
        if e.errno == errno.EEXIST:
            pass
        else:
            raise

    file_name = file_name + '.txt'
    
    try:
        included_spectra_str = json.dumps(included_spectra, sort_keys=True)
    except:
        print "Source: save_parameter_sweep\nJSON FAILED"
        return
    origin_import1 = param_name
    origin_import2 = unit
    for order in sb_included:
        origin_import1 += ",Sideband,Frequency,error,Amplitude,error,Linewidth,error"
        origin_import2 += ",order,eV,,arb. u.,,eV,"
    origin_total = origin_import1 + "\n" + origin_import2 + "\n"
    header = '#' + included_spectra_str + '\n' + origin_total
    #print "Spec header: ", spec_header

    np.savetxt(os.path.join(folder_str, file_name), param_array, delimiter=',',
               header=header, comments='', fmt='%f')

    print "Saved the file.\nDirectory: {}".format(os.path.join(folder_str, file_name))
    '''
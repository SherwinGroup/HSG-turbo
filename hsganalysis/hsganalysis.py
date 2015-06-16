# -*- coding: utf-8 -*-
"""
Created on Sun Feb 15 15:22:30 2015

@author: dreadnought

"Brevity required, prurience preferred"
"""

from __future__ import division
import os
import glob
import errno
import copy
import json
import numpy as np
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt

####################
# Objects 
####################

class CCD(object):
    """
    This class will be for loading CCD spectra.  I think future versions should
    be able to combine PMT and EMCCD data, as well as stitching scans.  It will
    contain all the methods useful for loading this kind of data.
    """
    
    def __init__(self, fname):
        """
        This will read the appropriate file.  The header needs to be fixed to
        reflect the changes to the output header from the Andor file.  Because
        another helper file will do the cleaning and background subtraction,
        those are no longer part of this init.  This also turns all wavelengths
        from nm (NIR ones) or cm-1 (THz ones) into eV.
        
        Input:
        fname = file name of the hsg spectrum
        
        Internal:
        self.fname = the filename
        self.parameters = string with all the relevant experimental perameters
        self.description = the description we added to the file as the data
                           was being taken
        self.raw_data = the unprocessed data that is wavelength vs signal
        self.hsg_data = processed data that has gone is frequency vs signal/pulse
        self.dark_stdev = the standard deviation of the dark noise using the
                          last 200 pixels since they probably(?) don't have 
                          strong signal
        self.std_error = Nothing yet, maybe this isn't useful
        self.addenda = the list of things that have been added to the file, in
                       form of [constant, *spectra_added]
        self.subtrahenda = the list of spectra that have been subtracted from
                           the file.  Constant subtraction is dealt with with
                           self.addenda
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
        '''
        ###
        # This section is for old code compatability
        ###
        self.parameters["nir_lambda"] = self.parameters["NIR_lambda"] 
        self.parameters["fel_lambda"] = self.parameters["FEL_lambda"]
        self.parameters["fel_pulses"] = self.parameters["FEL_pulses"]
        ### End section
        '''
        self.parameters["NIR_freq"] = 1239.84 / float(self.parameters["nir_lambda"])
        self.parameters["THz_freq"] = 0.000123984 * float(self.parameters["fel_lambda"])
        self.parameters["nir_power"] = float(self.parameters["nir_power"])
        self.parameters["thz_power"] = float(self.parameters["fel_power"])
        self.raw_data = np.flipud(np.genfromtxt(fname, comments='#', delimiter=','))
        # I used flipup so that the x-axis is an increasing function of frequency
        
        self.hsg_data = np.array(self.raw_data[:1600,:]) # By slicing up to 1600,
                                                         # we cut out the text 
                                                         # header
        self.hsg_data[:, 0] = 1239.84 / self.hsg_data[:, 0]
        # Need to do this next line AFTER adding all the spectra together
        self.hsg_data[:, 1] = self.hsg_data[:, 1] / self.parameters['fel_pulses']
        
        self.dark_stdev = self.parameters["background_darkcount_std"] / self.parameters['fel_pulses'] # What do I do with this now that hsg_data is not normalized?
        self.std_error = None
        # These will keep track of what operations have been performed
        
        self.addenda = self.parameters['addenda']
        self.subtrahenda = self.parameters['subtrahenda']
        #self.sb_results = None
    
    def __add__(self, other):
        """
        Add together the image data from self.hsg_data, or add a constant to 
        that np.array.  
        """
        ret = copy.deepcopy(self)

        # Add a constant offset to the data
        if type(other) in (int, float):
            ret.hsg_data[:, 1] = self.hsg_data[:, 1] + other # What shall the name be?
            ret.addenda[0] = ret.addenda[0] + other
        
        # or add the data of two hsg_spectra together
        else:
            if np.isclose(ret.parameters['center_lambda'], other.parameters['center_lambda']):
                ret.hsg_data[:, 1] = self.hsg_data[:, 1] + other.hsg_data[:, 1] # Again, need to choose a name
                ret.addenda[0] = ret.addenda[0] + other.addenda[0]
                ret.addenda.extend(other.addenda[1:])
                ret.subtrahenda.extend(other.subtrahenda)
                ret.parameters['fel_pulses'] += other.parameters['fel_pulses']
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
        ret = copy.deepcopy(self)
        
        # Subtract a constant offset to the data
        if type(other) in (int, float):
            ret.hsg_data[:,1] = self.hsg_data[:,1] - other # Need to choose a name
            ret.addenda[0] = ret.addenda[0] - other
            
        # Subtract the data of two hsg_spectra from each other
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

#    def get_dark_stdev(self):
#        """
#        This method will divide everything by the number of FEL shots that contributed
#        :return:
#        """
#        self.dark_stdev = np.std(self.hsg_data[1400:1600, 1])
        
    def add_std_error(self, std_errors):
        """
        This adds a numpy array of standard errors to the self.hsg_data array.
        It's actually just an append_column method, but there's no reason for 
        that (I think)
        """
        try:
            
            self.hsg_data = np.hstack((self.hsg_data, std_errors.reshape((1600, 1))))
        except:
            print "Spectrum.add_std_error fucked up.  What's wrong?"
            print std_errors.shape
    
    def calc_approx_sb_order(self, freq):
        """
        This simple method will simply return a float approximating the order
        of the frequency input.
        """
        nir_freq = self.parameters['NIR_freq']
        thz_freq = self.parameters['THz_freq']
        approx_order = (freq - nir_freq) / thz_freq
        return approx_order
    
    def image_normalize(self, num_images):
        """
        This method will divide the hsg_data by the number of images that were
        used in the sum_spectra function.
        """
        self.hsg_data[:, 1] = self.hsg_data[:, 1] / num_images
        self.hsg_data[:, 2] = self.hsg_data[:, 2] / num_images

    
    def guess_better(self, cutoff=4.5):
        """
        This method will replace the previous one for guessing sidebands.
        """
        
        x_axis = np.array(self.hsg_data[:, 0])
        y_axis = np.array(self.hsg_data[:, 1])
        error = np.array(self.hsg_data[:, 2])
        #window = 20
        
        min_sb = int(self.calc_approx_sb_order(x_axis[0])) + 1
        max_sb = int(self.calc_approx_sb_order(x_axis[-1]))
        
        #pre = np.empty(window)
        #post = np.empty(window)
        
        #pre[:] = -np.Inf
        #post[:] = -np.Inf
        #y_axis_temp = np.concatenate((pre, y_axis, post))
        
        nir_freq = self.parameters["NIR_freq"]
        thz_freq = self.parameters["THz_freq"]
        
        # Find max strength sideband and it's order
        global_max = np.argmax(y_axis)
        order_init = int(round(self.calc_approx_sb_order(x_axis[global_max])))

        check_y = y_axis[global_max - 15:global_max + 15]
        print "Global max checking:", check_y
        check_max_area = np.sum(y_axis[global_max - 1:global_max + 2])
        check_ave = np.mean(check_y)
        check_stdev = np.std(check_y)
        print "\ncheck_max_area is", check_max_area
        print "check_ave is", check_ave
        print "check_stdev is", check_stdev
        check_ratio = (check_max_area - 3 * check_ave) / check_stdev
        print "check_ratio is", check_ratio
        if check_ratio > cutoff:
            self.sb_list = [order_init]
            self.sb_index = [global_max]
            sb_freq_guess = [x_axis[global_max]]
            sb_amp_guess = [y_axis[global_max]]
            sb_error_est = [np.sqrt(sum([i**2 for i in error[global_max - 1:global_max + 2]])) / (check_max_area - 3 * check_ave)]
        else:
            print "There are no sidebands in", self.fname
            assert False
        
        
        # Look for lower order sidebands
        
        last_sb = sb_freq_guess[0]
        index_guess = global_max
        consecutive_null_sb = 0
        consecutive_null_odd = 0
        no_more_odds = False
        break_condition = False
        for order in xrange(order_init - 1, min_sb - 1, -1):
            if no_more_odds == True and order % 2 == 1:
                last_sb = last_sb - thz_freq
                continue
            lo_freq_bound = last_sb - thz_freq * (1 + 0.22) # Not sure what to do about these
            hi_freq_bound = last_sb - thz_freq * (1 - 0.22)
            start_index = False
            end_index = False
            print "\nSideband", order, "\n"            
            for i in xrange(index_guess, 0, -1):
                if end_index == False and i == 1:
                    break_condition = True
                    break
                if end_index == False and x_axis[i] < hi_freq_bound:
                    print "end_index is", i
                    end_index = i
                elif i == 1:
                    start_index = 0
                    print "hit end of data, start_index is 0"
                elif start_index == False and x_axis[i] < lo_freq_bound:
                    start_index = i
                    print "start_index is", i
                    index_guess = i
                    break
            
            if break_condition:
                break
            check_y = y_axis[start_index:end_index]
            print "check_y is", check_y
            #check_max = check_y.max()
            #print "check_max is", check_max
            check_max_index = np.argmax(check_y) # This assumes that two floats won't be identical
            check_max_area = np.sum(check_y[check_max_index - 1:check_max_index + 2])
            check_ave = np.mean(check_y)
            check_stdev = np.std(check_y)
            check_ratio = (check_max_area - 3 * check_ave) / check_stdev
            print "\ncheck_max_area is", check_max_area
            print "check_ave is", check_ave
            print "check_stdev is", check_stdev
            print "check_ratio is", check_ratio
            
            if check_ratio > cutoff:
                found_index = check_max_index + start_index
                self.sb_index.append(found_index)
                last_sb = x_axis[found_index]
                print "I just found", last_sb
                
                sb_freq_guess.append(x_axis[found_index])
                sb_amp_guess.append(check_max_area - 3 * check_ave)
                error_est = np.sqrt(sum([i**2 for i in error[found_index - 1:found_index + 2]])) / (check_max_area - 3 * check_ave)
                print "My error estimate is:", error_est
                sb_error_est.append(error_est)
                self.sb_list.append(order)
                consecutive_null_sb = 0
                if order % 2 == 1:
                    consecutive_null_odd = 0
            else:
                print "I could not find sideband with order", order
                last_sb = last_sb + thz_freq
                consecutive_null_sb += 1
                if order % 2 == 1:
                    consecutive_null_odd += 1
            if consecutive_null_odd == 1 and no_more_odds == False:
                print "I'm done looking for odd sidebands"
                no_more_odds = True
            if consecutive_null_sb == 2:
                print "I can't find any more sidebands"
                break  
        
        # Look for higher sidebands
        
        last_sb = sb_freq_guess[0]
        index_guess = global_max
        consecutive_null_sb = 0
        consecutive_null_odd = 0
        no_more_odds = False
        break_condition = False
        for order in xrange(order_init + 1, max_sb + 1):
            if no_more_odds == True and order % 2 == 1:
                last_sb = last_sb + thz_freq
                continue
            lo_freq_bound = last_sb + thz_freq * (1 - 0.22) # Not sure what to do about these
            hi_freq_bound = last_sb + thz_freq * (1 + 0.22)
            start_index = False
            end_index = False
            print "\nSideband", order, "\n"            
            for i in xrange(index_guess, 1600):
                if start_index == False and i == 1599:
                    print "I'm all out of space, captain!"
                    break_condition = True
                    break
                elif start_index == False and x_axis[i] > lo_freq_bound:
                    print "start_index is", i
                    start_index = i
                elif i == 1599:
                    end_index = 1599
                    print "hit end of data, end_index is 1599"
                elif end_index == False and x_axis[i] > hi_freq_bound:
                    end_index = i
                    print "end_index is", i
                    index_guess = i
                    break
            if break_condition:
                break
            check_y = y_axis[start_index:end_index]
            print "check_y is", check_y
            #check_max = check_y.max()
            #print "check_max is", check_max
            check_max_index = np.argmax(check_y) # This assumes that two floats won't be identical
            check_max_area = np.sum(check_y[check_max_index - 1:check_max_index + 2])
            check_ave = np.mean(check_y)
            check_stdev = np.std(check_y)
            check_ratio = (check_max_area - 3 * check_ave) / check_stdev
            print "\ncheck_max_area is", check_max_area
            print "check_ave is", check_ave
            print "check_stdev is", check_stdev
            print "check_ratio is", check_ratio
            
            if check_ratio > cutoff:
                found_index = check_max_index + start_index
                self.sb_index.append(found_index)
                last_sb = x_axis[found_index]
                print "I just found", last_sb
                
                sb_freq_guess.append(x_axis[found_index])
                sb_amp_guess.append(check_max_area - 3 * check_ave)
                error_est = np.sqrt(sum([i**2 for i in error[found_index - 1:found_index + 2]])) / (check_max_area - 3 * check_ave)
                print "My error estimate is:", error_est
                sb_error_est.append(error_est)
                self.sb_list.append(order)
                consecutive_null_sb = 0
                if order % 2 == 1:
                    consecutive_null_odd = 0
            else:
                print "I could not find sideband with order", order
                last_sb = last_sb + thz_freq
                consecutive_null_sb += 1
                if order % 2 == 1:
                    consecutive_null_odd += 1
            if consecutive_null_odd == 1 and no_more_odds == False:
                print "I'm done looking for odd sidebands"
                no_more_odds = True
            if consecutive_null_sb == 2:
                print "I can't find any more sidebands"
                break  
        
        print "I found these sidebands:", self.sb_list
        self.sb_guess = np.array([np.asarray(sb_freq_guess), np.asarray(sb_amp_guess), np.asarray(sb_error_est)]).T

        
    def guess_sidebands(self, cutoff=4.5):
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
                called a sideband.  3 is too insensitive, 4 might be too much.  
                3.5 might be slightly to insensitive.  I've seen definite 
                sidebands at 3.7 and definitely bad sidebands at 3.9.
        
        important local variables:
        x_axis: a np.array consisting of the frequency axis of self.hsg_data.
        y_axis: a np.array consisting of the signal axis of self.hsg_data.
        sb_init: the smallest-order sideband that the frequency axis could see.
        found_anything: boolean that keeps the method looking for sidebands 
                        if they aren't within two orders of the start.
        break_condition: boolean that will stop the method from looking past 
                         the edge of the spectrum.
        check_max_index: index of maximum value in the 30% region
        check_max: amplitude of the sideband candidate, important for fitting
        check_max_area: approximate integral of the sideband candidate
        
        attribute outputs:
        self.sb_list: list of sideband orders detected by this method
        self.sb_index: list of the index of the sideband
        self.sb_guess: 2D array that contains the frequency and amplitude of
                       the found sidebands.  It now also contains an estimate
                       of the error of area of the sidebands by adding the 
                       three standard errors nearest the peak.  Seems to work?
        """
        x_axis = np.array(self.hsg_data[:, 0])
        y_axis = np.array(self.hsg_data[:, 1])
        error = np.array(self.hsg_data[:, 2])
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
        sb_error_estimate = []
        
        last_sb = NIR_freq + THz_freq * (sb_init - 1)
        index_guess = 0
        consecutive_null_sb = 0
        consecutive_null_odd = 0
        break_condition = False
        no_more_odds = False
        found_anything = False
        
        for order in xrange(sb_init, 50):
            if no_more_odds == True and order % 2 == 1:
                last_sb = last_sb + THz_freq
                continue
            lo_freq_bound = last_sb + THz_freq * (1 - 0.22) # Not sure what to do about these
            hi_freq_bound = last_sb + THz_freq * (1 + 0.22)
            start_index = False
            end_index = False
            print "\nSideband", order, "\n"            
            for i in xrange(index_guess, 1600):
                if start_index == False and x_axis[i] > lo_freq_bound:
                    print "start_index is", i
                    start_index = i
                elif i == 1599:
                    break_condition = True
                    print "hit end of data, end_index is", i
                    break
                elif end_index == False and x_axis[i] > hi_freq_bound:
                    end_index = i
                    print "end_index is", i
                    index_guess = i
                    break
            if break_condition:
                break
            
            check_y = y_axis_temp[start_index + window:end_index + window]
            print "check_y is", check_y
            #check_max = check_y.max()
            #print "check_max is", check_max
            check_max_index = np.argmax(check_y) # This assumes that two floats won't be identical
            print "points in check_max_area:", check_y[check_max_index - 1:check_max_index + 2]
            check_max_area = np.sum(check_y[check_max_index - 1:check_max_index + 2])
            check_ave = np.mean(check_y[np.isfinite(check_y)])
            check_stdev = np.std(check_y[np.isfinite(check_y)])
            print "\ncheck_max_area is", check_max_area
            print "check_ave is", check_ave
            print "check_stdev is", check_stdev
            print "check_ratio is", (check_max_area - 3 * check_ave) / check_stdev
            
            if (check_max_area - 3 * check_ave) > cutoff * check_stdev:
                found_index = np.argmax(check_y) + start_index
                self.sb_index.append(found_index)
                last_sb = x_axis[found_index]
                print "I just found", last_sb
                
                sb_freq_guess.append(x_axis[found_index])
#                sb_amp_guess.append(y_axis[found_index])
                sb_amp_guess.append(check_max_area - 3 * check_ave)
                error_est = np.sqrt(sum([i**2 for i in error[found_index - 1:found_index + 2]])) / (check_max_area - 3 * check_ave)
                # I think that should be a "+ 2" in the second index of slicing?
                print "My error estimate is:", error_est
                sb_error_estimate.append(error_est)
                self.sb_list.append(order)
                consecutive_null_sb = 0
                found_anything = True
                if order % 2 == 1:
                    consecutive_null_odd = 0
            else:
                print "I could not find sideband with order", order
                last_sb = last_sb + THz_freq
                if found_anything == True:
                    consecutive_null_sb += 1
                if found_anything == True and order % 2 == 1:
                    consecutive_null_odd += 1
            if consecutive_null_odd == 2 and no_more_odds == False:
                print "I'm done looking for odd sidebands"
                no_more_odds = True
            if consecutive_null_sb == 2:
                print "I can't find any more sidebands"
                break
        
        print "I found these sidebands:", self.sb_list
        self.sb_guess = np.array([np.asarray(sb_freq_guess), np.asarray(sb_amp_guess), np.asarray(sb_error_estimate)]).T
    
    def fit_sidebands(self, plot=False):
        """
        This takes self.sb_guess and fits to each maxima to get the details of
        each sideband.  It's really ugly, but it works.  The error of the 
        sideband area is approximated from the data, not the curve fit.  All
        else is from the curve fit.
        
        Inputs:
        plot = if you want to plot the fits on the same plot as before
        
        Temporary stuff:
        sb_fits = holder of the fitting results until all spectra have been fit
        
        Attributes created:
        self.sb_results = the money maker
        """
        print "Trying to fit these"
        sb_fits = []
        for elem, num in enumerate(self.sb_index): # Have to do this because guess_sidebands doesn't out put data in the most optimized way
            data_temp = self.hsg_data[self.sb_index[elem] - 25:self.sb_index[elem] + 25, :]
            p0 = [self.sb_guess[elem, 0], self.sb_guess[elem, 1] / 30000, 0.0005, 1.0]
            #print "Let's fit this shit!"
            try:
                coeff, var_list = curve_fit(gauss, data_temp[:, 0], data_temp[:, 1], p0=p0)
                coeff[1] = abs(coeff[1])
                coeff[2] = abs(coeff[2]) # The linewidth shouldn't be negative
                #print "coeffs:", coeff
#                print coeff[1] / coeff[2], " vs. ", np.max(data/temp[:, 1]), "of", self.sb_list[elem]
#                print "sigma for {}: {}".format(self.sb_list[elem], coeff[2])
                if 10e-4 > coeff[2] > 10e-6:
                    sb_fits.append(np.hstack((self.sb_list[elem], coeff, np.sqrt(np.diag(var_list)))))
                    sb_fits[-1][6] = self.sb_guess[elem, 2] * sb_fits[-1][2] # the var_list wasn't approximating the error well enough, even when using sigma and absoluteSigma
                    # And had to scale by the area?
                if plot:
                    x_vals = np.linspace(data_temp[0, 0], data_temp[-1, 0], num=500)
                    plt.plot(x_vals, gauss(x_vals, *coeff), 
                             plt.gca().get_lines()[-1].get_color()+'--' # I don't really know. Mostly
                                                         # just looked around at what functions
                                                         # matplotlib has...
                             , linewidth = 3)
            except:
                print "I couldn't fit", elem
                self.sb_list[elem] = None
        sb_fits_temp = np.asarray(sb_fits)
        reorder = [0, 1, 5, 2, 6, 3, 7, 4, 8]
        try:
            sb_fits = sb_fits_temp[:, reorder]
        except:
            print "The file is:", self.fname
            print "\n!!!!!\nSHIT WENT WRONG\n!!!!!"
                
        # Going to label the appropriate row with the sideband
        self.sb_list = sorted(list([x for x in self.sb_list if x is not None]))
        sb_names = np.vstack(self.sb_list)
#        print "sb_names:", sb_names
#        print "sb_fits:", sb_fits
        sorter = np.argsort(sb_fits[:, 0])
        
        self.sb_results = np.array(sb_fits[sorter, :7])
#        print "sb_results:", self.sb_results

    def fit_sidebands_for_NIR_freq(self, sensitivity=2.5, plot=False):
        """
        This takes self.sb_guess and fits to each maxima to get the details of
        each sideband.
        
        sensitivity - the multiplier of the noise stdev required to save a peak
        plot - if you want to plot the fits on the same plot as before
        """
        self.sb_fits = []
        
        for elem in xrange(len(self.sb_index)):
            data_temp = self.hsg_data[self.sb_index[elem] - 50:self.sb_index[elem] + 50, :]
            p0 = [self.sb_guess[elem, 0], self.sb_guess[elem, 1], 0.0005, 1.0]
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
            print "\n!!!!!!!!!\nSome shit went wrong here!!!\n!!!!!!!!\n"
        
        # Going to label the appropriate row with the sideband
        sb_names = np.vstack(self.sb_list)
        self.sb_results = np.hstack((sb_names, self.sb_fits[:,:6]))

    
    def save_processing(self, file_name, folder_str, marker='', index=''):
        """
        This will save all of the self.hsg_data and the results from the 
        fitting of this individual file.
        
        Inputs:
        file_name = the beginning of the file name to be saved
        folder_str = the location of the folder where the file will be saved, 
                     will create the folder, if necessary.
        marker = I...I don't know what this was originally for
        index = used to keep these files from overwriting themselves when in a
                list
        
        Outputs:
        Two files, one that is self.hsg_data, the other is self.sb_results
        """
        try:
            os.mkdir(folder_str)
        except OSError, e:
            if e.errno == errno.EEXIST:
                pass
            else:
                raise
        
        self.sb_results[:, 5:] = self.sb_results[:, 5:] * 1000 # For meV linewidths
        
        spectra_fname = file_name + '_' + marker + '_' + str(index) + '.txt'
        fit_fname = file_name + '_' + marker + '_' + str(index) + '_fits.txt'
        self.save_name = spectra_fname
        
        self.parameters['addenda'] = self.addenda
        self.parameters['subtrahenda'] = self.subtrahenda
        try:
            parameter_str = json.dumps(self.parameters, sort_keys=True)
        except:
            print "Source: EMCCD_image.save_images\nJSON FAILED"
            print "Here is the dictionary that broke JSON:\n", self.parameters
            return

        origin_import_spec = '\nNIR frequency,Signal,Standard error\neV,arb. u.,arb. u.'
        spec_header = '#' + parameter_str + '\n#' + self.description[:-2] + origin_import_spec
        
        origin_import_fits = '\nSideband,Center energy,error,Sideband strength,error,Linewidth,error\norder,eV,,arb. u.,,meV,,arb. u.,\n' + marker
        fits_header = '#' + parameter_str + '\n#' + self.description[:-2] + origin_import_fits
        
        np.savetxt(os.path.join(folder_str, spectra_fname), self.hsg_data, delimiter=',',
                   header=spec_header, comments='', fmt='%f')
        np.savetxt(os.path.join(folder_str, fit_fname), self.sb_results, delimiter=',',
                   header=fits_header, comments='', fmt='%f')

        print "Save image.\nDirectory: {}".format(os.path.join(folder_str, spectra_fname))

    def stitch_spectra(self):
        """
        I really hope this is important later!
        """        
        raise NotImplementedError

class PMT(object):
    """
    This class will handle a bunch of files related to one HSG spectrum.  I'm 
    not sure currently how to handle the two HSG sources together.
    """
    def __init__(self, folder_path):
        """
        Initializes a SPEX spectrum.  It'll open a folder, and bring in all of
        the individual sidebands into this object. To work, all the individual
        sidebands need to be in that folder, not in separate folders
        
        attributes:
            self.parameters - dictionary of important experimental parameters
            self.description - string of the description of the file(s)
            self.sb_dict - keys are sideband order, values are PMT data arrays
            self.sb_list - sorted
        """
        print "This started"
        file_list = glob.glob(os.path.join(folder_path, '*.txt'))
        self.sb_dict = {}
        # in __main__.py, look for the method "genSaveHeader"
        #                 look for method "initSettings" to see what parameters are saved
        f = open(file_list[0],'rU')
        throw_away = f.readline() # Just need to get things down to the next line
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
        
        for sb_file in file_list:
            f = open(sb_file, 'rU')
            sb_num = int(f.readline().split(' ')[-1])
            print "Sideband number is", sb_num
            f.close()
            raw_temp = np.genfromtxt(sb_file, comments='#')#, delimiter=',') 
            # Hopefully saving will be comma delimited soon!
            frequencies = set(raw_temp[:, 0])
            fire_condition = np.mean(raw_temp[:, 2]) / 2 # Say FEL fired if the cavity dump signal is more than half the mean of the cavity dump signal
            print "The fire condition is", fire_condition
            temp = None
            for freq in frequencies:
                data_temp = np.array([])
                for raw_point in raw_temp:
                    if raw_point[0] == freq and raw_point[2] > fire_condition:
                        #print "I'm going to add this", raw_point[0], raw_point[3]
                        data_temp = np.hstack((data_temp, raw_point[3])) # I don't know why hstack works here and not concatenate
                print "The data temp is", data_temp
                print len(data_temp)
                try:
                    temp = np.vstack((temp, np.array([freq, np.mean(data_temp), np.std(data_temp) / np.sqrt(len(data_temp))])))
                except:
                    temp = np.array([freq, np.mean(data_temp), np.std(data_temp) / np.sqrt(len(data_temp))])
            temp[:, 0] = temp[:, 0] / 8065.6 # turn NIR freq into eV
            temp = temp[temp[:, 0].argsort()]
            self.sb_dict[sb_num] = np.array(temp)
        
        self.sb_list = sorted(self.sb_dict.keys())
        
    def fit_sidebands(self, plot=False):
        """
        This method will fit a gaussian to each of the sidebands provided in 
        the self.sb_dict and make a list just like in the EMCCD version.  It 
        will also use the standard error of the integral of the PMT peak as the
        error of the gaussian area instead of that element from the covariance
        matrix.  Seems more legit.  
        
        attributes:
        self.sb_results: the numpy array that contains all of the fit info just
                         like it does in the CCD class.
        """
        sb_fits = {}
        for sideband in self.sb_dict.items():
            print "Sideband number", sideband[0]
            index = np.argmax(sideband[1][:, 1])
            nir_frequency = sideband[1][index, 0]
            peak = sideband[1][index, 1]
            p0 = [nir_frequency, peak / 3000, 0.00006, 0.00001]
            print "p0:", p0
            try: 
                coeff, var_list = curve_fit(gauss, sideband[1][:, 0], sideband[1][:, 1], p0=p0)#, sigma=10*sideband[1][:, 2], absolute_sigma=True)
                coeff[1] = abs(coeff[1])
                coeff[2] = abs(coeff[2])
                print "coeffs:", coeff
                if np.sqrt(np.diag(var_list))[0] < 0.001: # The error on where the sideband is should be small
                    sb_fits[sideband[0]] = np.concatenate((np.array([sideband[0]]), coeff, np.sqrt(np.diag(var_list))))
                    #print "error then:", sb_fits[sideband[0]][6]
                    relative_error = np.sqrt(sum([x**2 for x in sideband[1][index - 1:index + 2, 2]])) / np.sum(sideband[1][index - 1:index + 2, 1])
                    print "relative error:", relative_error
                    sb_fits[sideband[0]][6] = coeff[1] * relative_error
                    #print "error now:", sb_fits[sideband[0]][6]                
                    if plot:
                        x_vals = np.linspace(np.amin(sideband[1][:, 0]), np.amax(sideband[1][:, 0]), num=50)
                        plt.plot(x_vals, gauss(x_vals, *coeff))
            except:
                print "God damn it, Leroy.\nYou couldn't fit this."
                sb_fits[sideband[0]] = None
            
        for result in sorted(sb_fits.keys()):
            try:
                self.sb_results = np.vstack((self.sb_results, sb_fits[result]))
            except:
                self.sb_results = np.array(sb_fits[result])
        
        self.sb_results = self.sb_results[:, [0, 1, 5, 2, 6, 3, 7, 4, 8]]
        self.sb_results = self.sb_results[:, :7]
        print "And the results, please:", self.sb_results 
    
    def save_processing(self, file_name, folder_str, marker='', index=''):
        """
        This will save all of the self.hsg_data and the results from the 
        fitting of this individual file.
        
        Inputs:
        file_name = the beginning of the file name to be saved
        folder_str = the location of the folder where the file will be saved, 
                     will create the folder, if necessary.
        marker = I...I don't know what this was originally for
        index = used to keep these files from overwriting themselves when in a
                list
        
        Outputs:
        Two files, one that is self.hsg_data, the other is self.sb_results
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
        self.save_name = spectra_fname
        
        self.parameters['addenda'] = self.addenda
        self.parameters['subtrahenda'] = self.subtrahenda
        try:
            parameter_str = json.dumps(self.parameters, sort_keys=True)
        except:
            print "Source: EMCCD_image.save_images\nJSON FAILED"
            print "Here is the dictionary that broke JSON:\n", self.parameters
            return

        origin_import_spec = '\nNIR frequency,Signal,Standard error\neV,arb. u.,arb. u.'
        spec_header = '#' + parameter_str + '\n#' + self.description[:-2] + origin_import_spec
        
        origin_import_fits = '\nCenter energy,error,Amplitude,error,Linewidth,error\neV,,arb. u.,,eV,,\n,,'# + marker
        fits_header = '#' + parameter_str + '\n#' + self.description[:-2] + origin_import_fits
        
        for sideband in sorted(self.sb_dict.keys()):
            try:
                complete = np.vstack((complete, self.sb_dict[sideband]))
            except:
                complete = np.array([self.sb_dict[sideband]])
        
        np.savetxt(os.path.join(folder_str, spectra_fname), complete, delimiter=',',
                   header=spec_header, comments='', fmt='%f')
        np.savetxt(os.path.join(folder_str, fit_fname), self.sb_results, delimiter=',',
                   header=fits_header, comments='', fmt='%f')

        print "Save image.\nDirectory: {}".format(os.path.join(folder_str, spectra_fname))
            
class Spectrum(object):
    """
    This class will collect the sb_results arrays from PMT and CCD arrays and
    combine them to make a complete spectrum.  Details are forthcoming.
    """
    def __init__(self, PMT_spectrum, CCD_spectrum):
        """
        This init will currently take a PMT object and a CCD object, grab the
        sb_results attributes, then overwrite the sidebands too low-order to 
        see accurately with the CCD with PMT results.  I have a feeling the CCD
        data will always be better, which is definitely true from a resolution
        standpoint.
        
        As it stands, we'll have to pair up the appropriate measurements 
        beforehand.  I'm not sure how to pass this a list of pmt spectra and
        ccd spectra and make it sort everything out.  I don't think that's how
        we'd want to do it.
        
        Hopefully we will have to add multiple CCD spectra to this class in the
        future.  I think the reason to make this a class is that this 
        functionality could be best implemented with methods.  
        """
        self.pmt_results = PMT_spectrum.sb_results
        self.ccd_results = CCD_spectrum.sb_results
        self.parameters = CCD_spectrum.parameters
        # Some functionality will be much easier with dictionaries, I think
        self.pmt_dict = {}
        self.full_dict = {}
        for sb in self.pmt_results:
            self.pmt_dict[sb[0]] = np.asarray(sb[1:])
        for sb in self.ccd_results:
            self.full_dict[sb[0]] = np.asarray(sb[1:])
        print "Full dictionary:", self.full_dict
        
        self.full_dict = stitch_dicts(self.full_dict, self.pmt_dict, bad_order=1)
        
    def plot_prep(self):
        """
        This method is called once the full_dict has been completed.  It will
        turn full_dict into a list like the results lists.
        """
        for sb in sorted(self.full_dict):
            temp = np.concatenate((np.array([sb]), self.full_dict[sb]))
            try:
                self.full_results = np.vstack((self.full_results, temp))
            except:
                self.full_results = np.array(temp)
        
        self.full_results[:, 5:] = self.full_results[:, 5:] * 1000
    
    def add_sidebands(self, CCD_spectrum):
        """
        This method will add more sidebands to the end of the current spectrum,
        to increase the effective dynamic range.
        """
        self.ccd2_results = CCD_spectrum.sb_results
        
        self.ccd2_dict = {}
        for sb in self.ccd2_results:
            self.ccd2_dict[sb[0]] = np.asarray(sb[1:])
        
        self.full_dict = stitch_dicts(self.ccd2_dict, self.full_dict, bad_order=1)      
        
    def save_processing(self, file_name, folder_str, marker='', index=''):
        """
        This will save all of the self.hsg_data and the results from the 
        fitting of this individual file.
        
        Inputs:
        file_name = the beginning of the file name to be saved
        folder_str = the location of the folder where the file will be saved, 
                     will create the folder, if necessary.
        marker = I...I don't know what this was originally for
        index = used to keep these files from overwriting themselves when in a
                list
        
        Outputs:
        Two files, one that is self.hsg_data, the other is self.sb_results
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
        self.save_name = spectra_fname
        '''
        self.parameters['addenda'] = self.addenda
        self.parameters['subtrahenda'] = self.subtrahenda
        '''
        try:
            parameter_str = json.dumps(self.parameters, sort_keys=True)
        except:
            print "Source: EMCCD_image.save_images\nJSON FAILED"
            print "Here is the dictionary that broke JSON:\n", self.parameters
            return

        origin_import_fits = '\nSideband,Center energy,error,Sideband strength,error,Linewidth,error\norder,eV,,arb. u.,,meV,,\n,,' + marker
        fits_header = '#' + parameter_str + origin_import_fits
        
        np.savetxt(os.path.join(folder_str, fit_fname), self.full_results, delimiter=',',
                   header=fits_header, comments='', fmt='%f')

        print "Save image.\nDirectory: {}".format(os.path.join(folder_str, spectra_fname))

####################
# Fitting functions 
####################


"""
Hey nutsack

It doesn't seem to work, obviously. Sigma is seeded with 3e-5, which just so happens to be
nearly the exact value that all of the sigmas seem to settle around. If I had to guess,
the algorithm does something weird where it steps outside the allowed band, sees it can't go there,
thinks it fucked up hard, and doesn't move again. maybe the gradient is too large, so it assumes
that it can't move. Maybe lower the walls, or have them be a function of the deviation from allowed, 
instead of a sharp wall. Or see if you can look into the stepping algorithm to figure out
another way to stop it.

"""


def makeGauss(sbNum):
    # The width of the sidebands is vaguely a function of sb number
    # I fit a sample set without bounding and found the general trend
    # of width vs. sbnum and rough bounds. I assumed it was linear and
    # just picked bounds which fairly well fit everything and cut out
    # the clearly incorrect ones
    m = 2.6e-6
    b1 = 3e-5 + 4e-5
    b2 = 3e-5 - 4e-5
    def gaussian(x, *p):
        mu, A, sigma, y0 = p
        penalty = 0
#        if m*sbNum+b1 < sigma:
#            penalty += 1e3 * (m*sbNum+b1 - sigma)/sigma
#        elif sigma< m*sbNum+b2:
#            penalty += 1e3 * (sigma - m*sbNum+b2)/sigma
#        if abs(y0)>1e2: penalty+=1e60
        return (A / sigma) * np.exp(-(x - mu)**2 / (2. * sigma**2)) + y0 + penalty
    return gaussian

def gauss(x, *p):
    mu, A, sigma, y0 = p
    return (A / sigma) * np.exp(-(x - mu)**2 / (2. * sigma**2)) + y0
    
def lingauss(x, *p):
    mu, A, sigma, y0, m = p
    return (A / sigma) * np.exp(-(x - mu)**2 / (2. * sigma**2)) + y0 + m*x

def lorentzian(x, *p):
    mu, A, gamma, y0 = p
    return (A * gamma**2) / ((x - mu)**2 + gamma**2) + y0

def stitch_dicts(main, new, bad_order=5):
    """
    This function will make the Spectrum class more readable.
    """
    overlap = []
    for new_sb in sorted(new.keys()):
        if new_sb in main.keys():
            overlap.append(new_sb)
    print "overlap:", overlap
    # Cut out the likely-bad ones from the CCD data
    overlap = [x for x in overlap if x > (bad_order + 0.5)]
#    overlap = overlap[-2:]
        
    # Calculate the ratio of the signals
    ratio_list = []
    for sb in overlap:
        ratio_list.append(main[sb][2] / new[sb][2])
    print "ratio_list:", ratio_list
    ratio = np.mean(ratio_list)
    ratio_err = np.std(ratio_list) / np.sqrt(len(ratio_list))
    print "ratio:", ratio
    print "ratio error:", ratio_err
    print "New data:", new
        
    # Scale up PMT data and add it to the full_dict, and propagate errors
    for sb in sorted(new.keys()):
        old = new[sb][2]
#        print "sideband:", sb
#        print "data relative error:", new[sb][3] / old
#        print "ratio relative error:", ratio_err / ratio
        new[sb][2] = ratio * new[sb][2]
        new[sb][3] = new[sb][2] * np.sqrt((ratio_err / ratio)**2 + (new[sb][3] / old)**2)
        if sb not in overlap: # Because I want to get rid of the lowest order guys attenuated by SP filter
            main[sb] = new[sb]
    return main

####################
# Collection functions
####################

def sum_spectra(object_list):
    """
    This function will add all the things that should be added.  Obvs.  It will
    also calculate the standard error of the mean for every NIR frequency.  The
    standard error is the sum of the dark noise and the "shot" noise.

    object_list: A list of spectrum objects
    :param object_list:
    :return:
    """
    print "I'm trying!"
    
    good_list = []
    for index in xrange(len(object_list)):
        dark_var = 0
        num_images = 0
        try:
            temp = object_list.pop(0)
            stderr_holder = np.array(temp.hsg_data[:, 1]).reshape((1600, 1))
#            print "Standard error holder shape 1:", stderr_holder.shape
        except:
#            print "God damn it, Leroy"
            break
#        print "temp has series: {}.\ttemp has cl: {}.\ttemp has series: {}".format(temp.parameters['series'], temp.parameters['center_lambda'], temp.parameters['series'])
        for spec in list(object_list):
#            print "\tspec has series: {}.\tspec has cl: {}.\tspec has fn: {}".format(spec.parameters['series'], spec.parameters['center_lambda'], spec.fname[-16:-13])
#            print "I am trying to add", temp.parameters['FELP'], spec.parameters['FELP']
            if temp.parameters['series'] == spec.parameters['series']:
                if temp.parameters['center_lambda'] == spec.parameters['center_lambda']:
                    temp += spec
                    num_images += 1
                    stderr_holder = np.hstack((stderr_holder, temp.hsg_data[:, 1].reshape((1600,1))))
                    print "Individual dark_stdev:", spec.dark_stdev
                    dark_var += (spec.dark_stdev)**2
#                    print "Standard error holder shape 2:", stderr_holder.shape
#                    print "\t\tadded"
                    #print "I ADDED", temp.parameters['FELP'], spec.parameters['FELP']
                    object_list.remove(spec)
        spec_number = stderr_holder.shape[1]
        std_error = np.sqrt(np.var(stderr_holder, axis=1, dtype=np.float64) + dark_var) / np.sqrt(spec_number) # Checking some sigma stuff from curve_fit
        # This standard error is for every point.  I think it actually overestimates
        # the error at places with no signal because we add the dark variance
        # effectively twice.
        print "final dark_stdev:", np.sqrt(dark_var)
        temp.add_std_error(std_error)
        temp.image_normalize(num_images)
        good_list.append(temp)
    return good_list

def stitch_hsg(object_list):
    """
    This will take the list of processed spectra and return a new one that has
    all hsg spectra that were taken in pieces and have them all stitched together.
    """
    good_list = []
    for index in xrange(len(object_list)):
        try:
            temp = object_list.pop(0)
        except:
            break
        
        for test in list(object_list):
            if temp.parameters['spec_section'] + 1 = test.parameters['spec_section']:
                
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
    included_spectra = dict()
    param_array = None
    sb_included = []
    
    for spec in spectrum_list:

        sb_included = sorted(list(set(sb_included + spec.sb_list)))
        included_spectra[spec.fname.split('/')[-1]] = spec.parameters[param_name]
        # If these are from summed spectra, then only the the first file name
        # from that sum will show up here, which should be fine?
    print "full name:", spectrum_list[0].fname
    print "included names:", included_spectra
    print "sb_included:", sb_included
    
    
    for spec in spectrum_list:
        temp_dict = {}
        print "the sb_results:", spec.sb_results
        for index in xrange(len(spec.sb_results[:, 0])):
            print "my array slice:", spec.sb_results[index, :]
            temp_dict[int(round(spec.sb_results[index, 0]))] = np.array(spec.sb_results[index, :])
        print temp_dict
        
        for sb in sb_included:
            blank = np.zeros(7)
            blank[0] = float(sb)
            print "checking sideband order:", sb
            print "blank", blank
            if not temp_dict.has_key(sb):
                print "\nNeed to add sideband order:", sb
                temp_dict[sb] = blank
        try:
            spec_data = np.array([float(spec.parameters[param_name]), spec.dark_stdev])
        except:
            spec_data = np.array([float(spec.parameters[param_name][:2]), spec.dark_stdev])
        for key in sorted(temp_dict.keys()):
            print "I am going to hstack this:", temp_dict[key]
            spec_data = np.hstack((spec_data, temp_dict[key]))
            
        try:
            param_array = np.vstack((param_array, spec_data))
        except:
            param_array = np.array(spec_data)
        print "The shape of the param_array is:", param_array.shape
        #print "The param_array itself is:", param_array
    
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
    origin_import1 = param_name + ",dark_stdev"
    origin_import2 = unit + ",post shot norm"
    origin_import3 = ","
    for order in sb_included:
        origin_import1 += ",Sideband,Frequency,error,Sideband strength,error,Linewidth,error"
        origin_import2 += ",order,eV,,arb. u.,,eV,"
        origin_import3 += ",,{0},,{0},,{0},".format(order)
    origin_total = origin_import1 + "\n" + origin_import2 + "\n" + origin_import3
    header = '#' + included_spectra_str + '\n' + origin_total
    #print "Spec header: ", spec_header
    print "the param_array is:", param_array
    np.savetxt(os.path.join(folder_str, file_name), param_array, delimiter=',', 
               header=header, comments='', fmt='%f')

    print "Saved the file.\nDirectory: {}".format(os.path.join(folder_str, file_name))
    
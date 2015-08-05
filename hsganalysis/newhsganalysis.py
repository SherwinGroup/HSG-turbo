# -*- coding: utf-8 -*-
"""
Created on Fri Jun 19 13:16:25 2015

@author: hbanks

Brevity required, prurience preferred
"""

from __future__ import division
import os
import glob
import errno
import copy
import json
import numpy as np
from scipy.optimize import curve_fit
import scipy.interpolate as spi
import scipy.optimize as spo
import scipy.fftpack as fft
import matplotlib.pyplot as plt

####################
# Objects 
####################

class CCD(object):
    def __init__(self, fname):
        """
        This will read the appropriate file and make a basic CCD object.  Fancier
        things will be handled with the sub classes.

        input:
        fname = file name where the data is saved

        creates:
        self.parameters = JSON dictionary holding all of the information from the
                          data file.
        self.description = string that is the text box from data taking GUI
        self.raw_data = raw data output by measurement software, wavelength vs.
                        data.  It may not entirely be numbers?
        self.ccd_data = semi-processed 1600 x 2 array of photon energy vs. data
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
        try:
            self.parameters["spec_step"] = int(self.parameters["spec_step"])
        except KeyError:
            pass
        self.raw_data = np.flipud(np.genfromtxt(fname, comments='#', delimiter=','))
        # I used flipup so that the x-axis is an increasing function of frequency
        
        self.ccd_data = np.array(self.raw_data[:1600,:]) # By slicing up to 1600,
                                                         # we cut out the text 
                                                         # header
        self.ccd_data[:, 0] = 1239.84 / self.ccd_data[:, 0]

    def __str__(self):
        return self.description

class Photoluminescence(CCD):
    def __init__(self, fname):
        """
        This object handles PL-type data.

        creates:
        self.proc_data = self.ccd_data divided by the exposure time
                         units: PL counts / second
        """
        super(Photoluminescence, self).__init__(fname)

        self.proc_data = np.array(self.ccd_data) # Does this work the way I want it to?
        self.proc_data[:, 1] = self.proc_data[:, 1] / self.parameters['exposure']
#     def save_processing(self):
        
class Absorbance(CCD):
    def __init__(self, fname):
        """
        There are several ways Absorbance data can be loaded
        You could try to load the abs data output from data collection directly,
        which has the wavelength, raw, blank and actual absorbance data itself

        Alternatively, you could want to load the raw transmission/reference
        data, ignoring (or maybe not even having) the abs calculated
        from the data collection software.

        I'm not sure the best way to input the files you want...
        Maybe either the abs data itself, or a list of the [ref, trans]?
        :param fname: either an absorbance filename, or a length 2 list of filenames
        :return:
        """

        if "abs_" in fname:
            CCD.__init__(self, fname)
            # Separate into the separate data sets
            #   The raw counts of the reference data
            self.ref_data = np.array(self.ccd_data[:, [0, 1]])
            #   Raw counts of the sample
            self.raw_data = np.array(self.ccd_data[:, [0, 2]])
            #   The calculated absorbance data (-log10(raw/ref))
            self.proc_data = np.array(self.ccd_data[:, [0, 3]])
        else:
            # Should be here if you pass the reference/trans filenames
            try:
                super(Absorbance, self).__init__(fname[0])
                self.ref_data = np.array(self.ccd_data)

                super(Absorbance, self).__init__(fname[1])
                self.raw_data = np.array(self.ccd_data)
            except ValueError:
                # ValueError gets thrown when importing older data
                # which had more headers than data columns. Enforce
                # only loading first two columns to avoid

                self.ref_data = np.flipud(np.genfromtxt(fname[0], comments='#',
                                                        delimiter=',', usecols=(0, 1)))

                self.ref_data = np.array(self.ref_data[:1600,:])
                self.ref_data[:, 0] = 1239.84 / self.ref_data[:, 0]

                self.raw_data = np.flipud(np.genfromtxt(fname[1], comments='#',
                                                        delimiter=',', usecols=(0, 1)))

                self.raw_data = np.array(self.raw_data[:1600,:])
                self.raw_data[:, 0] = 1239.84 / self.raw_data[:, 0]
            except Exception as e:
                print "Exception opening,", e
            self.proc_data = np.empty_like(self.ref_data)
            self.proc_data[:,0] = self.ref_data[:,0]
            self.proc_data[:,1] = np.log10(self.raw_data[:,1]/self.ref_data[:,1])

    def abs_per_QW(self, qw_number):
        """
        This method turns the absorption to the absorbance per quantum well.  Is
        that how this data should be reported?

        Also, I'm not sure if columns 1 and 2 are correct.
        """
        temp_abs = -np.log(self.proc_data[:, 1] / self.proc_data[:, 2]) / qw_number
        self.proc_data = np.hstack((self.proc_data, temp_abs))

    def save_processing(self, file_name, folder_str, marker='', index=''):
        try:
            os.mkdir(folder_str)
        except OSError, e:
            if e.errno == errno.EEXIST:
                pass
            else:
                raise
        
        spectra_fname = file_name + '_' + marker + '_' + str(index) + '.txt'
        self.save_name = spectra_fname

        try:
            parameter_str = json.dumps(self.parameters, sort_keys=True)
        except:
            print "Source: EMCCD_image.save_images\nJSON FAILED"
            print "Here is the dictionary that broke JSON:\n", self.parameters
            return

        origin_import_spec = '\nNIR frequency,Signal,Standard error\neV,arb. u.,arb. u.'
        spec_header = '#' + parameter_str + '\n#' + self.description[:-2] + origin_import_spec
        
        np.savetxt(os.path.join(folder_str, spectra_fname), self.proc_data, delimiter=',',
                   header=spec_header, comments='', fmt='%f')

        print "Save image.\nDirectory: {}".format(os.path.join(folder_str, spectra_fname))

class HighSidebandCCD(CCD):
    def __init__(self, fname):
        """
        This will read the appropriate file.  The header needs to be fixed to
        reflect the changes to the output header from the Andor file.  Because
        another helper file will do the cleaning and background subtraction,
        those are no longer part of this init.  This also turns all wavelengths
        from nm (NIR ones) or cm-1 (THz ones) into eV.
        
        Input:
        fname = file name of the hsg spectrum from CCD superclass
        
        Internal:
        self.fname = the filename
        self.parameters = string with all the relevant experimental perameters
        self.description = the description we added to the file as the data
                           was being taken
        self.proc_data = processed data that has gone is frequency vs counts/pulse
        self.dark_stdev = the standard deviation of the dark noise taken from
                          the background file when data was taken (I hope).  It
                          is divided by the number of FEL pulses
        self.addenda = the list of things that have been added to the file, in
                       form of [constant, *spectra_added]
        self.subtrahenda = the list of spectra that have been subtracted from
                           the file.  Constant subtraction is dealt with with
                           self.addenda
        """
        super(HighSidebandCCD, self).__init__(fname)

        self.proc_data = np.array(self.ccd_data) # Does this work the way I want it to?
        self.proc_data[:, 1] = self.proc_data[:, 1] / self.parameters['fel_pulses']

        self.parameters["nir_freq"] = 1239.84 / float(self.parameters["nir_lambda"])
        self.parameters["thz_freq"] = 0.000123984 * float(self.parameters["fel_lambda"])
        self.parameters["nir_power"] = float(self.parameters["nir_power"])
        self.parameters["thz_power"] = float(self.parameters["fel_power"])

        # Need to do this next line AFTER adding all the spectra together

        
        self.dark_stdev = self.parameters["background_darkcount_std"] / self.parameters['fel_pulses'] 
        # I think this error stuff is working now

        
        self.addenda = self.parameters['addenda']
        self.subtrahenda = self.parameters['subtrahenda']
        #self.sb_results = None

    def __add__(self, other):
        """
        Add together the image data from self.proc_data, or add a constant to 
        that np.array.  It will then combine the addenda and subtrahenda lists,
        as well as add the fel_pulses together.

        Input:
        self = CCD-like object
        other = int, float or CCD object

        Internal:
        ret.proc_data = the self.proc_data + other(.proc_data)
        ret.addenda = combination of two input addenda lists
        """
        ret = copy.deepcopy(self)

        # Add a constant offset to the data
        if type(other) in (int, float):
            ret.proc_data[:, 1] = self.proc_data[:, 1] + other 
            ret.addenda[0] = ret.addenda[0] + other
        
        # or add the data of two hsg_spectra together
        else:
            if np.isclose(ret.parameters['center_lambda'], other.parameters['center_lambda']):
                ret.proc_data[:, 1] = self.proc_data[:, 1] + other.proc_data[:, 1] 
                ret.addenda[0] = ret.addenda[0] + other.addenda[0]
                ret.addenda.extend(other.addenda[1:])
                ret.subtrahenda.extend(other.subtrahenda)
                ret.parameters['fel_pulses'] += other.parameters['fel_pulses']
            else:
                raise Exception('Source: Spectrum.__add__:\nThese are not from the same grating settings')
        return ret

    def __sub__(self, other):
        """
        This subtracts constants or other data sets between self.proc_data.  I 
        think it even keeps track of what data sets are in the file and how 
        they got there.

        See how __add__ works for more information.
        """
        ret = copy.deepcopy(self)
        
        # Subtract a constant offset to the data
        if type(other) in (int, float):
            ret.proc_data[:,1] = self.proc_data[:,1] - other # Need to choose a name
            ret.addenda[0] = ret.addenda[0] - other
            
        # Subtract the data of two hsg_spectra from each other
        else:
            if np.isclose(ret.proc_data[0,0], other.proc_data[0,0]):
                ret.proc_data[:,1] = self.proc_data[:,1] - other.proc_data[:,1]
                ret.subtrahenda.extend(other.addenda[1:])
                ret.addenda.extend(other.subtrahenda)
            else:
                raise Exception('Source: Spectrum.__sub__:\nThese are not from the same grating settings')
        return ret

    def add_std_error(self, std_errors):
        """
        This adds a numpy array of standard errors to the self.proc_data array.
        It's actually just an append_column method, but there's no reason for 
        that (I think)

        Input:
        std_errors = 1600x1 array that includes the standard error for every 
                     individual point

        Internal:
        self.proc_data = 1600x3 array that has [photon energy, signal, error]
        """
        try:
            self.proc_data = np.hstack((self.proc_data, std_errors.reshape((1600, 1))))
        except:
            print "Spectrum.add_std_error fucked up.  What's wrong?"
            print std_errors.shape
    
    def calc_approx_sb_order(self, test_nir_freq):
        """
        This simple method will simply return a float approximating the order
        of the frequency input.  We need this because the CCD wavelength 
        calibration is not even close to perfect.  And it shifts by half a nm 
        sometimes.

        Input: 
        test_nir_freq = NIR frequency that we're guessing the order for

        Returns:
        approx_order = the approximate sideband order of the test frequency
        """
        nir_freq = self.parameters['nir_freq']
        thz_freq = self.parameters['thz_freq']
        approx_order = (test_nir_freq - nir_freq) / thz_freq
        return approx_order
    
    def image_normalize(self, num_images):
        """
        This method will divide the proc_data by the number of images that were
        used in the hsg_sum_spectra function.  The proc_data array contains CCD
        signal per FEL pulse already.  

        Input:
        num_images = number of images proc_data came from

        Internal:
        self.proc_data = normalizes the data and error columns to be signal/pulse
                         again
        """
        self.proc_data[:, 1] = self.proc_data[:, 1] / num_images
        self.proc_data[:, 2] = self.proc_data[:, 2] / num_images
    
    def guess_sidebands(self, cutoff=4, verbose=False):
        """
        Finds the locations of all the sidebands in the proc_data array to be 
        able to seed the fitting method.  This works by finding the maximum data
        value in the array and guessing what sideband it is.  It creates an array
        that includes this information.  It will then step down, initially by one
        THz frequency, then by twos after it hasn't found any odd ones.  It then
        goes up from the max and finds everything above in much the same way.

        Input:
        cutoff = signal-to-noise threshold to count a sideband candidate.

        Internal:
        self.sb_list = List of all of the orders the method found
        self.sb_index = index of all of the peaks of the sidebands
        self.sb_guess = three-part list including the frequency, amplitude and
                        error guesses for each sideband
        """
        x_axis = np.array(self.proc_data[:, 0])
        y_axis = np.array(self.proc_data[:, 1])
        error = np.array(self.proc_data[:, 2])
        
        min_sb = int(self.calc_approx_sb_order(x_axis[0])) + 1
        max_sb = int(self.calc_approx_sb_order(x_axis[-1]))
        
        nir_freq = self.parameters["nir_freq"]
        thz_freq = self.parameters["thz_freq"]
        
        # Find max strength sideband and it's order
        global_max = np.argmax(y_axis)
        order_init = int(round(self.calc_approx_sb_order(x_axis[global_max])))

        check_y = y_axis[global_max - 15:global_max + 15]
        check_max_area = np.sum(y_axis[global_max - 1:global_max + 2])
        check_ave = np.mean(check_y)
        check_stdev = np.std(check_y)

        check_ratio = (check_max_area - 3 * check_ave) / check_stdev
        
        
        if verbose:
            print "Global max checking:", check_y
            print "\ncheck_max_area is", check_max_area
            print "check_ave is", check_ave
            print "check_stdev is", check_stdev
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
        if verbose:
            print "Now I'll look for sidebands with frequency less than", sb_freq_guess
        last_sb = sb_freq_guess[0]
        index_guess = global_max
        consecutive_null_sb = 0
        consecutive_null_odd = 0
        no_more_odds = False
        break_condition = False
        for order in xrange(order_init - 1, min_sb - 1, -1):
            if no_more_odds == True and order % 2 == 1:
                last_sb = last_sb - thz_freq
                if verbose:
                    print "I skipped", order
                continue
            lo_freq_bound = last_sb - thz_freq * (1 + 0.22) # Not sure what to do about these
            hi_freq_bound = last_sb - thz_freq * (1 - 0.22)

            start_index = False
            end_index = False
            if verbose:
                print "\nSideband", order, "\n"
                print "The low frequency bound is", lo_freq_bound
                print "The high frequency bound is", hi_freq_bound
            for i in xrange(index_guess, 0, -1):
                if end_index == False and i == 1:
                    break_condition = True
                    break
                if end_index == False and x_axis[i] < hi_freq_bound:
                    #print "end_index is", i
                    end_index = i
                elif i == 1:
                    start_index = 0
                    #print "hit end of data, start_index is 0"
                elif start_index == False and x_axis[i] < lo_freq_bound:
                    start_index = i
                    #print "start_index is", i
                    index_guess = i
                    break
            
            if break_condition:
                break
            check_y = y_axis[start_index:end_index]

            #check_max = check_y.max()
            #print "check_max is", check_max
            check_max_index = np.argmax(check_y) # This assumes that two floats won't be identical
            check_max_area = np.sum(check_y[check_max_index - 1:check_max_index + 2])
            check_ave = np.mean(check_y)
            check_stdev = np.std(check_y)
            check_ratio = (check_max_area - 3 * check_ave) / check_stdev

            if verbose:
                print "check_y is", check_y
                print "\ncheck_max_area is", check_max_area
                print "check_ave is", check_ave
                print "check_stdev is", check_stdev
                print "check_ratio is", check_ratio
            
            if check_ratio > cutoff:
                found_index = check_max_index + start_index
                self.sb_index.append(found_index)
                last_sb = x_axis[found_index]
                
                if verbose:
                    print "I just found", last_sb
                
                sb_freq_guess.append(x_axis[found_index])
                sb_amp_guess.append(check_max_area - 3 * check_ave)
                error_est = np.sqrt(sum([i**2 for i in error[found_index - 1:found_index + 2]])) / (check_max_area - 3 * check_ave)
                if verbose:
                    print "My error estimate is:", error_est
                sb_error_est.append(error_est)
                self.sb_list.append(order)
                consecutive_null_sb = 0
                if order % 2 == 1:
                    consecutive_null_odd = 0
            else:
                #print "I could not find sideband with order", order
                last_sb = last_sb - thz_freq
                consecutive_null_sb += 1
                if order % 2 == 1:
                    consecutive_null_odd += 1
            if consecutive_null_odd == 1 and no_more_odds == False:
                #print "I'm done looking for odd sidebands"
                no_more_odds = True
            if consecutive_null_sb == 2:
                #print "I can't find any more sidebands"
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
            window_size = 0.22 + 0.0006 * last_sb
            lo_freq_bound = last_sb + thz_freq * (1 - window_size) # Not sure what to do about these
            hi_freq_bound = last_sb + thz_freq * (1 + window_size)
            start_index = False
            end_index = False

            if verbose:
                print "\nSideband", order, "\n"            
            for i in xrange(index_guess, 1600):
                if start_index == False and i == 1599:
                    #print "I'm all out of space, captain!"
                    break_condition = True
                    break
                elif start_index == False and x_axis[i] > lo_freq_bound:
                    #print "start_index is", i
                    start_index = i
                elif i == 1599:
                    end_index = 1599
                    #print "hit end of data, end_index is 1599"
                elif end_index == False and x_axis[i] > hi_freq_bound:
                    end_index = i
                    #print "end_index is", i
                    index_guess = i
                    break
            if break_condition:
                break
            check_y = y_axis[start_index:end_index]

            #check_max = check_y.max()
            #print "check_max is", check_max
            check_max_index = np.argmax(check_y) # This assumes that two floats won't be identical
            check_max_area = np.sum(check_y[check_max_index - 1:check_max_index + 2])
            check_ave = np.mean(check_y)
            check_stdev = np.std(check_y)
            check_ratio = (check_max_area - 3 * check_ave) / check_stdev

            if verbose:
                print "check_y is", check_y
                print "\ncheck_max_area is", check_max_area
                print "check_ave is", check_ave
                print "check_stdev is", check_stdev
                print "check_ratio is", check_ratio
            
            if check_ratio > cutoff:
                found_index = check_max_index + start_index
                self.sb_index.append(found_index)
                last_sb = x_axis[found_index]
                
                if verbose:
                    print "I just found", last_sb
                
                sb_freq_guess.append(x_axis[found_index])
                sb_amp_guess.append(check_max_area - 3 * check_ave)
                error_est = np.sqrt(sum([i**2 for i in error[found_index - 1:found_index + 2]])) / (check_max_area - 3 * check_ave)
                if verbose:
                    print "My error estimate is:", error_est
                sb_error_est.append(error_est)
                self.sb_list.append(order)
                consecutive_null_sb = 0
                if order % 2 == 1:
                    consecutive_null_odd = 0
            else:
                #print "I could not find sideband with order", order
                last_sb = last_sb + thz_freq
                consecutive_null_sb += 1
                if order % 2 == 1:
                    consecutive_null_odd += 1
            if consecutive_null_odd == 1 and no_more_odds == False:
                #print "I'm done looking for odd sidebands"
                no_more_odds = True
            if consecutive_null_sb == 2:
                #print "I can't find any more sidebands"
                break  
        
        if verbose:
            print "I found these sidebands:", self.sb_list
        self.sb_guess = np.array([np.asarray(sb_freq_guess), np.asarray(sb_amp_guess), np.asarray(sb_error_est)]).T

    def fit_sidebands(self, plot=False, verbose=False):
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
        self.full_dict = a dictionary similar to sb_results, but now the keys 
                         are the sideband orders
        """
        #print "Trying to fit these"
        sb_fits = []
        for elem, num in enumerate(self.sb_index): # Have to do this because guess_sidebands 
                                                   # doesn't out put data in the most optimized way
            data_temp = self.proc_data[self.sb_index[elem] - 25:self.sb_index[elem] + 25, :]
            width_guess = 0.0001 + 0.000001*num # so the width guess gets wider as order goes up
            p0 = [self.sb_guess[elem, 0], self.sb_guess[elem, 1] * width_guess, width_guess, 1.0]
            #print "Let's fit this shit!"
            try:
                coeff, var_list = curve_fit(gauss, data_temp[:, 0], data_temp[:, 1], p0=p0)
                coeff[1] = abs(coeff[1])
                coeff[2] = abs(coeff[2]) # The linewidth shouldn't be negative
                if verbose:
                    print "coeffs:", coeff
                    print "sigma for {}: {}".format(self.sb_list[elem], coeff[2])
                if 10e-4 > coeff[2] > 10e-6:
                    sb_fits.append(np.hstack((self.sb_list[elem], coeff, np.sqrt(np.diag(var_list)))))
                    sb_fits[-1][6] = self.sb_guess[elem, 2] * sb_fits[-1][2] # the var_list wasn't approximating the error well enough, even when using sigma and absoluteSigma
                    # And had to scale by the area?
                if plot:
                    linewidth = 5
                    x_vals = np.linspace(data_temp[0, 0], data_temp[-1, 0], num=500)
                    if elem!=0:
                        try:
                            plt.plot(x_vals, gauss(x_vals, *coeff), 
                                     plt.gca().get_lines()[-1].get_color()+'--' # I don't really know. Mostly
                                                                 # just looked around at what functions
                                                                 # matplotlib has...
                                     , linewidth = linewidth)
                        except: # to prevent weird mac issues with the matplotlib things?
                            plt.plot(x_vals, gauss(x_vals, *coeff), '--', linewidth=linewidth)
                                         
                    else:
                        plt.plot(x_vals, gauss(x_vals, *coeff), '--', linewidth=linewidth)
            except:
                print "I couldn't fit", elem
                print "In file", self.fname
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

        sorter = np.argsort(sb_fits[:, 0])
        
        self.sb_results = np.array(sb_fits[sorter, :7])
        if verbose:
            print "sb_names:", sb_names
            print "sb_fits:", sb_fits
            print "sb_results:", self.sb_results
        self.full_dict = {}
        for sb in self.sb_results:
            self.full_dict[sb[0]] = np.asarray(sb[1:])
        
    def save_processing(self, file_name, folder_str, marker='', index=''):
        """
        This will save all of the self.proc_data and the results from the 
        fitting of this individual file.
        
        Inputs:
        file_name = the beginning of the file name to be saved
        folder_str = the location of the folder where the file will be saved, 
                     will create the folder, if necessary.
        marker = I...I don't know what this was originally for
        index = used to keep these files from overwriting themselves when in a
                list
        
        Outputs:
        Two files, one that is self.proc_data, the other is self.sb_results
        """
        try:
            os.mkdir(folder_str)
        except OSError, e:
            if e.errno == errno.EEXIST:
                pass
            else:
                raise
        self.sb_results[:, 5:7] = self.sb_results[:, 5:7] * 1000 # For meV linewidths
        area = np.array([self.sb_results[:, 3] / self.sb_results[:, 5]])
        print "sb_results", self.sb_results.shape
        print "area", area.shape
        save_results = np.hstack((self.sb_results, area.T))
        
        
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
        
        origin_import_fits = '\nSideband,Center energy,error,Sideband strength,error,Linewidth,error,Area\norder,eV,,arb. u.,,meV,,arb. u.\n' + marker
        fits_header = '#' + parameter_str + '\n#' + self.description[:-2] + origin_import_fits
        
        np.savetxt(os.path.join(folder_str, spectra_fname), self.proc_data, delimiter=',',
                   header=spec_header, comments='', fmt='%f')
        np.savetxt(os.path.join(folder_str, fit_fname), save_results, delimiter=',',
                   header=fits_header, comments='', fmt='%f')

        print "Save image.\nDirectory: {}".format(os.path.join(folder_str, spectra_fname))

class PMT(object):
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
            fire_condition = np.mean(raw_temp[:, 2]) / 2 # Say FEL fired if the 
                                                         # cavity dump signal is
                                                         # more than half the mean 
                                                         # of the cavity dump signal
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
    
class HighSidebandPMT(PMT):
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
            fire_condition = np.mean(raw_temp[:, 2]) / 2 # Say FEL fired if the 
                                                         # cavity dump signal is
                                                         # more than half the mean 
                                                         #of the cavity dump signal
            print "The fire condition is", fire_condition
            temp = None
            for freq in frequencies:
                data_temp = np.array([])
                for raw_point in raw_temp:
                    if raw_point[0] == freq and raw_point[2] > fire_condition:
                        #print "I'm going to add this", raw_point[0], raw_point[3]
                        data_temp = np.hstack((data_temp, raw_point[3])) 
                                    # I don't know why hstack works here and not concatenate
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
    
    def laser_line(self):
        """
        This method is designed to scale everything in the PMT to the conversion
        efficiency based on our measurement of the laser line with a fixed 
        attenuation.
        """
        if 0 not in self.sb_list:
            return
        else:
            return

    def fit_sidebands(self, plot=False, verbose=False):
        """
        This method will fit a gaussian to each of the sidebands provided in 
        the self.sb_dict and make a list just like in the EMCCD version.  It 
        will also use the standard error of the integral of the PMT peak as the
        error of the gaussian area instead of that element from the covariance
        matrix.  Seems more legit.  
        
        attributes:
        self.sb_results: the numpy array that contains all of the fit info just
                         like it does in the CCD class.
        self.full_dict = A dictionary version of self.sb_results
        """
        sb_fits = {}
        for sideband in self.sb_dict.items():
            if verbose:
                print "Sideband number", sideband[0]
            index = np.argmax(sideband[1][:, 1])
            nir_frequency = sideband[1][index, 0]
            peak = sideband[1][index, 1]
            width_guess = 0.0003 # Yep, another magic number
            p0 = [nir_frequency, peak * width_guess, width_guess, 0.00001]
            if verbose:
                print "p0:", p0
            try: 
                coeff, var_list = curve_fit(gauss, sideband[1][:, 0], sideband[1][:, 1], p0=p0)
                coeff[1] = abs(coeff[1])
                coeff[2] = abs(coeff[2])
                if verbose:
                    print "coeffs:", coeff
                if np.sqrt(np.diag(var_list))[0] < 0.001: # The error on where the sideband is should be small
                    sb_fits[sideband[0]] = np.concatenate((np.array([sideband[0]]), coeff, np.sqrt(np.diag(var_list))))
                    #print "error then:", sb_fits[sideband[0]][6]
                    relative_error = np.sqrt(sum([x**2 for x in sideband[1][index - 1:index + 2, 2]])) / np.sum(sideband[1][index - 1:index + 2, 1])
                    if verbose:
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
        if verbose:
            print "And the results, please:", self.sb_results

        self.full_dict = {}
        for sb in self.sb_results:
            self.full_dict[sb[0]] = np.asarray(sb[1:])
    
    def save_processing(self, file_name, folder_str, marker='', index=''):
        """
        This will save all of the self.proc_data and the results from the 
        fitting of this individual file.
        
        Inputs:
        file_name = the beginning of the file name to be saved
        folder_str = the location of the folder where the file will be saved, 
                     will create the folder, if necessary.
        marker = I...I don't know what this was originally for
        index = used to keep these files from overwriting themselves when in a
                list
        
        Outputs:
        Two files, one that is self.proc_data, the other is self.sb_results
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

class FullSpectrum(object):
    def __init__(self):
        pass

class FullAbsorbance(FullSpectrum):
    """
    I'm imagining this will sew up absorption spectra, but I'm not at all sure
    how to do that at the moment.
    """
    def __init__(self):
        pass

class FullHighSideband(FullSpectrum):
    """
    I'm imagining this class is created with a base CCD file, then gobbles up
    other spectra that belong with it, then grabs the PMT object to normalize
    everything, assuming that PMT object exists.
    """
    def __init__(self, initial_CCD_piece):
        """
        Initialize a full HSG spectrum.  Starts with a single CCD image, then
        adds more on to itself using stitch_hsg_dicts.
        """
        self.fname = initial_CCD_piece.fname
        self.description = initial_CCD_piece.description
        self.sb_results = initial_CCD_piece.sb_results
        self.parameters = initial_CCD_piece.parameters
        self.parameters['files_here'] = [initial_CCD_piece.fname.split('/')[-1]]
        self.full_dict = {}
        for sb in self.sb_results:
            self.full_dict[sb[0]] = np.asarray(sb[1:])

    def add_CCD(self, ccd_object):
        """
        This method will be called by the stitch_hsg_results function to add another
        CCD image to the spectrum.
        """
        if self.parameters["gain"] == ccd_object.parameters["gain"]:
            calc = False
        else:
            print "!!!"
            print "What happened"
            print "!!!"

            calc = True
        self.full_dict = stitch_hsg_dicts(self.full_dict, ccd_object.full_dict, need_ratio=calc)
        self.parameters['files_here'].append(ccd_object.fname.split('/')[-1])

    def add_PMT(self, pmt_object):
        """
        This method will be called by the stitch_hsg_results function to add the PMT
        data to the spectrum.
        """
        self.full_dict = stitch_hsg_dicts(pmt_object.full_dict, self.full_dict)
        self.parameters['files_here'].append(pmt_object.fname.split('/')[-1])

    def make_results_array(self):
        """
        The idea behind this method is to create the sb_results array from the 
        finished full_dict dictionary.
        """
        self.sb_results = None
        #print "I'm making the results array:", sorted(self.full_dict.keys())
        for sb in sorted(self.full_dict.keys()):
            #print "Going to add this", sb
            try:
                self.sb_results = np.vstack((self.sb_results, np.hstack((sb, self.full_dict[sb]))))
            except ValueError:
                #print "It didn't exist yet!"
                self.sb_results = np.hstack((sb, self.full_dict[sb]))
        #print "and I made this array:", self.sb_results[:, 0]

    def save_processing(self, file_name, folder_str, marker='', index=''):
        """
        This will save all of the self.proc_data and the results from the 
        fitting of this individual file.
        
        Inputs:
        file_name = the beginning of the file name to be saved
        folder_str = the location of the folder where the file will be saved, 
                     will create the folder, if necessary.
        marker = I...I don't know what this was originally for
        index = used to keep these files from overwriting themselves when in a
                list
        
        Outputs:
        Two files, one that is self.proc_data, the other is self.sb_results
        """
        try:
            os.mkdir(folder_str)
        except OSError, e:
            if e.errno == errno.EEXIST:
                pass
            else:
                raise
        
        area = np.array([self.sb_results[:, 3] / self.sb_results[:, 5]]) * 1000
        print "sb_results", self.sb_results.shape
        print "area", area.shape
        save_results = np.hstack((self.sb_results, area.T))
        
        
        #spectra_fname = file_name + '_' + marker + '_' + str(index) + '.txt'
        fit_fname = file_name + '_' + marker + '_' + str(index) + '_full.txt'
        #self.save_name = spectra_fname
        
        #self.parameters['addenda'] = self.addenda
        #self.parameters['subtrahenda'] = self.subtrahenda
        try:
            parameter_str = json.dumps(self.parameters, sort_keys=True)
        except:
            print "Source: EMCCD_image.save_images\nJSON FAILED"
            print "Here is the dictionary that broke JSON:\n", self.parameters
            return

        #origin_import_spec = '\nNIR frequency,Signal,Standard error\neV,arb. u.,arb. u.'
        #spec_header = '#' + parameter_str + '\n#' + self.description[:-2] + origin_import_spec
        
        origin_import_fits = '\nSideband,Center energy,error,Sideband strength,error,Linewidth,error,Area\norder,eV,,arb. u.,,meV,,arb. u.\n' + marker
        fits_header = '#' + parameter_str + '\n#' + self.description[:-2] + origin_import_fits
        
        #np.savetxt(os.path.join(folder_str, spectra_fname), self.proc_data, delimiter=',',
        #           header=spec_header, comments='', fmt='%f')
        np.savetxt(os.path.join(folder_str, fit_fname), save_results, delimiter=',',
                   header=fits_header, comments='', fmt='%f')

        print "Save image.\nDirectory: {}".format(os.path.join(folder_str, fit_fname))

####################
# Fitting functions 
####################

def gauss(x, *p):
    mu, A, sigma, y0 = p
    return (A / sigma) * np.exp(-(x - mu)**2 / (2. * sigma**2)) + y0
    
def lingauss(x, *p):
    mu, A, sigma, y0, m = p
    return (A / sigma) * np.exp(-(x - mu)**2 / (2. * sigma**2)) + y0 + m*x

def lorentzian(x, *p):
    mu, A, gamma, y0 = p
    return (A * gamma**2) / ((x - mu)**2 + gamma**2) + y0

def background(x, *p):
    """
    Arbitrary model background data for absorbance FFT
    for the intention of replacing a peak in the FFT
    with the background
    """
    a, b = p
    return a * (1/x)**b

def gaussWithBackground(x, *p):
    pGauss = p[:4]
    a, b = p[4:]
    return gauss(x, *pGauss) + background(x, a, b)

####################
# Collection functions 
####################

def hsg_sum_spectra(object_list):
    """
    This function will add all the things that should be added.  Obvs.  It will
    also calculate the standard error of the mean for every NIR frequency.  The
    standard error is the sum of the dark noise and the "shot" noise.

    Also, remember, we're adding together and averaging the counts per pulse.  
    Hence why we use self.proc_data after it has been divided by the number of 
    fel pulses.  

    object_list: A list of spectrum objects
    """
    print "I'm trying!"
    
    good_list = []
    for index in xrange(len(object_list)):
        dark_var = 0
        num_images = 0
        try:
            temp = object_list.pop(0)
            stderr_holder = np.array(temp.proc_data[:, 1]).reshape((1600, 1))
            # print "Standard error holder shape 1:", stderr_holder.shape
        except:
            # print "God damn it, Leroy"
            break
        #print "temp has series: {}.\ttemp has cl: {}.\ttemp has series: {}".format(temp.parameters['series'], temp.parameters['center_lambda'], temp.parameters['series'])
        for spec in list(object_list):
            #print "\tspec has series: {}.\tspec has cl: {}.\tspec has fn: {}".format(spec.parameters['series'], spec.parameters['center_lambda'], spec.fname[-16:-13])
            #print "I am trying to add", temp.parameters['FELP'], spec.parameters['FELP']
            if temp.parameters['series'] == spec.parameters['series']:
                if temp.parameters['center_lambda'] == spec.parameters['center_lambda']:
                    temp += spec
                    num_images += 1
                    stderr_holder = np.hstack((stderr_holder, temp.proc_data[:, 1].reshape((1600,1))))
                    #print "Individual dark_stdev:", spec.dark_stdev
                    dark_var += (spec.dark_stdev)**2
                    #print "Standard error holder shape 2:", stderr_holder.shape
                    #print "\t\tadded"
                    #print "I ADDED", temp.parameters['FELP'], spec.parameters['FELP']
                    object_list.remove(spec)
        spec_number = stderr_holder.shape[1]
        std_error = np.sqrt(np.var(stderr_holder, axis=1, dtype=np.float64) + dark_var) / np.sqrt(spec_number) # Checking some sigma stuff from curve_fit
        # This standard error is for every point.  I think it actually overestimates
        # the error at places with no signal because we add the dark variance
        # effectively twice.
        #print "final dark_stdev:", np.sqrt(dark_var)
        temp.add_std_error(std_error)
        temp.image_normalize(num_images)
        good_list.append(temp)
    return good_list

def hsg_combine_spectra(spectra_list):
    """
    This function is all about smooshing different parts of the same hsg 
    spectrum together.  It takes a list of HighSidebandCCD spectra and turns the
    zeroth spec_step into a FullHighSideband object.  It then uses the function
    stitch_hsg_dicts over and over again for the smooshing.  

    Input:
    spectra_list = list of HighSidebandCCD objects that have sideband spectra
                   larger than the spectrometer can see.

    Returns:
    good_list = A list of FullHighSideband objects that have been combined as 
                much as can be.  
    """
    good_list = []
    spectra_list.sort(key=lambda x: x.parameters["spec_step"])
    for index in xrange(len(spectra_list)):
        try:
            temp = spectra_list.pop(0)
        except:
            break

        good_list.append(FullHighSideband(temp))

        counter = temp.parameters["spec_step"] + 1

        for piece in spectra_list:
            if temp.parameters["series"] == piece.parameters["series"]:
                if piece.parameters["spec_step"] == counter:
                    good_list[-1].add_CCD(piece)
                    spectra_list.remove(piece)
                    counter += 1
        good_list[-1].make_results_array()
    return good_list

def stitch_abs_results(main, new):
    raise NotImplementedError

####################
# Helper functions 
####################

def fft_filter(data, cutoffFrequency = 1520, inspectPlots = False, tryFitting = False, freqSigma = 50, ftol = 1e-4, isInteractive=False):
    """
    Performs an FFT, then fits a peak in frequency around the
    input with the input width.

    If only data is given, it will cut off all frequencies above the default value.

    inspectPlots = True will plot the FFT and the filtering at each step, as well as the results

    tryFitting = True will try to fit the peak in frequency space centered at the cutoffFrequency
    and with a width of freqSigma, using the background function above. Will replace
    the peak with the background function. Feature not very well tested

    isInteractive: Will pop up interactive windows to move the cutoff frequency and view the
    FFT in real time. Requires pyqtgraph and PyQt4 installed (pyqt4 is standard with
    anaconda/winpython, but pyqtgraph is not)
    """
    # Make a copy so we can return the same thing
    retData = np.array(data)
    x = np.array(data[:,0])
    y = np.array(data[:,-1])
    # Let's you place with zero padding.
    zeroPadding = len(x)
    N = len(x)

    if isInteractive:
        try:
            import pyqtgraph as pg
            from PyQt4 import QtCore, QtGui
        except:
            raise ImportError("Cannot do interactive plotting without pyqtgraph installed")
        # Need to make some basic classes fir signals and slots to make things simple
        class FFTWin(pg.PlotWindow):
            sigCutoffChanged = QtCore.pyqtSignal(object)
            sigClosed = QtCore.pyqtSignal()
            def __init__(self, x, y):
                super(FFTWin, self).__init__()
                # Plot the log of the data,
                # it breaks text boxes to do semilogy
                self.plotItem.plot(x, np.log10(y))
                # The line for picking the cutoff
                # Connect signals so the textbox updates and the
                # realspace window can recalcualte the FFT
                self.line = pg.InfiniteLine(cutoffFrequency, movable=True)
                self.line.sigPositionChanged.connect(lambda x:self.sigCutoffChanged.emit(x.value()))
                self.line.sigPositionChanged.connect(self.updateText)
                self.addItem(self.line)
                # Set up the textbox so user knows the frequency
                # If this ends up being useful, may need
                # a way to set the cutoff manually
                self.text = pg.TextItem("{:.4f}".format(cutoffFrequency))
                self.addItem(self.text)
                self.text.setPos(min(x), max(np.log10(y)))

                # Cheap magic to get the close event
                # of the main window. Need to keep a reference
                # to the old function so that we can call it
                # to properly clean up afterwards
                self.oldCloseEvent = self.win.closeEvent
                self.win.closeEvent = self.closeEvent
            def updateText(self, val):
                self.text.setText("{:.4f}".format(val.value()))

            def closeEvent(self, ev):
                # Just emit that we've been closed and
                # pass it along to the window closer
                self.sigClosed.emit()
                self.oldCloseEvent(ev)

        class RealWin(pg.PlotWindow):
            sigClosed = QtCore.pyqtSignal()
            def __init__(self, data, fftWin):
                super(RealWin, self).__init__()
                # To connect signals from it
                self.fftWin = fftWin
                self.data = data

                # Start off with the FFT given by the original
                # inputted cutoff
                self.updatePlot(cutoffFrequency)

                # See above comments
                self.oldClose = self.win.closeEvent
                self.win.closeEvent = self.closeEvent
                fftWin.sigCutoffChanged.connect(self.updatePlot)
                # Close self if other window is closed
                fftWin.sigClosed.connect(self.win.close)

            def updatePlot(self, val):
                self.plotItem.clear()
                self.plotItem.plot(*self.data.T, pen=pg.mkPen(width=3))
                # Recursion! Call this same function to do the FFT
                newData = fft_filter(self.data, cutoffFrequency=val)
                self.plotItem.plot(*newData.T, pen=pg.mkPen('r', width=3))

            def closeEvent(self, ev):
                self.sigClosed.emit()
                try:
                    self.fftWin.win.close()
                except:
                    pass
                self.oldClose(ev)


        k = fft.fftfreq(zeroPadding, x[1]-x[0])
        Y = fft.fft(y, n=zeroPadding)
        # Make the windows
        fftWin = FFTWin(k, np.abs(Y))
        realWin = RealWin(np.array(retData), fftWin)
        realWin.show()
        # Need to pause the program until the frequency is selected
        # Done with this qeventloop.
        loop = QtCore.QEventLoop()
        realWin.sigClosed.connect(loop.exit)
        loop.exec_()
        # Return with the desired output value
        return fft_filter(retData, fftWin.line.value())


    if inspectPlots:
        plt.figure("Real Space")
        plt.plot(x, y, label="Input Data")


    # Replicate origin directy
    # http://www.originlab.com/doc/Origin-Help/Smooth-Algorithm
    # "rotate" the data set so it ends at 0,
    # enforcing a periodicity in the data. Otherwise
    # oscillatory artifacts result at the ends
    onePerc = int(0.01*N)
    x1 = np.mean(x[:onePerc])
    x2 = np.mean(x[-onePerc:])
    y1 = np.mean(y[:onePerc])
    y2 = np.mean(y[-onePerc:])

    m = (y1-y2)/(x1-x2)
    b = y1 - m * x1

    flattenLine = m * x + b
    y -= flattenLine

    if inspectPlots:
        plt.plot(x, y, label="Rotated Data")

    # Perform the FFT and find the appropriate frequency spacing
    k = fft.fftfreq(zeroPadding, x[1]-x[0])
    Y = fft.fft(y, n=zeroPadding)
    if inspectPlots:
        plt.figure("Frequency Space")
        plt.semilogy(k, np.abs(Y), label="Raw FFT")

    if tryFitting:
        try:
            # take +/- 4 sigma points around peak to fit to
            sl = np.abs(k-cutoffFrequency).argmin() + np.array([-1, 1]) * 10 * freqSigma / np.abs(k[0]-k[1])
            sl = slice(*[int(j) for j in sl])
            p0 = [cutoffFrequency,
                  np.abs(Y)[sl].max() * freqSigma, # estimate the height baased on the max in the set
                  freqSigma,
                  0.14, 2e3, 1.1] # magic test numbers, they fit the background well


            if inspectPlots:
                plt.semilogy(k[sl], gaussWithBackground(k[sl], *p0), label="Peak with initial values")
            p, _ = curve_fit(gaussWithBackground, k[sl], np.abs(Y)[sl], p0=p0, ftol=ftol)
            if inspectPlots:
                plt.semilogy(k[sl], gaussWithBackground(k[sl], *p), label="Fitted Peak")


            # Want to remove data within 5 sigma ( arb value... )
            st = int(p[0] - 5*p[2])
            en = int(p[0] + 5*p[2])

            # Find get the indices to remove.
            refitRangeIdx = np.argwhere((k>st) & (k<en))
            refitRangeIdxNeg = np.argwhere((k<-st) & (k>-en))

            # Replace the data with the backgroudn
            # Note: abuses the symmetry of the FFT of a real function
            # to get the negative side of the data
            Y[refitRangeIdx] = background(k[refitRangeIdx], *p[-2:])
            Y[refitRangeIdxNeg] = background(k[refitRangeIdx], *p[-2:])[::-1]
        except:
            print "ERROR: Trouble fitting the peak in frequency space.\n\t Defaulting to cutting off"

            # Assume cutoffFrequency was the peak, not the actual cutoff
            # Leaving it alone means half the peak would remain and the data
            # wouldn't really be smoothed
            cutoffFrequency -= 5 * freqSigma

            # Reset this so the next part gets called
            tryFitting = False

    # "if not" instead of "else" because if the above
    # fitting fails, we can default to the sharp cutoff
    if not tryFitting:
        # Define where to remove the data
        st = cutoffFrequency
        en = int(max(k))+1

        # Find the indices to remove the data
        refitRangeIdx = np.argwhere((k>st) & (k<en))
        refitRangeIdxNeg = np.argwhere((k<-st) & (k>-en))

        # Kill it all after the cutoff
        Y[refitRangeIdx] = 0
        Y[refitRangeIdxNeg] = 0

        smoothIdx = np.argwhere((-st<k) & (k< st))
        smoothr = -1./cutoffFrequency**2 * k[smoothIdx]**2 + 1

        Y[smoothIdx] *= smoothr

    if inspectPlots:
        plt.plot(k, np.abs(Y), label="FFT with removed parts")
        a = plt.legend()
        a.draggable(True)

    # invert the FFT
    y = fft.ifft(Y, n=zeroPadding)

    # unshift the data
    y += flattenLine

    # using fft, not rfft, so data may have some
    # complex parts. But we can assume they'll be negligible and
    # remove them
    # ( Safer to use np.real, not np.abs? )
    # Need the [:len] to remove zero-padded stuff
    y = np.abs(y)[:len(x)]

    if inspectPlots:
        plt.figure("Real Space")
        print x.size, y.size
        plt.plot(x, y, label="Smoothed Data")
        a = plt.legend()
        a.draggable(True)

    retData[:,0] = x
    retData[:,-1] = y
    return retData

def stitchData(dataList, plot=False):
    """
    Attempt to stitch together absorbance data. Will translate the second data set
    to minimize leastsq between the two data sets.
    :param dataList: Iterable of the data sets to be fit. Currently
            it only takes the first two elements of the list, but should be fairly
            straightforward to recursivly handle a list>2. Shifts the second
            data set to overlap the first

             elements of dataList can be either np.arrays or Absorbance class,
              where it will take the proc_data itself
    :param plot: bool whether or not you want the fit iterations to be plotted
            (for debugging)
    :return: a, a (2,) np.array of the shift
    """

    # Data coercsion, make sure we know what we're working wtih
    first = dataList[0]
    if isinstance(first, Absorbance):
        first = first.proc_data
        second = dataList[1]
    if isinstance(second, Absorbance):
        second = second.proc_data
    if plot:
        # Keep a reference to whatever plot is open at call-time
        # Useful if the calling script has plots before and after, as
        # omitting this will cause future plots to be added to figures here
        firstFig = plt.gcf()
        plt.figure("Stitcher")
        # Plot the raw input data
        plt.plot(*first.T)
        plt.plot(*second.T)

    # Algorithm is set up such that the "second" data set spans the
    # higher domain than first. Need to enforce this, and remember it
    # so the correct shift is applied
    flipped = False
    if max(first[:,0])>max(second[:,0]):
        flipped = True
        first, second = second, first

def fitter(p, shiftable, immutable):
    # Function for leastsq to minimize

    # Get the shifts
    dx = p[0]
    dy = p[1]

    # Don't want pass-by-reference nonsense, recast our own refs
    shiftable = np.array(shiftable)
    immutable = np.array(immutable)

    # Shift the data set
    shiftable[:,1]+=dy
    shiftable[:,0]+=dx

    # Create an interpolator. We want a
    # direct comparision for subtracting the two functions
    # Different spec grating positions have different wavelengths
    # so they're not directly comparable.
    shiftF = spi.interp1d(*shiftable.T)

    # Find the bounds of where the two data sets overlap
    overlap = (min(shiftable[:,0]), max(immutable[:,0]))
        #print "overlap", overlap

    # Determine the indices of the immutable function
    # where it overlaps. argwhere returns 2-d thing,
    # requiring the [0] at the end of each call
    fOlIdx = (min(np.argwhere(immutable[:,0]>=overlap[0]))[0],
              max(np.argwhere(immutable[:,0]<=overlap[1]))[0])

    # Get the interpolated values of the shiftable function at the same
    # x-coordinates as the immutable case
    newShift = shiftF(immutable[fOlIdx[0]:fOlIdx[1],0])

    if plot:
        plt.plot(*immutable[fOlIdx[0]:fOlIdx[1],:].T, marker='o', label="imm", markersize=10)
        plt.plot(immutable[fOlIdx[0]:fOlIdx[1], 0], newShift, marker='o', label="shift")
    imm = immutable[fOlIdx[0]:fOlIdx[1],1]
    shift = newShift
    return imm-shift

    a, _, _, msg, err= spo.leastsq(fitter, [0.0001, 0.01*max(first[:,1])], args=(second, first), full_output = 1)
    # print "a", a
    if plot:
        # Revert back to the original figure, as per top comments
        plt.figure(firstFig.number)

    # Need to invert the shift if we flipped which
    # model we're supposed to move
    if flipped: a*=-1

    return a

def stitch_hsg_dicts(full, new_dict, need_ratio=False, verbose=False):
    """
    This helper function takes a FullHighSideband.full_dict attribute and a sideband 
    object, either CCD or PMT and smushes the new sb_results into the full_dict.

    The first input doesn't change, so f there's a PMT set of data involved, it 
    should be in the full variable to keep the laser normalization intact.

    Inputs:
    full = full_dict from FullHighSideband, or HighSidebandPMT.  It's important 
           that it contains lower orders than the new_dict.
    new_dict = another full_dict.
    need_ratio = If gain or other parameters aren't equal and must resort to 
                 calculating the ratio instead of the measurements being equivalent.
                 Changing integration time still means N photons made M counts, 
                 but changing gain or using PMT or whatever does affect things.

    Returns:
    full = extended version of the input full.  Overlapping sidebands are 
           averaged because that makes sense?
    """
    #print "I'm adding these sidebands", sorted(new_dict.keys())
    overlap = []
    for new_sb in sorted(new_dict.keys()):
        if new_sb in full.keys():
            overlap.append(new_sb)
    if verbose:
        print "overlap:", overlap

    new_starter = overlap[-1]
    overlap = [x for x in overlap if (x%2 == 0) and (x != min(overlap) and (x != max(overlap)))]

    if need_ratio:
    # Calculate the appropriate ratio to multiply the new sidebands by.
    # I'm not entirely sure what to do with the error of this guy.
        ratio_list = []
        for sb in overlap:
            ratio_list.append(full[sb][2] / new_dict[sb][2])
        ratio = np.mean(ratio_list)
        error = np.std(ratio_list) / np.sqrt(len(ratio_list))
        print "Ratio list", ratio_list
        print "Ratio", ratio
        print "Error", error
    # Adding the new sidebands to the full set and moving errors around.
    # I don't know exactly what to do about the other aspecs of the sidebands
    # besides the strength and its error.
        for sb in overlap:
            full[sb][2] = ratio * new_dict[sb][2]
            full[sb][3] = full[sb][2] * np.sqrt((error / ratio)**2 + (new_dict[3] / new_dict[2])**2)
            
            # Now for linewidths
            lw_error = np.sqrt(full[sb][5]**(-2) + new_dict[sb][5]**(-2))**(-1)
            lw_avg = (full[sb][4] / (full[sb][5]**2) + new_dict[sb][4] / (new_dict[sb][5]**2)) / (full[sb][5]**(-2) + new_dict[sb][3]**(-2))
            full[sb][4] = lw_avg
            full[sb][5] = lw_error

    else:
        for sb in overlap:
            error = np.sqrt(full[sb][3]**(-2) + new_dict[sb][3]**(-2))**(-1)
            avg = (full[sb][2] / (full[sb][3]**2) + new_dict[sb][2] / (new_dict[sb][3]**2)) / (full[sb][3]**(-2) + new_dict[sb][3]**(-2))
            full[sb][2] = avg
            full[sb][3] = error

            lw_error = np.sqrt(full[sb][5]**(-2) + new_dict[sb][5]**(-2))**(-1)
            lw_avg = (full[sb][4] / (full[sb][5]**2) + new_dict[sb][4] / (new_dict[sb][5]**2)) / (full[sb][5]**(-2) + new_dict[sb][3]**(-2))
            full[sb][4] = lw_avg
            full[sb][5] = lw_error

    for sb in [x for x in new_dict.keys() if (x >= new_starter)]:
        full[sb] = new_dict[sb]
    #print "I made this dictionary", sorted(full.keys())
    return full

def save_parameter_sweep(spectrum_list, file_name, folder_str, param_name, unit, verbose=False):
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
    spectrum_list.sort(key=lambda x: x.parameters[param_name])
    included_spectra = dict()
    param_array = None
    sb_included = []
    
    for spec in spectrum_list:
        sb_included = sorted(list(set(sb_included + spec.full_dict.keys())))
        included_spectra[spec.fname.split('/')[-1]] = spec.parameters[param_name]
        # If these are from summed spectra, then only the the first file name
        # from that sum will show up here, which should be fine?
    if verbose:
        #print "full name:", spectrum_list[0].fname
        print "included names:", included_spectra
        print "sb_included:", sb_included
    
    
    for spec in spectrum_list:
        temp_dict = {} # This is different from full_dict in that the list has the
                       # sideband order as the zeroth element.
        if verbose:
            print "the sb_results:", spec.sb_results
        for index in xrange(len(spec.sb_results[:, 0])):
            if verbose:
                print "my array slice:", spec.sb_results[index, :]
            temp_dict[int(round(spec.sb_results[index, 0]))] = np.array(spec.sb_results[index, :])
        
        if verbose:
            print temp_dict
        
        for sb in sb_included:
            blank = np.zeros(7)
            blank[0] = float(sb)
            #print "checking sideband order:", sb
            #print "blank", blank
            if not temp_dict.has_key(sb):
                #print "\nNeed to add sideband order:", sb
                temp_dict[sb] = blank
        try: # Why is this try-except here?
            spec_data = np.array([float(spec.parameters[param_name])])
        except:
            spec_data = np.array([float(spec.parameters[param_name][:2])])
        for key in sorted(temp_dict.keys()):
            #print "I am going to hstack this:", temp_dict[key]
            spec_data = np.hstack((spec_data, temp_dict[key]))
            
        try:
            param_array = np.vstack((param_array, spec_data))
        except:
            param_array = np.array(spec_data)
        if verbose:
            print "The shape of the param_array is:", param_array.shape
        #print "The param_array itself is:", param_array
    
    param_array_norm = np.array(param_array).T # python iterates over rows
    for elem in [x for x in xrange(len(param_array_norm)) if (x-1)%7 == 3]:
        temp_max = np.max(param_array_norm[elem])
        param_array_norm[elem] = param_array_norm[elem] / temp_max
        param_array_norm[elem + 1] = param_array_norm[elem + 1] / temp_max


    param_array_norm = np.array(param_array_norm.T)


    try:
        os.mkdir(folder_str)
    except OSError, e:
        if e.errno == errno.EEXIST:
            pass
        else:
            raise
    norm_name = file_name + '_norm.txt'
    file_name = file_name + '.txt'
    
    try:
        included_spectra_str = json.dumps(included_spectra, sort_keys=True)
    except:
        print "Source: save_parameter_sweep\nJSON FAILED"
        return
    origin_import1 = param_name
    origin_import2 = unit
    origin_import3 = ""
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
    np.savetxt(os.path.join(folder_str, norm_name), param_array_norm, delimiter=',', 
               header=header, comments='', fmt='%f')
    print "Saved the file.\nDirectory: {}".format(os.path.join(folder_str, file_name))




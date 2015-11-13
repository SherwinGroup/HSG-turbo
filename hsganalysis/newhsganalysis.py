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
import scipy.ndimage as ndimage

####################
# Objects 
####################

class CCD(object):
    def __init__(self, fname, spectrometer_offset=None):
        """
        This will read the appropriate file and make a basic CCD object.  Fancier
        things will be handled with the sub classes.

        input:
        fname = file name where the data is saved
        spectrometer_offset = number of nanometers that the spectrometer is 
                              shifted by

        creates:
        self.parameters = JSON dictionary holding all of the information from the
                          data file.
        self.description = string that is the text box from data taking GUI
        self.raw_data = raw data output by measurement software, wavelength vs.
                        data.  There may be text for some of the entries
        self.ccd_data = semi-processed 1600 x 3 array of photon energy vs. data with error
        """
        self.fname = fname
        
        print "I'm going to open", fname
        #f = open(fname,'rU')
        with open(fname,'rU') as f:
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
        try:
            self.parameters["spec_step"] = int(self.parameters["spec_step"])
        except ValueError:
            self.parameters["spec_step"] = 0
        except KeyError:
            pass
        self.raw_data = np.flipud(np.genfromtxt(fname, comments='#', delimiter=','))
        # I used flipup so that the x-axis is an increasing function of frequency
        
        self.ccd_data = np.array(self.raw_data[:1600,:]) # By slicing up to 1600,
                                                         # we cut out the text 
                                                         # header
        if spectrometer_offset is not None:
            try:
                # print "Doing offset", self.parameters["offset"]
                self.ccd_data[:, 0] += float(self.parameters["offset"])
                # print "it worked"
            except:
                self.ccd_data[:, 0] += spectrometer_offset
                # print "it didn't work"
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
                   header=spec_header, comments='', fmt='%0.6e')

        print "Save image.\nDirectory: {}".format(os.path.join(folder_str, spectra_fname))

class NeonNoiseAnalysis(CCD):
    def __init__(self, fname, spectrometer_offset=None):
        super(HighSidebandCCD, self).__init__(fname, spectrometer_offset=spectrometer_offset)

        self.addenda = self.parameters['addenda']
        self.subtrahenda = self.parameters['subtrahenda']

        self.noise_and_signal()
        self.process_stuff()

    def noise_and_signal(self):
        """
        This bad boy calculates the standard deviation of the space between the 
        neon lines.

        The noise regions are, in nm:
        high: 784-792
        low1: 795-806
        low2: 815-823
        low3: 831-834

        the peaks are located at, in nm:
        #1, weak: 793.6
        #2, medium: 794.3 
        #3, medium: 808.2
        #4, weak: 825.9
        #5, strong: 830.0
        """
        self.high_noise_region = np.array(self.ccd_data[30:230, :])
        self.low_noise_region1 = np.array(self.ccd_data[380:700, :])
        self.low_noise_region2 = np.array(self.ccd_data[950:1200, :])
        self.low_noise_region3 = np.array(self.ccd_data[1446:1546, :])

        self.high_noise = np.std(self.high_noise_region[:, 1])
        self.low_noise1 = np.std(self.low_noise_region1[:, 1])
        self.low_noise2 = np.std(self.low_noise_region2[:, 1])
        self.low_noise3 = np.std(self.low_noise_region3[:, 1])

        self.noise_list = [self.high_noise, self.low_noise1, self.low_noise2, self.low_noise3]
        print "Noise list:", self.noise_list

        self.peak1 = np.array(self.ccd_data[303:323, :])
        self.peak2 = np.array(self.ccd_data[319:339, :])
        self.peak3 = np.array(self.ccd_data[736:746, :])
        self.peak4 = np.array(self.ccd_data[1268:1288, :])
        self.peak5 = np.array(self.ccd_data[1381:1421, :])

        self.signal1 = np.sum(self.peak1[:, 1])
        self.error1 = np.sqrt(np.sum(self.peak1[:, 2]**2))
        self.signal2 = np.sum(self.peak2[:, 1])
        self.error2 = np.sqrt(np.sum(self.peak2[:, 2]**2))
        self.signal3 = np.sum(self.peak3[:, 1])
        self.error3 = np.sqrt(np.sum(self.peak3[:, 2]**2))
        self.signal4 = np.sum(self.peak4[:, 1])
        self.error4 = np.sqrt(np.sum(self.peak4[:, 2]**2))
        self.signal5 = np.sum(self.peak5[:, 1])
        self.error5 = np.sqrt(np.sum(self.peak5[:, 2]**2))

        self.signal_list = [self.signal1, self.signal2, self.signal3, self.signal4, self.signal5]
        self.error_list = [self.error1, self.error2, self.error3, self.error4, self.error5]
        print "Signal list:", self.signal_list

    def process_stuff(self):
        """
        This one puts high_noise, low_noise1, signal2, and error2 in a nice horizontal array
        """
        self.results = np.array([self.high_noise, self.low_noise1, self.signal2, self.error2])

def collect_noise(neon_list, param_name, folder_name, file_name):
    """
    This function acts like save parameter sweep.

    param_name = string that we're gonna save!
    """
    
    for elem in neon_list:
        temp = np.column_stack((elem.parameters[param_name], elem.results))
        try:
            param_array = np.row_stack(param_array, elem.results)
        except:
            param_array = np.array(elem.results)
    
    try:
        os.mkdir(folder_name)
    except OSError, e:
        if e.errno == errno.EEXIST:
            pass
        else:
            raise

    file_name = file_name + '.txt'

    origin_import1 = param_name + ",Noise,Noise,Signal,error"
    origin_import2 = ",counts,counts,counts,counts"
    origin_import3 = ",High noise region,Low noise region,808nm peak,808nm peak error"

    header_total = origin_import1 + "\n" + origin_import2 + "\n" + origin_import3

    #print "Spec header: ", spec_header
    print "the param_array is:", param_array
    np.savetxt(os.path.join(folder_name, file_name), param_array, delimiter=',', 
               header=header_total, comments='', fmt='%0.6e')
    print "Saved the file.\nDirectory: {}".format(os.path.join(folder_name, file_name))

class HighSidebandCCD(CCD):
    def __init__(self, fname, spectrometer_offset=None):
        """
        This will read the appropriate file.  The header needs to be fixed to
        reflect the changes to the output header from the Andor file.  Because
        another helper file will do the cleaning and background subtraction,
        those are no longer part of this init.  This also turns all wavelengths
        from nm (NIR ones) or cm-1 (THz ones) into eV.
        
        Input:
        fname = file name of the hsg spectrum from CCD superclass
        spectrometer_offset = number of nanometers the spectrometer is off by, 
                              should be 0.0...but can be 0.2 or 1.0
        
        Internal:
        self.fname = the filename
        self.parameters = string with all the relevant experimental perameters
        self.description = the description we added to the file as the data
                           was being taken
        self.proc_data = processed data that has gone is frequency vs counts/pulse
        self.dark_stdev = this is not currently handled appropriately
        self.addenda = the list of things that have been added to the file, in
                       form of [constant, *spectra_added]
        self.subtrahenda = the list of spectra that have been subtracted from
                           the file.  Constant subtraction is dealt with with
                           self.addenda
        """
        super(HighSidebandCCD, self).__init__(fname, spectrometer_offset=spectrometer_offset)

        self.proc_data = np.array(self.ccd_data) # Does this work the way I want it to?

        self.parameters["nir_freq"] = 1239.84 / float(self.parameters["nir_lambda"])
        self.parameters["thz_freq"] = 0.000123984 * float(self.parameters["fel_lambda"])
        self.parameters["nir_power"] = float(self.parameters["nir_power"])
        self.parameters["thz_power"] = float(self.parameters["fel_power"])

        self.addenda = self.parameters['addenda']
        self.subtrahenda = self.parameters['subtrahenda']

    def __add__(self, other):
        """
        Add together the image data from self.proc_data, or add a constant to 
        that np.array.  It will then combine the addenda and subtrahenda lists,
        as well as add the fel_pulses together.  If type(other) is a CCD object,
        then it will add the errors as well.

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
                ret.proc_data[:, 2] = np.sqrt(self.proc_data[:, 1]**2 + other.proc_data[:, 1]**2)
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
            ret.proc_data[:, 1] = self.proc_data[:, 1] - other # Need to choose a name
            ret.addenda[0] = ret.addenda[0] - other
            
        # Subtract the data of two hsg_spectra from each other
        else:
            if np.isclose(ret.proc_data[0, 0], other.proc_data[0, 0]):
                ret.proc_data[:, 1] = self.proc_data[:, 1] - other.proc_data[:, 1]
                ret.proc_data[:, 2] = np.sqrt(self.proc_data[:, 1]**2 + other.proc_data[:, 1]**2)
                ret.subtrahenda.extend(other.addenda[1:])
                ret.addenda.extend(other.subtrahenda)
            else:
                raise Exception('Source: Spectrum.__sub__:\nThese are not from the same grating settings')
        return ret
    
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
    
    def guess_sidebands(self, cutoff=5, verbose=False):
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

        if global_max < 15:
            check_y = y_axis[:global_max + 15]
        elif global_max > 1585:
            check_y = y_axis[global_max - 15:]
        else:
            check_y = y_axis[global_max - 15:global_max + 15]

        check_max_area = np.sum(y_axis[global_max - 1:global_max + 2])
        check_ave = np.mean(check_y)
        check_stdev = np.std(check_y)

        check_ratio = (check_max_area - 3 * check_ave) / check_stdev
        
        
        if verbose:
            print "\nI'm checking", self.fname, "\n"
            print "Global max checking:", global_max, check_y
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
            window_size = 0.28 + 0.0004 * order # used to be last_sb?
            lo_freq_bound = last_sb - thz_freq * (1 + window_size) # Not sure what to do about these
            hi_freq_bound = last_sb - thz_freq * (1 - window_size)

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

            check_max_index = np.argmax(check_y) # This assumes that two floats won't be identical
            check_max_area = np.sum(check_y[check_max_index - 1:check_max_index + 2])
            check_ave = np.mean(check_y)
            check_stdev = np.std(check_y[[0,1,2,3,-1,-2,-3,-4]]) 
                # So the sideband isn't included in the noise calculation
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
            window_size = 0.28 + 0.0004 * order # used to be last_sb?
            lo_freq_bound = last_sb + thz_freq * (1 - window_size) # Not sure what to do about these
            hi_freq_bound = last_sb + thz_freq * (1 + window_size)
            start_index = False
            end_index = False

            if verbose:
                print "\nSideband", order, "\n"
                print "lower bound", lo_freq_bound
                print "upper bound", hi_freq_bound ,"\n"        
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

            check_max_index = np.argmax(check_y) # This assumes that two floats won't be identical
            check_max_area = np.sum(check_y[check_max_index - 1:check_max_index + 2])
            check_ave = np.mean(check_y)
            check_stdev = np.std(check_y[[0,1,2,3,-1,-2,-3,-4]])
            check_ratio = (check_max_area - 3 * check_ave) / check_stdev

            if verbose:
                print "Start and end index:", start_index, end_index
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
            if self.sb_index[elem] < 30:
                data_temp = self.proc_data[:self.sb_index[elem] + 30, :]
            elif (1600 - self.sb_index[elem]) < 30:
                data_temp = self.proc_data[self.sb_index[elem] - 30:, :]
            else:
                data_temp = self.proc_data[self.sb_index[elem] - 30:self.sb_index[elem] + 30, :]
            width_guess = 0.0001 + 0.000001*self.sb_list[elem] # so the width guess gets wider as order goes up
            p0 = [self.sb_guess[elem, 0], self.sb_guess[elem, 1] * width_guess, width_guess, 0.1]
            #print "Let's fit this shit!"
            if verbose:
                print "number:", elem, num
                #print "data_temp:", data_temp
                print "p0:", p0
            if plot:
                plt.figure('CCD data')
                linewidth = 3
                x_vals = np.linspace(data_temp[0, 0], data_temp[-1, 0], num=500)
                if elem!=0:
                    try:
                        plt.plot(x_vals, gauss(x_vals, *p0), 
                                 plt.gca().get_lines()[-1].get_color()+'--' # I don't really know. Mostly
                                                             # just looked around at what functions
                                                             # matplotlib has...
                                 , linewidth=linewidth)
                    except: # to prevent weird mac issues with the matplotlib things?
                        plt.plot(x_vals, gauss(x_vals, *p0), '--', linewidth=linewidth)
                                     
                else:
                    plt.plot(x_vals, gauss(x_vals, *p0), '--', linewidth=linewidth)

            try:
                if verbose:
                    print "I'm going to try to fit", self.sb_list[elem]
                coeff, var_list = curve_fit(gauss, data_temp[:, 0], data_temp[:, 1], p0=p0)
                if verbose:
                    print "the fit worked"
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
                print "It's sideband", num
                print "In file", self.fname
                self.sb_list[elem] = None
        sb_fits_temp = np.asarray(sb_fits)
        reorder = [0, 1, 5, 2, 6, 3, 7, 4, 8]
        try:
            sb_fits = sb_fits_temp[:, reorder]
        except:
            print "The file is:", self.fname
            print "\n!!!!!\nSHIT WENT WRONG\n!!!!!\n"
                
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
        
        Format:
        spectra_fname = file_name + '_' + marker + '_' + str(index) + '.txt'
        fit_fname = file_name + '_' + marker + '_' + str(index) + '_fits.txt'

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
        temp = np.array(self.sb_results)
        
        ampli = np.array([temp[:, 3] / temp[:, 5]]) # But [:, 3] is already area?
                                                    # (The old name was area)
                                                    # I think it must be amplitude
        temp[:, 5:7] = temp[:, 5:7] * 1000 # For meV linewidths
        print "sb_results", self.sb_results.shape
        print "ampli", ampli.shape
        save_results = np.hstack((temp, ampli.T))
        
        
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
        
        origin_import_fits = '\nSideband,Center energy,error,Sideband strength,error,Linewidth,error,Amplitude\norder,eV,,arb. u.,,meV,,arb. u.\n' + marker
        fits_header = '#' + parameter_str + '\n#' + self.description[:-2] + origin_import_fits
        
        np.savetxt(os.path.join(folder_str, spectra_fname), self.proc_data, delimiter=',',
                   header=spec_header, comments='', fmt='%0.6e')
        np.savetxt(os.path.join(folder_str, fit_fname), save_results, delimiter=',',
                   header=fits_header, comments='', fmt='%0.6e')

        print "Save image.\nDirectory: {}".format(os.path.join(folder_str, spectra_fname))



class PMT(object):
    def __init__(self, file_path):
        """
        Initializes a SPEX spectrum.  It'll open a file, and bring in the details
        of a sideband spectrum into the object.  There isn't currently any reason
        to use inheritance here, but it could be extended later to include PLE or
        something of the sort.
        
        attributes:
            self.parameters - dictionary of important experimental parameters
                              this will not necessarily be the same for each
                              file in the object
            self.description - string of the description of the file(s)
            self.sb_dict - keys are sideband order, values are PMT data arrays
            self.sb_list - sorted list of included sidebands
            self.files - included file paths
        """
        print "This started"

        self.files = [file_path]

        with open(file_path, 'rU') as f:
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

class HighSidebandPMT(PMT):
    def __init__(self, file_path, verbose=False):
        """
        Initializes a SPEX spectrum.  It'll open a single file, then read
        the data from that file using .add_sideband().  The super's init will handle the parameters
        and the description.  
        
        attributes:
            self.parameters - dictionary of important experimental parameters
            self.description - string of the description of the file(s)
            self.sb_dict - keys are sideband order, values are PMT data arrays
            self.sb_list - sorted list of included sidebands
        """
        super(HighSidebandPMT, self).__init__(file_path)
        self.sb_dict = {}
        self.sb_list = []
        self.add_sideband(file_path)

    def add_sideband(self, file_name):
        """
        This bad boy will add a sideband to the sideband spectrum of this object
        """
        self.files.append(file_name)
        with open(file_name, 'rU') as f:
            sb_num = int(f.readline().split(' ')[-1])
        raw_temp = np.genfromtxt(file_name, comments='#')#, delimiter=',')
        # Make things comma delimited?
        try:
            self.sb_dict[sb_num].vstack((raw_temp))
        except:
            self.sb_dict[sb_num] = np.array(raw_temp)

    def process_sidebands(self, verbose=False):
        """
        This bad boy will clean up the garbled mess that is the object before hand, 
        including clearing out misfired shots and doing the averaging.

        """
        for sb_deets in list(self.sb_dict.items()):
            sb_num, sb = sb_deets[0], sb_deets[1]
            frequencies = sorted(list(set(sb[:, 0])))
            fire_condition = np.mean(sb[:, 2]) / 2 # Say FEL fired if the 
                                                   # cavity dump signal is
                                                   # more than half the mean 
                                                   # of the cavity dump signal
            temp = None
            for freq in frequencies:
                data_temp = np.array([])
                for point in sb:
                    if point[0] == freq and point[2] > fire_condition:
                        data_temp = np.hstack((data_temp, point[3]))
                try:
                    temp = np.vstack((temp, np.array([freq, np.mean(data_temp), np.std(data_temp) / np.sqrt(len(data_temp))])))
                except:
                    temp = np.array([freq, np.mean(data_temp), np.std(data_temp) / np.sqrt(len(data_temp))])
            temp[:, 0] = temp[:, 0] / 8065.6 # turn NIR freq into eV
            temp = temp[temp[:, 0].argsort()]
            self.sb_dict[sb_num] = np.array(temp)
        self.sb_list = sorted(self.sb_dict.keys())
        if verbose:
            print "Sidebands included", self.sb_list
    
    def laser_line(self):
        """
        This method is designed to scale everything in the PMT to the conversion
        efficiency based on our measurement of the laser line with a fixed 
        attenuation.
        """

        if 0 not in self.sb_list:
            self.parameters['normalized?'] = False
            return
        else:
            laser_index = np.where(self.sb_results[:, 0] == 0)[0][0]
            print "sb_results", self.sb_results[:, 0]
            print "laser_index", laser_index

            laser_strength = np.array(self.sb_results[laser_index, 3:5])
            
            for sb in self.sb_results:
                print "Laser_strength", laser_strength
                sb[4] = (sb[3] / laser_strength[0]) * np.sqrt((sb[4] / sb[3])**2 + (laser_strength[1] / laser_strength[0])**2)
                sb[3] = sb[3] / laser_strength[0]
            for sb in self.full_dict.values():
                sb[3] = (sb[2] / laser_strength[0]) * np.sqrt((sb[3] / sb[2])**2 + (laser_strength[1] / laser_strength[0])**2)
                sb[2] = sb[2] / laser_strength[0]
            self.parameters['normalized?'] = True

    def integrate_sidebands(self, verbose=False):
        """
        This method will integrate the sidebands to find their strengths, and then
        use a magic number to get the width, since they are currently so utterly
        undersampled for fitting.
        """
        self.full_dict = {}
        for sideband in self.sb_dict.items():
            index = np.argmax(sideband[1][:, 1])
            nir_frequency = sideband[1][index, 0]
            area = np.trapz(np.nan_to_num(sideband[1][:, 1]), sideband[1][:, 0])
            error = np.sqrt(np.sum(np.nan_to_num(sideband[1][:, 2])**2)) / 8065.6 # Divide by the step size?
            if verbose:
                print "order", sideband[0]
                print "area", area
                print "error", error
            details = np.array([sideband[0], nir_frequency, 1/8065.6, area, error, 2/8065.6, 1/8065.6])
            if area < 0:
                continue
            try:
                self.sb_results = np.vstack((self.sb_results, details))
            except:
                self.sb_results = np.array(details)
            self.full_dict[sideband[0]] = details[1:]
        self.sb_results = self.sb_results[self.sb_results[:, 0].argsort()]
        
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
                print "Sideband data:\n", sideband[1]
            index = np.argmax(sideband[1][:, 1])
            nir_frequency = sideband[1][index, 0]
            peak = sideband[1][index, 1]
            width_guess = 0.0001 # Yep, another magic number
            p0 = [nir_frequency, peak * width_guess, width_guess, 0.00001]
            
            if verbose:
                x_vals = np.linspace(np.amin(sideband[1][:, 0]), np.amax(sideband[1][:, 0]), num=50)
                plt.plot(x_vals, gauss(x_vals, *p0))
                print "p0:", p0
            try: 
                coeff, var_list = curve_fit(gauss, sideband[1][:, 0], sideband[1][:, 1], sigma=sideband[1][:, 2], p0=p0)
                coeff[1] = abs(coeff[1])
                coeff[2] = abs(coeff[2])
                if verbose:
                    print "coeffs:", coeff
                    print "stdevs:", np.sqrt(np.diag(var_list))
                    print "integral", np.trapz(sideband[1][:, 1], sideband[1][:, 0])
                if np.sqrt(np.diag(var_list))[0] / coeff[0] < 0.5: # The error on where the sideband is should be small
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
                        #plt.plot(x_vals, gauss(x_vals, *p0))
                else:
                    print "what happened?"
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
            print "And the results, please:\n", self.sb_results

        self.full_dict = {}
        for sb in self.sb_results:
            self.full_dict[sb[0]] = np.asarray(sb[1:])
    
    def save_processing(self, file_name, folder_str, marker='', index=''):
        """
        This will save all of the self.proc_data and the results from the 
        fitting of this individual file.
        
        Format:
        spectra_fname = file_name + '_' + marker + '_' + str(index) + '.txt'
        fit_fname = file_name + '_' + marker + '_' + str(index) + '_fits.txt'

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
        self.parameters['included_files'] = list(self.files)
        try:
            parameter_str = json.dumps(self.parameters, sort_keys=True)
        except:
            print "Source: PMT.save_images\nJSON FAILED"
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
                complete = np.array(self.sb_dict[sideband])
        
        np.savetxt(os.path.join(folder_str, spectra_fname), complete, delimiter=',',
                   header=spec_header, comments='', fmt='%0.6e')
        
        np.savetxt(os.path.join(folder_str, fit_fname), self.sb_results, delimiter=',',
                   header=fits_header, comments='', fmt='%0.6e')
        
        print "Saved PMT spectrum.\nDirectory: {}".format(os.path.join(folder_str, spectra_fname))

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
        self.full_dict = stitch_hsg_dicts(pmt_object.full_dict, self.full_dict, need_ratio=True, verbose=True)
        self.parameters['files_here'].append(pmt_object.files)
        self.make_results_array()

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
        
        Format:
        fit_fname = file_name + '_' + marker + '_' + str(index) + '_full.txt'

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

        temp = np.array(self.sb_results)
        
        ampli = np.array([temp[:, 3] / temp[:, 5]]) # I'm pretty sure this is
                                                    # amplitude, not area
        temp[:, 5:7] = temp[:, 5:7] * 1000 # For meV linewidths
        print "sb_results", self.sb_results.shape
        print "ampli", ampli.shape
        save_results = np.hstack((temp, ampli.T))
        
        #spectra_fname = file_name + '_' + marker + '_' + str(index) + '.txt'
        fit_fname = file_name + '_' + marker + '_' + str(index) + '_full.txt'
        #self.save_name = spectra_fname
        
        #self.parameters['addenda'] = self.addenda
        #self.parameters['subtrahenda'] = self.subtrahenda
        try:
            parameter_str = json.dumps(self.parameters, sort_keys=True)
        except Exception as e:
            print e
            print "Source: EMCCD_image.save_images\nJSON FAILED"
            print "Here is the dictionary that broke JSON:\n", self.parameters
            return

        #origin_import_spec = '\nNIR frequency,Signal,Standard error\neV,arb. u.,arb. u.'
        #spec_header = '#' + parameter_str + '\n#' + self.description[:-2] + origin_import_spec
        
        origin_import_fits = '\nSideband,Center energy,error,Sideband strength,error,Linewidth,error,Amplitude\norder,eV,,arb. u.,,meV,,arb. u.\n' + marker
        fits_header = '#' + parameter_str + '\n#' + self.description[:-2] + origin_import_fits
        
        #np.savetxt(os.path.join(folder_str, spectra_fname), self.proc_data, delimiter=',',
        #           header=spec_header, comments='', fmt='%f')
        np.savetxt(os.path.join(folder_str, fit_fname), save_results, delimiter=',',
                   header=fits_header, comments='', fmt='%0.6e')

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
    return (A / np.pi) * (gamma / ((x - mu)**2 + gamma**2)) + y0

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

def hsg_sum_spectra(object_list, do_fvb_crr=False):
    """
    This function will add all the things that should be added.  Obvs.  It will
    also calculate the standard error of the mean for every NIR frequency.  The
    standard error is the sum of the dark noise and the "shot" noise.

    Also, remember, we're adding together and averaging the counts per pulse.  
    Hence why we use self.proc_data after it has been divided by the number of 
    fel pulses.  

    This function should not need to be called anymore.  This should be handled
    in the GUI.

    object_list: A list of spectrum objects
    """
    print "I'm trying!"

    good_list = []
    for index in xrange(len(object_list)):
        # dark_var = 0
        num_images = 1
        try:
            temp = object_list.pop(0)
            # var_holder = np.array(temp.proc_data[:, 2])**2
            crr_holder = np.array(temp.proc_data[:, 1]).reshape((1600, 1))

        except Exception as E:
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

                    crr_holder = np.hstack((crr_holder, spec.proc_data[:, 1].reshape((1600,1))))

                    object_list.remove(spec)
        if do_fvb_crr and num_images > 1:
            print "\nI am doing cosmic ray removal!!\n"
            crr_holder = fvb_crr(crr_holder, debugging=False)
            temp.proc_data[:, 1] = np.mean(crr_holder, axis=1)
        else:
            print "\nI did not do anything bad, I think\n"
        # std_error = np.sqrt(np.var(stderr_holder, axis=1, dtype=np.float64) + dark_var) / np.sqrt(spec_number) # Checking some sigma stuff from curve_fit
        # This standard error is for every point.  I think it actually overestimates
        # the error at places with no signal because we add the dark variance
        # effectively twice.
        #print "final dark_stdev:", np.sqrt(dark_var)
        # temp.add_std_error(std_error)
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

    for elem in spectra_list:
        print "Spec_step is", elem.parameters["spec_step"]
    for index in xrange(len(spectra_list)):
        try:
            temp = spectra_list.pop(0)
            print "\nStarting with this guy", temp.parameters["spec_step"], "\n"
        except:
            break

        good_list.append(FullHighSideband(temp))

        counter = temp.parameters["spec_step"] + 1
        print "Starting counter is", counter
        temp_list = list(spectra_list)
        for piece in temp_list:
            print "checking this spec_step", piece.parameters["spec_step"]
            print "The counter is", counter
            if temp.parameters["series"] == piece.parameters["series"]:
                if piece.parameters["spec_step"] == counter:
                    print "I found this one", piece.parameters["series"]
                    counter += 1
                    good_list[-1].add_CCD(piece)
                    spectra_list.remove(piece)                    
        good_list[-1].make_results_array()
    return good_list

def pmt_sorter(folder_path):
    """
    This function will be fed a folder with a bunch of PMT data files in it.  
    The folder should contain a bunch of spectra with at least one sideband in
    them, each differing by the series entry in the parameters dictionary.  

    This function will return a list of HighSidebandPMT objects.  
    """
    file_list = glob.glob(os.path.join(folder_path, '*[0-9].txt'))

    pmt_list = []
    for sb_file in file_list:
        with open(sb_file, 'rU') as f:
            sb_num = int(f.readline().split(' ')[-1])
            parameters_str = f.readline()
            parameters = json.loads(parameters_str[1:])
        try:
            for pmt_spectrum in pmt_list: # pmt_spectrum is a pmt object?
                if parameters['series'] == pmt_spectrum.parameters['series']:
                    pmt_spectrum.add_sideband(sb_file)
                    break
            else: # this will execute IF the break was NOT called
                pmt_list.append(HighSidebandPMT(sb_file))
        except:
            pmt_list.append(HighSidebandPMT(sb_file))

    for pmt_spectrum in pmt_list:
        pmt_spectrum.process_sidebands()
    return pmt_list

def stitch_abs_results(main, new):
    raise NotImplementedError

####################
# Helper functions 
####################

def fvb_crr(raw_array, offset=0, medianRatio=1, noiseCoeff=5, debugging=False):
    """

        Remove cosmic rays from a sequency of identical exposures
        :param raw_array: The array to be cleaned. Successive spectra should
                be the columns (i.e. 1600 x n) of the raw_array
        :param offset: baseline to add to raw_array.
               Not used, but here if it's needed in the future
        :param medianRatio: Multiplier to the median when deciding a cutoff
        :param noiseCoeff: Multiplier to the noise on the median
                    May need changing for noisy data
        :return:
    """

    d = np.array(raw_array)

    med = ndimage.filters.median_filter(d, size=(1, d.shape[1]), mode='wrap')
    med = np.median(d, axis=1).reshape(d.shape[0], 1)
    if debugging:
        print "shape of median filter:", med.shape
    meanMedian  = med.mean(axis=1)
    # meanMedian = med.copy()
    if debugging:
        print "shape of meaned median filter:", meanMedian.shape
    # Construct a cutoff for each pixel. It was kind of guess and
    # check
    cutoff = meanMedian * medianRatio + noiseCoeff * np.std(meanMedian[-100:])
    if debugging:
        print "shape of cutoff criteria:", cutoff.shape
        import pyqtgraph as pg

        winlist = []
        app = pg.QtGui.QApplication([])

        win = pg.GraphicsLayoutWidget()
        win.setWindowTitle("Raw Image")
        p1 = win.addPlot()

        img = pg.ImageItem()
        img.setImage(d.copy().T)
        p1.addItem(img)

        hist = pg.HistogramLUTItem()
        hist.setImageItem(img)
        win.addItem(hist)

        win.nextRow()
        p2 = win.addPlot(colspan=2)
        p2.setMaximumHeight(250)
        p2.addLegend()
        for i, v in enumerate(d.T):
            p2.plot(v, pen=(i, d.shape[1]), name=str(i))
        p2.plot(np.sum(d, axis=1), pen=pg.mkPen('w', width=3))
        win.show()
        winlist.append(win)

        win2 = pg.GraphicsLayoutWidget()
        win2.setWindowTitle("Median Image")
        p1 = win2.addPlot()

        img = pg.ImageItem()
        img.setImage(med.T)
        p1.addItem(img)

        hist = pg.HistogramLUTItem()
        hist.setImageItem(img)
        win2.addItem(hist)

        win2.nextRow()
        p2 = win2.addPlot(colspan=2)
        p2.setMaximumHeight(250)

        p2.plot(np.sum(med, axis=1)/d.shape[1])
        win2.show()
        winlist.append(win2)



        win2 = pg.GraphicsLayoutWidget()
        win2.setWindowTitle("d-m")
        p1 = win2.addPlot()

        img = pg.ImageItem()
        img.setImage((d - med).T)
        p1.addItem(img)

        hist = pg.HistogramLUTItem()
        hist.setImageItem(img)
        win2.addItem(hist)

        win2.nextRow()
        p2 = win2.addPlot(colspan=2)
        p2.setMaximumHeight(250)
        p2.addLegend()
        for i, v in enumerate((d-med).T):
            p2.plot(v, pen=(i, d.shape[1]), name=str(i))
        p2.plot(cutoff, pen=pg.mkPen('w', width=3))
        win2.show()
        winlist.append(win2)



    # Find the bad pixel positions
    # Note the [:, None] - needed to cast the correct shapes
    badPixs = np.argwhere((d - med)>(cutoff.reshape(len(cutoff), 1)))

    for pix in badPixs:
        # get the other pixels in the row which aren't the cosmic
        if debugging:
            print "cleaning pixel", pix
        p = d[pix[0], [i for i in range(d.shape[1]) if not i==pix[1]]]
        if debugging:
            print "\tRemaining pixels in row are", p
        # Replace the cosmic by the average of the others
        # Could get hairy if more than one cosmic per row.
        # Maybe when doing many exposures?
        d[pix[0], pix[1]] = np.mean(p)

    if debugging:
        win = pg.GraphicsLayoutWidget()
        win.setWindowTitle("Clean Image")
        p1 = win.addPlot()

        img = pg.ImageItem()
        img.setImage(d.copy().T)
        p1.addItem(img)

        hist = pg.HistogramLUTItem()
        hist.setImageItem(img)
        win.addItem(hist)

        win.nextRow()
        p2 = win.addPlot(colspan=2)
        p2.setMaximumHeight(250)
        p2.plot(np.sum(d, axis=1))
        win.show()
        winlist.append(win)
        app.exec_()

    return np.array(d)

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
    if need_ratio:
    # Calculate the appropriate ratio to multiply the new sidebands by.
    # I'm not entirely sure what to do with the error of this guy.
        ratio_list = []
        print '\n1979\nfull[2]', full[0][2]
        new_starter = overlap[-1]
        if len(overlap) > 2:
            overlap = [x for x in overlap if (x%2 == 0) and (x != min(overlap) and (x != max(overlap)))]
        for sb in overlap:
            ratio_list.append(full[sb][2] / new_dict[sb][2])
        ratio = np.mean(ratio_list)
        error = np.std(ratio_list) / np.sqrt(len(ratio_list))
        print "Ratio list", ratio_list
        print "Ratio", ratio
        print "Error", error
        print '\n1990\nfull[2]', full[0][2]
    # Adding the new sidebands to the full set and moving errors around.
    # I don't know exactly what to do about the other aspecs of the sidebands
    # besides the strength and its error.
        for sb in overlap:
            full[sb][2] = ratio * new_dict[sb][2]
            full[sb][3] = full[sb][2] * np.sqrt((error / ratio)**2 + (new_dict[sb][3] / new_dict[sb][2])**2)
            print '\n1997\nfull[2]', full[0][3]
            # Now for linewidths
            lw_error = np.sqrt(full[sb][5]**(-2) + new_dict[sb][5]**(-2))**(-1)
            lw_avg = (full[sb][4] / (full[sb][5]**2) + new_dict[sb][4] / (new_dict[sb][5]**2)) / (full[sb][5]**(-2) + new_dict[sb][5]**(-2))
            full[sb][4] = lw_avg
            full[sb][5] = lw_error
        print '\n2003\nfull[2]', full[0][2]
    else:
        new_starter = overlap[-1]
        overlap = [x for x in overlap if (x%2 == 0) and (x != min(overlap) and (x != max(overlap)))]
        for sb in overlap:
            if verbose:
                print "The sideband", sb
                print "Old value", full[sb][4]*1000
                print "Add value", new_dict[sb][4]*1000
            error = np.sqrt(full[sb][3]**(-2) + new_dict[sb][3]**(-2))**(-1)
            avg = (full[sb][2] / (full[sb][3]**2) + new_dict[sb][2] / (new_dict[sb][3]**2)) / (full[sb][3]**(-2) + new_dict[sb][3]**(-2))
            full[sb][2] = avg
            full[sb][3] = error
            
            lw_error = np.sqrt(full[sb][5]**(-2) + new_dict[sb][5]**(-2))**(-1)
            lw_avg = (full[sb][4] / (full[sb][5]**2) + new_dict[sb][4] / (new_dict[sb][5]**2)) / (full[sb][5]**(-2) + new_dict[sb][5]**(-2))
            full[sb][4] = lw_avg
            full[sb][5] = lw_error # This may not be the exactly right way to calculate the error
            if verbose:
                print "New value", lw_avg * 1000
    if need_ratio:
        print '\n2024\nfull[2]', full[0][2]
    for sb in [x for x in new_dict.keys() if (x >= new_starter)]:
        full[sb] = new_dict[sb]
        if need_ratio:
            full[sb][2] = ratio * full[sb][2]
            full[sb][3] = full[sb][2] * np.sqrt((error / ratio)**2 + (ratio * full[sb][3] / full[sb][2])**2)
            print '\n2030\nfull[2]', full[0][2]
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
    '''
    param_array_norm = np.array(param_array).T # python iterates over rows
    for elem in [x for x in xrange(len(param_array_norm)) if (x-1)%7 == 3]:
        temp_max = np.max(param_array_norm[elem])
        param_array_norm[elem] = param_array_norm[elem] / temp_max
        param_array_norm[elem + 1] = param_array_norm[elem + 1] / temp_max
    '''
    snipped_array = param_array[:, 0]
    norm_array = param_array[:, 0]
    print "Snipped_array is", snipped_array
    for ii in xrange(len(param_array.T)):
        if (ii - 1) % 7 == 0:
            print "param_array shape", param_array[:, ii]
            snipped_array = np.vstack((snipped_array, param_array[:, ii]))
            norm_array = np.vstack((norm_array, param_array[:, ii]))
        elif (ii - 1) % 7 == 1:
            snipped_array = np.vstack((snipped_array, param_array[:, ii]))
            norm_array = np.vstack((norm_array, param_array[:, ii]))
        elif (ii - 1) % 7 == 3:
            snipped_array = np.vstack((snipped_array, param_array[:, ii]))

            temp_max = np.max(param_array[:, ii])
            norm_array = np.vstack((norm_array, param_array[:, ii] / temp_max))
        elif (ii - 1) % 7 == 4:
            snipped_array = np.vstack((snipped_array, param_array[:, ii]))
            norm_array = np.vstack((norm_array, param_array[:, ii] / temp_max))

    snipped_array = snipped_array.T
    norm_array = norm_array.T


    try:
        os.mkdir(folder_str)
    except OSError, e:
        if e.errno == errno.EEXIST:
            pass
        else:
            raise
    norm_name = file_name + '_norm.txt'
    snip_name = file_name + '_snip.txt'
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
        origin_import2 += ",order,eV,,arb. u.,,meV,"
        origin_import3 += ",,{0},,{0},,{0},".format(order)
    origin_total = origin_import1 + "\n" + origin_import2 + "\n" + origin_import3

    origin_import1 = param_name
    origin_import2 = unit
    origin_import3 = ""
    for order in sb_included:
        origin_import1 += ",Sideband,Frequency,Sideband strength,error"
        origin_import2 += ",order,eV,arb. u.,"
        origin_import3 += ",,{0},{0},".format(order)
    origin_snip = origin_import1 + "\n" + origin_import2 + "\n" + origin_import3

    header_total = '#' + included_spectra_str + '\n' + origin_total
    header_snip = '#' + included_spectra_str + '\n' + origin_snip

    #print "Spec header: ", spec_header
    print "the param_array is:", param_array
    np.savetxt(os.path.join(folder_str, file_name), param_array, delimiter=',', 
               header=header_total, comments='', fmt='%0.6e')
    np.savetxt(os.path.join(folder_str, snip_name), snipped_array, delimiter=',', 
               header=header_snip, comments='', fmt='%0.6e')
    np.savetxt(os.path.join(folder_str, norm_name), norm_array, delimiter=',', 
               header=header_snip, comments='', fmt='%0.6e')
    print "Saved the file.\nDirectory: {}".format(os.path.join(folder_str, file_name))

####################
# Complete functions 
####################

def proc_n_plotPMT(folder_path, plot=False, save=None, verbose=False):
    """
    This function will take a pmt object, process it completely.
    """
    pmt_data = pmt_sorter(folder_path)

    for spectrum in pmt_data:
        spectrum.integrate_sidebands()
        spectrum.laser_line()
        print spectrum.full_dict
        index = 0
        if plot:
            plt.figure('PMT data')
            for elem in spectrum.sb_dict.values():
                plt.errorbar(elem[:, 0], elem[:, 1], elem[:, 2], marker='o')
            plt.figure('Sideband strengths')
            plt.errorbar(spectrum.sb_results[:, 0], spectrum.sb_results[:, 3], spectrum.sb_results[:, 4], label=spectrum.parameters['series'], marker='o')
        if type(save) is tuple:
            spectrum.save_processing(save[0], save[1], index=index)
    plt.legend()
    return pmt_data

def proc_n_plotCCD(folder_path, cutoff=5, offset=None, plot=False, save=None, verbose=False):
    """
    This function will take a list of ccd files and process it completely.
    save_name is a tuple (file_base, folder_path)
    """
    file_list = glob.glob(os.path.join(folder_path, '*seq_spectrum.txt'))
    raw_list = []
    for fname in file_list:
        raw_list.append(HighSidebandCCD(fname, spectrometer_offset=offset))

    index = 0
    spectra_list = hsg_sum_spectra(raw_list, do_fvb_crr=True)
    for spectrum in spectra_list:
        spectrum.guess_sidebands(cutoff=cutoff, verbose=verbose)
        spectrum.fit_sidebands(plot=plot, verbose=verbose)
        if plot:
            plt.figure('CCD data')
            plt.errorbar(spectrum.proc_data[:, 0], spectrum.proc_data[:, 1], spectrum.proc_data[:, 2], label=spectrum.parameters['series'])
            plt.legend()
            #plt.yscale('log')
            plt.figure('Sideband strengths')
            plt.errorbar(spectrum.sb_results[:, 0], spectrum.sb_results[:, 3], spectrum.sb_results[:, 4], label=spectrum.parameters['series'], marker='o')
            plt.legend()    
            plt.yscale('log')
        if type(save) is tuple:
            spectrum.save_processing(save[0], save[1], marker=spectrum.parameters["series"] + '_' + str(spectrum.parameters["spec_step"]), index=index)
            index += 1
    return spectra_list


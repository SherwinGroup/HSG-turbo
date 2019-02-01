# -*- coding: utf-8 -*-
"""
Created on Fri Jun 19 13:16:25 2015

@author: hbanks

Brevity required, prurience preferred
"""


import os
import io
import glob
import errno
import copy
import json
import warnings
import numpy as np
from scipy.optimize import curve_fit
import scipy.interpolate as spi
import scipy.optimize as spo
import scipy.fftpack as fft
import matplotlib.pyplot as plt
import scipy.ndimage as ndimage
import itertools as itt
np.set_printoptions(linewidth=500)

# One of the main results is the HighSidebandCCD.sb_results array. These are the
# various mappings between index and real value
# I deally, this code should be converted to pandas to avoid this issue,
# but that's outside the scope of current work.
# [sb number, Freq (eV), Freq error (eV), Gauss area (arb.), Area error, Gauss linewidth (eV), Linewidth error (eV)]
# [    0    ,      1   ,        2,      ,        3         ,      4    ,         5           ,        6            ]
class sbarr(object):
    SBNUM = 0
    CENFREQ = 1
    CENFREQERR = 2
    AREA = 3
    AREAERR = 4
    WIDTH = 5
    WIDTHERR = 6


####################
# Objects
####################

class CCD(object):
    def __init__(self, fname, spectrometer_offset=None):
        """
        This will read the appropriate file and make a basic CCD object.  Fancier
        things will be handled with the sub classes.

        Creates:
        self.parameters = Dictionary holding all of the information from the
                          data file, which comes from the JSON encoded header in the data
                          file
        self.description = string that is the text box from data taking GUI
        self.raw_data = raw data output by measurement software, wavelength vs.
                        data, errors.  There may be text for some of the entries
                        corresponding to text used for Origin imports, but they
                        should appear as np.nan
        self.ccd_data = semi-processed 1600 x 3 array of photon energy vs. data with standard error of mean at that pixel
                        calculated by taking multiple images. Standard error is calculated from
                        the data collection software

        Most subclasses should make a self.proc_data, which will do whatever
        processing is required to the ccd_data, such as normalizing, taking ratios,
        etc.

        :param fname: file name where the data is saved
        :type fname: str
        :param spectrometer_offset: if the spectrometer won't go where it's told, use this to correct the wavelengths (nm)
        :type spectrometer_offset: float
        """

        self.fname = fname

        # Checking restrictions from Windows path length limits. Check if you can
        # open the file:
        try:
            with open(fname) as f: pass
        except FileNotFoundError:
            # Couldn't find the file. Could be you passed the wrong one, but I'm
            # finding with a large number of subfolders for polarimetry stuff,
            # you end up exceeding Windows'  filelength limit.
            # Haven't tested on Mac or UNC moutned drives (e.g \\128.x.x.x\Sherwin\)
            fname = r"\\?\\" + os.path.abspath(fname)

        # Read in the JSON-formatted parameter string.
        # The lines are all prepended by '#' for easy numpy importing
        # so loop over all those lines
        with open(fname, 'r') as f:
            param_str = ''
            line = f.readline()
            while line[0] == '#':
                ### changed 09/17/18
                # This line assumed there was a single '#'
                # param_str += line[1:]
                # while this one handles everal (because I found old files
                # which had '## <text>...'
                param_str += line.replace("#", "")
                line = f.readline()
            # Parse the JSON string
            try:
                self.parameters = json.loads(param_str)
            except json.JSONDecodeError:
                # error from _really_ old data where comments were dumped after a
                # single-line json dumps
                self.parameters=json.loads(param_str.splitlines()[0])

        # Spec[trometer] steps are set to define the same physical data, but taken at
        # different spectrometer center wavelengths. This value is used later
        # for stitching these scans together
        try:
            self.parameters["spec_step"] = int(self.parameters["spec_step"])
        except (ValueError, KeyError):
            # If there isn't a spe
            self.parameters["spec_step"] = 0

        # Slice through 3 to get rid of comments/origin info.
        # Would likely be better to check np.isnan() and slicing out those nans.
        # I used flipup so that the x-axis is an increasing function of frequency
        self.raw_data = np.flipud(np.genfromtxt(fname, comments='#', delimiter=',')[3:])


        # The camera chip is 1600 pixels wide. This line was redudent with the [3:]
        # slice above and served to make sure there weren't extra stray bad lines
        # hanging around.
        #
        # This should also be updated some day to compensate for any horizontal bining
        # on the chip, or masking out points that are bad (cosmic ray making it
        # through processing, room lights or monitor lines interfering with signal)
        self.ccd_data = np.array(self.raw_data[:1600, :])

        # Check to see if the spectrometer offset is set. This isn't specified
        # during data collection. This is a value that can be appended
        # when processing if it's realized the data is offset.
        # This allows the offset to be specified and kept with the data file itself,
        # instead of trying to do it in individual processing scripts
        #
        # It's allowed as a kwarg parameter in this script for trying to determine
        # what the correct offset should be
        if spectrometer_offset is not None or "offset" in self.parameters:
            try:
                self.ccd_data[:, 0] += float(self.parameters["offset"])
            except:
                self.ccd_data[:, 0] += spectrometer_offset

        # Convert from nm to eV
        # self.ccd_data[:, 0] = 1239.84 / self.ccd_data[:, 0]
        self.ccd_data[:, 0] = photon_converter["nm"]["eV"](self.ccd_data[:, 0])


class Photoluminescence(CCD):
    def __init__(self, fname):
        """
        This object handles PL-type data. The only distinction from the parent class 
        is that the CCD data gets normalized to the exposure time to make different 
        exposures directly comparable.

        creates:
        self.proc_data = self.ccd_data divided by the exposure time
                         units: PL counts / second
        :param fname: name of the file
        :type fname: str

        """
        super(Photoluminescence, self).__init__(fname)

        # Create a copy of the array , and then normalize the signal and the errors
        # by the exposure time
        self.proc_data = np.array(self.ccd_data) 
        self.proc_data[:, 1] = self.proc_data[:, 1] / self.parameters['exposure']
        self.proc_data[:, 2] = self.proc_data[:, 2] / self.parameters['exposure']


class Absorbance(CCD):
    def __init__(self, fname):
        """
        There are several ways Absorbance data can be loaded
        You could try to load the abs data output from data collection directly,
        which has the wavelength, raw, blank and actual absorbance data itself.
        This is best way to do it.

        Alternatively, you could want to load the raw transmission/reference
        data, ignoring (or maybe not even having) the abs calculated
        from the data collection software. If you want to do it this way,
        you should pass fname as a list where the first element is the
        file name for the reference data, and the second is the absorbance data
        At first, it didn't really seem to make sense to let you pass just the
        raw reference or raw abs data,


        Creates:
        self.ref_data = np array of the reference,
            freq (eV) vs. reference (counts)
        self.raw_data = np.array of the raw absorption spectrum,
            freq (eV) vs. reference (counts)
        self.proc_data = np.array of the absorption spectrum
            freq (eV) vs. "absorbance" (dB)

        Note, the error bars for this data haven't been defined.

        :param fname: either an absorbance filename, or a length 2 list of filenames
        :type fname: str
        :return: None
        """
        if "abs_" in fname:
            super(Absorbance, self).__init__(fname)
            # Separate into the separate data sets
            #   The raw counts of the reference data
            self.ref_data = np.array(self.ccd_data[:, [0, 1]])
            #   Raw counts of the sample
            self.raw_data = np.array(self.ccd_data[:, [0, 2]])
            #   The calculated absorbance data (-10*log10(raw/ref))
            self.proc_data = np.array(self.ccd_data[:, [0, 3]]) # Already in dB's

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
                # only loading first two columns to avoid numpy trying
                # to parse all of the data

                # See CCD.__init__ for what's going on.

                self.ref_data = np.flipud(np.genfromtxt(fname[0], comments='#',
                                                        delimiter=',', usecols=(0, 1)))

                self.ref_data = np.array(self.ref_data[:1600, :])
                self.ref_data[:, 0] = 1239.84 / self.ref_data[:, 0]

                self.raw_data = np.flipud(np.genfromtxt(fname[1], comments='#',
                                                        delimiter=',', usecols=(0, 1)))

                self.raw_data = np.array(self.raw_data[:1600, :])
                self.raw_data[:, 0] = 1239.84 / self.raw_data[:, 0]
            except Exception as e:
                print("Exception opening absorbance data,", e)

            # Calculate the absorbance from the raw camera counts.
            self.proc_data = np.empty_like(self.ref_data)
            self.proc_data[:, 0] = self.ref_data[:, 0]
            self.proc_data[:, 1] = -10*np.log10(self.raw_data[:, 1] / self.ref_data[:,
                                                                     1])

    def abs_per_QW(self, qw_number):
        """

        :param qw_number: number of quantum wells in the sample.
        :type qw_number: int
        :return: None
        """
        """
        This method turns the absorption to the absorbance per quantum well.  Is
        that how this data should be reported?

        Also, I'm not sure if columns 1 and 2 are correct.
        """
        temp_abs = -np.log(self.proc_data[:, 1] / self.proc_data[:, 2]) / qw_number
        self.proc_data = np.hstack((self.proc_data, temp_abs))

    def fft_smooth(self, cutoff, inspectPlots=False):
        """
        This function removes the Fabry-Perot that affects the absorption data

        creates:
        self.clean = np.array of the Fourier-filtered absorption data, freq (eV) vs. absorbance (dB!)
        self.parameters['fourier cutoff'] = the low pass cutoff frequency, in eV**(-1)
        :param cutoff: Fourier frequency of the cut off for the low pass filter
        :type cutoff: int or float
        :param inspectPlots: Do you want to see the results?
        :type inspectPlots: bool
        :return: None
        """
        # self.fixed = -np.log10(abs(self.raw_data[:, 1]) / abs(self.ref_data[:, 1]))
        # self.fixed = np.nan_to_num(self.proc_data[:, 1])
        # self.fixed = np.column_stack((self.raw_data[:, 0], self.fixed))
        self.parameters['fourier cutoff'] = cutoff
        self.clean = low_pass_filter(self.proc_data[:, 0], self.proc_data[:, 1], cutoff, inspectPlots)

    def save_processing(self, file_name, folder_str, marker='', index=''):
        """
        This bad boy saves the absorption spectrum that has been manipulated.

        Saves 100 lines of comments.

        :param file_name: The base name of the file to be saved
        :type file_name: str
        :param folder_str: The name of the folder where the file will be saved
        :type folder_str: str
        :param marker: A further label that might be the series tag or something
        :type marker: str
        :param index: If multiple files are being saved with the same name, include an integer to append to the end of the file
        :type index: int
        :return: None
        """
        try:
            os.mkdir(folder_str)
        except OSError as e:
            if e.errno == errno.EEXIST:
                pass
            else:
                raise

        spectra_fname = file_name + '_' + marker + '_' + str(index) + '.txt'
        self.save_name = spectra_fname

        try:
            parameter_str = json.dumps(self.parameters, sort_keys=True, indent=4, separators=(',', ': '))
        except:
            print("Source: EMCCD_image.save_images\nJSON FAILED")
            print("Here is the dictionary that broke JSON:\n", self.parameters)
            return
        parameter_str = parameter_str.replace('\n', '\n#')

        num_lines = parameter_str.count('#')  # Make the number of lines constant so importing into Origin is easier
        # for num in range(99 - num_lines): parameter_str += '\n#'
        parameter_str += '\n#' * (99 - num_lines)

        origin_import_spec = '\nNIR frequency,Signal,Standard error\neV,arb. u.,arb. u.'
        spec_header = '#' + parameter_str + origin_import_spec
        # spec_header = '#' + parameter_str + '\n#' + self.description[:-2] + origin_import_spec

        np.savetxt(os.path.join(folder_str, spectra_fname), self.proc_data, delimiter=',',
                   header=spec_header, comments='', fmt='%0.6e')
        spectra_fname = 'clean ' + spectra_fname
        np.savetxt(os.path.join(folder_str, spectra_fname), self.clean, delimiter=',',
                   header=spec_header, comments='', fmt='%0.6e')
        print("Save image.\nDirectory: {}".format(os.path.join(folder_str, spectra_fname)))

# class LaserLineCCD(HighSidebandCCD):
#     """
#     Class for use when doing alinging/testing by sending the laser
#     directly into the CCD. Modifies how "sidebands" and guess and fit,
#     simply looking at the max signal.
#     """
#     def guess_sidebands(self, cutoff=8, verbose=False, plot=False):
#         pass

class NeonNoiseAnalysis(CCD):
    """
    This class is used to make handling neon calibration lines easier.  It's not great.
    """
    def __init__(self, fname, spectrometer_offset=None):
        # print 'opening', fname
        super(NeonNoiseAnalysis, self).__init__(fname, spectrometer_offset=spectrometer_offset)

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
        print('\n\n')

        self.ccd_data = np.flipud(self.ccd_data)
        # self.high_noise_region = np.array(self.ccd_data[30:230, :])
        self.high_noise_region = np.array(self.ccd_data[80:180, :])  # for dark current measurements
        self.low_noise_region1 = np.array(self.ccd_data[380:700, :])
        self.low_noise_region2 = np.array(self.ccd_data[950:1200, :])
        self.low_noise_region3 = np.array(self.ccd_data[1446:1546, :])

        # self.high_noise = np.std(self.high_noise_region[:, 1])
        self.high_noise_std = np.std(self.high_noise_region[:, 1])
        self.high_noise_sig = np.mean(self.high_noise_region[:, 1])
        self.low_noise1 = np.std(self.low_noise_region1[:, 1])
        self.low_noise2 = np.std(self.low_noise_region2[:, 1])
        self.low_noise_std = np.std(self.low_noise_region2[:, 1])
        self.low_noise_sig = np.mean(self.low_noise_region2[:, 1])
        self.low_noise3 = np.std(self.low_noise_region3[:, 1])

        # self.noise_list = [self.high_noise, self.low_noise1, self.low_noise2, self.low_noise3]

        self.peak1 = np.array(self.ccd_data[303:323, :])
        self.peak2 = np.array(self.ccd_data[319:339, :])
        self.peak3 = np.array(self.ccd_data[736:746, :])
        self.peak4 = np.array(self.ccd_data[1268:1288, :])
        self.peak5 = np.array(self.ccd_data[1381:1421, :])

        temp_max = np.argmax(self.peak1[:, 1])
        self.signal1 = np.sum(self.peak1[temp_max - 1:temp_max + 2, 1])
        self.error1 = np.sqrt(np.sum(self.peak1[temp_max - 1:temp_max + 2, 2] ** 2))

        temp_max = np.argmax(self.peak2[:, 1])
        self.signal2 = np.sum(self.peak2[temp_max - 1:temp_max + 2, 1])
        self.error2 = np.sqrt(np.sum(self.peak2[temp_max - 1:temp_max + 2, 2] ** 2))

        temp_max = np.argmax(self.peak3[:, 1])
        self.signal3 = np.sum(self.peak3[temp_max - 1:temp_max + 2, 1])
        self.error3 = np.sqrt(np.sum(self.peak3[temp_max - 1:temp_max + 2, 2] ** 2))

        temp_max = np.argmax(self.peak4[:, 1])
        self.signal4 = np.sum(self.peak4[temp_max - 1:temp_max + 2, 1])
        self.error4 = np.sqrt(np.sum(self.peak4[temp_max - 1:temp_max + 2, 2] ** 2))

        temp_max = np.argmax(self.peak5[:, 1])
        self.signal5 = np.sum(self.peak5[temp_max - 1:temp_max + 2, 1])
        self.error5 = np.sqrt(np.sum(self.peak5[temp_max - 1:temp_max + 2, 2] ** 2))

        self.signal_list = [self.signal1, self.signal2, self.signal3, self.signal4, self.signal5]
        self.error_list = [self.error1, self.error2, self.error3, self.error4, self.error5]
        print("Signal list:", self.signal_list)
        self.ccd_data = np.flipud(self.ccd_data)

    def process_stuff(self):
        """
        This one puts high_noise, low_noise1, signal2, and error2 in a nice horizontal array
        """
        # self.results = np.array([self.high_noise, self.low_noise1, self.signal5, self.error5])
        # average = np.mean([self.low_noise1, self.low_noise2, self.low_noise3])
        # self.results = np.array([self.high_noise, self.low_noise1, self.low_noise2, self.low_noise3, self.high_noise/average])
        self.results = np.array([self.high_noise_sig, self.high_noise_std, self.low_noise_sig, self.low_noise_std])

def collect_noise(neon_list, param_name, folder_name, file_name, name='Signal'):
    """
    This function acts like save parameter sweep.

    param_name = string that we're gonna save!
    """
    # param_array = None
    for elem in neon_list:
        print("pname: {}".format(elem.parameters[param_name]))
        print("results:", elem.results)
        temp = np.insert(elem.results, 0, elem.parameters[param_name])
        try:
            param_array = np.row_stack((param_array, temp))
        except UnboundLocalError:
            param_array = np.array(temp)

    if len(param_array.shape) == 1:
        print("I don't think you want this file")
        return
        # append the relative peak error

    print('\n', param_array, '\n')

    param_array = np.column_stack((param_array, param_array[:, 4] / param_array[:, 3]))
    # append the snr
    param_array = np.column_stack((param_array, param_array[:, 3] / param_array[:, 2]))

    try:
        param_array = param_array[param_array[:, 0].argsort()]
    except:
        print("param_array shape", param_array.shape)
        raise

    try:
        os.mkdir(folder_name)
    except OSError as e:
        if e.errno == errno.EEXIST:
            pass
        else:
            raise

    file_name = file_name + '.txt'

    origin_import1 = param_name + ",Noise,Noise,Signal,error,rel peak error,peak signal-to-noise"
    # origin_import1 = param_name + ",Noise,Noise,Noise,Noise,Ratio"
    origin_import2 = ",counts,counts,counts,counts,,"
    # origin_import2 = ",counts,counts,counts,,"
    origin_import3 = ",High noise region,Low noise region,{},{} error,{} rel error, {}".format(name, name, name, name)
    # origin_import3 = ",High noise region,Low noise region 1,Low noise region 2,Low noise region 3,High/low"
    header_total = origin_import1 + "\n" + origin_import2 + "\n" + origin_import3

    # print "Spec header: ", spec_header
    print("the param_array is:", param_array)
    np.savetxt(os.path.join(folder_name, file_name), param_array, delimiter=',',
               header=header_total, comments='', fmt='%0.6e')
    print("Saved the file.\nDirectory: {}".format(os.path.join(folder_name, file_name)))

class HighSidebandCCD(CCD):
    def __init__(self, hsg_thing, parameter_dict=None, spectrometer_offset=None):
        """
        This will read the appropriate file.  The header needs to be fixed to
        reflect the changes to the output header from the Andor file.  Because
        another helper file will do the cleaning and background subtraction,
        those are no longer part of this init.  This also turns all wavelengths
        from nm (NIR ones) or cm-1 (THz ones) into eV.

        OR, if an array is thrown in there, it'll handle the array and dict

        Input:
            For post-processing analysis:
                hsg_thing = file name of the hsg spectrum from CCD superclass
                spectrometer_offset = number of nanometers the spectrometer is off by,
                                      should be 0.0...but can be 0.2 or 1.0
            For Live-software:
                hsg_thing = np array of spectrum from camera
                parameter_dict = equipment dict generated by software

        Internal:
        self.hsg_thing = the filename
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

        :param hsg_thing: file name for the file to be opened.  OR the actually hsg np.ndarray.  Fun!
        :type hsg_thing: str OR np.ndarray
        :param parameter_dict: If being loaded through the data acquisition GUI, throw the dict in here
        :type parameter_dict: dict
        :param spectrometer_offset: Number of nm the spectrometer is off by
        :type spectrometer_offset: float
        :return: None, technically
        """
        if isinstance(hsg_thing, str):
            super(HighSidebandCCD, self).__init__(hsg_thing, spectrometer_offset=spectrometer_offset)
            # TODO: fix addenda bullshit
            self.addenda = []
            self.subtrahenda = []
        elif isinstance(hsg_thing, np.ndarray):
            self.parameters = parameter_dict.copy()  # Probably shouldn't shoehorn this in this way
            self.addenda = []
            self.subtrahenda = []
            self.ccd_data = np.array(hsg_thing)
            self.ccd_data[:, 0] = 1239.84 / self.ccd_data[:, 0]
            # This data won't have an error column, so attached a column of ones
            self.ccd_data = np.column_stack((self.ccd_data, np.ones_like(self.ccd_data[:,1])))
            self.ccd_data = np.flipud(self.ccd_data) # Because turning into eV switches direction
            self.fname = "Live Data"
        else:
            raise Exception("I don't know what this file type is {}, type: {}".format(
                hsg_thing, type(hsg_thing)
            ))
        self.proc_data = np.array(self.ccd_data)
        # proc_data is now a 1600 long array with [frequency (eV), signal (counts / FEL pulse), S.E. of signal mean]

        # self.parameters["nir_freq"] = 1239.84 / float(self.parameters["nir_lambda"])
        self.parameters["nir_freq"] = 1239.84 / float(self.parameters.get("nir_lambda", -1))
        # self.parameters["thz_freq"] = 0.000123984 * float(self.parameters["fel_lambda"])
        self.parameters["thz_freq"] = 0.000123984 * float(self.parameters.get("fel_lambda", -1))
        # self.parameters["nir_power"] = float(self.parameters["nir_power"])
        self.parameters["nir_power"] = float(self.parameters.get("nir_power", -1))
        try: # This is the new way of doing things.  Also, now it's power
            self.parameters["thz_energy"] = float(self.parameters["pulseEnergies"]["mean"])
            self.parameters["thz_energy_std"] = float(self.parameters["pulseEnergies"]["std"])
        except: # This is the old way TODO: DEPRECATE THIS
            self.parameters["thz_energy"] = float(self.parameters.get("fel_power", -1))

        # things used in fitting/guessing
        self.sb_list = np.array([])
        self.sb_index = np.array([])
        self.sb_dict = {}
        self.sb_results = np.array([])
        self.full_dict = {}

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

        This raises a FutureWarning because these were designed early on and
        haven't been used much.

        :param other: The thing to be added, it's either a int/float or a HighSidebandCCD object
        :type other: int/float or HighSidebandCCD
        :return: Sum of self and other
        :rtype: HighSidebandCCD
        """
        raise FutureWarning
        ret = copy.deepcopy(self)
        # Add a constant offset to the data
        if type(other) in (int, float):
            ret.proc_data[:, 1] = self.proc_data[:, 1] + other
            ret.addenda[0] = ret.addenda[0] + other

        # or add the data of two hsg_spectra together
        else:
            if np.isclose(ret.parameters['center_lambda'], other.parameters['center_lambda']):
                ret.proc_data[:, 1] = self.proc_data[:, 1] + other.proc_data[:, 1]
                ret.proc_data[:, 2] = np.sqrt(self.proc_data[:, 1] ** 2 + other.proc_data[:, 1] ** 2)
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

        This raises a FutureWarning because these were designed early on and
        haven't been used much.

        :param other: The thing to be subtracted, it's either a int/float or a HighSidebandCCD object
        :type other: int/float or HighSidebandCCD
        :return: Sum of self and other
        :rtype: HighSidebandCCD
        """
        raise FutureWarning
        ret = copy.deepcopy(self)
        # Subtract a constant offset to the data
        if type(other) in (int, float):
            ret.proc_data[:, 1] = self.proc_data[:, 1] - other  # Need to choose a name
            ret.addenda[0] = ret.addenda[0] - other

        # Subtract the data of two hsg_spectra from each other
        else:
            if np.isclose(ret.proc_data[0, 0], other.proc_data[0, 0]):
                ret.proc_data[:, 1] = self.proc_data[:, 1] - other.proc_data[:, 1]
                ret.proc_data[:, 2] = np.sqrt(self.proc_data[:, 1] ** 2 + other.proc_data[:, 1] ** 2)
                ret.subtrahenda.extend(other.addenda[1:])
                ret.addenda.extend(other.subtrahenda)
            else:
                raise Exception('Source: Spectrum.__sub__:\nThese are not from the same grating settings')
        return ret

    def __repr__(self):
        base = """
            fname: {},
            Series: {series},
            spec_step: {spec_step},
            fel_lambda: {fel_lambda},
            nir_lambda: {nir_lambda}""".format(os.path.basename(self.fname),**self.parameters)
        return base

    __str__ = __repr__

    def calc_approx_sb_order(self, test_nir_freq):
        """
        This simple method will simply return a float approximating the order
        of the frequency input.  We need this because the CCD wavelength
        calibration is not even close to perfect.  And it shifts by half a nm
        sometimes.

        :param test_nir_freq: the frequency guess of the nth sideband
        :type test_nir_freq: float
        :return: The approximate order of the sideband in question
        :rtype: float
        """
        nir_freq = self.parameters['nir_freq']
        thz_freq = self.parameters['thz_freq']
        # If thz = 0, prevent error
        if not thz_freq: thz_freq = 1
        approx_order = (test_nir_freq - nir_freq) / thz_freq
        return approx_order

    def guess_sidebands(self, cutoff=4.5, verbose=False, plot=False, **kwargs):
        """
        Update 05/24/18:
        Hunter had two different loops for negative order sidebands,
        then positive order sidebands. They're done pretty much identically,
        so I've finally merged them into one.

        Finds the locations of all the sidebands in the proc_data array to be
        able to seed the fitting method.  This works by finding the maximum data
        value in the array and guessing what sideband it is.  It creates an array
        that includes this information.  It will then step down, initially by one
        THz frequency, then by twos after it hasn't found any odd ones.  It then
        goes up from the max and finds everything above in much the same way.

        There is currently no rhyme or reason to a cutoff of 8.  I don't know what
        it should be changed to, though.

        Input:
        cutoff = signal-to-noise threshold to count a sideband candidate.

        kwargs:
           window_size: how big of a window (in pixels) to use for checking for
                sidebands. Specified in half-width
              default: 15


        Internal:
        self.sb_list = List of all of the orders the method found
        self.sb_index = index of all of the peaks of the sidebands
        self.sb_guess = three-part list including the frequency, amplitude and
                        error guesses for each sideband
        """
        # TODO: this isn't commented appropriately.  Will it be made more readable first?

        if "cutoff" in self.parameters:
            cutoff = self.parameters["cutoff"]
        else:
            self.parameters['cutoff for guess_sidebands'] = cutoff

        if verbose:
            print("=" * 15)
            print()
            print("Guessing CCD Sideband parameters")
            print(os.path.basename(self.fname))
            print("\tCutoff = {}".format(cutoff))
            print()
            print("=" * 15)
        x_axis = np.array(self.proc_data[:, 0])
        y_axis = np.array(self.proc_data[:, 1])
        try:
            error = np.array(self.proc_data[:, 2])
        except IndexError:
            # Happens on old data where spectra weren't calculated in the live
            # software.
            error = np.ones_like(x_axis)

        min_sb = int(self.calc_approx_sb_order(x_axis[0])) + 1
        try:
            max_sb = int(self.calc_approx_sb_order(x_axis[-1]))
        except ValueError:
            print(x_axis)

        nir_freq = self.parameters["nir_freq"]
        thz_freq = self.parameters["thz_freq"]

        if verbose:
            print("min_sb: {} | max_sb: {}".format(min_sb, max_sb))

        # Find max strength sideband and it's order
        global_max = np.argmax(y_axis)
        order_init = int(round(self.calc_approx_sb_order(x_axis[global_max])))
        # if verbose:
        #     print "The global max is at index", global_max
        if global_max < 15:
            check_y = y_axis[:global_max + 15]
            check_y = np.concatenate((np.zeros(15 - global_max), check_y))
        elif global_max > 1585:
            check_y = y_axis[global_max - 15:]
            check_y = np.concatenate((check_y, np.zeros(global_max - 1585)))
        else:
            check_y = y_axis[global_max - 15:global_max + 15]

        check_max_index = np.argmax(check_y)
        check_max_area = np.sum(check_y[check_max_index - 2:check_max_index + 3])

        check_ave = np.mean(check_y[[0, 1, 2, 3, 4, -1, -2, -3, -4, -5]])
        check_stdev = np.std(check_y[[0, 1, 2, 3, 4, -1, -2, -3, -4, -5]])
        check_ratio = (check_max_area - 3 * check_ave) / check_stdev

        if verbose:
            print(("{:^16}" * 5).format(
                "global_max idx", "check_max_area", "check_ave", "check_stdev",
                "check_ratio"))
            print(("{:^16.5g}" * 5).format(
                global_max, check_max_area, check_ave, check_stdev, check_ratio))

        if check_ratio > cutoff:
            self.sb_list = [order_init]
            self.sb_index = [global_max]
            sb_freq_guess = [x_axis[global_max]]
            sb_amp_guess = [y_axis[global_max]]
            sb_error_est = [
                np.sqrt(sum([i ** 2 for i in error[global_max - 2:global_max + 3]])) / (
                            check_max_area - 5 * check_ave)]
        else:
            print("There are no sidebands in", self.fname)
            raise RuntimeError

        if verbose:
            print("\t Looking for sidebands with f < {:.6f}".format(sb_freq_guess[0]))
        last_sb = sb_freq_guess[0]
        index_guess = global_max
        # keep track of how many consecutive sidebands we've skipped. Sometimes one's
        #  noisy or something, so we want to keep looking after skipping one
        consecutive_null_sb = 0
        consecutive_null_odd = 0
        no_more_odds = False
        break_condition = False
        for order in range(order_init - 1, min_sb - 1, -1):
            # Check to make sure we're not looking at an odd when
            # we've decided to skip them.
            if no_more_odds == True and order % 2 == 1:
                last_sb = last_sb - thz_freq
                if verbose:
                    print("I skipped", order)
                continue

            # Window size to look for next sideband. Needs to be order dependent
            # because higher orders get wider, so we need to look at more.
            # Values are arbitrary.
            window_size = 0.45 + 0.0004 * order  # used to be last_sb?
            lo_freq_bound = last_sb - thz_freq * (
                        1 + window_size)  # Not sure what to do about these
            hi_freq_bound = last_sb - thz_freq * (1 - window_size)

            if verbose:
                print("\nSideband", order)
                print("\t{:.4f} < f_{} < {:.4f}".format(lo_freq_bound, order,
                                                        hi_freq_bound))

            # Get the indices where the energies lie within the bounds for this SB
            sliced_indices = \
            np.where((x_axis > lo_freq_bound) & (x_axis < hi_freq_bound))[0]
            start_index, end_index = sliced_indices.min(), sliced_indices.max()

            # Get a slice of the y_data which is only in the region of interest
            check_y = y_axis[sliced_indices]

            check_max_index = np.argmax(
                check_y)  # This assumes that two floats won't be identical
            # Calculate the "area" of the sideband by looking at the peak value
            # within the range, and the pixel above/below it
            check_max_area = np.sum(check_y[check_max_index - 1:check_max_index + 2])

            if verbose and plot:
                plt.figure("CCD data")
                plt.plot([lo_freq_bound] * 2, [0, check_y[check_max_index]], 'b')
                plt.plot([hi_freq_bound] * 2, [0, check_y[check_max_index]], 'b')
                plt.plot([lo_freq_bound, hi_freq_bound], [check_y[check_max_index]] *
                         2, 'b', label="{} Box".format(order))
                plt.text((lo_freq_bound + hi_freq_bound) / 2, check_y[check_max_index],
                         order)

            # get the slice that doesn't have the peak in it to compare statistics
            check_region = np.append(check_y[:check_max_index - 1],
                                     check_y[check_max_index + 2:])
            check_ave = check_region.mean()
            check_stdev = check_region.std()

            # Calculate an effective SNR, where check_ave is roughly the
            # background level
            check_ratio = (check_max_area - 3 * check_ave) / check_stdev

            if order % 2 == 1:  # This raises the barrier for odd sideband detection
                check_ratio = check_ratio / 1.5
            if verbose:
                print("\t" + ("{:^14}" * 4).format(
                    "check_max_area", "check_ave", "check_stdev", "check_ratio"))
                print("\t" + ("{:^14.5g}" * 4).format(
                    check_max_area, check_ave, check_stdev, check_ratio))

            if check_ratio > cutoff:
                found_index = check_max_index + start_index
                self.sb_index.append(found_index)
                last_sb = x_axis[found_index]

                if verbose:
                    print("I just found", last_sb)

                sb_freq_guess.append(x_axis[found_index])
                sb_amp_guess.append(check_max_area - 3 * check_ave)
                error_est = np.sqrt(
                    sum(
                        [i ** 2 for i in error[found_index - 1:found_index + 2]]
                    )) / (check_max_area - 3 * check_ave)
                if verbose:
                    print("My error estimate is:", error_est)
                sb_error_est.append(error_est)
                self.sb_list.append(order)
                consecutive_null_sb = 0
                if order % 2 == 1:
                    consecutive_null_odd = 0
            else:
                # print "I could not find sideband with order", order
                last_sb = last_sb - thz_freq
                consecutive_null_sb += 1
                if order % 2 == 1:
                    consecutive_null_odd += 1
            if consecutive_null_odd == 1 and no_more_odds == False:
                # print "I'm done looking for odd sidebands"
                no_more_odds = True
            if consecutive_null_sb == 2:
                # print "I can't find any more sidebands"
                break

        # Look for higher sidebands
        if verbose: print("\nLooking for higher energy sidebands")

        last_sb = sb_freq_guess[0]
        index_guess = global_max
        consecutive_null_sb = 0
        consecutive_null_odd = 0
        no_more_odds = False
        break_condition = False
        for order in range(order_init + 1, max_sb + 1):
            if no_more_odds == True and order % 2 == 1:
                last_sb = last_sb + thz_freq
                continue
            window_size = 0.45 + 0.001 * order  # used to be 0.28 and 0.0004
            lo_freq_bound = last_sb + thz_freq * (
                        1 - window_size)  # Not sure what to do about these
            hi_freq_bound = last_sb + thz_freq * (1 + window_size)

            start_index = False
            end_index = False

            if verbose:
                print("\nSideband", order)
                # print "The low frequency bound is", lo_freq_bound
                # print "The high frequency bound is", hi_freq_bound
                print("\t{:.4f} < f_{} < {:.4f}".format(lo_freq_bound, order,
                                                        hi_freq_bound))
            for i in range(index_guess, 1600):
                if start_index == False and i == 1599:
                    # print "I'm all out of space, captain!"
                    break_condition = True
                    break
                elif start_index == False and x_axis[i] > lo_freq_bound:
                    # print "start_index is", i
                    start_index = i
                elif i == 1599:
                    end_index = 1599
                    # print "hit end of data, end_index is 1599"
                elif end_index == False and x_axis[i] > hi_freq_bound:
                    end_index = i
                    # print "end_index is", i
                    index_guess = i
                    break
            if break_condition:
                break
            check_y = y_axis[start_index:end_index]

            check_max_index = np.argmax(
                check_y)  # This assumes that two floats won't be identical
            octant = len(check_y) // 8  # To be able to break down check_y into eighths
            if octant < 1:
                octant = 1

            check_max_area = np.sum(
                check_y[check_max_index - octant - 1:check_max_index + octant + 1])

            if verbose and plot:
                plt.figure("CCD data")
                plt.plot([lo_freq_bound] * 2, [0, check_y[check_max_index]], 'b')
                plt.plot([hi_freq_bound] * 2, [0, check_y[check_max_index]], 'b')
                plt.plot([lo_freq_bound, hi_freq_bound], [check_y[check_max_index]] *
                         2, 'b', label=order)
                plt.text((lo_freq_bound + hi_freq_bound) / 2, check_y[check_max_index],
                         order)

            no_peak = (2 * len(
                check_y)) // 6  # The denominator is in flux, used to be 5
            # if verbose: print "\tcheck_y length", len(check_y)

            check_ave = np.mean(np.take(check_y, np.concatenate(
                (np.arange(no_peak), np.arange(-no_peak, 0)))))
            check_stdev = np.std(np.take(check_y, np.concatenate(
                (np.arange(no_peak), np.arange(-no_peak, 0)))))

            check_ratio = (check_max_area - (2 * octant + 1) * check_ave) / check_stdev

            if verbose:
                print("\tIndices: {}->{} (d={})".format(start_index, end_index,
                                                        len(check_y)))
                # print "check_y is", check_y
                # print "\ncheck_max_area is", check_max_area
                # print "check_ave is", check_ave
                # print "check_stdev is", check_stdev
                # print "check_ratio is", check_ratio

                print("\t" + ("{:^14}" * 4).format(
                    "check_max_area", "check_ave", "check_stdev", "check_ratio"))
                print("\t" + ("{:^14.6g}" * 4).format(
                    check_max_area, check_ave, check_stdev, check_ratio))

            if order % 2 == 1:  # This raises the barrier for odd sideband detection
                check_ratio = check_ratio / 2
            if check_ratio > cutoff:
                found_index = check_max_index + start_index
                self.sb_index.append(found_index)
                last_sb = x_axis[found_index]

                if verbose:
                    print("\tI'm counting this SB at index {} (f={:.4f})".format(
                        found_index, last_sb), end=' ')
                    # print "\tI found", order, "at index", found_index, "at freq", last_sb

                sb_freq_guess.append(x_axis[found_index])
                sb_amp_guess.append(check_max_area - (2 * octant + 1) * check_ave)
                error_est = np.sqrt(sum([i ** 2 for i in error[
                                                         found_index - octant:found_index + octant]])) / (
                                    check_max_area - (2 * octant + 1) * check_ave)
                # This error is a relative error.
                if verbose:
                    print(". Err = {:.3g}".format(error_est))
                    # print "\tMy error estimate is:", error_est
                # print "My relative error is:", error_est / sb_amp_guess
                sb_error_est.append(error_est)
                self.sb_list.append(order)
                consecutive_null_sb = 0
                if order % 2 == 1:
                    consecutive_null_odd = 0
            else:
                # print "I could not find sideband with order", order
                last_sb = last_sb + thz_freq
                consecutive_null_sb += 1
                if order % 2 == 1:
                    consecutive_null_odd += 1
                if verbose:
                    print("\t\tI did not count this sideband")
            if consecutive_null_odd == 1 and no_more_odds == False:
                # print "I'm done looking for odd sidebands"
                no_more_odds = True
            if consecutive_null_sb == 2:
                # print "I can't find any more sidebands"
                break

        if verbose:
            print("I found these sidebands:", self.sb_list)
            print('-' * 15)
            print()
            print()
        self.sb_guess = np.array([np.asarray(sb_freq_guess), np.asarray(sb_amp_guess),
                                  np.asarray(sb_error_est)]).T
        # self.sb_guess = [frequency guess, amplitude guess, relative error of amplitude] for each sideband.

    def guess_sidebandsOld(self, cutoff=4.5, verbose=False, plot=False, **kwargs):
        """
        05/24/18
        Old code from Hunter's days (or nearly, I've already started cleaning some
        stuff up). keeping it around in case I break too much stuff

        Finds the locations of all the sidebands in the proc_data array to be
        able to seed the fitting method.  This works by finding the maximum data
        value in the array and guessing what sideband it is.  It creates an array
        that includes this information.  It will then step down, initially by one
        THz frequency, then by twos after it hasn't found any odd ones.  It then
        goes up from the max and finds everything above in much the same way.

        There is currently no rhyme or reason to a cutoff of 8.  I don't know what
        it should be changed to, though.

        Input:
        cutoff = signal-to-noise threshold to count a sideband candidate.

        kwargs:
           window_size: how big of a window (in pixels) to use for checking for
                sidebands. Specified in half-width
              default: 15


        Internal:
        self.sb_list = List of all of the orders the method found
        self.sb_index = index of all of the peaks of the sidebands
        self.sb_guess = three-part list including the frequency, amplitude and
                        error guesses for each sideband
        """
        # TODO: this isn't commented appropriately.  Will it be made more readable first?

        if "cutoff" in self.parameters:
            cutoff = self.parameters["cutoff"]
        else:
            self.parameters['cutoff for guess_sidebands'] = cutoff

        if verbose:
            print("=" * 15)
            print()
            print("Guessing CCD Sideband parameters")
            print(os.path.basename(self.fname))
            print("\tCutoff = {}".format(cutoff))
            print()
            print("=" * 15)
        x_axis = np.array(self.proc_data[:, 0])
        y_axis = np.array(self.proc_data[:, 1])
        error = np.array(self.proc_data[:, 2])

        min_sb = int(self.calc_approx_sb_order(x_axis[0])) + 1
        try:
            max_sb = int(self.calc_approx_sb_order(x_axis[-1]))
        except ValueError:
            print(x_axis)

        nir_freq = self.parameters["nir_freq"]
        thz_freq = self.parameters["thz_freq"]

        if verbose:
            print("min_sb: {} | max_sb: {}".format(min_sb, max_sb))

        # Find max strength sideband and it's order
        global_max = np.argmax(y_axis)
        order_init = int(round(self.calc_approx_sb_order(x_axis[global_max])))
        # if verbose:
        #     print "The global max is at index", global_max
        if global_max < 15:
            check_y = y_axis[:global_max + 15]
            check_y = np.concatenate((np.zeros(15 - global_max), check_y))
        elif global_max > 1585:
            check_y = y_axis[global_max - 15:]
            check_y = np.concatenate((check_y, np.zeros(global_max - 1585)))
        else:
            check_y = y_axis[global_max - 15:global_max + 15]

        check_max_index = np.argmax(check_y)
        check_max_area = np.sum(check_y[check_max_index - 2:check_max_index + 3])

        check_ave = np.mean(check_y[[0, 1, 2, 3, 4, -1, -2, -3, -4, -5]])
        check_stdev = np.std(check_y[[0, 1, 2, 3, 4, -1, -2, -3, -4, -5]])
        check_ratio = (check_max_area - 3 * check_ave) / check_stdev

        if verbose:
            print(("{:^16}" * 5).format(
                "global_max idx", "check_max_area", "check_ave", "check_stdev",
                "check_ratio"))
            print(("{:^16.5g}" * 5).format(
                global_max, check_max_area, check_ave, check_stdev, check_ratio))

        if check_ratio > cutoff:
            self.sb_list = [order_init]
            self.sb_index = [global_max]
            sb_freq_guess = [x_axis[global_max]]
            sb_amp_guess = [y_axis[global_max]]
            sb_error_est = [
                np.sqrt(sum([i ** 2 for i in error[global_max - 2:global_max + 3]])) / (
                            check_max_area - 5 * check_ave)]
        else:
            print("There are no sidebands in", self.fname)
            raise RuntimeError

        if verbose:
            print("\t Looking for sidebands with f < {:.6f}".format(sb_freq_guess[0]))
        last_sb = sb_freq_guess[0]
        index_guess = global_max
        # keep track of how many consecutive sidebands we've skipped. Sometimes one's
        #  noisy or something, so we want to keep looking after skipping one
        consecutive_null_sb = 0
        consecutive_null_odd = 0
        no_more_odds = False
        break_condition = False
        for order in range(order_init - 1, min_sb - 1, -1):
            # Check to make sure we're not looking at an odd when
            # we've decided to skip them.
            if no_more_odds == True and order % 2 == 1:
                last_sb = last_sb - thz_freq
                if verbose:
                    print("I skipped", order)
                continue

            # Window size to look for next sideband. Needs to be order dependent
            # because higher orders get wider, so we need to look at more.
            # Values are arbitrary.
            window_size = 0.45 + 0.0004 * order  # used to be last_sb?
            lo_freq_bound = last_sb - thz_freq * (
                        1 + window_size)  # Not sure what to do about these
            hi_freq_bound = last_sb - thz_freq * (1 - window_size)

            if verbose:
                print("\nSideband", order)
                print("\t{:.4f} < f_{} < {:.4f}".format(lo_freq_bound, order,
                                                        hi_freq_bound))

            # Get the indices where the energies lie within the bounds for this SB
            sliced_indices = \
            np.where((x_axis > lo_freq_bound) & (x_axis < hi_freq_bound))[0]
            start_index, end_index = sliced_indices.min(), sliced_indices.max()

            # Get a slice of the y_data which is only in the region of interest
            check_y = y_axis[sliced_indices]

            check_max_index = np.argmax(
                check_y)  # This assumes that two floats won't be identical
            # Calculate the "area" of the sideband by looking at the peak value
            # within the range, and the pixel above/below it
            check_max_area = np.sum(check_y[check_max_index - 1:check_max_index + 2])

            if verbose and plot:
                plt.figure("CCD data")
                plt.plot([lo_freq_bound] * 2, [0, check_y[check_max_index]], 'b')
                plt.plot([hi_freq_bound] * 2, [0, check_y[check_max_index]], 'b')
                plt.plot([lo_freq_bound, hi_freq_bound], [check_y[check_max_index]] *
                         2, 'b', label="{} Box".format(order))
                plt.text((lo_freq_bound + hi_freq_bound) / 2, check_y[check_max_index],
                         order)

            # get the slice that doesn't have the peak in it to compare statistics
            check_region = np.append(check_y[:check_max_index - 1],
                                     check_y[check_max_index + 2:])
            check_ave = check_region.mean()
            check_stdev = check_region.std()

            # Calculate an effective SNR, where check_ave is roughly the
            # background level
            check_ratio = (check_max_area - 3 * check_ave) / check_stdev

            if order % 2 == 1:  # This raises the barrier for odd sideband detection
                check_ratio = check_ratio / 1.5
            if verbose:
                print("\t" + ("{:^14}" * 4).format(
                    "check_max_area", "check_ave", "check_stdev", "check_ratio"))
                print("\t" + ("{:^14.5g}" * 4).format(
                    check_max_area, check_ave, check_stdev, check_ratio))

            if check_ratio > cutoff:
                found_index = check_max_index + start_index
                self.sb_index.append(found_index)
                last_sb = x_axis[found_index]

                if verbose:
                    print("I just found", last_sb)

                sb_freq_guess.append(x_axis[found_index])
                sb_amp_guess.append(check_max_area - 3 * check_ave)
                error_est = np.sqrt(
                    sum(
                        [i ** 2 for i in error[found_index - 1:found_index + 2]]
                    )) / (check_max_area - 3 * check_ave)
                if verbose:
                    print("My error estimate is:", error_est)
                sb_error_est.append(error_est)
                self.sb_list.append(order)
                consecutive_null_sb = 0
                if order % 2 == 1:
                    consecutive_null_odd = 0
            else:
                # print "I could not find sideband with order", order
                last_sb = last_sb - thz_freq
                consecutive_null_sb += 1
                if order % 2 == 1:
                    consecutive_null_odd += 1
            if consecutive_null_odd == 1 and no_more_odds == False:
                # print "I'm done looking for odd sidebands"
                no_more_odds = True
            if consecutive_null_sb == 2:
                # print "I can't find any more sidebands"
                break

        # Look for higher sidebands
        if verbose: print("\nLooking for higher energy sidebands")

        last_sb = sb_freq_guess[0]
        index_guess = global_max
        consecutive_null_sb = 0
        consecutive_null_odd = 0
        no_more_odds = False
        break_condition = False
        for order in range(order_init + 1, max_sb + 1):
            if no_more_odds == True and order % 2 == 1:
                last_sb = last_sb + thz_freq
                continue
            window_size = 0.45 + 0.001 * order  # used to be 0.28 and 0.0004
            lo_freq_bound = last_sb + thz_freq * (
                        1 - window_size)  # Not sure what to do about these
            hi_freq_bound = last_sb + thz_freq * (1 + window_size)

            start_index = False
            end_index = False

            if verbose:
                print("\nSideband", order)
                # print "The low frequency bound is", lo_freq_bound
                # print "The high frequency bound is", hi_freq_bound
                print("\t{:.4f} < f_{} < {:.4f}".format(lo_freq_bound, order,
                                                        hi_freq_bound))
            for i in range(index_guess, 1600):
                if start_index == False and i == 1599:
                    # print "I'm all out of space, captain!"
                    break_condition = True
                    break
                elif start_index == False and x_axis[i] > lo_freq_bound:
                    # print "start_index is", i
                    start_index = i
                elif i == 1599:
                    end_index = 1599
                    # print "hit end of data, end_index is 1599"
                elif end_index == False and x_axis[i] > hi_freq_bound:
                    end_index = i
                    # print "end_index is", i
                    index_guess = i
                    break
            if break_condition:
                break
            check_y = y_axis[start_index:end_index]

            check_max_index = np.argmax(
                check_y)  # This assumes that two floats won't be identical
            octant = len(check_y) // 8  # To be able to break down check_y into eighths
            if octant < 1:
                octant = 1

            check_max_area = np.sum(
                check_y[check_max_index - octant - 1:check_max_index + octant + 1])

            if verbose and plot:
                plt.figure("CCD data")
                plt.plot([lo_freq_bound] * 2, [0, check_y[check_max_index]], 'b')
                plt.plot([hi_freq_bound] * 2, [0, check_y[check_max_index]], 'b')
                plt.plot([lo_freq_bound, hi_freq_bound], [check_y[check_max_index]] *
                         2, 'b', label=order)
                plt.text((lo_freq_bound + hi_freq_bound) / 2, check_y[check_max_index],
                         order)

            no_peak = (2 * len(
                check_y)) // 6  # The denominator is in flux, used to be 5
            # if verbose: print "\tcheck_y length", len(check_y)

            check_ave = np.mean(np.take(check_y, np.concatenate(
                (np.arange(no_peak), np.arange(-no_peak, 0)))))
            check_stdev = np.std(np.take(check_y, np.concatenate(
                (np.arange(no_peak), np.arange(-no_peak, 0)))))

            check_ratio = (check_max_area - (2 * octant + 1) * check_ave) / check_stdev

            if verbose:
                print("\tIndices: {}->{} (d={})".format(start_index, end_index,
                                                        len(check_y)))
                # print "check_y is", check_y
                # print "\ncheck_max_area is", check_max_area
                # print "check_ave is", check_ave
                # print "check_stdev is", check_stdev
                # print "check_ratio is", check_ratio

                print("\t" + ("{:^14}" * 4).format(
                    "check_max_area", "check_ave", "check_stdev", "check_ratio"))
                print("\t" + ("{:^14.6g}" * 4).format(
                    check_max_area, check_ave, check_stdev, check_ratio))

            if order % 2 == 1:  # This raises the barrier for odd sideband detection
                check_ratio = check_ratio / 2
            if check_ratio > cutoff:
                found_index = check_max_index + start_index
                self.sb_index.append(found_index)
                last_sb = x_axis[found_index]

                if verbose:
                    print("\tI'm counting this SB at index {} (f={:.4f})".format(
                        found_index, last_sb), end=' ')
                    # print "\tI found", order, "at index", found_index, "at freq", last_sb

                sb_freq_guess.append(x_axis[found_index])
                sb_amp_guess.append(check_max_area - (2 * octant + 1) * check_ave)
                error_est = np.sqrt(sum([i ** 2 for i in error[
                                                         found_index - octant:found_index + octant]])) / (
                                    check_max_area - (2 * octant + 1) * check_ave)
                # This error is a relative error.
                if verbose:
                    print(". Err = {:.3g}".format(error_est))
                    # print "\tMy error estimate is:", error_est
                # print "My relative error is:", error_est / sb_amp_guess
                sb_error_est.append(error_est)
                self.sb_list.append(order)
                consecutive_null_sb = 0
                if order % 2 == 1:
                    consecutive_null_odd = 0
            else:
                # print "I could not find sideband with order", order
                last_sb = last_sb + thz_freq
                consecutive_null_sb += 1
                if order % 2 == 1:
                    consecutive_null_odd += 1
                if verbose:
                    print("\t\tI did not count this sideband")
            if consecutive_null_odd == 1 and no_more_odds == False:
                # print "I'm done looking for odd sidebands"
                no_more_odds = True
            if consecutive_null_sb == 2:
                # print "I can't find any more sidebands"
                break

        if verbose:
            print("I found these sidebands:", self.sb_list)
            print('-' * 15)
            print()
            print()
        self.sb_guess = np.array([np.asarray(sb_freq_guess), np.asarray(sb_amp_guess),
                                  np.asarray(sb_error_est)]).T
        # self.sb_guess = [frequency guess, amplitude guess, relative error of amplitude] for each sideband.

    def fit_sidebands(self, plot=False, verbose=False):
        """
        This takes self.sb_guess and fits to each maxima to get the details of
        each sideband.  It's really ugly, but it works.  The error of the
        sideband area is approximated from the data, not the curve fit.  All
        else is from the curve fit.  Which is definitely underestimating the
        error, but we don't care too much about those errors (at this point).

        self.sb_guess = [frequency guess, amplitude guess, relative error of amplitude] for each sideband.

        Temporary stuff:
        sb_fits = holder of the fitting results until all spectra have been fit
        window = an integer that determines the "radius" of the fit window, proportional to thz_freq.

        Attributes created:
        self.sb_results = the money maker.  Column order:
                          [sb number, Freq (eV), Freq error (eV), Gauss area (arb.), Area error, Gauss linewidth (eV), Linewidth error (eV)]
                          [    0    ,      1   ,        2,      ,        3         ,      4    ,         5           ,        6            ]
        self.full_dict = a dictionary similar to sb_results, but now the keys
                         are the sideband orders.  Column ordering is otherwise the same.
        :param plot: Do you want to see the fits plotted with the data?
        :type plot: bool
        :param verbose: Do you want to see the details AND the initial guess fits?
        :type verbose: bool
        :return: None
        """
        # print "Trying to fit these"
        sb_fits = []

        if verbose:
            print("=" * 15)
            print()
            print("Fitting CCD Sidebands")
            print(os.path.basename(self.fname))
            print()
            print("=" * 15)
        # pretty sure you want this up here so things don't break
        # when no sidebands found
        self.full_dict = {}
        thz_freq = self.parameters["thz_freq"]
        window = 15 + int(15 * thz_freq / 0.0022) # Adjust the fit window based on the sideband spacing
                                                  # The 15's are based on empirical knowledge that for
                                                  # 540 GHz (2.23 meV), the best window size is 30 and
                                                  # that it seems like the window size should grow slowly?
        for elem, peakIdx in enumerate(self.sb_index):  # Have to do this because guess_sidebands
            # doesn't out put data in the most optimized way
            if peakIdx < window:
                data_temp = self.proc_data[:peakIdx + window, :]
            elif (1600 - peakIdx) < window:
                data_temp = self.proc_data[peakIdx - window:, :]
            else:
                data_temp = self.proc_data[peakIdx - window:peakIdx + window, :]
            width_guess = 0.0001 + 0.000001 * self.sb_list[elem]  # so the width guess gets wider as order goes up
            p0 = np.array([self.sb_guess[elem, 0],
                           self.sb_guess[elem, 1] * width_guess,
                           width_guess,
                           0.1])
            # print "Let's fit this shit!"
            if verbose:
                print("Fitting SB {}. Peak index: {}, {}th peak in spectra".format(
                    self.sb_list[elem], peakIdx, elem
                ))
                # print "\nnumber:", elem, num
                # print "data_temp:", data_temp
                # print "p0:", p0
                print(' '*20 +"p0 = " + np.array_str(p0, precision=4))
            # plot_guess = True  # This is to disable plotting the guess function
            if verbose and plot:
                plt.figure('CCD data')
                linewidth = 3
                x_vals = np.linspace(data_temp[0, 0], data_temp[-1, 0], num=500)
                if elem != 0:
                    try:
                        plt.plot(x_vals, gauss(x_vals, *p0),
                                 plt.gca().get_lines()[-1].get_color() + '--'  # I don't really know. Mostly
                                 # just looked around at what functions
                                 # matplotlib has...
                                 , linewidth=linewidth)
                    except:  # to prevent weird mac issues with the matplotlib things?
                        plt.plot(x_vals, gauss(x_vals, *p0), '--', linewidth=linewidth)

                else:
                    plt.plot(x_vals, gauss(x_vals, *p0), '--', linewidth=linewidth)

            try:
                # 11/1/16
                # needed to bump maxfev up to 2k because a sideband wasn't being fit
                # Fix for sb 106
                # 05-23 Loren 10nm\hsg_640_Perp352seq_spectrum.txt
                coeff, var_list = curve_fit(
                    gauss, data_temp[:, 0], data_temp[:, 1], p0=p0, maxfev = 2000)
            except Exception as e:
                if verbose:
                    print("\tThe fit failed:")
                    print("\t\t", e)
                    print("\tFitting region: {}->{}".format(peakIdx-window, peakIdx+window))
                    # print "I couldn't fit", elem
                    # print "It's sideband", num
                    # print "In file", self.fname
                    # print "because", e
                    # print "wanted to fit xindx", peakIdx, "+-", window
                self.sb_list[elem] = None
                continue # This will ensure the rest of the loop is not run without an actual fit.

            coeff[1] = abs(coeff[1])  # The amplitude could be negative if the linewidth is negative
            coeff[2] = abs(coeff[2])  # The linewidth shouldn't be negative
            if verbose:
                print("\tFit successful: ", end=' ')
                print("p = " + np.array_str(coeff, precision=4))
                # print "coeffs:", coeff
                # print "sigma for {}: {}".format(self.sb_list[elem], coeff[2])
            if 10e-4 > coeff[2] > 10e-6:
                try:
                    sb_fits.append(np.hstack((self.sb_list[elem], coeff, np.sqrt(np.diag(var_list)))))
                except RuntimeWarning:
                    sb_fits.append(np.hstack((self.sb_list[elem], coeff, np.sqrt(np.abs(np.diag(var_list))))))

                # the var_list wasn't approximating the error well enough, even when using sigma and absoluteSigma
                # self.sb_guess[elem, 2] is the relative error as calculated by the guess_sidebands method
                # coeff[1] is the area from the fit.  Therefore, the product should be the absolute error
                # of the integrated area of the sideband.  The other errors are still underestimated.
                #
                # 1/12/18 note: So it looks like what hunter did is calculate an error estimate
                # for the strength/area by the quadrature sum of errors of the points in the peak
                # (from like 813 in guess_sidebands:
                #    error_est = np.sqrt(sum([i ** 2 for i in error[found_index - 1:found_index + 2]])) / (
                # Where the error is what comes from the CCD by averaging 4 spectra. As far as I can tell,
                # it doesn't currently pull in the dark counts or anything like that, except maybe
                # indirectly since it'll cause the variations in the peaks
                sb_fits[-1][6] = self.sb_guess[elem, 2] * coeff[1]
                if verbose:
                    print("\tRel.Err: {:.4e}  |  Abs.Err: {:.4e}".format(
                        self.sb_guess[elem, 2], coeff[1] * self.sb_guess[elem, 2]
                    ))
                    print()
                    # print "The rel. error guess is", self.sb_guess[elem, 2]
                    # print "The abs. error guess is", coeff[1] * self.sb_guess[elem, 2]

                # The error from self.sb_guess[elem, 2] is a relative error
            if plot and verbose:
                plt.figure('CCD data')
                linewidth = 5
                x_vals = np.linspace(data_temp[0, 0], data_temp[-1, 0], num=500)
                if elem != 0:
                    try:
                        plt.plot(x_vals, gauss(x_vals, *coeff),
                                 plt.gca().get_lines()[-1].get_color() + '--'  # I don't really know. Mostly
                                 # just looked around at what functions
                                 # matplotlib has...
                                 , linewidth=linewidth)
                    except:  # to prevent weird mac issues with the matplotlib things?
                        plt.plot(x_vals, gauss(x_vals, *coeff), '--', linewidth=linewidth)

                else:
                    plt.plot(x_vals, gauss(x_vals, *coeff), '--', linewidth=linewidth)
        sb_fits_temp = np.asarray(sb_fits)
        reorder = [0, 1, 5, 2, 6, 3, 7, 4, 8]
        # Reorder the list to put the error of the i-th parameter as the i+1th.
        try:
            sb_fits = sb_fits_temp[:, reorder]
            # if verbose: print "The abs. error guess is", sb_fits[:, 0:5]
        except:
            raise RuntimeError("No sidebands to fit?")

        # Going to label the appropriate row with the sideband
        self.sb_list = sorted(list([x for x in self.sb_list if x is not None]))
        sb_names = np.vstack(self.sb_list)

        # Sort by SB order
        sorter = np.argsort(sb_fits[:, 0])
        self.sb_results = np.array(sb_fits[sorter, :7])

        if verbose:
            print("\tsb_results:")
            print("\t\t" + ("{:^5s}" + ("{:^12s}")*(self.sb_results.shape[1]-1)).format(
                "SB", "Cen.En.", "", "Area", "", "Width",""))
            for line in self.sb_results:
                print('\t\t[' + ("{:^5.0f}"+ "{:<12.4g}"*(line.size-1)).format(*line) + ']')
            print('-'*19)
        self.full_dict = {}
        for sb in self.sb_results:
            self.full_dict[sb[0]] = np.asarray(sb[1:])

    def infer_frequencies(self, nir_units="wavenumber", thz_units="GHz", bad_points=-2):
        """
        This guy tries to fit the results from fit_sidebands to a line to get the relevant frequencies
        :param nir_units: What units do you want this to output?
        :type nir_units: 'nm', 'wavenumber', 'eV', 'THz'
        :param thz_units: What units do you want this to output for the THz?
        :type thz_units: 'GHz', 'wavenumber', 'meV'
        :param bad_points: How many more-positive order sidebands shall this ignore?
        :type bad_points: int
        :return: freqNIR, freqTHz, the frequencies in the appropriate units
        """
        # force same units for in dict
        freqNIR, freqTHz = calc_laser_frequencies(self, "wavenumber", "wavenumber", bad_points)

        self.parameters["calculated NIR freq (cm-1)"] = "{}".format(freqNIR, nir_units)
        self.parameters["calculated THz freq (cm-1)"] = "{}".format(freqTHz, freqTHz)
        freqNIR, freqTHz = calc_laser_frequencies(self, nir_units, thz_units, bad_points)
        return freqNIR, freqTHz

    def save_processing(self, file_name, folder_str, marker='', index='', verbose=''):
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
        Two files:
            self.proc_data = the continuous spectrum
            self.sb_results = the individual sideband details

        :param file_name: The base name for the saved file
        :type file_name: str
        :param folder_str: The full name for the folder hte file is saved it.  Folder can be created
        :type folder_str: str
        :param marker: Marker for the file, appended to file_name, often the self.parameters['series']
        :type marker: str
        :param index: used to keep these files from overwriting themselves when marker is the same
        :type index: str or int
        :return: None
        """
        try:
            os.mkdir(folder_str)
        except OSError as e:
            if e.errno == errno.EEXIST:
                pass
            else:
                raise
        temp = np.array(self.sb_results)

        ampli = np.array([temp[:, 3] / temp[:, 5]])  # But [:, 3] is already area?
        # (The old name was area)
        # I think it must be amplitude
        temp[:, 5:7] = temp[:, 5:7] * 1000  # For meV linewidths
        if verbose:
            print("sb_results", self.sb_results.shape)
            print("ampli", ampli.shape)
        save_results = np.hstack((temp, ampli.T))

        spectra_fname = file_name + '_' + marker + '_' + str(index) + '.txt'
        fit_fname = file_name + '_' + marker + '_' + str(index) + '_fits.txt'
        self.save_name = spectra_fname

        self.parameters['addenda'] = self.addenda
        self.parameters['subtrahenda'] = self.subtrahenda
        try:
            parameter_str = json.dumps(self.parameters, sort_keys=True, indent=4, separators=(',', ': '))
        except:
            print("Source: EMCCD_image.save_images\nJSON FAILED")
            print("Here is the dictionary that broke JSON:\n", self.parameters)
            return
        parameter_str = parameter_str.replace('\n', '\n#')

        num_lines = parameter_str.count('#')  # Make the number of lines constant so importing is easier
        # for num in range(99 - num_lines): parameter_str += '\n#'
        parameter_str += '\n#' * (99 - num_lines)
        origin_import_spec = '\nNIR frequency,Signal,Standard error\neV,arb. u.,arb. u.'
        spec_header = '#' + parameter_str + origin_import_spec

        origin_import_fits  = '\nSideband,Center energy,error,Sideband strength,error,Linewidth,error,Amplitude'
        origin_import_fits += '\norder,eV,,arb. u.,,meV,,arb. u.'
        origin_import_fits += "\n{},,,{},,,".format(marker, marker)
        fits_header = '#' + parameter_str + origin_import_fits

        # print "DEBUG: in saving", folder_str, ",", spectra_fname

        np.savetxt(os.path.join(folder_str, spectra_fname), self.proc_data, delimiter=',',
                   header=spec_header, comments='', fmt='%0.6e')
        np.savetxt(os.path.join(folder_str, fit_fname), save_results, delimiter=',',
                   header=fits_header, comments='', fmt='%0.6e')
        if verbose:
            print("Save image.\nDirectory: {}".format(os.path.join(folder_str, spectra_fname)))

class HighSidebandCCDRaw(HighSidebandCCD):
    """
    This class is meant for passing in an image file (currently supports a 2x1600)
    Which it does all the processing on.
    """
    def __init__(self, hsg_thing, parameter_dict=None, spectrometer_offset=None):
        # let the supers do the hard work of importing the json dict and all that jazz
        super(HighSidebandCCDRaw, self).__init__(hsg_thing, parameter_dict=None, spectrometer_offset=None)
        self.ccd_data = np.genfromtxt(hsg_thing, delimiter=',').T
        self.proc_data = np.column_stack((
            self.gen_wavelengths(self.parameters["center_lambda"], self.parameters["grating"]),
            np.array(self.ccd_data[:,1], dtype=float)-np.median(self.ccd_data[:,1]),
            np.ones_like(self.ccd_data[:,1], dtype=float)
                                         ))


        self.proc_data[:, 0] = 1239.84 / self.proc_data[:, 0]
        self.proc_data = np.flipud(self.proc_data)

    @staticmethod
    def gen_wavelengths(center_lambda, grating):
        '''
        This returns a 1600 element list of wavelengths for each pixel in the EMCCD based on grating and center wavelength

        grating = which grating, 1 or 2
        center = center wavelength in nanometers
        '''
        b = 0.75  # length of spectrometer, in m
        k = -1.0  # order looking at
        r = 16.0e-6  # distance between pixles on CCD

        if grating == 1:
            d = 1. / 1800000.
            gamma = 0.213258508834
            delta = 1.46389935365
        elif grating == 2:
            d = 1. / 1200000.
            gamma = 0.207412628027
            delta = 1.44998344749
        elif grating == 3:
            d = 1. / 600000.
            gamma = 0.213428934011
            delta = 1.34584754696
        else:
            print("What a dick, that's not a valid grating")
            return None

        center = center_lambda * 10 ** -9
        wavelength_list = np.arange(-799.0, 801.0)

        output = d * k ** (-1) * ((-1) * np.cos(delta + gamma + (-1) * np.arccos(
            (-1 / 4) * (1 / np.cos((1 / 2) * gamma)) ** 2 * (
            2 * (np.cos((1 / 2) * gamma) ** 4 * (2 + (-1) * d ** (-2) * k ** 2 * center ** 2 + 2 * np.cos(gamma))) ** (
            1 / 2) + d ** (-1) * k * center * np.sin(gamma))) + np.arctan(
            b ** (-1) * (r * wavelength_list + b * np.cos(delta + gamma)) * (1 / np.sin(delta + gamma)))) + (
                                  1 + (-1 / 16) * (1 / np.cos((1 / 2) * gamma)) ** 4 * (2 * (
                                  np.cos((1 / 2) * gamma) ** 4 * (
                                  2 + (-1) * d ** (-2) * k ** 2 * center ** 2 + 2 * np.cos(gamma))) ** (1 / 2) + d ** (
                                                                                        -1) * k * center * np.sin(
                                      gamma)) ** 2) ** (1 / 2))

        output = (output + center) * 10 ** 9
        return output

class PMT(object):
    def __init__(self, file_name):
        """
        Initializes a SPEX spectrum.  It'll open a file, and bring in the details
        of a sideband spectrum into the object.  There isn't currently any reason
        to use inheritance here, but it could be extended later to include PLE or
        something of the sort.

        attributes:
            self.parameters - dictionary of important experimental parameters
                              this will not necessarily be the same for each
                              file in the object
            self.fname - the current file path

        :param file_name: The name of the PMT file
        :type file_name: str
        :return: None
        """
        # print "This started"
        self.fname = file_name
        # self.files_included = [file_name]
        with open(file_name, 'r') as f:
            param_str = ''
            line = f.readline()  # Needed to move past the first line, which is the sideband order.  Not generally useful
            line = f.readline()
            while line[0] == '#':
                param_str += line[1:]
                line = f.readline()

            self.parameters = json.loads(param_str)

class HighSidebandPMT(PMT):
    def __init__(self, file_path, verbose=False):
        """
        Initializes a SPEX spectrum.  It'll open a single file, then read
        the data from that file using .add_sideband().  The super's init will handle the parameters
        and the description.

        attributes:
            self.parameters - dictionary of important experimental parameters, created in PMT
            self.sb_dict - keys are sideband order, values are PMT data arrays
            self.sb_list - sorted list of included sidebands

        :param file_path: path to the current file
        :type file_path: str
        :param verbose: Flag to see the nitty gritty details
        :type verbose: bool
        :return:
        """
        super(HighSidebandPMT, self).__init__(
            file_path)  # Creates the json parameters dictionary
        self.fname = file_path
        self.parameters["files included"] = [file_path]
        with open(file_path, 'r') as f:
            sb_num = int(f.readline()[1:])
        raw_temp = np.genfromtxt(file_path, comments='#', delimiter=',')[3:, :]

        if self.parameters.get("photon counted", False):
            # The scale factor for photon counting to generic
            # PMT data depends on... things. It's different each
            # day. Unfortunately, the overlap in dynamic range between
            # the two is small, and generally only one sideband
            # can been seen by both methods. I don't really have
            # the motivation to automatically calculate the
            # appropriate factor, so this is your reminder to find
            # it yourself.
            import time
            # assert time.strftime("%x") == "03/15/17"
            assert self.parameters.get("pc ratio", -1) != -1, self.fname
            raw_temp[:,3] *= self.parameters["pc ratio"]
            pass
        raw_temp[:, 0] = raw_temp[:, 0] / 8065.6  # turn NIR freq into eV
        self.parameters["thz_freq"] = 0.000123984 * float(
            self.parameters.get("fel_lambda", -1))
        self.parameters["nir_freq"] =  float(
            self.parameters.get("nir_lambda", -1))/8065.6


        self.initial_sb = sb_num
        self.initial_data = np.array(raw_temp)
        self.sb_dict = {sb_num: np.array(raw_temp)}
        self.sb_list = [sb_num]

    def add_sideband(self, other):
        """
        This bad boy will add another PMT sideband object to the sideband spectrum of this object.  It handles
        when you measure the same sideband twice.  It assumes both are equally "good"
        NOTE: This means that if both aren't equally "good" (taking a second scan with higher
        gain/photon counting because you didn't see it), you need to not add the file
        (remove/rename the file, etc.)
        I'd love to overhall the data collection/analysis so this can be more intelligent
        (Effectively offload a lot of the processing (especially not saving 10 arbitrary
        points to process later) onto the live software and add sideband strengths alone,
        like the CCD works. But this would be a bigger change that I can seem to find
        time for).

        It currently doesn't do any sort of job combining dictionaries or anything, but it definitely could, if
        you have two incomplete dictionaries

        :param other: the new sideband data to add to the larger spectrum.  Add means append, no additino is performed
        :type other: HighSidebandPMT
        :return:
        """
        """
        This bad boy will add another PMT sideband object to the sideband spectrum of this object

        It currently doesn't do any sort of job combining dictionaries or anything, but it definitely could
        """
        self.parameters["files included"].append(other.fname)

        if other.initial_sb in self.sb_list:
            self.sb_list.append(other.initial_sb)

        # Make things comma delimited?
        try:
            self.sb_dict[other.initial_sb] = np.row_stack(
                (self.sb_dict[other.initial_sb], other.initial_data)
            )
        except KeyError:
            self.sb_dict[other.initial_sb] = np.array(other.initial_data)
        except Exception as e:
            print("THIS IS THE OTHER ERROR", e)
            raise

    def process_sidebands(self, verbose=False, baselineCorr = False):
        """
        This bad boy will clean up the garbled mess that is the object before hand,
        including clearing out misfired shots and doing the averaging.

        Affects:
            self.sb_dict = Averages over sidebands

        Creates:
            self.sb_list = The sideband orders included in this object.

        :param verbose: Flag to see the nitty gritty details.
        :type verbose: bool
        :param baselineCorr: Whether to subtract the average across
        the two endpoints
        :return: None
        """

        for sb_num, sb in list(self.sb_dict.items()):
            if sb_num == 0:
                fire_condition = -np.inf  # This way the FEL doesn't need to be on during laser line measurement
            else:
                fire_condition = np.mean(sb[:, 2]) / 2  # Say FEL fired if the
                # cavity dump signal is
                # more than half the mean
                # of the cavity dump signal
            frequencies = sorted(list(set(sb[:, 0])))

            temp = None
            for freq in frequencies:
                data_temp = np.array([])
                for point in sb:
                    if point[0] == freq and point[2] > fire_condition:
                        data_temp = np.hstack((data_temp, point[3]))
                try:
                    temp = np.vstack(
                        (temp, np.array([freq, np.mean(data_temp),
                                         np.std(data_temp) / np.sqrt(len(data_temp))])))
                except:
                    temp = np.array([freq, np.mean(data_temp),
                                     np.std(data_temp) / np.sqrt(len(data_temp))])
            # temp[:, 0] = temp[:, 0] / 8065.6  # turn NIR freq into eV
            temp = temp[temp[:, 0].argsort()]
            if baselineCorr:
                x = temp[[0, -1], 0]
                y = temp[[0, -1], 1]
                p = np.polyfit(x, y, 1)
                temp[:, 1] -= np.polyval(p, temp[:,0])


            self.sb_dict[sb_num] = np.array(temp)
        self.sb_list = sorted(self.sb_dict.keys())
        if verbose:
            print("Sidebands included", self.sb_list)

    def integrate_sidebands(self, verbose=False, cutoff=1.5, **kwargs):
        """
        This method will integrate the sidebands to find their strengths, and then
        use a magic number to define the width, since they are currently so utterly
        undersampled for fitting.

        cutoff is the ratio of area/error which must be exceeded to count

        It is currently the preferred method for calculating sideband strengths.
        self.fit_sidebands is probably better with better-sampled lines.

        Creates:
        self.sb_results = full list of integrated data. Column order is:
                          [sb order, Freq (eV), "error" (eV), Integrate area (arb.), area error, "Linewidth" (eV), "Linewidth error" (eV)
        self.full_dict = Dictionary where the SB order column is removed and turned into the keys.  The values
                         are the rest of that sideband's results.

        :param verbose: Flag to see the nitty gritty details
        :type verbose: bool
        :return: None
        """

        if verbose:
            print("="*15)
            print()
            print("Integrating PMT Sidebands")
            print("Cutoff: {}".format(cutoff))
            print(os.path.basename(self.fname))
            print()
            print("=" * 15)

        self.full_dict = {}
        for sideband in list(self.sb_dict.items()):
            index = np.argmax(sideband[1][:, 1])
            nir_frequency = sideband[1][index, 0]

            # stroff = np.nan_to_num(sideband[1][[0,1,-2,1], 1]).sum()/4.

            area = np.trapz(np.nan_to_num(sideband[1][:, 1]), sideband[1][:, 0])
            error = np.sqrt(np.sum(np.nan_to_num(
                sideband[1][:, 2]) ** 2)) / 8065.6  # Divide by the step size?
            if verbose:
                print("\torder: {}, area: {:.3g}, error: {:.3g}, ratio: {:.3f}".format(
                    sideband[0], area, error, area/error
                ))
            details = np.array(
                [sideband[0], nir_frequency, 1 / 8065.6, area, error, 2 / 8065.6,
                 1 / 8065.6])
            if area < 0:
                if verbose:
                    print("\t\tarea < 0")
                continue
            elif area < cutoff * error:  # Two seems like a good cutoff?
                if verbose:
                    print("\t\tI did not keep sideband")
                continue
            try:
                self.sb_results = np.vstack((self.sb_results, details))
            except:
                self.sb_results = np.array(details)
            self.full_dict[sideband[0]] = details[1:]
        try:
            self.sb_results = self.sb_results[self.sb_results[:, 0].argsort()]

        except (IndexError, AttributeError):
            # IndexError where there's only one sideband
            # AttributeError when there aren't any (one sb which wasn't fit)
            pass

        if verbose:
            print('-'*19)

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

        :param plot: Flag to see the results plotted
        :type plot: bool
        :param verbose: Flag to see the nitty gritty details
        :type verbose: bool
        :return: None
        """
        sb_fits = {}
        for sideband in list(self.sb_dict.items()):
            if verbose:
                print("Sideband number", sideband[0])
                print("Sideband data:\n", sideband[1])
            index = np.argmax(sideband[1][:, 1])
            nir_frequency = sideband[1][index, 0]
            peak = sideband[1][index, 1]
            width_guess = 0.0001  # Yep, another magic number
            p0 = [nir_frequency, peak * width_guess, width_guess, 0.00001]

            if verbose:
                x_vals = np.linspace(np.amin(sideband[1][:, 0]),
                                     np.amax(sideband[1][:, 0]), num=50)
                plt.plot(x_vals, gauss(x_vals, *p0),
                         label="fit :{}".format(sideband[1]))
                print("p0:", p0)
            try:
                coeff, var_list = curve_fit(gauss, sideband[1][:, 0], sideband[1][:, 1],
                                            sigma=sideband[1][:, 2], p0=p0)
                coeff[1] = abs(coeff[1])
                coeff[2] = abs(coeff[2])
                if verbose:
                    print("coeffs:", coeff)
                    print("stdevs:", np.sqrt(np.diag(var_list)))
                    print("integral", np.trapz(sideband[1][:, 1], sideband[1][:, 0]))
                if np.sqrt(np.diag(var_list))[0] / coeff[
                    0] < 0.5:  # The error on where the sideband is should be small
                    sb_fits[sideband[0]] = np.concatenate(
                        (np.array([sideband[0]]), coeff, np.sqrt(np.diag(var_list))))
                    # print "error then:", sb_fits[sideband[0]][6]
                    relative_error = np.sqrt(sum([x ** 2 for x in
                                                  sideband[1][index - 1:index + 2,
                                                  2]])) / np.sum(
                        sideband[1][index - 1:index + 2, 1])
                    if verbose:
                        print("relative error:", relative_error)
                    sb_fits[sideband[0]][6] = coeff[1] * relative_error
                    # print "error now:", sb_fits[sideband[0]][6]
                    if plot:
                        x_vals = np.linspace(np.amin(sideband[1][:, 0]),
                                             np.amax(sideband[1][:, 0]), num=50)
                        plt.plot(x_vals, gauss(x_vals, *coeff))
                        # plt.plot(x_vals, gauss(x_vals, *p0))
                else:
                    print("what happened?")
            except:
                print("God damn it, Leroy.\nYou couldn't fit this.")
                sb_fits[sideband[0]] = None

        for result in sorted(sb_fits.keys()):
            try:
                self.sb_results = np.vstack((self.sb_results, sb_fits[result]))
            except:
                self.sb_results = np.array(sb_fits[result])

        self.sb_results = self.sb_results[:, [0, 1, 5, 2, 6, 3, 7, 4, 8]]
        self.sb_results = self.sb_results[:, :7]
        if verbose:
            print("And the results, please:\n", self.sb_results)

        self.full_dict = {}
        for sb in self.sb_results:
            self.full_dict[sb[0]] = np.asarray(sb[1:])

    def laser_line(self, verbose=False, **kwargs):
        """
        This method is designed to scale everything in the PMT to the conversion
        efficiency based on our measurement of the laser line with a fixed
        attenuation.

        Creates:
            self.parameters['normalized?'] = Flag to specify if the laser has been
            accounted for.

        :return: None
        """

        if 0 not in self.sb_list:
            self.parameters['normalized?'] = False
            return
        else:
            laser_index = np.where(self.sb_results[:, 0] == 0)[0][0]
            if verbose:
                print("sb_results", self.sb_results[laser_index, :])
                print("laser_index", laser_index)

            laser_strength = np.array(self.sb_results[laser_index, 3:5])

            if verbose:
                print("Laser_strength", laser_strength)

            for sb in self.sb_results:
                sb[4] = (sb[3] / laser_strength[0]) * np.sqrt(
                    (sb[4] / sb[3]) ** 2 + (laser_strength[1] / laser_strength[0]) ** 2)
                sb[3] = sb[3] / laser_strength[0]
            for sb in list(self.full_dict.values()):
                sb[3] = (sb[2] / laser_strength[0]) * np.sqrt(
                    (sb[3] / sb[2]) ** 2 + (laser_strength[1] / laser_strength[0]) ** 2)
                sb[2] = sb[2] / laser_strength[0]
            self.parameters['normalized?'] = True

    def save_processing(self, file_name, folder_str, marker='', index='', verbose=False):
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
        Two files:
            self.proc_data = the continuous spectrum
            self.sb_results = the individual sideband details

        :param file_name: The base name for the saved file
        :type file_name: str
        :param folder_str: The full name for the folder hte file is saved it.  Folder can be created
        :type folder_str: str
        :param marker: Marker for the file, appended to file_name, often the self.parameters['series']
        :type marker: str
        :param index: used to keep these files from overwriting themselves when marker is the same
        :type index: str or int
        :return: None
        """
        try:
            os.mkdir(folder_str)
        except OSError as e:
            if e.errno == errno.EEXIST:
                pass
            else:
                raise

        spectra_fname = file_name + '_' + marker + '_' + str(index) + '.txt'
        fit_fname = file_name + '_' + marker + '_' + str(index) + '_fits.txt'
        self.save_name = spectra_fname
        # self.parameters["files included"] = list(self.files)
        try:
            parameter_str = json.dumps(self.parameters, sort_keys=True, indent=4,
                                       separators=(',', ': '))
        except:
            print("Source: PMT.save_images\nJSON FAILED")
            print("Here is the dictionary that broke JSON:\n", self.parameters)
            return
        parameter_str = parameter_str.replace('\n', '\n#')

        num_lines = parameter_str.count(
            '#')  # Make the number of lines constant so importing is easier
        # for num in range(99 - num_lines): parameter_str += '\n#'
        parameter_str += '\n#' * (99 - num_lines)

        origin_import_spec = '\nNIR frequency,Signal,Standard error\neV,arb. u.,arb. u.\n,{:.3f},'.format(
            self.parameters["fieldStrength"]["mean"])
        spec_header = '#' + parameter_str + origin_import_spec

        origin_import_fits = '\nCenter energy,error,Amplitude,error,Linewidth,error\neV,,arb. u.,,eV,,\n,,'  # + marker
        fits_header = '#' + parameter_str + origin_import_fits

        for sideband in sorted(self.sb_dict.keys()):
            try:
                complete = np.vstack((complete, self.sb_dict[sideband]))
            except:
                complete = np.array(self.sb_dict[sideband])

        np.savetxt(os.path.join(folder_str, spectra_fname), complete, delimiter=',',
                   header=spec_header, comments='', fmt='%0.6e')

        try:
            np.savetxt(os.path.join(folder_str, fit_fname), self.sb_results,
                       delimiter=',',
                       header=fits_header, comments='', fmt='%0.6e')
        except AttributeError:
            # Catch the error that happens if you save something without files
            print("warning, couldn't save fit file (no sidebands found?)")

        if verbose:
            print("Saved PMT spectrum.\nDirectory: {}".format(
                os.path.join(folder_str, spectra_fname)))

class HighSidebandPMTOld(PMT):
    """
    Old version: Replaced March 01, 2017

    Class initialized by loading in data set.

    Multiple copies of the same sideband were stacked as raw data and combined,
    effectively causing (2) 10-pt scans to be treated the same as (1) 20pt scan.
    This works well until you have photon counted pulses.
    """
    def __init__(self, file_path, verbose=False):
        """
        Initializes a SPEX spectrum.  It'll open a single file, then read
        the data from that file using .add_sideband().  The super's init will handle the parameters
        and the description.

        attributes:
            self.parameters - dictionary of important experimental parameters, created in PMT
            self.sb_dict - keys are sideband order, values are PMT data arrays
            self.sb_list - sorted list of included sidebands

        :param file_path: path to the current file
        :type file_path: str
        :param verbose: Flag to see the nitty gritty details
        :type verbose: bool
        :return:
        """
        super(HighSidebandPMT, self).__init__(
            file_path)  # Creates the json parameters dictionary
        self.fname = file_path
        self.parameters["files included"] = [file_path]
        with open(file_path, 'r') as f:
            sb_num = int(f.readline()[1:])
        raw_temp = np.genfromtxt(file_path, comments='#', delimiter=',')[3:, :]
        self.initial_sb = sb_num
        self.initial_data = np.array(raw_temp)
        self.sb_dict = {sb_num: np.array(raw_temp)}
        self.sb_list = [sb_num]

    def add_sideband(self, other):
        """
        This bad boy will add another PMT sideband object to the sideband spectrum of this object.  It handles
        when you measure the same sideband twice.  It assumes both are equally "good"

        It currently doesn't do any sort of job combining dictionaries or anything, but it definitely could, if
        you have two incomplete dictionaries

        :param other: the new sideband data to add to the larger spectrum.  Add means append, no additino is performed
        :type other: HighSidebandPMT
        :return:
        """
        """
        This bad boy will add another PMT sideband object to the sideband spectrum of this object

        It currently doesn't do any sort of job combining dictionaries or anything, but it definitely could
        """
        self.parameters["files included"].append(other.fname)

        if other.initial_sb in self.sb_list:
            self.sb_list.append(other.initial_sb)

        # Make things comma delimited?
        try:
            self.sb_dict[other.initial_sb].vstack((other.initial_data))
        except:
            self.sb_dict[other.initial_sb] = np.array(other.initial_data)

    def process_sidebands(self, verbose=False):
        """
        This bad boy will clean up the garbled mess that is the object before hand,
        including clearing out misfired shots and doing the averaging.

        Affects:
            self.sb_dict = Averages over sidebands

        Creates:
            self.sb_list = The sideband orders included in this object.

        :param verbose: Flag to see the nitty gritty details.
        :type verbose: bool
        :return: None
        """

        for sb_num, sb in list(self.sb_dict.items()):
            if sb_num == 0:
                fire_condition = -np.inf  # This way the FEL doesn't need to be on during laser line measurement
            else:
                fire_condition = np.mean(sb[:, 2]) / 2  # Say FEL fired if the
                # cavity dump signal is
                # more than half the mean
                # of the cavity dump signal
            frequencies = sorted(list(set(sb[:, 0])))

            temp = None
            for freq in frequencies:
                data_temp = np.array([])
                for point in sb:
                    if point[0] == freq and point[2] > fire_condition:
                        data_temp = np.hstack((data_temp, point[3]))
                try:
                    temp = np.vstack(
                        (temp, np.array([freq, np.mean(data_temp),
                                         np.std(data_temp) / np.sqrt(len(data_temp))])))
                except:
                    temp = np.array([freq, np.mean(data_temp),
                                     np.std(data_temp) / np.sqrt(len(data_temp))])
            temp[:, 0] = temp[:, 0] / 8065.6  # turn NIR freq into eV
            temp = temp[temp[:, 0].argsort()]
            self.sb_dict[sb_num] = np.array(temp)
        self.sb_list = sorted(self.sb_dict.keys())
        if verbose:
            print("Sidebands included", self.sb_list)

    def integrate_sidebands(self, verbose=False):
        """
        This method will integrate the sidebands to find their strengths, and then
        use a magic number to define the width, since they are currently so utterly
        undersampled for fitting.

        It is currently the preferred method for calculating sideband strengths.
        self.fit_sidebands is probably better with better-sampled lines.

        Creates:
        self.sb_results = full list of integrated data. Column order is:
                          [sb order, Freq (eV), "error" (eV), Integrate area (arb.), area error, "Linewidth" (eV), "Linewidth error" (eV)
        self.full_dict = Dictionary where the SB order column is removed and turned into the keys.  The values
                         are the rest of that sideband's results.

        :param verbose: Flag to see the nitty gritty details
        :type verbose: bool
        :return: None
        """
        self.full_dict = {}
        for sideband in list(self.sb_dict.items()):
            index = np.argmax(sideband[1][:, 1])
            nir_frequency = sideband[1][index, 0]
            area = np.trapz(np.nan_to_num(sideband[1][:, 1]), sideband[1][:, 0])
            error = np.sqrt(np.sum(np.nan_to_num(
                sideband[1][:, 2]) ** 2)) / 8065.6  # Divide by the step size?
            if verbose:
                print("order", sideband[0])
                print("area", area)
                print("error", error)
                print("ratio", area / error)
            details = np.array(
                [sideband[0], nir_frequency, 1 / 8065.6, area, error, 2 / 8065.6,
                 1 / 8065.6])
            if area < 0:
                if verbose:
                    print("area less than 0", sideband[0])
                continue
            elif area < 1.5 * error:  # Two seems like a good cutoff?
                if verbose:
                    print("I did not keep sideband ", sideband[0])
                continue
            try:
                self.sb_results = np.vstack((self.sb_results, details))
            except:
                self.sb_results = np.array(details)
            self.full_dict[sideband[0]] = details[1:]
        try:
            self.sb_results = self.sb_results[self.sb_results[:, 0].argsort()]

        except (IndexError, AttributeError):
            # IndexError where there's only one sideband
            # AttributeError when there aren't any (one sb which wasn't fit)
            pass

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

        :param plot: Flag to see the results plotted
        :type plot: bool
        :param verbose: Flag to see the nitty gritty details
        :type verbose: bool
        :return: None
        """
        sb_fits = {}
        for sideband in list(self.sb_dict.items()):
            if verbose:
                print("Sideband number", sideband[0])
                print("Sideband data:\n", sideband[1])
            index = np.argmax(sideband[1][:, 1])
            nir_frequency = sideband[1][index, 0]
            peak = sideband[1][index, 1]
            width_guess = 0.0001  # Yep, another magic number
            p0 = [nir_frequency, peak * width_guess, width_guess, 0.00001]

            if verbose:
                x_vals = np.linspace(np.amin(sideband[1][:, 0]),
                                     np.amax(sideband[1][:, 0]), num=50)
                plt.plot(x_vals, gauss(x_vals, *p0),
                         label="fit :{}".format(sideband[1]))
                print("p0:", p0)
            try:
                coeff, var_list = curve_fit(gauss, sideband[1][:, 0], sideband[1][:, 1],
                                            sigma=sideband[1][:, 2], p0=p0)
                coeff[1] = abs(coeff[1])
                coeff[2] = abs(coeff[2])
                if verbose:
                    print("coeffs:", coeff)
                    print("stdevs:", np.sqrt(np.diag(var_list)))
                    print("integral", np.trapz(sideband[1][:, 1], sideband[1][:, 0]))
                if np.sqrt(np.diag(var_list))[0] / coeff[
                    0] < 0.5:  # The error on where the sideband is should be small
                    sb_fits[sideband[0]] = np.concatenate(
                        (np.array([sideband[0]]), coeff, np.sqrt(np.diag(var_list))))
                    # print "error then:", sb_fits[sideband[0]][6]
                    relative_error = np.sqrt(sum([x ** 2 for x in
                                                  sideband[1][index - 1:index + 2,
                                                  2]])) / np.sum(
                        sideband[1][index - 1:index + 2, 1])
                    if verbose:
                        print("relative error:", relative_error)
                    sb_fits[sideband[0]][6] = coeff[1] * relative_error
                    # print "error now:", sb_fits[sideband[0]][6]
                    if plot:
                        x_vals = np.linspace(np.amin(sideband[1][:, 0]),
                                             np.amax(sideband[1][:, 0]), num=50)
                        plt.plot(x_vals, gauss(x_vals, *coeff))
                        # plt.plot(x_vals, gauss(x_vals, *p0))
                else:
                    print("what happened?")
            except:
                print("God damn it, Leroy.\nYou couldn't fit this.")
                sb_fits[sideband[0]] = None

        for result in sorted(sb_fits.keys()):
            try:
                self.sb_results = np.vstack((self.sb_results, sb_fits[result]))
            except:
                self.sb_results = np.array(sb_fits[result])

        self.sb_results = self.sb_results[:, [0, 1, 5, 2, 6, 3, 7, 4, 8]]
        self.sb_results = self.sb_results[:, :7]
        if verbose:
            print("And the results, please:\n", self.sb_results)

        self.full_dict = {}
        for sb in self.sb_results:
            self.full_dict[sb[0]] = np.asarray(sb[1:])

    def laser_line(self, verbose=False):
        """
        This method is designed to scale everything in the PMT to the conversion
        efficiency based on our measurement of the laser line with a fixed
        attenuation.

        Creates:
            self.parameters['normalized?'] = Flag to specify if the laser has been
            accounted for.

        :return: None
        """

        if 0 not in self.sb_list:
            self.parameters['normalized?'] = False
            return
        else:
            laser_index = np.where(self.sb_results[:, 0] == 0)[0][0]
            if verbose:
                print("sb_results", self.sb_results[laser_index, :])
                print("laser_index", laser_index)

            laser_strength = np.array(self.sb_results[laser_index, 3:5])

            if verbose:
                print("Laser_strength", laser_strength)

            for sb in self.sb_results:
                sb[4] = (sb[3] / laser_strength[0]) * np.sqrt(
                    (sb[4] / sb[3]) ** 2 + (laser_strength[1] / laser_strength[0]) ** 2)
                sb[3] = sb[3] / laser_strength[0]
            for sb in list(self.full_dict.values()):
                sb[3] = (sb[2] / laser_strength[0]) * np.sqrt(
                    (sb[3] / sb[2]) ** 2 + (laser_strength[1] / laser_strength[0]) ** 2)
                sb[2] = sb[2] / laser_strength[0]
            self.parameters['normalized?'] = True

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
        Two files:
            self.proc_data = the continuous spectrum
            self.sb_results = the individual sideband details

        :param file_name: The base name for the saved file
        :type file_name: str
        :param folder_str: The full name for the folder hte file is saved it.  Folder can be created
        :type folder_str: str
        :param marker: Marker for the file, appended to file_name, often the self.parameters['series']
        :type marker: str
        :param index: used to keep these files from overwriting themselves when marker is the same
        :type index: str or int
        :return: None
        """
        try:
            os.mkdir(folder_str)
        except OSError as e:
            if e.errno == errno.EEXIST:
                pass
            else:
                raise

        spectra_fname = file_name + '_' + marker + '_' + str(index) + '.txt'
        fit_fname = file_name + '_' + marker + '_' + str(index) + '_fits.txt'
        self.save_name = spectra_fname
        # self.parameters["files included"] = list(self.files)
        try:
            parameter_str = json.dumps(self.parameters, sort_keys=True, indent=4,
                                       separators=(',', ': '))
        except:
            print("Source: PMT.save_images\nJSON FAILED")
            print("Here is the dictionary that broke JSON:\n", self.parameters)
            return
        parameter_str = parameter_str.replace('\n', '\n#')

        num_lines = parameter_str.count(
            '#')  # Make the number of lines constant so importing is easier
        # for num in range(99 - num_lines): parameter_str += '\n#'
        parameter_str += '\n#' * (99 - num_lines)

        origin_import_spec = '\nNIR frequency,Signal,Standard error\neV,arb. u.,arb. u.\n,{:.3f},'.format(
            self.parameters["fieldStrength"]["mean"])
        spec_header = '#' + parameter_str + origin_import_spec

        origin_import_fits = '\nCenter energy,error,Amplitude,error,Linewidth,error\neV,,arb. u.,,eV,,\n,,'  # + marker
        fits_header = '#' + parameter_str + origin_import_fits

        for sideband in sorted(self.sb_dict.keys()):
            try:
                complete = np.vstack((complete, self.sb_dict[sideband]))
            except:
                complete = np.array(self.sb_dict[sideband])

        np.savetxt(os.path.join(folder_str, spectra_fname), complete, delimiter=',',
                   header=spec_header, comments='', fmt='%0.6e')

        try:
            np.savetxt(os.path.join(folder_str, fit_fname), self.sb_results,
                       delimiter=',',
                       header=fits_header, comments='', fmt='%0.6e')
        except AttributeError:
            # Catch the error that happens if you save something without files
            print("warning, couldn't save fit file (no sidebands found?)")

        print("Saved PMT spectrum.\nDirectory: {}".format(
            os.path.join(folder_str, spectra_fname)))

class TimeTrace(PMT):
    """
    This class will be able to handle time traces output by the PMT softare.
    """
    def __init__(self, file_path):
        super(HighSidebandPMT, self).__init__(file_path)

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

        Creates:
        self.fname = file name of the initial_CCD_piece
        self.sb_results = The sideband details from the initializing data
        self.parameters = The parameter dictionary of the initializing data.  May
                          not have all details of spectrum pieces added later.
        self.full_dict = a copy of the sb_results without the zeroth column, which
                         is SB order

        :param initial_CCD_piece: The starting part of the spectrum, often the lowest orders seen by CCD
        :type initial_CCD_piece: HighSidebandCCD
        :return: None
        """
        self.fname = initial_CCD_piece.fname
        try:
            self.sb_results = initial_CCD_piece.sb_results
        except AttributeError:
            print(initial_CCD_piece.full_dict)
            raise
        self.parameters = initial_CCD_piece.parameters
        self.parameters['files_here'] = [initial_CCD_piece.fname.split('/')[-1]]
        self.full_dict = {}
        for sb in self.sb_results:
            self.full_dict[sb[0]] = np.asarray(sb[1:])


    @staticmethod
    def parse_sb_array(arr):
        """
        Check to make sure the first even order sideband in an array is not weaker
        than the second even order. If this happens, it's likely because the SB was in
        the short pass filter and isn't work counting.

        We cut it out to prevent it from itnerfering with calculating overlaps
        :param arr:
        :return:
        """
        arr = np.array(arr)

        if (arr[0, sbarr.SBNUM]>0 and arr[1, sbarr.SBNUM]>0 and # make sure they're both pos
               arr[0, sbarr.AREA] < arr[1, sbarr.AREA]): # and the fact the area is less
            # print "REMOVING FIRST SIDEBAND FROM FULLSIDEBAND"
            # print arr[0]
            # print arr[1]
            arr = arr[1:]


        full_dict = {}
        for sb in arr:
            full_dict[sb[0]] = np.asarray(sb[1:])
        return full_dict, arr

    def add_CCD(self, ccd_object, verbose=False, force_calc=None, **kwargs):
        """
        This method will be called by the stitch_hsg_results function to add another
        CCD image to the spectrum.

        :param ccd_object: The CCD object that will be stiched into the current FullHighSideband object
        :type ccd_object: HighSidebandCCD
        :return: None
        """
        if self.parameters["gain"] == ccd_object.parameters["gain"]:
            calc = False
        else:
            calc = True
        if force_calc is not None:
            calc = force_calc
        if "need_ratio" in kwargs: #cascading it through, starting to think
            # everything should be in a kwarg
            calc = kwargs.pop("need_ratio")
        try:
            # self.full_dict = stitch_hsg_dicts(self.full_dict, ccd_object.full_dict,
            #                                   need_ratio=calc, verbose=verbose)
            self.full_dict = stitch_hsg_dicts(self, ccd_object, need_ratio=calc,
                                              verbose=verbose, **kwargs)
            self.parameters['files_here'].append(ccd_object.fname.split('/')[-1])
            # update sb_results, too
            sb_results = [[k]+list(v) for k, v in list(self.full_dict.items())]
            sb_results = np.array(sb_results)
            self.sb_results = sb_results[sb_results[:,0].argsort()]
        except AttributeError:
            print('Error, not enough sidebands to fit here! {}, {}, {}, {}'.format(
                self.parameters["series"], self.parameters["spec_step"],
                ccd_object.parameters["series"], ccd_object.parameters["spec_step"]
            ))

    def add_PMT(self, pmt_object, verbose=False):
        """
        This method will be called by the stitch_hsg_results function to add the PMT
        data to the spectrum.
        """
        # print "I'm adding PMT once"
        # self.full_dict = stitch_hsg_dicts(pmt_object.full_dict, self.full_dict,
                                          # need_ratio=True, verbose=False)
        self.full_dict = stitch_hsg_dicts(pmt_object, self,
                                          need_ratio=True, verbose=verbose)
        # if verbose:
        #     self.full_dict, ratio = self.full_dict
        # print "I'm done adding PMT data"
        self.parameters['files_here'].append(pmt_object.parameters['files included'])
        self.make_results_array()
        # if verbose:
        #     return ratio

    def make_results_array(self):
        """
        The idea behind this method is to create the sb_results array from the
        finished full_dict dictionary.
        """
        self.sb_results = None
        # print "I'm making the results array:", sorted(self.full_dict.keys())
        for sb in sorted(self.full_dict.keys()):
            # print "Going to add this", sb
            try:
                self.sb_results = np.vstack((self.sb_results, np.hstack((sb, self.full_dict[sb]))))
            except ValueError:
                # print "It didn't exist yet!"
                self.sb_results = np.hstack((sb, self.full_dict[sb]))
                # print "and I made this array:", self.sb_results[:, 0]

    def save_processing(self, file_name, folder_str, marker='', index='', verbose=''):
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
        except OSError as e:
            if e.errno == errno.EEXIST:
                pass
            else:
                raise

        temp = np.array(self.sb_results)

        ampli = np.array([temp[:, 3] / temp[:, 5]])  # I'm pretty sure this is
        # amplitude, not area
        temp[:, 5:7] = temp[:, 5:7] * 1000  # For meV linewidths
        if verbose:
            print("sb_results", self.sb_results.shape)
            print("ampli", ampli.shape)
        save_results = np.hstack((temp, ampli.T))

        # spectra_fname = file_name + '_' + marker + '_' + str(index) + '.txt'
        fit_fname = file_name + '_' + marker + '_' + str(index) + '_full.txt'
        # self.save_name = spectra_fname

        # self.parameters['addenda'] = self.addenda
        # self.parameters['subtrahenda'] = self.subtrahenda
        try:
            # PMT files add unnecessary number of lines, dump it into one line
            # by casting it to a string.
            reduced = self.parameters.copy()
            reduced["files_here"] = str(reduced["files_here"])
            parameter_str = json.dumps(reduced, sort_keys=True, indent=4, separators=(',', ': '))
        except Exception as e:
            print(e)
            print("Source: EMCCD_image.save_images\nJSON FAILED")
            print("Here is the dictionary that broke JSON:\n", self.parameters)
            return
        parameter_str = parameter_str.replace('\n', '\n#')

        num_lines = parameter_str.count('#')  # Make the number of lines constant so importing is easier
        # for num in range(99 - num_lines): parameter_str += '\n#'
        parameter_str += '\n#' * (99 - num_lines)
        # origin_import_spec = '\nNIR frequency,Signal,Standard error\neV,arb. u.,arb. u.'
        # spec_header = '#' + parameter_str + '\n#' + self.description[:-2] + origin_import_spec

        origin_import_fits = '\nSideband,Center energy,error,Sideband strength,error,Linewidth,error,Amplitude'+\
                             '\norder,eV,,arb. u.,,meV,,arb. u.\n' + ','.join([marker]*8)
        fits_header = '#' + parameter_str + origin_import_fits

        # np.savetxt(os.path.join(folder_str, spectra_fname), self.proc_data, delimiter=',',
        #           header=spec_header, comments='', fmt='%f')
        np.savetxt(os.path.join(folder_str, fit_fname), save_results, delimiter=',',
                   header=fits_header, comments='', fmt='%0.6e')

        if verbose:
            print("Save image.\nDirectory: {}".format(os.path.join(folder_str, fit_fname)))

####################
# Fitting functions
####################
def gauss(x, *p):
    """
    Gaussian fit function.

    :param x: The independent variable
    :type x: np.array, or int or float
    :param p: [mean, area, width, y offset] to be unpacked
    :type p: list of floats or ints
    :return: Depends on x, returns another np.array or float or int
    :rtype: type(x)
    """
    mu, A, sigma, y0 = p
    return (A / sigma) * np.exp(-(x - mu) ** 2 / (2. * sigma ** 2)) + y0

def lingauss(x, *p):
    """
    Gaussian fit function with a linear offset

    :param x: The independent variable
    :type x: np.array, or int or float
    :param p: [mean, area, width, constant offset of background, slope of background] to be unpacked
    :type p: list of floats or ints
    :return: Depends on x, returns another np.array or float or int
    :rtype: type(x)
    """
    mu, A, sigma, y0, m = p
    return (A / sigma) * np.exp(-(x - mu) ** 2 / (2. * sigma ** 2)) + y0 + m * x

def lorentzian(x, *p):
    """
    Lorentzian fit with constant offset

    :param x: The independent variable
    :type x: np.array, or int or float
    :param p: [mean, area, width, constant offset of background, slope of background] to be unpacked
    :type p: list of floats or ints
    :return: Depends on x, returns another np.array or float or int
    :rtype: type(x)
    """
    mu, A, gamma, y0 = p
    return (A / np.pi) * (gamma / ((x - mu) ** 2 + gamma ** 2)) + y0

def background(x, *p):
    """
    Arbitrary pink-noise model background data for absorbance FFT
    for the intention of replacing a peak in the FFT
    with the background

    :param x: The independent variable
    :type x: np.array, or int or float
    :param p: [proportionality factor, exponent of power law]
    :type p: list of floats or ints
    :return: Depends on x
    :rtype: type(x)
    """
    a, b = p
    return a * (1 / x) ** b

def gaussWithBackground(x, *p):
    """
    Gaussian with pink-noise background function

    :param x: independent variable
    :type x: np.array, or int or float
    :param p: [mean, area, width, constant background, proportionality of power law, exponent of power law]
    :type p: list of floats or ints
    :return: Depends on x
    :rtype: type(x)
    """
    pGauss = p[:4]
    a, b = p[4:]
    return gauss(x, *pGauss) + background(x, a, b)

####################
# Collection functions
####################
def hsg_combine_spectra(spectra_list, verbose = False, **kwargs):
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

    :param spectra_list: randomly-ordered list of HSG spectra, some of which can be stitched together
    :type spectra_list: List of HighSidebandCCD objects
    kwargs gets passed onto add_item
    :return: fully combined list of full hsg spectra.  No PMT business yet.
    :rtype: list of FullHighSideband
    """
    good_list = []
    spectra_list = spectra_list.copy()
    spectra_list.sort(key=lambda x: x.parameters["spec_step"])

    # keep a dict for each series' spec step
    # This allows you to combine spectra whose spec steps
    # change by values other than 1 (2, if you skip, or 0.5 if you
    # decide to insert things, or arbitary strings)
    spec_steps = {}

    for elem in spectra_list:
        # if verbose:
        #     print "Spec_step is", elem.parameters["spec_step"]
        current_steps = spec_steps.get(elem.parameters["series"], [])
        current_steps.append(elem.parameters["spec_step"])
        spec_steps[elem.parameters["series"]] = current_steps
    if verbose:
        print("I found these spec steps for each series:")
        print("\n\t".join("{}: {}".format(*ii) for ii in spec_steps.items()))

    # sort the list of spec steps
    for series in spec_steps:
        spec_steps[series].sort()

    same_freq = lambda x,y: x.parameters["fel_lambda"] == y.parameters["fel_lambda"]

    for index in range(len(spectra_list)):
        try:
            temp = spectra_list.pop(0)
            if verbose:
                print("\nStarting with this guy", temp, "\n")
        except:
            break

        good_list.append(FullHighSideband(temp))

        counter = 1
        temp_list = list(spectra_list)
        for piece in temp_list:
            if verbose:
                print("\tchecking this spec_step", piece.parameters["spec_step"], end=' ')
                print(", the counter is", counter)
            if not same_freq(piece, temp):
                if verbose:
                    print("\t\tnot the same fel frequencies ({} vs {})".format(piece.parameters["fel_lambda"], temp.parameters["fel_lambda"]))
                continue
            if temp.parameters["series"] == piece.parameters["series"]:
                if piece.parameters["spec_step"] == spec_steps[temp.parameters["series"]][counter]:
                    if verbose:
                        print("I found this one", piece)
                    counter += 1
                    good_list[-1].add_CCD(piece, verbose=verbose, **kwargs)
                    spectra_list.remove(piece)
                else:
                    print("\t\tNot the right spec step?", type(piece.parameters["spec_step"]))

            else:
                if verbose:
                    print("\t\tNot the same series ({} vs {}".format(
                        piece.parameters["series"],temp.parameters["series"]))
        good_list[-1].make_results_array()
    return good_list

def hsg_combine_spectra_arb_param(spectra_list, param_name="series", verbose = False):
    """
    This function is all about smooshing different parts of the same hsg
    spectrum together.  It takes a list of HighSidebandCCD spectra and turns the
    zeroth spec_step into a FullHighSideband object.  It then uses the function
    stitch_hsg_dicts over and over again for the smooshing.

    This is different than hsg_combine_spectra in that you pass which
    criteria distinguishes the files to be the "same". Since it can be any arbitrary
    value, things won't be exactly the same (field strength will never be identical
    between images). It will start with the first (lowest) spec step, then compare the
    number of images in the next step. Whichever has

    Input:
    spectra_list = list of HighSidebandCCD objects that have sideband spectra
                   larger than the spectrometer can see.

    Returns:
    good_list = A list of FullHighSideband objects that have been combined as
                much as can be.

    :param spectra_list: randomly-ordered list of HSG spectra, some of which can be stitched together
    :type spectra_list: list of HighSidebandCCD
    :return: fully combined list of full hsg spectra.  No PMT business yet.
    :rtype: list of FullHighSideband
    """
    if not spectra_list:
        raise RuntimeError("Passed an empty spectra list!")
    if isinstance(param_name, list):
        # if you pass two things because the param you want
        # is in a dict (e.g. field strength has mean/std)
        # do it that way
        param_name_list = list(param_name)
        paramGetter = lambda x: x.parameters[param_name_list[0]][param_name_list[1]]
        param_name = param_name[0]
    elif isinstance(spectra_list[0].parameters[param_name], dict):
        paramGetter = lambda x: x.parameters[param_name]["mean"]
    else:
        paramGetter = lambda x: x.parameters[param_name]

    good_list = []
    spectra_list.sort(key=lambda x: x.parameters["spec_step"])

    # keep a dict for each spec step.
    spec_steps = {}
    for elem in spectra_list:
        if verbose:
            print("Spec_step is", elem.parameters["spec_step"])
        current_steps = spec_steps.get(elem.parameters["spec_step"], [])
        current_steps.append(elem)
        spec_steps[elem.parameters["spec_step"]] = current_steps


    # Next, loop over all of the elements. For each element, if it has not
    # already been added to a spectra, look at all of the combinations from
    # other spec steps to figure out which has the smallest overall deviation
    # to make a new full spectrum
    good_list = []
    already_added = set()
    for elem in spectra_list:
        if elem in already_added: continue
        already_added.add(elem)
        good_list.append(FullHighSideband(elem))

        other_spec_steps = [v for k, v in list(spec_steps.items()) if
                            k != good_list[-1].parameters["spec_step"]]
        min_distance = np.inf
        cur_value = paramGetter(good_list[-1])
        best_match = None
        for comb in itt.product(*other_spec_steps):
            new_values = list(map(paramGetter, comb))
            all_values = new_values + [cur_value]

            if np.std(all_values) < min_distance:
                min_distance = np.std(all_values)
                best_match = list(comb)

        if best_match is None:
            raise RuntimeError("No matches found. Empty lists passed?")

        best_values = list(map(paramGetter, best_match))
        for spec in best_match:
            print("Adding new spec step\n\tStarted with spec={},series={}".format(
                good_list[-1].parameters["spec_step"],good_list[-1].parameters["series"]
            ))
            print("\tAdding with spec={},series={}\n".format(
                spec.parameters["spec_step"],
                spec.parameters["series"]
            ))
            print("\n\nfirst SBs:\n", good_list[-1].sb_results)
            print("\n\nsecond SBs:\n", spec.sb_results)
            good_list[-1].add_CCD(spec, True)
            print("\n\nEnding SBs:\n", good_list[-1].sb_results)

            already_added.add(spec)
        best_match.append(good_list[-1])
        best_values.append(cur_value)
        new_value = np.mean(best_values)
        new_std = np.std(best_values)

        if isinstance(good_list[-1].parameters[param_name], dict):
            best_values = np.array([x.parameters[param_name]["mean"] for x in best_match])
            best_std = np.array([x.parameters[param_name]["std"] for x in best_match])
            new_value = np.average(best_values, weights = best_std)
            new_std = np.sqrt(np.average((best_values-new_value)**2, weights=best_std))

        good_list[-1].parameters[param_name] = {
            "mean": new_value,
            "std": new_std
        }
    return good_list

def pmt_sorter(folder_path, plot_individual = True):
    """
    This function will be fed a folder with a bunch of PMT data files in it.
    The folder should contain a bunch of spectra with at least one sideband in
    them, each differing by the series entry in the parameters dictionary.

    This function will return a list of HighSidebandPMT objects.

    :param folder_path: Path to a folder containing a bunch of PMT data, can be
                        part of a parameter sweep
    :type folder_path: str
    :param plot_individual: Whether to plot each sideband itself
    :return: A list of all the possible hsg pmt spectra, organized by series tag
    :rtype: list of HighSidebandPMT
    """
    file_list = glob.glob(os.path.join(folder_path, '*[0-9].txt'))

    pmt_list = []

    plot_sb = lambda x: None

    if plot_individual:
        plt.figure("PMT data")
        def plot_sb(spec):
            spec = copy.deepcopy(spec)
            spec.process_sidebands()
            elem = spec.sb_dict[spec.initial_sb]
            plt.errorbar(elem[:, 0], elem[:, 1], elem[:, 2],
                     marker='o',
                     label="{} {}, {}.{} ".format(
                         spec.parameters["series"], spec.initial_sb,
                         spec.parameters["pm_hv"],
                         't' if spec.parameters.get("photon counted", False) else 'f')
                         )

    for sb_file in file_list:
        temp = HighSidebandPMT(sb_file)
        plot_sb(temp)
        try:
            for pmt_spectrum in pmt_list:  # pmt_spectrum is a pmt object
                if temp.parameters['series'] == pmt_spectrum.parameters['series']:
                    pmt_spectrum.add_sideband(temp)
                    break
            else:  # this will execute IF the break was NOT called
                pmt_list.append(temp)
        except:
            pmt_list.append(temp)
    # for sb_file in file_list:
    #     with open(sb_file,'rU') as f:
    #         param_str = ''
    #         line = f.readline()
    #         line = f.readline()
    #         while line[0] == '#':
    #             param_str += line[1:]
    #             line = f.readline()
    #
    #         parameters = json.loads(param_str)
    #     try:
    #         for pmt_spectrum in pmt_list: # pmt_spectrum is a pmt object?
    #             if parameters['series'] == pmt_spectrum.parameters['series']:
    #                 pmt_spectrum.add_sideband(sb_file)
    #                 break
    #         else: # this will execute IF the break was NOT called
    #             pmt_list.append(HighSidebandPMT(sb_file))
    #     except:
    #         pmt_list.append(HighSidebandPMT(sb_file))

    for pmt_spectrum in pmt_list:
        pmt_spectrum.process_sidebands()
    return pmt_list

def stitch_abs_results(main, new):
    raise NotImplementedError

def hsg_combine_qwp_sweep(path, loadNorm = True, save = False, verbose=False,
                          skipOdds = True):
    """
    Given a path to data taken from rotating the QWP (doing polarimetry),
    process the data (fit peaks), and parse it into a matrix of sb strength vs
    QWP angle vs sb number.

    By default, saves the file into "Processed QWP Dependence"

    Return should be passed directly into fitting

         -1     |     SB1     |   SB1  |     SB2     |   SB2  |    ...    |   ...  |     SBn     |   SBn  |
      angle1    | SB Strength | SB err | SB Strength | SB Err |
      angle2    |     ...     |    .   |
      .
      .
      .

    :param path: Path to load
    :param loadNorm: if true, load the normalized data
    :param save: Save the processed file or not
    :param verbose:
    :param skipOdds: Passed on to save sweep; determine whether or not to save
            odd orders. Generally, odds are artifacts and I don't want
            them messing up the data, so default to True.
    :return:
    """
    def getData(fname):
        """
        Helper function for loading the data and getting the header information for incident NIR stuff
        :param fname:
        :return:
        """
        if isinstance(fname, str):
            if loadNorm:
                ending = "_norm.txt"
            else:
                ending = "_snip.txt"
            header = ''
            with open(os.path.join("Processed QWP Dependence", fname + ending)) as fh:
                ln = fh.readline()
                while ln[0] == '#':
                    header += ln[1:]
                    ln = fh.readline()
            data = np.genfromtxt(os.path.join("Processed QWP Dependence", fname + ending),
                                 delimiter=',', dtype=str)
        if isinstance(fname, io.BytesIO):
            header = b''
            ln = fname.readline()
            while ln.decode()[0] == '#':
                header += ln[1:]
                ln = fname.readline()
            fname.seek(0)
            data = np.genfromtxt(fname,
                                 delimiter=',', dtype=str)

        header = json.loads(header)
        return data, float(header["lAlpha"]), float(header["lGamma"]), float(header["nir"]), float(header["thz"])
        ######### End getData

    try:
        sbData, lAlpha, lGamma, nir, thz = getData(path)
    except:
        # Do the processing on all the files
        specs = proc_n_plotCCD(path, keep_empties=True, verbose=verbose)

        for sp in specs:
            try:
                sp.parameters["series"] = round(float(sp.parameters["rotatorAngle"]), 2)
            except KeyError:
                # Old style of formatting
                sp.parameters["series"] = round(float(sp.parameters["detectorHWP"]), 2)
        specs = hsg_combine_spectra(specs, ignore_weaker_lowers=False)
        if not save:
            # If you don't want to save them, set everything up for doing Bytes objects
            # to replacing saving files
            full, snip, norm = io.BytesIO(), io.BytesIO(), io.BytesIO()
            if "nir_pola" not in specs[0].parameters:
                # in the olden days. Force them. Hopefully making them outside of ±360
                # makes it obvious
                specs[0].parameters["nir_pola"] = 361
                specs[0].parameters["nir_polg"] = 361
            keyName = "rotatorAngle"
            if keyName not in specs[0].parameters:
                # from back before I changed the name
                keyName = "detectorHWP"

            save_parameter_sweep(specs, [full, snip, norm], None,
                                     keyName, "deg", wanted_indices=[3, 4],
                                     header_dict={
                                         "lAlpha": specs[0].parameters["nir_pola"],
                                         "lGamma": specs[0].parameters["nir_polg"],
                                         "nir": specs[0].parameters["nir_lambda"],
                                         "thz": specs[0].parameters["fel_lambda"], },
                                 only_even=skipOdds)

            if loadNorm:
                sbData, lAlpha, lGamma, nir, thz = getData(norm)
            else:
                sbData, lAlpha, lGamma, nir, thz = getData(snip)
        else:
            save_parameter_sweep(specs, os.path.basename(path), "Processed QWP Dependence",
                                 "rotatorAngle", "deg", wanted_indices=[3, 4],
                                 header_dict={
                                     "lAlpha": specs[0].parameters["nir_pola"],
                                     "lGamma": specs[0].parameters["nir_polg"],
                                     "nir": specs[0].parameters["nir_lambda"],
                                     "thz": specs[0].parameters["fel_lambda"], },
                                 only_even=skipOdds)
            sbData, lAlpha, lGamma, nir, thz = getData(os.path.basename(path))

    laserParams = {
        "lAlpha": lAlpha,
        "lGamma": lGamma,
        "nir": nir,
        "thz": thz
    }

    # get which sidebands were found in this data set
    # first two rows are origin header, second is sideband number
    # (and empty strings, which is why the "if ii" below, to prevent
    # ValueErrors on int('').
    foundSidebands = np.array(sorted([float(ii) for ii in set(sbData[2]) if ii]))

    # Remove first 3 rows, which are strings for origin header, and cast it to floats
    sbData = sbData[3:].astype(float)

    # double the sb numbers (to account for sb strength/error) and add a dummy
    # number so the array is the same shape
    foundSidebands = np.insert(foundSidebands, range(len(foundSidebands)), foundSidebands)
    foundSidebands = np.insert(foundSidebands, 0, -1)
    return laserParams, np.row_stack((foundSidebands, sbData))

def makeCurve(eta, isVertical):
    """

    :param eta: QWP retardance at the wavelength
    :return:
    """
    cosd = lambda x: np.cos(x * np.pi / 180)
    sind = lambda x: np.sin(x * np.pi / 180)
    eta = eta * 2 * np.pi
    if isVertical:
        # vertical polarizer
        def analyzerCurve(x, *S):
            S0, S1, S2, S3 = S
            return S0-S1/2*(1+np.cos(eta)) \
                   + S3*np.sin(eta)*sind(2*x) \
                   + S1/2*(np.cos(eta)-1)*cosd(4*x) \
                   + S2/2*(np.cos(eta)-1)*sind(4*x)
    else:
        # vertical polarizer
        def analyzerCurve(x, *S):
            S0, S1, S2, S3 = S
            return S0+S1/2*(1+np.cos(eta)) \
                   - S3*np.sin(eta)*sind(2*x) \
                   + S1/2*(1-np.cos(eta))*cosd(4*x) \
                   + S2/2*(1-np.cos(eta))*sind(4*x)
    return analyzerCurve

def proc_n_fit_qwp_data(data, laserParams = dict(), wantedSBs = None, vertAnaDir = True, plot=False,
                        save = False, plotRaw = lambda sbidx, sbnum: False, series = '', eta=None,
                        **kwargs):
    """
    Fit a set of sideband data vs QWP angle to get the stoke's parameters
    :param data: data in the form of the return of hsg_combine_qwp_sweep
    :param laserParams: dictionary of the parameters of the laser, the angles and frequencies. See function for
                expected keys. I don't think the errors are used (except for plotting?), or the wavelengths (but
                left in for potential future use (wavelength dependent stuff?))
    :param wantedSBs: List of the wanted sidebands to fit out.
    :param vertAnaDir: direction of the analzyer. True if vertical, false if horizontal.
    :param plot: True/False to plot alpha/gamma/dop. Alternatively, a list of "a", "g", "d" to only plot selected ones
    :param save: filename to save the files. Accepts BytesIO
    :param plotRaw: callable that takes an index of the sb and sb number, returns true to plot the raw curve
    :param series: a string to be put in the header for the origin files
    :param eta: a function to call to calculate the desired retardance. Input will be the SB order.

    if saveStokes is in kwargs and False, it will not save the stokes parameters, since I rarely actually use them.
    :return:
    """
    defaultLaserParams = {
        "lAlpha": 90,
        "ldAlpha": 0.2,
        "lGamma": 0.0,
        "ldGamma": 0.2,
        "lDOP": 1,
        "ldDOP": 0.02,
        "nir": 765.7155,
        "thz": 21.1
    }
    defaultLaserParams.update(laserParams)
    lAlpha, ldAlpha, lGamma, ldGamma, lDOP, ldDOP = defaultLaserParams["lAlpha"], \
                                                    defaultLaserParams["ldAlpha"], \
                                                    defaultLaserParams["lGamma"], \
                                                    defaultLaserParams["ldGamma"], \
                                                    defaultLaserParams["lDOP"], \
                                                    defaultLaserParams["ldDOP"]
    allSbData = data
    angles = allSbData[1:, 0]

    # angles += -5
    # print("="*20)
    # print("\n"*3)
    # print("             WARNING")
    # print("\n"*3)
    # print("ANGLES HAVE BEEN MANUALLY OFFEST IN proc_n_fit_qwp_data")
    # print("\n"*3)
    # print("="*20)

    allSbData = allSbData[:, 1:] # trim out the angles

    if wantedSBs is None:
        # set to get rid of duplicates, 1: to get rid of the -1 used for
        # getting arrays the right shape
        wantedSBs = set(allSbData[0, 1:])

    if eta is None:
        """
        It might be easier for the end user to do this by passing eta(wavelength) instead of eta(sborder), 
        but then this function would need to carry around wavelengths, which is extra work. It could convert
        between NIR/THz wavelengths to SB order, but it's currently unclear whether you'd rather use what the WS6 
        claims, or what the sidebands say, and you'd probably want to take the extra step to ensure the SB fit rseults
        if using the spectromter wavelengths. In general, if you have a function as etal(wavelength), you'd probably 
        want to pass this as 
        eta = lambda x: etal(1239.84/(nirEv + x*THzEv))
        assuming nirEv/THzEv are the photon energies of the NIR/THz. 
        """
        eta = lambda x: 0.25

    # allow pasing a flag it ignore odds. I think I generally do, so set it to
    # default to True
    skipOdds = kwargs.get("skip_odds", True)

    # Make an array to keep all of the sideband information.
    # Start it off by keeping the NIR information (makes for easier plotting into origin)
    sbFits = [[0] + [-1] * 8 + [lAlpha, ldAlpha, lGamma, ldGamma, lDOP, ldDOP]]
    # Also, for convenience, keep a dictionary of the information.
    # This is when I feel like someone should look at porting this over to pandas
    sbFitsDict = {}
    sbFitsDict["S0"] = [[0, -1, -1]]
    sbFitsDict["S1"] = [[0, -1, -1]]
    sbFitsDict["S2"] = [[0, -1, -1]]
    sbFitsDict["S3"] = [[0, -1, -1]]
    sbFitsDict["alpha"] = [[0, lAlpha, ldAlpha]]
    sbFitsDict["gamma"] = [[0, lGamma, ldGamma]]
    sbFitsDict["DOP"] = [[0, lDOP, ldDOP]]

    # Iterate over all sb data. Skip by 2 because error bars are included
    for sbIdx in range(0, allSbData.shape[1], 2):
        sbNum = allSbData[0, sbIdx]
        if sbNum not in wantedSBs: continue
        if skipOdds and sbNum%2: continue
        # if verbose:
        #     print("\tlooking at sideband", sbNum)
        sbData = allSbData[1:, sbIdx]
        sbDataErr = allSbData[1:, sbIdx + 1]

        # try:
        #     p0 = sbFits[-1][1:8:2]
        # except:
        #     p0 = [1, 1, 0, 0]
        p0 = [1, 1, 0, 0]

        etan = eta(sbNum)
        try:
            p, pcov = curve_fit(makeCurve(etan, vertAnaDir), angles, sbData, p0=p0)
        except ValueError:
            # This is getting tossed around, especially when looking at noisy data,
            # especially with the laser line, and it's fitting erroneous values.
            # Ideally, I should be cutting this out and not even returning them,
            # but that's immedaitely causing
            p = np.nan*np.array(p0)
            pcov = np.eye(len(p))


        if plot and plotRaw(sbIdx, sbNum):
            # pg.figure("{}: sb {}".format(dataName, sbNum))
            plt.figure("All Curves")
            plt.errorbar(angles, sbData, sbDataErr, 'o-', name=f"{series}, {sbNum}")
            # plt.plot(angles, sbData,'o-', label="Data")
            fineAngles = np.linspace(angles.min(), angles.max(), 300)
            # plt.plot(fineAngles,
            #         makeCurve(eta, "V" in dataName)(fineAngles, *p0), name="p0")
            plt.plot(fineAngles,
                    makeCurve(etan, vertAnaDir)(fineAngles, *p))
            # plt.show()
            plt.ylim(0, 1)
            plt.xlim(0, 360)
            plt.ylabel("Normalized Intensity")
            plt.xlabel("QWP Angle (&theta;)")
            print(f"\t{series} {sbNum}, p={p}")


        # get the errors
        d = np.sqrt(np.diag(pcov))
        thisData = [sbNum] + list(p) + list(d)
        d0, d1, d2, d3 = d
        S0, S1, S2, S3 = p
        # reorder so errors are after values
        thisData = [thisData[i] for i in [0, 1, 5, 2, 6, 3, 7, 4, 8]]

        sbFitsDict["S0"].append([sbNum, S0, d0])
        sbFitsDict["S1"].append([sbNum, S1, d1])
        sbFitsDict["S2"].append([sbNum, S2, d2])
        sbFitsDict["S3"].append([sbNum, S3, d3])

        # append alpha value
        thisData.append(np.arctan2(S2, S1) / 2 * 180. / np.pi)
        # append alpha error
        variance = (d2 ** 2 * S1 ** 2 + d1 ** 2 * S2 ** 2) / (S1 ** 2 + S2 ** 2) ** 2
        thisData.append(np.sqrt(variance) * 180. / np.pi)

        sbFitsDict["alpha"].append([sbNum, thisData[-2], thisData[-1]])

        # append gamma value
        thisData.append(np.arctan2(S3, np.sqrt(S1 ** 2 + S2 ** 2)) / 2 * 180. / np.pi)
        # append gamma error
        variance = (d3 ** 2 * (S1 ** 2 + S2 ** 2) ** 2 + (d1 ** 2 * S1 ** 2 + d2 ** 2 * S2 ** 2) * S3 ** 2) / (
        (S1 ** 2 + S2 ** 2) * (S1 ** 2 + S2 ** 2 + S3 ** 2) ** 2)
        thisData.append(np.sqrt(variance) * 180. / np.pi)
        sbFitsDict["gamma"].append([sbNum, thisData[-2], thisData[-1]])

        # append degree of polarization
        thisData.append(np.sqrt(S1 ** 2 + S2 ** 2 + S3 ** 2) / S0)
        variance = ((d1 ** 2 * S0 ** 2 * S1 ** 2 + d0 ** 2 * (S1 ** 2 + S2 ** 2 + S3 ** 2) ** 2 + S0 ** 2 * (
        d2 ** 2 * S2 ** 2 + d3 ** 2 * S3 ** 2)) / (S0 ** 4 * (S1 ** 2 + S2 ** 2 + S3 ** 2)))
        thisData.append(np.sqrt(variance))
        sbFitsDict["DOP"].append([sbNum, thisData[-2], thisData[-1]])

        sbFits.append(thisData)

    sbFits = np.array(sbFits)
    sbFitsDict = {k: np.array(v) for k, v in sbFitsDict.items()}
    # This chunk used to insert the "alpha deviation", the difference between the angles and the
    # nir. I don't think I use this anymore, so stop saving it
                # origin_header = 'Sideband,S0,S0 err,S1,S1 err,S2,S2 err,S3,S3 err,alpha,alpha deviation,alpha err,gamma,gamma err,DOP,DOP err\n'
                # origin_header += 'Order,arb.u,arb.u,arb.u,arb.u,arb.u,arb.u,arb.u,arb.u,deg,deg,deg,deg,deg,arb.u.,arb.u.\n'
                # origin_header += 'Sideband,{},{},{},{},{},{},{},{},{},{},{},{},{},{},{}'.format(*["{}".format(series)] * 15)
                # sbFits = np.array(sbFits)
                # sbFits = np.insert(sbFits, 10, sbFits[:, 9] - lAlpha, axis=1)
                # sbFits = sbFits[sbFits[:, 0].argsort()]

    origin_header = "#\n"*100 # to fit all other files for easy origin importing
    origin_header += 'Sideband,S0,S0 err,S1,S1 err,S2,S2 err,S3,S3 err,alpha,alpha err,gamma,gamma err,DOP,DOP err\n'
    origin_header += 'Order,arb.u,arb.u,arb.u,arb.u,arb.u,arb.u,arb.u,arb.u,deg,deg,deg,deg,arb.u.,arb.u.\n'
    origin_header += 'Sideband,{},{},{},{},{},{},{},{},{},{},{},{},{},{}'.format(*["{}".format(series)] * 14)
    sbFits = sbFits[sbFits[:, 0].argsort()]

    if isinstance(save, str):
        sbFitsSave = sbFits
        if not kwargs.get("saveStokes", True):
            headerlines = origin_header.splitlines()
            ln, units, coms = headerlines[-3:]
            ln = ','.join([ln.split(',')[0]] + ln.split(',')[9:])
            units = ','.join([units.split(',')[0]] + units.split(',')[9:])
            coms = ','.join([coms.split(',')[0]] + coms.split(',')[9:])
            headerlines[-3:] = ln, units, coms
            # remove them from the save data
            origin_header = '\n'.join(headerlines)
            sbFitsSave = np.delete(sbFits, range(1, 9), axis=1)

        if not os.path.exists(os.path.dirname(save)):
            os.mkdir(os.path.dirname(save))
        np.savetxt(save, np.array(sbFitsSave), delimiter=',', header=origin_header,
                   comments='', fmt='%.6e')

    # print("a = {:.2f} ± {:.2f}".format(sbFits[1, 9], sbFits[1, 10]))
    # print("g = {:.2f} ± {:.2f}".format(sbFits[1, 11], sbFits[1, 12]))

    if plot:
        plt.figure("alpha")
        plt.errorbar(sbFitsDict["alpha"][:, 0],
                     sbFitsDict["alpha"][:, 1],
                     sbFitsDict["alpha"][:, 2],
                     'o-', name = series
                     )
        plt.figure("gamma")
        plt.errorbar(sbFitsDict["gamma"][:, 0],
                     sbFitsDict["gamma"][:, 1],
                     sbFitsDict["gamma"][:, 2],
                     'o-', name=series
                     )
    return sbFits, sbFitsDict

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
        print("shape of median filter:", med.shape)
    meanMedian = med.mean(axis=1)
    # meanMedian = med.copy()
    if debugging:
        print("shape of meaned median filter:", meanMedian.shape)
    # Construct a cutoff for each pixel. It was kind of guess and
    # check
    cutoff = meanMedian * medianRatio + noiseCoeff * np.std(meanMedian[-100:])
    if debugging:
        print("shape of cutoff criteria:", cutoff.shape)
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

        p2.plot(np.sum(med, axis=1) / d.shape[1])
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
        for i, v in enumerate((d - med).T):
            p2.plot(v, pen=(i, d.shape[1]), name=str(i))
        p2.plot(cutoff, pen=pg.mkPen('w', width=3))
        win2.show()
        winlist.append(win2)

    # Find the bad pixel positions
    # Note the [:, None] - needed to cast the correct shapes
    badPixs = np.argwhere((d - med) > (cutoff.reshape(len(cutoff), 1)))

    for pix in badPixs:
        # get the other pixels in the row which aren't the cosmic
        if debugging:
            print("cleaning pixel", pix)
        p = d[pix[0], [i for i in range(d.shape[1]) if not i == pix[1]]]
        if debugging:
            print("\tRemaining pixels in row are", p)
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
    if max(first[:, 0]) > max(second[:, 0]):
        flipped = True
        first, second = second, first

def stitch_hsg_dicts(full_obj, new_obj, need_ratio=False, verbose=False, ratios=[1,1],
                     override_ratio = False, ignore_weaker_lowers = True):
    """
    This helper function takes a FullHighSideband and a sideband
    object, either CCD or PMT and smushes the new sb_results into the full_dict.

    The first input doesn't change, so f there's a PMT set of data involved, it
    should be in the full variable to keep the laser normalization intact.

    This function almost certainly does not work for stitching many negative orders
    in it's current state

    11/14/16
    --------
    This function has been updated to take the CCD objects themselves to be more
    intelligent about stitching. Consider two scans, (a) spec step 0 with 1 gain, spec
    step 2 with 110 gain and (b) spec step 0 with 50 gain and spec step 1 with 110 gain.
    The old version would always take spec step 0 to scale to, so while comparisons
    between spec step 0 and 1 for either case is valid, comparison between (a) and (b)
    were not, since they were scaled to different gain parameters. This new code will
    check what the gain values are and scale to the 110 data set, if present. This seems
    valid because we currently always have a 110 gain exposure for higher order
    sidebands.
    The exception is if the laser is present (sideband 0), as that is an absolute
    measure to which all else should be related.
    TODO: run some test cases to test this.

    06/11/18
    --------
    That sometimes was breaking if there were only 3-4 sidebands to fit with poor
    SNR. I've added the override_ratio to be passed to set a specific ratio to scale
    by. From data on 06/03/18, the 50gain to 110gain is a ~3.6 ratio. I haven't done
    a clean way of specifying which data set it should be scaled. Right now,
    it leaves the laser line data, or the 110 gain data alone.


    Inputs:
    full = full_dict from FullHighSideband, or HighSidebandPMT.  It's important
           that it contains lower orders than the new_dict.
    new_dict = another full_dict.
    need_ratio = If gain or other parameters aren't equal and must resort to
                 calculating the ratio instead of the measurements being equivalent.
                 Changing integration time still means N photons made M counts,
                 but changing gain or using PMT or whatever does affect things.
    ratios: Will update with the values to the ratios needed to scale the data.
            ratios[0] is the ratio for the "full_obj"
            ratios[1] is the ratio for the "new_obj"
            one of them will be one, one will be the appropriate scale, since one of
            them is unscaled. This is strictly speaking an output
    override_ratio: Pass a float to specify the ratio that should be used.
    ignore_weaker_lowers: Sometimes, a SB is in the short pass filter so a lower
        order is weaker than the next highest. If True, causes script to ignore all
        sidebands which are weaker and lower order.

    Returns:
    full = extended version of the input full.  Overlapping sidebands are
           averaged because that makes sense?
    """
    if isinstance(full_obj, dict) and isinstance(new_obj, dict):
        return stitch_hsg_dicts_old(full_obj, new_obj, need_ratio, verbose)

    if verbose:
        print("=" * 15)
        print()
        print("Stitching HSG dicts")
        print()
        print("=" * 15)

    # remove potentially offensive SBs, i.e. a 6th order SB being in the SPF for more
    #  data, but being meaningless to pull intensity information from.
    # Note: this might not be the best if you get to higher order stitches where it's
    #  possible that the sidebands might not be monotonic (from noise?)
    if ignore_weaker_lowers:
        full_obj.full_dict, full_obj.sb_results = FullHighSideband.parse_sb_array(full_obj.sb_results)
        new_obj.new_dict, new_obj.sb_results = FullHighSideband.parse_sb_array(new_obj.sb_results)

    # was fucking around with references and causing updates to arrays when it shouldn't
    # be
    full = copy.deepcopy(full_obj.full_dict)
    new_dict = copy.deepcopy(new_obj.full_dict)

    # Force a rescaling if you've passed a specified parameter
    # if isinstance(override_ratio, float):
    #     need_ratio = True

    # Do some testing to see which dict should be scaled to the other
    # I honestly forget why I prioritized the PMT first like this. But the third
    # check looks to make a gain 110 prioritize non-110, unless the non-110 includes
    # a laser line
    scaleTo = ""
    if need_ratio:
        if isinstance(new_obj, HighSidebandPMT):
            scaleTo = "new"
        elif isinstance(full_obj, HighSidebandPMT):
            scaleTo = "full"
        elif new_obj.parameters["gain"] == 110 and full_obj.parameters["gain"] != 110 \
            and 0 not in full:
            scaleTo = "new"
        else:
            scaleTo = "full"

    if verbose:
        print("\tI'm adding these sidebands", sorted(new_dict.keys()))
        print("\t  With these:", sorted(full.keys()))
    overlap = [] # The list that hold which orders are in both dictionaries
    missing = [] # How to deal with sidebands that are missing from full but in new.
    for new_sb in sorted(new_dict.keys()):
        full_sbs = sorted(full.keys())
        if new_sb in full_sbs:
            overlap.append(new_sb)
        elif new_sb not in full_sbs and new_sb < full_sbs[-1]:
            # This probably doesn't work with bunches of negative orders
            missing.append(new_sb)

    if verbose:
        print("\t  ( overlap:", overlap, ")")
        print("\t  ( missing:", missing, ")")


    # This if-else clause handles how to average together overlapping sidebands
    # which are seen in both spectra,
    if need_ratio:
        # Calculate the appropriate ratio to multiply the new sidebands by.
        # I'm not entirely sure what to do with the error of this guy.
        ratio_list = []
        try:
            new_starter = overlap[-1]
            if verbose:
                print("\n\tadding these ratios,", end=' ')
            if len(overlap) > 2:
                overlap = [x for x in overlap if (x % 2 == 0)
                           ]# and (x != min(overlap) and (x != max(overlap)))]
            if scaleTo == "new":
                if verbose:
                    print("scaling to new :")
                for sb in overlap:
                    ratio_list.append(new_dict[sb][2]/full[sb][2])
                    if verbose:
                        print("\t\t{:2.0f}: {:.3e}/{:.3e} ~ {:.3e},".format(sb, new_dict[sb][2],
                                                               full[sb][2], ratio_list[-1]))
                # new_ratio = 1 06/11/18 Not sure what these were used for
                ratio = np.mean(ratio_list)
            else:
                if verbose:
                    print("scaling to full:")
                for sb in overlap:
                    ratio_list.append(full[sb][2] / new_dict[sb][2])
                    if verbose:
                        print("\t\t{:2.0f}: {:.3e}/{:.3e} ~ {:.3e},".format(sb, full[sb][2],
                                                               new_dict[sb][2], ratio_list[-1]))
                # new_ratio = np.mean(ratio_list) 06/11/18 Not sure what these were used for

                ratio = np.mean(ratio_list)
            # Maybe not the best way to do it, performance wise, since you still
            # iterate through the list, even though you'll override it.
            if isinstance(override_ratio, float):
                ratio = override_ratio
                if verbose:
                    print("overriding calculated ratio with user inputted")
            error = np.std(ratio_list) / np.sqrt(len(ratio_list))

        except IndexError:
            # If there's no overlap (which you shouldn't let happen), hardcode a ratio
            # and error. I looked at all the ratios for the overlaps from 6/15/16
            # (540ghz para) to get the rough average. Hopefully they hold for all data.
            if not overlap:
                ratio = 0.1695
                error = 0.02
                # no overlap, so make sure it grabs all the sidebands
                new_starter = min(new_dict.keys())
            else:
                raise
        if verbose:
            # print "Ratio list\n\t", ("{:.3g}, "*len(ratio_list))[:-2].format(*ratio_list)
            # print "Overlap   \n\t", [round(ii, 3) for ii in overlap]
            print("\t Ratio: {:.3g} +- {:.3g} ({:.2f}%)\n".format(ratio, error, error/ratio*100))
        # Adding the new sidebands to the full set and moving errors around.
        # I don't know exactly what to do about the other aspects of the sidebands
        # besides the strength and its error.
        if scaleTo == "full":
            ratios[1] = ratio
            for sb in overlap:
                if verbose:
                    print("For SB {:02d}, original strength is {:.3g} +- {:.3g} ({:.3f}%)".format(int(sb), new_dict[sb][2], new_dict[sb][3],
                                                                        new_dict[sb][3]/new_dict[sb][2]*100
                            ))

                new_dict[sb][3] = ratio * new_dict[sb][2] * np.sqrt((error / ratio) ** 2 + (new_dict[sb][3] / new_dict[sb][2]) ** 2)
                new_dict[sb][2] = ratio * new_dict[sb][2]
                if verbose:
                    print("\t\t   scaled\t\t\t\t{:.3g} +- {:.3g} ({:.3f}%)".format(new_dict[sb][2],
                                                                        new_dict[sb][3],
                                                                        new_dict[sb][3]/new_dict[sb][2]*100))
                    print("\t\t   full\t\t\t\t\t{:.3g} +- {:.3g} ({:.3f}%)".format(full[sb][2],
                                                                        full[sb][3],
                                                                        full[sb][3]/full[sb][2]*100))


                sb_error = np.sqrt(full[sb][3] ** (-2) + new_dict[sb][3] ** (-2)) ** (-1)

                avg = (full[sb][2] / (full[sb][3] ** 2) + new_dict[sb][2] / (
                    new_dict[sb][3] ** 2)) / (full[sb][3] ** (-2) + new_dict[sb][3] ** (-2))
                full[sb][2] = avg
                full[sb][3] = sb_error
                if verbose:
                    print("\t\t   replaced with \t\t{:.3g} +- {:.3g} ({:.3f}%)".format(full[sb][2],
                                                                        full[sb][3],
                                                                        full[sb][3]/full[sb][2]*100))
                    print()

                lw_error = np.sqrt(full[sb][5] ** (-2) + new_dict[sb][5] ** (-2)) ** (-1)
                lw_avg = (full[sb][4] / (full[sb][5] ** 2) + new_dict[sb][4] / (
                new_dict[sb][5] ** 2)) / (
                             full[sb][5] ** (-2) + new_dict[sb][5] ** (-2))
                full[sb][4] = lw_avg
                full[sb][5] = lw_error  # This may not be the exactly right way to calculate the error
        else:
            ratios[0] = ratio
            for sb in overlap:
                full[sb][3] = ratio * full[sb][2] * np.sqrt((error / ratio) ** 2 + (full[sb][3] / full[sb][2]) ** 2)
                full[sb][2] = ratio * full[sb][2]

                sberror = np.sqrt(full[sb][3] ** (-2) + new_dict[sb][3] ** (-2)) ** (-1)
                avg = (full[sb][2] / (full[sb][3] ** 2) + new_dict[sb][2] / (
                    new_dict[sb][3] ** 2)) / (full[sb][3] ** (-2) + new_dict[sb][3] ** (-2))
                full[sb][2] = avg
                full[sb][3] = sberror

                lw_error = np.sqrt(full[sb][5] ** (-2) + new_dict[sb][5] ** (-2)) ** (-1)
                lw_avg = (full[sb][4] / (full[sb][5] ** 2) + new_dict[sb][4] / (
                new_dict[sb][5] ** 2)) / (
                             full[sb][5] ** (-2) + new_dict[sb][5] ** (-2))
                full[sb][4] = lw_avg
                full[sb][5] = lw_error  # This may not be the exactly right way to calculate the error


    else: # not needing a new ratio
        try:
            new_starter = overlap[-1] # This grabs the sideband order where only the new dictionary has
                                      # sideband information.  It's not clear why it necessarily has to be
                                      # at this line.
            overlap = [x for x in overlap if (x % 2 == 0)
                       ] # and (x != min(overlap) and (x != max(overlap)))]
            # This cuts out the lowest order sideband in the overlap for mysterious reasons
            for sb in overlap: # This for loop average two data points weighted by their relative errors
                if verbose:
                    print("The sideband", sb)
                    print("Old value", full[sb][4] * 1000)
                    print("Add value", new_dict[sb][4] * 1000)
                try:
                    error = np.sqrt(full[sb][3] ** (-2) + new_dict[sb][3] ** (-2)) ** (-1)
                    avg = (full[sb][2] / (full[sb][3] ** 2) + new_dict[sb][2] / (new_dict[sb][3] ** 2)) / (
                        full[sb][3] ** (-2) + new_dict[sb][3] ** (-2))
                    full[sb][2] = avg
                    full[sb][3] = error
                except RuntimeWarning:
                    raise IOError()

                lw_error = np.sqrt(full[sb][5] ** (-2) + new_dict[sb][5] ** (-2)) ** (-1)
                lw_avg = (full[sb][4] / (full[sb][5] ** 2) + new_dict[sb][4] / (new_dict[sb][5] ** 2)) / (
                full[sb][5] ** (-2) + new_dict[sb][5] ** (-2))
                full[sb][4] = lw_avg
                full[sb][5] = lw_error  # This may not be the exactly right way to calculate the error
                if verbose:
                    print("New value", lw_avg * 1000)
        except:
            new_starter = 0  # I think this makes things work when there's no overlap
    if verbose:
        print("appending new elements. new_starter={}".format(new_starter))


    for sb in [x for x in list(new_dict.keys()) if ((x > new_starter) or (x in missing))]:
        full[sb] = new_dict[sb]
        if scaleTo == "full":
            full[sb][2] = ratio * full[sb][2]
            full[sb][3] = full[sb][2] * np.sqrt((error / ratio) ** 2 + (ratio * full[sb][3] / full[sb][2]) ** 2)
    if scaleTo == "new":
        for sb in set(full.keys()) - set(sorted(new_dict.keys())[:]):
            full[sb][2] *= ratio
            # TODO: I think this is an invalid error
            # propagation (since ratio has error associated with it
            full[sb][3] *= ratio
    if verbose:
        print("I made this dictionary", sorted(full.keys()))
        print('-'*19)
        return full
        return full, ratio #the fuck? Why was this here?

    return full

def stitch_hsg_dicts_old(full, new_dict, need_ratio=False, verbose=False):
    """
    This helper function takes a FullHighSideband.full_dict attribute and a sideband
    object, either CCD or PMT and smushes the new sb_results into the full_dict.

    The first input doesn't change, so f there's a PMT set of data involved, it
    should be in the full variable to keep the laser normalization intact.

    This function almost certainly does not work for stitching many negative orders
    in it's current state

    11/14/16
    --------
    The original function has been updated to take the full object (instead of
    the dicts alone) to better handle calculating ratios when stitching. This is called
    once things have been parsed in the original function (or legacy code where dicts
    are passed instead of the object)

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
    if verbose:
        print("I'm adding these sidebands in old stitcher", sorted(new_dict.keys()))
    overlap = [] # The list that hold which orders are in both dictionaries
    missing = [] # How to deal with sidebands that are missing from full but in new.
    for new_sb in sorted(new_dict.keys()):
        full_sbs = sorted(full.keys())
        if new_sb in full_sbs:
            overlap.append(new_sb)
        elif new_sb not in full_sbs and new_sb < full_sbs[-1]:
            # This probably doesn't work with bunches of negative orders
            missing.append(new_sb)

    if verbose:
        print("overlap:", overlap)
        print("missing:", missing)

    # This if-else clause handles how to average together overlapping sidebands
    # which are seen in both spectra,
    if need_ratio:
        # Calculate the appropriate ratio to multiply the new sidebands by.
        # I'm not entirely sure what to do with the error of this guy.
        ratio_list = []
        #print '\n1979\nfull[2]', full[0][2]
        try:
            new_starter = overlap[-1]
            if len(overlap) > 2:
                overlap = [x for x in overlap if (x % 2 == 0)
                           ]#and (x != min(overlap) and (x != max(overlap)))]
            for sb in overlap:
                ratio_list.append(full[sb][2] / new_dict[sb][2])
            ratio = np.mean(ratio_list)
            # print
            # print '-'*15
            # print "ratio for {}: {}".format()
            error = np.std(ratio_list) / np.sqrt(len(ratio_list))
        except IndexError:
            # If there's no overlap (which you shouldn't let happen),
            # hardcode a ratio and error.
            # I looked at all the ratios for the overlaps from 6/15/16
            # (540ghz para) to get the rough average. Hopefully they hold
            # for all data.
            if not overlap:
                ratio = 0.1695
                error = 0.02
                # no overlap, so make sure it grabs
                # all the sidebands
                new_starter = min(new_dict.keys())
            else:
                raise
        if verbose:
            print("Ratio list","\n", [round(ii, 3) for ii in ratio_list])
            print("Overlap   ","\n", [round(ii, 3) for ii in overlap])
            print("Ratio", ratio)
            print("Error", error)
        #print '\n2118\nfull[2]', full[0][2]
        # Adding the new sidebands to the full set and moving errors around.
        # I don't know exactly what to do about the other aspects of the sidebands
        # besides the strength and its error.
        for sb in overlap:
            full[sb][2] = ratio * new_dict[sb][2]
            full[sb][3] = full[sb][2] * np.sqrt((error / ratio) ** 2 + (new_dict[sb][3] / new_dict[sb][2]) ** 2)
            #print '\n2125\nfull[2]', full[0][3]
            # Now for linewidths
            lw_error = np.sqrt(full[sb][5] ** (-2) + new_dict[sb][5] ** (-2)) ** (-1)
            lw_avg = (full[sb][4] / (full[sb][5] ** 2) + new_dict[sb][4] / (new_dict[sb][5] ** 2)) / (
            full[sb][5] ** (-2) + new_dict[sb][5] ** (-2))
            full[sb][4] = lw_avg
            full[sb][5] = lw_error
        #print '\n2132\nfull[2]', full[0][2]
    else:
        try:
            new_starter = overlap[-1] # This grabs the sideband order where only the new dictionary has
                                      # sideband information.  It's not clear why it necessarily has to be
                                      # at this line.
            overlap = [x for x in overlap if (x % 2 == 0) and (x != min(overlap) and (x != max(overlap)))]
            # This cuts out the lowest order sideband in the overlap for mysterious reasons
            for sb in overlap: # This for loop average two data points weighted by their relative errors
                if verbose:
                    print("The sideband", sb)
                    print("Old value", full[sb][4] * 1000)
                    print("Add value", new_dict[sb][4] * 1000)
                error = np.sqrt(full[sb][3] ** (-2) + new_dict[sb][3] ** (-2)) ** (-1)
                avg = (full[sb][2] / (full[sb][3] ** 2) + new_dict[sb][2] / (new_dict[sb][3] ** 2)) / (
                    full[sb][3] ** (-2) + new_dict[sb][3] ** (-2))
                full[sb][2] = avg
                full[sb][3] = error

                lw_error = np.sqrt(full[sb][5] ** (-2) + new_dict[sb][5] ** (-2)) ** (-1)
                lw_avg = (full[sb][4] / (full[sb][5] ** 2) + new_dict[sb][4] / (new_dict[sb][5] ** 2)) / (
                full[sb][5] ** (-2) + new_dict[sb][5] ** (-2))
                full[sb][4] = lw_avg
                full[sb][5] = lw_error  # This may not be the exactly right way to calculate the error
                if verbose:
                    print("New value", lw_avg * 1000)
        except:

            new_starter = 0  # I think this makes things work when there's no overlap
    if verbose:
        print("appending new elements. new_starter={}".format(new_starter))

    # This loop will add the sidebands which were only seen in the second step
    for sb in [x for x in list(new_dict.keys()) if ((x >= new_starter) or (x in missing))]:
        full[sb] = new_dict[sb]
        if need_ratio:
            full[sb][2] = ratio * full[sb][2]
            full[sb][3] = full[sb][2] * np.sqrt((error / ratio) ** 2 + (ratio * full[sb][3] / full[sb][2]) ** 2)
            #print '\n2164\nfull[2]', full[0][2]
    if verbose:
        print("I made this dictionary", sorted(full.keys()))
    return full

def save_parameter_sweep_no_sb(spectrum_list, file_name, folder_str, param_name, unit,
                         verbose=False):
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
        sb_included = sorted(list(set(sb_included + list(spec.full_dict.keys()))))
        included_spectra[spec.fname.split('/')[-1]] = spec.parameters[param_name]
        # If these are from summed spectra, then only the the first file name
        # from that sum will show up here, which should be fine?
    if verbose:
        # print "full name:", spectrum_list[0].fname
        print("included names:", included_spectra)
        print("sb_included:", sb_included)

    for spec in spectrum_list:
        temp_dict = {}  # This is different from full_dict in that the list has the
        # sideband order as the zeroth element.
        if verbose:
            print("the sb_results:", spec.sb_results)
        if spec.sb_results.ndim == 1: continue
        for index in range(len(spec.sb_results[:, 0])):
            if verbose:
                print("my array slice:", spec.sb_results[index, :])
            temp_dict[int(round(spec.sb_results[index, 0]))] = np.array(
                spec.sb_results[index, 1:])

        if verbose:
            print(temp_dict)

        for sb in sb_included:
            blank = np.zeros(6)
            # print "checking sideband order:", sb
            # print "blank", blank
            if sb not in temp_dict:
                # print "\nNeed to add sideband order:", sb
                temp_dict[sb] = blank
        try:  # Why is this try-except here?
            spec_data = np.array([float(spec.parameters[param_name])])
        except:
            spec_data = np.array([float(spec.parameters[param_name][:2])])
        for key in sorted(temp_dict.keys()):
            # print "I am going to hstack this:", temp_dict[key]
            spec_data = np.hstack((spec_data, temp_dict[key]))

        try:
            param_array = np.vstack((param_array, spec_data))
        except:
            param_array = np.array(spec_data)
        if verbose:
            print("The shape of the param_array is:", param_array.shape)
            # print "The param_array itself is:", param_array
    '''
    param_array_norm = np.array(param_array).T # python iterates over rows
    for elem in [x for x in xrange(len(param_array_norm)) if (x-1)%7 == 3]:
        temp_max = np.max(param_array_norm[elem])
        param_array_norm[elem] = param_array_norm[elem] / temp_max
        param_array_norm[elem + 1] = param_array_norm[elem + 1] / temp_max
    '''
    snipped_array = param_array[:, 0]
    norm_array = param_array[:, 0]
    if verbose:
        print("Snipped_array is", snipped_array)
    for ii in range(len(param_array.T)):
        if (ii - 1) % 6 == 0:
            if verbose:
                print("param_array shape", param_array[:, ii])
            snipped_array = np.vstack((snipped_array, param_array[:, ii]))
            norm_array = np.vstack((norm_array, param_array[:, ii]))
        elif (ii - 1) % 6 == 2:
            snipped_array = np.vstack((snipped_array, param_array[:, ii]))

            temp_max = np.max(param_array[:, ii])
            norm_array = np.vstack((norm_array, param_array[:, ii] / temp_max))
        elif (ii - 1) % 6 == 3:
            snipped_array = np.vstack((snipped_array, param_array[:, ii]))
            norm_array = np.vstack((norm_array, param_array[:, ii] / temp_max))

    snipped_array = snipped_array.T
    norm_array = norm_array.T

    try:
        os.mkdir(folder_str)
    except OSError as e:
        if e.errno == errno.EEXIST:
            pass
        else:
            raise
    norm_name = file_name + '_norm.txt'
    snip_name = file_name + '_snip.txt'
    file_name = file_name + '.txt'

    try:
        included_spectra_str = json.dumps(included_spectra, sort_keys=True, indent=4,
                                          separators=(',', ': '))
    except:
        print("Source: save_parameter_sweep\nJSON FAILED")
        return
    included_spectra_str = included_spectra_str.replace('\n', '\n#')

    included_spectra_str += '\n#' * (99 - included_spectra_str.count('\n'))
    origin_import1 = param_name
    origin_import2 = unit
    origin_import3 = ""
    for order in sb_included:
        origin_import1 += "Frequency,error,Sideband strength,error,Linewidth,error"
        origin_import2 += ",eV,,arb. u.,,meV,"
        origin_import3 += ",{0},,{0},,{0},".format(order)
    origin_total = origin_import1 + "\n" + origin_import2 + "\n" + origin_import3

    origin_import1 = param_name
    origin_import2 = unit
    origin_import3 = ""
    for order in sb_included:
        origin_import1 += ",Frequency,Sideband strength,error"
        origin_import2 += ",eV,arb. u.,"
        origin_import3 += ",{0},{0},".format(order)
    origin_snip = origin_import1 + "\n" + origin_import2 + "\n" + origin_import3

    header_total = '#' + included_spectra_str + '\n' + origin_total
    header_snip = '#' + included_spectra_str + '\n' + origin_snip

    # print "Spec header: ", spec_header
    if verbose:
        print("the param_array is:", param_array)
    np.savetxt(os.path.join(folder_str, file_name), param_array, delimiter=',',
               header=header_total, comments='', fmt='%0.6e')
    np.savetxt(os.path.join(folder_str, snip_name), snipped_array, delimiter=',',
               header=header_snip, comments='', fmt='%0.6e')
    np.savetxt(os.path.join(folder_str, norm_name), norm_array, delimiter=',',
               header=header_snip, comments='', fmt='%0.6e')
    if verbose:
        print("Saved the file.\nDirectory: {}".format(
            os.path.join(folder_str, file_name)))

def save_parameter_sweep(spectrum_list, file_name, folder_str, param_name, unit,
                         wanted_indices = [1, 3, 4], skip_empties = False, verbose=False,
                         header_dict = {}, only_even=False):
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


    Thus function has been update to pass a list of indices to slice for the return
    values

    skip_empties: If False, will add a row of zeroes for the parameter even if no sidebands
    are found. If True, will not add a line for that parameter

    only_even: don't include odd orders in the saved sweep

    [sb number, Freq (eV), Freq error (eV), Gauss area (arb.), Area error, Gauss linewidth (eV), Linewidth error (eV)]
    [    0    ,      1   ,        2,      ,        3         ,      4    ,         5           ,        6            ]
    """
    if isinstance(param_name, list):
        # if you pass two things because the param you want
        # is in a dict (e.g. field strength has mean/std)
        # do it that way
        param_name_list = list(param_name) # keep reference to old one
        paramGetter = lambda x: x.parameters[param_name_list[0]][param_name_list[1]]
        # Keep the name for labeling things later on
        param_name = param_name[0]
    else:
        paramGetter = lambda x: x.parameters[param_name]

    # Sort all of the spectra based on the desired key
    spectrum_list.sort(key=paramGetter)

    # keep track of which file name corresponds to which parameter which gets put in
    included_spectra = dict()

    # The big array which will be stacked up to keep all of the sideband details vs desired parameter
    param_array = None

    # list of which sidebands are seen throughout.
    sb_included = []
    # how many parameters (area, strength, linewidth, pos, etc.) are there?
    # Here incase software changes and more things are kept in
    # sb results. Needed to handle how to slice the arrays
    try:
        num_params = spectrum_list[0].sb_results.shape[1]
    except IndexError:
        # There's a file with only 1 sb and it happens to be first
        # in the list.
        num_params = spectrum_list[0].sb_results.shape[0]
    except AttributeError:
        # The first file has no sidebands, so just hardcode it, as stated below.
        num_params=0

    # Rarely, there's an issue where I'm doing some testing and there's a set
    # where the first file has no sidebands in it, so the above thing returns 0
    # It seems really silly to do a bunch of testing to try and correct for that, so
    # I'm going to hardcode the number of parameters.
    if num_params == 0:
        num_params = 7

    # loop through all of them once to figure out which sidebands are seen in all spectra
    for spec in spectrum_list:
        try:
            # use sets to keep track of only unique sidebands
            sb_included = sorted(list(set(sb_included + list(spec.full_dict.keys()))))
        except AttributeError:
            print("No full dict?", spec.fname)
            print(spec.sb_list)
        # If these are from summed spectra, then only the the first file name
        # from that sum will show up here, which should be fine?
        included_spectra[spec.fname.split('/')[-1]] = paramGetter(spec)

    if only_even:
        sb_included = [ii for ii in sb_included if not ii%2]
    if verbose:
        print("included names:", included_spectra)
        print("sb_included:", sb_included)

    for spec in spectrum_list:
        # Flag to keep whethere there are no sidebands or not. Used to skip
        # issues when trying to index on empty arrays
        noSidebands = False
        if verbose:
            print("the sb_results:", spec.sb_results)

        # if no sidebands were found, skip this one
        try:
            # TODO: (08/14/18) the .ndim==1 isn't the correct check, since it fails
            # when looking at the laser line. Need to test this with a real
            # empty data set, vs data set with 1 sb
            #
            #
            # (08/28/18) I'm not sure what the "not spec" is trying to handle
            #      spec.sb_results is None occurs when _no_ sidebands were fit
            #     spec.sb_results.ndim == 1 happens when only one sideband is found
            if not spec or spec.sb_results is None or spec.sb_results.ndim == 1:
                if spec.sb_results is None:
                    # Flag no sidebands are afound
                    noSidebands = True
                elif spec.sb_results[0] == 0:
                    # Cast it to 2d to allow slicing later on. Not sure hwy this is
                    # only done if the laser line is the one found.
                    spec.sb_results = np.atleast_2d(spec.sb_results)
                elif skip_empties:
                    continue
                else:
                    noSidebands = True
        except (AttributeError, TypeError):
            # continue
            raise

        # Make an sb_results of all zeroes where we'll fill
        # in the sideband info we found
        new_spec = np.zeros((len(sb_included), num_params))
        if not noSidebands:
            sb_results = spec.sb_results.copy()
            saw_sbs = sb_results[:, 0]
            found_sb = sorted(list(set(sb_included) & set(saw_sbs)))
            found_idx = [sb_included.index(ii) for ii in found_sb]
            try:
                new_spec[:, 0] = sb_included
            except:
                print("new_spec", new_spec)
                raise
            try:
                if only_even:
                    new_spec[found_idx, :] = sb_results[sb_results[:,0]%2==0]
                else:
                    new_spec[found_idx, :] = sb_results
            except ValueError:
                print(spec.fname)
                print("included:", sb_included)
                print("found:", found_sb, found_idx)
                print(new_spec.shape, sb_results.shape)
                print(sb_results)
                print(new_spec)
                raise

        spec_data = np.insert(new_spec.flatten(), 0, float(paramGetter(spec)))

        try:
            param_array = np.row_stack((param_array, spec_data))
        except:
            param_array = np.array(spec_data)

    if param_array.ndim == 1: # if you only pass one spectra
        param_array = param_array[None, :] # recast it to 2D for slicing
    # the indices we want from the param array from the passed argument
    snip = wanted_indices
    N = len(sb_included)
    # run it out across all of the points across the param_array
    snipped_indices = [0] + list(
        1+np.array(snip * N) + num_params * np.array(sorted(list(range(N)) * len(snip))))
    snipped_array = param_array[:, snipped_indices]
    norm_array = snipped_array.copy()
    # normalize the area if it's requested
    if 3 in snip:
        num_snip = len(snip)
        strength_idx = snip.index(3)
        if 4 in snip:
            #normalize error first if it was requested
            idx = snip.index(4)
            norm_array[:, 1 + idx + np.arange(N) * num_snip] /= norm_array[:,1 + strength_idx + np.arange(N) * num_snip].max(axis=0)
        strength_idx = snip.index(3)
        norm_array[:, 1+strength_idx+np.arange(N)*num_snip]/=norm_array[:, 1+strength_idx+np.arange(N)*num_snip].max(axis=0)

    try:
        os.mkdir(folder_str)
    except TypeError:
        pass # if you pass None as folder_str (for using byteIO)
    except OSError as e:
        if e.errno == errno.EEXIST:
            pass
        else:
            raise

    included_spectra.update(header_dict)
    try:
        included_spectra_str = json.dumps(included_spectra, sort_keys=True, indent=4,
                                          separators=(',', ': '))
    except:
        print("Source: save_parameter_sweep\nJSON FAILED")
        return
    included_spectra_str = included_spectra_str.replace('\n', '\n#')

    included_spectra_str += '\n#' * (99 - included_spectra_str.count('\n'))

    # this will make the header chunk for the full, un-sliced data set
    # TODO: fix naming so you aren't looping twice
    ### 1/9/18 This isn't needed, right? Why isn't it deleted?
    origin_import1 = param_name
    origin_import2 = unit
    origin_import3 = ""
    for order in sb_included:
        origin_import1 += ",sideband,Frequency,error,Sideband strength,error,Linewidth,error"
        origin_import2 += ",order,eV,eV,arb. u.,arb.u.,meV,meV"
        origin_import3 += ",,{0},,{0},,{0},".format(order)
    origin_total = origin_import1 + "\n" + origin_import2 + "\n" + origin_import3


    # This little chunk will make a chunk block of header strings for the sliced
    # data set which can be looped over
    origin_import1 = param_name
    origin_import2 = unit
    origin_import3 = ""
    wanted_titles = ["Sideband", "Frequency", "error", "Sideband strength","error","Linewidth","error"]
    wanted_units  = ["order", "eV", "eV", "arb. u.", "arb. u.", "eV", "eV"]
    wanted_comments = ["", "{0}", "", "{0}", "", "{0}", ""]
    wanted_titles = ",".join([wanted_titles[ii] for ii in wanted_indices])
    wanted_units = ",".join([wanted_units[ii] for ii in wanted_indices])
    wanted_comments = ",".join([wanted_comments[ii] for ii in wanted_indices])

    for order in sb_included:
        origin_import1 += ","+wanted_titles
        origin_import2 += ","+wanted_units
        origin_import3 += ","+wanted_comments.format(order)
    origin_snip = origin_import1 + "\n" + origin_import2 + "\n" + origin_import3

    header_total = '#' + included_spectra_str + '\n' + origin_total
    header_snip = '#' + included_spectra_str + '\n' + origin_snip

    # print "Spec header: ", spec_header
    if verbose:
        print("the param_array is:", param_array)
    if isinstance(file_name, list):
        if isinstance(file_name[0], io.BytesIO):
            np.savetxt(file_name[0], param_array, delimiter=',',
                       header=header_total, comments='', fmt='%0.6e')
            np.savetxt(file_name[1], snipped_array, delimiter=',',
                       header=header_snip, comments='', fmt='%0.6e')
            np.savetxt(file_name[2], norm_array, delimiter=',',
                       header=header_snip, comments='', fmt='%0.6e')
            # Need to reset the file position if you want to read them immediately
            # Is it better to do that here, or assume you'll do it later?
            # I'm gonna assume here, because I can't currently think of a time when I'd want
            # to be at the end of the file
            [ii.seek(0) for ii in file_name]
            if verbose:
                print("Saved the file to bytes objects")
    else:
        if file_name:
            norm_name = file_name + '_norm.txt'
            snip_name = file_name + '_snip.txt'
            file_name = file_name + '.txt'
            np.savetxt(os.path.join(folder_str, file_name), param_array, delimiter=',',
                       header=header_total, comments='', fmt='%0.6e')
            np.savetxt(os.path.join(folder_str, snip_name), snipped_array, delimiter=',',
                       header=header_snip, comments='', fmt='%0.6e')
            np.savetxt(os.path.join(folder_str, norm_name), norm_array, delimiter=',',
                       header=header_snip, comments='', fmt='%0.6e')
            if verbose:
                print("Saved the file.\nDirectory: {}".format(os.path.join(folder_str, file_name)))
        else:
            if verbose:
                print("Didn't save")

    return sb_included, param_array, snipped_array, norm_array

def save_parameter_sweep_vs_sideband(spectrum_list, file_name,
                                     folder_str, param_name, unit, verbose=False,
                                     wanted_indices = [1, 3, 4]):
    """
    Similar to save_parameter_sweep, but the data[:,0] column is sideband number instead of
    series, and each set of columns correspond to a series step. Pretty much compiles
    all of the fit parameters from the files that are already saved and puts it into
    one file to keep from polluting the Origin folder
    :param spectrum_list:
    :param file_name:
    :param folder_str:
    :param param_name:
    :param unit:
    :param verbose:

    sb number is automatically prepended, so do not include in slicing list

    [sb number, Freq (eV), Freq error (eV), Gauss area (arb.), Area error, Gauss linewidth (eV), Linewidth error (eV)]
    [    0    ,      1   ,        2,      ,        3         ,      4    ,         5           ,        6            ]

    :return:
    """
    spectrum_list.sort(key=lambda x: x.parameters[param_name])
    included_spectra = dict()
    param_array = None
    sb_included = []

    # what parameters were included (for headers)
    params = sorted([x.parameters[param_name] for x in spectrum_list])

    for spec in spectrum_list:
        sb_included = sorted(list(set(sb_included + list(spec.full_dict.keys()))))
        included_spectra[spec.fname.split('/')[-1]] = spec.parameters[param_name]
        # If these are from summed spectra, then only the the first file name
        # from that sum will show up here, which should be fine?
    if verbose:
        # print "full name:", spectrum_list[0].fname
        print("included names:", included_spectra)
        print("sb_included:", sb_included)

    param_array = np.array(sb_included)

    for spec in spectrum_list:
        temp_dict = spec.full_dict.copy()

        #prevent breaking if no sidebands in spectrum
        if not temp_dict:
            if verbose:
                print("No sidebands here? {}, {}".format(spec.parameters["series"],
                                                         spec.parameters["spec_step"]))
            continue

        if verbose:
            print(temp_dict)

        # matrix for holding all of the sb information
        # for a given spectrum
        spec_matrix = None
        for sb in sb_included:
            blank = np.zeros(6)
            # print "checking sideband order:", sb
            # print "blank", blank
            sb_data = temp_dict.get(sb, blank)
            try:
                spec_matrix = np.row_stack((spec_matrix, sb_data))
            except:
                spec_matrix = sb_data
        param_array = np.column_stack((param_array, spec_matrix))

    # the indices we want from the param array
    # 1- freq, 3-area, 4-area error
    snip = wanted_indices
    N = len(spectrum_list)
    # run it out across all of the points across the param_array
    snipped_indices = [0] + list( np.array(snip*N) + 6*np.array(sorted(list(range(N))*len(snip))) )
    snipped_array = param_array[:, snipped_indices]

    try:
        os.mkdir(folder_str)
    except OSError as e:
        if e.errno == errno.EEXIST:
            pass
        else:
            raise
    snip_name = file_name + '_snip.txt'
    file_name = file_name + '.txt'

    try:
        included_spectra_str = json.dumps(included_spectra, sort_keys=True, indent=4, separators=(',', ': '))
    except:
        print("Source: save_parameter_sweep\nJSON FAILED")
        return
    included_spectra_str = included_spectra_str.replace('\n', '\n#')

    included_spectra_str += '\n#' * (99 - included_spectra_str.count('\n'))
    origin_import1 = "Sideband"
    origin_import2 = "Order"
    origin_import3 = "SB"
    for param in params:
        origin_import1 += ",Frequency,error,Sideband strength,error,Linewidth,error"
        origin_import2 += ",eV,,arb. u.,,meV,"
        origin_import3 += ",{0},,{0},,{0},".format(param)
    origin_total = origin_import1 + "\n" + origin_import2 + "\n" + origin_import3

    # This little chunk will make a chunk block of header strings for the sliced
    # data set which can be looped over
    origin_import1 = "Sideband"
    origin_import2 = "Order"
    origin_import3 = "SB"
    wanted_titles = ["Sideband", "Frequency", "error", "Sideband strength", "error",
                     "Linewidth", "error"]
    wanted_units = ["order", "eV", "eV", "arb. u.", "arb. u.", "eV", "eV"]
    wanted_comments = ["", "{0}", "", "{0}", "", "{0}", ""]
    wanted_titles = ",".join([wanted_titles[ii] for ii in wanted_indices])
    wanted_units = ",".join([wanted_units[ii] for ii in wanted_indices])
    wanted_comments = ",".join([wanted_comments[ii] for ii in wanted_indices])

    for param in params:
        origin_import1 += "," + wanted_titles
        origin_import2 += "," + wanted_units
        origin_import3 += "," + wanted_comments.format(param)
    origin_snip = origin_import1 + "\n" + origin_import2 + "\n" + origin_import3

    header_total = '#' + included_spectra_str + '\n' + origin_total
    header_snip = '#' + included_spectra_str + '\n' + origin_snip

    # print "Spec header: ", spec_header
    if verbose:
        print("the param_array is:", param_array)
    if file_name: # allow passing false (or empty string) to prevent saving
        np.savetxt(os.path.join(folder_str, file_name), param_array, delimiter=',',
                   header=header_total, comments='', fmt='%0.6e')
        np.savetxt(os.path.join(folder_str, snip_name), snipped_array, delimiter=',',
                   header=header_snip, comments='', fmt='%0.6e')
    if verbose:
        print("Saved the file.\nDirectory: {}".format(os.path.join(folder_str, file_name)))
    return None

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
    if max(first[:, 0]) > max(second[:, 0]):
        flipped = True
        first, second = second, first

    def fitter(p, shiftable, immutable):
        # designed to over

        # Get the shifts
        dx = p[0]
        dy = p[1]

        # Don't want pass-by-reference nonsense, recast our own refs
        shiftable = np.array(shiftable)
        immutable = np.array(immutable)

        # Shift the data set
        shiftable[:, 1] += dy
        shiftable[:, 0] += dx

        # Create an interpolator. We want a
        # direct comparision for subtracting the two functions
        # Different spec grating positions have different wavelengths
        # so they're not directly comparable.
        shiftF = spi.interp1d(*shiftable.T)

        # Find the bounds of where the two data sets overlap
        overlap = (min(shiftable[:, 0]), max(immutable[:, 0]))
        print("overlap", overlap)

        # Determine the indices of the immutable function
        # where it overlaps. argwhere returns 2-d thing,
        # requiring the [0] at the end of each call
        fOlIdx = (min(np.argwhere(immutable[:, 0] >= overlap[0]))[0],
                  max(np.argwhere(immutable[:, 0] <= overlap[1]))[0])
        print("fOlIdx", fOlIdx)

        # Get the interpolated values of the shiftable function at the same
        # x-coordinates as the immutable case
        newShift = shiftF(immutable[fOlIdx[0]:fOlIdx[1], 0])

        if plot:
            plt.plot(*immutable[fOlIdx[0]:fOlIdx[1], :].T, marker='o', label="imm", markersize=10)
            plt.plot(immutable[fOlIdx[0]:fOlIdx[1], 0], newShift, marker='o', label="shift")
        imm = immutable[fOlIdx[0]:fOlIdx[1], 1]
        shift = newShift
        return imm - shift

    a, _, _, msg, err = spo.leastsq(fitter, [0.0001, 0.01 * max(first[:, 1])], args=(second, first), full_output=1)
    # print "a", a
    if plot:
        # Revert back to the original figure, as per top comments
        plt.figure(firstFig.number)

    # Need to invert the shift if we flipped which
    # model we're supposed to move
    if flipped: a *= -1

    return a


def integrateData(data, t1, t2, ave=False):
    """
    Integrate a discrete data set for a
    given time period. Sums the data between
    the given bounds and divides by dt. Optional
    argument to divide by T = t2-t1 for calculating
    averages.

    data = 2D array. data[:,0] = t, data[:,1] = y
    t1 = start of integration
    t2 = end of integration


    if data is a NxM, with M>=3, it will take the
    third column to be the errors of the points,
    and return the error as the quadrature sum
    """
    t = data[:, 0]
    y = data[:, 1]
    if data.shape[0] >= 3:
        errors = data[:, 2]
    else:
        errors = np.ones_like(y) * np.nan

    gt = set(np.where(t > t1)[0])
    lt = set(np.where(t < t2)[0])

    # find the intersection of the sets
    vals = list(gt & lt)

    # Calculate the average
    tot = np.sum(y[vals])
    error = np.sqrt(np.sum(errors[vals] ** 2))

    # Multiply by sampling
    tot *= (t[1] - t[0])
    error *= (t[1] - t[0])

    if ave:
        # Normalize by total width if you want an average
        tot /= (t2 - t1)
        errors /= (t2 - t1)
    if not np.isnan(error):
        return tot, error
    return tot


def fourier_prep(x_vals, y_vals, num=None):
    """
    This function will take a Nx2 array with unevenly spaced x-values and make
    them evenly spaced for use in fft-related things.

    And remove nans!
    """
    y_vals = handle_nans(y_vals)
    spline = spi.interp1d(x_vals, y_vals,
                          kind='linear')  # for some reason kind='quadratic' doesn't work? returns all nans
    if num is None:
        num = len(x_vals)
    even_x = np.linspace(x_vals[0], x_vals[-1], num=num)
    even_y = spline(even_x)
    # even_y = handle_nans(even_y)
    return even_x, even_y


def handle_nans(y_vals):
    """
    This function removes nans and replaces them with linearly interpolated
    values.  It requires that the array maps from equally spaced x-values.
    Taken from Stack Overflow: "Interpolate NaN values in a numpy array"
    """
    nan_idx = np.isnan(y_vals)
    my_lambda = lambda x: x.nonzero()[0]  # Returns the indices where Trues reside
    y_vals[nan_idx] = np.interp(my_lambda(nan_idx), my_lambda(~nan_idx), y_vals[~nan_idx])
    return y_vals


def calc_laser_frequencies(spec, nir_units="eV", thz_units="eV",
                           bad_points=-2, inspect_plots=False):
    """
    Calculate the NIR and FEL frequency for a spectrum
    :param spec: HSGCCD object to fit
    :type spec: HighSidebandCCD
    :param nir_units: str of desired units.
        Options: wavenumber, eV, meV, THz, GHz, nm
    :param thz_units: str of desired units.
        Options: wavenumber, eV, meV, THz, GHz, nm
    :param bad_points: How many bad points which shouldn't be used
        to calculate the frequencies (generally because the last
        few points are noisy and unreliable)
    :return: <NIR freq>, <THz freq>
    """
    if not hasattr(spec, "sb_results"):
        spec.guess_sidebands()
        spec.fit_sidebands()

    sidebands = spec.sb_results[:, 0]
    locations = spec.sb_results[:, 1]
    errors = spec.sb_results[:, 2]
    try:
        p = np.polyfit(sidebands[1:bad_points],
                       # This is 1 because the peak picker function was calling the 10th order the 9th
                       locations[1:bad_points], deg=1)
    except TypeError:
        # if there aren't enough sidebands to fit, give -1
        p = [-1, -1]

    NIRfreq = p[1]
    THzfreq = p[0]

    if inspect_plots:
        plt.figure("Frequency Fit")
        plt.errorbar(sidebands, locations, errors, marker='o')
        plt.errorbar(sidebands[:bad_points], locations[:bad_points],
                     errors[:bad_points], marker='o')
        plt.plot(sidebands, np.polyval(p, sidebands))

    converter = {
        "eV": lambda x: x,
        "meV": lambda x: 1000. * x,
        "wavenumber": lambda x: 8065.6 * x,
        "THz": lambda x: 241.80060 * x,
        "GHz": lambda x: 241.80060 * 1e3 * x,
        "nm": lambda x: 1239.83 / x
    }

    freqNIR = converter.get(nir_units, converter["eV"])(NIRfreq)
    freqTHz = converter.get(thz_units, converter["eV"])(THzfreq)

    return freqNIR, freqTHz

def get_data_and_header(fname, returnOrigin = False):
    """
    Given a file to a raw data file, returns the data
    and the json decoded header.

    Can choose to return the origin header as well
    :param fname: Filename to open
    :return: data, header (dict)
    """
    with open(fname) as fh:
        line = fh.readline()
        header_string = ''
        while line[0]=='#':
            header_string += line[1:]
            line = fh.readline()

        # image files don't have an origin header
        if not "Images" in fname:
            oh = line
            # last readline in loop removes first line in Origin Header
            # strip the remaining two
            oh += fh.readline()
            oh += fh.readline()[:-1] #remove final \n

        # data = np.genfromtxt(fh, delimiter=',')
    data = np.genfromtxt(fname, delimiter=',')

    header = json.loads(header_string)

    if returnOrigin:
        return data, header, oh
    return data, header

def natural_glob(*args):
    # glob/python sort alphabetically, so 1, 10, 11, .., 2, 21,
    # but I sometimes wnat "natural" sorting: 1, 2, 3, ..., 10, 11, 12, ..., 20, 21, 21 ...
    # There's tons of stack overflows, so I grabbed one of them. I put it in here
    # because I use it all the damned time. I also almost always use it when
    # glob.glob'ing, so just internally do it that way
    #
    # This is taken from
    # https://stackoverflow.com/questions/5967500/how-to-correctly-sort-a-string-with-a-number-inside

    import re
    def atoi(text):
        try:
            return int(text)
        except ValueError:
            return text
        # return int(text) if text.isdigit() else text

    def natural_keys(text):
        '''
        alist.sort(key=natural_keys) sorts in human order
        http://nedbatchelder.com/blog/200712/human_sorting.html
        (See Toothy's implementation in the comments)
        '''
        return [atoi(c) for c in re.split('(-?\d+)', text)]

    return sorted(glob.glob(os.path.join(*args)), key=natural_keys)

def convertTime(timeStr):
    """
    The data file headers have the timestamp of data collection. Sometimes you want to
    convert that to numbers for data's sake, but I constantly forget the functions
    to convert it from the time-stamp string. So here you go
    :param timeStr: the time as a string from the data file
    :return: int of the time since the epoch
    """
    import time
    return time.mktime(time.strptime(timeStr, "%x %X%p"))


# photonConverter[A][B](x):
#    convert x from A to B.
photon_converter = {
    "nm":         {"nm": lambda x: x,           "eV": lambda x:1239.84/x,            "wavenumber": lambda x: 10000000./x},
    "eV":         {"nm": lambda x: 1239.84/x,   "eV": lambda x: x,                   "wavenumber":lambda x: 8065.56 * x},
    "wavenumber": {"nm": lambda x: 10000000./x, "eV": lambda x: x/8065.56, "wavenumber": lambda x: x}
}

####################
# Smoothing functions
####################

def savitzky_golay(y, window_size, order, deriv=0, rate=1):
    r"""Smooth (and optionally differentiate) data with a Savitzky-Golay filter.
    The Savitzky-Golay filter removes high frequency noise from data.
    It has the advantage of preserving the original shape and
    features of the signal better than other types of filtering
    approaches, such as moving averages techniques.
    Parameters
    ----------
    y : array_like, shape (N,)
        the values of the time history of the signal.
    window_size : int
        the length of the window. Must be an odd integer number.
    order : int
        the order of the polynomial used in the filtering.
        Must be less then `window_size` - 1.
    deriv: int
        the order of the derivative to compute (default = 0 means only smoothing)
    Returns
    -------
    ys : ndarray, shape (N)
        the smoothed signal (or it's n-th derivative).
    Notes
    -----
    The Savitzky-Golay is a type of low-pass filter, particularly
    suited for smoothing noisy data. The main idea behind this
    approach is to make for each point a least-square fit with a
    polynomial of high order over a odd-sized window centered at
    the point.
    Examples
    --------
    t = np.linspace(-4, 4, 500)
    y = np.exp( -t**2 ) + np.random.normal(0, 0.05, t.shape)
    ysg = savitzky_golay(y, window_size=31, order=4)
    import matplotlib.pyplot as plt
    plt.plot(t, y, label='Noisy signal')
    plt.plot(t, np.exp(-t**2), 'k', lw=1.5, label='Original signal')
    plt.plot(t, ysg, 'r', label='Filtered signal')
    plt.legend()
    plt.show()
    References
    ----------
    .. [1] A. Savitzky, M. J. E. Golay, Smoothing and Differentiation of
       Data by Simplified Least Squares Procedures. Analytical
       Chemistry, 1964, 36 (8), pp 1627-1639.
    .. [2] Numerical Recipes 3rd Edition: The Art of Scientific Computing
       W.H. Press, S.A. Teukolsky, W.T. Vetterling, B.P. Flannery
       Cambridge University Press ISBN-13: 9780521880688

    source:
    http://scipy.github.io/old-wiki/pages/Cookbook/SavitzkyGolay
    """
    import numpy as np
    from math import factorial

    try:
        window_size = np.abs(np.int(window_size))
        order = np.abs(np.int(order))
    except ValueError as msg:
        raise ValueError("window_size and order have to be of type int")
    if window_size % 2 != 1 or window_size < 1:
        raise TypeError("window_size size must be a positive odd number")
    if window_size < order + 2:
        raise TypeError("window_size is too small for the polynomials order")
    order_range = list(range(order + 1))
    half_window = (window_size - 1) // 2
    # precompute coefficients
    b = np.mat([[k ** i for i in order_range] for k in range(-half_window, half_window + 1)])
    m = np.linalg.pinv(b).A[deriv] * rate ** deriv * factorial(deriv)
    # pad the signal at the extremes with
    # values taken from the signal itself
    firstvals = y[0] - np.abs(y[1:half_window + 1][::-1] - y[0])
    lastvals = y[-1] + np.abs(y[-half_window - 1:-1][::-1] - y[-1])
    y = np.concatenate((firstvals, y, lastvals))
    return np.convolve(m[::-1], y, mode='valid')


def fft_filter(data, cutoffFrequency=1520, inspectPlots=False, tryFitting=False, freqSigma=50, ftol=1e-4,
               isInteractive=False):
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
    x = np.array(retData[:, 0])
    y = np.array(retData[:, -1])
    # Let's you place with zero padding.
    zeroPadding = len(x)
    N = len(x)

    if isInteractive:
        try:
            import pyqtgraph as pg
            from PyQt5 import QtCore, QtWidgets
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
                self.plotItem.plot(x, np.log10(y), pen='k')
                # The line for picking the cutoff
                # Connect signals so the textbox updates and the
                # realspace window can recalcualte the FFT
                self.line = pg.InfiniteLine(cutoffFrequency, movable=True)
                self.line.sigPositionChanged.connect(lambda x: self.sigCutoffChanged.emit(x.value()))
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
                self.plotItem.plot(*self.data.T, pen=pg.mkPen('k', width=3))
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

        k = fft.fftfreq(zeroPadding, x[1] - x[0])
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
    onePerc = int(0.01 * N)
    x1 = np.mean(x[:onePerc])
    x2 = np.mean(x[-onePerc:])
    y1 = np.mean(y[:onePerc])
    y2 = np.mean(y[-onePerc:])

    m = (y1 - y2) / (x1 - x2)
    b = y1 - m * x1

    flattenLine = m * x + b
    y -= flattenLine

    if inspectPlots:
        plt.plot(x, y, label="Rotated Data")

    # Perform the FFT and find the appropriate frequency spacing
    k = fft.fftfreq(zeroPadding, x[1] - x[0])
    Y = fft.fft(y, n=zeroPadding)
    if inspectPlots:
        plt.figure("Frequency Space")
        plt.semilogy(k, np.abs(Y), label="Raw FFT")

    if tryFitting:
        try:
            # take +/- 4 sigma points around peak to fit to
            sl = np.abs(k - cutoffFrequency).argmin() + np.array([-1, 1]) * 10 * freqSigma / np.abs(k[0] - k[1])
            sl = slice(*[int(j) for j in sl])
            p0 = [cutoffFrequency,
                  np.abs(Y)[sl].max() * freqSigma,  # estimate the height baased on the max in the set
                  freqSigma,
                  0.14, 2e3, 1.1]  # magic test numbers, they fit the background well

            if inspectPlots:
                plt.semilogy(k[sl], gaussWithBackground(k[sl], *p0), label="Peak with initial values")
            p, _ = curve_fit(gaussWithBackground, k[sl], np.abs(Y)[sl], p0=p0, ftol=ftol)
            if inspectPlots:
                plt.semilogy(k[sl], gaussWithBackground(k[sl], *p), label="Fitted Peak")

            # Want to remove data within 5 sigma ( arb value... )
            st = int(p[0] - 5 * p[2])
            en = int(p[0] + 5 * p[2])

            # Find get the indices to remove.
            refitRangeIdx = np.argwhere((k > st) & (k < en))
            refitRangeIdxNeg = np.argwhere((k < -st) & (k > -en))

            # Replace the data with the backgroudn
            # Note: abuses the symmetry of the FFT of a real function
            # to get the negative side of the data
            Y[refitRangeIdx] = background(k[refitRangeIdx], *p[-2:])
            Y[refitRangeIdxNeg] = background(k[refitRangeIdx], *p[-2:])[::-1]
        except:
            print("ERROR: Trouble fitting the peak in frequency space.\n\t Defaulting to cutting off")

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
        en = int(max(k)) + 1

        # Find the indices to remove the data
        refitRangeIdx = np.argwhere((k > st) & (k < en))
        refitRangeIdxNeg = np.argwhere((k < -st) & (k > -en))

        # Kill it all after the cutoff
        Y[refitRangeIdx] = 0
        Y[refitRangeIdxNeg] = 0

        smoothIdx = np.argwhere((-st < k) & (k < st))
        smoothr = -1. / cutoffFrequency ** 2 * k[smoothIdx] ** 2 + 1

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
        print(x.size, y.size)
        plt.plot(x, y, label="Smoothed Data")
        a = plt.legend()
        a.draggable(True)

    retData[:, 0] = x
    retData[:, -1] = y
    return retData


def low_pass_filter(x_vals, y_vals, cutoff, inspectPlots=True):
    """
    Replicate origin directy
    http://www.originlab.com/doc/Origin-Help/Smooth-Algorithm
    "rotate" the data set so it ends at 0,
    enforcing a periodicity in the data. Otherwise
    oscillatory artifacts result at the ends

    This uses a 50th order Butterworth filter.
    """
    x_vals, y_vals = fourier_prep(x_vals, y_vals)
    if inspectPlots:
        plt.figure("Real Space")
        plt.plot(x_vals, y_vals, label="Non-nan Data")

    zeroPadding = len(x_vals)
    # print "zero padding", zeroPadding  # This needs to be this way because truncation is bad and actually zero padding
    N = len(x_vals)
    onePerc = int(0.01 * N)
    x1 = np.mean(x_vals[:onePerc])
    x2 = np.mean(x_vals[-onePerc:])
    y1 = np.mean(y_vals[:onePerc])
    y2 = np.mean(y_vals[-onePerc:])

    m = (y1 - y2) / (x1 - x2)
    b = y1 - m * x1

    flattenLine = m * x_vals + b
    y_vals -= flattenLine

    if inspectPlots:
        plt.figure("Real Space")
        plt.plot(x_vals, y_vals, label="Rotated Data")

    # even_data = np.column_stack((x_vals, y_vals))
    # Perform the FFT and find the appropriate frequency spacing
    x_fourier = fft.fftfreq(zeroPadding, x_vals[1] - x_vals[0])
    y_fourier = fft.fft(y_vals)  # , n=zeroPadding)

    if inspectPlots:
        plt.figure("Frequency Space")
        plt.semilogy(x_fourier, np.abs(y_fourier), label="Raw FFT")

    # Define where to remove the data
    band_start = cutoff
    band_end = int(max(abs(x_fourier))) + 1

    '''
    # Find the indices to remove the data
    refitRangeIdx = np.argwhere((x_fourier > band_start) & (x_fourier <= band_end))
    refitRangeIdxNeg = np.argwhere((x_fourier < -band_start) & (x_fourier >= -band_end))

    #print "x_fourier", x_fourier[795:804]
    #print "max(x_fourier)", max(x_fourier)
    #print "refitRangeIdxNeg", refitRangeIdxNeg[:-400]

    # Kill it all after the cutoff
    y_fourier[refitRangeIdx] = 0
    y_fourier[refitRangeIdxNeg] = 0

    # This section does a square filter on the remaining code.
    smoothIdx = np.argwhere((-band_start < x_fourier) & (x_fourier < band_start))
    smoothr = -1 / band_start**2 * x_fourier[smoothIdx]**2 + 1

    y_fourier[smoothIdx] *= smoothr
    '''

    # print abs(y_fourier[-10:])
    butterworth = np.sqrt(1 / (1 + (x_fourier / cutoff) ** 100))
    y_fourier *= butterworth

    if inspectPlots:
        plt.plot(x_fourier, np.abs(y_fourier), label="FFT with removed parts")
        a = plt.legend()
        a.draggable(True)
        # print "y_fourier", len(y_fourier)

    # invert the FFT
    y_vals = fft.ifft(y_fourier, n=zeroPadding)

    # using fft, not rfft, so data may have some
    # complex parts. But we can assume they'll be negligible and
    # remove them
    # ( Safer to use np.real, not np.abs? )
    # Need the [:len] to remove zero-padded stuff
    y_vals = y_vals[:len(x_vals)]
    # unshift the data
    y_vals += flattenLine
    y_vals = np.abs(y_vals)

    if inspectPlots:
        plt.figure("Real Space")
        # print x_vals.size, y_vals.size
        plt.plot(x_vals, y_vals, linewidth=3, label="Smoothed Data")
        a = plt.legend()
        a.draggable(True)

    return np.column_stack((x_vals, y_vals))


def high_pass_filter(x_vals, y_vals, cutoff, inspectPlots=True):
    """
    Replicate origin directy
    http://www.originlab.com/doc/Origin-Help/Smooth-Algorithm
    "rotate" the data set so it ends at 0,
    enforcing a periodicity in the data. Otherwise
    oscillatory artifacts result at the ends

    This uses a 50th order Butterworth filter.
    """
    x_vals, y_vals = fourier_prep(x_vals, y_vals)
    if inspectPlots:
        plt.figure("Real Space")
        plt.plot(x_vals, y_vals, label="Non-nan Data")

    zeroPadding = len(x_vals)
    print("zero padding", zeroPadding)  # This needs to be this way because truncation is bad and actually zero padding
    N = len(x_vals)
    onePerc = int(0.01 * N)
    x1 = np.mean(x_vals[:onePerc])
    x2 = np.mean(x_vals[-onePerc:])
    y1 = np.mean(y_vals[:onePerc])
    y2 = np.mean(y_vals[-onePerc:])

    m = (y1 - y2) / (x1 - x2)
    b = y1 - m * x1

    flattenLine = m * x_vals + b
    y_vals -= flattenLine

    if inspectPlots:
        plt.figure("Real Space")
        plt.plot(x_vals, y_vals, label="Rotated Data")

    # even_data = np.column_stack((x_vals, y_vals))
    # Perform the FFT and find the appropriate frequency spacing
    x_fourier = fft.fftfreq(zeroPadding, x_vals[1] - x_vals[0])
    y_fourier = fft.fft(y_vals)  # , n=zeroPadding)

    if inspectPlots:
        plt.figure("Frequency Space")
        plt.semilogy(x_fourier, np.abs(y_fourier), label="Raw FFT")

    # Define where to remove the data
    band_start = cutoff
    band_end = int(max(abs(x_fourier))) + 1

    '''
    # Find the indices to remove the data
    refitRangeIdx = np.argwhere((x_fourier > band_start) & (x_fourier <= band_end))
    refitRangeIdxNeg = np.argwhere((x_fourier < -band_start) & (x_fourier >= -band_end))

    #print "x_fourier", x_fourier[795:804]
    #print "max(x_fourier)", max(x_fourier)
    #print "refitRangeIdxNeg", refitRangeIdxNeg[:-400]

    # Kill it all after the cutoff
    y_fourier[refitRangeIdx] = 0
    y_fourier[refitRangeIdxNeg] = 0

    # This section does a square filter on the remaining code.
    smoothIdx = np.argwhere((-band_start < x_fourier) & (x_fourier < band_start))
    smoothr = -1 / band_start**2 * x_fourier[smoothIdx]**2 + 1

    y_fourier[smoothIdx] *= smoothr
    '''

    print(abs(y_fourier[-10:]))
    butterworth = 1 - np.sqrt(1 / (1 + (x_fourier / cutoff) ** 50))
    y_fourier *= butterworth

    if inspectPlots:
        plt.plot(x_fourier, np.abs(y_fourier), label="FFT with removed parts")
        a = plt.legend()
        a.draggable(True)
        print("y_fourier", len(y_fourier))

    # invert the FFT
    y_vals = fft.ifft(y_fourier, n=zeroPadding)

    # using fft, not rfft, so data may have some
    # complex parts. But we can assume they'll be negligible and
    # remove them
    # ( Safer to use np.real, not np.abs? )
    # Need the [:len] to remove zero-padded stuff
    y_vals = y_vals[:len(x_vals)]
    # unshift the data
    y_vals += flattenLine
    y_vals = np.abs(y_vals)

    if inspectPlots:
        plt.figure("Real Space")
        print(x_vals.size, y_vals.size)
        plt.plot(x_vals, y_vals, label="Smoothed Data")
        a = plt.legend()
        a.draggable(True)

    return np.column_stack((x_vals, y_vals))


def band_pass_filter(x_vals, y_vals, cutoff, inspectPlots=True):
    """
    Replicate origin directy
    http://www.originlab.com/doc/Origin-Help/Smooth-Algorithm
    "rotate" the data set so it ends at 0,
    enforcing a periodicity in the data. Otherwise
    oscillatory artifacts result at the ends

    This uses a 50th order Butterworth filter.
    """
    x_vals, y_vals = fourier_prep(x_vals, y_vals)
    if inspectPlots:
        plt.figure("Real Space")
        plt.plot(x_vals, y_vals, label="Non-nan Data")

    zeroPadding = len(x_vals)
    print("zero padding", zeroPadding)  # This needs to be this way because truncation is bad and actually zero padding
    N = len(x_vals)
    onePerc = int(0.01 * N)
    x1 = np.mean(x_vals[:onePerc])
    x2 = np.mean(x_vals[-onePerc:])
    y1 = np.mean(y_vals[:onePerc])
    y2 = np.mean(y_vals[-onePerc:])

    m = (y1 - y2) / (x1 - x2)
    b = y1 - m * x1

    flattenLine = m * x_vals + b
    y_vals -= flattenLine

    if inspectPlots:
        plt.figure("Real Space")
        plt.plot(x_vals, y_vals, label="Rotated Data")

    # even_data = np.column_stack((x_vals, y_vals))
    # Perform the FFT and find the appropriate frequency spacing
    x_fourier = fft.fftfreq(zeroPadding, x_vals[1] - x_vals[0])
    y_fourier = fft.fft(y_vals)  # , n=zeroPadding)

    if inspectPlots:
        plt.figure("Frequency Space")
        plt.semilogy(x_fourier, np.abs(y_fourier), label="Raw FFT")

    # Define where to remove the data
    band_start = cutoff
    band_end = int(max(abs(x_fourier))) + 1

    '''
    # Find the indices to remove the data
    refitRangeIdx = np.argwhere((x_fourier > band_start) & (x_fourier <= band_end))
    refitRangeIdxNeg = np.argwhere((x_fourier < -band_start) & (x_fourier >= -band_end))

    #print "x_fourier", x_fourier[795:804]
    #print "max(x_fourier)", max(x_fourier)
    #print "refitRangeIdxNeg", refitRangeIdxNeg[:-400]

    # Kill it all after the cutoff
    y_fourier[refitRangeIdx] = 0
    y_fourier[refitRangeIdxNeg] = 0

    # This section does a square filter on the remaining code.
    smoothIdx = np.argwhere((-band_start < x_fourier) & (x_fourier < band_start))
    smoothr = -1 / band_start**2 * x_fourier[smoothIdx]**2 + 1

    y_fourier[smoothIdx] *= smoothr
    '''

    print(abs(y_fourier[-10:]))
    butterworth = 1 - np.sqrt(1 / (1 + (x_fourier / cutoff[0]) ** 50))
    butterworth *= np.sqrt(1 / (1 + (x_fourier / cutoff[1]) ** 50))
    y_fourier *= butterworth

    if inspectPlots:
        plt.plot(x_fourier, np.abs(y_fourier), label="FFT with removed parts")
        a = plt.legend()
        a.draggable(True)
        print("y_fourier", len(y_fourier))

    # invert the FFT
    y_vals = fft.ifft(y_fourier, n=zeroPadding)

    # using fft, not rfft, so data may have some
    # complex parts. But we can assume they'll be negligible and
    # remove them
    # ( Safer to use np.real, not np.abs? )
    # Need the [:len] to remove zero-padded stuff
    y_vals = y_vals[:len(x_vals)]
    # unshift the data
    y_vals += flattenLine
    y_vals = np.abs(y_vals)

    if inspectPlots:
        plt.figure("Real Space")
        print(x_vals.size, y_vals.size)
        plt.plot(x_vals, y_vals, label="Smoothed Data")
        a = plt.legend()
        a.draggable(True)

    return np.column_stack((x_vals, y_vals))


####################
# Complete functions
####################

def proc_n_plotPMT(folder_path, plot=False, confirm_fits=False, save=None, verbose=False, **kwargs):
    """
    This function will take a pmt object, process it completely.

    :rtype: list of HighSidebandPMT
    """
    pmt_data = pmt_sorter(folder_path, plot_individual=plot)

    index = 0
    for spectrum in pmt_data:
        spectrum.integrate_sidebands(verbose=verbose, **kwargs)
        spectrum.laser_line(verbose=verbose, **kwargs)  # This function is broken
        # because process sidebands can't handle the laser line
        # print spectrum.full_dict
        if plot:
            plt.figure('PMT data')
            for sb, elem in list(spectrum.sb_dict.items()):
                plt.errorbar(elem[:, 0], elem[:, 1], elem[:, 2],
                             marker='o', label="{} {}".format(spectrum.parameters["series"],sb))
            plt.figure('Sideband strengths')
            plt.yscale("log")
            plt.errorbar(spectrum.sb_results[:, 1], spectrum.sb_results[:, 3], spectrum.sb_results[:, 4],
                         label=spectrum.parameters['series'], marker='o')
        if plot and confirm_fits:
            plt.figure('PMT confirm fits')
            for elem in list(spectrum.sb_dict.values()):
                plt.errorbar(elem[:, 0], elem[:, 1], elem[:, 2], marker='o')
            plt.errorbar(spectrum.sb_results[:, 1], spectrum.sb_results[:, 3], spectrum.sb_results[:, 4],
                         label=spectrum.parameters['series'], marker='o')
            plt.ylim([-0.005, 0.025])
        if type(save) is tuple:
            spectrum.save_processing(save[0], save[1], index=index)
            index += 1
        elif isinstance(save, str):
            dirr = os.path.dirname(save) if os.path.dirname(save) else '.' # if you just pass a filename tos ave
            spectrum.save_processing(os.path.basename(save), dirr,
                                     index=index)
            index += 1
    if plot:
        plt.legend()
    return pmt_data


def proc_n_plotCCD(folder_path, offset=None, plot=False, confirm_fits=False,
                   save=None, keep_empties = False, verbose=False, **kwargs):
    """
    This function will take a list of ccd files and process it completely.
    save_name is a tuple (file_base, folder_path)

    keep_empties: If True, keep the HighSidebandCCD object in the list if no sidebands
    are found. Else, cut it off.

    The cutoff of 8 is too high, but I don't know what to change it to
    :rtype: list of HighSidebandCCD
    """
    if isinstance(folder_path, list):
        file_list = folder_path
    else:
        # if verbose:
            # print "Looking in:", os.path.join(folder_path, '*seq_spectrum.txt')
        # file_list = glob.glob(os.path.join(folder_path, '*seq_spectrum.txt'))
        file_list = natural_glob(folder_path, '*seq_spectrum.txt')
        # if verbose:
            # print "found these files:", "\n".join([os.path.basename(ii) for ii in file_list])
    raw_list = []
    for fname in file_list:
        raw_list.append(HighSidebandCCD(fname, spectrometer_offset=offset))

    index = 0
    for spectrum in raw_list:
        try:
            spectrum.guess_sidebands(verbose=verbose, plot=plot)
        except RuntimeError:
            print("\n\n\nNo sidebands??\n\n")
            # No sidebands, say it's empty
            if not keep_empties:
                raw_list.pop(raw_list.index(spectrum))
            continue
        try:
            spectrum.fit_sidebands(plot=plot, verbose=verbose)
        except RuntimeError:
            print("\n\n\nNo sidebands??\n\n")
            # No sidebands, say it's empty
            if not keep_empties:
                raw_list.pop(raw_list.index(spectrum))
            continue
        if "calculated NIR freq (cm-1)" not in list(spectrum.parameters.keys()):
            spectrum.infer_frequencies()
        if plot:
            plt.figure('CCD data')
            plt.errorbar(spectrum.proc_data[:, 0], spectrum.proc_data[:, 1], spectrum.proc_data[:, 2],
                         label=spectrum.parameters['series'])
            plt.legend()
            # plt.yscale('log')
            plt.figure('Sideband strengths')
            plt.errorbar(spectrum.sb_results[:, 1], spectrum.sb_results[:, 3], spectrum.sb_results[:, 4],
                         label=spectrum.parameters['series'], marker='o')
            plt.legend()
            plt.yscale('log')
        if plot and confirm_fits:
            plt.figure('CCD confirm fits')
            plt.plot(spectrum.proc_data[:, 0], spectrum.proc_data[:, 1],# spectrum.proc_data[:, 2],
                         label=spectrum.parameters['series'])
            plt.plot(spectrum.sb_results[:, 1], spectrum.sb_results[:, 3] / spectrum.sb_results[:, 5],# spectrum.sb_results[:, 4],
                         label=spectrum.parameters['series'], marker='o')
            plt.legend()
            plt.ylim([-0.1, 1])
        if type(save) is tuple:
            spectrum.save_processing(save[0], save[1],
                                     marker=spectrum.parameters["series"].replace(
                                         r"/", "p"),
                                     index=index)
            index += 1
        elif isinstance(save, str):
            # print "DEBUG: trying to save CCD with ", os.path.dirname(save),'_at_', os.path.basename(save)
            spectrum.save_processing(os.path.basename(save), os.path.dirname(save),
                                     marker=spectrum.parameters["series"].replace(
                                         r"/", "p"),
                                     index=index)
            index += 1
    return raw_list


def create_full_spectra(folder_path, skipLaser = True, *args, **kwargs):
    """
    Given the folder path of raw data (where the PMT data is held in the subfolder "PMT"),
    scale all the data to create a raw comb spectra.
    :param folder_path:
    :param args:
    :param kwargs:
    :return:
    """
    output = np.empty((0,2))
    # have proc_n_plot do all the integrating for the sbs
    pmt = proc_n_plotPMT(os.path.join(folder_path, "PMT"))

    ccd_file_list = glob.glob(os.path.join(folder_path, '*seq_spectrum.txt'))
    ccd_list = [HighSidebandCCD(fname) for fname in ccd_file_list]




    for pmtsb in sorted(pmt[0].sb_dict.keys()):
        if skipLaser and pmtsb == 0: continue
        data = pmt[0].sb_dict[pmtsb]
        try:
            print(pmtsb, pmt[0].full_dict[pmtsb])
        except:
            continue
        output = np.row_stack((output, np.abs(data[:,[0,1]])))
        output = np.row_stack((output, [np.nan, np.nan]))

    # insert the pmt so I can iterate over scaling consecutive pairs
    ccd_list.insert(0, pmt[0])

    # make sure all things get scaled down by the factors before them
    runningRatio = 1
    for idx, ccd in enumerate(ccd_list[1:]):
        ccd.guess_sidebands()
        ccd.fit_sidebands()
        ratio = [1, 1]

        stitch_hsg_dicts(ccd_list[idx], ccd, need_ratio = True, ratios=ratio)

        print("new ratio", ratio)
        runningRatio *= ratio[1]
        ccd.proc_data[:,1]*=runningRatio

        output = np.row_stack((output, np.abs(ccd.proc_data[:,[0,1]])))
        output = np.row_stack((output, [np.nan, np.nan]))

    offsetEnergy = (output[:,0] - pmt[0].full_dict[0][0])*1e3
    print(offsetEnergy.shape, output.shape)
    output = np.column_stack((output[:,0], offsetEnergy.T, output[:,1]))

    return output

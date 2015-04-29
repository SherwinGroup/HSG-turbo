# -*- coding: utf-8 -*-
"""
Created on Tue Apr 28 11:28:51 2015

@author: hbanks

Brevity required, prurience preferred
"""
from __future__ import division
import os
import glob
import numpy as np
import matplotlib.pyplot as plt
from hsganalysis import hsganalysis as hsg

my_path = 'C:/Users/hbanks/Dropbox/Python practice/04-23 Playing with signal and noise'
file_list = glob.glob(os.path.join(my_path, '*_spectrumb.txt'))
print file_list

object_list = []
for fname in file_list:
    object_list.append(hsg.Spectrum(fname))
print "object_list length:", len(object_list)
plt.close('all')
spectra_list = hsg.sum_spectra(object_list)
print "spectra_list length:", len(spectra_list)
for spectrum in spec_list:
    print "About to process file", spectrum.fname
    spectrum.initial_process()
    spectrum.guess_sidebands(cutoff=3.8)
    spectrum.fit_sidebands(plot=True)
    #spectrum.save_processing('300GHz', './Processed spectra2')
    plt.plot(spectrum.hsg_data[:, 0], spectrum.hsg_data[:, 1])

#hsg.save_parameter_sweep(object_list, "test_file", "moment_of_truth", "gain", "pain")
#plt.show()

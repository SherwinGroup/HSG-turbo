# -*- coding: utf-8 -*-
"""
Created on Fri Feb 20 13:46:10 2015

@author: hbanks
"""
import os
import glob
import numpy as np
import matplotlib.pyplot as plt
from hsganalysis import hsganalysis as hsg

my_path = '/Users/dreadnought/Dropbox/OTST stuff/300 GHz low NIR THz sweep'
file_list = glob.glob(os.path.join(my_path, '*'))
print file_list

object_list = []
for fname in file_list:
    object_list.append(hsg.Spectrum(fname))

print object_list

plt.close()
spectra_list = hsg.sum_spectra(object_list)
for spectrum in spectra_list[:2]:
    spectrum.initial_process()
    spectrum.guess_sidebands()
    spectrum.fit_sidebands()
    #spectrum.save_processing('300GHz', './Processed spectra2')
    plt.plot(spectrum.hsg_data[:, 0], spectrum.hsg_data[:, 1])

plt.show()
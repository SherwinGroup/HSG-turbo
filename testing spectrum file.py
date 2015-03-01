# -*- coding: utf-8 -*-
"""
Created on Fri Feb 20 13:46:10 2015

@author: hbanks
"""

import numpy as np
import matplotlib.pyplot as plt
from hsganalysis import hsganalysis as hsg

test1 = hsg.Spectrum('FEL_power_sweep87_spectrum.txt')

print test1.parameters['FEL_pulses']
test1.shot_normalize()
test1.guess_sidebands()

plt.plot(test1.hsg_data[:, 0], test1.hsg_data[:, 1])

test1.fit_sidebands(plot=True)
print "I'll try to save now"
test1.save_processing('test', 'Processed spectra')
print "I got past the saving method"
plt.show()

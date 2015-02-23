# -*- coding: utf-8 -*-
"""
Created on Fri Feb 20 13:46:10 2015

@author: hbanks
"""

import numpy as np
import matplotlib.pyplot as plt
import spectrum as spec

test1 = spec.Base_spectrum('checkingPolarization33_spectrum.txt')
test2 = spec.Base_spectrum('checkingPolarization34_spectrum.txt')
test3 = spec.Base_spectrum('checkingPolarization35_spectrum.txt')
test4 = spec.Base_spectrum('checkingPolarization36_spectrum.txt')
plt.close()

#plt.plot(test1.spectrum[:, 0], test1.spectrum[:, 1])
#plt.plot(test2.spectrum[:, 0], test2.spectrum[:, 1])
#plt.plot(test3.spectrum[:, 0], test3.spectrum[:, 1])
plt.plot(test4.spectrum[:, 0], test4.spectrum[:, 1])
plt.show()
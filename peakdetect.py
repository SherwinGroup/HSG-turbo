# -*- coding: utf-8 -*-
"""
Created on Sat Feb  7 11:11:28 2015

@author: dreadnought

"Brevity required, prurience preferred"
"""

from __future__ import division
import numpy as np
import matplotlib.pyplot as plt

test = np.genfromtxt("fullNIRpower000deg128.txt")

print test
plt.close('all')

plt.plot(test[:,0], test[:,1])



def peak_detect(x_axis, y_axis, window=10):
    '''
    Not sure
    '''
    pre = np.empty(window)
    post = np.empty(window)
    
    pre[:] = -np.Inf
    post[:] = -np.Inf
    y_axis_temp = np.concatenate((pre, y_axis, post))
    
    max_x = []
    max_y = []
    index = 0
    while index < len(x_axis):
        test_value = index + window
        
        check_y = y_axis_temp[test_value - window:test_value + window]
        check_max = check_y.max()
        check_ave = np.mean(abs(check_y[np.isfinite(check_y)]))
        check_value = y_axis_temp[test_value]
        if check_value == check_max and check_value > 4 * check_ave:
            print check_ave
            max_x.append(x_axis[index])
            max_y.append(y_axis[index])
            if check_value > 8 * check_ave:
                index += 20
            else:
                index += 1
        else:
            index += 1
        
    return max_x, max_y

x, y = peak_detect(test[:,0], test[:,1], 80)
plt.plot(x, y, 'o')

plt.show()

newx = [10000000/elem for elem in x]
plt.figure()
plt.plot(range(len(x)), newx, 'o')
plt.show()
print x, y

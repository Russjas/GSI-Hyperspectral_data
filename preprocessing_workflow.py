# -*- coding: utf-8 -*-
"""
Created on Tue Sep 27 09:21:37 2022

@author: russj
"""

import spectral as sp
import spectral.io.envi as envi
import numpy as np
import matplotlib.pyplot as plt
import matplotlib
from scipy.stats import chi2


#%% open data

'''requires either the hdr files to be in the environment variable directory or absolute pathnames for data, white ref and dark ref'''
box_16_2 = envi.open('16-2_55_230m90_235m40_2019-01-10_11-16-33.hdr')
white_ref = envi.open('WHITEREF_16-2_55_230m90_235m40_2019-01-10_11-16-33.hdr')
dark_ref =envi.open('DARKREF_16-2_55_230m90_235m40_2019-01-10_11-16-33.hdr')

#%%'''convert data to np arrays and reduce the number of bands'''
data = np.array(box_16_2.load())
white = np.array(white_ref.load())
dark = np.array(dark_ref.load())
data = data[:, :, 13:262]
white = white[:, :, 13:262]
dark = dark[:, :, 13:262]
bands = box_16_2.bands.centers[13:262]

#%% ''' calib files produces arrays with shape (0, width, bands), but every calib file is the same! and I have dark ref values over 65000, so this file is useless, or I dont understand it'''
calib = envi.open('D:/PYTHON_Messing/Hyperspectral_data/bpr_SWIR_461061.hdr', 'D:/PYTHON_Messing/Hyperspectral_data/bpr_SWIR_461061.bpr')
calibar = np.array(calib.load())

#%% '''instead use errors in the dark ref to identify bad bands'''
l, r, b = np.where(dark > 9000) #9000 is an arbitrary value I came up with by looking at the graphs 
listr = list(r)
listb = list(b)
listl = list(l)

def flatten_bad_band(A, l, r):
    resarr = A.copy()  
    loop_count = 0
    for idx in zip(l, r):
        loop_count += 1
        resarr[:,idx[0], idx[1]] = np.average(np.stack((resarr[:,idx[0], idx[1]], resarr[:,idx[0]+2, idx[1]])), axis=0)
    print(loop_count)    
    return resarr
    print(loop_count)


dark_clean = flatten_bad_band(dark, listr, listb)

white_clean = flatten_bad_band(white,listr, listb)

data_clean = flatten_bad_band(data, listr, listb)


#%% '''reflectance corrected data'''


def reflect_correct(x, w, d):
    wmax = np.mean(w, axis=0)
    dmin = np.mean(d, axis=0)
    result = np.divide(np.subtract(x, dmin), np.subtract(wmax, dmin))
    return result 

data_reflect = reflect_correct(data_clean, white_clean, dark_clean)



#%% ''' save reflectance image '''
envi.save_image(r'D:\PYTHON_Messing\Hyperspectral_data\REFLECTANCE\16_2_55_reflectance.hdr', data_reflect)






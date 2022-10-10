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

#%%
rgb_folder = 'C:/Users/Russell.Rogers/Documents/Hyperspectral/Datafiles/RGB'
server = r'Z:/SWIR/16-2/gsi_10-01-2019/16-2_21_83m00_88m90_2019-01-10_08-23-53/capture'
raw_folder = 'C:/Users/Russell.Rogers/Documents/Hyperspectral/Datafiles/Raw'
reflectance_folder = ['C:/Users/Russell.Rogers/Documents/Hyperspectral/Datafiles/Reflectance/Corrected_terracore', 'C:/Users/Russell.Rogers/Documents/Hyperspectral/Datafiles/Reflectance/Corrected_RR']
useful_folder = 'C:/Users/Russell.Rogers/Documents/Hyperspectral/Utility_files'
filename = '98_2236_05_1_6m00_10m80_2019-05-31_07-39-03.hdr'
white_name = 'WHITEREF_' + filename
dark_name = 'DARKREF_' + filename

'''open data

requires either the hdr files to be in the environment variable directory or absolute pathnames for data, white ref and dark ref'''

box = envi.open('Z:/SWIR/98_2236_05/pl2236/98_2236_05_1_6m00_10m80_2019-05-31_07-39-03/capture/98_2236_05_1_6m00_10m80_2019-05-31_07-39-03.hdr')
white_ref = envi.open('Z:/SWIR/98_2236_05/pl2236/98_2236_05_1_6m00_10m80_2019-05-31_07-39-03/capture/WHITEREF_98_2236_05_1_6m00_10m80_2019-05-31_07-39-03.hdr')
dark_ref = envi.open('Z:/SWIR/98_2236_05/pl2236/98_2236_05_1_6m00_10m80_2019-05-31_07-39-03/capture/DARKREF_98_2236_05_1_6m00_10m80_2019-05-31_07-39-03.hdr')



def local_open(filename):
    box = envi.open(raw_folder + r'/' + filename)
    white_ref = envi.open((raw_folder + r'/' + white_name))
    dark_ref =envi.open((raw_folder + r'/' + dark_name))
    return box, white_ref, dark_ref
def server_open(filename):
    box = envi.open(server + r'/' + filename)
    white_ref = envi.open((server + r'/' + white_name))
    dark_ref =envi.open((server + r'/' + dark_name))
    return box, white_ref, dark_ref

def rgb_open(filename):
    box = envi.open(rgb_folder + r'/' + filename)
    white_ref = envi.open((rgb_folder + r'/' + white_name))
    dark_ref =envi.open((rgb_folder + r'/' + dark_name))
    return box, white_ref, dark_ref

#box, white_ref, dark_ref = server_open(filename)
#box, white_ref, dark_ref = rgb_open(filename)
'''convert data to np arrays and reduce the number of bands'''
data = np.array(box.load())
white = np.array(white_ref.load())
dark = np.array(dark_ref.load())
data = data[:, :, 13:262]
white = white[:, :, 13:262]
dark = dark[:, :, 13:262]
bands = box.bands.centers[13:262]
#np.save('Clipped_Bands.npy', bands) #commented out because file already saved


''' 
calib files produces arrays with shape (0, width, bands), and are supplied by the instrument manufacturer to designate bad pixels.
A filter still needs to be applied to detect inconsistently bad pixels
instead use errors in the dark ref to identify bad bands
This filter is arbritrary, but works for now. Can look into median and gaussian filters
'''

l, r, b = np.where(dark > 9000) # arbritrary value from looking at the graphs
listr = list(r)
listb = list(b)
listl = list(l)

def flatten_bad_band(A, l, r):
    resarr = A.copy()  
    for idx in zip(l, r):       
        #resarr[:,idx[0], idx[1]] = np.average(np.stack((resarr[:,idx[0]+1, idx[1]], resarr[:,idx[0]-1, idx[1]])), axis=0)
        resarr[:,idx[0], idx[1]] = np.average(np.stack((resarr[:,idx[0], idx[1]], resarr[:,idx[0]+2, idx[1]])), axis=0)        
    return resarr 
    
dark_clean = flatten_bad_band(dark, listr, listb)

white_clean = flatten_bad_band(white,listr, listb)

data_clean = flatten_bad_band(data, listr, listb)




'''reflectance corrected data'''


def reflect_correct(x, w, d):
    wmax = np.mean(w, axis=0)
    dmin = np.mean(d, axis=0)
    result = np.divide(np.subtract(x, dmin), np.subtract(wmax, dmin))
    return result * 100

data_reflect = reflect_correct(data_clean, white_clean, dark_clean)



''' save reflectance image '''
envi.save_image((reflectance_folder[1] + r'/ref_' + filename[:7] +'.hdr'), data_reflect)

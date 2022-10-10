

# -*- coding: utf-8 -*-
"""
Created on Wed Sep 28 14:15:21 2022

@author: Russell.Rogers
"""

import spectral as sp
import spectral.io.envi as envi
import numpy as np
import matplotlib.pyplot as plt
import matplotlib
from scipy.stats import chi2
import os

#%%
raw_folder = 'C:/Users/Russell.Rogers/Documents/Hyperspectral/Datafiles/Raw'
reflectance_folder = ['C:/Users/Russell.Rogers/Documents/Hyperspectral/Datafiles/Reflectance/Corrected_terracore', 'C:/Users/Russell.Rogers/Documents/Hyperspectral/Datafiles/Reflectance/Corrected_RR']
useful_folder = 'C:/Users/Russell.Rogers/Documents/Hyperspectral/Utility_files'
###MUST INPUT FILENAME HERE ###
#filename = 'ref_16-2_55_230m90_235m40_2019-01-10_11-16-33.hdr' 
#outputs_folder = reflectance_folder[1] + '/' + filename[:-4]
#if not os.path.isdir(outputs_folder):
#    os.makedirs(outputs_folder)




'''open data

requires either the hdr files to be in the environment variable directory or absolute pathnames for data, white ref and dark ref'''
box = envi.open('C:/Users/Russell.Rogers/Documents/Hyperspectral/Datafiles/Reflectance/Corrected_RR/ref_16-2_22/denoisedLength_ref_16-2_22.hdr')
data = np.array(box.load())
bands = np.load('C:/Users/Russell.Rogers/Documents/Hyperspectral/Utility_files/Clipped_Bands.npy')

#%%
def normalize_cube(data):
    newarr = np.zeros_like(data)
    for i in range(data.shape[0]):
        for j in range(data.shape[1]):
            newarr[i,j] = data[i, j] / np.max(data[i, j])
    return newarr

norm = normalize_cube(data)

Cr_tight = sp.remove_continuum(norm[:, :, 220:], bands[220:])
Cr_16_2_22 = sp.remove_continuum(norm, bands)

tightminima = np.argmin(Cr_tight, axis=2)




#%%

envi.save_image('C:/Users/Russell.Rogers/Documents/Hyperspectral/Datafiles/Reflectance/Corrected_RR/ref_16-2_22/normalised.hdr' , norm)
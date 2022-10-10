# -*- coding: utf-8 -*-
"""
Created on Thu Oct  6 12:38:49 2022

@author: Russell.Rogers

Transcription of the SQM method for determining wavelenghth position of the minima of the spectral absorbtion feature. Transcribed from:
   Rodger, Andrew & Laukamp, Carsten & Haest, Maarten & Cudahy, Thomas. (2012). 
   A simple quadratic method of absorption feature wavelength estimation in continuum removed spectra. 
   Remote Sensing of Environment. 118. 273–283. 10.1016/j.rse.2011.11.025. 
input data must be continuum removed and recommend that bands are sliced to 
"""

import spectral as sp
import spectral.io.envi as envi
import numpy as np
import matplotlib.pyplot as plt


#%%

class SQM:
    def __init__(self, wavelengths, depths):
        self.wavelengths = wavelengths
        self.depths = depths
       

def get_SQM(data, bands):
        m, n = data.shape[:2]
        tru_arr = np.zeros((m,n))
        dep_arr = np.zeros((m,n))
        for i in range(m):
            for j in range(n):
                b = np.argmin(data[i, j])                
                D0 = data[i, j, b]
                D1 = data[i, j, (b+1)]
                Dx1 = data[i, j, (b-1)]
                W0 = bands[b]
                W1 = bands[b+1]
                Wx1 = bands[b-1]
                An = ((D0 - Dx1)  * (Wx1- W1)) + ((D1 - Dx1) * (W0 - Wx1))
                Ad = ((Wx1 - W1) * ((W0**2) - (Wx1**2))) + ((W0 - Wx1) * ((W1**2) - (Wx1**2)))
                A = An / Ad
                B = ((D0-Dx1) - (A *(W0**2 - Wx1**2))) / (W0 - Wx1)
                tru_arr[i, j] = (-1 * B) / (2 * A)
                C = Dx1 - (A * (Wx1**2)) - B
                dep_arr[i, j] = 1 - ((A * tru_arr[i, j]**2) + (B * (tru_arr[i, j])) + C)        
        return SQM(tru_arr, dep_arr)
#%%
if __name__ == "__main__":

    box = envi.open('C:/Users/Russell.Rogers/Documents/Hyperspectral/Datafiles/Reflectance/Corrected_RR/ref_16-2_22/Clip_CRtight.hdr')
    data = np.array(box.load())
    bandsa = np.load('C:/Users/Russell.Rogers/Documents/Hyperspectral/Utility_files/Clipped_Bands.npy')
    bands = bandsa[220:]    
    
    sqm = get_SQM(data, bands)
    wav = sqm.wavelengths
    dep = sqm.depths
#%%    
    bins = np.arange(np.min(wav), np.max(wav), ((np.max(wav) - np.min(wav)) / 10))
    bins2 = np.arange(np.min(dep), np.max(dep), ((np.max(wav) - np.min(wav)) / 10))
    x = np.digitize(wav, bins)
    x = np.swapaxes(x, 0, 1)
    y = np.digitize(dep, bins2)
#%%    
    rgb = envi.open('C:/Users/Russell.Rogers/Documents/Hyperspectral/Datafiles/Reflectance/Corrected_RR/ref_16-2_22/denoisedClip_ref_16-2_22.hdr')
    rgbData = np.array(rgb.load())
    Swav = np.swapaxes(wav, 0, 1)
    Sdep = np.swapaxes(dep, 0, 1)
    Srgb = np.swapaxes(rgbData, 0, 1)
    
    #sp.imshow(x, cmap='RdBu_r', aspect='auto', fignum=1)
    sp.imshow(Swav, cmap='RdBu_r')
    sp.imshow(Sdep, cmap='RdBu_r')
    sp.imshow(Srgb, classes=x)
#%%
    plt.figure()
    plt.subplot(211)
    plt.plot(bands, data[392, 16], bands, data[423, 16])
    plt.title('Carbonate absorbtion band in convex hull quotient')    
    plt.subplot(212)
    plt.plot(bandsa, Srgb[16, 423], bandsa, Srgb[16, 392])
    plt.title('Absorbtion bands in reflectance data')
    plt.suptitle('Pixes either side of contact')
    plt.tight_layout()
           



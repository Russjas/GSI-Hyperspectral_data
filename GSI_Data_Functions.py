# -*- coding: utf-8 -*-
"""
Created on Thu Oct 13 09:36:04 2022

@author: Russell.Rogers
"""

import spectral as sp
import spectral.io.envi as envi
import numpy as np
import scipy.signal as sig

''' Pre=processing'''
def data_open(filename, path):
    white_name = 'WHITEREF_' + filename
    dark_name = 'DARKREF_' + filename
    box = envi.open(path + r'/' + filename)
    white_ref = envi.open((path + r'/' + white_name))
    dark_ref =envi.open((path + r'/' + dark_name))
    return box, white_ref, dark_ref

def flatten_bad_band(A, dark):
    l, r, b = np.where(dark > 9000) # arbritrary value from looking at the graphs
    print(len(l)), print(len(r)), print(len(b))
    listr = list(r)
    listl = list(l)
    listb = list(b)
    resarr = A.copy()  
    for idx in zip(listr, listb):       
        resarr[:,idx[0], idx[1]] = np.average(np.stack((resarr[:,idx[0]+1, idx[1]], resarr[:,idx[0]-1, idx[1]])), axis=0)
        #resarr[:,idx[0], idx[1]] = np.average(np.stack((resarr[:,idx[0], idx[1]], resarr[:,idx[0]-1, idx[1]])), axis=0)        
    return resarr 

def reflect_correct(x, w, d):
    wmax = np.mean(w, axis=0)
    dmin = np.mean(d, axis=0)
    result = np.divide(np.subtract(x, dmin), np.subtract(wmax, dmin))
    return result * 100

''' useful functions for a variety of things'''
def normalize(data):
    newarr = np.zeros_like(data)
    for i in range(data.shape[0]):
        for j in range(data.shape[1]):
            newarr[i,j] = data[i, j] / np.max(data[i, j])
    return newarr

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
    
''' re-descretize the library spectra to have the same bands as the data'''
def get_matching_spectrum(bands, lib):
    x = np.array(lib.x) * 1000
    y = np.array(lib.y)
    match_spectra = []
    for i in range(len(bands)):
        a = np.interp(bands[i], x, y)
        match_spectra.append(a)
    return np.array(match_spectra)

def get_peaks_cube(data):
    M, N, B = data.shape
    A = np.zeros((M, N, 4))
   
    for i in range(M):
        for j in range(N):
            x, props = sig.find_peaks(-data[i,j], prominence= 0)
            y = np.sort(props['prominences'])[-5]
            z = []
            for b in range(len(x)):
                if props['prominences'][b] > y:
                    z.append(x[b])
            if len(z) == 3:
                z.append(-999)
            elif len(z) == 2:
                z.append([-999, -999])
            elif len(z) == 1:
                z.append([-999, -999, -999])
            elif len(z) == 0:
                z.append([-999, -999, -999, -999])
            A[i, j] = np.array(z)
    return A
def get_deepest_peak(data, bands, peak):
    M, N, B = data.shape
    A = np.zeros((M, N, 1))
   
    for i in range(M):
        for j in range(N):
            x, props = sig.find_peaks(-data[i,j], height= -0.99)
            #print(props)
            y = np.sort(props['peak_heights'])[-peak]
            for b in range(len(x)):
                if props['peak_heights'][b] == y:
                    A[i, j] = np.array(bands[x[b]])            
    return A

def class_by_peak(data):
    '''data must be a M x N x 1 dataset where the one is the wavelength band of
    deepest absorption feature'''
    M, N, B = data.shape
    wavelengths = np.sort(np.unique(data))
    a = data.copy()
    for i in range(M):
        for j in range(N):
            for x in range(len(wavelengths)):
                if data[i, j] == wavelengths[x]:
                    a[i, j] = x
    return a
                
            
# -*- coding: utf-8 -*-
"""
Created on Thu Sep 29 08:38:50 2022

@author: russj
"""

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


Box = envi.open(r'D:\PYTHON_Messing\Hyperspectral_data\REFLECTANCE\16_2_55_reflectance2.hdr')
data = np.array(Box.load())
view = sp.imshow(data)

#    ''''be sure to look at the image produced and note the pixels you want to clip out in the nect cell. Bands were clipped in pre-processing''''
#%%

clip = data[29:1283, 127:256]

#%%
pc_clip = sp.principal_components(clip)

plt.figure()
plt.plot(pc_clip.eigenvalues,linestyle='--', marker='o', color='b')
plt.title('"eignevalues in clipped data')

denoised1 = pc_clip.denoise(clip, num=1)

veiw = sp.imshow(denoised1)

#%%
import scipy.signal as sc


filtered1 = sc.savgol_filter(clip, 5, 1)

filtered3 = sc.savgol_filter(clip, 10, 1)
filtered4 = sc.savgol_filter(clip, 15, 1)
filtered5 = sc.savgol_filter(clip, 5, 2)

filtered7 = sc.savgol_filter(clip, 10, 2)
filtered8 = sc.savgol_filter(clip, 15, 2)
filtered9 = sc.savgol_filter(clip, 5, 3)

filtered11 = sc.savgol_filter(clip, 10, 3)
filtered12 = sc.savgol_filter(clip, 15, 3)

#%%
view = sp.imshow(filtered12)
#%%
plt.figure()
plt.subplot(321)
plt.plot(clip[701, 60])
plt.title('Original reflectance')
plt.subplot(322)
plt.plot(denoised1[701, 60])
plt.title('PCA')
plt.subplot(323)
plt.plot(filtered1[701, 60])
plt.title('savgol filt 1')
plt.subplot(324)
plt.plot(filtered3[701, 60])
plt.title('savgol filt 2')
plt.subplot(325)
plt.plot(filtered4[701, 60])
plt.title('savgol filt 3')
plt.tight_layout()


#%%
'''kmeans clustering'''



(m, c) = sp.kmeans(filtered1, 10, 100)


view = sp.imshow(filtered1, classes=m, figsize=(4,9))



#%% ''' load database '''
db = sp.EcostressDatabase('D:/PYTHON_Messing/Hyperspectral_data/fine_carbonates.db')

db.print_query('SELECT COUNT() FROM Samples WHERE Name LIKE "%dolomite%"')
db.print_query('SELECT SampleID, Name FROM Samples WHERE Name LIKE "%dolomite%"')
db.print_query('SELECT COUNT() FROM Samples WHERE Name LIKE "%calcite%"')
db.print_query('SELECT SampleID, Name FROM Samples WHERE Name LIKE "%calcite%"')
#%%

bands = np.load('D:/PYTHON_Messing/Hyperspectral_data/Clipped_Bands.npy')

def get_matching_spectrum(bands, lib):
    x = np.array(lib.x) * 1000
    y = np.array(lib.y)
    match_spectra = []
    for i in range(len(bands)):
        a = np.interp(bands[i], x, y)
        match_spectra.append(a)
    return np.array(match_spectra)

dolo = db.get_signature(42)
mdolo = get_matching_spectrum(bands, dolo)
calc = db.get_signature(35)
mcalc = get_matching_spectrum(bands, calc)
endmembers = np.vstack((mdolo, mcalc))

#%%



test = sp.unmix(filtered1, endmembers)
test2 = np.argmin(test, axis=2)
test3 = np.argmax(test, axis=2)


plt.plot(test[50, :, 0])
view = sp.imshow(filtered1, classes=test2, figsize=(4,9))
view = sp.imshow(filtered1, classes=test3, figsize=(4,9))
#%%

filtered_ppi = sp.ppi(filtered1,20 )




view = sp.imshow(filtered1, classes=filtered_ppi)


























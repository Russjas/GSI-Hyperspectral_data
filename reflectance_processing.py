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


Box_55 = envi.open(r'C:/Users/Russell.Rogers/Documents/Hyperspectral/reflectance/16_2_55_reflectance_meanSD.hdr')
data = np.array(Box_55.load())
clip= data[48:1292, 135:256]
view = sp.imshow(clip)

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
'''sav gol filter goes here, but I am not re-writing the code when I have it at home'''



(m, c) = sp.kmeans(filtered1, 10, 100)


view = sp.imshow(filtered1, classes=m, figsize=(4,9))


#%%
for i in range (c.shape[0]):
    plt.plot(c[i])
plt.grid()

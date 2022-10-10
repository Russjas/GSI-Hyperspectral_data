

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
filename = 'ref_RBG_16-2_22.hdr' 
outputs_folder = reflectance_folder[1] + '/' + filename[:-4]
if not os.path.isdir(outputs_folder):
    os.makedirs(outputs_folder)
'''open data

requires either the hdr files to be in the environment variable directory or absolute pathnames for data, white ref and dark ref'''
box = envi.open(reflectance_folder[1] + r'/' + filename)
data = np.array(box.load())
view = sp.imshow(data)

'''be sure to look at the image produced and note the pixels you want to clip out in the nect cell. Bands were clipped in pre-processing'''
#%%

clip = data[333:13653, 1536:2909] #full box

#%%
run1above = data[333:13653, 2515:2909] # rightmost run

run1below = data[224:1283, 232:266] # rightmost run
run1 = np.concatenate((run1above, run1below))
run2 = data[333:13653, 2025:2419] #middle run

run3above = data[333:13653, 1536:1930] #leftmost run
run3below = data[358:1300, 144:178] #leftmost run
run3 = np.concatenate((run3above, run3below))
length = np.concatenate((run1above, run2, run3above))


#%%
pc_clip = sp.principal_components(clip)
#%%
pc_run1 = sp.principal_components(run1)
pc_run2 = sp.principal_components(run2)
pc_run3 = sp.principal_components(run3)
pc_length = sp.principal_components(length)

plt.figure()
plt.plot(pc_clip.eigenvalues,linestyle='--', marker='o', color='b')
plt.title('"eignevalues in clipped data')

denoised_clip = pc_clip.denoise(clip, num=1)
denoised_run1 = pc_clip.denoise(run1, num=1)
denoised_run2 = pc_clip.denoise(run2, num=1)
denoised_run3 = pc_clip.denoise(run3, num=1)
denoised_length = pc_clip.denoise(length, num=1)



#veiw = sp.imshow(denoised1)

envi.save_image((outputs_folder + '/denoisedClip_' + filename) , denoised_clip)
envi.save_image((outputs_folder + '/denoisedRun1_' + filename) , denoised_run1)
envi.save_image((outputs_folder + '/denoisedRun2_' + filename) , denoised_run2)
envi.save_image((outputs_folder + '/denoisedRun3_' + filename) , denoised_run3)                
envi.save_image((outputs_folder + '/denoisedLength_' + filename) , denoised_length)
                

                
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

#envi.save_image(r'D:\PYTHON_Messing\Hyperspectral_data\REFLECTANCE\16_2_55_reflectance2_filtered.hdr', filtered1)




























#!/usr/bin/env python3
# -*- coding: utf-8 -*-
'''
   Notes
    -----
    -   BCARS Preprocessing Code provided by Ryan Muddiman and Timothy McNamara
    -   Based on the work of  Camp Jr. et al. [1]

    References
    ----------
    [1] C. H. Camp Jr, Y. J. Lee, and M. T. Cicerone, "Quantitative, \
    Comparable Coherent Anti-Stokes Raman Scattering (CARS) \
    Spectroscopy: Correcting Errors in Phase Retrieval," Journal of Raman \
    Spectroscopy 47, 408-415 (2016). arXiv:1507.06543.
'''

#%%
import os
file_path = "piezo flow 9.asc"
python_import_path = os.path.expanduser("FILEPATH")
os.chdir(python_import_path)


#%% Set working directory to source anscombe transform and other tools
os.chdir(python_import_path)

#%% Imports and load dataset
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import time as time
from numpy.linalg import svd as svd
from anscombe_transform import *
import h5py
from scipy.linalg import diagsvd as diagsvd
from scipy.signal import savgol_filter
from crikit.cri.kk import KramersKronig
from crikit.cri.error_correction import PhaseErrCorrectALS as PEC
from crikit.cri.error_correction import ScaleErrCorrectSG as SEC
import matplotlib as mpl
from brokenaxes import brokenaxes
import matplotlib.gridspec as gridspec
from matplotlib.colors import LogNorm, Normalize

import math
mpl.rcParams.update(mpl.rcParamsDefault)
%matplotlib widget

# Load data
raw_data = pd.read_csv(file_path, delimiter="\t")
raw_data = raw_data.to_numpy()[:,:]
raw_data=raw_data.T
raw_data = np.concatenate([raw_data[0:7812,:], raw_data[7833:,:]])
data_shape = raw_data.shape
print(data_shape)

# Store wavelength axis for later
wl_loaded = raw_data[0,:]
# First we crop the spectral length as sometimes we have laser present on the blue side of the spectrum
right_cutoff = 1010
#  Import pre-calibrated wavelength spectrum (this needs updating when plotting for a publication)
wl = np.squeeze(np.asarray(pd.read_csv(r'wl 1.csv',header=None)))
wl=wl[10:right_cutoff]
wn = 1e7*(1/wl-1/770.75)

# Remove wavelength axis row
raw_data = raw_data[1:-1,:]
fig = plt.figure()
plt.title("First spectra")
plt.plot(wl_loaded,raw_data[0])

df = raw_data

#%% Plot a single spectrum
figure=plt.figure()
plt.plot(df[999,:])
plt.show()


#%% Plot a section of spectra
figure=plt.figure()
for i in range(0,np.size(df,axis=0)):
    plt.plot(wl_loaded,df[i,:])
plt.show()

#%%  Remove air bubble
df  = np.concatenate([raw_data[0:7810,:], raw_data[8356:,:,]],axis=0)
fig = plt.figure()
plt.title("First spectra")
plt.plot(wl_loaded,raw_data[0])

#%% Ploting unprocessed spectra in colour mapS
row_min = np.min(np.min(df[:,:]))
row_max = np.max(np.max(df[:,:]))
figure=plt.figure()
plt.imshow(df[:,:].T, aspect='auto', cmap='Spectral', interpolation=None)
plt.show()


#%% Determining locations of plastic spectra
for i, spec in enumerate (df):
        if spec[854]-spec[775] > 420:
            print(f"{i}" )

figure=plt.figure()
plt.plot(df[175,:])
plt.show()

#%%
fig = plt.figure(figsize=(10,6),constrained_layout=False)
spec = gridspec.GridSpec(ncols=1, nrows=1, figure=fig)

bax1 = brokenaxes(xlims=((500, 1800),(2700, 3300)), subplot_spec=spec[0], wspace = 0.1)
print(np.min(wn))
bax1.plot(wn,df[610,:1000],label = 'glass',linewidth=1)

bax1.set_xlabel('Raman Shift ($\mathrm{cm}^{-1}$)',labelpad = 20)
bax1.set_ylabel('Intensity (A.U.)',labelpad = 40)
plt.show()

#%%
###
###         UPDATE SIGNAL VALUE & BACKGROUND RANGE
###
signal_i = 175
bg_range = (1000,1025)

#%% Singular Value Decomposition - Economy SVD
# We calculate the image mean and standard deviation in a background region - required for SVD reconstruction
mean = np.mean(np.mean(df[bg_range[0]:bg_range[1],450:500],1),0)

#%%
std = np.mean(np.std(df[bg_range[0]:bg_range[1],450:500],1),0)
new_im = df[bg_range[0]: bg_range[1], 0]
print(new_im) 
print(df[0])
# df_ans = Anscombe transformed data for removing heteroscedasticity
df_ans = np.zeros((np.size(df,0), np.size(df,1)))

for i in range(np.size(df,0)):
    ansc = generalized_anscombe(df[i,:], mu=mean, sigma=std)      
    df_ans[i,:] = ansc

#%%
# Calculate singular value decomposition of array
U,S,V = svd(df_ans, full_matrices=False)
N = np.size(df,0)
sigma= std
cutoff = (4/np.sqrt(3))*np.sqrt(N)*sigma
#r = np.sum((S>cutoff))

#% Display highest singular values eigen-spectra
fig, ax= plt.subplots(nrows = 5, ncols = 5,figsize = (10,8))
ax = ax.flatten()
for i in range(25):
    r = i
    denoised = V[r,:]  
    label = r'$SV_{{{}}}$'.format((str(i+1)))
    ax[i].plot(denoised,label=label)
    ax[i].set_title(label)
fig.tight_layout()
plt.show()
plt.plot(denoised)
plt.show()

#%% Show eigen-maps of singular vectors
fig, ax= plt.subplots(nrows = 4, ncols = 5,figsize = (11,11))
ax = ax.flatten()
for i in range(20):
    r = i
    U1 = U[:,r]
    denoised = U1
    label = r'$SV_{{{}}}$'.format((str(i+1)))
    #z1= ax[i].imshow(denoised,label=label,cmap = 'hsv')
    ax[i].set_title(label)
    ax[i].set_aspect('equal')
    ax[i].set_xticks([])
    ax[i].set_yticks([])
plt.show()

#%% Calculate compact SVD reconstruction based on custom singular values
# These are the singular values which will be used for low-rank reconstruction (must use 0-20 always)
ucr = [1,2,3,4,5,6,7,8,9,10,11,12,13,114,15,16,17,18,19,20]
r = []
for val in ucr:
    r.append(val-1)
r=range(20)

#%%
spec_len = np.size(df,1)
side_len = np.size(df,0)

s_select = np.zeros(S.size)
s_select[r] = S[r]

M = U.shape[-1]
N = V.shape[0]
Sr = diagsvd(s_select,M,N)
denoised = U @ (Sr) @ V
df_ans = np.reshape(denoised,(side_len,spec_len))
df_denoised = np.zeros((side_len,spec_len))

# Calculating the inverse Anscombe transform on the SVD reconstructed data
for i in range(side_len):
        inv_ans = inverse_generalized_anscombe(df_ans[i,:],mu=mean, sigma=std)     
        df_denoised[i,:] = inv_ans

#%% This is for checking the quality of denoising - set x and y to a signal_i pixel based off the image
fig = plt.figure()
plt.plot(df[signal_i,:],label = 'before SVD',c = 'k',linewidth = 1)
plt.plot(df_denoised[signal_i,:]+100,label = 'After SVD',c = 'r',linewidth = 1)
plt.legend()
plt.show()


#%% This is for optimizing the Kramers-Kronig NRB removal, we first perform it on a single spectrum and evaluate results
# Setting NRB as an average of a 'background' pixel area
nrb = np.mean(df[bg_range[0]:bg_range[1]],0)

fig = plt.figure()
plt.plot(nrb)
plt.show()


#%%
# Cropping signals
signal = df_denoised[:,10:right_cutoff]
nrb = nrb[10:right_cutoff]

nrb_min = nrb
nrb_signal = nrb_min
spec_len = np.size(signal,1)
new_signal = np.ones((side_len,spec_len))

# can be used for pixel normalisation (usually ignore)
for i in range(side_len):
    new_signal[i,:] = signal[i,:]#-signal[i,k,:].min()
    new_signal[i,:] = new_signal[i,:]#/new_signal[i,k,:].max()
        
cars_spectrum = new_signal[signal_i,:]       


fix_rng = np.hstack((np.arange(4,6))) #370-600
fix_rng=fix_rng.astype(int)

rng = np.arange(0,np.size(cars_spectrum))
#idx = np.s_[0:0]
#rng[idx]=0
rng=rng.astype(int)

asym_param = 1e-3*np.ones(spec_len)
asym_param[380:800] = 1e-3
asym_param[0:380] = 1e-3

kk = KramersKronig(conjugate=1, norm_to_nrb=True,pad_factor=1)

retrieved = kk.calculate(cars_spectrum, nrb_signal)

pec = PEC(wavenumber_increasing=False, smoothness_param=1e3, asym_param=asym_param, redux=1, 
          fix_end_points=False, rng=rng, fix_rng=fix_rng, fix_const=0)

sec = SEC()
phase_corrected = pec.calculate(retrieved)

result = savgol_filter(sec.calculate(phase_corrected).imag,7,5)

fig, ax= plt.subplots(nrows = 2, ncols = 1,figsize = (9,8))
ax = ax.flatten()

ax[0].plot(np.angle(retrieved),label = 'raw phase')
#plt.plot(np.angle(pec.calculate(phase_corrected)),label = 'phase error')
ax[0].plot(np.angle(retrieved)-np.angle(pec.calculate(phase_corrected)),label = ' error')
ax[0].legend()
ax[1].plot(wn,result,label = 'Retrieved Im[$\chi^{(3)}$]')
ax[1].legend()
plt.show()

#%% Kramers-Kronig Retrieval on whole dataset using above single-pixel parameters
KK_signal = kk.calculate(new_signal, nrb_signal)    
phase_corrected = pec.calculate(KK_signal)
result = savgol_filter(sec.calculate(phase_corrected).imag,11,5)

#%% Plot Signal and NRB example
fig = plt.figure(figsize=(10,6),constrained_layout=False)
spec = gridspec.GridSpec(ncols=1, nrows=1, figure=fig)
wn = 1e7*(1/wl-1/765)
bax1 = brokenaxes(xlims=((500, 1800),(2700, 3300)), subplot_spec=spec[0], wspace = 0.1)
bax1.plot(wn,result[signal_i,:],label = 'signal',c='#1874CD',linewidth=1)
bax1.plot(wn,result[bg_range[0],:],label = 'water',c='tab:green',linewidth=1)

#bax1.set_ylim(-0.0,0.1)
bax1.set_xlabel('Raman Shift ($\mathrm{cm}^{-1}$)',labelpad = 20)
bax1.set_ylabel('Intensity (A.U.)',labelpad = 40)
#bax1.legend(loc='upper left')
bax1.set_title('Retrieved BCARS spectra')

plt.show()
bax1.legend(loc=2)
plt.savefig('retrieved.png',dpi=600)

#%% This is for choosing the wavelengths to color
fig = plt.figure()
plt.plot(result[signal_i,:],label = 'signal',linewidth = 1)
plt.plot(result[bg_range[0],:],label = 'glass',linewidth = 1)
plt.legend()
plt.show()


#%%
plastic = []
for i, spec in enumerate(result):
    if result[i,840] > 0.01:
        plastic.append(result[i,:])
        print(i)
print(len(plastic))


#%%
fig = plt.figure(figsize=(10,6),constrained_layout=False)
spec = gridspec.GridSpec(ncols=1, nrows=1, figure=fig)
bax1 = brokenaxes(xlims=((500, 1800),(2700, 3300)), subplot_spec=spec[0], wspace = 0.1)

print(len(plastic))
for i in plastic:
    bax1.plot(wn, i)
bax1.set_xlabel('Raman Shift ($\mathrm{cm}^{-1}$)',labelpad = 20)
bax1.set_ylabel('Intensity (A.U.)',labelpad = 40)
bax1.set_title(f'Retrieved BCARS spectra')
bax1.legend(loc=2)
plt.show()



#%% Plot of Spectrum
index = 11
fig = plt.figure(figsize=(10,6),constrained_layout=False)
spec = gridspec.GridSpec(ncols=1, nrows=1, figure=fig)
bax1 = brokenaxes(xlims=((500, 1800),(2700, 3300)), subplot_spec=spec[0], wspace = 0.1)

result_normalised = np.zeros((np.size(result,0),np.size(result,1)))
X_min = np.min(result,axis=1)
X_max = np.max(result,axis=1)
# Perform max-min normalization
for i in range(0, np.size(result,0)):
    result_normalised[i,:] = (result[i,:] - X_min[i]) / (X_max[i] - X_min[i])

bax1.plot(wn,result_normalised[index,:],label = 'signal',c='#1874CD',linewidth=1.5)
for ax in bax1.axs:
    ax.yaxis.set_major_formatter('')

bax1.set_xlabel('Raman Shift ($\mathrm{cm}^{-1}$)',labelpad = 20)
bax1.set_ylabel('Intensity (A.U.)',labelpad = 10 )
bax1.set_title(f' BCARS Spectrum of PMMA')
plt.show()


#%% Plot of Spectrum
index =  175
fig = plt.figure(figsize=(10,6),constrained_layout=False)
spec = gridspec.GridSpec(ncols=1, nrows=1, figure=fig)
bax1 = brokenaxes(xlims=((500, 1800),(2700, 3300)), subplot_spec=spec[0], wspace = 0.1)

result_normalised = np.zeros((np.size(result,0),np.size(result,1)))
X_min = np.min(result,axis=1)
X_max = np.max(result,axis=1)
# Perform max-min normalization
for i in range(0, np.size(result,0)):
    result_normalised[i,:] = (result[i,:] - X_min[i]) / (X_max[i] - X_min[i])

bax1.plot(wn,result_normalised[index,:],label = 'signal',c='#1874CD',linewidth=1.5)
for ax in bax1.axs:
    ax.yaxis.set_major_formatter('')

bax1.set_xlabel('Raman Shift ($\mathrm{cm}^{-1}$)',labelpad = 20)
bax1.set_ylabel('Intensity (A.U.)',labelpad = 10 )
bax1.set_title(f' BCARS Spectrum of PS')
plt.show()

#%% Filtered Plastic Samplesand plot

fig=plt.figure()
plt.title("Filtered Plastic Samples")
plastic_index = []
for i, spec in enumerate (result):
        if spec[840]> 0.015:
            print(f"{i}" )
            plastic_index.append(i)
            plt.plot(result[i,:])
plt.show()


#%% Plot a range of spectra
fig = plt.figure(figsize=(10,6),constrained_layout=False)
spec = gridspec.GridSpec(ncols=1, nrows=1, figure=fig)

bax1 = brokenaxes(xlims=((500, 1800),(2700, 3300)), subplot_spec=spec[0], wspace = 0.1)
for i in range (9450,9499):
    bax1.plot(wn, result[i,:])    
    
bax1.set_xlabel('Raman Shift ($\mathrm{cm}^{-1}$)',labelpad = 20)
bax1.set_ylabel('Intensity (A.U.)',labelpad = 30)
bax1.set_title(f'Retrieved BCARS spectra')
bax1.legend(loc=2)
plt.show()


#%% Plot a spectrum
fig = plt.figure(figsize=(10,6),constrained_layout=False)
spec = gridspec.GridSpec(ncols=1, nrows=1, figure=fig)
bax1 = brokenaxes(xlims=((500, 1800),(2700, 3300)), subplot_spec=spec[0], wspace = 0.1)
bax1.set_xlabel('Raman Shift ($\mathrm{cm}^{-1}$)',labelpad = 20)
bax1.set_ylabel('Intensity (A.U.)',labelpad = 30)
#bax1.legend(loc='upper left')
bax1.set_title(f'Retrieved BCARS spectra')
bax1.plot(wn, result[signal_i,:])
plt.show()


#%% Export Result
export_file = h5py.File('denoised_plastic_datacube2.h5', 'w')
export_file.create_dataset("Spectra", data=result)
export_file.close()


#%% Plot of Complete Spectrum
norm_result = np.zeros((result.shape[0], result.shape[1]))

for i in range(result.shape[0]):
   
    for j in range (result.shape[1]):
        if i in plastic_index:

            min_val = np.min(np.min(result[i,:]))
            max_val = np.max(np.max(result[i,:]))
            norm_result[i,j] = (result[i,j] - min_val)/(max_val - min_val )
        else: 
            norm_result[i,j] = result[i,j]

fig, ax = plt.subplots()
wn_ir = [210,1000]
raster = ax.imshow(norm_result[:,wn_ir[0]:wn_ir[1]].T,cmap="Spectral", interpolation=None) #Remove norm=LogNorm() to go back to normal scale
cbar = fig.colorbar(raster, ax=ax)
cbar.set_label('Spectral Intensity', rotation=270, labelpad=15)  
ax.set_ylabel('Raman Shift ($\mathrm{cm}^{-1}$)',labelpad = 10)
ax.set_xlabel('Time (sec)',labelpad = 10 )
ax.set_aspect('auto')
# Define the number of wavenumbers to show on the y-axis
div_yaxis = 50 
wn_new = wn[wn_ir[0]:wn_ir[1]]
wn_new2 = (wn_new[::div_yaxis])
y_positions = np.arange(0, len(wn_new), div_yaxis) 
ax.set_yticks(y_positions)
ax.set_yticklabels([int(wn) for wn in wn_new2])  

# Define the acquisition time and range for x-axis
acquisition_time = 0.002
div_xaxis = 1000
xaxis = np.arange(result.shape[0]) 
xaxis_new = acquisition_time * np.arange(0, result.shape[0], div_xaxis)
xaxis_new_positions = np.arange(0, len(xaxis), div_xaxis) 
ax.set_xticks(xaxis_new_positions)
ax.set_xticklabels(xaxis_new)
plt.show()


#%% Plot of Fingerprint Region
norm_result = np.zeros((result.shape[0], result.shape[1]))
for i in range(result.shape[0]):
   
    for j in range (result.shape[1]):
        if i in plastic_index:

            min_val = np.min(np.min(result[i,:]))
            max_val = np.max(np.max(result[i,:]))
            norm_result[i,j] = (result[i,j] - min_val)/(max_val - min_val )
        else: 
            norm_result[i,j] = result[i,j]


wn_ir = [589,1000]
fig, ax = plt.subplots()

raster = ax.imshow(norm_result[:,wn_ir[0]:wn_ir[1]].T,cmap="Spectral", interpolation=None) 
cbar = fig.colorbar(raster, ax=ax)
cbar.set_label('Spectral Intensity', rotation=270, labelpad=15)  
ax.set_ylabel('Raman Shift ($\mathrm{cm}^{-1}$)',labelpad = 10)
ax.set_xlabel('Time (sec)',labelpad = 10 )
ax.set_aspect(10)

# Define the number of wavenumbers to show on the y-axis
div_yaxis = 50 
wn_new = wn[wn_ir[0]:wn_ir[1]]
wn_new2 = (wn_new[::div_yaxis])
y_positions = np.arange(0, len(wn_new), div_yaxis) 
ax.set_yticks(y_positions)
ax.set_yticklabels([int(wn) for wn in wn_new2])  

# Define the acquisition time and range for x-axis
acquisition_time = 0.002
div_xaxis = 1000
xaxis = np.arange(result.shape[0]) 
xaxis_new = acquisition_time * np.arange(0, result.shape[0], div_xaxis)
xaxis_new_positions = np.arange(0, len(xaxis), div_xaxis) 
ax.set_xticks(xaxis_new_positions)
ax.set_xticklabels(xaxis_new)
plt.show()


#%% Plot of CH region
norm_result = np.zeros((result.shape[0], result.shape[1]))
for i in range(result.shape[0]):
   
    for j in range (result.shape[1]):
        if i in plastic_index:

            min_val = np.min(np.min(result[i,:]))
            max_val = np.max(np.max(result[i,:]))
            norm_result[i,j] = (result[i,j] - min_val)/(max_val - min_val )
        else: 
            norm_result[i,j] = result[i,j]

wn_ir = [200,400]
fig, ax = plt.subplots()
raster = ax.imshow(result[:,wn_ir[0]:wn_ir[1]].T,cmap="Spectral", interpolation=None, norm=LogNorm()) #Remove norm=LogNorm() to go back to normal scale
cbar = fig.colorbar(raster, ax=ax)
cbar.set_label('Spectral Intensity', rotation=270, labelpad=15)
ax.set_ylabel('Raman Shift ($\mathrm{cm}^{-1}$)',labelpad = 10)
ax.set_xlabel('Time (sec)',labelpad = 10 ) 
ax.set_aspect(50)

# Define the number of wavenumbers to show on the y-axis
div_yaxis = 50 
wn_new = wn[wn_ir[0]:wn_ir[1]]
wn_new2 = (wn_new[::div_yaxis])
y_positions = np.arange(0, len(wn_new), div_yaxis) 
ax.set_yticks(y_positions)
ax.set_yticklabels([int(wn) for wn in wn_new2])  

# Define the acquisition time and range for x-axis ticks
acquisition_time = 0.002
div_xaxis = 1000
xaxis = np.arange(result.shape[0]) 
xaxis_new = acquisition_time * np.arange(0, result.shape[0], div_xaxis)
xaxis_new_positions = np.arange(0, len(xaxis), div_xaxis) 
ax.set_xticks(xaxis_new_positions)
ax.set_xticklabels(xaxis_new)
plt.show()


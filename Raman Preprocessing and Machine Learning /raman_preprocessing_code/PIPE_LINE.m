%% Load your dataset - should be single species with 4-ace recorded at the start
clear all;
close all;
load ../ML_Data/ML_Data_60.mat %should contain 'spectra_raw' (matrix eg 1024 x 59) and 'Four_Ace' (1024 x 1) - can optionally contain 'glass' %Acquisition time of glasss should be same as for spectra
figure;plot(spectra_raw);hold;plot(glass);
figure;plot(Four_Ace);
%% INTENSITY CALIBRATION
load int_cal.mat %should contain 'neon' (1024 x 1), 'White_Light' (1024 x 1), 'White_Light_Background' (1024 x 1)
WAVELENGTH_CAL_NEON %Make sure the peak positons and neon wavelengths are correct - see first two lines in m file
INTENSITY_CAL %This will produce a correction factor C, that you will then multiply all of your spectra by.
spectra_intensity_calibrated = spectra_calibrated;

%% WAVENUMBER CALIBRATION
WAVENUMBER_CALIBRATE %This requries Four_Ace to be availble and will generate wavenumber_axis which is relevant for the spectrum recorded soon after
p=348; %p should be the pixel position of the peak you want to use to calibrate all the spectra
spectra_calibrated = Peak_Calibrate(spectra_calibrated,p,'n',0); 

%% COSMIC RAY REMOVAL
REMOVE_COSMIC_RAYS
% This will produce a new dataset called spectra_cosmic_rays_removed

%% DENOISE
spectra_denoised = sgolayfilt(spectra_cosmic_ray_removed,3,11);

%% BASELINE REMOVAL
BASELINE_REMOVAL

%% PLOT FIGURES
PLOT_FIGURES

%% SAVE DATA
save ../PVC_6_processed.mat spectra_baseline_removed wavenumber_axis

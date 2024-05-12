%% EE499 - Individual Processing of Spectra in a mixed dataset
clear all
close all
clc
addpath Pipeline

load('ML_Data/ML_Data_60.mat');
figure;plot(spectra_raw);hold;plot(glass);
figure;plot(Four_Ace);

%% Step 1 - Instensity Calibration

load int_cal.mat %should contain 'neon' (1024 x 1), 'White_Light' (1024 x 1), 'White_Light_Background' (1024 x 1)
WAVELENGTH_CAL_NEON %Make sure the peak positons and neon wavelengths are correct - see first two lines in m file
INTENSITY_CAL %This will produce a correction factor C, that you will then multiply all of your spectra by.
spectra_intensity_calibrated = spectra_calibrated;


%% Step 2 -  Wavenumber Calibration
% Generate wavenumber axis
WAVENUMBER_CALIBRATE 

%% Step 3 - Denoise
spectra_denoised = sgolayfilt(spectra_calibrated,3,11);

%% Step 4 - Glass Removal
% Parameter
p=0.9;
N=3; 
iterations = 200; 

n = size(spectra_denoised,2);
spectra_baseline_removed = [];

for i = 1:n
    [corrected_spectrum, backgrounds] = berger_background_reduction(spectra_denoised(:,i), glass , p, N, iterations, 'y','n');
    spectra_baseline_removed = [spectra_baseline_removed corrected_spectrum];
end

%% Step 4 - Min-Max Normalisation

n = size(spectra_baseline_removed,2);
spectra_normalised = [];
% Normalising the dataset
for i = 1:n
    min_value = min(spectra_baseline_removed(:,i));
    max_value = max(spectra_baseline_removed(:,i));
    normalised_spectrum  = (spectra_baseline_removed(:,i) - min_value)/(max_value-min_value);
    spectra_normalised = [spectra_normalised normalised_spectrum];
end

figure;
plot(wavenumber_axis, spectra_normalised);
xlabel('Wavenumber (cm^{-1})');
ylabel('Intensity (AU)');
set(gca,'XLim',[200 1800]);



%% Creating ML Dataset with Labels
labels = [zeros([20,1]); ones([20,1]); 2*ones([20,1])];
dataset = spectra_normalised;
save('ML_Data_60_Processed.mat', 'dataset','wavenumber_axis', 'labels');
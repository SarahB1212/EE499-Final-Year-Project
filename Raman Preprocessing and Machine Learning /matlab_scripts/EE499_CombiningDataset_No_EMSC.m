%% EE499 - Combining Datasets and Normalisation
clear all
close all
clc

addpath Pipeline

% Defining Plastic Type
plastic_type = 'PMMA';

filename = [plastic_type '_matfiles/' plastic_type '_'];
number_of_datasets = 5;
combined_data = [];
combined_wn_axis = [];

%% Combining all of the spectra into one dataset 
for i = 1:number_of_datasets
    data = load([filename int2str(i) '_processed.mat']);
    dataset = data.spectra_baseline_removed;
    wn_axis = data.wavenumber_axis;
    combined_data = [combined_data dataset];
    combined_wn_axis = [combined_wn_axis ; wn_axis];
end

[m,n] = size(combined_data);

for i = 1:n
    % Normalising the combined dataset
    min_value = min(combined_data(:,i));
    max_value = max(combined_data(:,i));
    combined_data(:,i)  = (combined_data(:,i) - min_value)/(max_value-min_value);
end
wavenumber_axis = mean(combined_wn_axis);

figure;
plot(combined_data);

%% Peak Calibrate combined data
%p is pixel position of the peak you want to use to calibrate all the spectra
p = 278;
combined_data = Peak_Calibrate(combined_data,p,'n',0); 

figure;
plot(wavenumber_axis, combined_data);
title([plastic_type ' Raman Spectra Dataset']);
xlabel('Wavenumber (cm^{-1})');
ylabel('Intensity (AU)');
set(gca,'XLim',[200 1800]);
set(gca,'YLim',[min_value inf]);
set(gca,'YTick',[]);

file_name = [plastic_type '_Calibrated_Dataset.mat'];
dataset = combined_data;
save(file_name, 'wavenumber_axis','dataset');


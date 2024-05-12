%% EE499 - Formating Data to .mat file of a single plastic type dataset
clear all

% Define Paths
wn_path = "ASCII/Wavenumber_Calibration/4_ace_28_3_24.asc";
neon_path = "ASCII/Wavenumber_Calibration/Neon.asc";
glass_path = 'BackgroundSubtractionData';
folder_name = 'SortedData/PMMA';
mat_filename = 'PMMA_1';
dataset_number = 1;

fileName = [mat_filename '.mat'];

% Get wavenumber spectrum e.g. 4-ace
wn_data = load(wn_path);
Four_Ace = wn_data(:,2);
save(fileName, 'Four_Ace');

% Get neon spectrum
neon_data = load(neon_path);
neon_spectrum = neon_data(:,2);
save(fileName, 'neon_spectrum','-append');

% Get glass spectrum
filename = [glass_path '/background' int2str(dataset_number) '_1.txt'];
file_data = load(filename);
glass = file_data(:,2);
save(fileName, 'glass','-append');

% Get sample spectra
data = dir([folder_name '/sample' int2str(dataset_number) '*.txt']);
n = numel(data);
spectra_raw = [];
for i = 1:n
    filename = [folder_name '/' data(i).name];
    file_data = load(filename);
    spectra_raw = [spectra_raw file_data(:,2)];
end
save(fileName, 'spectra_raw', '-append');


%% EE499 - Formating Data to .mat file for mixed plastic type dataset
clear all

% Define Paths
wn_path = "ASCII/Wavenumber_Calibration/4_ace_28_3_24.asc";
neon_path = "ASCII/Wavenumber_Calibration/Neon.asc";
glass_path = 'BackgroundSubtractionData';
folder_name = 'RawData';
mat_filename = 'ML_Data_60';
dataset_number = 17;

fileName = [mat_filename '.mat'];

% Get wavenumber spectrum e.g. 4-ace
wn_data = load(wn_path);
Four_Ace = wn_data(:,2);
save(fileName, 'Four_Ace');

% Get neon spectrum
neon_data = load(neon_path);
neon_spectrum = neon_data(:,2);
save(fileName, 'neon_spectrum','-append');

% Get glass spectrum
filename = [glass_path '/background' int2str(dataset_number) '_1.txt'];
file_data = load(filename);
glass = file_data(:,2);
save(fileName, 'glass','-append');

% Get sample spectra
data = dir(['SortedData/PMMA/sample' int2str(dataset_number) '*.txt']);
n = numel(data);
spectra_raw = [];
figure;
for i = 1:n
    filename = [folder_name '/' data(i).name];
    file_data = load(filename);
    plot(file_data(:,2));
    hold on;
    spectra_raw = [spectra_raw file_data(:,2)];
end

figure;
data = dir(['SortedData/PVC/sample' int2str(dataset_number) '*.txt']);
n = numel(data);
for i = 1:n
    filename = [folder_name '/' data(i).name];
    file_data = load(filename);
        plot(file_data(:,2));
    hold on;
    spectra_raw = [spectra_raw file_data(:,2)];
end

figure
data = dir(['SortedData/PS/sample' int2str(dataset_number) '*.txt']);
n = numel(data);
for i = 1:n
    filename = [folder_name '/' data(i).name];
    file_data = load(filename);
        plot(file_data(:,2));
    hold on;
    spectra_raw = [spectra_raw file_data(:,2)];
end


save(fileName, 'spectra_raw', '-append');
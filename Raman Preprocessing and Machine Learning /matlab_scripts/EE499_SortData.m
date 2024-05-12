%% EE499 - Sort Data m-file

folder_name = '21_11_23';

%% Plot all spectra to sort each Raman spectrum by plastic type
folder = dir(['RawData/sample17_*']);
n = numel(folder);
file_names = {folder.name};
background_folder = dir('BackgroundSubtractionData/Background');
background_names = {background_folder.name};

xlim([0 2500])
for i =1:n
    name = file_names{i};
    dataset_number = name(7);
    background_data = load(['BackgroundSubtractionData/Background' dataset_number '_1.txt']);
    data = load(['RawData/' name]);
    fprintf(['RawData/' name '\n']);
    raman_spectum = data(:,2) - background_data(:,2);
    plot(data(:,1),raman_spectum);
    pause;
end

%% Plot all spectra
file_path = 'SortedData/PMMA/';
folder = dir([file_path 'sample*']);
n = numel(folder);
file_names = {folder.name};
background_folder = dir('BackgroundSubtractionData/Background');
background_names = {background_folder.name};

xlim([0 2500])
for i =1:n
    name = file_names{i};
    dataset_number = name(7);
    background_data = load(['BackgroundSubtractionData/Background' dataset_number '_1.txt']);
    data = load([file_path name]);
    fprintf([file_path name '\n']);
    raman_spectum = data(:,2) - background_data(:,2);
    plot(data(:,1),raman_spectum);
    hold on;
end 

file_path = 'SortedData/PS/';
folder = dir([file_path 'sample*']);
n = numel(folder);
file_names = {folder.name};
background_folder = dir('BackgroundSubtractionData/Background');
background_names = {background_folder.name};

xlim([0 2500])
for i =1:n
    name = file_names{i};
    dataset_number = name(7);
    background_data = load(['BackgroundSubtractionData/Background' dataset_number '_1.txt']);
    data = load([file_path name]);
    fprintf([file_path name '\n']);
    raman_spectum = data(:,2) - background_data(:,2);
    plot(data(:,1),raman_spectum);
    hold on;
end 
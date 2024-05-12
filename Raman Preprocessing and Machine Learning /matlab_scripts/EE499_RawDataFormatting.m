%% EE499 - Formating Data to .txt file
clear all
close all
clc

% Name of folder within the ASCII folder containing the .asc files
% In each folder the .asc files are named 
%   - background1.asc, background2.asc etc.
%   - sample1.asc, sample2.asc etc.
folder_name = '28_3_24';

% Assigning a dataset number to each recorded dataset
dataset_number = '17';

% Extracting the background spectra
data = dir(['ASCII/' folder_name '/background*.asc']);
n = numel(data);
for i = 1:n
    filename = ['ASCII/' folder_name '/background' int2str(i) '.asc'];
    file_data = load(filename);
    pixel_axis = file_data(:,1);
    Raman_data = file_data(:,2);
    formatted_data = [pixel_axis file_data];

    % Saving .txt file to the folder 'BackgroundSubtractionData'
    string = ['BackgroundSubtractionData/background' dataset_number '_' int2str(i) '.txt'];
    writematrix(formatted_data,string,'delimiter',' ');
end

% Extracting the recorded Raman spectra
data = dir(['ASCII/' folder_name '/sample*.asc']);
n = numel(data);
for i = 1:n
    filename = ['ASCII/' folder_name '/sample' int2str(i) '.asc'];
    file_data = load(filename);
    pixel_axis = file_data(:,1);
    Raman_data = file_data(:,2);
    formatted_data = [pixel_axis file_data];
    
    % Saving .txt file to the folder 'BackgroundSubtractionData'
    string = ['RawData/sample' dataset_number '_' int2str(i) '.txt'];
    writematrix(formatted_data,string,'delimiter',' ');
end

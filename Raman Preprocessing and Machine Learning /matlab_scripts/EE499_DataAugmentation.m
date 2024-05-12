%% EE499 - Data Aumentation
clear all 
clc

plastic_type = 'PS';
filename = [plastic_type '_matfiles/' plastic_type '_Calibrated_Dataset.mat'];
data = load(filename);

wavenumber_axis = data.wavenumber_axis;
spectrum_data = data.dataset;
[m, n] = size(spectrum_data);

augmented_dataset = [];
augmented_spectrum = zeros(m,1);

%% (1) Shifting each spectrum left or right randomly
for i = 1:n
    wn_shift = randi([1 5],1);
    shift = randi([0 1]);
    if shift == 0
        augmented_spectrum = [spectrum_data(1:wn_shift,i); spectrum_data(1:end-wn_shift,i)];
    else
        augmented_spectrum = [spectrum_data(1:end-wn_shift,i); spectrum_data(end-wn_shift+1:end,i)];
    end
   
    augmented_dataset = [augmented_dataset augmented_spectrum];
end

plot(augmented_dataset);
save([plastic_type '_Augmented_Dataset_1.mat'],'augmented_dataset','wavenumber_axis');

%% (2) Added a random noise 

for i = 1:n
    augmented_spectrum = awgn(spectrum_data(:,i),40);
    augmented_dataset = [augmented_dataset augmented_spectrum];
end

plot(augmented_dataset);

save([plastic_type '_Augmented_Dataset_2.mat'],'augmented_dataset','wavenumber_axis');


%% (3) Linear combinations of all spectra belonging to the same substance as augmented data. 
% Linear combination of 5 spectra
% The coefficients in the linear combination were chosen at random.
for j = 1:n
    spectra_index = randi([1 n],[1 5]);
    coefficient = randi(10,5,1);
    for i=1:5
        index = spectra_index(1,i);
        augmented_spectrum = augmented_spectrum+(coefficient(i).*spectrum_data(:,index));
    end
    augmented_spectrum = (augmented_spectrum - min(augmented_spectrum))/(max(augmented_spectrum)-min(augmented_spectrum));
    augmented_dataset = [augmented_dataset augmented_spectrum];
    augmented_spectrum = zeros(m,1);
end
plot(augmented_dataset)
save([plastic_type '_Augmented_Dataset_3.mat'],'augmented_dataset','wavenumber_axis');
%% PLEASE REFERENCE THE FOLLOWING PAPER FOR COSMIC RAY REMOVAL
% Barton, Sinead J., and Bryan M. Hennelly. "An algorithm for the removal of cosmic ray artifacts in spectral data sets." Applied spectroscopy 73.8 (2019): 893-901.

ref = mean(spectra_calibrated,2);
N=7;

if exist('glass', 'var') == 1
    [s, back, C_R, C_G] = EMSC_Matrix_Glass(spectra_calibrated,ref',glass',N);
else
    [s, back, C_R] = EMSC_Matrix_No_Glass(spectra_calibrated,ref,N);
end

temp = CRARreboot(s');

spectra_cosmic_ray_removed = temp.*(C_R') + back';
spectra_cosmic_ray_removed = spectra_cosmic_ray_removed';

clear ref N s back C_R C_N temp

if exist('glass', 'var') == 1
    clear C_G
end

%figure;plot(spectra);
%figure;plot(spectra_cosmic_ray_removed);



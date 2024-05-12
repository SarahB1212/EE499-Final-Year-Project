%% REMOVING BASELINE HAS FOUR OPTIONS
% Select the one you want - be sure to referene appropriately
option = 4;

%% 1 HAVE glass - HAVE REFERENCE FROM CAF2
% PLEASE REFERENCE THE FOLLOWING PAPERS
% Kerr, L. T., and B. M. Hennelly. "A multivariate statistical investigation of background subtraction algorithms for Raman spectra of cytology samples recorded on glass slides." Chemometrics and Intelligent Laboratory Systems 158 (2016): 61-68.
% Skogholt, Joakim, Kristian Hovde Liland, and Ulf Geir Indahl. "Preprocessing of spectral data in the extended multiplicative signal correction framework using multiple reference spectra." Journal of Raman Spectroscopy 50.3 (2019): 407-417.

if option == 1
load caf2 %should  wavenumber corrrect
ref = caf2; %Simply use the mean of the dataset as the REF and remove baseline/normalise against this REF
N=7; % select your polynomial order - preferably low
[spectra_baseline_removed, back, C_R, C_G] = EMSC_Matrix_Glass(spectra_denoised,ref,glass', N);
clear ref back C_R C_G caf2 N
end

%% 2 glass NOT AN ISSUE - SPECTRUM NOT ON A HILL - SIMPLIFIED EMSC
% PLEASE REFERENCE THE FOLLOWING PAPERS
% Kerr, L. T., and B. M. Hennelly. "A multivariate statistical investigation of background subtraction algorithms for Raman spectra of cytology samples recorded on glass slides." Chemometrics and Intelligent Laboratory Systems 158 (2016): 61-68.
% Skogholt, Joakim, Kristian Hovde Liland, and Ulf Geir Indahl. "Preprocessing of spectral data in the extended multiplicative signal correction framework using multiple reference spectra." Journal of Raman Spectroscopy 50.3 (2019): 407-417.

if option == 2
ref = mean(spectra_denoised,2); %Simply use the mean of the dataset as the REF and remove baseline/normalise against this REF
N=3; % select your polynomial order - preferably low
[spectra_baseline_removed, back, C_R] = EMSC_Matrix_No_Glass(spectra_denoised,ref,N);
clear ref back C_R N
end

%% 3 glass NOT AN ISSUE - SPECTRUM ON A HILL - RUBBER BAND METHOD - FOLLOWED BY EMSC
% PLEASE REFERENCE THE FOLLOWING PAPERS
% FOR RUBBER BAND
% S. Wartewig, IR and Raman Spectroscopy: Fundamental Processing, Spectroscopic Techniques, Wiley-VCH: Weinheim, 2003.
% FOR EMSC
% Kerr, L. T., and B. M. Hennelly. "A multivariate statistical investigation of background subtraction algorithms for Raman spectra of cytology samples recorded on glass slides." Chemometrics and Intelligent Laboratory Systems 158 (2016): 61-68.

if option == 3
N=5;% select your polynomial order - preferably low
iterations = 50;
ref = spectra_denoised(:,1);
%figure;plot(ref);
[ref,p] = rubber_band(ref, N, iterations);
%figure;plot(ref);
[spectra_baseline_removed, back, C_R] = EMSC_Matrix_No_Glass(spectra_denoised,ref,N);
clear ref back C_R N iterations p
end

%% 4 HAVE glass - HAVE NO REF - SPECTRUM ON A HILL WITH GLASS - BERGER METHOD - FOLLOWED BY EMSC
if option == 4
p=0.9; %guess at concentration
N=3; iterations = 200; 
%figure;hold;plot(spectra(:,1));plot(glass);hold off;
[ref, backgrounds] = berger_background_reduction(spectra_denoised(:,1), glass, p, N, iterations, 'y','y');
%figure; plot(glass);
%instead of using berger you might prefer to manually test different
%weights og glass being subtracted from spectra(:,1) and taking the one you
%like the most.
%figure;plot(ref);
N=3;
[spectra_baseline_removed, back, C_R, C_G] = EMSC_Matrix_Glass(spectra_denoised,ref,glass',N);
clear ref back C_R C_G N p iterations backgrounds
end

figure;plot(spectra_denoised);
figure;plot(spectra_baseline_removed);
clear option

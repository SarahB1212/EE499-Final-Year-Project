% PLEASE RERFERENCE THE FOLLOWING PAPER:
% Hutsebaut, Didier, Peter Vandenabeele, and Luc Moens. "Evaluation of an accurate calibration and spectral standardization procedure for Raman spectroscopy." Analyst 130.8 (2005): 1204-1214.

%figure; plot(White_Light); hold; plot(White_Light_Background); hold off;
X = White_Light - White_Light_Background;
%figure;plot(X);hold;
%figure; plot(X);
X=sgolayfilt(X,3,31);
X=X/max(X);
%plot(X);hold off;
white_ref=dlmread('7003P2027_HL-3 plus-INT-CAL_int_20190614_VIS(e).lmp','\t');%525-650
temp_a=polyfit(white_ref(1:21,1),white_ref(1:21,2),7);
TRUE_NIST=polyval(temp_a,wavelength_axis.*1e9); 
TURE_NIST=TRUE_NIST/max(TRUE_NIST);
%figure;plot(TRUE_NIST);hold;plot(X);
C=TRUE_NIST./X;
C=C/min(C);
C(1:200)=1;
%figure;plot(C);
% save C_NIST_532_600.mat C_NIST_532_600;


% NOW CALIBRATE THE DATA USING C
temp = sgolayfilt(White_Light_Background,3,41);%first subtract the background light as we dont want to correct (enhance!!) this
background = temp*min(spectra_raw)/min(temp);
spectra_calibrated = (spectra_raw-background).*C;
%figure;plot(background);

if exist('glass', 'var') == 1
    figure;plot(glass);hold;
    White_Light_Background=sgolayfilt(White_Light_Background,3,51);
    background = White_Light_Background*min(glass)/min(White_Light_Background);
    glass = (glass-background).*C; glass = sgolayfilt(glass,3,51); %may need to be commented out if not relevant - if you dont need glass
    plot(background);plot(glass);hold off;
end

save temp.mat glass Four_Ace spectra_raw spectra_calibrated
clear;
load temp.mat

%Author Bryan Hennelly
%Last Modified 23/07/15

%This function takes in a matrix of spectra and finds the best fit of the
%spectrum minus the backgorund and a polynomial .

%The method is based on using a modified version of the algorithm presented
%in B. Beier and A. Berger, Analyst, 2009, 134, 1198–1202.

%OUTPUTS
%The returned data is the spectrum minus the backgorund as well as the
%calculed backgorund itself. data is of dimension P x Q, see INPUTS 
%Note no normalisation is performed

%INPUTS:
%Spectrum is the spectrum data and can be a N x M matrix made up of N
%spectra of size M
%backgorund should be of size M
%N is the order of the polynomial e.g. 5 or 7
%iterations is the number of iteartions to be used to find the ebst fit
%under the spectrum. 
%area normalise should be set to 'y' or 'n' depending on whether area
%normalisation is preferred 
%print should be set to 'y' or 'n' depending on whether you wish to print
%the results 

function [data, backgrounds] = modified_berger_background_reduction(spectrum, background, N, iterations, area_normalise, print) 

%Set up output matrices
data = spectrum;
backgrounds = spectrum;

% Backgorund data - smoothe
Background = sgolayfilt(background,3,41);


%Perform reduction

    %First find a good estimate for the backgorund subtraction by applying
    %an N order polynomial + backgorund least squares fit to the first half
    %of the spectrum where the backgorund has its biggest effect. 
    [dummy,p] = rubber_band_adapted_BH(spectrum,Background, N, 0, 30, 'n');
    %The estimated concentration is given by p(order + 2. Next aply
    %fminsearch to fine tune the choice of concetration 
    x = fminsearch(@(x) find_best_fit(spectrum,Background ,x, 7, iterations),p(N+2));
    %finally apply runnerbad method to spectrum minus the estimate backgorund signal
    [data,p] = rubber_band(spectrum-x*Background, 7, iterations);
    
    backgrounds = spectrum-data;
    
    %Print results if required 
    if print == 'y' 
        figure;plot(spectrum,'b');hold;plot(data,'r');plot(backgrounds,'g');hold off;
    end
    
    
    
    %Perform area normalisation if required
    if area_normalise == 'y' 
        spectrum=spectrum/sum(spectrum);
        data=data/sum(data);
        backgrounds=backgrounds/sum(backgrounds);
    end
    


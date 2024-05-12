%Author Bryan Hennelly
%Last Modified 23/07/15

%This function takes in a matrix of spectra and finds a least squares best fit of the
%spectrum minus the backgorund and a polynomial .

%OUTPUTS
%The returned data is the spectrum minus the backgorund as well as the
%calculed backgorund itself. data is of dimension P x Q, see INPUTS 
%Note no normalisation is performed

%INPUTS:
%Spectrum is the spectrum data and can be a N x M matrix made up of N
%spectra of size M
%backgorund should be of size M
%N is the order of the polynomial e.g. 5 or 7
%nn is the order for the modulation of the polynomial modulating the
%backgorund signal, by default this should be set to zero
%itearations is the number of times tyhe algorithm will be run to fit the
%minimum of the signal and the previous result.

%print should be set to 'y' or 'n' depending on whether you wish to print
%the results 

function [back_spectrum,p] = rubber_band_adapted_BH(spectrum, background, N, nn, iterations, print)


%%
%THIS PART SHOULD ONLY BE USED IF WE ARE INTERESTED ONLY IN THE FIRST
%SECTION OF THE SPECTRUM WHERE THE BACKGORUND IS MORE DOMINANT. THIS MAY
%NEED TO BE ADAPTED DEPENDING ON THE WAVENUMBER AXIS OF A GIVEN SPECTRUM
spectrum = spectrum(1:258);
background = background(1:258);

%%
%Set up axis. Essentially we make sure that the axis used here
%matches the one we will use later in rubber_band where we use the entire
%spectrum

axis = (1:1:1024);
%%
%Set up matrix for retuned spectra
back_spectrum = spectrum;


    data = spectrum;

    %optikonal print
    if print == 'y' figure;hold;plot(spectrum,'r');  end 

    for m  = 1 : iterations
            % Apply least squares fit using polynomial and background
            [solution,p] = adapted_polyfit(axis, data, background, N, nn); % compute polynomial
            
            %Optional print
            if print == 'y' plot(solution, 'b'); pause(0.15); end
            
            %Set data to be the minimum of calculated backgorund and
            %previous calculation 
            temp(1,:) = spectrum;
            temp(2,:) = solution;
            data = min(temp); 
    end
    
    %Return the corrected data
    back_spectrum=spectrum - data;
    
    hold off;





end
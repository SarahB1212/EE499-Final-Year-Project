%Author: Bryan Hennelly NUIM 01/12/14
%Fucntion Wavelength_Calibrate will take in a datset Data made up of NxM
%points where N is the number of spectra of length M
%p is an estimate of the pixel number where the phenylalanine peak can be
%found. The function will return a new Dataset of same size as Data where
%the exact positions of the phenylalanine peaks are not all in the same position,
%i.e. at the  pixel point that was initially specified. The peak positions 
%have been calculated by
%analysing using cubic spline interpolation. Following this sub_pixel
%shifting is implemented using the DFT in order to ensure that the
%phenylalanine peak does fall on the pixel number p for each of the
%spectra. See the find_maximum for a discussion on how this is implemented.
%The other inputs are see_peaks which should be 'y' or 'n' and defines whetehr or
%not to display the cubic spline interpolation function for the
%phenylalanine peak in each of the spectra both before and after subpixel
%shifting. The delay parameter is the time in seconds that the function
%will delay between displaying each of these figures.

function new_Data = Peak_Calibrate(Data,p,see_peaks,delay)
Data = Data';
%Data = new_Data;

S = size(Data);

window = 4;

if(see_peaks == 'y') 
    figure;hold; 
    title('cubic spline interpolation in region of specified peaks');xlabel('pixels');ylabel('Intensity AU'); 
end
out = zeros(S(1),1);
for(n=1:1:S(1))
    pause(delay);
    out(n) = find_maximum(Data(n,p-window:p+window),p,see_peaks);
end

new_Data = Data;
f = 0:1:S(2)-1+20;
for(n=1:1:S(1))
    shift = out(n)/(length(Data(n,:)));
    lin_phase = exp(i*f*2*pi*shift);
    
    step = 21; diff = (Data(n,end) - Data(n,1))/step; a = Data(n,1); 
    FT = fftshift(fft([a+10*diff,a+9*diff,a+8*diff,a+7*diff,a+6*diff,a+5*diff,a+4*diff,a+3*diff,a+2*diff,a+diff,Data(n,:),a+20*diff,a+19*diff,a+18*diff,a+17*diff,a+16*diff,a+15*diff,a+14*diff,a+13*diff,a+12*diff,a+11*diff])); %we have zeropadded the input data to allow for 5 pixels of wrap around
  
    FT = FT.*lin_phase;
    FT = (ifft(FT));
    new_Data(n,:) = abs(FT(11:end-10));
    %new_Data(n,:) = new_Data(n,:)/sum(new_Data(n,:));
end

new_Data=new_Data;

out_new = out;
if(see_peaks == 'y') 
    figure;hold; 
    title('cubic spline interpolation in region of specified peaks after subpixel shifting');xlabel('pixels');ylabel('Intensity AU');
end

for(n=1:1:S(1))
    pause(delay);
    out_new(n) = find_maximum(new_Data(n,p-window:p+window),p,see_peaks);
end

if(S(1)>1)
    figure;title('shift of phenalanaline relative to the specified pixel (blue before, red after subpixel shifting)');xlabel('spectrum index');ylabel('deviation from specified pixel, units pixels');
    hold
    f=1:S(1);
    plot(f,out,'b',f,out_new,'r');
end

new_Data = new_Data';
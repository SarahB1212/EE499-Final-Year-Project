%% PLEASE REFERENCE THE FOLLOWING PAPER FOR NEON WAVELENGTH CALIBRATION
%Liu, Dongyue, and Bryan M. Hennelly. "Improved wavelength calibration by modeling the spectrometer." Applied Spectroscopy 76.11 (2022): 1283-1299.

%% BE CAREFUL TO ENTER THE RIGHT POSITIONS AND WAVELENGTHS
%You must ensure these first two lines are correct for your spectrum
neon_peak_positons=[641 680 765 808 886 953 987]; % You must enter the approximate pixel at which you can see the four Neon peaks
wavelength_NIST=[585.249 588.190 594.483 597.553 603.000 607.434 609.616]; %...and the corresponding wavelengths in neon
neon=neon';


wavelength_NIST=wavelength_NIST.*1e-9;
peak_length=size(wavelength_NIST,2);
wind1=ones(1,peak_length).*3; wind2=ones(1,peak_length).*3; 
for i=1:1
    for j=1:peak_length
        Min1(i,j)=1e14; %check
        for k=0.1:1:100
            p=neon_peak_positons(i,j)-wind1(j);
            q=neon_peak_positons(i,j)+wind2(j);    
            y=neon(p:q); 
            x=p:1:q;
            P2=neon_peak_positons(i,j);P3=k;C=0;
            [height,pos]= max(y);
            P1=P3*height;
            P0=[P1,P2,P3,C];
            [yprime2 P resnorm2 residual2] = lorentzfit(x,y,P0);
            [peak_positions_subpixel(i,j)]=P(2);
            xx=p:0.01:q;
            YPRIME = P(1)./((xx - P(2)).^2 + P(3)) + P(4);
            line=[xx',YPRIME'];
            [index,temp1]=find(x==line(:,1));
            test=sum(abs(line(index,2)-y'));
            if test<=Min1(i,j)
                
               Min1(i,j)=test;
               chosen_k(i,j)=k;
               chosen_P(i,j)={P};
               [peak_positions_subpixel2(i,j)]=P(2);
            end
        end
    end
end

temp5=polyfit(peak_positions_subpixel2,wavelength_NIST,2);
temp6=polyval(temp5,peak_positions_subpixel2);
error4=mean(sqrt((temp6-wavelength_NIST).^2));
wavelength_axis=polyval(temp5,(1:1:1024));
%plot(wavelength_axis)

save wavelength_axis

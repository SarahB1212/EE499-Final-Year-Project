%% PLEASE REFERENCE THE FOLLOWING PAPER FOR Four_Ace CALIBRATION
% Hutsebaut, Didier, Peter Vandenabeele, and Luc Moens. "Evaluation of an accurate calibration and spectral standardization procedure for Raman spectroscopy." Analyst 130.8 (2005): 1204-1214.

%% BE CAREFUL TO ENTER THE RIGHT POSITIONS AND WAVELENGTHS
%You must ensure these first two lines are correct for your spectrum
ace_peak_positons=[132 221 295 336 387 412 439 474 493 611];
%ace_peak_positons=[132 221 294 336 387 412 439 474 493 611];
wavenumber_NIST=[390.2 651.6 857.9 968.7 1105.5 1168.5 1236.8 1323.9 1371.5 1648.4];
Four_Ace=Four_Ace';

peak_length=size(wavenumber_NIST,2);
wind1=ones(1,peak_length).*3; wind2=ones(1,peak_length).*3; 
for i=1:1
    for j=1:peak_length
        Min1(i,j)=1e14; %check
        for k=0.1:1:100
            p=ace_peak_positons(i,j)-wind1(j);
            q=ace_peak_positons(i,j)+wind2(j);    
            y=Four_Ace(p:q); 
            x=p:1:q;
            P2=ace_peak_positons(i,j);P3=k;C=0;
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

temp5=polyfit(peak_positions_subpixel2,wavenumber_NIST,2);
temp6=polyval(temp5,peak_positions_subpixel2);
error4=mean(sqrt((temp6-wavenumber_NIST).^2));
wavenumber_axis=polyval(temp5,(1:1:1024));
%plot(wavenumber_axis)

clear temp5 wavenumber_NIST temp6 error4 peak_positions_subpixel2 chosen_P chosen_k Min1 test k P index y x index temp1
clear line xx YPRIME p q peak_positions_subpixel yprime2 resnorm2 residual2 P0 P1 P2 P3 i j k wind1 peak_length wind2 pos height C

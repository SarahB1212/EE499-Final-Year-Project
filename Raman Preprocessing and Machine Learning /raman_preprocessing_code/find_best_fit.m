%Author Bryan Hennelly
%Last Modified 23/07/15
%Perfomrs a minimisation on x (being the concentration of the substrate
%backgorund) on some crieterium as described by the line f = ....
%This criterium can be changed to meet a given requiremnt 
%Berger B. Beier and A. Berger, Analyst, 2009, 134, 1198–1202. bases it on 
%minimising the area between the estimated backgorund and the spectrum 
%here (10 times) more weight is given to those regions in which the backgorund signal
%is strongest..this will need to be adapte to suit a given spectrum/background 

function f = find_best_fit(spectrum,background,x, N, iterations)

[data,p] = rubber_band(spectrum-x*background, N, iterations);
f = sum(abs(data(1,:)));% + 10*(sum(abs(data(1:120)))+sum(abs(data(200:250))));

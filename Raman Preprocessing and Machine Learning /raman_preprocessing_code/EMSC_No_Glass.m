%This code uses the same notation as in Section 2.
%Assumes that the inputs r and g have been smoothed before inputting using
%g=sgolayfilt(raw_glass_spectrum,3,11) or some other smoothing filter.
%Outputs the corrected signal as well as the total background and other
%info. Assumes S,r,and g are all 1D vectors (just spectral irradiance
%values)
%%INPUTS
%S is the spectrum to be corrected (no need for smoothing in advance)
%r is a reference cell spectrum that is free of the glass signal (should be
%smooted in advance)
%g is the glass spectrum (also should be smoothed)
%N is the order to the background polymonial to be included (%Even N = 1,
%i.e.  straight line often works well)
%%OUTPUTS
%s is the returned signal with some weight of glass removed and an N order
%polymonial removed
%background is the total backgorund signal (glass + poly)
%c_r and c_g are the wights of r and g 
%B_N is the isolated polynomial signal in the backgorund. 
function [s,background,c_r,B_N] = EMSC_No_Glass(S,r,N)
%dummy axis
x = 1:length(S);
mu = [mean(x); std(x)];
x = (x - mu(1))/mu(2);
x=x';
% Construct Vandermonde matrix.
V(:,N+1) = ones(length(x),1,class(x));
for j = N:-1:1
   V(:,j) = x.*V(:,j+1);
end
V(:,N+2) = r;
% Solve least squares problem.
[Q,R] = qr(V,0);
p = R\(Q'*S'); 
% Calculate outputs
c_p = p(1:N+1);B_N=polyval(c_p,x);
c_r=p(N+2);
background=B_N;
s=(S-background')/c_r;



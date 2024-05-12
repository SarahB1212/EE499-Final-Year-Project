%%S is now a matrix of spectra. r and g are still 1D vectors

function [s, back, C_R, C_G] = EMSC_Matrix_Glass(S,r,g,N)

s=S;
back=S;

SIZE = size(S);



C_R = zeros(1,SIZE(2));
C_G = zeros(1,SIZE(2));
C_N = zeros(1,SIZE(2));

for n = 1:1:SIZE(2)
    X = S(:,n);[temp,background,c_r,c_g, B_N] = EMSC_Glass(X',r',g',N); %watch out for transposing r and g
    s(:,n)=temp';
    back(:,n)=background';
    C_R(n)=c_r;
    C_G(n)=c_g;
end



%%S is now a matrix of spectra. r and g are still 1D vectors

function [s, back, C_R, C_G] = EMSC_Matrix_No_Glass(S,r,N)

s=S;
back=S;

SIZE = size(S);



C_R = zeros(1,SIZE(2));
C_G = zeros(1,SIZE(2));
C_N = zeros(1,SIZE(2));

for n = 1:1:SIZE(2)
    X = S(:,n);
    [temp,background,c_r,B_N] = EMSC_No_Glass(X',r',N); %watch out for transposing r and g
    s(:,n)=temp';
    back(:,n)=background';
    C_R(n)=c_r;
end



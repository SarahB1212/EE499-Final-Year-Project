%Author Bryan Hennelly
%Last Modified 23/07/15

%   This function finds the coefficients of a polynomial P(X) of
%   degree N as well as a second polynomial of degree nn that modulates the backgorund signal 
%   that fits the data Y best in a least-squares sense. P is a
%   row vector of length nn + N + 1 containing the polynomial coefficients in
%   descending powers, P(1)*X^N + P(2)*X^(N-1) +...+ P(N)*X + P(N+1)+
%   backgorund.*(P(N+2)+ P(N+3)*X^1 + P(N+4)*X^(2) +...+ P(N+1+nn)*X^(nn) +

%   This function takes in a matrix of spectra and finds a least squares best fit of the
%   spectrum minus the backgorund and a polynomial .

%OUTPUTS
%The returned data is the best fit as well as the polynomial coefficients 

%INPUTS:
%y is the spectrum data and can be a 1 x M matrix
%x is the axis to calculate the polynomials on
%backgorund should be of size M
%n is the order of the polynomial e.g. 5 or 7
%nn is the order for the modulation of the polynomial modulating the
%backgorund signal, by default this should be set to zero

function [fit,p] = adapted_polyfit(x,y,background,n,nn)

if ~isequal(size(x),size(y))
    error('MATLAB:polyfit:XYSizeMismatch',...
          'X and Y vectors must be the same size.')
end

x = x(:);
y = y(:);

if nargout > 2
   mu = [mean(x); std(x)];
   x = (x - mu(1))/mu(2);
end

% Construct Vandermonde matrix.

V(:,n+1) = ones(length(x),1,class(x));
for j = n:-1:1
   V(:,j) = x.*V(:,j+1);
end

V(:,n+2) = background;

for j = n+3:1:n+nn+2
   V(:,j) = x.*V(:,j-1);
end



% Solve least squares problem.
[Q,R] = qr(V,0);
ws = warning('off','all'); 
p = R\(Q'*y);    % Same as p = V\y;

fit = V*p;


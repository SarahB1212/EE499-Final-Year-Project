%Author: Bryan Hennelly NUIM 01/12/14
%Fucntion find_maximum will take in a vector y and will find the maximum
%point as found by create a cubic spline interpolation of the dataset,
%differentiating that and finding the global maxima by investigating each
%if the cubic splines derivatives individually
%the input parameter see_peaks shold be 'y' or 'n' and defines whetehr or
%not to display the cubic spline interpolation function
% p is the pixel value of the estimated position of the phenylanaline in
% the original dataset

function out = find_maximum(y,p,see_peaks)
y=y-min(y);y = y/max(y); %normalise
L = length(y);%length should be odd
x = p-(L-1)/2:p+(L-1)/2;
cs = spline(x,[0 y 0]);
%cs = pchip(x, y);
xx = linspace(p-(L-1)/2,p+(L-1)/2,10001);

if(see_peaks == 'y') plot(x,y,'o',xx,ppval(cs,xx),'-'); end

S = size(cs.coefs);

max_x = 0;  %this is the point at which we pind the peak of the phenalanaline peak
max_value = 0;  %and this is the value at that point


for n = 1:1:S(1)
        % find quadratic equation that is derivative of each spine
        % polynomial
        a = 3*cs.coefs(n,1);
        b = 2*cs.coefs(n,2);
        c = 1*cs.coefs(n,3);
        % solve quadratic to find any maxima/minima
        x0_1 = 0;
        x0_2 = 0;
        if (b^2 - (4*a*c)) > 0
            x0_1 = (-1*b + sqrt(b^2 - (4*a*c)))/(2*a);
            x0_2 = (-1*b - sqrt(b^2 - (4*a*c)))/(2*a);
        end
        % make sure the minima/maxima that are discovered are within the
        % boiund of that given polynomial and chch the value of that point
        % to see if it is the maximum value
        lower_bound = 0;
        upper_bound = 1;
        % first x0_1
        if (x0_1 >= lower_bound) && (x0_1 < upper_bound)
            value_at_maximum = (x0_1^3)*cs.coefs(n,1) + (x0_1^2)*cs.coefs(n,2) + (x0_1)*cs.coefs(n,3) + cs.coefs(n,4);
            if (value_at_maximum > max_value)
                max_value = value_at_maximum;
                max_x = x0_1 + n - 1;
            end
        end
        %then x0_2
        if (x0_2 >= lower_bound) && (x0_2 < upper_bound)
            value_at_maximum = (x0_2^3)*cs.coefs(n,1) + (x0_2^2)*cs.coefs(n,2) + (x0_2)*cs.coefs(n,3) + cs.coefs(n,4);
            if (value_at_maximum > max_value)
                max_value = value_at_maximum;
                max_x = x0_2 + n - 1;
            end
        end
end

out = max_x - (L-1)/2;
        
        
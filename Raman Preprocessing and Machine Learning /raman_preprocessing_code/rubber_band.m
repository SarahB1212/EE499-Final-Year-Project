function [back_spectrum,p] = rubber_band(spectrum, N, iterations)
%N is the order


data = spectrum;
   
    for m  = 1 : iterations
            % Apply polynomial
%m
            %plot(data,'r');
%            size((1:1:1024)')
%            size(data)
            p = polyfit((1:1:1024)', data, N ); % compute polynomial
            solution = polyval(p,(1:1:1024)');
            
            temp(1,:) = spectrum;
            temp(2,:) = solution;
            


 %           plot(solution, 'b');
           % pause(0.15); 
            data = min(temp); 
            
%            temp(1,:) = data;
%            temp(2,:) = zeros(size(data));
            
%            data = max(temp); 
            
            data = data';
            %data(1:200)=spectrum(1:200);%anchor
            %data(800:end)=spectrum(800:end);%anchor
%            
%              if m ~= iterations 
%                
%                 data(500:600) = spectrum(n,500:600);
%                 
% 
%             end
%             
            
    end
   
    
    back_spectrum=spectrum - data;
    
    


end
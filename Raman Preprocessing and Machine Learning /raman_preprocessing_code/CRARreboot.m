function[dataset] = CRARreboot(dataset)

%Determine dimensions of matrix for reference
dim = size(dataset);

%figure, plot(dataset.')
dim(1)
% Median filtering to prevent CRA matching in high CRA data
for h = 1 : 1 : dim(1)
    refData(h,:) = medfilt1(dataset(h,:),11);
end

%figure, plot(refData.')

%Ascertain the correlation between each spectrum
%For each row of data
for i = 1 : dim(1)
    %Correlate that spectrum with every other spectrum
    for j = 1 : dim(1)
        %Brekaing down the correlation calculation
        d1 = dot(refData(i,:), refData(j,:));
        d2 = dot(refData(i,:), refData(i,:));
        d3 = dot(refData(j,:), refData(j,:));
        corrVec(j) = (d1^2/(d2*d3));
        
        %If a spectrum is going to be correlated with itself, record a null
        %value
        if i == j
            corrVec(j) = 0;
        end
    end
    
    %Determine the optimum pairings of spectra (corrVec should have the same numbber of elements as rows in the dataset so dim is reused)
    for k = 1 : dim(1)
        %If a particular element is equal to the max of the vector then
        %this is the most similar spectrum and the position of it is recorded in pairVec
        if corrVec(k) == max(corrVec)
            pairVec(i) = k;
        end
    end
end


sigma = mean(std(dataset(:,700:1024)));

%For every spectrum
for h = 1 : dim(1)
    
    
    
    %Subtract the correlated spectrum from the original
    tempSp = dataset(h,:) - dataset(pairVec(h),:);
    %figure;plot(dataset(h,:));hold;plot(dataset(pairVec(h),:));plot(tempSp);

    
    %If there is a significant difference between them then copy over the
    %values
    for g = 1 : dim(2)
        if tempSp(g) > sigma*5
            dataset(h, g) = dataset(pairVec(h), g);
            
            %If a cosmic ray is detected then apply lower filters to the
            %surrounding area
            for p = -1 : 1 : 1
                if g+p < 1 || g+p > dim(2)
                    continue
                else
                    if tempSp(g+p) > sigma*2
                        dataset(h, g+p) = dataset(pairVec(h), g+p);
                    end 
                end
            end
        
        end
    end
end


%figure, plot(dataset.')

end
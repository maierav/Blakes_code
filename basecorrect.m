function [BL] = basecorrect(data)
%Basecorrect - baseline correct the indexed data
%   Detailed explanation goes here
basetp = data(25:75,:); %retrieves the reference range to which the data will be baseline corrected
mean_basetp = mean(basetp, 1); 
BL = data(1:351,:)- mean_basetp; 
end


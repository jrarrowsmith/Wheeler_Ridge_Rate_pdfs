function [increment,range_of_values]= compute_increment_range(maxvalue, minvalue, n)
%This function computes the increment for the sample range given the min
%and max and the number of samples
%JRA October 20, 2016
increment=(maxvalue-minvalue)/(n-1);
range_of_values=[minvalue:increment:maxvalue];
end


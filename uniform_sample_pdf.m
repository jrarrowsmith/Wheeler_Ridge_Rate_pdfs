function [b_range, b, b_histogram] = uniform_sample_pdf(b_max, b_min, n)
%This function computes the age pdf.
%It assumes a uniform distribution across the range min to max and is
%sampled n times
%JRA October 20, 2016 updated Dec 5, 2020
[b_inc,b_range]= compute_increment_range(b_max, b_min, n);
b = b_min + (b_max-b_min).*rand(length(b_range),1);
b_histogram = hist(b, length(b_range));
end


function [rates] = compute_rates(younger_pos, older_pos, younger_age, older_age)
%This function computes the incision rates given elevation - age data
%JRA June 30, 2017 adapted December 5, 2020
%it also plots the resulting histogram

delta_pos=(older_pos-younger_pos);
delta_t=(older_age-younger_age);
locs=find(delta_t>0 & delta_pos>0); %We don't allow negative rates
delta_t=abs(delta_t(locs));
rates = delta_pos(locs)./delta_t;
end



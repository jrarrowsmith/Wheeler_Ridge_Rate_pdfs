function [sample_points_pos, sample_points_age]=collect_points_on_trajectories_continuous(ages, pos,n, step_size)
%This function collects points on the sampled trajectories for later 2D
%histograms. It marches through all of the input files
%J R Arrowsmith, December 2020

sample_points_pos=[];
sample_points_age=[];

% vq = interp1(x,v,xq) returns interpolated values of a 1-D function at
% specific query points using linear interpolation. Vector x contains the
% sample points, and v contains the corresponding values, v(x). Vector xq
% contains the coordinates of the query points.
xq=min(ages(:,1)):step_size:max(ages(:,end));

for i=1:n
  rates=ages(i,:)./pos(i,:);  
  pos_rates=rates>=0;
  num_positive=sum(pos_rates);
  if num_positive<5
      sprintf('hello %0.1f',num_positive)
  end
vq = interp1(ages(i,:),pos(i,:),xq);
sample_points_age=[sample_points_age xq];
sample_points_pos=[sample_points_pos vq];
end


end


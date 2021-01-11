%Simple script to explore Monte Carlo modeling of erosion or uplift or propagation histories
%JRA October 14 and December 5, 2020 adapted from Nagra
%Nagra_Incision_model_2D_pdf_ZNOdata_20181021.m
%%these are updates from Emily: "WITH AND WITHOUT Q4b/Q5b surfaces"

%X = rand(n) returns an n-by-n matrix of random numbers uniformly distributed from 0-1
%X = randn(n) returns an n-by-n matrix of random numbers normally
%distributed mean 0 standard deviation 1
clear all
close all

%Set variables
n=1000000; %number of samples
num_std=2; %number of standard deviations to sample for the elevations. We assume normal distributions for the elevations and uniform distributions for the ages
truncate = 1; %truncate = 1 yes cut the distribution at the 2 sigma or no if = other
sigmas=0.9545; %0.6827 is one sigma and 0.9545 is two sigma range to cut pdf 
lineweight=0.01; %lineweight for connections in figure 3 plot
hist_bins=100; %number of bins for the histogram figures
%plot_background_color=[211/255 211/255 211/255]; %light grey
plot_background_color=[80/255 80/255 80/255]; %dark grey
portion_of_connections = floor(n/50); %the denominator is how many of the total connections to draw as lines
SU_step_size = 0.05; %step size in kyr for the trajectories for SU
propagation_step_size = 0.5; %step size in kyr for the trajectories for propagation
grid_steps=100; %Gridding for the 2D histogram
max_SU_rate=10;%max rate in m/kyr for plotting
max_propagation_rate=150;%max rate in m/kyr for plotting
print_save=1; %1 to save files and print out the diagrams, 0 or other number for no
tic
%%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%Age, elevations, and propagations for the surfaces
%deal with with and without
run_name='-withoutQ4bQ5bv2';
run_name=strcat(date,run_name,sprintf('-%.4f',sigmas));
%Q4_SU_min=139; %Q4_surface uplift minimum in m WITH
Q4_SU_min=170; %Q4_surface uplift minimum in m WITHOUT
%Q5_SU_min=304; %Q5_surface uplift minimum in m WITH
Q5_SU_min=313; %Q5_surface uplift minimum in m WITHOUT

%Q1--no change from 2020 runs; same for WITH and WITHOUT
Q1rgb=[255/255 255/255 0/255]; %hex is FFFF00 https://www.countingcharacters.com/hex-to-rgb
Q1_mean_age=7.37; %Q1 surface mean age in ka
Q1_SD=0.14; %Q1 surface mean age in kyr
Q1_SU_min=0; %Q1_surface uplift minimum in m
Q1_SU_max=1; %Q1_surface uplift maximum in m
Q1_propagation_min=7750; %Q1 propagation minimum in m
Q1_propagation_max=7751; %Q1 propagation maximum in m

%Q2--age is updated from 2020 runs; same for WITH and WITHOUT
Q2rgb=[255/255 179/255 0/255]; %hex is FFB300 https://www.countingcharacters.com/hex-to-rgb
Q2_mean_age=32; %Q2 surface mean age in ka
Q2_SD=4.6; %Q2 surface mean age in kyr
Q2_SU_min=18; %Q2_surface uplift minimum in m
Q2_SU_max=23; %Q2_surface uplift maximum in m
Q2_propagation_min=6500; %Q2 propagation minimum in m
Q2_propagation_max=6610; %Q2 propagation maximum in m

%Q3--age is updated from 2020 runs; same for WITH and WITHOUT
Q3rgb=[255/255 128/255 77/255]; %hex is FF804D https://www.countingcharacters.com/hex-to-rgb
Q3_mean_age=70.3; %Q3 surface mean age in ka
Q3_SD=9.6; %Q3 surface mean age in kyr
Q3_SU_min=35; %Q3_surface uplift minimum in m
Q3_SU_max=70; %Q3_surface uplift maximum in m
Q3_propagation_min=5600; %Q3 propagation minimum in m
Q3_propagation_max=6300; %Q3 propagation maximum in m

%Q4--age and SU min is updated from 2020 runs; switch for with and without
Q4rgb=[168/255 112/255 0/255]; %hex is A87000 https://www.countingcharacters.com/hex-to-rgb
Q4_mean_age=128.7; %Q4 surface mean age in ka
Q4_SD=22.5; %Q4 surface mean age in kyr
Q4_SU_max=200; %Q4_surface uplift maximum in m
Q4_propagation_min=3000; %Q4 propagation minimum in m
Q4_propagation_max=3600; %Q4 propagation maximum in m

%Q5--age and SU min is updated from 2020 runs; switch for with and without
Q5rgb=[115/255 38/255 0/255]; %hex is 732600 https://www.countingcharacters.com/hex-to-rgb
Q5_mean_age=153; %Q5 surface mean age in ka
Q5_SD=24; %Q5 surface mean age in kyr
Q5_SU_max=331; %Q5_surface uplift maximum in m
Q5_propagation_min=380; %Q5 propagation minimum in m
Q5_propagation_max=750; %Q5 propagation maximum in m


%%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%%Initialization
%Set up these variables:
sample_points_Age_SU=[];
sample_points_SU=[];
sample_points_Age_propagation=[];
sample_points_propagation=[];
%%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%Sample normally distributed ages
[Q1_age_range, Q1_age_samples, Q1_age_histogram] = normal_sample_pdf(Q1_mean_age, Q1_SD, num_std, n, truncate);
[Q2_age_range, Q2_age_samples, Q2_age_histogram] = normal_sample_pdf(Q2_mean_age, Q2_SD, num_std, n, truncate);
[Q3_age_range, Q3_age_samples, Q3_age_histogram] = normal_sample_pdf(Q3_mean_age, Q3_SD, num_std, n, truncate);
[Q4_age_range, Q4_age_samples, Q4_age_histogram] = normal_sample_pdf(Q4_mean_age, Q4_SD, num_std, n, truncate);
[Q5_age_range, Q5_age_samples, Q5_age_histogram] = normal_sample_pdf(Q5_mean_age, Q5_SD, num_std, n, truncate);

figure(1)
set(gca,'Color',plot_background_color)
hold on
histogram(Q1_age_samples, hist_bins, 'EdgeColor','none', 'FaceColor',Q1rgb)
histogram(Q2_age_samples, hist_bins, 'EdgeColor','none', 'FaceColor',Q2rgb)
histogram(Q3_age_samples, hist_bins, 'EdgeColor','none', 'FaceColor',Q3rgb)
histogram(Q4_age_samples, hist_bins, 'EdgeColor','none', 'FaceColor',Q4rgb)
histogram(Q5_age_samples, hist_bins, 'EdgeColor','none', 'FaceColor',Q5rgb)
xlabel('age (ka)')
ylabel('counts')
title({'Normally distributed age ranges for surfaces',run_name})

%Sample uniformly distributed Surface Uplift
[Q1_SU_range, Q1_SU_samples, Q1_SU_histogram] = uniform_sample_pdf(Q1_SU_max, Q1_SU_min, n);
[Q2_SU_range, Q2_SU_samples, Q2_SU_histogram] = uniform_sample_pdf(Q2_SU_max, Q2_SU_min, n);
[Q3_SU_range, Q3_SU_samples, Q3_SU_histogram] = uniform_sample_pdf(Q3_SU_max, Q3_SU_min, n);
[Q4_SU_range, Q4_SU_samples, Q4_SU_histogram] = uniform_sample_pdf(Q4_SU_max, Q4_SU_min, n);
[Q5_SU_range, Q5_SU_samples, Q5_SU_histogram] = uniform_sample_pdf(Q5_SU_max, Q5_SU_min, n);

figure(2)
set(gca,'Color',plot_background_color)
hold on
histogram(Q1_SU_samples, hist_bins, 'EdgeColor','none', 'FaceColor',Q1rgb)
histogram(Q2_SU_samples, hist_bins, 'EdgeColor','none', 'FaceColor',Q2rgb)
histogram(Q3_SU_samples, hist_bins, 'EdgeColor','none', 'FaceColor',Q3rgb)
histogram(Q4_SU_samples, hist_bins, 'EdgeColor','none', 'FaceColor',Q4rgb)
histogram(Q5_SU_samples, hist_bins, 'EdgeColor','none', 'FaceColor',Q5rgb)
xlabel('surface uplift (m)')
ylabel('counts')
title({'Uniformly distributed surface uplift ranges for surfaces',run_name})


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%start section on surface uplift
figure(4)
set(gca,'Color',plot_background_color)
set(gca, 'ydir', 'reverse' )
xlabel('Age (ka)')
ylabel('Surface uplift (m)')
hold on
plot(Q1_age_samples,Q1_SU_samples,'.','color',Q1rgb)
plot(Q2_age_samples,Q2_SU_samples,'.','color',Q2rgb)
plot_connections(Q2_age_samples,Q1_age_samples,Q2_SU_samples,Q1_SU_samples, portion_of_connections, n, lineweight)
plot(Q3_age_samples,Q3_SU_samples,'.','color',Q3rgb)
plot_connections(Q3_age_samples,Q2_age_samples,Q3_SU_samples,Q2_SU_samples, portion_of_connections, n, lineweight)
plot(Q4_age_samples,Q4_SU_samples,'.','color',Q4rgb)
plot_connections(Q4_age_samples,Q3_age_samples,Q4_SU_samples,Q3_SU_samples, portion_of_connections, n, lineweight)
plot(Q5_age_samples,Q5_SU_samples,'.','color',Q5rgb)
plot_connections(Q5_age_samples,Q4_age_samples,Q5_SU_samples,Q4_SU_samples, portion_of_connections, n, lineweight)
title(run_name)

[sample_points_SU, sample_points_Age_SU]=collect_points_on_trajectories_continuous([Q1_age_samples Q2_age_samples Q3_age_samples Q4_age_samples Q5_age_samples],...
    [Q1_SU_samples Q2_SU_samples Q3_SU_samples Q4_SU_samples Q5_SU_samples],n, SU_step_size);

%now set up the binning for the hist3
x = linspace(min(sample_points_Age_SU),max(sample_points_Age_SU),grid_steps);
dx = (max(sample_points_Age_SU)-min(sample_points_Age_SU))./(grid_steps+1);
y = linspace(min(sample_points_SU),max(sample_points_SU),grid_steps);
dy = (max(sample_points_SU)-min(sample_points_SU))./(grid_steps+1);
age_SU_pts = [sample_points_Age_SU; sample_points_SU]';
ctrs={x;y};
clear sample_points_SU sample_points_Age_SU
N=hist3(age_SU_pts,ctrs);
clear age_SU_pts  
N(find(~N))=NaN;
volume_of_each_bar=N.*dx.*dy;
v =sum(sum(volume_of_each_bar,'omitnan'),'omitnan');
N=N./v;

figure(5)
set(gca, 'ydir', 'reverse' )
hold on
h=pcolor(x,y,N');
set(h,'edgecolor','none')
xlabel('age (ka)')
ylabel('surface uplift (m)')
title({'Surface uplift paths probability',run_name})
colorbar

figure(6)
set ( gca, 'ydir', 'reverse' )
hold on
h=pcolor(x,y,log10(N'));
set(h,'edgecolor','none')
xlabel('age (ka)')
ylabel('surface uplift (m)')
title({'Surface uplift paths  (log10 probability)',run_name})
colorbar

%Calculate the resulting rates
[Q1toQ2_SU_rates] = compute_rates(Q1_SU_samples, Q2_SU_samples, Q1_age_samples, Q2_age_samples);
[Q2toQ3_SU_rates] = compute_rates(Q2_SU_samples, Q3_SU_samples, Q2_age_samples, Q3_age_samples);
[Q3toQ4_SU_rates] = compute_rates(Q3_SU_samples, Q4_SU_samples, Q3_age_samples, Q4_age_samples);
[Q4toQ5_SU_rates] = compute_rates(Q4_SU_samples, Q5_SU_samples, Q4_age_samples, Q5_age_samples);
edges=0:SU_step_size:max_SU_rate;
figure(7)
clf
hold on
subplot(4,1,1)
Q1toQ2_SU_h=histogram(Q1toQ2_SU_rates, edges)
axis([0 max_SU_rate 0 inf])
xlabel('Q1 to Q2 SU rates')
title(run_name)

subplot(4,1,2)
Q2toQ3_SU_h=histogram(Q2toQ3_SU_rates, edges)
axis([0 max_SU_rate 0 inf])
xlabel('Q2 to Q3 SU rates')

subplot(4,1,3)
Q3toQ4_SU_h=histogram(Q3toQ4_SU_rates, edges)
axis([0 max_SU_rate 0 inf])
xlabel('Q3 to Q4 SU rates')

subplot(4,1,4)
Q4toQ5_SU_h=histogram(Q4toQ5_SU_rates, edges)
axis([0 max_SU_rate 0 inf])
xlabel('Q4 to Q5 SU rates (mm/yr)')
ylabel('counts')

figure(16)
label_pos=8;
clf
hold on
subplot(4,1,1)
clean_and_process_pdf(Q1toQ2_SU_h, sigmas, run_name, label_pos)
axis([0 max_SU_rate 0 inf])
xlabel('Q1 to Q2 SU rates')
title(run_name)
subplot(4,1,2)
clean_and_process_pdf(Q2toQ3_SU_h, sigmas, run_name, label_pos)
axis([0 max_SU_rate 0 inf])
xlabel('Q2 to Q3 SU rates')
subplot(4,1,3)
clean_and_process_pdf(Q3toQ4_SU_h, sigmas, run_name, label_pos)
axis([0 max_SU_rate 0 inf])
xlabel('Q3 to Q4 SU rates')
subplot(4,1,4)
clean_and_process_pdf(Q4toQ5_SU_h, sigmas, run_name, label_pos)
axis([0 max_SU_rate 0 inf])
xlabel('Q4 to Q5 SU rates (mm/yr)')
ylabel('probability')


clear Q1_SU_range Q1_SU_samples Q1_SU_histogram
clear Q2_SU_range Q2_SU_samples Q2_SU_histogram
clear Q3_SU_range Q3_SU_samples Q3_SU_histogram
clear Q4_SU_range Q4_SU_samples Q4_SU_histogram
clear Q5_SU_range Q5_SU_samples Q5_SU_histogram

%working on the composite rates and sampling them
All_SU_rates= [Q1toQ2_SU_rates; Q2toQ3_SU_rates; Q3toQ4_SU_rates; Q4toQ5_SU_rates]; %Concatenate the rates

nn=length(All_SU_rates); %Length of the concatenated list of rates
numsamples=n; %Sample it the same number of times as the first sampling
m=ceil(rand(numsamples,1).*nn); %Choose randomly and evenly across the length of the concatenated list of rates
sampled_rates=[];
for i = 1:numsamples
    sampled_SU_rates(i)=All_SU_rates(m(i)); %Sample the composite pdf
end

figure(8) % this figure provides a check that we are getting the same composite rate pdf 
subplot(2,1,1)
h=histogram(All_SU_rates,edges) %Second, make a histogram of the concatenated rates using our bins; THIS IS WHAT WE USE GOING FORWARD
title({'All surface uplift rates--histogram of the concatenated rates',run_name})
subplot(2,1,2)
histogram(sampled_SU_rates,edges); %Third, histogram of the sampled rates; this is just to check
xlabel('surface uplift rate, mm/yr')
ylabel('counts')
title('All surface uplift rates--histogram of the sampled rates')
clear sampled_SU_rates 

figure(9)
[most_common_SU_rate, max_SU_rate, min_SU_rate] = clean_and_process_pdf(h, sigmas, run_name, 4.5);
title({'surface uplift rate distribution',run_name})
xlabel('surface uplift rate, mm/yr')

clear sample_points_Age_SU sampled_SU_rates

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%start working ont the horizontal propagation
%Sample uniformly distributed horizontal propagation
[Q1_propagation_range, Q1_propagation_samples, Q1_propagation_histogram] = uniform_sample_pdf(Q1_propagation_max, Q1_propagation_min, n);
[Q2_propagation_range, Q2_propagation_samples, Q2_propagation_histogram] = uniform_sample_pdf(Q2_propagation_max, Q2_propagation_min, n);
[Q3_propagation_range, Q3_propagation_samples, Q3_propagation_histogram] = uniform_sample_pdf(Q3_propagation_max, Q3_propagation_min, n);
[Q4_propagation_range, Q4_propagation_samples, Q4_propagation_histogram] = uniform_sample_pdf(Q4_propagation_max, Q4_propagation_min, n);
[Q5_propagation_range, Q5_propagation_samples, Q5_propagation_histogram] = uniform_sample_pdf(Q5_propagation_max, Q5_propagation_min, n);

figure(3)
set(gca,'Color',plot_background_color)
hold on
histogram(Q1_propagation_samples, hist_bins, 'EdgeColor','none', 'FaceColor',Q1rgb)
histogram(Q2_propagation_samples, hist_bins, 'EdgeColor','none', 'FaceColor',Q2rgb)
histogram(Q3_propagation_samples, hist_bins, 'EdgeColor','none', 'FaceColor',Q3rgb)
histogram(Q4_propagation_samples, hist_bins, 'EdgeColor','none', 'FaceColor',Q4rgb)
histogram(Q5_propagation_samples, hist_bins, 'EdgeColor','none', 'FaceColor',Q5rgb)
xlabel('horizontal propagation (m)')
ylabel('counts')
title({'Uniformly distributed horizontal propagation ranges for surfaces',run_name})

figure(10)
clf
title(run_name)
set(gca,'Color',plot_background_color)
%set(gca, 'ydir', 'reverse' )
xlabel('Age (ka)')
ylabel('Horizontal propagation (m)')
hold on
plot(Q1_age_samples,-Q1_propagation_samples,'.','color',Q1rgb)
plot(Q2_age_samples,-Q2_propagation_samples,'.','color',Q2rgb)
plot_connections(Q2_age_samples,Q1_age_samples,-Q2_propagation_samples,-Q1_propagation_samples, portion_of_connections, n, lineweight)
plot(Q3_age_samples,-Q3_propagation_samples,'.','color',Q3rgb)
plot_connections(Q3_age_samples,Q2_age_samples,-Q3_propagation_samples,-Q2_propagation_samples, portion_of_connections, n, lineweight)
plot(Q4_age_samples,-Q4_propagation_samples,'.','color',Q4rgb)
plot_connections(Q4_age_samples,Q3_age_samples,-Q4_propagation_samples,-Q3_propagation_samples, portion_of_connections, n, lineweight)
plot(Q5_age_samples,-Q5_propagation_samples,'.','color',Q5rgb)
plot_connections(Q5_age_samples,Q4_age_samples,-Q5_propagation_samples,-Q4_propagation_samples, portion_of_connections, n, lineweight)

[sample_points_propagation, sample_points_Age_propagation]=collect_points_on_trajectories_continuous([Q1_age_samples Q2_age_samples Q3_age_samples Q4_age_samples Q5_age_samples],...
    [Q1_propagation_samples Q2_propagation_samples Q3_propagation_samples Q4_propagation_samples Q5_propagation_samples],n, propagation_step_size);
 
%now set up the binning for the hist3
x = linspace(min(sample_points_Age_propagation),max(sample_points_Age_propagation),grid_steps);
dx = (max(sample_points_Age_propagation)-min(sample_points_Age_propagation))./(grid_steps+1);
y = linspace(min(sample_points_propagation),max(sample_points_propagation),grid_steps);
%dy = (max(sample_points_SU)-min(sample_points_SU))./(grid_steps+1);
dy = (max(sample_points_propagation)-min(sample_points_propagation))./(grid_steps+1);
age_propagation_pts = [sample_points_Age_propagation; sample_points_propagation]';
ctrs={x;y};
N=hist3(age_propagation_pts,ctrs);
clear age_propagation_pts
N(find(~N))=NaN;
volume_of_each_bar=N.*dx.*dy;
v =sum(sum(volume_of_each_bar,'omitnan'),'omitnan');
N=N./v;

figure(11)
clf
set(gca, 'ydir', 'reverse' )
hold on
h=pcolor(x,y,N');
set(h,'edgecolor','none')
xlabel('age (ka)')
ylabel('horizontal propagation (m)')
title({'horizontal propagation paths probability',run_name})
colorbar

figure(12)
clf
set(gca, 'ydir', 'reverse' )
hold on
h=pcolor(x,y,log10(N'));
set(h,'edgecolor','none')
xlabel('age (ka)')
ylabel('horizontal propagation (m)')
title({'Horizontal propagation paths (log10 probability)',run_name})
colorbar

%Calculate the resulting rates
[Q1toQ2_propagation_rates] = compute_rates(-Q1_propagation_samples, -Q2_propagation_samples, Q1_age_samples, Q2_age_samples);
[Q2toQ3_propagation_rates] = compute_rates(-Q2_propagation_samples, -Q3_propagation_samples, Q2_age_samples, Q3_age_samples);
[Q3toQ4_propagation_rates] = compute_rates(-Q3_propagation_samples, -Q4_propagation_samples, Q3_age_samples, Q4_age_samples);
[Q4toQ5_propagation_rates] = compute_rates(-Q4_propagation_samples, -Q5_propagation_samples, Q4_age_samples, Q5_age_samples);
edges=0:propagation_step_size:max_propagation_rate;
figure(13)
clf
hold on
subplot(4,1,1)
Q1toQ2_propagation_h=histogram(Q1toQ2_propagation_rates, edges);
axis([0 max_propagation_rate 0 inf])
xlabel('Q1 to Q2 propagation rates')
title(run_name)

subplot(4,1,2)
Q2toQ3_propagation_h=histogram(Q2toQ3_propagation_rates, edges);
axis([0 max_propagation_rate 0 inf])
xlabel('Q2 to Q3 propagation rates')

subplot(4,1,3)
Q3toQ4_propagation_h=histogram(Q3toQ4_propagation_rates, edges);
axis([0 max_propagation_rate 0 inf])
xlabel('Q3 to Q4 propagation rates')

subplot(4,1,4)
Q4toQ5_propagation_h=histogram(Q4toQ5_propagation_rates, edges);
axis([0 max_propagation_rate 0 inf])
xlabel('Q4 to Q5 propagation rates (mm/yr)')
ylabel('counts')

figure(17)
label_pos=90;
clf
hold on
subplot(4,1,1)
clean_and_process_pdf(Q1toQ2_propagation_h, sigmas, run_name, label_pos)
axis([0 max_propagation_rate 0 inf])
xlabel('Q1 to Q2 propagation rates')
title(run_name)
subplot(4,1,2)
clean_and_process_pdf(Q2toQ3_propagation_h, sigmas, run_name, label_pos)
axis([0 max_propagation_rate 0 inf])
xlabel('Q2 to Q3 propagation rates')
subplot(4,1,3)
clean_and_process_pdf(Q3toQ4_propagation_h, sigmas, run_name, label_pos)
axis([0 max_propagation_rate 0 inf])
xlabel('Q3 to Q4 propagation rates')
subplot(4,1,4)
clean_and_process_pdf(Q4toQ5_propagation_h, sigmas, run_name, label_pos)
axis([0 max_propagation_rate 0 inf])
xlabel('Q4 to Q5 propagation rates (mm/yr)')
ylabel('probability')

%working on the composite rates and sampling them
All_propagation_rates= [Q1toQ2_propagation_rates; Q2toQ3_propagation_rates; Q3toQ4_propagation_rates; Q4toQ5_propagation_rates]; %Concatenate the rates

nn=length(All_propagation_rates); %Length of the concatenated list of rates
numsamples=n; %Sample it the same number of times as the first sampling
m=ceil(rand(numsamples,1).*nn); %Choose randomly and evenly across the length of the concatenated list of rates
sampled_rates=[];
for i = 1:numsamples
    sampled_propagation_rates(i)=All_propagation_rates(m(i)); %Sample the composite pdf
end

figure(14) % this figure provides a check that we are getting the same composite rate pdf 
subplot(2,1,1)
h=histogram(All_propagation_rates,edges); %Second, make a histogram of the concatenated rates using our bins; THIS IS WHAT WE USE GOING FORWARD
title({'All horizontal propagation rates--histogram of the concatenated rates',run_name})
subplot(2,1,2)
histogram(sampled_propagation_rates,edges); %Third, histogram of the sampled rates;l this is just to check
title('All horizontal propagation rates--histogram of the sampled rates')
xlabel('horizontal propagation rate, mm/yr')
ylabel('counts')
clear sampled_propagation_rates

figure(15)
[most_common_propagation_rate, max_propagation_rate, min_propagation_rate] = clean_and_process_pdf(h, sigmas, run_name, 90);
title({'horizontal propagation rate distribution',run_name})
xlabel('horizontal propagation rate, mm/yr')

clear sample_points_Age_propagation sample_points_propagation

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Saving and printing
if print_save==1
mkdir(run_name)
cd(run_name)

figure(1)
    s=strcat(run_name,'_1_age_pdfs.png');
    print(s,'-dpng')
figure(2)
    s=strcat(run_name,'_2_SU_pdfs.png');
    print(s,'-dpng')
figure(3)
    s=strcat(run_name,'_3_propagation_pdfs.png');
    print(s,'-dpng')
figure(4)
    s=strcat(run_name,'_4_SU_paths.png');
    print(s,'-dpng')
figure(5)
    strcat(run_name,'_5_SU_2D_pdf.png');
    print(s,'-dpng')
figure(6)
    s=strcat(run_name,'_6_log_SU_2D_pdf.png');
    print(s,'-dpng')
figure(16)
    s=strcat(run_name,'_7_SU_rates.png');
    print(s,'-dpng')
figure(8)
    s=strcat(run_name,'_8_SU_concatenated_rates.png');
    print(s,'-dpng')
figure(9)
    s=strcat(run_name,'_9_SU_normalized_rates.png');
    print(s,'-dpng')
figure(10)
    s=strcat(run_name,'_10_horizontal_paths.png');
    print(s,'-dpng')
figure(11)
    s=strcat(run_name,'_11_horizontal_2D_pdf.png');
    print(s,'-dpng')
figure(12)
    s=strcat(run_name,'_12_log_horizontal_2D_pdf.png');
    print(s,'-dpng')
figure(17)
    s=strcat(run_name,'_13_horizontal_rates.png');
    print(s,'-dpng')
figure(14)
    s=strcat(run_name,'_14_horizontal_concatenated_rates.png');
    print(s,'-dpng')
figure(15)
    s=strcat(run_name,'_15_horizontal_normalized_rates.png');
    print(s,'-dpng')
cd ..
end

toc



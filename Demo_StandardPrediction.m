%% Example using the standard prediction
clear; close all; clc

%load the TEX data
load ModelOutput/tex_testdata.mat
%this loads three data structures: 
%castaneda2010
%lopes_santos2010
%shevenell2011
% each contains the following feilds: lat, lon, age, tex86

%select an example timeseries here:

%stuct_tp=castaneda2010;
stuct_tp=lopes_santos2010;
%stuct_tp=shevenell2011;

%% Set the inputs for the prediction code
dats=stuct_tp.tex86;
lon=stuct_tp.lon;
lat=stuct_tp.lat;
prior_std=6;

%select which model to use:
%SST
runname='subT';

%predict:
Output_Struct = bayspar_tex(dats, lon, lat, prior_std, runname);

%and need a timeline:
timeline=stuct_tp.age;
timeline=timeline(:);
%% Now make a plot of the predictions:
figure(1), clf
set(gca, 'fontsize', 16)

fill_trip=[0.7, 0.7, 0.95];
fill([timeline; flipud(timeline)], [Output_Struct.Preds(:,3); flipud(Output_Struct.Preds(:,1))], fill_trip, 'edgecolo', fill_trip), hold on
plot(timeline, Output_Struct.Preds(:,2), 'k-', 'linewidth', 3), hold on
plot(timeline, ones(size(timeline))*Output_Struct.PriorMean, 'r--', 'linewidth', 2)
axis tight
ylabel('Temperature in C')
xlabel('Age')
legend('90% Uncertainty', 'Mean', 'Prior Mean')
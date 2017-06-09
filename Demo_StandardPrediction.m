%% Example using the standard prediction
clear all; close all; clc

%load the TEX data
load ModelOutput/tex_testdata.mat
%this loads three data structures: 
%castaneda2010
%lopes_santos2010
%shevenell2011
% each contains the following feilds: lat, lon, age, tex86

%select an example timeseries here:

%stuct_tp=castaneda2010;
%stuct_tp=lopes_santos2010;
stuct_tp=shevenell2011;

%% Set the inputs for the prediction code
% 
% Idea: import your own data, get it into the form below, and then use the
% code. 
%

dats=stuct_tp.tex86;
lon=stuct_tp.lon;
lat=stuct_tp.lat;
prior_std=6;
min_num=1;
max_dist=500;
%optional inputs
Nsamps=200;
ens_sel=1;

%select which model to use:

%SST
runname=char('SST');

%subT
%runname=char('subT');

%predict:
Output_Struct = bayspar_tex(dats, lon, lat, prior_std, runname, Nsamps, ens_sel);

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

%% and if available, plot 10 of the ensemble members:
if ens_sel==1
    inders=floor(linspace(1, Nsamps, 10));
    for kk=1:1:length(inders)
        plot(timeline, Output_Struct.PredsEns(:,inders(kk)), 'b-'), hold on
    end    
end


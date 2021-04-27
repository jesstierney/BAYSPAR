%% Example using the analog prediction
clear; close all; clc

%%%%%%%%%
%%%%%%%%%
%load example TEX data from the PETM
load ModelOutput/wilsonlake
%structure with fields:
%depth
%tex86
%lat
%lon
%In this case the lat and lon are the paleo estimates,  not the sampling locaitns. 
%relabel for consistency:
wilsonlake.paleolat=wilsonlake.lat;
wilsonlake.paleolon=wilsonlake.lon;
%Set to the standard name to play nice with the code below.
tex_data=wilsonlake; clear wilsonlake

%% Set the inputs for the prediction code:
%Data
dats=tex_data.tex86;
%Prior mean
prior_mean=30;
%Prior standard deviatio
prior_std=20;
%search tolerance
search_tol=std(dats)*2; %set search tolerance to data timeseries std*2
%optional inputs:
%Nsamps=10;
%ens_sel=0; %do not save ensemble


%select which model to use:

%SST
runname='subT';

% Run the prediction:
Output_Struct=bayspar_tex_analog(dats, prior_mean, prior_std, search_tol, runname,1000,1);

%in this case, there are 
N_analogs=length(Output_Struct.AnLocs(:,1))

%% plot the paleo location of the data time series and the analog locations 

figure(1), clf
set(gca, 'fontsize', 14)
set(gcf, 'color', 'w')
worldmap('World')
load coast
plotm(lat, long)
geoshow(tex_data.paleolat, tex_data.paleolon, 'DisplayType', 'point', 'marker', '+', 'linewidth', 4, 'markersize', 12, 'MarkerEdgeColor', 'r'), hold on
% load the original input and plot the modern tex locations:
load(['ModelOutput/Data_Input_SpatAg_', runname])
lt_s=[Data_Input.ObsLocs{:,:}];
tex_lon=lt_s(1:2:end);
tex_lat=lt_s(2:2:end);
geoshow(tex_lat, tex_lon, 'DisplayType', 'point', 'marker', '+', 'linewidth', 1.5, 'markersize', 3, 'MarkerEdgeColor', 'k'), hold on
% plot the outlines of the large gridboxes:
for kk=1:1:length(Output_Struct.AnLocs(:,2))
     lat_vec=[Output_Struct.AnLocs(kk,2)-10, Output_Struct.AnLocs(kk,2)+10, Output_Struct.AnLocs(kk,2)+10, Output_Struct.AnLocs(kk,2)-10, Output_Struct.AnLocs(kk,2)-10];
     lon_vec=[Output_Struct.AnLocs(kk,1)-10, Output_Struct.AnLocs(kk,1)-10, Output_Struct.AnLocs(kk,1)+10, Output_Struct.AnLocs(kk,1)+10, Output_Struct.AnLocs(kk,1)-10];
     geoshow(lat_vec, lon_vec, 'linewidth', 1, 'color', 'k'), hold on 
end
title('Red: Data paleo location. Black cross: modern TEX86 locations.  Black boxes: analog locations.')
%% Now make a plot of the predictions:
%need a depth scale:
depth=tex_data.depth;
depth=depth(:);

figure(2), clf
set(gca, 'fontsize', 16)

fill_trip=[0.7, 0.7, 0.95];
fill([depth; flipud(depth)], [Output_Struct.Preds(:,3); flipud(Output_Struct.Preds(:,1))], fill_trip, 'edgecolo', fill_trip), hold on
plot(depth, Output_Struct.Preds(:,2), 'k-', 'linewidth', 3), hold on
plot(depth, ones(size(depth))*Output_Struct.PriorMean, 'r--', 'linewidth', 2)
axis tight
ylabel('Temperature in C')
xlabel('Depth')
legend('90% Uncertainty', 'Mean', 'Prior Mean')
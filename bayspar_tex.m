function Output_Struct = bayspar_tex(dats, lon, lat, prior_std, runname, varargin)
%
% function Output_Struct = bayspar_tex(dats, lon, lat, prior_std, min_num, max_dist, runname, varargin)
%
% INPUTS:
% dats      - N by 1 vector of TEX86 observations from a single location
% lon       - longitude, from -180 to 180
% lat       - latitude, from -90 to 90
% prior_std - prior standard deviation.
% runname   - specify which model to use. Enter either 'SST' or 'subT'
%
% varargin  - used to set the number of posterior draws of parameters to
% mix across, and whether to save the 5th/50th/95th percentiles of the
% predictions, or the whole ensemble as well.
% - if left empty, then 5000 draws are used and the ensemble of predictions
% is not saved. 
% - if only one argument, it is the number of draws (cannot exceed 15000),
% and the ensemble is not saved. Note that samples are thinned to use the
% full span of the ensmeble. 
% - if two arguments, the first gives the number of draws, while the second
% is an indicator:
% 0: (default) save only the 5th/50th/95th percentiles of the predictions. 
% 1: save the whole ensemble as well. 
%
% in all of the above, the total number of ensemble members is capped by
% the number of draws in the model output. 
%
%prior mean is set as the mean over all instrumental SST observations
%within max_dist (distance is chordal), or the closest min_num points:
%whichever has the larger number of observations.  
%To use the closest K points, set min_num=K and max_dist=0. 
%to use all points within L km, set min_num=1 (ensure at least one) and
%max_dist=L. Results from the paper are with min_num=1; max_dist=500; 
%
%
%Output structure:
%Output_Struct.
%Preds      - Nd by 3 array, where Nd=length(dats), giving the 5/50/95
%percentiles of the predictions. 
%.SiteLoc   - Location of dats, as the inputted [lon, lat];
%.GridLoc 	- Location of the grid centroid closest to dats, used to pull the alpha and beta
%samples. 
%.PriorMean - The prior mean as calculated above. 
%.PriorStd  - THe prior std as input. 
%.PredsEns  - Nd by Nsamps array of predictions. Only include if
%ens_sel==1. 

%% deal with nargin
ng=nargin;
if ng==7
    Nsamps=varargin{1};
    ens_sel=varargin{2};
elseif ng==6
    Nsamps=varargin{1};
    ens_sel=0;
elseif ng==5;
    Nsamps=5000;
    ens_sel=0;
end

%% Load data files needed in the analysis:
load(['ModelOutput/', 'Output_SpatAg_', runname, '/alpha_samples_comp'])
load(['ModelOutput/', 'Output_SpatAg_', runname, '/beta_samples_comp'])
load(['ModelOutput/', 'Output_SpatAg_', runname, '/tau2_samples'])
load(['ModelOutput/', 'Output_SpatAg_', runname, '/Locs_Comp'])
% 

% % and the SST data and locs: THIS IS BRITTLE
if length(strfind(runname, char('SST')))==1
    load ModelOutput/locs_woa_1degree_asvec_SST
    load ModelOutput/st_woa_1degree_asvec_SST
elseif length(strfind(runname, char('subT')))==1
    load ModelOutput/locs_woa_1degree_asvec_subT
    load ModelOutput/st_woa_1degree_asvec_subT
end
%


% grid spacing is hard-coded here:
grid_half_space=10;

% mininum number of grid cells used to calculate modern prior SST
min_num=1;

% maximum distance to search over for modern prior SST
max_dist=500;

%% make sure input is column:
dats=dats(:);
%and that Nsamps is less than 15000:

% trim tau^2 as it may not have had burnin removed:
Ntk=length(alpha_samples_comp(1,:));
tau2_samples=tau2_samples(end-Ntk+1:1:end);
Nsamps=min([Nsamps, length(tau2_samples)]);
%thin the samples to the right number (so as to use the full span of the
%ensemble even if few samples are used.)
ind_s=round(linspace(1, length(tau2_samples), Nsamps));
alpha_samples_comp=alpha_samples_comp(:, ind_s);
beta_samples_comp=beta_samples_comp(:, ind_s);
tau2_samples=tau2_samples(ind_s);


%get the number of obs:
Nd=length(dats);
% get the three pers:
pers3=round([.05, .50, .95]*Nsamps);

% build output structure. 
Output_Struct.Preds=NaN(Nd, 3);
Output_Struct.SiteLoc=[lon, lat];
Output_Struct.GridLoc=NaN(2,1);
Output_Struct.PriorMean=NaN(1,1);
Output_Struct.PriorStd=prior_std;

if ens_sel==1
    Output_Struct.PredsEns=NaN(Nd, Nsamps);
end


%% get the prior mean: constant across the calibration cases, as we just
%average min_num sst obs, on the original 1 degree scale, that are closest to the
%timeseries, or all that are within a max_dist search radius. 

%get the distance to each of the sst gridcells, on the original
%resolution:
dists_prior=EarthChordDistances_2([lon, lat], locs_st_obs);
%order by distance:
[vals_dist_prior, inds_dist_prior]=sort(dists_prior);
%get the number that are below the distance cutoff:
num_below_dist=find(vals_dist_prior<max_dist, 1, 'last' );
   
%if this is larger than the min number, use it to select them:
if num_below_dist>min_num
  prior_mean_val=mean(st_obs_ave_vec(inds_dist_prior(1:1:num_below_dist)));
%otherwise use the smalles min_num:
else
   prior_mean_val=mean(st_obs_ave_vec(inds_dist_prior(1:1:min_num)));
end
% fill in the prior mean value in the output:
Output_Struct.PriorMean=prior_mean_val;

%% figure out the alpha and beta series to draw from.       

% just find the index of the Locs_Comp entry that contains the inputted
% location:

inder_g=find(abs(Locs_Comp(:,1)-lon)<=grid_half_space & abs(Locs_Comp(:,2)-lat)<=grid_half_space);

% Extract the alpha, beta series:
alpha_samples_comp=alpha_samples_comp(inder_g, :)';
beta_samples_comp=beta_samples_comp(inder_g, :)';

% Fill in the outputs:
Output_Struct.GridLoc=Locs_Comp(inder_g, :);


%% set the priors in vector form:
    
Prior_Pars.mu=ones(Nd,1)*prior_mean_val;
Prior_Pars.inv_cov=eye(Nd)*prior_std^(-2); %recall that input is std not var. 
   

%% set the blank array to fill
Preds=NaN(Nd, Nsamps);

tic;
%cycle through to get the predictions in each case:
for JJ=1:1:Nsamps
    %NoNorth
    Preds(:,JJ)= Target_All_Predict(alpha_samples_comp(JJ), beta_samples_comp(JJ), tau2_samples(JJ), dats, Prior_Pars);

    if floor(JJ/1000)-JJ/1000==0
        tt=toc;
        display(['Finished iteration ', num2str(JJ), ' of ', num2str(Nsamps)])
        trem=round((Nsamps-JJ)*tt/1000);
        display(['Approximately ', num2str(trem), ' seconds remaining'])
        tic;
    end
end

%% sort, take the percentiles

Preds_S=sort(Preds, 2);
Output_Struct.Preds=Preds_S(:, pers3);

%% if necessary, save the ensemble as well. 
if ens_sel==1
    Output_Struct.PredsEns=Preds;
end





function Output_Struct = bayspar_tex(dats, lon, lat, prior_std, runname, varargin)
%
% function Output_Struct = bayspar_tex(dats, lon, lat, prior_std, runname, varargin)
%
% INPUTS:
% dats      - N by 1 vector of TEX86 observations from a single location
% lon       - longitude, from -180 to 180
% lat       - latitude, from -90 to 90
% prior_std - prior standard deviation.
% runname   - specify which model to use. Enter either SST or subT
%
% varargin  - used to set the number of posterior draws of parameters to
% mix across, and whether to save the 5th/50th/95th percentiles of the
% predictions, or the whole ensemble as well.
% - if left empty, then 1000 draws are used and the ensemble of predictions
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
elseif ng==5
    Nsamps=1000;
    ens_sel=0;
end

%% Load data files needed in the analysis:
load(['ModelOutput/', 'Output_SpatAg_', runname, '/params_standard'],...
        'alpha_samples_comp','beta_samples_comp','tau2_samples','Locs_Comp');

% and the SST data and locs
if runname=="SST"
    load('ModelOutput/obsSST','locs_st_obs','st_obs_ave_vec');
elseif runname=="subT"
    load('ModelOutput/obssubT','locs_st_obs','st_obs_ave_vec');
end

% grid spacing is hard-coded here:
grid_half_space=10;

% mininum number of grid cells used to calculate modern prior SST
min_num=1;

% maximum distance to search over for modern prior SST
max_dist=500;

%% make sure input is column:
dats=dats(:);

%thin the samples to the right number (so as to use the full span of the
%ensemble even if few samples are used.)
ind_s=round(linspace(1, length(tau2_samples), Nsamps));
alpha_samples_comp=alpha_samples_comp(:, ind_s);
beta_samples_comp=beta_samples_comp(:, ind_s);
tau2_samples=tau2_samples(ind_s);

%get the number of obs:
Nd=length(dats);

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
  prior_mean=mean(st_obs_ave_vec(inds_dist_prior(1:1:num_below_dist)));
%otherwise use the smalles min_num:
else
   prior_mean=mean(st_obs_ave_vec(inds_dist_prior(1:1:min_num)));
end
% fill in the prior mean value in the output:
Output_Struct.PriorMean=prior_mean;

%% figure out the alpha and beta series to draw from.       

% just find the index of the Locs_Comp entry that contains the inputted
% location:

inder_g=find(abs(Locs_Comp(:,1)-lon)<=grid_half_space & abs(Locs_Comp(:,2)-lat)<=grid_half_space);

% Extract the alpha, beta series:
alpha_samples_comp=alpha_samples_comp(inder_g, :);
beta_samples_comp=beta_samples_comp(inder_g, :);

% Fill in the outputs:
Output_Struct.GridLoc=Locs_Comp(inder_g, :);

%% solve
% Prior mean and inverse covariance matrix
    pmu = repmat(ones(Nd, 1) * prior_mean,1,Nsamps);
    pinv_cov = repmat(prior_std,Nd,Nsamps).^-2;
    sigmaS = sqrt(tau2_samples);
    
    % Posterior calculations
    post_mean_num = pinv_cov .* pmu + repmat(sigmaS,Nd,1).^-2 .* repmat(beta_samples_comp,Nd,1) .* (dats - repmat(alpha_samples_comp,Nd,1));
    post_mean_den = pinv_cov + repmat(beta_samples_comp,Nd,1).^2 .* repmat(sigmaS,Nd,1).^-2;
    post_mean = post_mean_num ./ post_mean_den;
    post_sig = sqrt(post_mean_den.^-1);
    Preds = post_mean + randn(Nd,Nsamps).*post_sig;

%% take the percentiles
Output_Struct.Preds = prctile(sort(Preds,2),[5 50 95],2);

%% if asked save the ensemble as well. 
if ens_sel==1
    Output_Struct.PredsEns=Preds;
end

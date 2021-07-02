function Output_Struct = bayspar_tex_analog(dats, prior_mean, prior_std, search_tol, runname, varargin)
%
% function Output_Struct = bayspar_tex_analog(dats, prior_mean, prior_std, search_tol, runname, varargin)
%
% INPUTS:
% dats    - past TEX value or time series. Note that analog selection is
% done by the mean of these, if there are more than one. 
% prior_mean - prior mean on temperature. assumed same for all if dats is a
% vector
% prior_std  - prior std on temperature. assumed same for all. 
% search_tol - tolerance for finding analog locations. comparison is
% between the mean of dats and the mean tex value within each large
% gridcell. 
% runname   - enter 'SST' for sea-surface temperature or 'subT' for subsurface T
%


% varargin  - used to set the number of posterior draws of parameters to
% mix across, and whether to save the 5th/50th/95th percentiles of the
% predictions, or the whole ensemble as well.
% - if left empty, then 1000 draws are used and the ensemble of predictions
% is not saved. 
% - if only one argument, it is the number of draws (cannot exceed 15000),
% and the ensemble is not saved. Note that the first N_samps are used.
% - if two arguments, the first gives the number of draws, while the second
% is an indicator:
% 0: (default) save only the 5th/50th/95th percentiles of the predictions. 
% 1: save the whole ensemble as well. 




%Output structure:
%Output_Struct.
%Preds      - Nd by 3 array, where Nd=length(dats), giving the 5/50/95
%percentiles of the predictions. 
%.AnLocs=- Centroids of the grid cells selected as analogs, as [lon,
%lat] pairs 
%.PriorMean - The prior mean as calculated above. 
%.PriorStd  - The prior std as input. 
%.PredsEns  - Nd by No. analog locatioins by Nsamps array of predictions. Only included if
%ens_sel==1. Note that the second dimension corresponds to the locations
%in .AnLocs. 


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


%% Load the datas
load(['ModelOutput/', 'Output_SpatAg_', runname, '/params_analog'],...
        'alpha_samples','beta_samples','tau2_samples');
load(['ModelOutput/Data_Input_SpatAg_', runname],'Data_Input');
%% make sure input is column:
dats=dats(:);
%thin the samples to the right number (so as to use the full span of the
%ensemble even if few samples are used.)
ind_s=round(linspace(1, length(tau2_samples), Nsamps));
alpha_samples=alpha_samples(:, ind_s);
beta_samples=beta_samples(:, ind_s);
tau2_samples=tau2_samples(ind_s);

%get the number of obs:
Nd=length(dats);

% build output structure. 
Output_Struct.Preds=NaN(Nd, 3);
Output_Struct.AnLocs=[]; %don't yet know size, as it depends on # of analogs. 
Output_Struct.PriorMean=prior_mean;
Output_Struct.PriorStd=prior_std;

if ens_sel==1
    Output_Struct.PredsEns=[]; % don't yet know size. 
end

%% Find the analogs
%cycle through the alpha/beta grid cells, find those that feature mean
%modern TEX obs within the tolerance.

%number of big grids:
N_bg = length(Data_Input.Locs(:,1));

%NEW: calculate mean SSTs across spatial grid cells
spatialMean = NaN(N_bg,1);
for i=1:N_bg
    spatialMean(i)=mean(Data_Input.Obs_Stack(Data_Input.Inds_Stack==i));
end
%identify mean values within the tolerance
inder_g = spatialMean >= (mean(dats)-search_tol) & spatialMean <= (mean(dats)+search_tol);
if sum(inder_g)==0
    error('Your search tolerance is too narrow')
else
end
alpha_samples=alpha_samples(inder_g, :);
beta_samples=beta_samples(inder_g, :);
%tile tau2 to match
tau2_samples=repmat(tau2_samples,1,size(alpha_samples,1));
%reshape
alpha_samples=reshape(alpha_samples,1,size(alpha_samples,1)*size(alpha_samples,2));
beta_samples=reshape(beta_samples,1,size(beta_samples,1)*size(beta_samples,2));
%save the analog locations
Output_Struct.AnLocs=Data_Input.Locs(inder_g,:);

%% solve
% Prior mean and inverse covariance matrix
    pmu = repmat(ones(Nd, 1) * prior_mean,1,size(alpha_samples,2));
    pinv_cov = repmat(prior_std,Nd,size(alpha_samples,2)).^-2;
    sigmaS = sqrt(tau2_samples);
    
    % Posterior calculations
    post_mean_num = pinv_cov .* pmu + repmat(sigmaS,Nd,1).^-2 .* repmat(beta_samples,Nd,1) .* (dats - repmat(alpha_samples,Nd,1));
    post_mean_den = pinv_cov + repmat(beta_samples,Nd,1).^2 .* repmat(sigmaS,Nd,1).^-2;
    post_mean = post_mean_num ./ post_mean_den;
    post_sig = sqrt(post_mean_den.^-1);
    Preds = post_mean + randn(Nd,size(alpha_samples,2)).*post_sig;

%% get the three pers:
Output_Struct.Preds = prctile(sort(Preds,2),[5 50 95],2);

%% if requested, save the ensemble as well:
if ens_sel==1
    Output_Struct.PredsEns=Preds;
end

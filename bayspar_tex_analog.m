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
% - if left empty, then 5000 draws are used and the ensemble of predictions
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
elseif ng==5;
    Nsamps=5000;
    ens_sel=0;
end


%% Load the datas
load(['ModelOutput/', 'Output_SpatAg_', runname, '/alpha_samples'])
load(['ModelOutput/', 'Output_SpatAg_', runname, '/beta_samples'])
load(['ModelOutput/', 'Output_SpatAg_', runname, '/tau2_samples'])
load(['ModelOutput/Data_Input_SpatAg_', runname])

alpha_samples=[alpha_samples.field];
beta_samples=[beta_samples.field];
%% make sure input is column:
dats=dats(:);
%and that Nsamps is less than 15000:
Nsamps=min([Nsamps, length(tau2_samples)]);
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

%% get the prior in the right form:
Prior_Pars.mu=ones(Nd,1)*prior_mean;
Prior_Pars.inv_cov=eye(Nd)*prior_std^(-2);

%% Find the analogs

%cycle through the alpha/beta grid cells, find those that feature mean
%modern TEX obs within the tolerance: 

%number of big grids:
N_bg=length(Data_Input.Locs(:,1));
inder_g=[];

%%%% to here. 
for kk=1:1:N_bg
    %find the Tex obs corresponding to this index location, using the
    %stacked obs:
    inds_temp=find(Data_Input.Inds_Stack==kk);
    vals=Data_Input.Obs_Stack(inds_temp);

    %if the mean of vals in the big gridis within tolerance, add it to inder_g
    if mean(vals)>=(mean(dats)-search_tol) && mean(vals)<=(mean(dats)+search_tol)
        inder_g=[inder_g; kk];
    end    
    
end
N_Locs_g=length(inder_g);
Output_Struct.AnLocs=Data_Input.Locs(inder_g,:);

%cycle through to get the predictions
preds_T=NaN(size(dats,1),N_Locs_g, Nsamps);

tic;
ct=0;
for jj=1:1:N_Locs_g
    for kk=1:1:Nsamps
        preds_T(:,jj,kk) = Target_All_Predict(alpha_samples(inder_g(jj),kk), beta_samples(inder_g(jj),kk), tau2_samples(kk), dats, Prior_Pars);
    
        %deal with timing. 
        ct=ct+1;
        if floor(kk/1000)-kk/1000==0
            tt=toc;
            display(['Finished iteration ', num2str(kk), ' of ', num2str(Nsamps), ' for spatial analog ', num2str(jj), ' of ', num2str(N_Locs_g)])
            trem=round((Nsamps*N_Locs_g-ct)*tt/1000);
            display(['Approximately ', num2str(trem), ' seconds remaining'])
            tic;
        end
    end
    
end

%% get the three pers:
% recall that there are Nsamps predictions for each of the N_Locs_g
% analog locations. 
pers3=round([.05, .50, .95]*Nsamps*N_Locs_g);

%stack the predictions across the analog locations to enable sorting across them: 
preds_T_stack=reshape(preds_T, Nd, Nsamps*N_Locs_g);

%sort
preds_T_stack_S=sort(preds_T_stack,2);

%get the three percentiles:
Output_Struct.Preds=preds_T_stack_S(:,pers3);

%% if necessary, save the ensemble as well:
if ens_sel==1
    Output_Struct.PredsEns=preds_T;
end





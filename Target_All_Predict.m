function Target_timeseries_pred = Target_All_Predict(alpha_now, beta_now, tau2_now, Proxy_TS, prior_pars)
%
%
%inputs:
% alpha_now -> scalar
% beta_now -> scalar
% tau2_now -> current value of the residual variance; scalar 
% Proxy_TS -> time series of the proxy. assume no temporal strucutre for
% now, as timing is not equal. 
% prior_pars -> prior_pars.mu: vector of prior means for each element of
% the time series; prior_pars.inv_cov: INVERSE of the prior covariance
% matrix for the timeseries. 
% time 
% 
% output -> sample of target time series vector conditional on the rest. 

%% TEST
%JJ=1;
%kk=1;
%alpha_now=alphas_Y(JJ);
%beta_now=betas_Y(JJ);
%tau2_now=tau2_YesRS(JJ); 
%Proxy_TS=tex_timeseries(kk).Tex86;
%prior_pars=Prior_Pars_Y;

%% 
N_Ts=length(Proxy_TS);

% get the inverse posterior covariance matrix:
inv_post_cov=prior_pars.inv_cov + beta_now^2/tau2_now*eye(N_Ts);

%USE CHOLESKY HERE TO SPEED THINGS UP!
%first get the post_cov, which requires an inverse.
%use linsolve with options
opts.SYM=true; opts.POSDEF=true;
post_cov=linsolve(inv_post_cov, eye(N_Ts), opts);

%get the square root using the cholesky factor:
sqrt_post_cov=chol(post_cov)';

%now get the first factor for the mean:
Mean_first_factor = prior_pars.inv_cov *prior_pars.mu + (1/tau2_now)* beta_now *(Proxy_TS - alpha_now);


Mean_full=post_cov*Mean_first_factor;


%
Target_timeseries_pred = Mean_full + sqrt_post_cov*randn(N_Ts,1);



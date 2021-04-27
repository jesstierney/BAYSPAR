# BAYSPAR

BAYSPAR is a Bayesian, spatially-varying model for the TEX<sub>86</sub> paleotemperature proxy. When using this model, please cite the original paper describing the methodology:

Tierney, J.E. & Tingley, M.P. (2014) A Bayesian, spatially-varying calibration model for the TEX<sub>86</sub> proxy. Geochimica et Cosmochimica Acta, 127, 83-106. https://doi.org/10.1016/j.gca.2013.11.026.

This package includes functions to both forward-model TEX<sub>86</sub> from ocean temperatures and predict ocean temperatures from TEX<sub>86</sub>. These files make use of the regression parameters fit to the latest TEX<sub>86</sub> core top dataset (n = 1095, Tierney and Tingley, 2015, Scientific Data, https://doi.org/10.1038/sdata.2015.29). The regressions exclude TEX<sub>86</sub> data north of 70˚N due to the poor predictability of TEX<sub>86</sub> in the Arctic region (see discussion in TT14).

Two models are available: The SST model and the subT model. The latter uses a gamma-weighted average of temperatures between 0-200 meters water depth (see TT15). In addition, two modes of prediction are available, "standard" for Late Quaternary applications and "analogue" for deep time applications (see TT14 for details).

Below is a quick guide to using these functions.

## Predict TEX<sub>86</sub> from T (`TEX_forward.m`)

This function forward models TEX<sub>86</sub> from either SSTs or subTs. Usage is:

`tex = TEX_forward(lat, lon, t, varargin)`

INPUTS:
lon		- longitude of core site, decimal degrees from -180 to 180\
lat		- latitude of core site, decimal degrees from -90 to 90\
t		- a vector (N x 1) of SST or subT\
varargin        - if empty, the function uses the 'SST' model in standard mode\
		- use one argument to specify 'SST' vs 'subT'\
		- use three arguments to specify 'SST' vs 'subT' and use analog mode. The second argument is a string specifying "standard" or "analog" and the third argument is the search tolerance in temperature units. For example,

`tex = TEX_forward(lat, lon, t, 'SST', "analog", 10)`

would forward-model TEX<sub>86</sub> in analog mode with a search tolerance of 10 degrees C around the input t values.

OUTPUT: A N x 1000 ensemble of possible temperature values.

NOTE: tex_forward can now handle multiple lat/lon locations!

## Predict T from TEX<sub>86</sub> - Standard mode (`bayspar_tex.m`)

The assumption in this case is that oceanographic conditions are sufficiently similar to today so that it is reasonable to use the spatial distribution of regression parameters as fitted from the core tops.

The function is `bayspar_tex.m` and usage is demonstrated in `Demo_StandardPrediction.m`.

The inputs to and outputs from bayspar_tex.m are:

`Output_Struct = bayspar_tex(dats, lon, lat, prior_std, runname, varargin)`

INPUTS:

dats		- N by 1 vector of TEX<sub>86</sub> observations from a single location\
lon		- longitude of core site, decimal degrees from -180 to 180\
lat		- latitude of core site, decimal degrees from -90 to 90\
prior_std 	- prior standard deviation in degrees Celsius (scalar). Try a value like 10 to start. Don't make this too narrow.\
runname - defines which model to use. Enter 'SST' or 'subT'.\
varargin  	- optional arguments used to set the number of posterior draws of parameters to mix across, and whether to save only the 5th/50th/95th percentiles of the predictions or the whole ensemble as well.\
	- if left empty, then 1000 draws are used and the ensemble of predictions is not saved.\
	- if only one argument, it is the number of draws (cannot exceed 15000), and the ensemble is not saved.\
	- if two arguments, the first gives the number of draws, while the second is an indicator:\
		0: (default) save only the 5th/50th/95th percentiles of the predictions.\
		1: save the whole ensemble as well.

OUTPUT:

A structure with the following fields:
.Preds      		- N x 3 array, where N=length(dats), giving the 5/50/95 percentiles of the predictions in degrees C.\
.SiteLoc   		- Location of dats, as the inputted [lon, lat].\
.GridLoc   		- Location of the grid centroid closest to dats, used to pull the alpha and beta samples.\
.PriorMean 	- The prior mean (see below).\
.PriorStd  		- The prior std - same value as specified in the inputs.\
.PredsEns  	- N x 1000 array of predictions. Only included if the second optional input is set to 1.

The prior mean is set as the mean over all instrumental SST observations within 500 km of the observations, or the closest data point: whichever has the LARGER number of observations.

The alpha and beta values used in predictions are selected based on the location of the inputted TEX<sub>86</sub> series, and the possibility arises that the sample location will not be within one of the 20 by 20 degree grid boxes that contain core top calibration data. As the alpha and beta fields are modeled as a Gaussian spatial process, we can draw samples at any location conditional on the scalar parameters (mu, phi, sigma) and the values of alpha and beta at the centroids of the grids that contain core top calibration observations. We do so for each of the retained 15000 samples from the posterior distribution of alpha, beta, and the scalar parameters so as to build up the conditional distributions of alpha and beta and the centroids of all 20 by 20 grid boxes. This step is performed off-line and results simply looked up. Note that the distributions of alpha and beta at interpolated locations are necessarily wider, as are the predictive uncertainties.

## Predict T from TEX<sub>86</sub> - Analogue mode (`bayspar_tex_analog.m`)

The assumption in this case is that oceanographic conditions are sufficiently dissimilar to the modern configuration (e.g. due to tectonic changes) so as to preclude the use of the regression parameters that correspond to the current location of the study site (see Tierney and Tingley, 2014, GCA for details). The regression parameters used to make predictions are instead mixed across locations deemed analogous to the inputted TEX<sub>86</sub> values. To do this, the function finds all calibration grid boxes with average TEX<sub>86</sub> values within a user-set tolerance of the mean of inputted TEX<sub>86</sub> values.

The function is `bayspar_tex_analog.m` and usage is demonstrated in `Demo_AnalogPrediction.m`.

The inputs to and outputs from bayspar_tex_analog.m are:

`Output_Struct = bayspar_tex_analog(dats, prior_mean, prior_std, search_tol, runname, varargin)`

INPUTS:
dats			- N x 1 vector of TEX<sub>86</sub> values or time series. Note that analogue selection is performed using the mean of these values.\
prior_mean 	- prior mean on temperature. (scalar, in degrees Celsius)\
prior_std  		- prior std on temperature. (scalar, in degrees Celsuis). Don't make this too narrow.\
search_tol 	- tolerance for finding analog locations in TEX<sub>86</sub> units. Comparison is between the mean of dats and the mean TEX<sub>86</sub> value within each 20 by 20 grid box. In the demo code, search_tol is specified in terms of standard deviations of the input TEX<sub>86</sub> series.\
runname           - defines which model to use. Enter 'SST' or 'subT'.\
varargin  		- used to set the number of posterior draws of parameters to mix across, and whether to save only the 5th/50th/95th percentiles of the predictions or the whole ensemble as well.\
	- if left empty, then 1000 draws are used FOR EACH ANALOG LOCATION and the ensemble of predictions is not saved. \ 
	- if only one argument, it is the number of draws (cannot exceed 15000) FOR EACH ANALOG LOCATION, and the ensemble is not saved. Note that samples are thinned to use the full span of the ensemble.\
	- if two arguments, the first gives the number of draws, while the second is an indicator:\
		0: (default) save only the 5th/50th/95th percentiles of the predictions. \
		1: save the whole ensemble as well.

OUTPUT:

A structure with the following fields:

Preds      		- Nd by 3 array, where Nd=length(dats), giving the 5/50/95 percentiles of the predictions. \
.AnLocs		- Centroids of the grid cells selected as analogs, as [lon, lat] pairs \
.PriorMean 	- The prior mean as calculated above. \
.PriorStd  		- The prior std as input. \
.PredsEns  	- Nd by No. analog locatioins by Nsamps array of predictions. Only included if ens_sel==1. Note that the second dimension corresponds to the location in .AnLocs.

NOTE:

All necessary model outputs, as well as some sample TEX<sub>86</sub> series used by the demo scripts, are in the ModelOutput folder, and need not be touched.

# EXT-Tseries-ANITI
This repository has two branches: 1)the main branch which deals with udpated works on the simulation of extreme time series for wind speed and (2) the multivariate branch which deals with the multivariate adaptation.
You are here in the multivariate branch.

2) Simulation of extreme multivariate time series accounting for temporal dependence (MV case)

For better modeling surge-induced coastal flooding, we analyse extreme multivariate time series and build a simulator of extreme time series. 
Our main contributions are the following:

-(a) Accounting for temporal dependence and short-tailed behavior. Use of regularly varying functions, extreme value model and polar decomposition.

-(b) Polar decomposition of extreme multivariate time series with the modeling of polar coordinates. First, we model the distribution of normalised multivariate time series by approaching the distribution of their coordinates in a finite-dimensional basis. 

Then, we model the distribution of the radial elements with a CEV model, which enables to take into account the asymptotic independence between some pairs of the radial vector.

-(c) Tunable aspects with the choice of the compound risk function $g$ and the level of extremeness. 

-(d) Several methods proposed to validate the simulation method, using extreme value theory (univariate and bivariate), PCA decomposition and classification two-samples test. 


To run the codes, first run [Detrending_Season](./Detrending_Season.Rmd) for detrending the time series then run [Whitening](./Whitening.Rmd) for whitening them.
Then, run [Conversion_MV_data](./Conversion_MV_data.Rmd) to convert in the convenient all the marginal of the multivariate time series. This normalisation enables to simulate new extreme residuals in [Simul_ext_Residuals](./Simul_ext_Residuals.Rmd).
Finally, using the residuals displayed in folder (./residuals_MV), run  and [Simul_Ext_Forcing](./Simul_Ext_Forcing.Rmd) for simulating respectively extreme residuals and extreme time series of interest. consistency results are displayed in [Consistency_analysis](./Consistency_analysis.Rmd). 

The workflows displayed in (./Worflows) summarises the simulation method but also detail the effect of the choice of function $g$. The [preprocessing steps](Workflows/Preprocessing_steps.png) displays the whitening and the marginal transformations. These preliminar steps enable us to apply a [polar decomposition](Workflows/polar_decomp.png), which is the backbone of our probabilistic model. 
Once the distribution of the polar coordinates have been estimated, we simulate new extreme multivariate time series by first sampling from the [models](Workflows/simulation_polar_step.png), and, then, apply the [reverse transformations](Workflows/simulation_original.png). 

All details of the methods are provided in [Gorse et al. (2025)](https://arxiv.org/abs/2508.13687). The methods are applied to [surge](./Data/Data_Surge.csv), [significant wave height](./Data/Data_Hs.csv) and [wind speed data](./Data/Data_U.csv) of Gavres site (French Atlantic coast). If you use these datasets, please refer to [Idier et al. (2020)](https://link.springer.com/article/10.1007/s11069-020-03882-4).

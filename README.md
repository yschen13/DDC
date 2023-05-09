# Dynamical differential covariance (DDC)
Reference: https://www.pnas.org/doi/abs/10.1073/pnas.2117234119  <br />
All scripts were carefully annotated. Below is a quick summary.
## Connectivity estimation:
  * estimators.m: Common estimators including covariance <x,x>, Precision matrix, nonlinear averaging matrix <R(x),x>, dCov computed by symmetric difference quotient 
  * dCov_numerical.m: dCov calculated by different numerical estimations of the time derivative
  * derivative_123.m: Supporting file for dCov_numerical.m
  * DDC calculation: matrix product of dCov and  <x,x>^{-1} or <R(x),x>^{-1}
  * L2 regularized DDC calculation: dCov_linear_Reg.m
## PSD Bootstrapping (evaluate significance of connections)
 * PSDBootstrap_example.m: generate null timeseries under the null hypothesis that each time series is independent; calculate the null FC matrices
 * PSDBootstrap_Getmdl.m: fit the best AR model to the order of q for each time series
## Network simulations:
  * Linear_simulation.m: simulation of linear stochastic systems
  * Nonlinear_simulation.m: simulation of nonlinear stochastic systems
  * LIF_network_YC.m: Leaky-integrate-and-Fire network simulation
  * simulate_reduced_wong_wang.ipynb: Reduced Wong-Wang simulation
## Performance evaluation: 
* FC_GT_ROC.m: classification sensitivity for LIF network recovery
* c_sensitivity_YC.m: c-sensitivity calculation
## HCP dataset
* rs-fMRI time series data available through https://www.humanconnectome.org
* dMRI dataset available through https://www.eneuro.org/content/8/1/ENEURO.0416-20.2020/tab-article-info
* ATLAS_HCP_dMRIparcel_to_IC100node_connectivity.m: convert parcel-level dMRI strength to IC-level dMRI strength
## Demo (MATLAB syntax): 
For HCP ICA preprocessed time series: 
     
     TR = 0.72; 
     [T, N] = size(V); % number of timepoints x number of nodes
     V_obs = zscore(V)
     [dCov1, dCov2,~,~] = dCov_numerical(V_obs,TR);
     [Cov,Precision,B,~] = estimators(V_obs,prctile(V_obs(:),50),TR);
     Delta_L = dCov2*Precision
     Delta_ReLU = dCov2*B


# Dynamical differential covariance (DDC)
All scripts were carefully annotated. Below is a quick summary.
## Connectivity estimation:
  * estimators.m: Common estimators including covariance, precision matrix <x,x>, nonlinear averaging matrix <R(x),x>, first-order derivative computed by symmetric difference quotient 
  * dCov_numerical.m: numerical estimations of the time derivative
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
* ATLAS_HCP_dMRIparcel_to_IC100node_connectivity.m: convert parcel level dMRI strength to IC-level dMRI strength

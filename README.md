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

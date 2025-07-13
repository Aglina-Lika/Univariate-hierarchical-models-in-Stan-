# Univariate-hierarchical-models-in-Stan-
This repository contains a simulation for fitting Hierarchical models with different parametrizations and variable transformations. The two parametrizations are the hierarchical centering parametrization (HCP) and hierarchical non centering-parametrization.  We fit the models with raw, centred and standardized variables. 

In the folder Stan model you can find the Stan codes for fitting the models with the use of the commant rstan in R. There are aslo 3 extra files:
* The extractFrames_new  is a function that everytime provides us the tables that are required for a specific parametrization and variable transformation.
* The simExp_linear contains the code for the simulation of 100 data set for every combination of the parameters, parametrizations and variable transformations.
* The box plot for every parametrization and variable transformation is creating the box plots of the effective sample size devided by the convergence time (EES(t)/time) for every parameter and for every parametrization and variable transformation.  

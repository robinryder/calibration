
%README file - this directory contains the matlab scripts we used to produce 
%the results in section 5 of the paper

%Jeong Eun Lee, Geoff K. Nicholls, Robin J. Ryder (2018)
%"Calibration procedures for approximate Bayesian credible sets"
%https://arxiv.org/abs/1810.06433

%The two top-level scripts are runcal.m and recal.m
%runcal.m takes a few hours to run (on a laptop, with no parallelisation)

%The file ExactFreApproxCylIce.mat contains the output of the run used 
%to generate the figures in the paper: load ExactFreApproxCylIce.mat
%and run the output-analysis section of runcal.m, and recal.m

%Some R-code is included in the comments at the end of runcal.m -
%this was used to carry out the logistic-GAM regression in section 5.
%The R-code loads coverage data from isif.csv and fits the GAM.
%The coverage data in isif.csv comes from ExactFreApproxCylIce.mat.

%The file containing the function evaluating the Ising partition function
%for periodic boundary conditions is logZ.m

%Geoff Nicholls 2/1/19
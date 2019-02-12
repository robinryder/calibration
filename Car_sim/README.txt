README file - this directory contains the R scripts we used to produce the results in section 6 of the paper

Jeong Eun Lee, Geoff K. Nicholls, Robin J. Ryder (2018)
"Calibration procedures for approximate Bayesian credible sets"
https://arxiv.org/abs/1810.06433

a) carsim_result.RData contains 40500 particles and their KS distances and coverage indicators.
   Figure 6 is based on this result.

b) codes2.r contains all source codes.

c) sim.r contains the three parts and each part is separated by %%%%%. 
	First part : Data preparation and prior setup.
	Second part : Coverage estimation using Importance sampling (Algorithm 3) and exact coverage estimation (Algorithm 1). 
	Third part : Read 'carsim_result.RData' to load 40500 particle simulations and replicate Figure 6.


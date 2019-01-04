function [ theta_min, theta_max ] = getCI( p, theta, level )
%[theta_min,theta_max]=getCI(p,theta,level) - compute equal tailed credible
%set using approximate posterior density p (step-function approx of pi-tilde)
%INPUT
%   p - vector of density values at theta-values in theta vector
%   theta - vector of values of Ising smoothing parameter/Inverse Temp
%   level - the coverage level of the credible set we are computing
%OUTPUT
%   theta_min, theta_max - bottom and top ends of credible interval

%%
cdf=cumsum(p)*(theta(2)-theta(1));

%%
a2=(1-level)/2; %equal tailed
iL=find(cdf>a2,1,'first');
theta_min=theta(iL);
iU=find(cdf<1-a2,1,'last');
theta_max=theta(iU);
%%

end


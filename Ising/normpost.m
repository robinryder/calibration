function [np]=normpost(data,theta,m,n,nume)
%[np]=normpost(data,theta,m,n,nume) - compute normalised posterior density at pi(theta|data)
%INPUT
%   data - #x value
%   theta - inverse temperature/smoothing parameter values - a linspace
%     vector of values - spanning the whole space of theta
%   m,n - side lengths
%   nume - number of edges in nbrhood structure
%OUTPUT 
%   np - a vector of posterior density values (in histogram style approx) 

%%
lp=loglkd(data,theta,m,n,nume)+logprior(theta);

%%
up=exp(lp-max(lp));
%%
np=up/(sum(up)*(theta(2)-theta(1)));

end


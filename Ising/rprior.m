function [ theta ] = rprior( n )
%rprior() - simulate prior
%   prior is U(0,2), n = number of realisations

theta=2*rand([n,1]);

end


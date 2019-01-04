function [ sdat,LL,UL,c,d,e ] = Calibrate(S,phi,L,LSS,M,nbrsEXACT,th,numeC,postDT,postDTCDF)
%Calibrate() - simulate synthetic data y'=sdat for each phi and then for
%each y' compute the approx post. Calculate credible interval and record
%if covers phi. Also record the KS dist between post at y and post at y'.
%[sdat,LL,UL,c,d,e]=Calibrate(S,phi,L,LSS,M,nbrsEXACT,th,numeC,postDT,postDTCDF)
%INPUT
%   S - number of simulated y-primes
%   phi - vector of 'true' values
%   L,LSS - Ising MCMC run length and subsampling
%   M - Lattice size
%   nbrsEXACT - cell array giving nbrs of each pixel.
%   th - vector of theta values where posteriors are evaluated
%   numeC - number of edges in the lattice
%   postDT - (usually approx) posterior at data (pdf)
%   postCDF - (usually approx) posterior at data (CDF)
%OUTPUT
%   sdat - S sampled #x values (correspond to s(y^{(i)}) in paper)
%       these are sampled at the phi-values in the vector phi
%   LL - lower limits for credible intervals - LL(i) is the LL for the
%       credible interval computed using the approx post
%   UU - upper limit as LL
%   c - vector of S 0/1 values indicator for i'th CI covering phi(i)
%   d - vector of S KS-distances
%   e - vector of homebrew distances

c=ones(S,1);
sdat=c;
d=c;
e=c;
UL=c;
LL=c;

%simulating from smallest phi to largest so initialise for phi=0.
x=reshape(randsample(0:1,M^2,true),[M,M]);

for (s=1:S)
    s
    %simulate synthetic data y'~p(y'|phi) with y'=Ds here
    [Ds,x,junk]=ising(phi(s),L,LSS,M,nbrsEXACT,x); %done using open bconds
    %figure(4); imagesc(x); colormap(gray); axis square; drawnow;
    sdat(s)=Ds;
    %normalised posterior ie pi-tilde(theta|y') - analytical not sampled
    post=normpost(sdat(s),th,M,M,numeC);
    postCDF=cumsum(post)*(th(2)-th(1));
    %upper and lower limits of credible interval
    [LL(s),UL(s)]=getCI(post,th,0.95);
    %did we cover? 0/1
    c(s)=( (phi(s)<UL(s)) && (phi(s)>LL(s)) );
    %the max-distance between the CDF's (ie KS-stat)
    d(s)=max(abs(postDTCDF-postCDF));
    %the max distance between the pdf's
    e(s)=max(abs((postDT)-(post)));
end

end


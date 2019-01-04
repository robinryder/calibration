function [ llk ] = loglkd(data,theta,m,n,nume)
%loglkd() - evaluate the ising log likelihood for param theta
%  exact for periodic boundary conditions 
%  data=#x, theta=inv temp, mxn lattice, nume number of edges

%%
lt=length(theta);
llk=zeros(1,lt);
for i=1:lt
  llk(i)=-theta(i)*(data-nume/2)-logZ(m,n,theta(i)/2);
end
%%

end


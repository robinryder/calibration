function [ phi ] = rpost( n, cdf, theta )
%rpost() - n iid samples from the posterior
%   n=#samples, p=pdf at the values in the vector theta

r=rand(n,1);
phi=zeros(n,1);
for k=1:n
    phi(k)=theta(find(cdf>r(k),1,'first'));
end

end


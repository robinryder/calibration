function lZf=logZ(n,m,K)


%% lZf=logZ(n,m,K)
% log of the Ising partition function with K=J/K_BT
% on m x n square lattic with periodic boundary conditions
% Uses result from P. Beale, Phys. Rev. 76:78-81 (1996)
% due originally to B. Kaufman, Phys Rev, 76:1232-1243 (1949)
%
%% sample call
% n=512;m=n; c=m*n;
% s=1000;
% KC=log(1+sqrt(2))/2; %critical value
% b=linspace(KC-0.1,KC+0.1,s); lZf=zeros(1,s);
% for i=1:s
%     lZf(i)=logZ(m,n,b(i));
% end
% %% plot the heat capacity K^2 d^log(Z)/dK^2 
% plot(b(3:end)/KC,b(3:end).^2.*diff(diff(lZf/c))/(b(2)-b(1))^2);

%%
if (K==0) 
    lZf=m*n*log(2);
else
    c=n*m;
    k=0:(n-1);
    g=[log(exp(2*K)*tanh(K)),acosh(cosh(2*K)^2/sinh(2*K)-cos(pi*(1:(2*n))/n))];
    Y1=sum(log(2*cosh(m*g(2*k+2)/2)*sinh(2*K)^(m/2)));
    Y2=sum(log(2*sinh(m*g(2*k+2)/2)*sinh(2*K)^(m/2)));
    Y3=sum(log(2*cosh(m*g(2*k+1)/2)*sinh(2*K)^(m/2)));
    Y4=sum(log(2*sinh(m*g(2*k+1)/2)*sinh(2*K)^(m/2)));
    lZf=real( (c/2-1)*log(2)+Y1+log(1+exp(Y2-Y1)+exp(Y3-Y1)+exp(Y4-Y1)) );
end

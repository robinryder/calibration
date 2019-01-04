%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Coverage tests for Ising model smoothing parameter
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Common material - neighbourhood structure

REALDATA=true;
if REALDATA
    %the data for this example comes from 
    %Bornn, L., Jacob, P., Moral, P., and Doucet, A. (2013). 
    %"An Adaptive Interacting Wang-Landau Algorithm for Automatic 
    %Density Exploration" 
    %Journal of Computational and Graphical Statistics, 22(3):749-773. 
    xT=csvread('IceFloe.csv');
    M=size(xT,1); %this will be 40 x 40
else
    M=32;         %square lattice side, rectangles would need a few edits
end

%Get nbrs on square lattice with cylindrical boundary conditions
nbrsC=GetNbrs(M,M,'cylindrical');
numeC=length([nbrsC{:}])/2; %num edges - need this for lkd on theta-side

%Get nbrs on square lattice with free boundary conditions
nbrsF=GetNbrs(M,M,'free');

%exact is the phi-side, approx is the theta-side main example uses 
%free bnd conds on the phi-side and cylindrical on the theta side 
nbrsEXACT=nbrsF;
%nbrsEXACT=nbrsC; % if use this then no approximation phi/thetaside same

%% Common material - summary statistics for data, and approx posterior

Lsynth=M^2*10000;     %number of MCMC steps (one pixel flip per step)
LsynthSS=Lsynth;      %if subsample

if REALDATA
    DT=hashX(xT,nbrsEXACT);
else
    %simulate synthetic "true" data y=DT - y~p(y|thetaT) - we assume this
    %simulation is exact or at least the error is << error from bnd choice
    thetaT=0.3;           %true theta used to simulate y
    
    %simulate Ising state on MxM lattice
    [DT,xT,nbrsT]=ising(thetaT,Lsynth,LsynthSS,M,nbrsEXACT);
end
figure(1); imagesc(xT);colormap(gray);axis square;drawnow;

%now compute the theta-dbn at the data, theta|Y=DT
nt=1000;
th=linspace(0,2,nt);

%this is an accurate approximation to the approximate posterior - simply
%evaluated at nt theta-values and normalised
postDT=normpost(DT,th,M,M,numeC); 
postDTCDF=cumsum(postDT)*(th(2)-th(1));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Logistic Regression Approach 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%LR - Generate y', phi and c vectors
L=Lsynth/10;     %run length for each simulated data y'=Ds
LSS=L;           %sub-run subsample if required (not here)

S=1000;          %M from the paper - number of (phi,y') pairs to simulate
[phi]=sort(rprior(S));     %sorting helps a bit with the Ising simulation

[sdat,LL,UL,c,d,e]=Calibrate(S,phi,L,LSS,M,nbrsEXACT,th,numeC,postDT,postDTCDF);

%LR - Use LR/GAM to estimate coverage at data

%done in R - see note below.


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% IS Approach
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%
%IS - Generate y', phi and c vectors

Li=Lsynth/10;      %run length for each simulated data y'=Ds
LSSi=Li;           %sub-run subsample if required (not here)

Si=1000;           %number of phi's and synthetic data sims y'
[phii]=sort(rpost(Si,postDTCDF,th)); %sorting helps a bit with the Ising simulation
lkdi=loglkd(DT,phii,M,M,numeC);

[sdati,LLi,ULi,ci,di,ei]=Calibrate(Si,phii,Li,LSSi,M,nbrsEXACT,th,numeC,postDT,postDTCDF);

%%
%IS - compute IS estimate for coverage (output analysis)

fl=exp(-lkdi+max(lkdi));
rho=0.5;
fl=fl(di<=rho);
w=fl./sum(fl);
ESS=1/sum(w.^2) %ESS
ESS_sig=sum(w.^2)^2/sum(w.^4) %ESS for the variance estimate
cc=ci(di<=rho); ncc=length(cc);
cTis=w*cc   %estimated IS cover
cTisD=sqrt(sum(w.^2*(cc-cTis).^2))  %estimated standard error
plot(di(di<=rho),log(w),'*',[0.5,0.5],[min(log(w)),0]);

%% Notes on GAM regression in R
%In the end I wrote out the c,sdat,d table and did a regression using R and
%gam() - I dont think GAM's are available in MatLab?
%
% In MatLab ...
%outt=table(c,sdat,d)
%writetable(outt,'isif2.csv')
%
% And now in R...
% > cov.dat=read.csv("isif.csv")
% > 
% > library(mgcv)
% > summary(cov.gam<-gam(c~s(sdat),family=binomial,data=cov.dat))
% 
% Family: binomial 
% Link function: logit 
% 
% Formula:
% c ~ s(sdat)
% 
% Parametric coefficients:
%             Estimate Std. Error z value Pr(>|z|)    
% (Intercept)   1.9952     0.1026   19.45   <2e-16 ***
% ---
% Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
% 
% Approximate significance of smooth terms:
%           edf Ref.df Chi.sq p-value    
% s(sdat) 6.225  7.356  45.15   2e-07 ***
% ---
% Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
% 
% R-sq.(adj) =  0.0474   Deviance explained = 6.46%
% UBRE = -0.24864  Scale est. = 1         n = 1000
% > lp=predict(cov.gam)
% > plot(cov.dat$sdat,exp(lp)/(1+exp(lp)),xlab="Sufficient Statistic g(y';E_F)",ylab="Estimated Coverage")
% > et=predict(cov.gam,data.frame(sdat=503,d=0)) #the true value of the summary stat is 503
% > exp(et)/(1+exp(et))
%        1 
% 0.7951228 

%%

%uncomment this if you want to save the whole environment
%save ExactFreApproxCylIce


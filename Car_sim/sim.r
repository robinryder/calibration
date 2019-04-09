rm(list=ls())
library(rpart); 
require(MASS); require(mvtnorm); require(mixtools); require(lme4)
source('codes2.r')

#data preparation 
data(car90)
car=data.frame(list(Price=car90$Price/1000,Disp=car90$Disp,Hp.revs=car90$HP.revs,Type=car90$Type,Country=car90$Country)); 
car=car[rowSums(is.na(car))==0,]; n=dim(car)[1]
Z=matrix(1,n,1); X=Z; 
W=cbind(rep(1,n),car$Disp,car$Hp.revs)
y=as.matrix(car$Price); yName=car$Type; Name=levels(car$Type)
r = dim(Z)[2]; p = dim(X)[2]; q = dim(W)[2]; nL=length(Name)

#prior setup
hyper_condC=list(Csig=10,b0=rep(0,dim(X)[2]),gam0=rep(0,dim(W)[2]),alpha=1) #hyperparmeter setup for a posterior approximation 
hyper_MC=list(b0=rep(0,dim(X)[2]),gam0=rep(0,dim(W)[2]),alpha=1) #hyperparmeter setup for an exact posterior 

##########################################################
# T: length of Markov chain (this should be greater than 5e3)
# M: number of particles 

#Exact coverage estimation. Algorithm0
exact=MCestimation(y,X,Z,W,T=1e4,M=3e3,hyper_condC,hyper_MC)
exact_coverage=mean(exact$cover)

#Coverage estimation using an improtance sampling
est=ISestimation(y,X,Z,W,T=1.1e4,M=3.5e3,hyper_condC,hyper_MC)
#est$lam0 : 1/approx likelihood
#est$cover : coverage indicator
#est$ksd : Ks-distance


##########################################################
#Figure 6 implementation
load('result_1e5.RData')

ks0=seq(from=0.11,to=1,length.out=30); #vector of KS distances
cp=errks=ess=nn=rep(0,length(ks0))
for(i in 1:length(ks0)){
  cp[i]=sum(LAM[KSD<ks0[i]]/sum(LAM[KSD<ks0[i]])*COVER[KSD<ks0[i]])
  LAM1=LAM[KSD<ks0[i]]/sum(LAM[KSD<ks0[i]])
  errks[i]=sqrt(sum(LAM1^2*(COVER[KSD<ks0[i]]-cp[i])^2))
  ess[i]=1/sum(LAM1^2); nn[i]=sum(KSD<ks0[i])
}

exact0=0.9794999; sd0=sqrt(exact0*(1-exact0))

dev.new(); 
plot(ks0,cp,xlab='Threshold',ylim=c(0,1),ylab='Coverage probability',cex.axis=1.5,cex.lab=1.5); 
arrows(ks0, cp-2*errks, ks0, cp+2*errks, length=0.05, angle=90, code=3); 
lines(ks0,rep(exact0,length(ks0))); #lines(ks0,rep(exact0-sd0,length(ks0)),col='grey'); 

dev.new(); plot(ks0,nn,xlab='Threshold',log='y',ylab='Number of particles',cex.axis=1.5,cex.lab=1.5);
lines(ks0,ess,type='p',pch=3)

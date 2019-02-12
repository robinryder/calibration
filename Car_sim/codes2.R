#effect matrix
effmx<-function(rf){
  l=levels(rf); Z=matrix(0,length(rf),length(l))
  for(i in 1:length(l)){Z[,i]=as.numeric(rf==l[i])}
  colnames(Z)=l; return(Z); }

################################################
# data simulation 
datasim <- function(y,yName,Name,S,X,Z,W,sig_beta,sig_eta,sig_gam,sig,hyper){
	n=length(y); mu=s=rep(0,n);
	mu=

  simdada <- W%*%hyper$gam0+matrix(rnorm(dim(W)[1],0,sqrt(sig_gam*sig)),ncol=1) + rnorm(dim(W)[1],0,sqrt(sig))
  nL=length(Name)
  for (i in 1:nL){
    simdada=simdada+(yName==Name[i])*(Z%*%rnorm(dim(Z)[2],0,sqrt(sig_eta*sig)))
  }
  for(i in 1:max(S)){
    simdada=simdada+(yName %in% Name[S==i])*(X%*%rnorm(dim(X)[2],0,sqrt(sig_beta*sig)))
  } 
  return(matrix(simdada,ncol=1))
}

################################################
# likelihood
lnlike <- function(y,yName,Name,S,X,Z,W,beta,eta,gam,sig_beta,sig_eta,sig_gam,sig){
  n=length(y); mu=s=rep(0,n); gam=as.matrix(gam); #beta=as.matrix(beta); eta=as.matrix(eta); 
  mu=W%*%gam; s=sig*(sig_beta+sig_eta+sig_gam);
  for (i in 1:nL){mu=mu+(yName==Name[i])*eta[i];}
  for(i in 1:max(S)){mu=mu+(yName %in% Name[S==i])*beta[i];}
  return(sum(dnorm(y,mu,sqrt(s),log=TRUE)))
}


################################################
#posterior sample for sigma

sigpost <- function(CondML_flevel_Model3, flevelM3, y, yName, Name,X, Z, W, A, sig_beta, sig_eta, sig_gam){
  nL = length(Name)
  n = dim(y)[1]
  r = dim(Z)[2]
  p = dim(X)[2]
  q = dim(W)[2]
  b0 = matrix(rep(0, p), ncol = 1); gam0 = matrix(rep(0, q), ncol = 1)
  N=rep(0,nL); for (i in 1:nL) {N[i] = sum(A == i) }
  flevel = flevelM3(y, yName, Name ,X, Z, W, p, q, r, sig_eta)
  ff = CondML_flevel_Model3(flevel, y, yName, Name, X, Z, W, p, q, r, A, N, sig_beta, sig_eta, sig_gam, b0, gam0)
  return(1/rgamma(1,shape=n/2-1/2,rate=ff$ff))
}

################################################
# log-Marlike for effects by levels
flevelM3 <- function(y, yName, Name,X, Z, W, p, q, r, sig_eta) {
	nL=length(Name)
  Sig_eta = diag(r) * sig_eta; iSig_eta = solve(Sig_eta)
	U = matrix(0, q, nL)
	invQ = array(matrix(0, q, q), dim = c(q, q, nL))
	iGk = array(matrix(0, p, p), dim = c(p,p, nL))
	ldetVk = P = rep(0, nL)
	lk1 = array(matrix(0, p, q), dim = c(p, q, nL))
	lk2 = matrix(0, p, nL)

	for (i in 1:nL) {
		subset = which(yName == Name[i])
		Xsub = matrix(X[subset, ],nrow=length(subset),ncol=p); Zsub = matrix(Z[subset, ],nrow=length(subset),ncol=r); 
		ysub = matrix(y[subset, ],ncol = 1); Wsub = matrix(W[subset, ],nrow=length(subset),ncol=q)
	
		Vk = solve(t(Zsub) %*% Zsub + iSig_eta)
		ldetVk[i] = 0.5 * log(det(Vk))
		iGk[, , i] = t(Xsub) %*% Xsub - t(Xsub) %*% Zsub %*% Vk %*% t(Zsub) %*% Xsub
		lk1[, , i] = t(Xsub) %*% Zsub %*% Vk %*% t(Zsub) %*% Wsub - t(Xsub) %*% Wsub
		lk2[, i] = t(Xsub) %*% ysub - t(Xsub) %*% Zsub %*% Vk %*% t(Zsub) %*% ysub
		U[, i] = (-t(Wsub) %*% Zsub %*% Vk %*% t(Zsub) %*% ysub + t(Wsub) %*% ysub)
		invQ[, , i] = -t(Wsub) %*% Zsub %*% Vk %*% t(Zsub) %*% Wsub + t(Wsub) %*% Wsub
		P[i] = 0.5 * t(ysub) %*% Zsub %*% Vk %*% t(Zsub) %*% ysub - 0.5 * t(ysub) %*% ysub
	}
	f <- list(lk1 = lk1, lk2 = lk2, U = U, invQ = invQ, P = P, ldetVk = ldetVk, iGk = iGk)
	return(f)
}

################################################
# log-MarLike for a partition and hyperparameters 

CondML_flevel_Model3 <- function(flevel, y, yName, Name,X, Z, W, p, q, r, A, N, sig_beta, sig_eta, sig_gam, b0, gam0) {

	n = length(y); nL=length(Name)
	Sig_beta = diag(p) * sig_beta
	iSig_beta = solve(Sig_beta)
	Sig_eta = diag(r) * sig_eta
	iSig_eta = solve(Sig_eta)
	Sig_gam = diag(q) * sig_gam
	iSig_gam = solve(Sig_gam)
	
	logf = 0
	U = matrix(0, q, 1)
	invQ = matrix(0, q, q)
	P = 0
	for (j in c(which(N > 0))) {
		ind = which(A == j)
		iGk = matrix(0, p, p)
		ldetVk = 0
		lk1 = matrix(0, p, q)
		lk2 = matrix(0, p, 1)
		for (i in 1:length(ind)) {
			iGk = iGk + flevel$iGk[, , ind[i]]
			lk1 = lk1 + flevel$lk1[, , ind[i]]
			lk2 = lk2 + flevel$lk2[, ind[i]]
			U = U + flevel$U[, ind[i]]
			invQ = invQ + flevel$invQ[, , ind[i]]
			P = P + flevel$P[ind[i]]
			ldetVk = ldetVk + flevel$ldetVk[ind[i]]
		}
		lk2 = lk2 + iSig_beta %*% b0
		G = solve(iGk + iSig_beta)
		U = U + t(lk1) %*% t(G) %*% lk2
		invQ = invQ - t(lk1) %*% t(G) %*% lk1
		P = P - 0.5 * t(b0) %*% iSig_beta %*% b0 + 0.5 * t(lk2) %*% t(G) %*% lk2
		logf = logf + ldetVk + 0.5 * log(det(G)) - p/2 * log(sig_beta)
	}
	U = U + iSig_gam %*% gam0
	invQ = invQ + iSig_gam
	Q = solve(invQ)
	P = P - 0.5 * t(gam0) %*% iSig_gam %*% gam0
	ff = -0.5 * t(U) %*% Q %*% U - P
	logf = logf + 0.5 * log(det(Q)) - 0.5 * q * log(sig_gam) - 0.5 * nL * r * log(sig_eta) + lgamma(n/2) - (n/2) * log(ff)
	return(list(logf=logf, ff=ff))
}

#################################################
# MCMC - run

Model3 <- function(CondML_flevel_Model3, flevelM3, y, yName, Name,X, Z, W, T, C, Alph, hyper) {
	
	dbeta = deta = dgam = 0.5
	update = seq(from = 1, to = T, by = 49)
	nL = length(Name)
	n = dim(y)[1]
	r = dim(Z)[2]
	p = dim(X)[2]
	q = dim(W)[2]

	# hyper prior, beta~N(b0,sig_beta); eta~N(0,sig_eta); epsilon~N(0,sig); gam~N(gam0,sig_gam)
	b0 = hyper$b0; # matrix(rep(0, p), ncol = 1)
	gam0 = hyper$gam0 #matrix(rep(0, q), ncol = 1)
	sig_beta = sig_eta = sig_gam = hyper$Csig
	alpha = hyper$alpha

	A = rep(1, nL)
	#A = c(1:nL)
	N = N1 = rep(0, nL)
	for (i in 1:nL) {N[i] = sum(A == i) }
	k = sum(N > 0)
	flevel = flevelM3(y, yName, Name ,X, Z, W, p, q, r, sig_eta)
	logf0 = CondML_flevel_Model3(flevel, y, yName, Name, X, Z, W, p, q, r, A, N, sig_beta, sig_eta, sig_gam, b0, gam0)$logf
	outS = matrix(, T, nL)
	outk = lnML = Cb = Ceta = Cgam = Calpha = rep(0, T)

	for (t in 1:T) {

		#### partition update ####
		for (s in 1:nL) {
			A1 = A
			A1[s] = sample(c(1:min(nL, (1 + k))), 1)
			for (i in 1:nL) {N1[i] = sum(A1 == i) }
			logf = CondML_flevel_Model3(flevel, y, yName, Name,X, Z, W, p, q, r, A1, N1, sig_beta, sig_eta, sig_gam, b0, gam0)$logf
			h1 = logf + sum(N1 > 0) * log(alpha) + sum(lgamma(N1[N1 > 0])) - logf0 - sum(N > 0) * log(alpha) - sum(lgamma(N[N > 0]))
			if (runif(1) < exp(h1)) {
				logf0 = logf
				ind = which(N1 > 0)
				k = length(ind)
				for (i in 1:k) {
					A[A1 == ind[i]] = i
					N[i] = sum(A1 == ind[i])
				}
				if (k<nL){N[(k + 1):nL] = 0}
			}
		}

		######### hyperparameter update ##########
		
		if (C == "T") { #adaptive MH
			if (t > 49 & t %in% update) {
				dbeta = max(0.001, 1.3 * sd((Cb[(t - 49):(t - 1)])))
				deta = max(0.001, 1.3 * sd((Ceta[(t - 49):(t - 1)])))
				dgam = max(0.001, 1.3 * sd((Cgam[(t - 49):(t - 1)])))
			}

			#### beta update ####       
			sig_beta1 = exp(rnorm(1, log(sig_beta), dbeta))
			#sig_beta1 = rnorm(1,sig_beta,dbeta); if (sig_beta1<0){sig_beta1=sig_beta}
			logf = CondML_flevel_Model3(flevel, y, yName, Name,X, Z, W, p, q, r, A, N, sig_beta1, sig_eta, sig_gam, b0, gam0)$logf
			h1 = logf - logf0 + sum(dgamma(sig_beta1, 2, 0.5, log = T)) - sum(dgamma(sig_beta, 2, 0.5, log = T))
			#h1 = logf - logf0 + (dchisq(1/sig_beta1, 1, log = T)) - (dchisq(1/sig_beta, 1, log = T))
			
			if (runif(1) < exp(h1)) {
				logf0 = logf; sig_beta = sig_beta1
			}

			#### gam update ####
			sig_gam1 = exp(rnorm(1, log(sig_gam), dgam))
			#sig_gam1 = rnorm(1,sig_gam,dgam); if (sig_gam1<0){sig_gam1=sig_gam}
			logf = CondML_flevel_Model3(flevel, y, yName, Name,X, Z, W, p, q, r, A, N, sig_beta, sig_eta, sig_gam1, b0, gam0)$logf
			
			h1 = logf - logf0 + sum(dgamma(sig_gam1, 2, 0.5, log = T)) - sum(dgamma(sig_gam, 2, 0.5, log = T))
			#h1 = logf - logf0 + (dchisq(1/sig_gam1, 1, log = T)) - (dchisq(1/sig_gam, 1, log = T))
			if (runif(1) < exp(h1)) {
				logf0 = logf; sig_gam = sig_gam1
			}

			#### eta update ###
			sig_eta1 = exp(rnorm(1, log(sig_eta), deta))
			#sig_eta1 = rnorm(1,sig_eta,deta); if (sig_eta1<0){sig_eta1=sig_eta}

			flevel_new = flevelM3(y, yName, Name ,X, Z, W, p, q, r, sig_eta1)
			logf = CondML_flevel_Model3(flevel_new, y, yName, Name,X, Z, W, p, q, r, A, N, sig_beta, sig_eta1, sig_gam, b0, gam0)$logf
			h1 = logf - logf0 + sum(dgamma(sig_eta1, 2, 0.5, log = T)) - sum(dgamma(sig_eta, 2, 0.5, log = T))
			#h1 = logf - logf0 + (dchisq(1/sig_eta1, 1, log = T)) - (dchisq(1/sig_eta, 1, log = T))

			if (runif(1) < exp(h1)) {
				logf0 = logf; flevel = flevel_new; sig_eta = sig_eta1
			}
		}

		if (Alph == "T") { #### alpha update ###
			alpha1 = exp(rnorm(1, log(alpha), 3))
			h1 <- k * (log(alpha1/alpha)) + lgamma(alpha1) - lgamma(alpha1 + nL) - lgamma(alpha) + lgamma(alpha + nL) + dgamma(alpha1, shape = 2, scale = 3, log = T) - dgamma(alpha, shape = 2, scale = 3, log = T)
			if (runif(1) < exp(h1)) {alpha = alpha1; }
		}
		
		############# make outputs ###############
		outS[t, ] = as.vector(A)
		outk[t] = k
		lnML[t] = logf0 - n/2 * log(2 * pi)
		Cb[t] = sig_beta
		Ceta[t] = sig_eta
		Cgam[t] = sig_gam
		Calpha[t] = alpha
	}
	mod <- list(outS = outS, outk = outk, lnML = lnML, Cb = Cb, Ceta = Ceta, Cgam = Cgam, Calpha = Calpha)

	return(mod)
}

#################################################


#################################################
# MCMC - given hyperparameters 

MCMC_condC <- function(CondML_flevel_Model3, flevelM3, y, yName, Name,X, Z, W, T, hyper_condC) {
	
	nL = length(Name)
	n = dim(y)[1]
	r = dim(Z)[2]
	p = dim(X)[2]
	q = dim(W)[2]

	# hyper prior, beta~N(b0,sig_beta); eta~N(0,sig_eta); epsilon~N(0,sig); gam~N(gam0,sig_gam)
	b0 = hyper_condC$b0; 
	gam0 = hyper_condC$gam0 
	sig_beta = sig_eta = sig_gam = hyper_condC$Csig
	alpha = hyper_condC$alpha

	A = c(1:nL)
	N = N1 = rep(0, nL)
	for (i in 1:nL) {N[i] = sum(A == i) }
	k = sum(N > 0)
	flevel = flevelM3(y, yName, Name ,X, Z, W, p, q, r, sig_eta)
	logf0 = CondML_flevel_Model3(flevel, y, yName, Name, X, Z, W, p, q, r, A, N, sig_beta, sig_eta, sig_gam, b0, gam0)$logf
	outS = matrix(, T, nL)
	outk = lnML = Cb = Ceta = Cgam = Calpha = rep(0, T)

	for (t in 1:T) {

		#### partition update ####
		for (s in 1:nL) {
		  AA=A; AA[s] = sample(c(1:min(nL, (1 + k))), 1);
		  ss=sort(as.numeric(table(AA)),decreasing=TRUE,index.return=TRUE)
		  AAname=as.numeric(attributes(table(AA))$dimnames$AA)
		  j=0; N1=A1=rep(0,nL)
		  for (i in AAname[ss$ix]) {j=j+1; A1[AA == i] = j; }
		  N1[1:j]=as.numeric(table(A1)); 
		  
			logf = CondML_flevel_Model3(flevel, y, yName, Name,X, Z, W, p, q, r, A1, N1, sig_beta, sig_eta, sig_gam, b0, gam0)$logf
			h1 = logf + sum(N1 > 0) * log(alpha) + sum(lgamma(N1[N1 > 0])) - logf0 - sum(N > 0) * log(alpha) - sum(lgamma(N[N > 0]))
			if (runif(1) < exp(h1)) {
				logf0 = logf; A=A1; N=N1; 
				k = sum(N>0)
			}
		}
		############# make outputs ###############
		outS[t, ] = as.vector(A)
		outk[t] = k
		lnML[t] = logf0 - n/2 * log(2 * pi)
	}
	mod <- list(outS = outS, outk = outk, lnML = lnML)

	return(mod)
}
#################################################
# Marlikelihood MC method #
L_MC <- function(CondML_flevel_Model3, flevelM3, y, yName, Name,X, Z, W, A, N, nv, hyper_MC) {
  r = dim(Z)[2]
  p = dim(X)[2]
  q = dim(W)[2]
  v=rep(0,nv)
  for(i in 1:nv){sig_eta=1/rchisq(1,1); sig_beta=1/rchisq(1,1); sig_gam=1/rchisq(1,1);
  flevel = flevelM3(y, yName, Name ,X, Z, W, p, q, r, sig_eta)
  v[i] = CondML_flevel_Model3(flevel, y, yName, Name, X, Z, W, p, q, r, A, N, sig_beta, sig_eta, sig_gam, b0 = hyper_MC$b0, hyper_MC$gam0)$logf
  }
  logf0 = mean(exp(v-max(v)))+max(v); return(logf0)
}
#################################################
# Posterrior sample, sig_eta, sig_beta, sib_gam #
Sample_IS <- function(CondML_flevel_Model3, flevelM3, y, yName, Name,X, Z, W, A, N, nv, hyper_MC,TT) {
  r = dim(Z)[2]
  p = dim(X)[2]
  q = dim(W)[2]
  mx=matrix(,TT,3); v=rep(0,TT);
  for(i in 1:TT){sig_eta=1/rchisq(1,1); sig_beta=1/rchisq(1,1); sig_gam=1/rchisq(1,1);
  flevel = flevelM3(y, yName, Name ,X, Z, W, p, q, r, sig_eta)
  v[i] = CondML_flevel_Model3(flevel, y, yName, Name, X, Z, W, p, q, r, A, N, sig_beta, sig_eta, sig_gam, b0 = hyper_MC$b0, hyper_MC$gam0)$logf
  mx[i,]=c(sig_eta,sig_beta,sig_gam) }
  lam = exp(v-max(v))/sum(exp(v-max(v))); 
  #ind=sample(TT,nv,prob=lam,replace=TRUE)
  
  #return(matrix(mx[sample(i,nv,prob=lam,replace=TRUE),],nrow=nv))
  return(matrix(mx[sample(i,nv,prob=lam,replace=TRUE),],nrow=nv))
}
#################################################
#################################################
# MCMC - MCmethod

MCMC_MC <- function(CondML_flevel_Model3, L_MC,flevelM3, y, yName, Name,X, Z, W, T, hyper_MC) {
  
  nL = length(Name)
  n = dim(y)[1]
  r = dim(Z)[2]
  p = dim(X)[2]
  q = dim(W)[2]
  
  # hyper prior, beta~N(b0,sig_beta); eta~N(0,sig_eta); epsilon~N(0,sig); gam~N(gam0,sig_gam)
  b0 = hyper_MC$b0; 
  gam0 = hyper_MC$gam0 
  alpha = hyper_MC$alpha
  
  A = c(1:nL)
  N = N1 = rep(0, nL)
  for (i in 1:nL) {N[i] = sum(A == i) }
  k = sum(N > 0); nv=30
  
  logf0=L_MC(CondML_flevel_Model3, flevelM3, y, yName, Name,X, Z, W, A, N, nv, hyper_MC)
  outS = matrix(, T, nL)
  outk = lnML = Calpha = rep(0, T)
  S_table=matrix(A,nrow=1); logf0_table=logf0; 
  
  for (t in 1:T) {
    
    #### partition update ####
    for (s in 1:nL) {
      
      AA=A; AA[s] = sample(c(1:min(nL, (1 + k))), 1);
      ss=sort(as.numeric(table(AA)),decreasing=TRUE,index.return=TRUE)
      AAname=as.numeric(attributes(table(AA))$dimnames$AA)
      j=0; N1=A1=rep(0,nL)
      for (i in AAname[ss$ix]) {j=j+1; A1[AA == i] = j; }
      N1[1:j]=as.numeric(table(A1));   
      
      check=which(apply(S_table,1,function(x) sum(x==A1))==nL)
      if(length(check)>0){logf=logf0_table[check]
      }else{ 
        logf=L_MC(CondML_flevel_Model3, flevelM3, y, yName, Name,X, Z, W, A1, N1, nv, hyper_MC)
            S_table=rbind(S_table,A1);  logf0_table=c(logf0_table,logf) }
      h1 = logf + sum(N1 > 0) * log(alpha) + sum(lgamma(N1[N1 > 0])) - logf0 - sum(N > 0) * log(alpha) - sum(lgamma(N[N > 0]))
      if (runif(1) < exp(h1)){ logf0 = logf; A=A1; N=N1; k=sum(N>0) }
    }
    ############# make outputs ###############
    outS[t, ] = as.vector(A)
    outk[t] = k
    lnML[t] = logf0 - n/2 * log(2 * pi)
  }
  mod <- list(outS = outS, outk = outk, lnML = lnML)
  
  return(mod)
}

################################################
# MLE for effects by levels  
flevelMLE <- function(y, yName, Name,X, Z, W, p, q, r, A, N, b, gam, eta) {
  nL=length(Name)
  yy=rep(0,length(y))
  for (i in 1:nL) {
    subset = which(yName == Name[i])
    Xsub = matrix(X[subset, ],nrow=length(subset),ncol=p); Zsub = matrix(Z[subset, ],nrow=length(subset),ncol=r); 
    ysub = matrix(y[subset, ],ncol = 1); Wsub = matrix(W[subset, ],nrow=length(subset),ncol=q)
    eta=ginv(t(Zsub)%*%Zsub)%*%t(Zsub)%*%matrix(ysub-Xsub%*%b-Wsub%*%gam,ncol=1)
    yy[subset]=Zsub%*%eta
  }
  for (i in c(which(N > 0))) {
    subset = which(yName %in% Name[A==i]);
    Xsub = matrix(X[subset,],nrow=length(subset),ncol=p); Zsub = matrix(Z[subset, ],nrow=length(subset),ncol=r); 
    ysub = matrix(y[subset, ],ncol = 1); Wsub = matrix(W[subset, ],nrow=length(subset),ncol=q)
    b=ginv(t(Xsub)%*%Xsub)%*%t(Xsub)%*%matrix(ysub-Zsub%*%eta-Wsub%*%gam,ncol=1)
    yy[subset]=yy[subset]+Xsub%*%b
  }
  gam=ginv(t(W)%*%W)%*%t(W)%*%matrix(yy-Z%*%eta-X%*%b,ncol=1)
  logf=sum(dnorm(y-yy+W%*%matrix(gam,ncol=1),0,mean((y-yy+W%*%matrix(gam,ncol=1))^2),log=T))
  return(logf=logf)
}
#################################################

#################################################
# overall log-marginal likelihood approximation
ML_is <- function(flevelM3, CondML_flevel_Model3,y, yName,Name,X,Z,W,outS,lnML,Cb,Ceta,Cgam,Calpha,Nis,M){
	
	n = length(y);
	r = dim(Z)[2]; p = dim(X)[2]; q = dim(W)[2]
	# hyper prior, beta~N(b0,sig_beta); eta~N(0,sig_eta); epsilon~N(0,sig); gam~N(gam0,sig_gam)\n
	b0 = matrix(rep(0, p), ncol = 1)
	gam0 = matrix(rep(0, q), ncol = 1)
	
	# Nis : number of samples
	# M : number of posterior samples used to construct an importance function
	indM=sample(c(1:dim(outS)[1]),M,replace=TRUE)
	s_alpha=sd(log(Calpha))*1.06/M^5; s_b=sd(log(Cb))*1.06/M^5
	s_eta=sd(log(Ceta))*1.06/M^5; s_gam=sd(log(Cgam))*1.06/M^5
	m_alpha=log(Calpha[sample(indM, Nis,replace=TRUE)])
	m_b=log(Cb[sample(indM,Nis,replace=TRUE)])
	m_eta=log(Ceta[sample(indM,Nis,replace=TRUE)])
	m_gam=log(Cgam[sample(indM,Nis,replace=TRUE)])
	
	alpha1=exp(m_alpha+rnorm(Nis,0,s_alpha))
	b1=exp(m_b+rnorm(Nis,0,s_b))
	eta1=exp(m_eta+rnorm(Nis,0,s_eta))
	gam1=exp(m_gam+rnorm(Nis,0,s_gam))
	logf=qq=pr=rep(0,Nis); nL=length(Name)

	for (i in 1: Nis){
		A=outS[sample(indM,1), ]; 
		N=rep(0,nL); N[1:length(table(A))]=table(A)
		flevel = flevelM3(y, yName, Name ,X, Z, W, p, q, r, eta1[i])
		logf[i] = CondML_flevel_Model3(flevel,y, yName, Name,X, Z, W, p, q, r, A, N, b1[i], eta1[i],gam1[i], b0, gam0) - n/2 * log(2 * pi)
		qq[i] = log(mean(dnorm(log(alpha1[i]),m_alpha,s_alpha)))+log(mean(dnorm(log(eta1[i]),m_eta,s_eta)))+log(mean(dnorm(log(gam1[i]),m_gam,s_gam)))+log(mean(dnorm(log(b1[i]),m_b,s_b))) - log(M)
		pr[i] = sum(dgamma(c(eta1[i],gam1[i],b1[i]),2,0.5,log=TRUE))+sum(N>0)*log(alpha1[i]) + lgamma(alpha1[i]) - lgamma(alpha1[i] + nL) + dgamma(alpha1[i], shape = 2, scale = 3, log = T)	
	}
	w=logf+pr-qq; 
	ML=log(mean(exp(w-max(w))))+max(w); return(ML)
}


#################################################
#standardising partitions
gm2mg<-function(GM,M) {
  #convert groups listing members to list of members indicating groups
  group.indices=lapply(GM,function(x){which(is.element(as.matrix(M),x))})
  MG=rep(NA,length(M))
  group.sizes=unlist(lapply(group.indices,length))
  MG[unlist(group.indices)]=rep(1:length(GM),group.sizes)
  return(MG)
}

mg2gm<-function(MG,M) {
  #convert list of members indicating groups to groups listing members
  GM=lapply(as.list(unique(MG)), function(x) {M[MG==x]})
  return(GM)
}

prettyPart<-function(xx,AA) {
  #partition to string, pretty format
  ts=mg2gm(gm2mg(xx,AA),AA); 
  paste(lapply(ts,function(x){sprintf("(%s)",paste(x,collapse=","))}),collapse=",")
}
#################################################
# posterior for new partitions
hpdpartition_new <- function(mg2gm,gm2mg,S,p,Name){ T=dim(S)[1]; ss=rep(0,T)
for (t in 1:T){ x=gm2mg(mg2gm(S[t,],Name),Name)
ss[t]=paste(x,collapse="-") }
g=table(ss)
partnames=attributes(g)$dimnames[[1]]
sg=sort(as.numeric(g),decreasing=TRUE,index.return=TRUE)
pp=sg$x/sum(sg$x); ppnames=partnames[sg$ix]
if(pp[1]>=p){hpdnames=ppnames[1]; hpd=pp[1]
}else{
  hpd=pp[which(cumsum(pp)<p)]; hpdnames=ppnames[which(cumsum(pp)<p)] }
return(list(hpd=hpd,hpdnames=hpdnames))
}
#################################################
# coverage
Acoverage <- function(A,Alist){
  A=as.numeric(unlist(strsplit(A,"-",fixed=TRUE)))
  nL=length(A); mx=rep(0,nL);
  for(i in 1:length(Alist)){
        mx=rbind(mx,as.numeric(unlist(strsplit(Alist[i],"-",fixed=TRUE))))
      }
  mx=matrix(mx[-1,],ncol =nL); k=max(A); ps=perm(k,k,1:k); cover=0; A1=A;
  for(i in 1:dim(ps)[1]){ for(j in 1:k){A1[A==ps[i,j]]=j}
    match=apply(mx,1,function(x) sum(x==A1))
    cover=cover+sum(match==nL)
    #if(sum(match==nL)){print(A1)}
  }
  return(sum(cover>0))
}

#################################################
KS.dist <- function(s0,p0,s1,p1){ 
  #s0,p0 - partition and prob for the obs in the decreasing order
  #s1,p1 - partition and prob for the sim data in the decreasing order
  
  p0=cumsum(p0); ord=match(s0,s1); 
  if(sum(is.na(ord)==TRUE)>0){
    p2=rep(0,length(p0)); p2[is.na(ord)==FALSE]=p1[ord[is.na(ord)==FALSE]]; p2=cumsum(p2)
    rho=p0-p2;
    }else{rho=p0[is.na(ord)==FALSE]-cumsum(s1[ord[is.na(ord)==FALSE]]) }
  return(abs(rho))
}


##########################################################

post.samp2<- function(lnlike,y, yName, Name, X, Z, W, A, sig_beta, sig_eta, sig_gam, sig, T){
  nL = length(Name)
  n = dim(y)[1]
  r = dim(Z)[2]
  p = dim(X)[2]
  q = dim(W)[2]
  
  ZZ=matrix(0,dim(Z)[1],dim(Z)[2]*nL); for(i in 1:nL){ZZ[yName==Name[i],((i-1)*r+1):(i*r)]=Z[yName==Name[i],] }
  w0=rep(0,T); #importance weight, prior and post approx 
  
  Mx.eta=matrix(,T,nL); Mx.gam=matrix(,T,q); Mx.b=array(,dim=c(max(A),T,p)); ysim=matrix(0,dim(y)[1],T)
  Mx.gam[1,]=solve(t(W)%*%W)%*%t(W)%*%y; Mx.eta[1,]=rnorm(nL,0,sqrt(sig_eta*sig));
  y1=y-W%*%matrix(Mx.gam[1,],ncol=1)-ZZ%*%matrix(Mx.eta[1,],ncol=1)
  for(i in 1:max(A)){ ind=yName%in%Name[A==i]
  ss=ginv(t(X[ind,])%*%X[ind,]/sig+diag(p)/sig_beta*sig)
  Mx.b[i,1,] = mvrnorm(1,ss%*%t(X[ind,])%*%y1[ind]/sig,ss ) 
  ysim[ind,1]=X[ind,]%*%matrix(Mx.b[i,1,],ncol=1)}
  ysim[,1]=ysim[,1]+W%*%matrix(Mx.gam[1,],ncol=1)+ZZ%*%matrix(Mx.eta[1,],ncol=1)
  w0[1] = lnlike(matrix(ysim[,1],ncol=1),yName,Name,A,X,Z,W,Mx.b[,1,],Mx.eta[1,],Mx.gam[1,],sig_beta,sig_eta,sig_gam,sig)
  
  for(t in 2:T){
    #gam generate
    yy=rep(0,dim(y)[1])
    for(i in 1:max(A)){ ind=yName%in%Name[A==i]; 
    yy[ind]=matrix(y[ind]-X[ind,]%*%matrix(Mx.b[i,(t-1),],ncol=1)-ZZ[ind,]%*%matrix(Mx.eta[t-1,],ncol=1),ncol=1)
    }
    ss=ginv(t(W)%*%W/sig+diag(q)/sig_gam*sig)
    Mx.gam[t,] = mvrnorm(1,ss%*%t(W)%*%matrix(yy,ncol=1)/sig,ss)
    
    #eta generate
    yy=rep(0,dim(y)[1])
    for(i in 1:max(A)){ ind=yName%in%Name[A==i]; 
    yy[ind]=matrix(y[ind]-X[ind,]%*%matrix(Mx.b[i,(t-1),],ncol=1)-W[ind,]%*%matrix(Mx.gam[t,],ncol=1),ncol=1)
    }
    ss=ginv(t(ZZ)%*%ZZ/sig+diag(dim(ZZ)[2])/sig_eta*sig)
    Mx.eta[t,] = mvrnorm(1,ss%*%t(ZZ)%*%matrix(yy,ncol=1)/sig,ss)
    
    #b generate
    yy=y-W%*%matrix(Mx.gam[t,],ncol=1)-ZZ%*%matrix(Mx.eta[t,],ncol=1)
    for(i in 1:max(A)){ ind=yName%in%Name[A==i];
    ss=ginv(t(X[ind,])%*%X[ind,]/sig+diag(p)/sig_beta*sig)
    Mx.b[i,t,] = mvrnorm(1,ss%*%t(X[ind,])%*%yy[ind]/sig,ss)
    ysim[ind,t]=X[ind,]%*%matrix(Mx.b[i,t,],ncol=1) }
    
    ysim[,t]=ysim[,t]+W%*%matrix(Mx.gam[t,],ncol=1)+ZZ%*%matrix(Mx.eta[t,],ncol=1)
    
    w0[t] = lnlike(matrix(ysim[,t],ncol=1),yName,Name,A,X,Z,W,Mx.b[,t,],Mx.eta[t,],Mx.gam[t,],sig_beta,sig_eta,sig_gam,sig)
  }
  w1=1/exp(w0-min(w0)); ind=sample(T,T/2,prob=w1/sum(w1),replace=TRUE)
  ysim0=rowMeans(ysim[,ind])
  
  return(list(Mx.b=Mx.b,Mx.eta=Mx.eta,Mx.gam=Mx.gam,ysim=ysim,ysim0=ysim0))
}
##########################################################

ISestimation <- function(y,X,Z,W,T=1e4,M=10,hyper_condC,hyper_MC){
  #MCMC for a posterior approximation
  result=MCMC_condC(CondML_flevel_Model3, flevelM3, y, yName, Name, X, Z, W, T, hyper_condC)
  #credible set estimation
  cred0=hpdpartition_new(mg2gm,gm2mg,result$outS[5e3:T,],p=0.95,Name)
  #Importance sampling 
  ind=sample(c(5001:T),M,replace=FALSE)
  Amx=result$outS[ind,]; Nmx=result$outk[ind]; 
  Amx_ch=apply(Amx,1,function(x) paste(x,collapse='-'))
  Amx_name=attributes(table(Amx_ch))$dimnames[[1]]; Amx_n=as.numeric(table(Amx_ch)); 
  l=0; lam=cover=ksd=rep(0,M); ysim_mx=matrix(,M,dim(y)[1]); Smx=array(,dim=c(M,5e3,length(Name))); 
  para=list(); para$eta=matrix(,T,6); para$gam=matrix(,T,3); para$beta=matrix(0,T,6); para$sig=matrix(,T,4);
  for(j in 1:length(Amx_n)){ 
    Aj=as.numeric(unlist(strsplit(Amx_name[j],"-",fixed=TRUE))); nvj=Amx_n[j];
    Nj=rep(0,length(Name)); Nj[1:max(Aj)]=as.numeric(table(Aj))
    Csig=Sample_IS(CondML_flevel_Model3, flevelM3, y, yName, Name,X, Z, W,Aj, Nj, nvj, hyper_MC, T=3e3)
    for(i in 1:nvj){ l=l+1; 
    sig=sigpost(CondML_flevel_Model3, flevelM3, y, yName, Name,X, Z, W, Aj,sig_beta=Csig[i,2],sig_eta=Csig[i,1],sig_gam=Csig[i,3])
    ps=post.samp2(lnlike,y, yName, Name, X, Z, W, Aj, sig_beta=Csig[i,2],sig_eta=Csig[i,1],sig_gam=Csig[i,3], sig, T=2e3)
    ysim=matrix(ps$ysim0,ncol=1);
    ysim_mx[l,]=ysim; #simulated data
    
    result0=MCMC_condC(CondML_flevel_Model3, flevelM3, ysim, yName, Name,X, Z, W, T, hyper_condC) 
    cred=hpdpartition_new(mg2gm,gm2mg,result0$outS[5e3:T,],p=0.95,Name)
    cover[l]=Acoverage(paste(Aj,collapse="-"),cred$hpdnames)
    Smx[l,,]=result0$outS[5001:T,]; ksd[l]=max(KS.dist(cred0$hpdnames,cred0$hpd,cred$hpdnames,cred$hpd))
    flevel = flevelM3(y, yName, Name ,X, Z, W, p=dim(X)[2], q=dim(W)[2], r=dim(Z)[2], sig_eta=hyper_condC$Csig)
    lam[l] = CondML_flevel_Model3(flevel, y, yName, Name, X, Z, W, p=dim(X)[2], q=dim(W)[2], r=dim(Z)[2], Aj, Nj, sig_beta=hyper_condC$Csig, sig_eta=hyper_condC$Csig, sig_gam=hyper_condC$Csig, b0=hyper_condC$b0, gam0=hyper_condC$gam0)$logf
    }}
  lam0=1/exp(lam)
  
  return(list(ysim_mx=ysim_mx, Smx=Smx, cover=cover, lam0=lam0))
}

##########################################################

MCestimation <- function(y,X,Z,W,T=1e4,M=10,hyper_condC,hyper_MC){
  
  result=MCMC_MC(CondML_flevel_Model3,L_MC,flevelM3, y, yName, Name,X, Z, W, T, hyper_MC)
  
  ind=sample(c(5e3:dim(result$outS)[1]),M,replace=TRUE)
  Amx=result$outS[ind,]; Nmx=result$outk[ind]; lam=cover=dist=ksdist=rep(0,M); 
  Amx_ch=apply(Amx,1,function(x) paste(x,collapse='-'))
  Amx_name=attributes(table(Amx_ch))$dimnames[[1]]; Amx_n=as.numeric(table(Amx_ch)); 
  
  l=0; 
  for(j in 1:length(Amx_n)){ 
    Aj=as.numeric(unlist(strsplit(Amx_name[j],"-",fixed=TRUE))); nvj=Amx_n[j];
    for(i in 1:nvj){ l=l+1; #print(l)
    result0=MCMC_condC(CondML_flevel_Model3, flevelM3, y, yName, Name, X, Z, W, T, hyper_condC)
    cred=hpdpartition_new(mg2gm,gm2mg,result0$outS[sample(c(5e3:T),100),],p=0.95,Name)
    cover[l]=Acoverage(paste(Aj,collapse="-"),cred$hpdnames)
    } }
  return(list(cover=cover))
}
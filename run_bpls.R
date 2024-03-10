rm(list = ls())

source("./extrafunctions.R")
source("./loadPackages.R")

set.seed(123)
P = 1000 # dim(X) 
Q = 10  # dim(Z) 
R = 4   # dim(Y) 
N = 500 # training samples 
Ntest = 1000# test samples 
Nfull = N+Ntest


Sig2true = 0.1
Psi2true = 0.1



Wtrue = matrix(sqrt(1/Q),P,Q) * rnorm(P*Q)
Ctrue = matrix(sqrt(1/Q),R,Q) * rnorm(R*Q)

SpikeAndSlab = 0
LASSO=0




p_prior = 0.5
if(SpikeAndSlab){
	btrue = rep(c(0,1),Q/2)

	Btrue = diag(btrue)
	Ctrue = Ctrue %*% Btrue *sqrt(2)


}

if(LASSO){
	qq = Q/2
	for(r in 1:R){
		Ctrue[r,sample(1:Q, qq, F)] = 0
	}
	Ctrue = Ctrue*sqrt(2)
}


isoX = 1
isoY = 1
sparsity = 0
shrinkage = 1

N_MC = 2e4; N_MC_out = 0.7*N_MC; N_pred = N_MC_out
nthin = N_MC_out/N_pred
storeIDX = (N_MC-N_MC_out)+1:N_pred * nthin


Ztrue = matrix(rnorm(N*Q),nrow = N, ncol = Q)

X = Ztrue %*%t(Wtrue) + sqrt(Sig2true) *rnorm(N*P)
Y = Ztrue %*% t(Ctrue) + sqrt(Psi2true) * rnorm(N*R)

Ztest = matrix(rnorm(Ntest*Q),nrow = Ntest, ncol = Q)

Xtest = Ztest %*% t(Wtrue) + sqrt(Sig2true) *rnorm(Ntest*P)
Ytest = Ztest %*% t(Ctrue) + sqrt(Psi2true) * rnorm(Ntest*R)
Xfull = rbind(X,Xtest)
Yfull = rbind(Y,Ytest)
Xfull = std(Xfull)
Yfull = std(Yfull)

X = Xfull[1:N,] ; Xtest = Xfull[-(1:N),]
Y = Yfull[1:N,] ; Ytest = Yfull[-(1:N),]


Asig = 2.5; Bsig = 0.1
Apsi = 2.5; Bpsi = 1.5



shrinkX = 1
truncW = 0.0
shrinkY = 1
truncC = 0

nuW = c(2,3);nuC = c(2,3);
alphaW = c(1.01,2.2)
betaW  = c(1,1)


# alphaC = c(1.01,2.1)
# betaC  = c(1,1)
Qs = min(N-1,15)


sparsity = 0
source("./bpls_mcmc.R")
Yhat = matrix(NA, nrow = Ntest, ncol = R)

RMSE = rep(NA,R)
for(r in 1:R){
  Yhat[,r] = apply(EYpred[,,r],2,mean)
  RMSE[r] = sqrt(mean((Yhat[,r] - Ytest[,r])^2 ))
}


sparsity = 1
source("./bpls_mcmc.R")
Yhatss = matrix(NA, nrow = Ntest, ncol = R)

RMSEss = rep(NA,R)
for(r in 1:R){
  Yhatss[s,r] = apply(EYpred[,,r],2,mean)
  RMSEss[r] = sqrt(mean((Yhatss[,r] - Ytest[,r])^2 ))
}



sparsity = 0
source("./bpls_mcmc_lasso.R")
YhatL = matrix(NA, nrow = Ntest, ncol = R)

RMSEL = rep(NA,R)
for(r in 1:R){
  YhatL[,r] = apply(EYpred[,,r],2,mean)
  RMSEss[r] = sqrt(mean((YhatL[,r] - Ytest[,r])^2 ))
}

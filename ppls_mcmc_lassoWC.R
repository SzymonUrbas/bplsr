# LASSO on C
logdinvgamma = function(x, a, b){
  dgamma(1/x, shape = a, rate = b, log = TRUE) - 2 *log(x)
}
########################## INITIALISE
set.seed(NULL)
# par(mfrow = c(1,2))
Alam = 1
Blam = Alam
sparsity =0
# Qs = 15
mcmc_time = rep(0,7)
if(1){
  truncWC = truncW
  


  # lukewarm start
  if(1){

    pcaX = prcomp(X)
    # plot(pcaX$sdev,type = 'b', xlim = c(1,15))

    # cx <- c(0,cumsum(pcaX$sdev))
    # n = 5
    # rsum <- (cx[(n+1):length(cx)] - cx[1:(length(cx) - n)]) / n
    # plot(rsum, xlim = c(1,30))
    if(!exists('Qs')){
      Qs = min(22,N-1)
    }
    
    Wcurr = pcaX$rotation[,1:Qs]
    # mplot(Wcurr, digits = NA, border = F)
    Sig2curr = rep(mean(pcaX$sdev[-(1:Qs)])^2,P)


    # Psi2curr = 1/rgamma(R,Apsi,Bpsi)

    if(isoY){
      Psi2curr=rep(1/rgamma(1,Apsi,Bsig),R)
    } else {
      Psi2curr=1/rgamma(R,Apsi,Bsig)
    }

    SIGMAinv = diag(1/Sig2curr)
    if(R>1){
      PSIinv   = diag(1/Psi2curr)
    } else{
      PSIinv = matrix(1/Psi2curr,1,1)
    }

    Sinv = Id(Qs) + t(Wcurr)%*%SIGMAinv %*%Wcurr
    S = solve(Sinv)
    # plot(t(Wcurr)%*%SIGMAinv %*%Wcurr)

    Zcurr =  X %*%t(S%*% t(Wcurr)%*%SIGMAinv)  
    Ccurr = t(solve(t(Zcurr)%*%Zcurr)%*%t(Zcurr)%*%Y)
    # Ccurr = matrix(0,nrow = R,ncol = Qs)

    phiW= rgamma(P*Qs,nuW[1],nuW[2])

    L_lasso = 1


    phiC = 1/rexp(R*Qs,L_lasso^2/2)
    
    deltaWC = rep(1,Qs)
    if(shrinkage){
      # deltaWC = c(1.0, rltrgamma(Qs-1,alphaW[2],betaW[2],truncWC))
       deltaWC = c(1.0, rep(alphaW[2],Qs-1))
    }

    tauWC = cumprod(deltaWC)

    DELTA = matrix(1/(phiW*tauWC), ncol = Qs, byrow = T)
    OMEGA = matrix(1/(phiC*tauWC), ncol = Qs, byrow = T)


  # Ccurr = sqrt(OMEGA) * rnorm(R*Qs);

    Bcurr = diag(c(rep(as.logical(!sparsity), Qs)))
  }
  

  ####

  sumBstore = rep(NA, N_MC)
  storeB = matrix(NA,nrow = N_MC,ncol = Qs)
  storeCB = array(NA, dim = c(N_MC,R,Qs))
  storeC = array(NA, dim = c(N_MC,R,Qs))
  storeW = array(NA, dim = c(N_MC,P,Qs))
  storeCBa = array(NA, dim = c(N_MC,R,P))
  # storeW = array(NA, dim = c(N_MC,P,Qs))

  storeSig2 = matrix(NA,nrow = N_MC,ncol = P)
  storePsi2 = matrix(NA,nrow = N_MC, ncol = R)
  storePhiW = array(NA, dim = c(N_MC,P,Qs))
  storePhiC = array(NA, dim = c(N_MC,R,Qs))
  storeTauWC = matrix(NA,nrow = N_MC,ncol = Qs)
  storedelta = matrix(NA,nrow = N_MC,ncol = Qs)
  storeL = rep(NA, N_MC)
  
  # storeTauC = matrix(NA,nrow = N_MC, ncol = Qs)

  EYpred = array(NA, dim = c(N_pred,Ntest,R))
  Ypred = array(NA, dim = c(N_pred,Ntest,R))
  storeVarY = array(NA, dim = c(N_pred, R, R))
  iter = 1

  lp_store = rep(NA, N_MC)
}
######################### RUN CHAIN 
# cat("\f")
pb <- progress_bar$new(                                                             
  format = "  [:bar] :percent eta: :eta",                               
  total = N_MC, clear = FALSE, width= 60)  
for(iter in 1:N_MC){
  pb$tick()
  # if(iter){print(paste('--------',iter, '/',N_MC))}

  # print('Z_update')
  ### Z update
  t_tmp = Sys.time()
  Sinv = Id(Qs) + t(Wcurr)%*%SIGMAinv %*%Wcurr +
            t(Bcurr)%*%t(Ccurr)%*%PSIinv %*%Ccurr%*%Bcurr

  S = solve(Sinv)
  Zcurr = (X%*%SIGMAinv%*%Wcurr + Y%*%PSIinv%*%Ccurr%*%Bcurr)%*%S + matrix(rnorm(N*Qs),nrow = N)%*%.chol(S)

  mcmc_time[1] = mcmc_time[1]+(Sys.time()-t_tmp)

  ### W update
  # print('W_update')
  t_tmp = Sys.time()
  ZtZ = t(Zcurr) %*%Zcurr
  # mplot(abs(ZtZ))
  for(j in 1:P){
    choltDj = .chol(diag(1/DELTA[j,]) + SIGMAinv[j,j] * ZtZ)
    # choltDj = chol(diag(1/DELTA[j,]) + SIGMAinv[j,j] * ZtZ)


    Wcurr[j,] = backsolve(choltDj,
      backsolve(choltDj,SIGMAinv[j,j]*t(Zcurr)%*%X[,j],k=Qs,transpose=1)+rnorm(Qs),
      k=Qs)

  }
  mcmc_time[2] = mcmc_time[2]+(Sys.time()-t_tmp)
  
  # print('C_update')
  t_tmp = Sys.time()
  for(j in 1:R){
    choltOj = .chol(diag(1/OMEGA[j,]) +PSIinv[j,j]* Bcurr %*% ZtZ %*% Bcurr)
    Ccurr[j,] = backsolve(choltOj,
      backsolve(choltOj,PSIinv[j,j]* Bcurr%*%t(Zcurr)%*%Y[,j],transpose=1,k=Qs)+rnorm(Qs),
      k=Qs)
  }
  mcmc_time[3] = mcmc_time[3]+(Sys.time()-t_tmp)
  # mplot(Ccurr,digits = NA)

  # Bcurr = diag(as.numeric(runif(Q)<0.5))
  # Bcurr = diag(rep(1,Q))
  # print('B_update')
  t_tmp = Sys.time()

  mcmc_time[4] = mcmc_time[4]+(Sys.time()-t_tmp)
  # update precisions

  # print('PRECISION_update')

  t_tmp = Sys.time()

  phiW = rgamma(P*Qs,nuW[1]+0.5,1)/(nuW[2]+0.5*Wcurr^2 %*%diag(tauWC))


  # phiC = rgamma(R*Qs,nuC[1]+0.5,1)/(nuC[2]+0.5*Ccurr^2 %*%diag(tauWC))
  # meanPhiC = L_lasso*abs(Ccurr%*%diag(tauWC))

  
  # phiC = matrix(rinvgauss(R*Qs,mean = as.vector(meanPhiC), shape  = L_lasso),nrow = R,ncol = Qs) 


  meanPhiC = L_lasso / ( abs(Ccurr)%*%diag(sqrt(tauWC)) )
  phiC = matrix(rinvgauss(R*Qs,mean = as.vector(meanPhiC), shape  = L_lasso^2),nrow = R,ncol = Qs)
  L_lasso = sqrt(rgamma(1,Alam + Qs*R,Blam + 0.5*sum(1/phiC)))
  # L_lasso = sqrt(rltrgamma(1,Alam + Qs*R,Blam + 0.5*sum(1/phiC),0.05))
  # L_lasso = 0.25






  if(shrinkage){

    sum_phiW_W2 = colSums(phiW*Wcurr^2)

    sum_phiC_C2 = colSums(phiC*Ccurr^2 )

    SumBoth = sum_phiW_W2 + sum_phiC_C2
    # plot(log(SumBoth),type = 'b', main = iter)
     # deltaWC[1] = rgamma(1,alphaW[1] + 0.5*P*Qs,betaW[1]+0.5*sum(tauWC/deltaWC[1]*SumBoth))
    for(j in 2:Qs){
      deltaWC[j] = rltrgamma(1,alphaW[2]+0.5*(P+R)*(Qs-j+1),
                                betaW[2]+ 0.5*sum( tauWC[j:Qs]/deltaWC[j] * SumBoth[j:Qs]),truncWC)  
      tauWC = cumprod(deltaWC)
    }
    # plot(log(tauWC),type = 'b', main = round((Alam + Qs*R)/(Blam + 0.5*sum(1/phiC)),2))
    # if(any(tauWC<1)){print(L_lasso)}

    DELTA = 1/ (phiW %*%diag(tauWC))
    OMEGA = 1/ (phiC %*%diag(tauWC))

    # DELTA = matrix(1/(phiW*tauWC), ncol = Qs, byrow = T)


    # OMEGA = matrix(1/(phiC*tauWC), ncol = Qs, byrow = T)
    # plot(OMEGA,digits = 2)
  }
  mcmc_time[5] = mcmc_time[5]+(Sys.time()-t_tmp)
  # update noise variances
  # sum(X[,j] - Zcurr %*% t(Wcurr[,j]))
  t_tmp = Sys.time()
  if(isoX){
    Sig2curr=rep(1/rgamma(1,Asig+N*P/2,Bsig + 0.5 *sum(( X - Zcurr %*% t(Wcurr))^2)),P)
    } else {
    Sig2curr=1/rgamma(P,Asig+N/2,Bsig + 0.5 *colSums(( X - Zcurr %*% t(Wcurr))^2))
  }
  if(isoY){
    Psi2curr=rep(1/rgamma(1,Apsi+N*R/2,Bpsi + 0.5 *sum(( Y - Zcurr %*%Bcurr%*% t(Ccurr))^2)),R)
  } else {
    Psi2curr=1/rgamma(R,Apsi+N/2,Bpsi + 0.5 *colSums(( Y - Zcurr %*%Bcurr%*% t(Ccurr))^2))
    # Psi2curr=1/rgamma(R,Apsi+N/2,Bsig + 0.5 *colSums(( Y - Zcurr %*%Bcurr%*% t(Ccurr)))^2)
  }
  mcmc_time[6] = mcmc_time[6]+(Sys.time()-t_tmp)


  SIGMAinv = diag(1/Sig2curr)

  if(R>1){
    PSIinv   = diag(1/Psi2curr)
  } else{
    PSIinv = matrix(1/Psi2curr,1,1)
  }
  


  ### PREDICT??

  t_tmp = Sys.time()
  CB = Ccurr%*%Bcurr
  Sinv = Id(Qs) + t(Wcurr)%*%SIGMAinv %*%Wcurr
  S = solve(Sinv)
  a = S%*% ( t(Wcurr)%*%SIGMAinv)
  CBa = CB %*% a
  if(iter %in% storeIDX){
    idx_pred = (iter - (N_MC - N_MC_out))/nthin


    # cholSinv = .chol(Sinv)
    # # sqrtS = chol(S)
    if(R>1){
      storeVarY[idx_pred,,] = CB %*% S %*% t(CB) + diag(Psi2curr)
      cholVarY = t(chol(CB %*% S %*% t(CB) + diag(Psi2curr)))
    } else{
      storeVarY[idx_pred,,] = CB %*% S %*% t(CB) + matrix(Psi2curr,1,1)
      cholVarY = t(chol(CB %*% S %*% t(CB) + matrix(Psi2curr,1,1)))
    }


    for(i in 1:Ntest){


      EYpred[idx_pred,i,] = CBa %*% Xtest[i,]

      Ypred[idx_pred,i,] = EYpred[idx_pred,i,] + cholVarY %*%rnorm(R) 
    }




  }


  mcmc_time[7] = mcmc_time[7]+(Sys.time()-t_tmp)

  ### STORAGE
  storeW[iter,,] = Wcurr
  storeC[iter,,] = Ccurr
  storePhiW[iter,,] = phiW
  storePhiC[iter,,] = phiC
  storeTauWC[iter,] = tauWC
  storedelta[iter,] = deltaWC


  storeCB[iter,,] = CB
  storeCBa[iter,,] = CBa

  # storeB[iter,] = diag(Bcurr)

  storeSig2[iter,] = Sig2curr
  storePsi2[iter,] = Psi2curr

  storeL[iter] = L_lasso
  # plot(storeL[1:iter], type = 'b', main = iter)

  # ## calculater log_pi ----
  # lp_curr = 0
  # # likelihood
  # # lp_curr = lp_curr + 0.5*N*log(det(SIGMAinv)*det(PSIinv))
  # lp_curr = lp_curr - 0.5*N* (sum(log(Sig2curr))+sum(log(Psi2curr)))

  # # lp_curr = lp_curr - 0.5 * sum((X - Zcurr %*% t(Wcurr)) %*% SIGMAinv 
  # #                                 %*% t(X - Zcurr %*% t(Wcurr)) + 
  # #                                 (Y - Zcurr %*% t(Ccurr)) %*% PSIinv 
  # #                                 %*% t(Y - Zcurr %*% t(Ccurr)) + 
  # #                                  Zcurr %*% t(Zcurr))

  # lp_curr = lp_curr - 0.5 *sum(sapply(1:N, function(n){
  #           tmp = 0
  #           tmp = tmp + t( X[n,] - Wcurr %*% Zcurr[n,] )%*%SIGMAinv%*%( X[n,] - Wcurr %*% Zcurr[n,] )
  #           tmp = tmp + t( Y[n,] - Ccurr %*% Zcurr[n,] )%*%PSIinv%*%( Y[n,] - Ccurr %*% Zcurr[n,] )
  #           tmp = tmp + sum(Zcurr[n,]^2)
  #           return(tmp)
  #           }))

  # #loading priors
  # lp_curr = lp_curr + sum(dnorm(Wcurr, sd = sqrt(DELTA), log = TRUE))

  # lp_curr = lp_curr + sum(dnorm(Ccurr, sd = sqrt(OMEGA), log = TRUE))

  # #loading hyperpriors
  # lp_curr = lp_curr + sum(dgamma(phiW, shape = nuW[1], rate = nuW[2], log = TRUE)) 
  # lp_curr = lp_curr + sum(logdinvgamma(phiC,a = 1, b = L_lasso*L_lasso/2))
  # lp_curr = lp_curr + dgamma(L_lasso*L_lasso,Alam, Blam, log = TRUE)

  # #shrinkage 
  # lp_curr = lp_curr +sum(dgamma(deltaWC[-1], alphaW[2], betaW[2], log = TRUE))

  # #uniqueness priors
  # lp_curr = lp_curr + sum(logdinvgamma(Sig2curr, a = Asig, b = Bsig))
  # lp_curr = lp_curr + sum(logdinvgamma(Psi2curr, a = Apsi, b = Bpsi))

  # lp_store[iter] = lp_curr
}
if(0){
  plot(storeL, type = 'l');abline(h = Alam/Blam, col = 'red', lwd = 2)


  var(storeL[-(1:5000)])

  Ltmp = storeL[-(1:5000)]
  hist(Ltmp, freq = 0); lines(density(rgamma(1e4,Alam,Blam)))
  plot(density(log(1/rexp(length(Ltmp), Ltmp^2/2))))
}

# plot(lp_store, type = 'l')


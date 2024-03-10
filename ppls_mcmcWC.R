# COMBINE THE TWO SHRINKAGES, i.e. DeltaC = DeltaW

########################## INITIALISE
set.seed(NULL)

# Qs = 1.5*Q
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


    Psi2curr = 1/rgamma(R,Apsi,Bpsi)
    Psi2curr = 1/rep(1e3,R)

    SIGMAinv = diag(1/Sig2curr)
    if(R>1){
      PSIinv   = diag(1/Psi2curr)
    } else{
      PSIinv = matrix(1/Psi2curr,1,1)
    }




    Sinv = Id(Qs) + t(Wcurr)%*%SIGMAinv %*%Wcurr
    S = solve(Sinv)



    Zcurr =  X %*%t(S%*% t(Wcurr)%*%SIGMAinv) 

    Ccurr = t(solve(t(Zcurr)%*%Zcurr)%*%t(Zcurr)%*%Y)
    # Ccurr = matrix(0,nrow = R,ncol = Qs)

    phiW= rgamma(P*Qs,nuW[1],nuW[2])
    # hist(rgamma(P*Qs,nuW[1],nuW[2]))
    phiC = rgamma(R*Qs,nuC[1],nuC[2])

    deltaWC = rep(1,Qs)
    if(shrinkage){
      deltaWC = c(1.0, rltrgamma(Qs-1,alphaW[2],betaW[2],truncWC))
    }

    tauWC = cumprod(deltaWC)

    DELTA = matrix(1/(phiW*tauWC), ncol = Qs, byrow = T)
    OMEGA = matrix(1/(phiC*tauWC), ncol = Qs, byrow = T)


  # Ccurr = sqrt(OMEGA) * rnorm(R*Qs);

    Bcurr = diag(c(rep(1, Qs)))#diag(c(rep(as.logical(!sparsity), Qs)))
  }

  

  ####

  sumBstore = rep(NA, N_MC)
  storeB = matrix(NA,nrow = N_MC,ncol = Qs)
  storeCB = array(NA, dim = c(N_MC,R,Qs))
  storeC = array(NA, dim = c(N_MC,R,Qs))
  storeCBa = array(NA, dim = c(N_MC,R,P))
  # storeW = array(NA, dim = c(N_MC,P,Qs))

  storeSig2 = matrix(NA,nrow = N_MC,ncol = P)
  storePsi2 = matrix(NA,nrow = N_MC, ncol = R)
  # storePhiW = array(NA, dim = c(N_MC,P,Qs))
  # storePhiC = array(NA, dim = c(N_MC,R,Qs))
  storeTauWC = matrix(NA,nrow = N_MC,ncol = Qs)
  # storeTauC = matrix(NA,nrow = N_MC, ncol = Qs)

  EYpred = array(NA, dim = c(N_pred,Ntest,R))
  Ypred = array(NA, dim = c(N_pred,Ntest,R))
  iter = 1

  probSS = 0.5
  LogInvOdds = log(1/probSS - 1)
}
######################### RUN CHAIN 
# cat("\f")

# print('initiation done')
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

  t_tmp = Sys.time()
  ZtZ = t(Zcurr) %*%Zcurr
  # mplot(abs(ZtZ))
  for(j in 1:P){
    choltDj = .chol(diag(1/DELTA[j,]) + SIGMAinv[j,j] * ZtZ)
    Wcurr[j,] = backsolve(choltDj,
      backsolve(choltDj,SIGMAinv[j,j]*t(Zcurr)%*%X[,j],k=Qs,transpose=1)+rnorm(Qs),
      k=Qs)

  }
  mcmc_time[2] = mcmc_time[2]+(Sys.time()-t_tmp)
  

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
  if(sparsity){

    # MH
    for(j in 1:Qs){
      gam = 1-2*Bcurr[j,j]

      tSUM_RESID = t(Y - Zcurr%*%Bcurr%*%t(Ccurr))

      delta = t(Ccurr[,j])%*%PSIinv%*% (Ccurr[,j]*ZtZ[j,j] - 2*gam*tSUM_RESID%*%Zcurr[,j])
      # 1/(1+exp(0.5*delta))
      # exp(-0.5*delta)

      if(runif(1)<as.numeric(exp(-0.5*delta + gam*LogInvOdds))){
        Bcurr[j,j] = Bcurr[j,j]+gam
      }
    }
    SSsucc = sum(Bcurr)
    probSS = rbeta(1, 1+SSsucc, 1+Qs-SSsucc)
    LogInvOdds = log(1/probSS - 1)
  }

  mcmc_time[4] = mcmc_time[4]+(Sys.time()-t_tmp)
  # update precisions



  t_tmp = Sys.time()
  if(T){
    # summary(as.numeric(Wcurr^2 %*%diag(tauW)))
    # par(mfrow = c(1,1))
    # mplot(nuW[2]+Ccurr^2 %*%diag(tauWC), digits = NA, border = 0)
    phiW = rgamma(P*Qs,nuW[1]+0.5,1)/(nuW[2]+0.5*Wcurr^2 %*%diag(tauWC))
    phiC = rgamma(R*Qs,nuC[1]+0.5,1)/(nuC[2]+0.5*Ccurr^2 %*%diag(tauWC))

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


      
    }
    DELTA = 1/ (phiW %*%diag(tauWC))
    OMEGA = 1/ (phiC %*%diag(tauWC))

  }
  mcmc_time[5] = mcmc_time[5]+(Sys.time()-t_tmp)
  # update noise variances

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

  # if(R>1){
  #   PSIinv   = diag(1/Psi2curr)
  # } else{
  #   PSIinv = matrix(1/Psi2curr,1,1)
  # }
  


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
      cholVarY = chol(CB %*% S %*% t(CB) + diag(Psi2curr))
    } else{
      cholVarY = chol(CB %*% S %*% t(CB) + matrix(Psi2curr,1,1))
    }


    for(i in 1:Ntest){


      EYpred[idx_pred,i,] = CBa %*% Xtest[i,]

      Ypred[idx_pred,i,] = EYpred[idx_pred,i,] + cholVarY %*%rnorm(R) 
    }
  }


  mcmc_time[7] = mcmc_time[7]+(Sys.time()-t_tmp)

  ### STORAGE
  # storeW[iter,,] = Wcurr
  storeC[iter,,] = Ccurr
  # storePhiW[iter,,] = phiW
  # storePhiC[iter,,] = phiC
  storeTauWC[iter,] = tauWC
  # storeTauC[iter,] = tauC

  storeCB[iter,,] = CB
  storeCBa[iter,,] = CBa

  storeB[iter,] = diag(Bcurr)

  storeSig2[iter,] = Sig2curr
  storePsi2[iter,] = Psi2curr
}

#************************************Luo JBS ******************************************;
#**************************************************************************************;
#** d: Input dataframe that must contain:                                            **;
#**      1) event.time:                             Event Time var (in days say)     **;
#**      2) censor.ind{1=Event, 0=Censored}:        Censoring Indicator              **;
#**      3) R{1=Intervention 0=Control}:            Treatment variable               **;
#**      4) switch.ind{1=Switched, 0=Not-switched}: Switching indicator called       **;
#**      5) Vt:                  Switch time (in days say), 0 if no switching        **;
#**      6) x1, x2:              covariates                                          **;
#**************************************************************************************; 
LuoJBS <- function(d, JBSmodInit, JBSmod, JBSfinal, pThreshM=0.9, pThreshI=0.9, MaxIter=50){

  # latent time after switching
  d$latent.time <- d$Survt0 <- with(d, event.time-Vt) 
  d$Rswitch.ind <- with(d, R*switch.ind)
  pMnew <- pInew <- 0
  Gamma <- NULL

  fit <- survreg(JBSmodInit, data=d, dist='weibull', scale=1)
  pMnew <- summary(fit)$table['switch.ind', 'p']
  if(pMnew > pThreshM) GammaEst <- rep(0,3)
  ind <- 1
  iter <- 1
  error <- FALSE
  if(pMnew <= pThreshM) {
    while (ind & iter<MaxIter & !error) {
      fit <- tryCatch(survreg(JBSmod, data=d, dist='weibull', scale=1, control = list(maxiter=MaxIter)), error=function(e) 1)
      error <- !is.list(fit)
      if(!error) error <- is.na(fit$scale)
      if(!error) {
        pMnew <- summary(fit)$table['switch.ind', 'p']
        pInew <- summary(fit)$table['Rswitch.ind', 'p']
        ind <- (pMnew <= pThreshM)*(pInew <= pThreshI)
        coefLoc <- match(c('Vt', 'switch.ind', 'Rswitch.ind'), names(fit$coefficients))
        d$latent.time[d$switch.ind==1] <- (d$Survt0*exp(-apply(d[,c('Vt', 'switch.ind', 'Rswitch.ind')]*matrix(fit$coef[coefLoc],nrow=n,ncol=3,byrow=T), 1, sum)))[d$switch.ind==1]
        Gamma <- rbind(Gamma, fit$coefficients[coefLoc])
        iter <- iter+1
        # print(fit$coefficients)
        # cat('iter=',iter, 'pMnew=', pMnew, 'pInew=', pInew, 'ind=', ind, '\n')
      }
    }
    if(ind==0) Gamma[nrow(Gamma),] <- 0
    GammaEst <- apply(Gamma, 2, sum)
  }
  if(iter<MaxIter & !error) {
    d$latent.time[d$switch.ind==1] <- (d$Survt0*exp(-apply(d[,c('Vt', 'switch.ind', 'Rswitch.ind')]*matrix(GammaEst,nrow=n,ncol=3,byrow=T), 1, sum))+d$Vt)[d$switch.ind==1]
    # survreg.out <- survreg(Surv(latent.time, event.ind) ~ R + x, data=d, dist='weibull', scale=1)
    survreg.out <- survreg(JBSfinal, data=d, dist='weibull', scale=1)
    converge <- TRUE
  } else{
    survreg.out <- NA
    converge <- FALSE
  }
    
  return(list(survreg.out=survreg.out, converge=converge))
}

###########################################################################################################################################

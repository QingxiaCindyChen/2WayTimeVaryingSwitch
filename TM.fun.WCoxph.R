TM.fun.WCoxph <- function(SdatCount, Sdat, Umod, TDmod, TUmod, TGmod, sd=FALSE, ID) 
{

  alpFit <- tryCatch(glm(Umod, family=binomial(link="logit"), data=Sdat), error=function(e){1}, warning=function(w){2})
  if(is.list(alpFit)) oldalpvec <- alpFit$coefficients
                        
  fit0 <- tryCatch(coxph(as.formula("TDmod"), data=subset(SdatCount, Ee<1), weight=1-Ee), error=function(e){1}, warning=function(w){2})
  if(is.list(fit0)) {
    oldbeta0vec <- fit0$coefficients
    tmp <- mybasehaz(fit0); oldh0 <- tmp$h0; oldtimeD <- tmp$time; oldParModD <- data.frame(timeD=oldtimeD, oldh0=oldh0)
  }
  
  fit1 <- tryCatch(coxph(as.formula("TUmod"), data=subset(Sdat, Ee>0), weight=Ee), error=function(e){1}, warning=function(w){2})
  if(is.list(fit1)) {
    oldbeta1vec <- fit1$coefficients
    tmp <- mybasehaz(fit1); oldh1 <- tmp$h0; oldtimeU <- tmp$time; oldParModU <- data.frame(timeU=oldtimeU, oldh1=oldh1)
  }
  
  fit2 <- tryCatch(coxph(as.formula("TGmod"), data=subset(SdatCount, Group%in%c(2,3))), error=function(e){1}, warning=function(w){2})
  if(is.list(fit2)) {
    oldbeta2vec <- fit2$coefficients
    tmp <- mybasehaz(fit2); oldh2 <- tmp$h0; oldtimeG <- tmp$time; oldParModG <- data.frame(timeG=oldtimeG, oldh2=oldh2)
  }
  
  error <- 0
  if(!is.list(alpFit)) error <- 1
  if(!is.list(fit0)) error <- 2
  if(!is.list(fit1)) error <- 3
  if(!is.list(fit2)) error <- 4
  
  if(error>0) return(list(newpara=NA, iter=NA, error=error))
    
  oldpara <- list(beta0vec=oldbeta0vec, h0=oldh0, timeD=oldtimeD, beta1vec=oldbeta1vec, h1=oldh1, timeU=oldtimeU, 
                  beta2vec=oldbeta2vec, h2=oldh2, timeG=oldtimeG, alpvec=oldalpvec)
  
  epsilon <- 0.001
  maxiter <- 200
  iter <- 0
  absdiff <- 1
  
  # EM algorithm
  while(absdiff > epsilon & iter < maxiter) {
    ret <- Estep3(oldpara=oldpara, SSdatCount=SdatCount, SSdat=Sdat, Umod, TDmod, TUmod)
    Sdat$Ee <- ret$Eed; SdatCount$Ee <- ret$EedC
    newpara <- Mstep3(oldpara=oldpara, SSdatCount=SdatCount, SSdat=Sdat, Umod, TDmod, TUmod, TGmod)
    # make sure the baseline hazard estimates are corresponding to the same times
    ParModD <- data.frame(timeD=newpara$timeD, h0new=newpara$h0)
    NewParModD <- merge(oldParModD, ParModD, all=T)
    ParModU <- data.frame(timeU=newpara$timeU, h1new=newpara$h1)
    NewParModU <- merge(oldParModU, ParModU, all=T)
    ParModG <- data.frame(timeG=newpara$timeG, h2new=newpara$h2)
    NewParModG <- merge(oldParModG, ParModG, all=T)
    
    diff <- 0
    for (k in c(1,4,7,10)) diff <- c(diff,max(abs(oldpara[[k]]-newpara[[k]])))
    diff <- with(NewParModD, c(diff, max(abs(h0new-oldh0),na.rm=T)))
    diff <- with(NewParModU, c(diff, max(abs(h1new-oldh1), na.rm=T)))
    diff <- with(NewParModG, c(diff, max(abs(h2new-oldh2), na.rm=T)))
    
    oldpara <- newpara
    oldParModD <- data.frame(timeD=NewParModD$timeD, oldh0=NewParModD$h0new)
    oldParModU <- data.frame(timeU=NewParModU$timeU, oldh1=NewParModU$h1new)
    oldParModG <- data.frame(timeG=NewParModG$timeG, oldh2=NewParModG$h2new)
    
    iter <- iter+1
    absdiff <- max(diff)
  }
  if(iter==maxiter) error <- 5
  
  oldpara$timeD <- oldpara$timeD[oldpara$h0>0] 
  oldpara$h0 <- oldpara$h0[oldpara$h0>0]
  oldpara$timeU <- oldpara$timeU[oldpara$h1>0] 
  oldpara$h1 <- oldpara$h1[oldpara$h1>0]
  oldpara$timeG <- oldpara$timeG[oldpara$h2>0] 
  oldpara$h2 <- oldpara$h2[oldpara$h2>0]
  oldpara$nbeta0 <- length(oldpara$beta0vec)
  oldpara$n0 <- length(oldpara$h0)
  oldpara$nbeta1 <- length(oldpara$beta1vec)
  oldpara$n1 <- length(oldpara$h1)
  oldpara$nbeta2 <- length(oldpara$beta2vec)
  oldpara$n2 <- length(oldpara$h2)
  oldpara$nalpha <- length(oldpara$alpvec)
  
  if (sd) {
    # variance estimation
    templist <- Covest(par=oldpara, Umod, TDmod, TUmod, TGmod, SdatCount, Sdat, ID=ID)
    Sigma <- templist[[1]]
    Influence <- templist[[2]]
    
    IND <- cumsum(c(oldpara$nbeta0, oldpara$n0, oldpara$nbeta1, oldpara$n1, oldpara$nbeta2, oldpara$n2, oldpara$nalpha))
    indbeta0 <- 1:IND[1];          indh0 <- (IND[1]+1):IND[2]
    indbeta1 <- (IND[2]+1):IND[3]; indh1 <- (IND[3]+1):IND[4]
    indbeta2 <- (IND[4]+1):IND[5]; indh2 <-(IND[5]+1):IND[6]
    indalpha <- (IND[6]+1):IND[7]
    
    Sigmah0 <- Sigma[indh0,indh0]
    Sigmah1 <- Sigma[indh1, indh1]
    Sigmah2 <- Sigma[indh2, indh2]
    
    beta0cov <- Sigma[indbeta0, indbeta0]
    beta1cov <- Sigma[indbeta1, indbeta1]
    beta2cov <- Sigma[indbeta2, indbeta2]
    alpcov <- Sigma[indalpha, indalpha]
    
    return(list(newpara=oldpara, iter=iter, error=error, EeCount=SdatCount$Ee, EeD=Sdat$Ee, beta0cov=beta0cov, beta1cov=beta1cov, beta2cov=beta2cov, alpcov=alpcov, Sigmah0=Sigmah0, Sigmah1=Sigmah1, Sigmah2=Sigmah2))
  }
  else return(list(newpara=oldpara, iter=iter, error=error, EeCount=SdatCount$Ee, EeD=Sdat$Ee))
  
}

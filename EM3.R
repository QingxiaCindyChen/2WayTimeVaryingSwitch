
##### EM ALGORITHMS using Cox PH procedure in M-step 
################################################################################
# E-step

Estep3 <- function(oldpara, SSdatCount, SSdat, Umod, TDmod, TUmod) {
  
  #oldpara: old parameters
  # datCount: whole dataset in counting process
  # dat: whole dataset with one record per patient
  # Umod: the model formular for U
  # TDmod: the model formular for TD
  # TUmod: the model formular for YU
  # (StartD, EndD, DeltaD): the counting process data for TD model (time-dependent); time interval (StartD, EndD]
  # (SurvU, DeltaU): the right censoring data format for TU model (time-independent)
  
  n <- dim(SSdatCount)[1]; nid <- dim(SSdat)[1]; 
  # X <- model.matrix(Umod, SSdat)
  X <- as.matrix(cbind(intercept=rep(1,nid), model.frame(Umod, data=SSdat, na.action="na.pass")[,-1]))
  Xd <- model.matrix(TDmod, SSdatCount)[,-1]
  Xe <- model.matrix(TUmod, SSdat)[,-1]
  temp <- as.matrix(model.frame(TDmod, data=SSdatCount, na.action="na.pass"))
  StartD <- temp[,1]; EndD <- temp[,2]; DeltaD <- temp[,3]
  temp <- as.matrix(model.frame(TUmod, data=SSdat, na.action="na.pass"))
  SurvU <- temp[,1]; DeltaU <- temp[,2]
  
  y0 <- oldpara$timeD # the jump time of h0 including censored 
  y1 <- oldpara$timeU # the jump time of h1 including censored
  
  n0 <- length(y0); n1 <- length(y1); 
  indH0 <- (matrix(EndD,n,n0)>=matrix(y0,n,n0,byrow=TRUE))*(matrix(StartD,n,n0)<matrix(y0,n,n0,byrow=TRUE))
  indH1 <- (matrix(SurvU,nid,n1)>=matrix(y1,nid,n1,byrow=TRUE))
  nv <- table(SSdatCount$PTID)
  Tmp <- rep(1, nv[1])
  for (i in 2:nid) Tmp <- bdiag(Tmp, rep(1, nv[i])) # create block matrix n*nid
  
  HExbeta0 <- rowSums(t(Tmp)%*%((as.vector(exp(Xd%*%oldpara$beta0vec))%o%oldpara$h0)*indH0)) #nid*1
  EXalp <- as.vector(exp(X%*%oldpara$alpvec)) # nid*1
  HExbeta1 <- as.vector((indH1%*%oldpara$h1)*(exp(Xe%*%oldpara$beta1vec))) # nid*1
    
  PE0 <- 1/(1+EXalp)
  PE1 <- 1-PE0
  PE1.g4 <- PE1*exp(-HExbeta1)/(PE0*exp(-HExbeta0)+PE1*exp(-HExbeta1))   
  SSdat$Ee <- 0*(SSdat$Group==1)+(SSdat$Group %in% c(2,3))+(SSdat$Group==4)*PE1.g4
  SSdatCount <- merge(SSdatCount[, -which(names(SSdatCount)=='Ee')], subset(SSdat, select=c("PTID", "Ee")), by='PTID')
  
  return(list(Eed=SSdat$Ee, EedC=SSdatCount$Ee))
}                                              

###################################################################################
# M-step

Mstep3 <- function(oldpara, SSdatCount, SSdat, Umod, TDmod, TUmod, TGmod) {

  # estimating alpha
  # Ee <- dat$Ee
  # PE1 <- 1/(1+exp(-X%*%oldpara$alpvec))  
  # temp2 <- solve((t(X)*matrix(PE1*(1-PE1),dim(X)[2],n,byrow=TRUE))%*%X)
  # oldnewalpvec <- oldpara$alpvec + temp2%*%(t(X)%*%(Ee-PE1))  
  datG4.1 <- datG4.2 <- subset(SSdat, Group==4)
  datG123 <- subset(SSdat, Group %in% c(1,2,3))
  datG4.1$U <- 0; datG4.1$uweight <- 1-datG4.1$Ee
  datG4.2$U <- 1; datG4.2$uweight <- datG4.2$Ee
  datG123$uweight <- 1
  datGlm <- as.data.frame(rbind(datG123, datG4.1, datG4.2))
  datGlm$U <- as.integer(datGlm$U)
  newalpvec <- glm(Umod, family=quasibinomial(link="logit"), data=datGlm, weights=uweight, start=oldpara$alpvec)$coefficients
  # use quasibinomial in stead of binomial to suppress 

  # estimating beta0
  # newest <- Mstepsub((Group==1), 1-Ee, Xd, oldpara$beta0vec, Y)
  # oldnewbeta0vec <- newest$beta
  # oldnewh0 <- newest$h
  # TDmod <- Surv(Y, DeltaD2) ~ Xd with DeltaD2=(Group==1)
  fit0 <- coxph(as.formula("TDmod"), data=subset(SSdatCount, Ee<1), weight=1-Ee, init=oldpara$beta0vec)
  newbeta0vec <- fit0$coefficients
  tmp <- mybasehaz(fit0); newh0 <- tmp$h0; newtimeD <- tmp$time

  # estimating beta1
  # newest <- Mstepsub((Group==2|Group==3), Ee, Xe, oldpara$beta1vec, W)
  # oldnewbeta1vec <- newest$beta
  # oldnewh1 <- newest$h
  # TUmod <- Surv(W, DeltaU) ~ Xe with DeltaU=(Group==2|Group==3)
  fit1 <- coxph(as.formula("TUmod"), data=subset(SSdat, Ee>0), weight=Ee, init=oldpara$beta1vec)
  newbeta1vec <- fit1$coefficients
  tmp <- mybasehaz(fit1); newh1 <- tmp$h0; newtimeU <- tmp$time
  
  # estimating beta2
  # newest <- Mstepsub((Group==2), (Group%in%c(2,3)), Xg, oldpara$beta2vec, Y-W)
  # oldnewbeta2vec <- newest$beta
  # oldnewh2 <- newest$h
  # TGmod <- Surv(Y-W, DeltaG) ~ Xg with DeltaG=(Group==2)
  fit2 <- coxph(as.formula("TGmod"), data=subset(SSdatCount, Group%in%c(2,3) & StartG<EndG), init=oldpara$beta2vec)
  newbeta2vec <- fit2$coefficients
  tmp <- mybasehaz(fit2); newh2 <- tmp$h0; newtimeG <- tmp$time
  
  newpara <- list(beta0vec=newbeta0vec, h0=newh0, timeD=newtimeD, beta1vec=newbeta1vec, h1=newh1, timeU=newtimeU, 
                  beta2vec=newbeta2vec, h2=newh2, timeG=newtimeG, alpvec=newalpvec)
  return(newpara)
}

#######################################################################
# mybasehaz <- function(cumhaz, time)
  # cumhaz: the matrix returned by basehaz() in survival 
  # with two columns: hazard (baseline cumulative hazard function) and time (unique times including censoring time)
  # time: the observed unique event time
# {
#  H0 <- cumhaz$hazard
#  h0 <- diff(c(0, H0))
#  h0t <- cumhaz$time[h0>0]
#  h0 <- h0[h0>0]
#  return(h0[rank(time)])
# }

mybasehaz <- function(fit) 
{
  tmp <- basehaz(fit, centered=FALSE)
  if(tmp$hazard[1]==0) {
    newh <- diff(tmp$hazard)
    newtime <- tmp$time[-1] 
  }
  else {
    newh <- diff(c(0,tmp$hazard))
    newtime <- tmp$time 
  }
  return(list(h0=newh, time=newtime))
}

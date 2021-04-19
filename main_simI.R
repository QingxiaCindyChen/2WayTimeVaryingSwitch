####Main PROGRAMS FOR SEMI-COMPETING CROSSOVER STUDY
### compare to existing approaches 
#################################################################################

rm(list = ls())
library(survival)
library(Matrix)

# setup your own working directory
# dir <- '~/Dropbox/Papers/UConn_Julie/Program/TD_model/'

source(paste(dir, 'CountDat.R', sep=''))
source(paste(dir, 'EM3.R', sep=''))
source(paste(dir, 'TM.fun.WCoxph.R', sep=''))
source(paste(dir, 'resid.WCoxph.R', sep=''))
source(paste(dir, 'Predict.Fast.Nosd.Rev.R', sep=''))
source(paste(dir, 'Covest.R', sep=''))
source(paste(dir, 'Simulation/TrueSurv4.R', sep=''))
source(paste(dir, 'Simulation/shaocox3.r', sep=''))
source(paste(dir, 'Simulation/LuoJBS.R', sep=''))
set.seed(321)

# set parameters
# Simulation II
beta0vec <- c(-0.5, 1, 0.2)
beta0Tvec <- c(-1.1, 0.73)  # coefficients for (V(t), R*V(t))
beta1vec <- c(-0.16, 1, 0.5)
beta2vec <- c(-0.1, 0.6, -0.5, 0.5, -0.4)
beta2Tvec <- c(-1, 0.5)
alpvec <- c(0, -0.5, 1, 0.5)
gam1vec <- c(-0.16, 0.75, 0.2)

# increase treatment effect to evaluate the power
beta0vec <- c(-1, 1, 0.2)
beta0Tvec <- c(-1.1, 0.73)  # coefficients for (V(t), R*V(t))
beta1vec <- c(-1, 1, 0.5)
beta2vec <- c(-1, 0.6, -0.5, 0.5, -0.4)
beta2Tvec <- c(-1, 0.5)
alpvec <- c(0, -0.5, 1, 0.5)
gam1vec <- c(-0.16, 0.75, 0.2)

tau <- 3
pred.t <- seq(0.01,tau,0.01)
TS<- TrueSurv(pred.t, n=100000, beta0vec, beta0Tvec=c(0,0), beta1vec, beta2vec, beta2Tvec=c(0,0), alpvec, gam1vec)

testpt0 <- pred.t
testpt1 <- pred.t
testpt2 <- pred.t
truecoeff <- c(beta0vec, beta0Tvec, testpt0/2, beta1vec, testpt1, beta2vec, beta2Tvec, testpt2/2, alpvec, TS[[1]], TS[[2]])
# truecoeff <- c(beta0vec, beta0Tvec, testpt0, beta1vec, exp(testpt1)-1, beta2vec, beta2Tvec, testpt2/2, alpvec, TS[[1]], TS[[2]])

# start simulations

nsim <- 100
n <- 400
B <- 2

parest <- matrix(NA, nsim, length(truecoeff))
sdest <-  parest
reject <- matrix(NA, nsim, 4)
g0 <- matrix(NA, B, length(pred.t))
g1 <- matrix(NA, B, length(pred.t))
errormat <- iter <- rep(NA, nsim)

ITTKM0 <- matrix(NA, nsim, length(pred.t))
ITTKM1 <- matrix(NA, nsim, length(pred.t))
ITTCox0 <- matrix(NA, nsim, length(pred.t))
ITTCox1 <- matrix(NA, nsim, length(pred.t))
NoCrosKM0 <- matrix(NA, nsim, length(pred.t))
NoCrosKM1 <- matrix(NA, nsim, length(pred.t))
Shaocox0 <- matrix(NA, nsim, length(pred.t))
JBS0 <- matrix(NA, nsim, length(pred.t))
JBS1 <- matrix(NA, nsim, length(pred.t))
JBSest <- JBS_converge <- rep(NA, nsim)

censE <- 0
censD <- 0
G <- rep(0,4)
cross.rate0 <- 0
cross.rate1 <- 0

# ZZ <- matrix(rnorm(1000*length(pred.t)), nrow=1000, ncol=length(pred.t))
for (m in 1:nsim) {
  cat("sim=", m)
  print(Sys.time())

  dat0 <- Simulate(n, beta0vec, beta0Tvec, beta1vec, beta2vec, beta2Tvec, alpvec, gam1vec, tau)
   
  ############ from now on, use only dat0 ################
  cross.rate0 <- cross.rate0 + with(subset(dat0, R==0), sum(is.finite(Vtime))/length(Vtime)) 
  cross.rate1 <- cross.rate1 + with(subset(dat0, R==1), sum(is.finite(Vtime))/length(Vtime)) 
  
  dat0$U <- with(dat0, ifelse(Group==1, 0, ifelse(Group %in% c(2,3), 1, NA)))  
  dat0$Ee <- with(dat0, ifelse(Group==1, 0, ifelse(Group %in% c(2,3), 1, 0.5))); 
  dat0$PTID <- dat0$ID
  dat0 <- subset(dat0, select=c("PTID", "R", "x1", "x2", "Y", "DeltaY", "W", "DeltaE", "Vtime", "Group", "U", "Ee", "z1"))
  epsilon1 <- min(with(dat0, ifelse(is.finite(Vtime) & Vtime>W, Vtime-W, Vtime)))
  epsilon2 <- min(with(dat0, ifelse(is.finite(Vtime) & Y>W, Y-W, Y)))
  epsilon <- min(with(dat0, min(Y, W, Vtime, na.rm=T)), epsilon1, epsilon2)*0.25
  # datCount0 <- CountDat(dat0$Y, dat0$DeltaD, dat0$W, dat0$DeltaE, as.matrix(dat0[, -which(names(dat0)=='Vtime')]), dat0$Vtime, epsilon)
  datCount0 <- CountDat(dat0$Y, dat0$DeltaY, dat0$W, dat0$DeltaE, as.matrix(dat0), dat0$Vtime, epsilon)
  datCount0$RVt <- with(datCount0, R*Vt)
  # pred.t <- seq(0, 1, length.out=200)*max(dat0$W)
  
  Umod <- U ~ R + x1 + x2
  TDmod <- Surv(StartD, EndD, DeltaD) ~ R + x1 + x2 + Vt + RVt # true model
  # TDmod <- Surv(StartD, EndD, DeltaD) ~ R + x2 + Vt + RVt # misspecified model 
  TUmod <- Surv(W, DeltaE) ~ R + x1 + x2
  TGmod <- Surv(StartG, EndG, DeltaG) ~ R + x1 + x2 + z1 + W + Vt + RVt 

  ret <- TM.fun.WCoxph(datCount0, dat0, Umod, TDmod, TUmod, TGmod, sd=T, ID='PTID')
    
  errormat[m] <- ret$error
  iter[m] <- ret$iter
  if(ret$error ==0) {
    par <- ret$newpara
    beta0est <- par$beta0vec; beta1est <- par$beta1vec; beta2est <- par$beta2vec; alpest<-par$alpvec  
    beta0sd <- sqrt(diag(ret$beta0cov)); beta1sd <- sqrt(diag(ret$beta1cov)); beta2sd <- sqrt(diag(ret$beta2cov)); alpsd <- sqrt(diag(ret$alpcov));
    n0 <- length(par$timeD); n1 <- length(par$timeU); n2 <- length(par$timeG)

    testy0 <- matrix(testpt0,length(testpt0),n0,byrow=F)>=matrix(par$timeD, length(testpt0), n0,byrow=T)
    estH0 <- testy0%*%par$h0
    sdH0 <- sqrt(diag(testy0%*%ret$Sigmah0%*%t(testy0)))

    testy1 <- matrix(testpt1,length(testpt1),n1,byrow=F)>=matrix(par$timeU, length(testpt1),n1,byrow=T)
    estH1 <- testy1%*%par$h1
    sdH1 <- sqrt(diag(testy1%*%ret$Sigmah1%*%t(testy1)))

    testy2 <- matrix(testpt2,length(testpt2),n2,byrow=F)>=matrix(par$timeG, length(testpt2),n2,byrow=T)
    estH2 <- testy2%*%par$h2
    sdH2 <- sqrt(diag(testy2%*%ret$Sigmah2%*%t(testy2)))

  n <- dim(dat0)[1]
  X <- as.matrix(cbind(intercept=rep(1,n), model.frame(Umod, data=dat0, na.action="na.pass")[,-1]))
  Xd <- as.matrix(model.frame(TDmod, data=datCount0, na.action="na.pass"))[datCount0$Last==1,-(1:3)]
  Xe <- model.matrix(TUmod, dat0)[,-1]
  Xg <- as.matrix(model.frame(TGmod, data=datCount0, na.action="na.pass"))[datCount0$Last==1,-(1:3)]
  Xbaseline <- X[, -c(1:2)]
  Xd[,'Vt'] <- Xd[,'RVt'] <- Xg[,'Vt'] <- Xg[,'RVt'] <- 0  # prediction assuming no switching

   g0e <- Predict.Fast.Nosd(pred.t=pred.t, status=0, para=par, A0=dat0$R, Group=dat0$Group, X=X, Xd=Xd, Xe=Xe, 
                           Xg=Xg, Xbaseline=Xbaseline)$surv.pred 
   g1e <- Predict.Fast.Nosd(pred.t=pred.t, status=1, para=par, A0=dat0$R, Group=dat0$Group, X=X, Xd=Xd, Xe=Xe, 
                           Xg=Xg, Xbaseline=Xbaseline)$surv.pred 

### start of bootstrap sampling to estimate standard deviation.      
 for (b in 1:B) {
   cat("m=", m, "b=", b)
   print(Sys.time())
   sap <- sample(1:n, n, replace=T)
   Sdat <- dat0[sap,]
   Sdat$PTID <- 1:n # pseudo unique id
   SdatCount <- with(Sdat, CountDat(Y, DeltaY, W, DeltaE, as.matrix(Sdat), Vtime, epsilon))
   SdatCount$RVt <- with(SdatCount, R*Vt)
   retB <- TM.fun.WCoxph(SdatCount=SdatCount, Sdat=Sdat, Umod, TDmod, TUmod, TGmod, sd=FALSE, ID='PTID')
   g0[b,] <- Predict.Fast.Nosd(pred.t=pred.t, status=0, para=retB$newpara, A0=Sdat$R, Group=Sdat$Group, X=X[sap,], Xd=Xd[sap,], Xe=Xe[sap,],
                               Xg=Xg[sap,], Xbaseline=Xbaseline[sap,], z=Sdat$z1)$surv.pred
   g1[b,] <- Predict.Fast.Nosd(pred.t=pred.t, status=1, para=retB$newpara, A0=Sdat$R, Group=Sdat$Group, X=X[sap,], Xd=Xd[sap,], Xe=Xe[sap,],
                               Xg=Xg[sap,], Xbaseline=Xbaseline[sap,], z=Sdat$z1)$surv.pred
 }
 g0sd <- apply(g0,2,sd)
 g1sd <- apply(g1,2,sd)
### End of Bootstrap Loop
    
   parest[m,] <- c(beta0est, estH0, beta1est, estH1, beta2est, estH2, alpest, g0e, g1e)
   sdest[m,] <- c(beta0sd, sdH0, beta1sd, sdH1, beta2sd, sdH2, alpsd, g0sd, g1sd)

   }else {
    parest[m,] <- sdest[m,] <- NA
  }
     # summary
     censE <- with(dat0, censE+sum(DeltaE[Group!=1])/n)  
     censD <- with(dat0, censD+sum(DeltaY)/n)
     G <- with(dat0, G+table(Group)/n)
     summary <- cbind(censE, censD, G[1], G[2], G[3], G[4], cross.rate0, cross.rate1) 
  
  # ITT with KM
  fitITT <- coef(coxph(Surv(Y, DeltaY) ~ R, data = dat0))[1]
  fit0 <- survfit(Surv(Y, DeltaY) ~ 1, data = subset(dat0, R==0))
  fit1 <- survfit(Surv(Y, DeltaY) ~ 1, data = subset(dat0, R==1))
  t0 <- fit0$time
  L0 <- length(t0)
  t1 <- fit1$time
  L1 <- length(t1)
  P <- length(pred.t)
  tmpmat0 <- (matrix(pred.t, nrow=P, ncol=L0, byrow=FALSE) >= matrix(t0, nrow=P, ncol=L0, byrow=TRUE)) 
  f0 <- diff(c(fit0$surv[L0:1],1), lag=1)[L0:1]
  ITTKM0[m,] <- 1-tmpmat0%*%f0
  tmpmat1 <- (matrix(pred.t, nrow=P, ncol=L1, byrow=FALSE) >= matrix(t1, nrow=P, ncol=L1, byrow=TRUE)) 
  f1 <- diff(c(fit1$surv[L1:1],1), lag=1)[L1:1]
  ITTKM1[m,] <- 1-tmpmat1%*%f1
  
  # No crossover with KM
  fit0 <- survfit(Surv(Y, DeltaY) ~ 1, data = subset(dat0, R==0 & !is.finite(Vtime)))
  fit1 <- survfit(Surv(Y, DeltaY) ~ 1, data = subset(dat0, R==1 & !is.finite(Vtime)))
  t0 <- fit0$time
  L0 <- length(t0)
  t1 <- fit1$time
  L1 <- length(t1)
  P <- length(pred.t)
  tmpmat0 <- (matrix(pred.t, nrow=P, ncol=L0, byrow=FALSE) >= matrix(t0, nrow=P, ncol=L0, byrow=TRUE)) 
  f0 <- diff(c(fit0$surv[L0:1],1), lag=1)[L0:1]
  NoCrosKM0[m,] <- 1-tmpmat0%*%f0
  tmpmat1 <- (matrix(pred.t, nrow=P, ncol=L1, byrow=FALSE) >= matrix(t1, nrow=P, ncol=L1, byrow=TRUE)) 
  f1 <- diff(c(fit1$surv[L1:1],1), lag=1)[L1:1]
  NoCrosKM1[m,] <- 1-tmpmat1%*%f1

  # Luo JBS method
  dat0$event.time <- dat0$Y; dat0$event.ind <- dat0$DeltaY; dat0$trt.ind <- dat0$R
  dat0$switch.ind <- with(dat0, ifelse(is.finite(Vtime), 1, 0))
  dat0$Vt <- with(dat0, ifelse(is.finite(Vtime), Vtime, 0))
  JBSmodInit <- Surv(latent.time, event.ind) ~ R + x1 + x2 + Vt + switch.ind
  JBSmod <- Surv(latent.time, event.ind) ~ R +  x1 + x2 + Vt + switch.ind + Rswitch.ind
  JBSfinal <- Surv(latent.time, event.ind) ~ R + x1 + x2
  pThreshM <- 0.80; pThreshI <- 0.80; MaxIter <- 200
  JBS <- LuoJBS(d=dat0, JBSmodInit, JBSmod, JBSfinal, pThreshM, pThreshI, MaxIter=MaxIter)
  JBS_result <- JBS$survreg.out
  JBS_converge[m] <- JBS$converge 
  if(JBS_converge[m]) {
    JBSest[m] <- JBS_result$coef['R']
    JBS0[m,] <- pweibull(pred.t, scale=exp(JBS_result$coef[1]), shape=1/JBS_result$scale, lower.tail=FALSE)
    JBS1[m,] <- pweibull(pred.t*exp(-JBS_result$coef['R']/JBS_result$scale), scale=exp(JBS_result$coef[1]), shape=1/JBS_result$scale, lower.tail=FALSE)
  } 
  else JBS0[m,] <- JBS1[m,] <- NA
  
  # save(summary, parest, file='simI_result_102016.RData')
  print(Sys.time())  
}
 
save(list=ls(all=TRUE), file="working.RData")
save(censE, censD, G, parest, file='simI_result_102016.RData')

summary/nsim
pt0 <- c(16, 43, 99)
pt1 <- c(11, 29, 63)
pt2 <- c(20, 49, 99)
Lp <- length(estH0) 
L0 <- length(beta0est)+length(beta1est)+length(beta2est)+3*length(estH0)+length(alpest)
sel <- which(pred.t %in% c(0.75, 1.5, 2.25, 3))
est <- apply(parest[, c(1:5, 5+pt0, (6+Lp):(8+Lp), (8+Lp+pt1), (9+2*Lp):(15+2*Lp), (15+Lp+pt2), 
                        (16+3*Lp):(19+3*Lp), L0+sel, L0+Lp+sel)], 2, mean, na.rm=T)
EmpSd <- apply(parest[, c(1:5, 5+pt0, (6+Lp):(8+Lp), (8+Lp+pt1), (9+2*Lp):(15+2*Lp), (15+Lp+pt2), 
                        (16+3*Lp):(19+3*Lp), L0+sel, L0+Lp+sel)], 2, sd, na.rm=T)
MeanSd <- apply(sdest[, c(1:5, 5+pt0, (6+Lp):(8+Lp), (8+Lp+pt1), (9+2*Lp):(15+2*Lp), (15+Lp+pt2), 
                          (16+3*Lp):(19+3*Lp), L0+sel, L0+Lp+sel)], 2, mean, na.rm=T)
truecoeff1 <- truecoeff[c(1:5, 5+pt0, (6+Lp):(8+Lp), (8+Lp+pt1), (9+2*Lp):(15+2*Lp), (15+Lp+pt2), 
                        (16+3*Lp):(19+3*Lp), L0+sel, L0+Lp+sel)]
result <- as.matrix(cbind(truecoeff1, est, est-truecoeff1, EmpSd, MeanSd))
colnames(result) <- c('true','est','bias', 'EmpSd', "MeanSd")
rownames(result) <- c(paste("beta0",1:5, sep=''), 'H0(0.16)', 'H0(0.43)', 'H0(0.99)', 
                      paste("beta1",1:3, sep=''), 'H1(0.11)', 'H1(0.29)', 'H1(0.63)',
                      paste("beta2",1:7, sep=''), 'H2(0.20)', 'H2(0.49)', 'H2(0.99)',
                      paste("alp0",1:4, sep=''), 'S0(0.75)', 'S0(1.5)', 'S0(2.25)', 'S0(3)', 
                      'S1(0.75)', 'S1(1.5)', 'S1(2.25)', 'S1(3)')
result

g0 <- parest[,(L0+1):(L0+Lp)]
g1 <- parest[,(L0+Lp+1):(L0+2*Lp)]

mse <- matrix(NA,2,5)
mse[1,1] <- mean(apply((ITTKM0-matrix(TS[[1]], nsim, length(pred.t), byrow=TRUE))^2, 1, mean))
mse[2,1] <- mean(apply((ITTKM1-matrix(TS[[2]], nsim, length(pred.t), byrow=TRUE))^2, 1, mean))
mse[1,2] <- mean(apply((ITTCox0-matrix(TS[[1]], nsim, length(pred.t), byrow=TRUE))^2, 1, mean))
mse[2,2] <- mean(apply((ITTCox1-matrix(TS[[2]], nsim, length(pred.t), byrow=TRUE))^2, 1, mean))
mse[1,3] <- mean(apply((NoCrosKM0-matrix(TS[[1]], nsim, length(pred.t), byrow=TRUE))^2, 1, mean))
mse[2,3] <- mean(apply((NoCrosKM1-matrix(TS[[2]], nsim, length(pred.t), byrow=TRUE))^2, 1, mean))
mse[1,4] <- mean(apply((JBS0-matrix(TS[[1]], nsim, length(pred.t), byrow=TRUE))^2, 1, mean, na.rm=T),na.rm=T)
mse[2,4] <- mean(apply((JBS1-matrix(TS[[2]], nsim, length(pred.t), byrow=TRUE))^2, 1, mean, na.rm=T), na.rm=T)
mse[1,5] <- mean(apply((g0-matrix(TS[[1]], nsim, length(pred.t), byrow=TRUE))^2, 1, mean, na.rm=T), na.rm=T)
mse[2,5] <- mean(apply((g1-matrix(TS[[2]], nsim, length(pred.t), byrow=TRUE))^2, 1, mean, na.rm=T), na.rm=T)
colnames(mse) <- c("ITTKM", "ITTCox", "NoCrossKM", "JBS", "TM")
rownames(mse) <- c("control", "treatment")
sqrt(mse)

mmad <- matrix(NA,2,5)
mmad[1,1] <- mean(apply(abs(ITTKM0-matrix(TS[[1]], nsim, length(pred.t), byrow=TRUE)), 1, max))
mmad[2,1] <- mean(apply(abs(ITTKM1-matrix(TS[[2]], nsim, length(pred.t), byrow=TRUE)), 1, max))
mmad[1,2] <- mean(apply(abs(ITTCox0-matrix(TS[[1]], nsim, length(pred.t), byrow=TRUE)), 1, max))
mmad[2,2] <- mean(apply(abs(ITTCox1-matrix(TS[[2]], nsim, length(pred.t), byrow=TRUE)), 1, max))
mmad[1,3] <- mean(apply(abs(NoCrosKM0-matrix(TS[[1]], nsim, length(pred.t), byrow=TRUE)), 1, max))
mmad[2,3] <- mean(apply(abs(NoCrosKM1-matrix(TS[[2]], nsim, length(pred.t), byrow=TRUE)), 1, max))
mmad[1,4] <- mean(apply(abs(JBS0-matrix(TS[[1]], nsim, length(pred.t), byrow=TRUE)), 1, max))
mmad[2,4] <- mean(apply(abs(JBS1-matrix(TS[[2]], nsim, length(pred.t), byrow=TRUE)), 1, max))
mmad[1,5] <- mean(apply(abs(g0-matrix(TS[[1]], nsim, length(pred.t), byrow=TRUE)), 1, max), na.rm=T)
mmad[2,5] <- mean(apply(abs(g1-matrix(TS[[2]], nsim, length(pred.t), byrow=TRUE)), 1, max), na.rm=T)
colnames(mmad) <- c("ITTKM", "ITTCox", "NoCrossKM", "JBS", "TM")
rownames(mmad) <- c("control", "treatment")
mmad

pred.tn <- c(0, pred.t)
surv.mat <- matrix(NA, length(pred.tn), 13)
surv.mat[,1] <- pred.tn
surv.mat[,2] <- c(1,TS[[1]])
surv.mat[,3] <- c(1,TS[[2]])
surv.mat[,4] <- c(1,apply(ITTKM0, 2, mean, na.rm=TRUE))
surv.mat[,5] <- c(1,apply(ITTKM1, 2, mean, na.rm=TRUE))
surv.mat[,6] <- c(1,apply(ITTCox0, 2, mean, na.rm=TRUE))
surv.mat[,7] <- c(1,apply(ITTCox1, 2, mean, na.rm=TRUE))
surv.mat[,8] <- c(1,apply(NoCrosKM0, 2, mean, na.rm=TRUE))
surv.mat[,9] <- c(1,apply(NoCrosKM1, 2, mean, na.rm=TRUE))
surv.mat[,10] <- c(1,apply(JBS0, 2, mean, na.rm=TRUE))
surv.mat[,11] <- c(1,apply(JBS1, 2, mean, na.rm=TRUE))
surv.mat[,12] <- c(1,apply(g0, 2, mean, na.rm=TRUE))
surv.mat[,13] <- c(1,apply(g1, 2, mean, na.rm=TRUE))
surv.mat <- as.data.frame(surv.mat)
names(surv.mat) <- c("time", "True0", "True1", "ITTKM0", "ITTKM1", "ITTCox0", "ITTCox1", "NoCrosKM0", "NoCrosKM1", 
                     "JBS0", "JBS1", "TM0", "TM1")
save(surv.mat, file="SIMI_surv2mat.Rdata")

pdf("Surv_All_09212016.pdf", height = 6, width = 9)
par(mfrow=c(1,1))
plot(x=surv.mat$time, y=surv.mat$True0, type='n', xlim=c(0,3), ylim=c(0,1), lty=1, col=1, xlab="Time", ylab="Survival Probability")
lines(x=surv.mat$time, y=surv.mat$True0, xlim=c(0,3), type="s", lty=1, col=1, lwd=1)
lines(x=surv.mat$time, y=surv.mat$True1, xlim=c(0,3), type="s",lty=2, col=1, lwd=1)
lines(x=surv.mat$time, y=surv.mat$ITTKM0, xlim=c(0,3), type="s", lty=1, col=2, lwd=0.5)
lines(x=surv.mat$time, y=surv.mat$ITTKM1, xlim=c(0,3), type="s",lty=2, col=2, lwd=0.5)
lines(x=surv.mat$time, y=surv.mat$JBS0, xlim=c(0,3), type="s", lty=1, col=5, lwd=0.5)
lines(x=surv.mat$time, y=surv.mat$JBS1, xlim=c(0,3), type="s",lty=2, col=5, lwd=0.5)
lines(x=surv.mat$time, y=surv.mat$TM0, xlim=c(0,3), type="s", lty=1, col=6, lwd=0.5)
lines(x=surv.mat$time, y=surv.mat$TM1, xlim=c(0,3), type="s",lty=2, col=6, lwd=0.5)
#legend("topright", c("True S(t|A=0)", "True S(t|A=1)", "ITT KM S(t|A=0)", "ITT KM S(t|A=1)", "ITT Cox S(t|A=0)", "ITT Cox S(t|A=1)", "IPE S(t|A=0)", "IPE S(t|A=1)", "ShaoCOX S(t|A=0)", "ShaoCOX S(t|A=1)", "TM S(t|A=0)", 
#    "TM S(t|A=1)"), lty=rep(c(1,2),times=6), col=rep(1:6, each=2), cex=0.7)
legend("topright", c("True S(t|R=0)", "True S(t|R=1)", "ITT KM S(t|R=0)", "ITT KM S(t|R=1)", "JBS S(t|R=0)", "JBS S(t|R=1)", "TM S(t|R=0)", 
                     "TM S(t|R=1)"), lty=rep(c(1,2),times=4), col=rep(c(1:2,5:6), each=2), cex=0.7)
dev.off()

pdf("Surv_All_09212016_col.pdf", height = 6, width = 9)
par(mfrow=c(2,2), mar=c(4,4,2,2)) 
plot(x=surv.mat$time, y=surv.mat$True0, type='n', xlim=c(0,3), ylim=c(0,1), lty=1, col=1, xlab="Time", ylab="Survival Probability")
lines(x=surv.mat$time, y=surv.mat$True0, xlim=c(0,3), type="s", lty=1, col=1, lwd=1)
lines(x=surv.mat$time, y=surv.mat$True1, xlim=c(0,3), type="s",lty=2, col=1, lwd=1)
lines(x=surv.mat$time, y=surv.mat$ITTKM0, xlim=c(0,3), type="s", lty=1, col=2, lwd=0.5)
lines(x=surv.mat$time, y=surv.mat$ITTKM1, xlim=c(0,3), type="s",lty=2, col=2, lwd=0.5)
legend("topright", c("True S(t|R=0)", "True S(t|R=1)", "ITT KM S(t|R=0)", "ITT KM S(t|R=1)"), bty='n', lty=rep(c(1,2),times=6), col=rep(c(1,2), each=2), cex=0.7)

plot(x=surv.mat$time, y=surv.mat$True0, type='n', xlim=c(0,3), ylim=c(0,1), lty=1, col=1, xlab="Time", ylab="Survival Probability")
lines(x=surv.mat$time, y=surv.mat$True0, xlim=c(0,3), type="s", lty=1, col=1, lwd=1)
lines(x=surv.mat$time, y=surv.mat$True1, xlim=c(0,3), type="s",lty=2, col=1, lwd=1)
lines(x=surv.mat$time, y=surv.mat$NoCrosKM0, xlim=c(0,3), type="s", lty=1, col=2, lwd=0.5)
lines(x=surv.mat$time, y=surv.mat$NoCrosKM1, xlim=c(0,3), type="s",lty=2, col=2, lwd=0.5)
legend("topright", c("True S(t|R=0)", "True S(t|R=1)", "No Crossover KM S(t|R=0)", "No Crossover KM S(t|R=1)"), bty='n', lty=rep(c(1,2),times=6), col=rep(c(1,2), each=2), cex=0.7)

plot(x=surv.mat$time, y=surv.mat$True0, type='n', xlim=c(0,3), ylim=c(0,1), lty=1, col=1, xlab="Time", ylab="Survival Probability")
lines(x=surv.mat$time, y=surv.mat$True0, xlim=c(0,3), type="s", lty=1, col=1, lwd=1)
lines(x=surv.mat$time, y=surv.mat$True1, xlim=c(0,3), type="s",lty=2, col=1, lwd=1)
lines(x=surv.mat$time, y=surv.mat$JBS0, xlim=c(0,3), type="s", lty=1, col=2, lwd=0.5)
lines(x=surv.mat$time, y=surv.mat$JBS1, xlim=c(0,3), type="s",lty=2, col=2, lwd=0.5)
legend("topright", c("True S(t|R=0)", "True S(t|R=1)", "JBS S(t|R=0)", "JBS S(t|R=1)"), bty='n', lty=rep(c(1,2),times=6), col=rep(c(1,2), each=2), cex=0.7)

plot(x=surv.mat$time, y=surv.mat$True0, type='n', xlim=c(0,3), ylim=c(0,1), lty=1, col=1, xlab="Time", ylab="Survival Probability")
lines(x=surv.mat$time, y=surv.mat$True0, xlim=c(0,3), type="s", lty=1, col=1, lwd=1)
lines(x=surv.mat$time, y=surv.mat$True1, xlim=c(0,3), type="s",lty=2, col=1, lwd=1)
lines(x=surv.mat$time, y=surv.mat$TM0, xlim=c(0,3), type="s", lty=1, col=2, lwd=0.5)
lines(x=surv.mat$time, y=surv.mat$TM1, xlim=c(0,3), type="s",lty=2, col=2, lwd=0.5)
legend("topright", c("True S(t|R=0)", "True S(t|R=1)", "TM S(t|R=0)", "TM S(t|R=1)"), bty='n', lty=rep(c(1,2),times=6), col=rep(c(1,2), each=2), cex=0.7)
dev.off()

pdf("Surv_All_09212016_bw.pdf", height = 6, width = 9)
par(mfrow=c(2,2), mar=c(4,4,2,2)) 
plot(x=surv.mat$time, y=surv.mat$True0, type='n', xlim=c(0,3), ylim=c(0,1), lty=1, col=1, xlab="Time", ylab="Survival Probability")
lines(x=surv.mat$time, y=surv.mat$True0, xlim=c(0,3), type="s", lty=1, col='grey', lwd=1)
lines(x=surv.mat$time, y=surv.mat$True1, xlim=c(0,3), type="s",lty=2, col='grey', lwd=1)
lines(x=surv.mat$time, y=surv.mat$ITTKM0, xlim=c(0,3), type="s", lty=1, col='black', lwd=0.5)
lines(x=surv.mat$time, y=surv.mat$ITTKM1, xlim=c(0,3), type="s",lty=2, col='black', lwd=0.5)
legend("topright", c("True S(t|R=0)", "True S(t|R=1)", "ITT KM S(t|R=0)", "ITT KM S(t|R=1)"), bty='n', lty=rep(c(1,2),times=6), col=rep(c('grey', 'black'), each=2), cex=0.7)

plot(x=surv.mat$time, y=surv.mat$True0, type='n', xlim=c(0,3), ylim=c(0,1), lty=1, col=1, xlab="Time", ylab="Survival Probability")
lines(x=surv.mat$time, y=surv.mat$True0, xlim=c(0,3), type="s", lty=1, col='grey', lwd=1)
lines(x=surv.mat$time, y=surv.mat$True1, xlim=c(0,3), type="s",lty=2, col='grey', lwd=1)
lines(x=surv.mat$time, y=surv.mat$NoCrosKM0, xlim=c(0,3), type="s", lty=1, col='black', lwd=0.5)
lines(x=surv.mat$time, y=surv.mat$NoCrosKM1, xlim=c(0,3), type="s",lty=2, col='black', lwd=0.5)
legend("topright", c("True S(t|R=0)", "True S(t|R=1)", "No Crossover KM S(t|R=0)", "No Crossover KM S(t|R=1)"), bty='n',lty=rep(c(1,2),times=6), col=rep(c('grey', 'black'), each=2), cex=0.7)

plot(x=surv.mat$time, y=surv.mat$True0, type='n', xlim=c(0,3), ylim=c(0,1), lty=1, col=1, xlab="Time", ylab="Survival Probability")
lines(x=surv.mat$time, y=surv.mat$True0, xlim=c(0,3), type="s", lty=1, col='grey', lwd=1)
lines(x=surv.mat$time, y=surv.mat$True1, xlim=c(0,3), type="s",lty=2, col='grey', lwd=1)
lines(x=surv.mat$time, y=surv.mat$JBS0, xlim=c(0,3), type="s", lty=1, col='black', lwd=0.5)
lines(x=surv.mat$time, y=surv.mat$JBS1, xlim=c(0,3), type="s",lty=2, col='black', lwd=0.5)
legend("topright", c("True S(t|R=0)", "True S(t|R=1)", "JBS S(t|R=0)", "JBS S(t|R=1)"), bty='n', lty=rep(c(1,2),times=6), col=rep(c('grey', 'black'), each=2), cex=0.7)

plot(x=surv.mat$time, y=surv.mat$True0, type='n', xlim=c(0,3), ylim=c(0,1), lty=1, col=1, xlab="Time", ylab="Survival Probability")
lines(x=surv.mat$time, y=surv.mat$True0, xlim=c(0,3), type="s", lty=1, col='grey', lwd=1)
lines(x=surv.mat$time, y=surv.mat$True1, xlim=c(0,3), type="s",lty=2, col='grey', lwd=1)
lines(x=surv.mat$time, y=surv.mat$TM0, xlim=c(0,3), type="s", lty=1, col='black', lwd=0.5)
lines(x=surv.mat$time, y=surv.mat$TM1, xlim=c(0,3), type="s",lty=2, col='black', lwd=0.5)
legend("topright", c("True S(t|R=0)", "True S(t|R=1)", "TM S(t|R=0)", "TM S(t|R=1)"), bty='n', lty=rep(c(1,2),times=6), col=rep(c('grey', 'black'), each=2), cex=0.7)
dev.off()

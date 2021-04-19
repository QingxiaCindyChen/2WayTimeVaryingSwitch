####Main PROGRAMS with true function simulated from JBS
### compare to existing approaches 
#################################################################################

rm(list = ls())
library(survival)
# library(rms)
# library(Hmisc)
library(Matrix)

dir <- '~/Dropbox/Papers/UConn_Julie/Program/TD_model/'
# dir <- 'C:/Users/chenq3/Dropbox/Papers/UConn_Julie/Program/TD_model/'

source(paste(dir, 'CountDat.R', sep=''))
source(paste(dir, 'EM3.R', sep=''))
source(paste(dir, 'TM.fun.WCoxph.R', sep=''))
source(paste(dir, 'resid.WCoxph.R', sep=''))
source(paste(dir, 'Predict.Fast.Nosd.Rev.R', sep=''))
source(paste(dir, 'Covest.R', sep=''))
source(paste(dir, 'SimulateJBS/TrueSurvJBS.R', sep=''))
source(paste(dir, 'Simulation/LuoJBS.R', sep=''))

set.seed(321)

# set parameters

# Simulation JBS
r<-0.5; v<-0.5; c<-0.5; a<-2; b<-4; xi<-0.1; lambda<--log(0.75); eta<-c(0.5,1.5); gamma<-c(-log(0.3),-log(1.0),-log(0.2))

tau <- b
pred.t <- seq(0.01,tau,0.01)
TS<- TrueSurv(pred.t, n=100000, r, v, xi, lambda)
testpt0 <- testpt1 <- testpt2 <- pred.t
#save(TS, file="TS.Rdata")

# start simulations
nsim <- 1000
n <- 1000
B <- 2

parest <- matrix(NA, nsim, length(pred.t)*5+16)
g0 <- matrix(NA, B, length(pred.t))
g1 <- matrix(NA, B, length(pred.t))
errormat <- iter <- rep(NA, nsim)

ITTKM0 <- matrix(NA, nsim, length(pred.t))
ITTKM1 <- matrix(NA, nsim, length(pred.t))
NoCrosKM0 <- matrix(NA, nsim, length(pred.t))
NoCrosKM1 <- matrix(NA, nsim, length(pred.t))
JBS0 <- matrix(NA, nsim, length(pred.t))
JBS1 <- matrix(NA, nsim, length(pred.t))
JBS_converge <- rep(NA, nsim); JBSest <- rep(NA, nsim)

censE <- 0
censD <- 0
G <- rep(0,4)
cross.rate0 <- 0
cross.rate1 <- 0

# ZZ <- matrix(rnorm(1000*length(pred.t)), nrow=1000, ncol=length(pred.t))
for (m in 1:nsim) {
  cat("sim=", m)
  print(Sys.time())

  dat0 <- Simulate(n, r, v, c, a, b, xi, lambda, eta, gamma) 
   
  ############ from now on, use only dat0 ################
  cross.rate0 <- cross.rate0 + with(subset(dat0, R==0), sum(is.finite(Vtime))/length(Vtime)) 
  cross.rate1 <- cross.rate1 + with(subset(dat0, R==1), sum(is.finite(Vtime))/length(Vtime)) 
  
  # TM method
  dat1 <- JBSTranTM(dat0)
  dat1$U <- with(dat1, ifelse(Group==1, 0, ifelse(Group %in% c(2,3), 1, NA)))  
  dat1$Ee <- with(dat1, ifelse(Group==1, 0, ifelse(Group %in% c(2,3), 1, 0.5))); 
  dat1$Vtime2 <- with(dat1, ifelse(is.finite(Vtime), Vtime, 0)) # used in the JBS modeling to allow for correct model specification
  epsilon1 <- min(with(dat1, ifelse(is.finite(Vtime) & Vtime>W, Vtime-W, Vtime)))
  epsilon2 <- min(with(dat1, ifelse(is.finite(Vtime) & Y>W, Y-W, Y)))
  epsilon <- min(with(dat1, min(Y, W, Vtime, na.rm=T)), epsilon1, epsilon2)*0.25
  datCount0 <- CountDat(dat1$Y, dat1$DeltaY, dat1$W, dat1$DeltaE, as.matrix(dat1), dat1$Vtime, epsilon)
  datCount0$RVt <- with(datCount0, R*Vt)
  Umod <- U ~ R + x
  TDmod <- Surv(StartD, EndD, DeltaD) ~ R + x + Vtime2 + Vt + RVt 
  TUmod <- Surv(W, DeltaE) ~ R + x
  TGmod <- Surv(StartG, EndG, DeltaG) ~ R + x + W + Vtime2 + Vt + RVt
  ret <- TM.fun.WCoxph(datCount0, dat1, Umod, TDmod, TUmod, TGmod, sd=T, ID='PTID')
  
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
    
    X <- as.matrix(cbind(intercept=rep(1,n), model.frame(Umod, data=dat1, na.action="na.pass")[,-1]))
    Xd <- as.matrix(model.frame(TDmod, data=datCount0, na.action="na.pass"))[datCount0$Last==1,-(1:3)]
    Xe <- model.matrix(TUmod, dat1)[,-1]
    Xg <- as.matrix(model.frame(TGmod, data=datCount0, na.action="na.pass"))[datCount0$Last==1,-(1:3)]
    Xbaseline <- matrix(X[, -c(1,2)], n, 1)
    Xd[,'Vtime2'] <- Xd[,'Vt'] <- Xd[,'RVt'] <- Xg[,'Vtime2'] <- Xg[,'Vt'] <- Xg[,'RVt'] <- 0  # prediction assuming no switching
    
    g0e <- Predict.Fast.Nosd(pred.t=pred.t, status=0, para=par, A0=dat1$R, Group=dat1$Group, X=X, Xd=Xd, Xe=Xe, 
                             Xg=Xg, Xbaseline=Xbaseline)$surv.pred 
    g1e <- Predict.Fast.Nosd(pred.t=pred.t, status=1, para=par, A0=dat1$R, Group=dat1$Group, X=X, Xd=Xd, Xe=Xe, 
                             Xg=Xg, Xbaseline=Xbaseline)$surv.pred 
    parest[m,] <- c(beta0est, estH0, beta1est, estH1, beta2est, estH2, alpest, g0e, g1e)
  }else parest[m,] <- NA
  
  # ITT with KM
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
  JBSmodInit <- Surv(latent.time, event.ind) ~ R + x + Vt + switch.ind
  JBSmod <- Surv(latent.time, event.ind) ~ R + x + Vt + switch.ind + Rswitch.ind
  JBSfinal <- Surv(latent.time, event.ind) ~ R + x
  pThreshM <- 0.98; pThreshI <- 0.98; MaxIter <- 200
  JBS <- LuoJBS(d=dat0, JBSmodInit, JBSmod, JBSfinal, pThreshM, pThreshI, MaxIter=MaxIter)
  JBS_result <- JBS$survreg.out
  JBS_converge[m] <- JBS$converge 
  if(JBS_converge[m]) {
    JBSest[m] <- JBS_result$coef['R'] # the true value is -log(0.75), fit exponential model
    JBS0[m,] <- pweibull(pred.t, scale=exp(JBS_result$coef[1]), shape=1/JBS_result$scale, lower.tail=FALSE)
    JBS1[m,] <- pweibull(pred.t*exp(-JBS_result$coef['R']/JBS_result$scale), scale=exp(JBS_result$coef[1]), shape=1/JBS_result$scale, lower.tail=FALSE)
  } else JBS0[m,] <- JBS1[m,] <- NA

  save(summary, file='sim_result_102016.RData')
  print(Sys.time())  
}

mean(exp(-JBSest), na.rm=T)
table(JBS_converge)

# save(list=ls(all=TRUE), file="working.RData")
# save(censE, censD, G, parest, file='sim_result_102016.RData')
Lp <- length(estH0) 
L0 <- length(beta0est)+length(beta1est)+length(beta2est)+3*length(estH0)+length(alpest)
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
#n=400
# ITTKM ITTCox  NoCrossKM        JBS         TM
# control   0.1104346     NA 0.10371643 0.02784449 0.03748133
# treatment 0.1121627     NA 0.06069371 0.03168602 0.03225910
#n=1000
# ITTKM ITTCox  NoCrossKM        JBS        TM
# control   0.1073274     NA 0.10117044 0.01739263 0.0310268
# treatment 0.1091909     NA 0.05521104 0.01947555 0.0215636

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
# n=400
# ITTKM ITTCox NoCrossKM        JBS         TM
# control   0.1540973     NA 0.1991384 0.03243242 0.07565301
# treatment 0.1790544     NA 0.1099988 0.03287447 0.05333033
# n=1000
# ITTKM ITTCox  NoCrossKM        JBS         TM
# control   0.1452306     NA 0.18430159 0.02132002 0.06310208
# treatment 0.1753483     NA 0.09463076 0.02082266 0.03598575


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
# save(surv.mat, file="SIMJBS_surv2mat.Rdata")
save(surv.mat, file="SIMJBS_surv2mat1000.Rdata")

#pdf("SimII_n400_col.pdf", height = 6, width = 9)
pdf("SimII_n1000_col.pdf", height = 6, width = 9)
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

# pdf("SimII_n400_bw.pdf", height = 6, width = 9)
pdf("SimII_n1000_bw.pdf", height = 6, width = 9)
par(mfrow=c(2,2), mar=c(4,4,2,2)) 
plot(x=surv.mat$time, y=surv.mat$True0, type='n', xlim=c(0,3), ylim=c(0,1), lty=1, col=1, xlab="Time", ylab="Survival Probability")
lines(x=surv.mat$time, y=surv.mat$True0, xlim=c(0,3), type="s", lty=1, lwd=1)
lines(x=surv.mat$time, y=surv.mat$True1, xlim=c(0,3), type="s",lty=2, lwd=1)
lines(x=surv.mat$time, y=surv.mat$ITTKM0, xlim=c(0,3), type="s", lty=3, lwd=0.5)
lines(x=surv.mat$time, y=surv.mat$ITTKM1, xlim=c(0,3), type="s",lty=4, lwd=0.5)
# legend("topright", c("True S(t|R=0)", "True S(t|R=1)", "ITT KM S(t|R=0)", "ITT KM S(t|R=1)"), bty='n', lty=rep(c(1,2),times=6), col=rep(c('grey', 'black'), each=2), cex=0.7)

plot(x=surv.mat$time, y=surv.mat$True0, type='n', xlim=c(0,3), ylim=c(0,1), lty=1, col=1, xlab="Time", ylab="Survival Probability")
lines(x=surv.mat$time, y=surv.mat$True0, xlim=c(0,3), type="s", lty=1, lwd=1)
lines(x=surv.mat$time, y=surv.mat$True1, xlim=c(0,3), type="s",lty=2, lwd=1)
lines(x=surv.mat$time, y=surv.mat$NoCrosKM0, xlim=c(0,3), type="s", lty=3, lwd=0.5)
lines(x=surv.mat$time, y=surv.mat$NoCrosKM1, xlim=c(0,3), type="s",lty=4, lwd=0.5)
# legend("topright", c("True S(t|R=0)", "True S(t|R=1)", "No Crossover KM S(t|R=0)", "No Crossover KM S(t|R=1)"), bty='n',lty=rep(c(1,2),times=6), col=rep(c('grey', 'black'), each=2), cex=0.7)

plot(x=surv.mat$time, y=surv.mat$True0, type='n', xlim=c(0,3), ylim=c(0,1), lty=1, col=1, xlab="Time", ylab="Survival Probability")
lines(x=surv.mat$time, y=surv.mat$True0, xlim=c(0,3), type="s", lty=1, lwd=1)
lines(x=surv.mat$time, y=surv.mat$True1, xlim=c(0,3), type="s",lty=2, lwd=1)
lines(x=surv.mat$time, y=surv.mat$JBS0, xlim=c(0,3), type="s", lty=3, lwd=0.5)
lines(x=surv.mat$time, y=surv.mat$JBS1, xlim=c(0,3), type="s",lty=4, lwd=0.5)
# legend("topright", c("True S(t|R=0)", "True S(t|R=1)", "JBS S(t|R=0)", "JBS S(t|R=1)"), bty='n', lty=rep(c(1,2),times=6), col=rep(c('grey', 'black'), each=2), cex=0.7)

plot(x=surv.mat$time, y=surv.mat$True0, type='n', xlim=c(0,3), ylim=c(0,1), lty=1, col=1, xlab="Time", ylab="Survival Probability")
lines(x=surv.mat$time, y=surv.mat$True0, xlim=c(0,3), type="s", lty=1, lwd=1)
lines(x=surv.mat$time, y=surv.mat$True1, xlim=c(0,3), type="s",lty=2, lwd=1)
lines(x=surv.mat$time, y=surv.mat$TM0, xlim=c(0,3), type="s", lty=1, lwd=0.5)
lines(x=surv.mat$time, y=surv.mat$TM1, xlim=c(0,3), type="s",lty=2, lwd=0.5)
# legend("topright", c("True S(t|R=0)", "True S(t|R=1)", "TM S(t|R=0)", "TM S(t|R=1)"), bty='n', lty=rep(c(1,2),times=6), col=rep(c('grey', 'black'), each=2), cex=0.7)
dev.off()

# generate the figure for LDA
library(ggplot2)
library(gridExtra)

load("SIMJBS_surv2mat.Rdata")
K <- length(surv.mat$time)
d2 <- data.frame(Time=rep(surv.mat$time, 16), Method=c(rep('a',4*K), rep('b',4*K), rep('c',4*K), rep('d',4*K)), 
                 Group=rep(rep(1:4, each=K), times=4), 
                 Surv=c(surv.mat$True0, surv.mat$True1, surv.mat$ITTKM0, surv.mat$ITTKM1, 
                        surv.mat$True0, surv.mat$True1, surv.mat$NoCrosKM0, surv.mat$NoCrosKM1,
                        surv.mat$True0, surv.mat$True1, surv.mat$JBS0, surv.mat$JBS1,
                        surv.mat$True0, surv.mat$True1, surv.mat$TM0, surv.mat$TM1))

p2 <- ggplot(d2,aes(x=Time)) + 
  geom_line(data=d2[d2$Group %in% c(1,3),], aes(y=Surv, linetype=as.factor(Group), colour=as.factor(Group)), size=0.5, alpha=0.9) +
  geom_line(data=d2[d2$Group %in% c(2,4),], aes(y=Surv, linetype=as.factor(Group), colour=as.factor(Group)), size=1, alpha=0.9) +
  scale_linetype_manual(values=c("solid", "dashed", "solid", "dashed")) +
  scale_color_manual(values=c("grey59", "grey59", "black",'black')) + 
  labs(x = "Time", y = "Simulation II", color=NULL) + 
  scale_x_continuous(limits = c(0, 3)) + 
  theme_bw() + 
  facet_grid(.~ Method) + 
  theme(legend.position="none")
p2

load("SIMI_surv2mat.Rdata")

K <- length(surv.mat$time)
d1 <- data.frame(Time=rep(surv.mat$time, 16), Method=c(rep('a',4*K), rep('b',4*K), rep('c',4*K), rep('d',4*K)), 
                 Group=rep(rep(1:4, each=K), times=4), 
                 Surv=c(surv.mat$True0, surv.mat$True1, surv.mat$ITTKM0, surv.mat$ITTKM1, 
                        surv.mat$True0, surv.mat$True1, surv.mat$NoCrosKM0, surv.mat$NoCrosKM1,
                        surv.mat$True0, surv.mat$True1, surv.mat$JBS0, surv.mat$JBS1,
                        surv.mat$True0, surv.mat$True1, surv.mat$TM0, surv.mat$TM1))

p1 <- ggplot(d1,aes(x=Time, y=Surv, group=Group)) + 
  geom_line(data=d1[d1$Group %in% c(1,3),], aes(y=Surv, linetype=as.factor(Group), colour=as.factor(Group)), size=0.5, alpha=0.9) +
  geom_line(data=d1[d1$Group %in% c(2,4),], aes(y=Surv, linetype=as.factor(Group), colour=as.factor(Group)), size=1, alpha=0.9) +
    scale_linetype_manual(values=c("solid", "dashed", "solid", "dashed")) +
    scale_color_manual(values=c("grey59", "grey59", "black",'black')) + 
  labs(x = "Time", y = "Simulation I", color=NULL) + 
  scale_x_continuous(limits = c(0, 3)) + 
  theme_bw() + 
  facet_grid(.~ Method) + 
  theme(legend.position="none")
p1

postscript("SimIandII_n400_bw.ps")
# pdf("SimIandII_n400_bw.pdf", height = 9, width = 9)
grid.arrange(p1, p2, nrow = 2)
dev.off()

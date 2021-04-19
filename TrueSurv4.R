Simulate <- function(n, beta0vec, beta0Tvec, beta1vec, beta2vec, beta2Tvec, alpvec, gam1vec, tau) { 

  # baseline covariates XbaseAtt tells whether/not discrete     
  R <-  c(rep(1,n/2), rep(0, n/2))
  # x1 <- 2*(runif(n)-0.5)
  x1 <- rnorm(n, 1, 1)
  x2 <- runif(n)<0.6

  # prognostic factors at the disease progression
  z1 <- runif(n)

  # generate susceptibility indicator
  X <- cbind(rep(1,n), R, x1, x2)
  E <- rbinom(n, 1, 1/(1+exp(-X%*%alpvec)))

  # generate time to disease progression
  # baseline hazard e^(t), cumulative baseline hazard e^(t)-1
  # Te <- log(-log(runif(n))*exp(-Xe%*%beta1vec)+1)
  # baseline hazard 2t, cumulative baseline hazard t^2
  # Te <- sqrt(-log(runif(n))*exp(-Xe%*%beta1vec))
  Xe <- cbind(R, x1, x2)
  # baseline hazard 1, cumulative baseline hazard t
  Te <- -log(runif(n))*exp(-Xe%*%beta1vec)

  # generate swithing indicator
  # Xv <- cbind(rep(1,n), R, x1, x2)
  # V <- rbinom(n, 1, 1/(1+exp(-Xv%*%gam0vec)))
  
  # generate switching time
  Xv <- cbind(R, x1, x2)
  # baseline hazard 1/3, cumulative baseline hazard t/3
  Vtime <- -3*log(runif(n))*exp(-Xv%*%gam1vec)
  Vtime <- ifelse(Vtime<Te & E==1, Te, Vtime)
  # Vtime <- ifelse(V==1, Vtime, Inf)

  # generate time from disease progression to death with baseline hazard h2(t)=1/2 
  Xg <- cbind(R, x1, x2, z1, Te)  
  Xgt <- cbind(rep(1,n), R)  
  t0 <- -2*log(runif(n))*exp(-Xg%*%beta2vec)
  TGap <- ifelse(t0<=Vtime-Te, t0, (t0+Te-Vtime)*exp(-Xgt%*%beta2Tvec)+Vtime-Te)
 
  # generate latent death time with baseline hazard h0(t)=1/2
  Xd <- cbind(R, x1, x2)
  Xdt <- cbind(rep(1,n), R)  
  t0 <- -2*log(runif(n))*exp(-Xd%*%beta0vec)
  Tdl <- ifelse(t0<=Vtime, t0, (t0-Vtime)*exp(-Xdt%*%beta0Tvec)+Vtime)
  
  # generate observed death time 
  Td <- ifelse(E==0, Tdl, Te+TGap)

  # generate censoring time
  Cens <- 1+6*runif(n)
  Cens[Cens>tau] <- tau
    
  # Observed data
  DeltaY <- ifelse(Td<=Cens,1,0)     # death indicator
  DeltaE <- ifelse(Te <= Cens, 1, 0) # progression indicator
  Y <- ifelse(DeltaY==1, Td, Cens)   # death
  W <- ifelse(DeltaE==1, Te, Cens)   # progression
  W[E==0] <- Y[E==0]
  DeltaE[E==0] <- 0
  z1[E==0]<- NA
  Vtime <- ifelse(Vtime>=Y, Inf, Vtime) 
  Vtime <- ifelse(Vtime<W & E==1, Inf, Vtime) 
  simdat <- data.frame(ID=1:n, R, x1, x2, z1, Vtime, DeltaE, W, DeltaY, Y) 

  # define 4 groups
  Group <- ifelse(DeltaE==0 & DeltaY==1, 1, 
             ifelse(DeltaE==1 & DeltaY==1, 2, 
                ifelse(DeltaE==1 & DeltaY==0, 3,
                  ifelse(DeltaE==0 & DeltaY==0, 4, NA))))
  simdat <- data.frame(simdat, Group)

  return(simdat)
  
}

################################################################################

TrueSurv <- function(pred.t, n=10000000, beta0vec, beta0Tvec, beta1vec, beta2vec, beta2Tvec, alpvec, gam1vec)
{
  # baseline covariates XbaseAtt tells whether/not discrete     
  R <-  c(rep(1,n/2), rep(0, n/2))
  # x1 <- 2*(runif(n)-0.5)
  x1 <- rnorm(n, 1, 1)
  x2 <- runif(n)<0.6
  
  # prognostic factors at the disease progression
  z1 <- runif(n)
  
  # generate susceptibility indicator
  X <- cbind(rep(1,n), R, x1, x2)
  E <- rbinom(n, 1, 1/(1+exp(-X%*%alpvec)))
  
  # generate time to disease progression
  # baseline hazard e^(t), cumulative baseline hazard e^(t)-1
  # Te <- log(-log(runif(n))*exp(-Xe%*%beta1vec)+1)
  # baseline hazard 2t, cumulative baseline hazard t^2
  # Te <- sqrt(-log(runif(n))*exp(-Xe%*%beta1vec))
  Xe <- cbind(R, x1, x2)
  # baseline hazard 1, cumulative baseline hazard t
  Te <- -log(runif(n))*exp(-Xe%*%beta1vec)
  
  # generate swithing indicator
  # Xv <- cbind(rep(1,n), R, x1, x2)
  # V <- rbinom(n, 1, 1/(1+exp(-Xv%*%gam0vec)))
  
  # generate switching time
  Xv <- cbind(R, x1, x2)
  # baseline hazard 1/3, cumulative baseline hazard t/3
  Vtime <- -3*log(runif(n))*exp(-Xv%*%gam1vec)
  Vtime <- ifelse(Vtime<Te & E==1, Te, Vtime)
  # Vtime <- ifelse(V==1, Vtime, Inf)
  
  # generate time from disease progression to death with baseline hazard h2(t)=1/2 
  Xg <- cbind(R, x1, x2, z1, Te)  
  Xgt <- cbind(rep(1,n), R)  
  t0 <- -2*log(runif(n))*exp(-Xg%*%beta2vec)
  TGap <- ifelse(t0<=Vtime-Te, t0, (t0+Te-Vtime)*exp(-Xgt%*%beta2Tvec)+Vtime-Te)
  
  # generate latent death time with baseline hazard h0(t)=1/2
  Xd <- cbind(R, x1, x2)
  Xdt <- cbind(rep(1,n), R)  
  t0 <- -2*log(runif(n))*exp(-Xd%*%beta0vec)
  Tdl <- ifelse(t0<=Vtime, t0, (t0-Vtime)*exp(-Xdt%*%beta0Tvec)+Vtime)
  
  # generate observed death time 
  Td <- ifelse(E==0, Tdl, Te+TGap)
  
  # true HR 
  tau <- max(pred.t)
  Tdobs <- ifelse(Td>tau, tau, Td)
  delta <- ifelse(Td>tau, 0, 1)
  print(coxph(Surv(Tdobs, delta) ~ R))
  
  surv0 <- matrix(0, length(pred.t),1)
  surv1 <- matrix(0, length(pred.t),1)
  for (k in 1:length(pred.t)){
    surv0[k] <- mean(Td[R==0]>pred.t[k])
    surv1[k] <- mean(Td[R==1]>pred.t[k])
  }

  return(list(surv0, surv1))
}

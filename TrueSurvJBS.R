Simulate <- function(n, r, v, c, a, b, xi, lambda, eta, gamma) 
{ 
  V <- rbinom(n, 1, v)
  Z <- rbinom(n, 1, r)
  Ts <- rexp(n, rate=exp(-xi*V-lambda*Z))
  Cs <- runif(n, min=a, max=b)
  Tau0s <- rexp(n, rate=exp(-eta[1]*(1-Z)-eta[2]*Z))  
  Taus <- Tau0s^c*Ts^(1-c)
  TTau <- rexp(n, rate=exp(-xi*V-lambda*Z-gamma[1]*Taus-gamma[2]-gamma[3]*Z))
  Tau <- ifelse(Taus<pmin(Ts, Cs), Taus, 0)
  W <- ifelse(Tau>0, 1, 0)
  Ss <- ifelse(W==0, Ts, TTau)
  Css <- ifelse(W==0, Cs, Cs-Tau)
  Ytilde <- pmin(Ss, Css)
  Y <- Tau + Ytilde
  DeltaY <- ifelse(Ss <= Css, 1, 0)

  # change switching time to be Inf for no switching  
  Tau[Tau==0] <- Inf
  simdat <- data.frame(x=V, R=Z, Vtime=Tau, Y, DeltaY)
  return(simdat)
  
}

################################################################################
TrueSurv <- function(pred.t, n=10000000, r, v, xi, lambda)
{
  V <- rbinom(n, 1, v)
  Z <- rbinom(n, 1, r)
  Ts <- rexp(n, rate=exp(-xi*V-lambda*Z))
  surv0 <- matrix(0, length(pred.t),1)
  surv1 <- matrix(0, length(pred.t),1)
  for (k in 1:length(pred.t)){
    surv0[k] <- mean(Ts[Z==0]>pred.t[k])
    surv1[k] <- mean(Ts[Z==1]>pred.t[k])
  }
  
  return(list(surv0, surv1))
}

################################################################################
JBSTranTMtmp <- function(d) 
{  
  W <- d$Vtime
  R <- d$R
  Y <- d$Y
  DeltaY <- d$DeltaY
  x <- d$x
  PTID <- 1:dim(d)[1]
  DeltaE <- ifelse(W < Y, 1, 0)
  W[DeltaE==0] <- Y[DeltaE==0]
  Group <- ifelse(DeltaE==0 & DeltaY==1, 1, 
                  ifelse(DeltaE==1 & DeltaY==1, 2, 
                         ifelse(DeltaE==1 & DeltaY==0, 3,
                                ifelse(DeltaE==0 & DeltaY==0, 4, NA))))
  dat <- data.frame(PTID, R, x, Y, DeltaY, W, DeltaE, Group) 
  return(dat)
}

################################################################################
JBSTranTM <- function(d) 
{  
  R <- d$R
  Y <- d$Y
  DeltaY <- d$DeltaY
  x <- d$x
  n <- dim(d)[1]
  PTID <- 1:n
  Te <- runif(n, 0, max(Y)*1.1)
  Vtime <- d$Vtime

  DeltaE <- ifelse(Te < Y, 1, 0)
  W <- ifelse(DeltaE==1, Te, Y)  
  Group <- ifelse(DeltaE==0 & DeltaY==1, 1, 
                  ifelse(DeltaE==1 & DeltaY==1, 2, 
                         ifelse(DeltaE==1 & DeltaY==0, 3,
                                ifelse(DeltaE==0 & DeltaY==0, 4, NA))))
  dat <- data.frame(PTID, R, x, Y, DeltaY, W, DeltaE, Group, Vtime) 
  return(dat)
}
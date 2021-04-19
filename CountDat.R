 
CountDat <- function(Survt, Event, Progt, ProgEvt, X, Ts, epsilon) {
# construct counting process data for data

  N  <- length(Survt)
  n <- sum(is.finite(Ts))
  t1 <- rep(0, N+n)  # initialize start time at 0 for Td model
  t2 <- rep(NA, N+n)    # build vector for end times for Td model
  t1g <- rep(0, N+n)  # initialize start time at 0 for Tg model
  t2g <- rep(NA, N+n)    # build vector for end times for Tg model
  last <- rep(NA, N+n)    
  d  <- rep(NA, N+n)    # event indicator
  dg  <- rep(NA, N+n)    # event indicator for Tg model
  Xv <- matrix(NA, N+n, ncol(X)) # time independent covariates matrix
  Vt <- rep(0, N+n) # time-dependent switching indicator

  j  <- 1
  for(ii in 1:N){
    if(!is.finite(Ts[ii])){      # no Treatment switching
      t1[j] <- 0
      t2[j] <- Survt[ii]
      d[j]  <- Event[ii]
      t1g[j] <- 0
      if(Survt[ii]==Progt[ii]) t2g[j] <- epsilon else t2g[j] <- Survt[ii]-Progt[ii]
      dg[j]  <- Event[ii]*ProgEvt[ii]
      Xv[j,] <- as.vector(X[ii,])
      Vt[j] <- 0
      last[j] <- 1
      j <- j+1
    } else {                 # switching, split records
      ## create counting process for Td (death time)
      Xv[j, ] <- Xv[j+1,] <- X[ii,]  # time-independent covariates are same for each time
      t2[j]    <- Ts[ii]-epsilon # end of the first time interval that Vt=0
      d[j]     <- 0          # no event observed during (0,Ts-]
      Vt[j]    <- 0          # Vt=0 during (0,Ts-]
      last[j] <- 2
      t1[j+1]  <- t2[j]      # start of time interval that Vt=1
      t2[j+1]  <- Survt[ii]  # end of time interval that Vt=1
      d[j+1]   <- Event[ii]  # event occur after Ts?
      Vt[j+1]  <- 1          # Vt=1 during (Ts-, t2]
      last[j+1] <- 1
      
      ## create counting process for Tg (gap time)
      if(Ts[ii]<=Progt[ii]) t2g[j] <- epsilon else t2g[j] <- Ts[ii]-Progt[ii]-epsilon # end of the first time interval that Vt=0
      dg[j]     <- 0          # no event observed during (0,(Ts-Tu)-]
      t1g[j+1]  <- t2g[j]     # start of time interval that Vt=1
      if(Survt[ii]==Progt[ii]) t2g[j+1]  <- Survt[ii]-Progt[ii]+epsilon*2 else  t2g[j+1]  <- Survt[ii]-Progt[ii] # end of time interval that Vt=1
      dg[j+1]   <- Event[ii]*ProgEvt[ii]  # event occur after Ts?
      j <- j+2               # two records added
    }
  }

  dat.count <- as.data.frame(cbind(t1, t2, d, t1g, t2g, dg, Vt, Xv, last))
  names(dat.count) <- c("StartD", "EndD", "DeltaD", "StartG", "EndG", "DeltaG", "Vt", colnames(X), 'Last')
  return(dat.count)
}

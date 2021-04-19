Predict.Fast.Nosd <- function(pred.t, status, para, A0, Group, X, Xd, Xe, Xg, Xbaseline, z) 
{
    n <- length(A0)
    na <- sum(A0==status)
    
    beta0<-para$beta0vec; h0<-para$h0; y0<-para$timeD; beta1<-para$beta1vec; h1<-para$h1; y1<-para$timeU;
    beta2<-para$beta2vec; h2<-para$h2; y2<-para$timeG; alpha<-para$alpvec
    
    nbeta0 <- length(beta0); n0 <- length(h0)
    nbeta1 <- length(beta1); n1 <- length(h1)
    nbeta2 <- length(beta2); n2 <- length(h2)
    nalpha <- length(alpha); nt<- length(pred.t)

    H0t <- as.vector((matrix(pred.t, nt, n0, byrow=FALSE) >= matrix(y0, nt, n0, byrow=TRUE))%*%h0)
    H1 <- as.vector((matrix(y1, n1, n1, byrow=FALSE) >= matrix(y1, n1, n1, byrow=TRUE))%*%h1)
    H1t <- cbind((matrix(pred.t, nt, n1, byrow=FALSE) >= matrix(y1, nt, n1, byrow=TRUE))%*%h1)
    PE0 <- 1/(1+exp(X%*%alpha))   
   
    surv.pred <- matrix(0, nt, 1)
  
    # first term
    
    TEXbeta0 <- as.numeric(exp(Xd[A0==status,]%*%beta0))
    TPE0 <- PE0[A0==status]
    temp <- cbind(H0t)%*%rbind(TEXbeta0)
    SDt <- exp(-temp)
    surv.pred <- surv.pred+SDt%*%TPE0/na
    term1 <- apply(SDt, 1, mean)
    
    # second term

   weight <- matrix(0,n,n)
   an <- n^(-1/6)
   for (i in 1:n){
      if (A0[i]==status){
      for (j in 1:n){
        if ((A0[j]==status)& (Group[j]%in%c(2,3))){
            weight[i,j] <- exp(-sum((Xbaseline[j,]-Xbaseline[i,])^2)/an)
        }
      }
     }
   }
   den <- weight%*%rep(1,n)
   den[den==0] <- 1
   StWeight <- t(weight)%*%((1-PE0)/den*(A0==status))/na
   
   Tindex <- ((Group%in%c(2,3)) & A0==status)
   Tn <- sum(Tindex)
   TXe <- as.matrix(Xe[Tindex,])
   TXg <- as.matrix(Xg[Tindex,])
   TEXbeta1 <- as.numeric(exp(TXe%*%beta1))
   TEXbeta2 <- as.numeric(exp(TXg%*%beta2))

   St <- matrix(0, nt, n) # original St <- matrix(0, nt, n) 
   for (k in 1:length(pred.t)) {
    tk <- pred.t[k]-y1
    ntk <- length(tk)
    indtky2 <- (matrix(tk, ntk, n2, byrow=FALSE) >= matrix(y2, ntk, n2, byrow=TRUE))
    H2tk <- as.vector(indtky2%*%h2)
    temp1 <- cbind(H1)%*%rbind(TEXbeta1)
    temp2 <- cbind(H2tk)%*%rbind(TEXbeta2)
    tempk <- (exp(-temp2)*exp(-temp1)*(cbind(h1)%*%rbind(TEXbeta1)))
    Zero <- rep(0, ntk)
    Zero[tk>=0] <- 1
    
    St[k,Tindex] <- as.vector(t(tempk)%*%Zero)+exp(-H1t[k]*rbind(TEXbeta1))
  }
   
   surv.pred <- surv.pred+St%*%StWeight
   term2 <- apply(St, 1, sum)/sum(Tindex)
   
   return(list(surv.pred=surv.pred, term1=term1, term2=term2))
}


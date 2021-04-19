Covest <- function(par, Umod, TDmod, TUmod, TGmod, SdatCount, Sdat, ID) {

    # y0/y1/y2 are the jump time of h0/h1/h2 exclude censored
    beta0 <- par$beta0vec; h0 <- par$h0; y0 <- par$timeD; nbeta0 <- par$nbeta0; n0 <- par$n0
    beta1 <- par$beta1vec; h1 <- par$h1; y1 <- par$timeU; nbeta1 <- par$nbeta1; n1 <- par$n1
    beta2 <- par$beta2vec; h2 <- par$h2; y2 <- par$timeG; nbeta2 <- par$nbeta2; n2 <- par$n2 
    alpha <- par$alpvec; nalpha <- par$nalpha 

    IND <- cumsum(c(nbeta0, n0, nbeta1, n1, nbeta2, n2, nalpha))
    indbeta0 <- 1:IND[1];          indh0 <- (IND[1]+1):IND[2]
    indbeta1 <- (IND[2]+1):IND[3]; indh1 <- (IND[3]+1):IND[4]
    indbeta2 <- (IND[4]+1):IND[5]; indh2 <-(IND[5]+1):IND[6]
    indalpha <- (IND[6]+1):IND[7]
    npar <- max(indalpha)

    n <- dim(SdatCount)[1]
    nid <- dim(Sdat)[1]; 
    d2L <- matrix(0, npar, npar)
    score1 <- matrix(0, nid, npar)
    score0 <- matrix(0, nid, npar)
    X <- as.matrix(cbind(intercept=rep(1,nid), model.frame(Umod, data=Sdat, na.action="na.pass")[,-1]))
    Xd <- model.matrix(TDmod, SdatCount)[,-1]
    Xe <- model.matrix(TUmod, Sdat)[,-1]
    temp <- as.matrix(model.frame(TDmod, data=SdatCount, na.action="na.pass"))
    StartD <- temp[,1]; EndD <- temp[,2]; DeltaD <- temp[,3]
    temp <- as.matrix(model.frame(TUmod, data=Sdat, na.action="na.pass"))
    SurvU <- temp[,1]; DeltaU <- temp[,2]
    temp <- as.matrix(model.frame(TGmod, data=SdatCount, na.action="na.pass"))
    Xg <- temp[,-(1:3)]
    StartG <- temp[,1]; EndG <- temp[,2]; DeltaG <- temp[,3]
    
    indID <- (matrix(SdatCount[,ID], n, nid, byrow=F)==matrix(Sdat[,ID], n, nid, byrow=T))
    
     # estimating alpha

     PE1 <- exp(X%*%alpha)/(1+exp(X%*%alpha))  
     score1[, indalpha] <- matrix(1-PE1, nid, nalpha, byrow=FALSE)*(X)
     score0[, indalpha] <- matrix(0-PE1, nid, nalpha, byrow=FALSE)*(X)
     d2L[indalpha, indalpha] <- -(t(X)*matrix(PE1*(1-PE1),nalpha,nid,byrow=TRUE))%*%X   
   
     #  beta0 h0
     out <- Covestsub(Xd, beta0, h0, y0, StartD, EndD, DeltaD, (SdatCount$Group%in%c(1,4)), 1-SdatCount$Ee)
     score0[, c(indbeta0, indh0)] <- t(indID)%*%out$score
     d2L[c(indbeta0, indh0), c(indbeta0, indh0)] <- out$dscore
     
     #  beta1 h1
     out <- Covestsub(Xe, beta1, h1, y1, rep(0,nid), SurvU, DeltaU, (Sdat$Group!=1), Sdat$Ee)
     score1[, c(indbeta1, indh1)] <- out$score
     d2L[c(indbeta1, indh1), c(indbeta1, indh1)] <- out$dscore

    #  beta2 h2
     out <- Covestsub(Xg, beta2, h2, y2, StartG, EndG, DeltaG, (SdatCount$Group%in% c(2,3)), SdatCount$Ee)
     score1[, c(indbeta2, indh2)] <- t(indID)%*%out$score
     d2L[c(indbeta2, indh2), c(indbeta2, indh2)] <- out$dscore

    # information matrix
     Ee <- Sdat$Ee
     EdLdL <- t(score1)%*%(as.vector(Ee)*score1)+t(score0)%*%(as.vector(1-Ee)*score0)
     EdL <- as.vector(Ee)*score1+as.vector(1-Ee)*score0
     EdLEdL <- t(EdL)%*%EdL

     cov.est <- solve(-d2L-(EdLdL-EdLEdL))
     Influence <- cov.est%*%t(EdL)
     return(list(cov.est=cov.est, Infl=Influence))
}


Covestsub <- function(X, beta, hh0, yy0, Start, End, Delta, w, Ew) {
     XX<- X[w>0,]
     TStart<- Start[w>0]
     TEnd<- End[w>0]
     TDelta <- Delta[w>0]
     Tw <- w[w>0]
     ETw <- Ew[w>0]

     n <- length(TStart)
     nbeta <- length(beta)
     nh <- length(hh0)
     score <- matrix(0, n, nbeta+nh)
     dscore <- matrix(0, nbeta+nh, nbeta+nh)

     EXbeta <- exp(XX%*%beta)
     indY <- (matrix(TEnd, n, nh, byrow=FALSE)>=matrix(yy0,n,nh,byrow=TRUE))*(matrix(TStart, n, nh, byrow=FALSE)<matrix(yy0,n,nh,byrow=TRUE))
     indY2 <- (matrix(TEnd, n, nh, byrow=FALSE)<matrix(c(yy0[-1],Inf),n,nh,byrow=TRUE))
     HY <- indY%*%hh0
     OneY <- (indY*indY2)%*%rep(1,nh)
     score[,1:nbeta] <- as.vector(TDelta*OneY)*XX-as.vector(Tw*HY*EXbeta)*XX
     score[,nbeta+(1:nh)]    <- as.vector(TDelta)*indY*indY2%*%diag(1/hh0)-as.vector(Tw*EXbeta)*indY

     dscore[1:nbeta, 1:nbeta] <- -t(XX)%*%(matrix(ETw*HY*EXbeta,n,nbeta, byrow=FALSE)*XX)
     dscore[1:nbeta, nbeta+(1:nh)] <- -t(XX)%*%(matrix(ETw*EXbeta,n,nh,byrow=FALSE)*indY)
     dscore[nbeta+(1:nh), 1:nbeta] <- t(dscore[1:nbeta, nbeta+(1:nh)])
     dscore[nbeta+(1:nh), nbeta+(1:nh)] <- -diag(apply(TDelta*indY*indY2,2,sum)/hh0^2)
     
     Tscore <- score
     score <- matrix(0, length(Start), nbeta+nh)
     score[w>0,] <- Tscore
      
     out=list(score=score, dscore=dscore)
     return(out)
}

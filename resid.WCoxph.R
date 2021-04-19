resid.WCoxph <- function(SdatCount, Sdat, Umod, TDmod, TUmod, TGmod, par, type) 
{
  ret <- Estep3(oldpara=par, SSdatCount=SdatCount, SSdat=Sdat, Umod, TDmod, TUmod)
  Sdat$Ee <- ret$Eed; SdatCount$Ee <- ret$EedC
  
  datG4.1 <- datG4.2 <- subset(Sdat, Group==4)
  datG123 <- subset(Sdat, Group %in% c(1,2,3))
  datG4.1$U <- 0; datG4.1$uweight <- 1-datG4.1$Ee
  datG4.2$U <- 1; datG4.2$uweight <- datG4.2$Ee
  datG123$uweight <- 1
  datGlm <- as.data.frame(rbind(datG123, datG4.1, datG4.2))
  datGlm$U <- as.integer(datGlm$U)
  
  # model for beta0
  fit0 <- coxph(as.formula("TDmod"), data=subset(SdatCount, Ee<1), weight=1-Ee, init=par$beta0vec)
  res0 <- resid(fit0, type=type)
  
  # model for beta1
  fit1 <- coxph(as.formula("TUmod"), data=subset(Sdat, Ee>0), weight=Ee, init=par$beta1vec)
  res1 <- resid(fit1, type=type)
  
  # model for beta2
  fit2 <- coxph(as.formula("TGmod"), data=subset(SdatCount, Group%in%c(2,3) & StartG<EndG), init=par$beta2vec)
  res2 <- resid(fit2, type=type)
  
  return(list(TDres=res0, TUres=res1, TGres=res2, SdatCount=SdatCount, Sdat=Sdat))
}  

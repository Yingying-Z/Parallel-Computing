#########################################################################
# DESCRIPTION: functions that will be called in CYD15VE_hosp_yr3.R
#
# CODED BY: Yingying Zhuang
# SAVED AS: CYD15VEfunctions
#
# HISTORY:
# Date      Programmer   	
# Nov22-2015  Y Zhuang      Version 1: exactly same as CYD14VEfunctions
# Dec18-2015  Y Zhuang      Version 2: added perturbation functions, deleted non-hinge function, cleaned up stuff
##########################################################################

expit<-function(x) 1/(1+exp(-x))

#For categorical variables, get dummy variable for AGE and COUNTRY
getDummy.New<-function(dat){
  
  dat$AGE2=floor(dat$W/100)==2
  #Country
  dat$C2=dat$W-floor(dat$W/100)*100==2
  dat$C3=dat$W-floor(dat$W/100)*100==3
  dat$C4=dat$W-floor(dat$W/100)*100==4
  dat$C5=dat$W-floor(dat$W/100)*100==5
  return(dat)
}

R.logit<-function(Z,S1,W,beta,short){
  ####### calculate the logit risk ############
  
  Risk<-NULL
  datin=data.frame(Z=Z,S1=S1,W=W)
  if (short=="noW"){
  X<-cbind(1,datin$Z,datin$S1,datin$S1*datin$Z)
  }
  if (short=="withW"){
  datin=getDummy.New(datin)
  X<-cbind(1,datin$Z,datin$S1,datin$S1*datin$Z,datin$AGE2,datin$AGE2*datin$Z,datin$C2)
  }
  if (short=="withW2"){
    datin=getDummy.New(datin)
    X<-cbind(1,datin$Z,datin$S1,datin$S1*datin$Z,datin$AGE2,datin$AGE2*datin$Z)
  }
  Risk<-expit(X%*%as.matrix(beta))
  return(Risk)
}

try.error<-function(expr, silent=TRUE){
  
  # make it to error
  op <- options("warn")
  on.exit(options(op))
  options(warn = -1)
  
  # catch using the try() function
  
  try(expr, silent=TRUE)
}

###############################################################################################################
getVE<-function(dat,maxit,epsilon,short,Su){
  #Calculate the PSN beta's and VE after integrate on W using multinotmial splines
  Z<-dat$Z
  S1<-dat$S1
  W<-dat$W
  Y<-dat$Y
  delta<-dat$delta.CPV
  Wu<-unique(sort(dat$W)) 

  
  ### Recall from the paper: For a subject with delta_i=0, W_i=w_i, construct a set of filled-in data with lenfth equal to the number
  ## of observations in V, wjere V is the set of validation subjects with delta=1 and W=w_i
  
  ##constructing S in the validation dataset V for each unique value of W
  SS<-ll<-Sweights<-list()
  for (i in 1:length(Wu)){
    SS[[i]]<-S1[delta==1&W==Wu[i]]
    ll[[i]]<-sum(delta==1&W==Wu[i])
  }
  
  Sval<-SS
  lval<-ll
  
  rm(SS,ll,Sweights,i)
  
  oo<-match(W[delta==0],Wu)
  ind<-rep(1:sum(delta==0), unlist(lval[oo])) # rep(x, times) when length(times)=length(x), it gives the number of times to repeat each element 
  YY<-rep(Y[delta==0],unlist(lval[oo]))
  WW<-rep(W[delta==0],unlist(lval[oo]))
  ZZ<-rep(Z[delta==0],unlist(lval[oo]))
  SS1<-unlist(Sval[oo])
  
  if (short=="noW"){beta.old<-rep(0,4)}
  if(short=="withW"){beta.old<-rep(0,7)}
  if(short=="withW2"){beta.old<-rep(0,6)}
  
  iter<-0
  repeat{
    
    ### probability of delta=1 conditional on S1,W
    ### if random sample S1 from cases and controls in vaccine arm
    ### Also Recall in teh paper " the sampling prob depends on Y and Z only, we apply a saturated model with 
    ### pi={pi(0,0), pi(0,1), pi(1,0),pi(1,1)}
    P.Z1<-mean(delta[Z==1 & Y==1])*R.logit(Z=1,S1=SS1,W=WW,beta=beta.old,short)+
      mean(delta[Z==1 & Y==0])*(1-R.logit(Z=1,S1=SS1,W=WW,beta=beta.old,short))
    P.Z0<-mean(delta[Z==0 & Y==0])*(1-R.logit(Z=0,S1=SS1,W=WW,beta=beta.old,short)) #because CPV component are from controls in Placebo
    PP<-P.Z1*mean(Z)+P.Z0*(1-mean(Z))
    
    
    ff<-R.logit(Z=ZZ,S1=SS1,W=WW,beta=beta.old,short)
    ff<-ifelse(YY==1,ff,1-ff)
    ff<-ff/PP
    
    
    if (length(ff)>0){
      denom<-tapply(ff,ind,sum)
      denom<-denom[ind]
      weights<-as.numeric(ff)/denom
      weights<-as.numeric(weights)
    } else weights<-ff
    
    
    dataout<-data.frame(Y=c(Y[delta==1],YY),Z=c(Z[delta==1],ZZ),W=c(W[delta==1],WW),
                        S1=c(S1[delta==1],SS1),weights=c(rep(1,sum(delta==1)),weights))
    
    
    dataout=getDummy.New(dataout)
    
    if (short=="noW"){
      fit<-glm(Y~Z+S1+I(Z*S1),weights=weights,data=dataout,family=binomial(link=logit),start=beta.old)
      }
    if(short=="withW"){
      fit<-glm(Y~Z+S1+I(Z*S1)+AGE2+I(AGE2*Z)+C2,weights=weights,data=dataout,family=binomial(link=logit),start=beta.old)
    }
    if(short=="withW2"){
      fit<-glm(Y~Z+S1+I(Z*S1)+AGE2+I(AGE2*Z),weights=weights,data=dataout,family=binomial(link=logit),start=beta.old)
    }

    out<-fit$coef
    
    beta.new<-out
    
    iter<-iter+1
    
    diffcount<-sum(abs(beta.new-beta.old))<epsilon
    
    print(c(iter,beta.new));
    if (is.na(diffcount)) break
    else{
      if (sum(abs(beta.new-beta.old))<epsilon & sum(abs(beta.new-beta.old))<sum(abs(beta.old))*epsilon) break
    }
    
    beta.old<-beta.new
    
    
    if (iter>maxit | is.na(diffcount)) break
  }
  loglike<-logLik(fit)
  if (iter>maxit | is.na(diffcount)) {
    beta.new<-rep(NA,length(beta.old))
    loglike<-NA
  }
  
  beta<-beta.new
  
  fX<-table(dat$W)/sum(!is.na(dat$W))
  
  datnew<-data.frame(S1=Su)
  
  dat.Z1=dat[dat$Z==1 & !is.na(dat$S1),]
  #weights of sampling S1 in VACC arm
  dat.Z1$weights<-rep(0,nrow(dat.Z1))
  dat.Z1$weights<-ifelse(dat.Z1$Y==1, sum(!is.na(dat$S1) & dat$Y==1 & dat$Z==1)/sum(dat$Y==1 & dat$Z==1),
                         sum(!is.na(dat$S1) & dat$Y==0 & dat$Z==1)/sum(dat$Y==0 & dat$Z==1))
  
  fit=multinom(W~ns(S1,knots=quantile(Su,c(.25,.5,.75)),Boundary.knots=range(Su)),data=dat.Z1,weights=1/dat.Z1$weights)
  
  fx.s=predict(fit,datnew,type='probs')
  PY1.SZ0<-PY1.SZ1<-rep(0,length(Su))
  for (i in 1:length(Wu)){
    PY1.SZ1=PY1.SZ1+R.logit(Z=1,S1=Su,W=Wu[i],beta=beta,short)*fx.s[,i]
    PY1.SZ0=PY1.SZ0+R.logit(Z=0,S1=Su,W=Wu[i],beta=beta,short)*fx.s[,i]
  }
  
  VE.S<-1-PY1.SZ1/PY1.SZ0
  
  return(list(beta=beta,VE=VE.S))
  
}

###############################################################################################################
getVE.Perturb<-function(dat,maxit,epsilon,short,xi){
  
  ### Recall from the paper: For a subject with delta_i=0, W_i=w_i, construct a set of filled-in data with lenfth equal to the number
  ## of observations in V, wjere V is the set of validation subjects with delta=1 and W=w_i
  
  ##constructing S in the validation dataset V for each unique value of W
  Z<-dat$Z
  S1<-dat$S1
  W<-dat$W
  Y<-dat$Y
  delta<-dat$delta.CPV
  Wu<-unique(sort(dat$W))  
  Su=sort(unique(dat$S1[!is.na(dat$S1)]))
  
  SS<-ll<-Sweights<-list()
  for (i in 1:length(Wu)){
    SS[[i]]<-S1[delta==1&W==Wu[i]]
    Sweights[[i]]<-xi[delta==1&W==Wu[i]] #purturbation weights for SS1
    ll[[i]]<-sum(delta==1&W==Wu[i])
    
  }
  
  Sval<-SS
  lval<-ll
  
  rm(SS,ll,i)
  
  oo<-match(W[delta==0],Wu)
  ind<-rep(1:sum(delta==0), unlist(lval[oo])) # rep(x, times) when length(times)=length(x), it gives the number of times to repeat each element 
  YY<-rep(Y[delta==0],unlist(lval[oo]))
  WW<-rep(W[delta==0],unlist(lval[oo]))
  ZZ<-rep(Z[delta==0],unlist(lval[oo]))
  SS1<-unlist(Sval[oo])
  SSweights<-unlist(Sweights[oo])  #purturbation weights for SS1
  
  if (short=="noW"){beta.old<-rep(0,4)}
  if(short=="withW"){beta.old<-rep(0,7)}
  if(short=="withW2"){beta.old<-rep(0,6)}
  
  
  ii<-0
  repeat{
    
    ### probability of delta=1 conditional on S1,W
    ### if random sample S1 from cases and controls in vaccine arm
    ### Also Recall in teh paper " the sampling prob depends on Y and Z only, we apply a saturated model with 
    ### pi={pi(0,0), pi(0,1), pi(1,0),pi(1,1)}
    P.Z1<-mean(delta[Z==1 & Y==1])*R.logit(Z=1,S1=SS1,W=WW,beta=beta.old,short)+
      mean(delta[Z==1 & Y==0])*(1-R.logit(Z=1,S1=SS1,W=WW,beta=beta.old,short))
    P.Z0<-mean(delta[Z==0 & Y==0])*(1-R.logit(Z=0,S1=SS1,W=WW,beta=beta.old,short)) #because CPV component are from controls in Placebo
    PP<-P.Z1*mean(Z)+P.Z0*(1-mean(Z))
    PP<-PP[,1]
    
    
    ff<-R.logit(Z=ZZ,S1=SS1,W=WW,beta=beta.old,short)
    ff<-ifelse(YY==1,ff,1-ff)
    ff<-ff/PP
    ff<-ff*SSweights #purturbation weights for SS1
    
    if (length(ff)>0){
      denom<-tapply(ff,ind,sum)
      denom<-denom[ind]
      weights<-as.numeric(ff)/denom
      weights<-as.numeric(weights)
    } else weights<-ff
    
    #add purturbation weights for Likelihood
    xxi<-rep(xi[delta==0],unlist(lval[oo]))
    weights<-weights*xxi
    dataout<-data.frame(Y=c(Y[delta==1],YY),Z=c(Z[delta==1],ZZ),W=c(W[delta==1],WW),
                        S1=c(S1[delta==1],SS1),weights=c(xi[delta==1],weights))
    
    
    dataout=getDummy.New(dataout)
    
    if (short=="noW"){
      fit<-glm(Y~Z+S1+I(Z*S1),weights=weights,data=dataout,family=binomial(link=logit),start=beta.old)
    }
    if(short=="withW"){
      fit<-glm(Y~Z+S1+I(Z*S1)+AGE2+I(AGE2*Z)+C2,weights=weights,data=dataout,family=binomial(link=logit),start=beta.old)
    }
    if(short=="withW2"){
      fit<-glm(Y~Z+S1+I(Z*S1)+AGE2+I(AGE2*Z),weights=weights,data=dataout,family=binomial(link=logit),start=beta.old)
    }

    out<-fit$coef
    
    beta.new<-out
    
    ii<-ii+1
    
    diffcount<-sum(abs(beta.new-beta.old))<epsilon
    
    print(c(ii,beta.new));
    if (is.na(diffcount)) break
    else{
      if (sum(abs(beta.new-beta.old))<epsilon & sum(abs(beta.new-beta.old))<sum(abs(beta.old))*epsilon) break
    }
    
    beta.old<-beta.new
    
    
    if (ii>maxit | is.na(diffcount)) break
  }
  loglike<-logLik(fit)
  if (ii>maxit | is.na(diffcount)) {
    beta.new<-rep(NA,length(beta.old))
    loglike<-NA
  }
  
  beta<-beta.new
  
  fX<-table(dat$W)/sum(!is.na(dat$W))
  
  datnew<-data.frame(S1=Su)
  
  dat.Z1=dat[dat$Z==1 & !is.na(dat$S1),]
  #weights of sampling S1 in VACC arm
  dat.Z1$weights<-rep(0,nrow(dat.Z1))
  dat.Z1$weights<-ifelse(dat.Z1$Y==1, sum(!is.na(dat$S1) & dat$Y==1 & dat$Z==1)/sum(dat$Y==1 & dat$Z==1),
                         sum(!is.na(dat$S1) & dat$Y==0 & dat$Z==1)/sum(dat$Y==0 & dat$Z==1))
  wts<-1/dat.Z1$weights
  wts2<-xi[dat$Z==1 & !is.na(dat$S1)] #purturbation for estimating multinomial parameters
  wts<-wts*wts2
  
  fit=multinom(W~ns(S1,knots=quantile(Su,c(.25,.5,.75)),Boundary.knots=range(Su)),data=dat.Z1,weights=wts)
  
  fx.s=predict(fit,datnew,type='probs')
  PY1.SZ0<-PY1.SZ1<-PS<-rep(0,length(Su))
  for (i in 1:length(Wu)){
    PY1.SZ1=PY1.SZ1+R.logit(Z=1,S1=Su,W=Wu[i],beta=beta,short)*fx.s[,i]
    PY1.SZ0=PY1.SZ0+R.logit(Z=0,S1=Su,W=Wu[i],beta=beta,short)*fx.s[,i]
  }
  
  VE.S<-1-PY1.SZ1/PY1.SZ0
  
  return(list(beta=beta,VE=VE.S))
}

########################################################################################

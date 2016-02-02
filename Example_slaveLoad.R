# saveLoad file, see master files for detailed description
# HISTORY:
# Date      Programmer   	
# Dec19-2015  Y Zhuang      Version 1
##########################################################################
library("Rmpi")
library("kyotil",lib.loc="~/R/library")
library(nnet)
library(splines)
library(MASS)
library("chngpt",lib.loc="~/R/library")

dirName <- "."
source(file.path(dirName,"myfunctions.R"))

init.seed <- mpi.comm.rank()*100

bVE.Nohinge<- function (iter) {
  
  oo1<-(1:nrow(dat))[dat$Y==0 & dat$Z==0]
  oo2<-(1:nrow(dat))[dat$Y==1 & dat$Z==0]
  oo3<-(1:nrow(dat))[dat$Y==0 & dat$Z==1]
  oo4<-(1:nrow(dat))[dat$Y==1 & dat$Z==1]
  
  Su<-sort(unique(dat$S1[!is.na(dat$S1)]))
  
  i<-1
  set.seed(init.seed+i)
  oo=c(sample(oo1,replace=T),sample(oo2,replace=T),sample(oo3,replace=T),sample(oo4,replace=T))
  out<-getVE(dat[oo,],maxit=500,epsilon=0.005,short,Su)
  out<-c(out$beta,out$VE)
  
  for (i in 2:iter){
    set.seed(init.seed+i)
    oo=c(sample(oo1,replace=T),sample(oo2,replace=T),sample(oo3,replace=T),sample(oo4,replace=T))
    temp<-getVE(dat[oo,],maxit=500,epsilon=0.005,short,Su)
    out<-rbind(out,c(temp$beta,temp$VE))
  }
  return(out)
}


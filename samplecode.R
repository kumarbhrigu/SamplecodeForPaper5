################################################################################
###         This file includes sample Rcode for the manuscript         #########
# A Study on Censored Data Imputation through Propensity Score Matching in     #
#                   Multi State Frailty Model                                  #
################################################################################
rm(list = ls())
library(SurviMChd)
library(survival)
library(msm)
library(mstate)
library(simsurv)
library(cmprsk)
library(riskRegression)
library(survminer)
library(dplyr)
################################################################################
#-------------------------------------------------------------------------------
# survMC function performs survival analysis using Cox Proportional Hazards with MCMC.
# description of various arguments are as follows:-
# m= Starting column number from where variables of high dimensional data will get selected.
# n= Ending column number till where variables of high dimensional data will get selected.
# Time= Variable/Column name containing the information on duration of survival
# Event= Variable/Column name containing the information of survival event
# chains= Number of chains to perform
# adapt= Number of adaptations to perform
# iter= Number of iterations to perform
# data= High dimensional data having survival duration and event.
#-------------------------------------------------------------------------------
survMC <- function(m,n,Time,Event,chains,adapt,iter,data)
{
  if(Time!="OS"){
    names(data)[names(data) == Time] <- "OS"
  }
  if(Event!="Death"){
    names(data)[names(data) == Event] <- "Death"
  }
  data<-data[order(data$OS),]
  var1 <- colnames(data)
  nr <- nrow(data)
  data1 <- subset(data, Death == 1) #subsetting data with death status = 1
  u <- unique(data1$OS) #creating a vector with unique values of OS
  #adding a condition for censoring time vector to include the last censored patient when censoring = 0
  if ((data$Death[nrow(data)])==0){
    u1<-c(u,data$OS[nrow(data)])
  } else {
    u1 <- u
  }
  u2 <- sort(u1)
  u2[length(u2)]<-u2[length(u2)]+1E-3
  t.len<-(length(u2)-1)
  model_jags <- "
  data{
    # Set up data
  for(i in 1:N) {
    for(j in 1:T) {
    Y[i,j] <- step(obs.t[i] - t[j] + eps)
    dN[i, j] <- Y[i, j] * step(t[j + 1] - obs.t[i] - eps) * fail[i]
    }
  }
  }

  # Model
  model{
  for(i in 1:N){
    betax[i,1] <- 0
    for(k in 2:(p+1)){
      betax[i,k] <- betax[i,k-1] + beta[k-1]*x[i,k-1]
    }
  }
  for(j in 1:T) {
    for(i in 1:N) {
    dN[i, j] ~ dpois(Idt[i, j]) # Likelihood
    Idt[i, j] <- Y[i, j] * exp(betax[i,p+1]) * dL0[j] # Intensity
    }
    dL0[j] ~ dgamma(mu[j], c)
    mu[j] <- dL0.star[j] * c # prior mean hazard
  }
  c <- 0.001
  r <- 0.1
  for (j in 1 : T) {
    dL0.star[j] <- r * (t[j + 1] - t[j])
  }
  for(k in 1:p){
    beta[k] ~ dnorm(0.0,0.000001)
  }
  }"
  
  params <- c("beta","dL0")
  inits <-  function(){list( beta = rep(0,p), dL0 = rep(0.0001,bigt))}
  x2 <- rep(0,nrow(data))
  q <- matrix(nrow=0,ncol=5)
  s <- matrix(nrow=0,ncol=2)
  di <- matrix(nrow=0,ncol=1)
  for(i in m:n){
    x1 <- data[(1:nrow(data)),i]
    x = t(rbind(x1,x2))
    datafi <- list(x=x,obs.t=data$OS,t=u2,T=t.len,N=nrow(data),fail=data$Death,eps=1E-10,p=2)
    jags <- jags.model(textConnection(model_jags),
                       data = datafi,
                       n.chains = chains,
                       n.adapt = adapt)
    samps <- coda.samples(jags, params, n.iter=iter)
    s1 <- summary(samps)
    stats <- s1$statistics[1,c(1:2)]
    s <- rbind(s,stats)
    quan <- s1$quantiles[1,]
    q <- rbind(q,quan)
    d = dic.samples(jags, n.iter=iter)
    meandeviance <- round(sum(d$deviance),2)
    di <- rbind(di,meandeviance)
  }
  results <- cbind(s,q)
  expresults <- exp(results)
  Variables <- names(data)[m:n]
  expresults <- data.frame(Variables,expresults,di)
  colnames(expresults)<-c("Variable","Posterior Means","SD","2.5%","25%","50%","75%","97.5%","DIC")
  rownames(expresults) <- NULL
  return(expresults)
}


################################################################################
#' Function to simulate survival data using coxph and weibull baseline hazard
#' @param parm parameters for distribution of survival time at baseline
#' @param dist distribution of survival time at baseline
#' @param ctime censored time 
#' @param beta regression parameter for cox ph
#' @param gparm frailty parameter
#' @param data covariate data
#' @return simulated survival data
################################################################################
invflsurv<-function(parm,dist,ctime,beta,gparm,data){
  s<-gparm
  n<-nrow(data)
  if(dist=="exponential"){
    u<-runif(n,0,1)
    lu<-log(u)
    data1<-as.matrix(data[colnames(data)!="id"])
    xbeta<-exp(data1%*%beta)
    #v<-rgamma(n,s)
    v<-rnorm(n,1,s)
    vxbeta<-as.vector(v%*%xbeta)
    st<--(lu)/(parm*vxbeta)
  }else if( dist=="weibull"){
    u<-runif(n,0,1)
    lu<-log(u)
    data1<-as.matrix(data[colnames(data)!="id"])
    xbeta<-exp(data1%*%beta)
    #v<-rgamma(n,s)
    v<-rnorm(n,1,s)
    vxbeta<-as.vector(v%*%xbeta)
    st<-(-(lu)/(parm[1]*vxbeta))^(1/parm[2])
  }else{
    u<-runif(n,0,1)
    lu<-log(u)
    data1<-as.matrix(data[colnames(data)!="id"])
    xbeta<-exp(data1%*%beta)
    #v<-rgamma(n,s)
    v<-rnorm(n,0,s)
    vxbeta<-as.vector(v%*%xbeta)
    st<-(1/parm[2])*log(((-parm[2]*lu)/(parm[1]*vxbeta))+1)
  } 
  dt<-list()
  ctime<-runif(n,0,1)
  dt$id<-data[colnames(data)=='id']
  dt$eventtime<-st
  dt$status<-as.numeric(st<ctime)
  dt<-cbind(data.frame(dt),data[colnames(data)!='id'])
  dt         
}
################################################################################
#              method argument can take different value like 
#   "euclidean", "maximum", "manhattan", "canberra", "binary" or "minkowski"
################################################################################
#' dscore: Function of Propensity Score Matching for censored survival data   
#' @param status status column name from the input dataset
#' @param data dataset
#' @param prob threshhold probability
#' @param m 1 term of the sequence of variable used in the PSM 
#' @param n last term of the sequence of variable to be used in PSM
#' @param method the distance metric to be used for Propensity score matching 
#' @return updated dataset with 
#' @examples
################################################################################
dscore<-function(status,data,prob,m,n,method="euclidean"){
  # dataset must contain two columns named status,time 
  mdata<-data
  mdata[is.na(mdata)]<-0
  k<-which(colnames(mdata)==status)
  mdata1<-mdata[,m:n]
  mdata2<-cbind(status=mdata$status,time=mdata$time,mdata[,m:n])
  d11<-subset(mdata2,mdata2$status==1)
  d00<-subset(mdata2,mdata2$status==0)
  coefficient<-c()
  for(i in 1:ncol(mdata2)-2){
    fit<-survMC(m=i+2,n=i+2,Time="time",Event="status",chains=2,adapt=100,iter=50,data=mdata2)
    if(fit$`Posterior Means`&fit$`2.5%`&fit$`97.5%`>1){
      coefficient[i]<-fit$`Posterior Means`
    }else if(fit$`Posterior Means`&fit$`2.5%`&fit$`97.5%`<1){
      coefficient[i]<-fit$`Posterior Means`
    } else {
      coefficient[i]<-0
    }
  } 
  #parm<-summary(model1)$coefficients[,1]
  d1<-d11[,-which(colnames(d11)=="status")]
  t<-which(colnames(d1)=="time")
  d1<-d1[,-t]
  d0<-d00[,-which(colnames(d11)=="status")]
  xmin<-c()
  k<-length(coefficient)
  
  xmin<-function(x,data){
    xmin<-c()
    for(i in 1:length(x)){
      if(x[i]>0){
        xmin[i]<-min(data[,i])
      }else{xmin[i]<-max(data[,i])}
    }
    xmin
  }
  xmin1<-xmin(x=coefficient,data=d1)
  
  
  Dscore<-c()
  for( i in 1:nrow(d0)){
    if(anyNA(xmin1)==T){
      k<-which(is.na(xmin1)==T)
    }
    d<-rbind(d0[i,][-k],xmin1[-k])
    Dscore[i]<-dist(d,method=method)
  }  
  rslt<-list() 
  rslt$dscore<-Dscore
  rslt$olddeath<-nrow(d11)
  risk<-ecdf(Dscore)
  proximity<-1-risk(Dscore)
  nd00<-cbind(d00,proximity)
  nstatus<-c()
  for(i in 1:nrow(nd00)){
    if(nd00$proximity[i]>prob){
      nstatus[i]<-1
    }else{
      nstatus[i]<-0
    }
  }
  d11<-subset(mdata,mdata$status==1)
  d00<-subset(mdata,mdata$status==0)
  d00$status<-as.factor(nstatus)
  newdata<-rbind(d11,d00)
  newdata<-newdata[order(newdata$id),]
  rslt$updateddeath<-nrow(subset(newdata,newdata$status=='1'))
  rslt$newdata<-newdata
  #-------------------------------------------------------------------------
  rslt
}
################################################################################
#' Simulate survival data with frailty terms
#' @param p1 lambda is the scale parameter of baseline distribution that is Weibul 
#' @param p2 sigma^2 is the variance for frailty distribution
#' @param p3 gamma is the shape parameter of baseline distribution that is Weibul
#' @param x covariates for the coxph
#' @param beta regression coefficient for the coxph
#' @param n number of datapoints to be simulated
#' @return simulated survival data using given parameter 
#' @examples
################################################################################
ftime<-function(p1,p2,p3,x,beta,n){
  #p1 is lambda
  #p2 is sigma^2
  #p3 is gamma
  n=n
  #p1<-0.029;p2<-0.254;p3<-0.951
  p1<-p1;p2<-p2;p3<-p3
  x<-as.matrix(x)
  beta<-beta
  u<-runif(nrow(x),0,1)
  xb<-x%*%beta
  exb<-exp(xb)
  u2<-(1/u)^{p2}-1
  t<-(u2/(p1*p2*exb))^{1/p3}
  st<-rbinom(n,1,0.8)
  fdata<-data.frame(id=1:n,time=t,status=st,x)
  return(fdata)
}
################################################################################
#################### Sample code for Simulation study ##########################
tmat <- transMat(x = list(c(2, 3), c(3), 
                          c()), names = c("Tx", "Rec", "Death"))
cph1gCoeff<-list();cph2gCoeff<-list();cph3gCoeff<-list();cph1Theta<-list();cph2Theta<-list();cph3Theta<-list()
frl1<-list();frl2<-list();frl3<-list()
cph1g0Coeff<-list();cph2g0Coeff<-list();cph3g0Coeff<-list()
cph10Theta<-list();cph20Theta<-list();cph30Theta<-list()

for(j in 1:5){
  covs<-data.frame(x1=rnorm(100,0,1),x2=rnorm(100,10,20),x3=rnorm(100,5,4),x4=rnorm(100,3,2))
  ssim13<-ftime(p1=0.029,p2=0.254,p3=0.951,beta=t(t(c(0.013,0.045,0.011,0.033))),n=100,x=covs)
  ssim23<-ftime(p1=0.020,p2=0.120,p3=0.6,beta=t(t(c(0.045,0.05,0.06,0.077))),n=100,x=covs)
  stime13<-ssim13$time
  stime23<-ssim23$time
  #sstatus12<-ssim12$status
  sstatus13<-ssim13$status
  sstatus23<-ssim23$status
  data1<-data.frame(id=1:100, stime13,stime23,sstatus13,sstatus23,covs)
  for(i in 1:nrow(data1)){
    if(data1[i,]$sstatus13==0&data1[i,]$sstatus23==0){
      data1[i,]$stime23=data1[i,]$stime13 
    }else if(data1[i,]$sstatus13==1&data1[i,]$sstatus23==0) {
      data1[i,]$stime23=data1[i,]$stime13+data1[i,]$stime23
    } else if(data1[i,]$sstatus13==1&data1[i,]$sstatus23==1){
      data1[i,]$stime23=data1[i,]$stime13+data1[i,]$stime23
    } else if(data1[i,]$sstatus13==0&data1[i,]$sstatus23==1){
      data1[i,]$stime13=data1[i,]$stime23
    }
  }
  udata1<-list()
  udata1$id<-data1$id;udata1$status<-data1$sstatus23;udata1$time<-data1$stime23;udata1$x1<-data1$x1;udata1$x2<-data1$x2
  udata1$x3<-data1$x3;udata1$x4<-data1$x4
  udata1<-as.data.frame(udata1)
  udata<-dscore(status="status",data=udata1,prob=0.9,m=4,n=7)
  ndata<-udata$newdata
  data1$arm<-rbinom(100,1,0.6)
  data2<-data1
  data2$sstatus23<-as.numeric(as.character(ndata$status))
  covs1<-c("x1", "x2", "x3","x4","arm")
  msbmt<-msprep(time = c(NA, "stime13", "stime23"), status = c(NA,
                                                                 "sstatus13", "sstatus23"), data =data1, trans = tmat, keep = covs1)
  
  msbmt1<-msprep(time = c(NA, "stime13", "stime23"), status = c(NA,
                                                                  "sstatus13", "sstatus23"), data =data2, trans = tmat, keep = covs1)
  
  msbmt<-expand.covs(msbmt, covs1, append = TRUE, longnames = FALSE)
  msbmt1<-expand.covs(msbmt1, covs1, append = TRUE, longnames = FALSE)
  
  cph1g<-coxph(Surv(time,status)~x1+x2+x3+x4+frailty(id,distribution = 'gaussian'),data=msbmt[msbmt$trans==1,])
  cph2g<-coxph(Surv(time,status)~x1+x2+x3+x4+frailty(id,distribution = 'gaussian'),data=msbmt[msbmt$trans==2,])
  cph3g<-coxph(Surv(time,status)~x1+x2+x3+x4+frailty(id,distribution = 'gaussian'),data=msbmt[msbmt$trans==3,])
  
  cph10g<-coxph(Surv(time,status)~x1+x2+x3+x4+frailty(id,distribution = 'gaussian'),data=msbmt1[msbmt1$trans==1,])
  cph20g<-coxph(Surv(time,status)~x1+x2+x3+x4+frailty(id,distribution = 'gaussian'),data=msbmt1[msbmt1$trans==2,])
  cph30g<-coxph(Surv(time,status)~x1+x2+x3+x4+frailty(id,distribution = 'gaussian'),data=msbmt1[msbmt1$trans==3,])
  
  cph1gCoeff[[j]]<-cph1g$coefficients
  cph2gCoeff[[j]]<-cph2g$coefficients
  cph3gCoeff[[j]]<-cph3g$coefficients
  cph1Theta[[j]]<-cph1g$history$f$theta
  cph2Theta[[j]]<-cph2g$history$f$theta
  cph3Theta[[j]]<-cph3g$history$f$theta
  
  cph1g0Coeff[[j]]<-cph10g$coefficients
  cph2g0Coeff[[j]]<-cph20g$coefficients
  cph3g0Coeff[[j]]<-cph30g$coefficients
  cph10Theta[[j]]<-cph10g$history$f$theta
  cph20Theta[[j]]<-cph20g$history$f$theta
  cph30Theta[[j]]<-cph30g$history$f$theta
  #############################################################################
  
}
simreslt1<-Reduce("cbind",cph1gCoeff)
simreslt2<-Reduce("cbind",cph2gCoeff)
simreslt3<-Reduce("cbind",cph3gCoeff)
simreslt11<-Reduce('cbind',cph1Theta)
simreslt22<-Reduce('cbind',cph2Theta)
simreslt33<-Reduce('cbind',cph3Theta)

simr2<-list()
simr2$Mean12<-apply(simreslt1,1,mean)
simr2$Mean13<-apply(simreslt2,1,mean)
simr2$Mean23<-apply(simreslt3,1,mean)
simr2$var12<-apply(simreslt1,1,var)
simr2$var13<-apply(simreslt2,1,var)
simr2$var23<-apply(simreslt3,1,var)

simr2$Theta12<-mean(simreslt11)
simr2$Theta13<-mean(simreslt22)
simr2$Theta23<-mean(simreslt33)
simr2$frlvar12<-mean(frlvar1,na.rm = TRUE)
simr2$frlvar13<-mean(frlvar2,na.rm=TRUE)
simr2$frlvar23<-mean(frlvar3,na.rm=TRUE)

simreslt110<-Reduce("cbind",cph1g0Coeff)
simreslt220<-Reduce("cbind",cph2g0Coeff)
simreslt330<-Reduce("cbind",cph3g0Coeff)
simreslt111<-Reduce('cbind',cph10Theta)
simreslt222<-Reduce('cbind',cph20Theta)
simreslt333<-Reduce('cbind',cph30Theta)

simr0<-list()
simr0$Mean12<-apply(simreslt11,1,mean)
simr0$Mean13<-apply(simreslt22,1,mean)
simr0$Mean23<-apply(simreslt33,1,mean)
simr0$var12<-apply(simreslt111,1,var)
simr0$var13<-apply(simreslt222,1,var)
simr0$var23<-apply(simreslt333,1,var)

simr2$Theta12<-mean(simreslt11)
simr2$Theta13<-mean(simreslt22)
simr2$Theta23<-mean(simreslt33)
simr2$frlvar12<-mean(frlvar1,na.rm = TRUE)
simr2$frlvar13<-mean(frlvar2,na.rm=TRUE)
simr2$frlvar23<-mean(frlvar3,na.rm=TRUE)

################################################################################
####################### ebmt3 data analysis ####################################
################################################################################ 
rm(list = ls())
library(SurviMChd)
library(survival)
library(msm)
library(mstate)
library(simsurv)
library(simsurv)
library(cmprsk)
library(riskRegression)
library(survminer)
library(dplyr)
data('ebmt3')
ebmt3<-ebmt3
ebmt3$arm<-rbinom(2204,1,0.7)
ebmt3$x1<-rnorm(2204,0,1)
ebmt3$x2<-rnorm(2204,3,5)
ebmt3$x3<-rnorm(2204,10,3)
ebmt3$x4<-rnorm(2204,5,6)
ebmt3new<-list()
ebmt3new$id<-1:2204
ebmt3new$status<-ebmt3$rfsstat;ebmt3new$time<-ebmt3$rfstime;ebmt3new$x1<-ebmt3$x1;
ebmt3new$x2<-ebmt3$x2;ebmt3new$x3<-ebmt3$x3;ebmt3new$x4<-ebmt3$x4
ebmt3new<-as.data.frame(ebmt3new)
udata1<-dscore(status="status",data=ebmt3new[1:500,],prob=0.7,m=4,n=7)
udata2<-dscore(status="status",data=ebmt3new[501:1000,],prob=0.7,m=4,n=7)
udata3<-dscore(status="status",data=ebmt3new[1001:1500,],prob=0.7,m=4,n=7)
udata4<-dscore(status="status",data=ebmt3new[1501:2204,],prob=0.7,m=4,n=7)



ndata<-rbind(udata1$newdata,udata2$newdata,udata3$newdata,udata4$newdata)
ebmt3d<-ebmt3
ebmt3d$rfsstat<-as.numeric(as.character(ndata$status))
tmat <- transMat(x = list(c(2, 3), c(3), c()), names = c("Tx", "Rec", "Death"))
covs <- c("dissub", "age", "drmatch", "tcd", "prtime","arm","x1","x2","x3","x4")
msbmt <- msprep(time = c(NA, "prtime", "rfstime"), status = c(NA,"prstat", "rfsstat"),
                data = ebmt3, trans = tmat, keep = covs)
msbmt1 <- msprep(time = c(NA, "prtime", "rfstime"), status = c(NA,"prstat", "rfsstat"),
                 data = ebmt3d, trans = tmat, keep = covs)
msbmt <- expand.covs(msbmt, covs, append = TRUE, longnames = FALSE)
msbmt1 <- expand.covs(msbmt1, covs, append = TRUE, longnames = FALSE)


msph1<-coxph(Surv(time,status)~dissub+age +drmatch+ tcd+x1+x2+x3+x4+    
               frailty(id,distribution = 'gaussian'),data=msbmt[msbmt$trans==1,])

msph2<-coxph(Surv(time,status)~dissub+age +drmatch+ tcd+ x1+x2+x3+x4+  
               frailty(id,distribution = 'gaussian'),data=msbmt[msbmt$trans==2,])
msph3<-coxph(Surv(time,status)~dissub+age +drmatch+ tcd+  x1+x2+x3+x4+  
               frailty(id,distribution = 'gaussian'),data=msbmt[msbmt$trans==3,])
msph4<-coxph(Surv(time,status)~ x1.1+x2.1+x3.1+x4.1+x1.2+x2.2+x3.2+x4.2+x1.3+x2.3+x3.3+x4.3+strata(trans)
             ,data=msbmt)

mdata<-list()
mdata$mscoeff1<-round(summary(msph1)$coefficients,4)[,c(1,3,6)]
mdata$mscoeff2<-round(summary(msph2)$coefficients,4)[,c(1,3,6)]
mdata$mscoeff3<-round(summary(msph3)$coefficients,4)[,c(1,3,6)]
mdata$msfrail1<-msph1$history$f$theta
mdata$msfrail2<-msph2$history$f$theta
mdata$msfrail3<-msph3$history$f$theta

msph11<-coxph(Surv(Tstart,Tstop,status)~dissub+age+x1+x2+x3+x4+
                +drmatch+ tcd+frailty(id,distribution = 'gaussian'),data=msbmt1[msbmt1$trans==1,])

msph22<-coxph(Surv(Tstart,Tstop,status)~dissub+age +drmatch+ tcd+ x1+x2+x3+x4+   
                frailty(id,distribution = 'gaussian'),data=msbmt1[msbmt1$trans==2,])
msph33<-coxph(Surv(Tstart,Tstop,status)~dissub+age +drmatch+ tcd+x1+x2+x3+x4+    
                frailty(id,distribution = 'gaussian'),data=msbmt1[msbmt1$trans==3,])

mdata1<-list()
mdata1$mscoeff1<-round(summary(msph11)$coefficients,4)[,c(1,3,6)]
mdata1$mscoeff2<-round(summary(msph22)$coefficients,4)[,c(1,3,6)]
mdata1$mscoeff3<-round(summary(msph33)$coefficients,4)[,c(1,3,6)]
mdata1$msfrail1<-msph11$history$f$theta
mdata1$msfrail2<-msph22$history$f$theta
mdata1$msfrail3<-msph33$history$f$theta

bmsph1<-basehaz(msph1,centered=F)
bmsph2<-basehaz(msph2,centered=F)
bmsph3<-basehaz(msph3,centered=F)
bmsph11<-basehaz(msph11,centered=F)
bmsph22<-basehaz(msph22,centered=F)
bmsph33<-basehaz(msph33,centered=F)
par(mfrow=c(3,1))
hist(exp(msph11$frail),main='1->2',xlab='Estimated Frailty')
lines(density(exp(msph11$frail)))
hist(exp(msph22$frail),main='1->3',xlab='Estimated Frailty')
lines(density(exp(msph22$frail)))
hist(exp(msph33$frail),main='2->3',xlab='Estimated Frailty')
lines(density(exp(msph33$frail)))

################################################################################
########### Predicted survival probability Plot before and after PSM ###########
########### for each transition 1->2,1->3,2->3 #################################
################################################################################
sProbMax1<-exp(-5*exp(as.vector(c(0,1,1,0,1,1,0.05,0.45,-0.76,0.8))%*%t(t(c(0.088	,-0.0694,	0.0511,	0.3118,	-0.1068,	0.2083,	-0.117,	0.0089,	-0.0065,	-0.0267
))))*bmsph1$hazard)
sProbMin1<-exp(-0.5*exp(as.vector(c(0,1,1,0,1,1,0.05,0.45,-0.76,0.8))%*%t(t(c(0.088	,-0.0694,	0.0511,	0.3118,	-0.1068,	0.2083,	-0.117,	0.0089,	-0.0065,	-0.0267))))*bmsph1$hazard)
sProbMax0<-exp(-5*exp(as.vector(c(0,1,1,0,1,1,0.05,0.45,-0.76,0.8))%*%t(t(c(0.1129,	0.0089,	0.0404,	0.4171,	0.2243,	-0.0679,	-0.0161,	0.0049,	0.0011,	-0.0454
))))*bmsph11$hazard)
sProbMin0<-exp(-0.5*exp(as.vector(c(0,1,1,0,1,1,0.05,0.45,-0.76,0.8))%*%t(t(c(0.1129,	0.0089,	0.0404,	0.4171,	0.2243,	-0.0679,	-0.0161,	0.0049,	0.0011,	-0.0454))))*bmsph11$hazard)
plot(bmsph22$time,sProbMax1,type = 'l',ylim=c(0,1),xlim=c(1,200),xlab='Time',ylab='Survival probability')
lines(bmsph22$time,sProbMin1,type = 'l',lty=2)
lines(bmsph22$time,sProbMax0,type = 'l',col='red')
lines(bmsph22$time,sProbMin0,type = 'l',col='red',lty=2)
legend("topright", legend=c("2->3", "1->3"),
       col=c('black', 'red'), lty=1:2, cex=0.8)

sProbMax1<-exp(-exp(as.vector(c(0,1,1))%*%t(t(c(.1762,0.0579,-.1019)))+1.2)*bmsph2$hazard)
sProbMin1<-exp(-exp(as.vector(c(0,1,1))%*%t(t(c(.1762,0.0579,-.1019)))-0.1099327)*bmsph2$hazard)
sProbMax0<-exp(-exp(as.vector(c(0,1,0))%*%t(t(c(.1762,0.0579,-.1019)))+1.2)*bmsph2$hazard)
sProbMin0<-exp(-exp(as.vector(c(0,1,0))%*%t(t(c(.1762,0.0579,-.1019)))-0.1099327)*bmsph2$hazard)
plot(bmsph2$time,sProbMax1,type = 'l',ylim=c(0,1),xlim=c(1,500),xlab='Time',ylab='Survival probability')
lines(bmsph2$time,sProbMin1,type = 'l',lty=2)
lines(bmsph2$time,sProbMax0,type = 'l',col='red')
lines(bmsph2$time,sProbMin0,type = 'l',col='red',lty=2)
legend("topright", legend=c("arm=1", "arm=0"),
       col=c('black', 'red'), lty=1:2, cex=0.8)

sProbMax1<-exp(-exp(as.vector(c(0,1,1))%*%t(t(c(.0166,0.3478,0.2644)))+1.2)*bmsph3$hazard)
sProbMin1<-exp(-exp(as.vector(c(0,1,1))%*%t(t(c(.0166,0.3478,0.2644)))-0.1099327)*bmsph3$hazard)
sProbMax0<-exp(-exp(as.vector(c(0,1,0))%*%t(t(c(.0166,0.3478,0.2644)))+1.2)*bmsph3$hazard)
sProbMin0<-exp(-exp(as.vector(c(0,1,0))%*%t(t(c(.0166,0.3478,0.2644)))-0.1099327)*bmsph3$hazard)
plot(bmsph3$time,sProbMax1,type = 'l',ylim=c(0,1),xlim=c(1,500),xlab='Time',ylab='Survival probability')
lines(bmsph3$time,sProbMin1,type = 'l',lty=2)
lines(bmsph3$time,sProbMax0,type = 'l',col='red')
lines(bmsph3$time,sProbMin0,type = 'l',col='red',lty=2)
legend("topright", legend=c("arm=1", "arm=0"),
       col=c('black', 'red'), lty=1:2, cex=0.8)




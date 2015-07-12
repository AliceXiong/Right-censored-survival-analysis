library(gamlss)
library(gamlss.dist)
library(gamlss.cens)
library(eha)
library(survival)
library(Hmisc)
gen.cens(WEI,type="right") #defines a right censored Weibull
gen.cens(EXP,type="right") #defines a right censored exponential
gen.cens(LOGNO,type="right") #defines a right censored lognomral
gen.cens(GA,type="right") #defines a right censored gamma

##############################################################################
effect.size<-function(e.size,samp.size1,samp.size2){ 
  #weibull,effect size=0
  mu<-5
  sigma<-2
  WEI.sd<-sqrt((mu^2)*(gamma((2/sigma)+1)-gamma((1/sigma)+1)^2))
  level1<-rWEI(n=samp.size1,mu,sigma) 
  level2<-rWEI(n=samp.size2,mu,sigma)
  #weibull,censor data,effect size=0
  cen1<-rWEI(n=samp.size1,mu,sigma) 
  cen2<-rWEI(n=samp.size2,mu,sigma)
  
  y1<-pmin(level1,cen1)
  c1<-as.numeric(level1<=cen1)
  y2<-pmin(level2,cen2)
  c2<-as.numeric(level2<=cen2)
  c<-c(c1,c2)
  g<-c(rep(1,samp.size1),rep(2,samp.size2))
  y<-c(y1+e.size*WEI.sd,y2)
  
  #effect size 0, creates a right-censored survival object
  right<-Surv(y,c)
  
  #three traditional method
  #logrank, nonparametric
  t2<-survdiff(right~g)
  p2<-pchisq(t2$chisq,df=1,lower.tail=F)
  
  #Cox regression, semiparametric
  t3<-coxph(right~g)
  t3ll<-t3$loglik[1]-t3$loglik[2] #use likelihood ratio test statistic
  p3<-pchisq(-2*t3ll,df=1,lower.tail=F)
  
  #Weibull regression, fully parametric
  t4<-weibreg(right~g)
  t4ll<-t4$loglik[1]-t4$loglik[2] #use likelihood ratio test statistic
  p4<-pchisq(-2*t4ll,df=1,lower.tail=F)
  
  #####the method that we try to compare,GAMLSS Model 1, change family
  t5<-gamlss(right~g,family=LOGNOrc) #alt hyp model
  t0<-gamlss(right~1,family=LOGNOrc) #null hyp model
  #when deviance(null)<deviance(alt), then problem as LR stat is negative, so set p-value=1
  #deviance=-2LL
  p5<-ifelse(deviance(t0)>deviance(t5),LR.test(t0,t5,print=FALSE)$p.val,1) #likelihood ratio
  
  #2nd GAMLSS model also models the variance term, change family
  t6<-gamlss(right~g,sigma.formula=~g,family=LOGNOrc)
  p6<-ifelse(deviance(t0)>deviance(t6),LR.test(t0,t6,print=FALSE)$p.val,1)
  ####return the p-value in the effect size of zero, which means that there is no difference between these two groups in the variance.(the same mean) 
  ####then, we need to repeat this situaiton for 5 times in the different variances by malnicate the mean to make it keep the same.  
  return(c(p2,p3,p4,p5,p6))}
###############################################################################


samp.size1<-25
samp.size2<-50

####################################################################################
NewSimu<-function(samp.size1,samp.size2){
  a<-effect.size(0,samp.size1,samp.size2)
  b<-effect.size(0.1,samp.size1,samp.size2)
  c<-effect.size(0.3,samp.size1,samp.size2)
  d<-effect.size(0.5,samp.size1,samp.size2)
  e<-effect.size(1,samp.size1,samp.size2)
  return(c(a,b,c,d,e))
}
all<-NewSimu(samp.size1,samp.size2)
all.p<-matrix(all,nrow=5,byrow=T)
all.p
#############################################################################


alpha<-0.05
effect<-c(0,0.1,0.3,0.5,1)
NewSimu(samp.size1,samp.size2)
#set.seed(800)#gamma model
#set.seed(900)#weibull model
set.seed(1000)#LOGNO model

result<-matrix(NA,nrow=1000,ncol=25)
result
for(i in 1:100) result[i,]<-try(NewSimu(samp.size1,samp.size2),silent=TRUE)
for(i in 101:200) result[i,]<-try(NewSimu(samp.size1,samp.size2),silent=TRUE)
for(i in 201:300) result[i,]<-try(NewSimu(samp.size1,samp.size2),silent=TRUE)
for(i in 301:400) result[i,]<-try(NewSimu(samp.size1,samp.size2),silent=TRUE)
for(i in 401:500) result[i,]<-try(NewSimu(samp.size1,samp.size2),silent=TRUE)
for(i in 501:600) result[i,]<-try(NewSimu(samp.size1,samp.size2),silent=TRUE)
for(i in 601:700) result[i,]<-try(NewSimu(samp.size1,samp.size2),silent=TRUE)
for(i in 701:800) result[i,]<-try(NewSimu(samp.size1,samp.size2),silent=TRUE)
for(i in 801:900) result[i,]<-try(NewSimu(samp.size1,samp.size2),silent=TRUE)
for(i in 901:1000) result[i,]<-try(NewSimu(samp.size1,samp.size2),silent=TRUE)
result2<-matrix(ncol=25,as.numeric(result))

for(i in 1:1000) if (is.na(result2[i,1]))
  result2[i,]<-try(NewSimu(samp.size), silent=TRUE)
#workaround to get rid of NAs when crash
power2<-sum(result2[,1]<=alpha,na.rm=T) #logrank
power3<-sum(result2[,2]<=alpha,na.rm=T) #Cox regression
power4<-sum(result2[,3]<=alpha,na.rm=T) #Weibull regression
power5<-sum(result2[,4]<=alpha,na.rm=T) #GAMLSS1 (models mean only)
power6<-sum(result2[,5]<=alpha,na.rm=T) #GAMLSS2 (models mean and variance)
power2l<-sum(result2[,6]<=alpha,na.rm=T) #logrank
power3l<-sum(result2[,7]<=alpha,na.rm=T) #Cox regression
power4l<-sum(result2[,8]<=alpha,na.rm=T) #Weibull regression
power5l<-sum(result2[,9]<=alpha,na.rm=T) #GAMLSS1 (models mean only)
power6l<-sum(result2[,10]<=alpha,na.rm=T) #GAMLSS2 (models mean and variance)
power2m<-sum(result2[,11]<=alpha,na.rm=T) #logrank
power3m<-sum(result2[,12]<=alpha,na.rm=T) #Cox regression
power4m<-sum(result2[,13]<=alpha,na.rm=T) #Weibull regression
power5m<-sum(result2[,14]<=alpha,na.rm=T) #GAMLSS1 (models mean only)
power6m<-sum(result2[,15]<=alpha,na.rm=T) #GAMLSS2 (models mean and variance)
power2h<-sum(result2[,16]<=alpha,na.rm=T) #logrank
power3h<-sum(result2[,17]<=alpha,na.rm=T) #Cox regression
power4h<-sum(result2[,18]<=alpha,na.rm=T) #Weibull regression
power5h<-sum(result2[,19]<=alpha,na.rm=T) #GAMLSS1 (models mean only)
power6h<-sum(result2[,20]<=alpha,na.rm=T) #GAMLSS2 (models mean and variance)
power2s<-sum(result2[,21]<=alpha,na.rm=T) #logrank
power3s<-sum(result2[,22]<=alpha,na.rm=T) #Cox regression
power4s<-sum(result2[,23]<=alpha,na.rm=T) #Weibull regression
power5s<-sum(result2[,24]<=alpha,na.rm=T) #GAMLSS1 (models mean only)
power6s<-sum(result2[,25]<=alpha,na.rm=T) #GAMLSS2 (models mean and variance)
N<-nrow(result2)
power.zero<-c(power2,power3,power4,power5,power6)
power.zero<-power.zero/N
power.low<-c(power2l,power3l,power4l,power5l,power6l)
power.low<-power.low/N
power.med<-c(power2m,power3m,power4m,power5m,power6m)
power.med<-power.med/N
power.high<-c(power2h,power3h,power4h,power5h,power6h)
power.high<-power.high/N
power.super<-c(power2s,power3s,power4s,power5s,power6s)
power.super<-power.super/N

final<-c(samp.size,alpha)
names(final)<-c("Distribution","mu","sigma","n","alpha")
power<-matrix(NA,nrow=5,ncol=6)
power[1,]<-c(effect[1],power.zero)
power[2,]<-c(effect[2],power.low)
power[3,]<-c(effect[3],power.med)
power[4,]<-c(effect[4],power.high)
power[5,]<-c(effect[5],power.super)
colnames(power)<-c("effect","logrank","Cox","Weibull","GAMLSS1","GAMLSS2")
power
#creates table of results
resultz<-list(t(final),power)
resultz

#check to make sure that all 1000 models were fit
N
#Edit this code to fit different distributions accordingly

#Figure (Figure5)
plot(power[,2]~power[,1],col=1,lty=1,lwd=2,type="l",xlab="effect",ylab="power")
points(power[,3]~power[,1],col=2,lty=2,lwd=2,type="l")
points(power[,4]~power[,1],col=3,lty=3,lwd=2,type="l")
points(power[,5]~power[,1],col=4,lty=4,lwd=2,type="l")
points(power[,6]~power[,1],col=5,lty=5,lwd=2,type="l")
title("0-1 LOGNO Model to Weibull Data, n=25,50, Unequal means")
legend("bottomright",legend=c("Logrank","Cox","Weibull","GAMLSS1","GAMLSS 2"),col=1:5,lty=1:5,lwd=c(2,2,2,2,2))




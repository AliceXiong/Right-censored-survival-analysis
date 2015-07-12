

library(gamlss)
library(gamlss.dist)
library(gamlss.cens)
library(eha)
library(survival)
library(Hmisc)
gen.cens(WEI,type="right") #defines a right censored Weibull

attach(veteran)
veteran
vet.surv<-with(veteran,Surv(time,status))

###logrank test, p-value
vet.logrank<-survdiff(vet.surv~trt, data=veteran)
summary(vet.logrank)
p.logrank<-pchisq(vet.logrank$chisq,df=1,lower.tail=F)
p.logrank

###cox regression, p-value
vet.cox<-coxph(vet.surv~trt, data=veteran)
summary(vet.cox)
coxll<-vet.cox$loglik[1]-vet.cox$loglik[2]
p.cox<-pchisq(-2*coxll,df=1,lower.tail=F)
p.cox

###gamlss, p-value
gen.cens(WEI,type="right")
vet.gamlss1<-gamlss(Surv(time,status)~trt, data=veteran, family=WEIrc)
summary(vet.gamlss)
vet.gamlss0<-gamlss(Surv(time,status)~1, data=veteran, family=WEIrc)
p.gamlss<-ifelse(deviance(vet.gamlss0)>deviance(vet.gamlss1), LR.test(vet.gamlss0,vet.gamlss1,print=FALSE)$p.val,1)
p.gamlss

###weibull regression, p-value
detach(veteran)

veteran$c1<-ifelse(veteran$trt=="1",1,0)
veteram$c2<-ifelse(veteran$trt=="2",1,0)

attach(veteran)
S<-Surv(time,status)

wei<-weibreg(S~c2)
summary(wei)

weill<-wei$loglik[1]-wei$loglik[2]
p.wei<-pchisq(-2*weill,df=1,lower.tail=F)
p.wei



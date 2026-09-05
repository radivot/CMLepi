## 3_workflow.R
graphics.off();rm(list=ls())#clear plots and environment 
library(tidyverse)
# library(SEERaBomb) #call funcs with ::
library(mgcv)
options(rgl.useNULL = TRUE)
options(rgl.printRglwidget = TRUE)
# library(rgl)  #call funcs with ::
load("~/data/CMLepi/cml20.RData") #made in mkSEER.R
(d=d20|>filter(histo3%in%c(9863,9875),agedx<90,surv<80,surv>0)) #43932
load("~/data/CMLepi/Gac.RData") #made in mkMorts.R
Sex="Both"
head(PYin<-d|>mutate(py=surv,age=agedx,year=yrdx)|>select(py,age,year))
head(PYin<-as.matrix(PYin))
yrs=1975:2023
ages=0.5:125.5
Z=matrix(0,ncol=length(yrs),nrow=length(ages))
colnames(Z)=yrs
rownames(Z)=ages
head(Z)# all zeros initially
PY=Z+0 # add 0 to make sure not a lazy copy of just the pointer 
SEERaBomb::fillPYM(PYin, PY)
(d=d|>mutate(lc=ifelse(COD2!="alive",1,0)))
N=1e4 ## shrink by N big then multiply by it later, so each PY of 1 lands in one bin
head(Od<-d|>mutate(py=lc/N,age=agedx+surv,year=yrdx+surv)|>select(py,age,year)) 
head(Oin<-as.matrix(Od))
O=Z+0
SEERaBomb::fillPYM(Oin, O)
O=O*N
head(dO<-reshape2:::melt(O,value.name="Obs"))
head(dP<-reshape2:::melt(PY,value.name="PY"))
head(dd<-left_join(dO,dP))
names(dd)[1:2]<-c("age","year")
(D=dd|>filter(PY>0))
D$sex=Sex
D=D|>mutate(denom=PY,num=Obs)
D=D|>filter(age<90,age>20)
D$Eus=exp(predict(Gac,D))  # Gac based on US rate

# Using rgl package, get it in  plot tab, rotate view, and save as png
with(D,rgl::plot3d(age,year,PY,xlab="",ylab="",zlab="",alpha=1,type="p"))
with(D,rgl::plot3d(age,year,Eus,xlab="",ylab="",zlab="",alpha=1,type="p")) 
with(D,rgl::plot3d(age,year,Obs,xlab="",ylab="",zlab="",alpha=1,type="h"))

EARbrks=c(0.02,0.03,0.1,0.3)
(Drr=D|>group_by(year)|>summarise(O=sum(Obs),E=sum(Eus),PY=sum(denom)))
(Drr=Drr|>mutate(EAR=(O-E)/PY,LL=EAR-1.96*sqrt(O)/PY,UL=EAR+1.96*sqrt(O)/PY))

library(splines)
DF=3
D12=Drr|>filter(year<=2013)
summary(lmo<-lm(log(EAR)~ns(year,df=DF),data=D12))
pD=data.frame(year=2000:2023,yG="2013")
pD$fit=exp(predict(lmo,newdata=pD))

geE=geom_errorbar(aes(ymin=LL,ymax=UL),width=0.2)#for absolute risks
tc=function(sz) theme_classic(base_size=sz)
ghp2=geom_hline(yintercept=0.02,col="gray")
ghp03=geom_hline(yintercept=0.03,col="gray")
leg=theme(legend.margin=margin(0,0,0,0),legend.position= c(0.75, 0.85))
cc1=coord_cartesian(ylim=c(0.015,NA))
(Drr=Drr|>mutate(yG=ifelse(year<=2013,"2013","2023")))
Drr|>ggplot(aes(x=year,y=EAR,col=yG))+geE+geom_point(size=0.7) + #scale_y_log10(breaks=EARbrks)+
  geom_line(aes(y=fit),data=pD)+cc1+
  scale_y_log10(breaks=EARbrks)+tc(13)+ghp2+ghp03+ylab("Excess Absolute Risk of Death")+leg+labs(x="Year",color="Data through")
ggsave(file="LE/outs/3_workflowEAR.pdf",height=3,width=3) # EAR

(D12=D12|>mutate(RR=O/E,rrL=qchisq(.025,2*O)/(2*E),rrU=qchisq(.975,2*O+2)/(2*E)))
summary(lm1<-lm(log(RR)~ns(year,df=DF),data=D12))
pD=data.frame(year=2000:2023,yG="2013")
pD$fit=exp(predict(lm1,newdata=pD))
(Drr=Drr|>mutate(RR=O/E,rrL=qchisq(.025,2*O)/(2*E),rrU=qchisq(.975,2*O+2)/(2*E)))
gh1=geom_hline(yintercept=1)
gh2=geom_hline(yintercept=2,col="gray")
geR=geom_errorbar(aes(ymin=rrL,ymax=rrU),width=.2,col="gray")
RRbrks=c(0,1,2,5,10,15)
Drr|>ggplot(aes(x=year,y=RR,col=yG))+geR+geom_point(size=0.7) + #scale_y_log10(breaks=EARbrks)+
  geom_line(aes(y=fit),data=pD)+
  scale_y_log10(breaks=RRbrks)+
  tc(13)+gh2+gh1+ylab("Relative Risk of Death")+leg+labs(x="Year",color="Data through")
ggsave(file="LE/outs/3_workflowRR.pdf",height=3,width=3) # EAR

(Ages=seq(min(D$age),max(D$age)))
(Years=seq(min(D$year),max(D$year)))
head(nD<-expand.grid(Ages,Years))
names(nD)<-c("age","year")
nD$denom=1
head(nD)
nD=nD|>mutate(sex=Sex)
nD$E=exp(predict(Gac,nD)) # Gh based on US rate
with(nD,rgl::plot3d(age,year,log10(E),xlab="",ylab="",zlab="",type="n"))
US=reshape2::acast(nD|>filter(sex==Sex)|>select(year,age,E), age~year, value.var="E")
rgl::surface3d(Ages,Years,log10(US),col="violet",alpha=0.5) # mort surface in middle
# save manually from rstudio window 

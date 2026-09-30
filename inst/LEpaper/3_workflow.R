## 3_workflow.R
graphics.off();rm(list=ls())#clear plots and environment 
library(tidyverse)
library(splines)
options(rgl.useNULL = TRUE) # rgl no longer works with xquartz on macs
options(rgl.printRglwidget = TRUE)  #so these two lines are now needed
load("~/data/CMLepi/grp_KK.RData") #run mkGrp.R first to create this 
(D=D|>select(age,year,t,sex,PY,O,E,m)|>filter(sex=="Male"))
D=D|>group_by(age,year,sex,m)|>summarize(O=sum(O),E=sum(E),PY=sum(PY))|>  
    mutate(EAR=(O-E)/PY,LL=EAR-1.96*sqrt(O)/PY,UL=EAR+1.96*sqrt(O)/PY)|>ungroup() 
(D=D|>filter(age<90))
# for rgl plots, from the plot tab, export, save as image, rotate view, save as png 
with(D,rgl::plot3d(age,year,PY,xlab="",ylab="",zlab="",alpha=1,type="p"))
with(D,rgl::plot3d(age,year,E,xlab="",ylab="",zlab="",alpha=1,type="p")) 
with(D,rgl::plot3d(age,year,O,xlab="",ylab="",zlab="",alpha=1,type="h"))

(Dr=D|>group_by(year)|>summarise(O=sum(O),E=sum(E),PY=sum(PY)))
(Dr=Dr|>mutate(EAR=(O-E)/PY,LL=EAR-1.96*sqrt(O)/PY,UL=EAR+1.96*sqrt(O)/PY))
D12=Dr|>filter(year<=2013)
summary(lmo<-lm(log(EAR)~ns(year,df=3),data=D12))
pD=data.frame(year=2000:2023,yG="2013")
pD$fit=exp(predict(lmo,newdata=pD))
geE=geom_errorbar(aes(ymin=LL,ymax=UL),width=0.2)#for absolute risks
tc=function(sz) theme_classic(base_size=sz)
EARbrks=c(0.02,0.03,0.1,0.3)
ghp2=geom_hline(yintercept=0.02,col="gray")
ghp03=geom_hline(yintercept=0.03,col="gray")
leg=theme(legend.margin=margin(0,0,0,0),legend.position= c(0.75, 0.85))
cc1=coord_cartesian(ylim=c(0.015,NA),xlim=c(2000,2025))
(Dr=Dr|>mutate(yG=ifelse(year<=2013,"2013","2023")))
Dr|>ggplot(aes(x=year,y=EAR,col=yG))+geE+geom_point(size=0.7) + #scale_y_log10(breaks=EARbrks)+
  geom_line(aes(y=fit),data=pD)+cc1+
  scale_y_log10(breaks=EARbrks)+tc(13)+ghp2+ghp03+ylab("Excess Absolute Risk of Death")+leg+labs(x="Year",color="Data through")
ggsave(file="LE/outs/3_workflowEARm.pdf",height=3,width=3) # EAR for males

(US=reshape2::acast(D|>select(year,age,m),age~year,value.var="m")) # get into matrix form
(Ages=seq(min(D$age),max(D$age)))
(Years=seq(min(D$year),max(D$year)))
with(D,rgl::plot3d(age,year,log10(m),xlab="",ylab="",zlab="",type="n")) #set up axes without points
rgl::surface3d(Ages,Years,log10(US),col="violet",alpha=0.5) # add AC surface

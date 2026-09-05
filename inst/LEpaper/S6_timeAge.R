## S6_timeAge.R
graphics.off();rm(list=ls())#clear plots and environment 
library(tidyverse)
load("~/data/CMLepi/cml20.RData") #made in mkSEER.R
(d20=d20|>filter(histo3%in%c(9863,9875),agedx<90,surv<80,surv>0)) #43932
d=d20|>select(yrdx,agedx,sex,surv,status) # 43,932 CML cases
(ages=seq(40,90,10))
leg=theme(legend.margin=margin(0,0,0,0),legend.position="top",legend.text=element_text(size=12))
gx=xlab("Years Since Dx")
gl=geom_line()
gp=geom_point()
gh0=geom_hline(yintercept=0)
gh1=geom_hline(yintercept=1)
gh2=geom_hline(yintercept=2,col="gray")
tc=function(sz) theme_classic(base_size=sz)
ccEAR=coord_cartesian(ylim=c(0,NA))
geE=geom_errorbar(aes(ymin=LL,ymax=UL),width=0.2)#for absolute risks
geRR=geom_errorbar(aes(ymin=rrL,ymax=rrU),width=.2)
yint=c(0.01,0.014,0.025,0.06,0.1)
gyE=ylab("Excess Absolute Risk of Death")
gh2=geom_hline(yintercept=2,col="gray")
for (i in 1:5) {
  dt=d|>filter(agedx>=ages[i],agedx<ages[i+1]) #9k
  (txt=paste0("Age at Dx in ",ages[i],"-",ages[i+1]-1))
  load("~/data/mrt/mrtUSA.RData")#loads US mortality data
  (D=SEERaBomb::msd(dt,mrt,brkst=c(0,1,2,3,4,5,6,8,10,12,15)))
  D=D|>rename(Group="sex")%>%select(Group,int,everything())
  (D=SEERaBomb::foldD(D,keep=c("int")))
  D|>ggplot(aes(x=t,y=RR)) +gp+gl+gx+gh1+gh2+tc(13)+leg+scale_y_continuous(breaks=c(1,2,5,10,15))+geRR+
    ggtitle(txt)+  theme(plot.title = element_text(size = 10))+
    scale_color_manual(values = c("gray","black"),name = "Year of Dx")+ylab("Relative Risk of Death")
  ggsave(paste0("LE/outs/supp/S6_RNR",ages[i],".pdf"),width=2.5,height=3)
  sy= scale_y_continuous(breaks=c(0,yint[i],0.025,0.05,0.075,0.1,0.2,0.3))
  ghp2=geom_hline(yintercept=yint[i],col="gray")
  if (i==4) {sy= scale_y_continuous(breaks=c(0,0.06,0.1,0.15))}
  if (i==5) {sy= scale_y_continuous(breaks=c(0,0.1,0.2,0.3))}
  D|>ggplot(aes(x=t,y=EAR))+gp+gl+gx+gyE+gh0+ghp2+tc(13)+leg+geE+ccEAR+
    ggtitle(txt)+  theme(plot.title = element_text(size = 10))+
    scale_x_continuous(breaks=c(0,5,10,15,20,25))+sy
  ggsave(paste0("LE/outs/supp/S6_EAR",ages[i],".pdf"),width=2.5,height=3)
}

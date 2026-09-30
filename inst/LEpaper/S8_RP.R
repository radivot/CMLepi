# Figure 7A RP fit
graphics.off();rm(list=ls())#clear plots and environment 
library(flexsurv) 
library(tidyverse)
load("~/data/CMLepi/cml20.RData") 
(d=d20|>filter(histo3%in%c(9863,9875),agedx<90,surv<80,surv>0)) #43932
range(d$agedx) # 90 is 90+
agedxUL=64  # not too old so we can get out past 10 years
agedxLL=54  # all above 65 so transplants not an issue
agedxGap=10  #
(d=d|>filter(agedx>=agedxLL,agedx<=agedxUL))
d=d|>mutate(ageC=cut(agedx,seq(agedxLL,agedxUL,agedxGap),include.lowest = TRUE, right = FALSE))
(L=split.data.frame(d,d$ageC))
# library(SEERaBomb)
(brks=seq(agedxLL,agedxUL,agedxGap))
(binS=levels(cut(brks+0.1,breaks=c(brks)))) #make a vector of intervals
SEERaBomb::getBinInfo(binS[1],binS)
load("~/data/mrt/mrtUSA.RData")#loads US mortality data mrt
(L=lapply(L,function(x) SEERaBomb::msd(x,mrt,brkst=c(0:21))|>SEERaBomb::foldD(keep=c("int"))))
(D=bind_rows(L,.id="ageGroup"))
D$agedx=sapply(D$ageGroup,function(x) mean(SEERaBomb::getBinInfo(x,binS)[1:2]))
(D=D|>mutate(age=agedx+t))
d=d|>mutate(s=as.character(sex))
d=d|>mutate(a=as.character(round(agedx+surv)))
d=d|>mutate(y=as.character(round(yrdx+surv)))
getMort=function(s,a,y) mrt[[s]][a,y]
d$h=mapply(getMort,d$s,d$a,d$y)
load("~/data/mrt/hBack.RData")  # made in S7_mkhBack.R
(D$y=approx(x=hb$time,y=hb$hazard,xout=D$t)$y)
D=D|>mutate(LLy=LL+y,ULy=UL+y,EARy=EAR+y)
# cc1=coord_cartesian(xlim=c(0,22),ylim=c(0,0.25))
(fsp <- flexsurvspline(Surv(surv, status) ~ 1, data = d,k=2,scale="hazard",bhazard=h))
(h <- as.data.frame(summary(fsp, type = "hazard", t = seq(0.1,25,0.1)))) 
(h$y=approx(x=hb$time,y=hb$hazard,xout=h$time)$y)
(h=h|>mutate(LL=lcl+y,UL=ucl+y,h=est+y))
# source("LE/setup.R")
gx=xlab("Years Since Dx")
gyE=ylab("Excess Absolute Risk of Death")
gp=geom_point()
ghp03=geom_hline(yintercept=c(0.03),col="gray")
gh0=geom_hline(yintercept=0)
tc=function(sz) theme_classic(base_size=sz)
cc1=coord_cartesian(ylim=c(0.0,NA))
leg=theme(legend.margin=margin(0,0,0,0),
          legend.position= c(0.7, 0.85),legend.text=element_text(size=12))
D|>ggplot(aes(x=t,y=EARy))+gp+gx+gyE+gh0+ghp03+tc(13)+leg+
  geom_ribbon(aes(x=time,y=h,ymin=LL, ymax=UL), data=h,alpha=0.2, colour=NA) +
  geom_line(aes(x=time,y=h), data=h) +
  geom_step(aes(x=time,y=hazard),data=hb)+
  theme_minimal() + xlab("Years after diagnosis") + ylab("Hazard") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  geom_errorbar(aes(x=t,y=EARy,ymin=LLy,ymax=ULy),data=D,width=0.2)+ #tc(13)+
  geom_vline(xintercept = exp(fsp$aux$knots), col="gray80", lty=2) + cc1
ggsave("LE/outs/supp/S7_RP.pdf",width=3.5,height=3.5)  

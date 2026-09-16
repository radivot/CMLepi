# testCOD2.R   
graphics.off();rm(list=ls())#clear plots and environment 
library(tidyverse)  
load("~/data/CMLepi/D2.RData") 
Dt=D2|>group_by(t)|>summarize(Olc=sum(LC),Elc=sum(Elc),Ooc=sum(OC),Eoc=sum(Eoc),PY=sum(PY))|>
 mutate(EARlc=(Olc-Elc)/PY,LLlc=EARlc-1.96*sqrt(Olc)/PY,ULlc=EARlc+1.96*sqrt(Olc)/PY,
        EARoc=(Ooc-Eoc)/PY,LLoc=EARoc-1.96*sqrt(Ooc)/PY,ULoc=EARoc+1.96*sqrt(Ooc)/PY) 

(Dl=Dt|>select(t,EARlc,LLlc,ULlc)|>mutate(cause="Leukemic Cause")|>rename(EAR=EARlc,LL=LLlc,UL=ULlc))
(Do=Dt|>select(t,EARoc,LLoc,ULoc)|>mutate(cause="Other Cause")|>rename(EAR=EARoc,LL=LLoc,UL=ULoc))
D=bind_rows(Do,Dl)

leg=theme(legend.margin=margin(0,0,0,0),legend.title=element_blank(),
          legend.position="top") #redfine leg without size increase
EARbrks=c(0,0.01,0.05,0.1)
ccE=coord_cartesian(ylim=c(0,0.08))
tc=function(sz) theme_classic(base_size=sz)
ghp01=geom_hline(yintercept=c(0.01),col="gray")
gh0=geom_hline(yintercept=0)
gx=xlab("Years Since CML Diagnosis")
gyE=ylab("Excess Absolute Risk of Death")
geE=geom_errorbar(aes(ymin=LL,ymax=UL),width=0.2)#for absolute risks
gl=geom_line()
D=D|>mutate(t=t+0.5)
D|>ggplot(aes(x=t,y=EAR,col = cause))+ghp01+geE+geom_point(size=0.7)+gyE +gl+ccE+
  scale_y_continuous(minor_breaks=NULL,breaks=EARbrks) +tc(12)+gh0+leg+gx+
  scale_color_manual(values = c("gray","black"),name = "cause")
ggsave(file=paste0("LExtra/outs/7A_OCtimeD2.pdf"),height=3,width=3.2) 

# testCOD7.R   
graphics.off();rm(list=ls())#clear plots and environment 
library(tidyverse)  
load("~/data/CMLepi/D7e.RData") 
D=D|>group_by(t)|>summarize(Olc=sum(LC),Elc=sum(Elc),
                             Oin=sum(IN),Ein=sum(Ein),Odk=sum(DK),Edk=sum(Edk),
                             Ocv=sum(CV),Ecv=sum(Ecv),Oca=sum(CA),Eca=sum(Eca),
                        Oash=sum(ASH),Eash=sum(Eash),Oyoc=sum(YOC),Eyoc=sum(Eyoc),PY=sum(PY))|>
 mutate(EARlc=(Olc-Elc)/PY,LLlc=EARlc-1.96*sqrt(Olc)/PY,ULlc=EARlc+1.96*sqrt(Olc)/PY,
        EARin=(Oin-Ein)/PY,LLin=EARin-1.96*sqrt(Oin)/PY,ULin=EARin+1.96*sqrt(Oin)/PY, 
        EARdk=(Odk-Edk)/PY,LLdk=EARdk-1.96*sqrt(Odk)/PY,ULdk=EARdk+1.96*sqrt(Odk)/PY, 
        EARcv=(Ocv-Ecv)/PY,LLcv=EARcv-1.96*sqrt(Ocv)/PY,ULcv=EARcv+1.96*sqrt(Ocv)/PY, 
        EARca=(Oca-Eca)/PY,LLca=EARca-1.96*sqrt(Oca)/PY,ULca=EARca+1.96*sqrt(Oca)/PY, 
        EARash=(Oash-Eash)/PY,LLash=EARash-1.96*sqrt(Oash)/PY,ULash=EARash+1.96*sqrt(Oash)/PY, 
        EARyoc=(Oyoc-Eyoc)/PY,LLyoc=EARyoc-1.96*sqrt(Oyoc)/PY,ULyoc=EARyoc+1.96*sqrt(Oyoc)/PY) 


(Dlc=D|>select(t,EARlc,LLlc,ULlc)|>mutate(cause="LC")|>rename(EAR=EARlc,LL=LLlc,UL=ULlc))
(Dyoc=D|>select(t,EARyoc,LLyoc,ULyoc)|>mutate(cause="YOC")|>rename(EAR=EARyoc,LL=LLyoc,UL=ULyoc))
(Dash=D|>select(t,EARash,LLash,ULash)|>mutate(cause="ASH")|>rename(EAR=EARash,LL=LLash,UL=ULash))
(Dca=D|>select(t,EARca,LLca,ULca)|>mutate(cause="CA")|>rename(EAR=EARca,LL=LLca,UL=ULca))
(Dcv=D|>select(t,EARcv,LLcv,ULcv)|>mutate(cause="CV")|>rename(EAR=EARcv,LL=LLcv,UL=ULcv))
(Ddk=D|>select(t,EARdk,LLdk,ULdk)|>mutate(cause="DK")|>rename(EAR=EARdk,LL=LLdk,UL=ULdk))
(Din=D|>select(t,EARin,LLin,ULin)|>mutate(cause="IN")|>rename(EAR=EARin,LL=LLin,UL=ULin))
Dt=bind_rows(Dlc,Dcv,Ddk,Dash,Dca,Din,Dyoc)|>mutate(t=t+0.5)

myt=theme(legend.text=element_text(size=12),strip.text=element_text(size=12))
geEAR=geom_errorbar(aes(ymin=LL,ymax=UL),width=.2,col="gray")
EARbrks=c(0,0.01,0.05,0.1)
ghE=geom_hline(yintercept=c(0.01),col="gray")
leg=theme(legend.margin=margin(0,0,0,0),legend.title=element_blank(),
          legend.position="top") #redfine leg without size increase
gyE=ylab("Excess Absolute Risk of Death")
geE=geom_errorbar(aes(ymin=LL,ymax=UL),width=0.2)#for absolute risks
jco=ggsci::scale_color_jco()
tc=function(sz) theme_classic(base_size=sz)
gh0=geom_hline(yintercept=0)
gx=xlab("Years Since CML Diagnosis")
ccE=coord_cartesian(ylim=c(-0.003,0.012))
Dt|>filter(cause%in%c("CV","DK","ASH"))|>mutate(cause=as_factor(cause))|>
  ggplot(aes(x=t,y=EAR,col=cause,shape=cause))+ghE+geEAR+geom_point(size=2)+geom_line()+
  gyE+scale_y_continuous(minor_breaks=NULL,breaks=EARbrks)+ jco+tc(13)+gh0 + leg + gx +ccE
ggsave(file=paste0("LExtra/outs/7_Cv2.pdf"),height=3,width=3)

Dt|>filter(cause%in%c("CA","IN","YOC"))|>
  ggplot(aes(x=t,y=EAR,col=cause,shape=cause))+ghE+geEAR+geom_point(size=2)+geom_line()+
  gyE+scale_y_continuous(minor_breaks=NULL,breaks=EARbrks)+ jco+tc(13)+gh0 + leg + gx # +ccE
ggsave(file=paste0("LExtra/outs/7_Dv2.pdf"),height=3,width=3)

Dt|>filter(cause%in%c("LC","YOC"))|>
  ggplot(aes(x=t,y=EAR,col=cause,shape=cause))+ghE+geEAR+geom_point(size=2)+geom_line()+
  gyE+scale_y_continuous(minor_breaks=NULL,breaks=EARbrks)+ jco+tc(13)+gh0 + leg + gx
ggsave(file=paste0("LExtra/outs/7_LC-YOC_inverse_correlations.pdf"),height=3,width=3)

Dt|>filter(cause%in%c("LC","CA"))|>
  ggplot(aes(x=t,y=EAR,col=cause,shape=cause))+ghE+geEAR+geom_point(size=2)+geom_line()+
  gyE+scale_y_continuous(minor_breaks=NULL,breaks=EARbrks)+ jco+tc(13)+gh0 + leg + gx
ggsave(file=paste0("LExtra/outs/7_LC-CA__correlations.pdf"),height=3,width=3)


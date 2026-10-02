# 7A_time.R   
graphics.off();rm(list=ls())#clear plots and environment 
library(tidyverse)  
load("~/data/CMLepi/grp_KK.RData") 
D=D|>filter(year>=2005,age>=60,age<=89) 
names(D) 
#  [1] "age"   "year"  "t"     "sex"   "PY"    "ASH"   "BEN"   "CA"    "COPD"  "CV"    "DK"    "ILL"   "IN"    "LC"    "MCA"   "MSPC" 
# [17] "OCD"   "Eash"  "Eben"  "Eca"   "Ecopd" "Ecv"   "Edk"   "Eill"  "Ein"   "Elc"   "Emca"  "Emspc" "Eocd" 
D=D|>group_by(t)|>summarize(O=sum(O),E=sum(E),Oash=sum(ASH),Eash=sum(Eash),Oben=sum(BEN),Eben=sum(Eben),
                             Oca=sum(CA),Eca=sum(Eca),Ocopd=sum(COPD),Ecopd=sum(Ecopd),
                             Ocv=sum(CV),Ecv=sum(Ecv),Odk=sum(DK),Edk=sum(Edk),
                             Ocvdk=sum(CV)+sum(DK),Ecvdk=sum(Ecv)+sum(Edk),
                             Oin=sum(IN),Ein=sum(Ein),Oill=sum(ILL),Eill=sum(Eill),
                             Olc=sum(LC),Elc=sum(Elc), 
                             Omca=sum(MCA),Emca=sum(Emca),Omspc=sum(MSPC),Emspc=sum(Emspc),
                            Oocd=sum(OCD),Eocd=sum(Eocd),PY=sum(PY))|>  
 mutate(EAR=(O-E)/PY,LL=EAR-1.96*sqrt(O)/PY,UL=EAR+1.96*sqrt(O)/PY, 
        EARash=(Oash-Eash)/PY,LLash=EARash-1.96*sqrt(Oash)/PY,ULash=EARash+1.96*sqrt(Oash)/PY, 
        EARben=(Oben-Eben)/PY,LLben=EARben-1.96*sqrt(Oben)/PY,ULben=EARben+1.96*sqrt(Oben)/PY, 
        EARca=(Oca-Eca)/PY,LLca=EARca-1.96*sqrt(Oca)/PY,ULca=EARca+1.96*sqrt(Oca)/PY, 
        EARcopd=(Ocopd-Ecopd)/PY,LLcopd=EARcopd-1.96*sqrt(Ocopd)/PY,ULcopd=EARcopd+1.96*sqrt(Ocopd)/PY, 
        EARcv=(Ocv-Ecv)/PY,LLcv=EARcv-1.96*sqrt(Ocv)/PY,ULcv=EARcv+1.96*sqrt(Ocv)/PY, 
        EARdk=(Odk-Edk)/PY,LLdk=EARdk-1.96*sqrt(Odk)/PY,ULdk=EARdk+1.96*sqrt(Odk)/PY, 
        EARcvdk=(Ocvdk-Ecvdk)/PY,LLcvdk=EARcvdk-1.96*sqrt(Ocvdk)/PY,ULcvdk=EARcvdk+1.96*sqrt(Ocvdk)/PY, 
        EARill=(Oill-Eill)/PY,LLill=EARill-1.96*sqrt(Oill)/PY,ULill=EARill+1.96*sqrt(Oill)/PY, 
        EARin=(Oin-Ein)/PY,LLin=EARin-1.96*sqrt(Oin)/PY,ULin=EARin+1.96*sqrt(Oin)/PY, 
        EARlc=(Olc-Elc)/PY,LLlc=EARlc-1.96*sqrt(Olc)/PY,ULlc=EARlc+1.96*sqrt(Olc)/PY,
        EARmca=(Omca-Emca)/PY,LLmca=EARmca-1.96*sqrt(Omca)/PY,ULmca=EARmca+1.96*sqrt(Omca)/PY, 
        EARmspc=(Omspc-Emspc)/PY,LLmspc=EARmspc-1.96*sqrt(Omspc)/PY,ULmspc=EARmspc+1.96*sqrt(Omspc)/PY, 
        EARocd=(Oocd-Eocd)/PY,LLocd=EARocd-1.96*sqrt(Oocd)/PY,ULocd=EARocd+1.96*sqrt(Oocd)/PY) 

(Dac=D|>select(t,EAR,LL,UL)|>mutate(cause="AC")) # all cause totals using HMD data (not a fit) for E
(Dash=D|>select(t,EARash,LLash,ULash)|>mutate(cause="ASH")|>rename(EAR=EARash,LL=LLash,UL=ULash))
(Dben=D|>select(t,EARben,LLben,ULben)|>mutate(cause="BEN")|>rename(EAR=EARben,LL=LLben,UL=ULben))
(Dca=D|>select(t,EARca,LLca,ULca)|>mutate(cause="CA")|>rename(EAR=EARca,LL=LLca,UL=ULca))
(Dcopd=D|>select(t,EARcopd,LLcopd,ULcopd)|>mutate(cause="COPD")|>rename(EAR=EARcopd,LL=LLcopd,UL=ULcopd))
(Dcv=D|>select(t,EARcv,LLcv,ULcv)|>mutate(cause="CV")|>rename(EAR=EARcv,LL=LLcv,UL=ULcv))
(Ddk=D|>select(t,EARdk,LLdk,ULdk)|>mutate(cause="DK")|>rename(EAR=EARdk,LL=LLdk,UL=ULdk))
(Dcvdk=D|>select(t,EARcvdk,LLcvdk,ULcvdk)|>mutate(cause="CVDK")|>rename(EAR=EARcvdk,LL=LLcvdk,UL=ULcvdk))
(Dill=D|>select(t,EARill,LLill,ULill)|>mutate(cause="ILL")|>rename(EAR=EARill,LL=LLill,UL=ULill))
(Din=D|>select(t,EARin,LLin,ULin)|>mutate(cause="IN")|>rename(EAR=EARin,LL=LLin,UL=ULin))
(Dlc=D|>select(t,EARlc,LLlc,ULlc)|>mutate(cause="LC")|>rename(EAR=EARlc,LL=LLlc,UL=ULlc))
(Dmca=D|>select(t,EARmca,LLmca,ULmca)|>mutate(cause="MCA")|>rename(EAR=EARmca,LL=LLmca,UL=ULmca))
(Dmspc=D|>select(t,EARmspc,LLmspc,ULmspc)|>mutate(cause="MSPC")|>rename(EAR=EARmspc,LL=LLmspc,UL=ULmspc))
(Docd=D|>select(t,EARocd,LLocd,ULocd)|>mutate(cause="OCD")|>rename(EAR=EARocd,LL=LLocd,UL=ULocd))
Dt=bind_rows(Dac,Dash,Dben,Dca,Dcopd,Dcv,Ddk,Dcvdk,Dlc,Dill,Din,Dmca,Dmspc,Docd)|>mutate(t=t+0.5) #in survRate, grouping uses left ends 

myt=theme(legend.text=element_text(size=12),strip.text=element_text(size=12))
geEAR=geom_errorbar(aes(ymin=LL,ymax=UL),width=.05,col="gray")
EARbrks=c(0,0.01,0.05,0.1)
ghE=geom_hline(yintercept=c(0.01),col="gray")
leg=theme(legend.margin=margin(0,0,0,0),legend.title=element_blank(),
          legend.position="top") #redfine leg without size increase
gyE=ylab("Excess Absolute Risk of Death")
geE=geom_errorbar(aes(ymin=LL,ymax=UL),width=0.2)#for absolute risks
jco=ggsci::scale_color_jco()
shps=scale_shape_manual(values = c("circle", "triangle", "square","diamond","square cross"))
tc=function(sz) theme_classic(base_size=sz)
gh0=geom_hline(yintercept=0)
gx=xlab("Time Since Diagnosis in Years")
sbb=theme(strip.background=element_blank(),strip.text=element_text(size=10)) 

codA=c("LC","CVDK","OCD","ASH")
codA=c("LC","CVDK","ASH")
codA=c("AC","LC","CVDK","OCD")
Dt|>filter(cause%in%codA)|>mutate(cause=factor(cause,levels=codA))|>
  ggplot(aes(x=t,y=EAR,col=cause,shape=cause))+gh0 +ghE+geEAR+geom_line()+geom_point(size=1)+ 
  gyE+scale_y_continuous(minor_breaks=NULL,breaks=EARbrks)+ 
  scale_x_continuous(breaks=c(0,10,20))+ 
  jco+tc(13)+ leg + gx
ggsave(file="LE/outs/7A_time.pdf",height=3.5,width=3.5) # LC + TKI driven excess risks

sbb=theme(strip.background=element_blank(),strip.text=element_text(size=10)) 
codA=c("LC","CVDK","ASH","CV","DK")
codA=c("LC","CVDK","OCD","CV","DK","AC")
codA=c("CVDK","AC")
Dt|>filter(!cause%in%codA)|>mutate(cause=factor(cause))|>
  ggplot(aes(x=t,y=EAR))+gh0 +ghE+geEAR+geom_line()+geom_point(size=1)+ 
  gyE+scale_y_continuous(minor_breaks=NULL,breaks=EARbrks)+ 
  scale_x_continuous(breaks=c(0,10,20))+ 
  jco+tc(13)+ leg + gx +facet_wrap(~cause,nrow=2)+sbb #+coord_cartesian(ylim=c(-0.005,0.05))
ggsave(file="LE/outs/7A_timePanel.pdf",height=4,width=8) 

d=D|>select(t,EARocd,EARcvdk,EAR)
d=d|>mutate(ratio=EARcvdk/EAR,t=t+0.5)  # 40%
d=d|>mutate(ratio=(EARcvdk+EARocd)/EAR,t=t+0.5) #>50% 
ghE=geom_hline(yintercept=c(0.01),col="gray")
gv=geom_vline(xintercept=c(15),col="gray")
d|>ggplot(aes(x=t,y=ratio))+gh0+ghE+gv+geom_line()+geom_point(size=0.5)+ gyE +tc(12)
# average over times >2 years, TKI account for roughly 25% of total EAR.  
# focus on times >10 years in 7B_time.R to get average ratio at steady state


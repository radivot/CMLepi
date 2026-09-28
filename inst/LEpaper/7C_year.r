# 7C_year.R   
# graphics.off();rm(list=ls())#clear plots and environment 
library(tidyverse)  
load("~/data/CMLepi/grp_KK.RData") 
D=D|>filter(age>=60,age<=89,t>=2,year>=2005) #no diff in EE vs KK since t>=5
D=D|>group_by(year)|>summarize(Oash=sum(ASH),Eash=sum(Eash),Oben=sum(BEN),Eben=sum(Eben),
                             Oca=sum(CA),Eca=sum(Eca),Ocopd=sum(COPD),Ecopd=sum(Ecopd),
                             Ocv=sum(CV),Ecv=sum(Ecv),Odk=sum(DK),Edk=sum(Edk),
                             Ocvdk=sum(CV)+sum(DK),Ecvdk=sum(Ecv)+sum(Edk),
                             Oin=sum(IN),Ein=sum(Ein),Oill=sum(ILL),Eill=sum(Eill),
                             Olc=sum(LC),Elc=sum(Elc), 
                             Omca=sum(MCA),Emca=sum(Emca),Omspc=sum(MSPC),Emspc=sum(Emspc),
                            Oocd=sum(OCD),Eocd=sum(Eocd),PY=sum(PY))|>  #have to stick UNK somewhere 
 mutate(EARash=(Oash-Eash)/PY,LLash=EARash-1.96*sqrt(Oash)/PY,ULash=EARash+1.96*sqrt(Oash)/PY, 
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


(Dash=D|>select(year,EARash,LLash,ULash)|>mutate(cause="ASH")|>rename(EAR=EARash,LL=LLash,UL=ULash))
(Dben=D|>select(year,EARben,LLben,ULben)|>mutate(cause="BEN")|>rename(EAR=EARben,LL=LLben,UL=ULben))
(Dca=D|>select(year,EARca,LLca,ULca)|>mutate(cause="CA")|>rename(EAR=EARca,LL=LLca,UL=ULca))
(Dcopd=D|>select(year,EARcopd,LLcopd,ULcopd)|>mutate(cause="COPD")|>rename(EAR=EARcopd,LL=LLcopd,UL=ULcopd))
(Dcv=D|>select(year,EARcv,LLcv,ULcv)|>mutate(cause="CV")|>rename(EAR=EARcv,LL=LLcv,UL=ULcv))
(Ddk=D|>select(year,EARdk,LLdk,ULdk)|>mutate(cause="DK")|>rename(EAR=EARdk,LL=LLdk,UL=ULdk))
(Dcvdk=D|>select(year,EARcvdk,LLcvdk,ULcvdk)|>mutate(cause="CVDK")|>rename(EAR=EARcvdk,LL=LLcvdk,UL=ULcvdk))
(Dill=D|>select(year,EARill,LLill,ULill)|>mutate(cause="ILL")|>rename(EAR=EARill,LL=LLill,UL=ULill))
(Din=D|>select(year,EARin,LLin,ULin)|>mutate(cause="IN")|>rename(EAR=EARin,LL=LLin,UL=ULin))
(Dlc=D|>select(year,EARlc,LLlc,ULlc)|>mutate(cause="LC")|>rename(EAR=EARlc,LL=LLlc,UL=ULlc))
(Dmca=D|>select(year,EARmca,LLmca,ULmca)|>mutate(cause="MCA")|>rename(EAR=EARmca,LL=LLmca,UL=ULmca))
(Dmspc=D|>select(year,EARmspc,LLmspc,ULmspc)|>mutate(cause="MSPC")|>rename(EAR=EARmspc,LL=LLmspc,UL=ULmspc))
(Docd=D|>select(year,EARocd,LLocd,ULocd)|>mutate(cause="OCD")|>rename(EAR=EARocd,LL=LLocd,UL=ULocd))
Dy=bind_rows(Dash,Dben,Dca,Dcopd,Dcv,Ddk,Dcvdk,Dlc,Dill,Din,Dmca,Dmspc,Docd)|>mutate(year=year+0.5) #in survRate, grouping uses left ends 

myt=theme(legend.text=element_text(size=12),strip.text=element_text(size=12))
geEAR=geom_errorbar(aes(ymin=LL,ymax=UL),width=.2,col="gray")
EARbrks=c(0,0.01,0.05,0.1)
ghE=geom_hline(yintercept=c(0.01,.1),col="gray")
leg=theme(legend.margin=margin(0,0,0,0),legend.title=element_blank(),
          legend.position="top") #redfine leg without size increase
gyE=ylab("Excess Absolute Risk of Death")
geE=geom_errorbar(aes(ymin=LL,ymax=UL),width=0.05)#for absolute risks
jco=ggsci::scale_color_jco()
shps=scale_shape_manual(values = c("circle", "triangle", "square","diamond","square cross"))
tc=function(sz) theme_classic(base_size=sz)
gh0=geom_hline(yintercept=0)
gv90=geom_vline(xintercept=90,col="gray50")
gx=xlab("Attained Year")
ccE=coord_cartesian(ylim=c(-0.005,0.06),xlim=c(2005,2025))
codA=c("LC","CVDK","OCD")
Dy|>filter(cause%in%codA)|>mutate(cause=factor(cause,levels=codA))|>
  ggplot(aes(x=year,y=EAR,col=cause,shape=cause))+gh0 +geEAR+geom_line()+geom_point(size=1)+shps+
  gyE+scale_y_continuous(minor_breaks=NULL,breaks=EARbrks)+ jco+tc(13)+ leg + gx + ccE
ggsave(file="LE/outs/7C_year.pdf",height=3,width=3.5) 

sbb=theme(strip.background=element_blank(),strip.text=element_text(size=10)) 
codA=c("LC","CVDK","OCD","CV","DK")
Dy|>filter(!cause%in%codA)|>mutate(cause=factor(cause))|>
  ggplot(aes(x=year,y=EAR))+gh0 +geEAR+geom_line()+geom_point(size=1)+ gyE+
  scale_y_continuous(minor_breaks=NULL,breaks=EARbrks)+
  scale_x_continuous(minor_breaks=NULL,breaks=c(2010,2015,2020,2025))+
  jco+tc(13)+ leg + gx +facet_wrap(~cause,ncol=4)+sbb +coord_cartesian(ylim=c(-0.005,0.01))+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave(file="LE/outs/7C_yearPanel.pdf",height=3,width=8) 



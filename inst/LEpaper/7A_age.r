# 7A_age.R   
# graphics.off();rm(list=ls())#clear plots and environment 
library(tidyverse)  
# load("~/data/CMLepi/grp_KKold.RData") # no diffs till 95 if t>2 
load("~/data/CMLepi/grp_KK.RData")
D=D|>filter(t>=2,year>=2005,age<90) #no diff in EE vs KK since t>=5
names(D) 
#  [1] "age"   "year"  "t"     "sex"   "PY"    "ASH"   "BEN"   "CA"    "COPD"  "CV"    "DK"    "ILL"   "IN"    "LC"    "LIV"   "MCA"  
# [17] "MSPC"  "OCD"   "UNK"   "Eash"  "Eben"  "Eca"   "Ecopd" "Ecv"   "Edk"   "Eill"  "Ein"   "Elc"   "Eliv"  "Emca"  "Emspc" "Eocd" 
D=D|>group_by(age)|>summarize(Oash=sum(ASH),Eash=sum(Eash),Oben=sum(BEN),Eben=sum(Eben),
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


(Dash=D|>select(age,EARash,LLash,ULash)|>mutate(cause="ASH")|>rename(EAR=EARash,LL=LLash,UL=ULash))
(Dben=D|>select(age,EARben,LLben,ULben)|>mutate(cause="BEN")|>rename(EAR=EARben,LL=LLben,UL=ULben))
(Dca=D|>select(age,EARca,LLca,ULca)|>mutate(cause="CA")|>rename(EAR=EARca,LL=LLca,UL=ULca))
(Dcopd=D|>select(age,EARcopd,LLcopd,ULcopd)|>mutate(cause="COPD")|>rename(EAR=EARcopd,LL=LLcopd,UL=ULcopd))
(Dcv=D|>select(age,EARcv,LLcv,ULcv)|>mutate(cause="CV")|>rename(EAR=EARcv,LL=LLcv,UL=ULcv))
(Ddk=D|>select(age,EARdk,LLdk,ULdk)|>mutate(cause="DK")|>rename(EAR=EARdk,LL=LLdk,UL=ULdk))
(Dcvdk=D|>select(age,EARcvdk,LLcvdk,ULcvdk)|>mutate(cause="CVDK")|>rename(EAR=EARcvdk,LL=LLcvdk,UL=ULcvdk))
(Dill=D|>select(age,EARill,LLill,ULill)|>mutate(cause="ILL")|>rename(EAR=EARill,LL=LLill,UL=ULill))
(Din=D|>select(age,EARin,LLin,ULin)|>mutate(cause="IN")|>rename(EAR=EARin,LL=LLin,UL=ULin))
(Dlc=D|>select(age,EARlc,LLlc,ULlc)|>mutate(cause="LC")|>rename(EAR=EARlc,LL=LLlc,UL=ULlc))
(Dmca=D|>select(age,EARmca,LLmca,ULmca)|>mutate(cause="MCA")|>rename(EAR=EARmca,LL=LLmca,UL=ULmca))
(Dmspc=D|>select(age,EARmspc,LLmspc,ULmspc)|>mutate(cause="MSPC")|>rename(EAR=EARmspc,LL=LLmspc,UL=ULmspc))
(Docd=D|>select(age,EARocd,LLocd,ULocd)|>mutate(cause="OCD")|>rename(EAR=EARocd,LL=LLocd,UL=ULocd))
Da=bind_rows(Dash,Dben,Dca,Dcopd,Dcv,Ddk,Dcvdk,Dlc,Dill,Din,Dmca,Dmspc,Docd)|>mutate(age=age+0.5) #in survRate, grouping uses left ends 

myt=theme(legend.text=element_text(size=12),strip.text=element_text(size=12))
geEAR=geom_errorbar(aes(ymin=LL,ymax=UL),width=.05,col="gray")
EARbrks=c(0,0.01,0.05,0.1)
ghE=geom_hline(yintercept=c(0.01,.1),col="gray")
leg=theme(legend.margin=margin(0,0,0,0),legend.title=element_blank(),
          legend.position="top") #redfine leg without size increase
gyE=ylab("Excess Absolute Risk of Death")
jco=ggsci::scale_color_jco()
shps=scale_shape_manual(values = c("circle", "triangle", "square","diamond","square cross"))
tc=function(sz) theme_classic(base_size=sz)
gh0=geom_hline(yintercept=0)
gv69=geom_vline(xintercept=c(60),col="gray")
gx=xlab("Attained Age")
ccE=coord_cartesian(ylim=c(-0.002,0.1))
codA=c("LC","CVDK","OCD")
Da|>filter(cause%in%codA,age<95)|>mutate(cause=factor(cause,levels=codA))|>
  ggplot(aes(x=age,y=EAR,col=cause,shape=cause))+gh0+ghE+geEAR+geom_line()+geom_point(size=0.5)+shps+
  scale_x_continuous(minor_breaks=NULL,breaks=seq(20,90,10))+ 
  gyE+scale_y_continuous(minor_breaks=NULL,breaks=EARbrks)+ jco+tc(13)+ leg + gx  + ccE
ggsave(file="LE/outs/7A_age.pdf",height=3,width=3.5) # TKI driven excess risks rise with aging

sbb=theme(strip.background=element_blank(),strip.text=element_text(size=10)) 
ghp01=geom_hline(yintercept=c(0.01),col="gray")
codA=c("LC","CVDK","OCD","CV","DK")
Da|>filter(!cause%in%codA,age<95)|>mutate(cause=factor(cause))|>
  ggplot(aes(x=age,y=EAR))+gh0+ghp01 +geEAR+geom_line()+geom_point(size=0.5)+ gyE+
  scale_y_continuous(minor_breaks=NULL,breaks=EARbrks)+
  scale_x_continuous(minor_breaks=NULL,breaks=seq(30,90,20))+ 
  jco+tc(13)+ leg + gx +facet_wrap(~cause,ncol=4)+sbb +coord_cartesian(ylim=c(-0.007,0.014))
  # jco+tc(13)+ leg + gx +facet_wrap(~cause,ncol=4,scales="free_y")+sbb 
ggsave(file="LE/outs/7A_agePanel.pdf",height=3,width=8) 


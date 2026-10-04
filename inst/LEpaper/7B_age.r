# 7B_age.R   
# graphics.off();rm(list=ls())#clear plots and environment 
library(tidyverse)  
load("~/data/CMLepi/grp_KK.RData")
D=D|>filter(t>=10,year>=2005,age<90) 
names(D) 
#  [1] "age"   "year"  "t"     "sex"   "PY"    "ASH"   "BEN"   "CA"    "COPD"  "CV"    "DK"    "ILL"   "IN"    "LC"    "MCA"   "MSPC" 
# [17] "OCD"   "Eash"  "Eben"  "Eca"   "Ecopd" "Ecv"   "Edk"   "Eill"  "Ein"   "Elc"   "Emca"  "Emspc" "Eocd"  "m"     "O"     "E"    
D=D|>group_by(age)|>summarize(O=sum(O),E=sum(E),Oash=sum(ASH),Eash=sum(Eash),Oben=sum(BEN),Eben=sum(Eben),
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


(Dac=D|>select(age,EAR,LL,UL)|>mutate(cause="AC")) # all cause totals using HMD data (not a fit) for E
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
Da=bind_rows(Dac,Dash,Dben,Dca,Dcopd,Dcv,Ddk,Dcvdk,Dlc,Dill,Din,Dmca,Dmspc,Docd)|>mutate(age=age+0.5) #in survRate, grouping uses left ends 

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
ccE=coord_cartesian(ylim=c(-0.002,0.10))
codA=c("AC","LC","CVDK","OCD")
Da|>filter(cause%in%codA,age<95)|>mutate(cause=factor(cause,levels=codA))|>
  ggplot(aes(x=age,y=EAR,col=cause,shape=cause))+gh0+ghE+geEAR+geom_line()+geom_point(size=0.5)+shps+
  scale_x_continuous(minor_breaks=NULL,breaks=seq(20,90,10))+ 
  gyE+scale_y_continuous(minor_breaks=NULL,breaks=EARbrks)+ jco+tc(13)+ leg + gx  + ccE
ggsave(file="LE/outs/7B_age.pdf",height=3.5,width=3.5) # TKI driven excess risks rise with aging

sbb=theme(strip.background=element_blank(),strip.text=element_text(size=10)) 
ghp01=geom_hline(yintercept=c(0.01),col="gray")
codA=c("LC","CVDK","OCD","CV","DK","AC")
codA=c("CVDK","AC")
Da|>filter(!cause%in%codA,age<95)|>mutate(cause=factor(cause))|>
  ggplot(aes(x=age,y=EAR))+gh0+ghp01 +geEAR+geom_line()+geom_point(size=0.5)+ gyE+
  scale_y_continuous(minor_breaks=NULL,breaks=EARbrks)+
  scale_x_continuous(minor_breaks=NULL,breaks=seq(30,90,20))+ 
  jco+tc(13)+ leg + gx +facet_wrap(~cause,ncol=6)+sbb +coord_cartesian(ylim=c(-0.007,0.07))
ggsave(file="LE/outs/7B_agePanel.pdf",height=4,width=8) 
Da|>filter(age<60,cause=="CVDK")|>mutate(ageG=cut(age,seq(20,60,10)))|>
  group_by(ageG)|>summarize(EAR=mean(EAR))
#   ageG         EAR
#   <fct>      <dbl>
# 1 (20,30] 0.000398
# 2 (30,40] 0.000622
# 3 (40,50] 0.000717
# 4 (50,60] 0.00104 
Da|>filter(cause=="CVDK")|>mutate(ageG=cut(age,seq(20,90,10)))|>
  group_by(ageG)|>summarize(EAR=mean(EAR))#2.7%
Da|>filter(cause=="OCD")|>mutate(ageG=cut(age,seq(20,90,10)))|>
  group_by(ageG)|>summarize(EAR=mean(EAR)) #0.5%
Da|>filter(cause=="ASH")|>mutate(ageG=cut(age,seq(20,90,10)))|>
  group_by(ageG)|>summarize(EAR=mean(EAR)) #0.4% 
# so they sum to 3.8%, which is roughly half of 6.9% below
Da|>filter(cause=="AC")|>mutate(ageG=cut(age,seq(20,90,10)))|>
  group_by(ageG)|>summarize(EAR=mean(EAR))#7%
# and higher than 3.5% via LC
Da|>filter(cause=="LC")|>mutate(ageG=cut(age,seq(20,90,10)))|>
  group_by(ageG)|>summarize(EAR=mean(EAR))#3.4%


# 
# d=D|>filter(age>60)|>select(age,EAR,EARcvdk,EARocd)
# d=d|>mutate(ratio=(EARcvdk+EARocd)/EAR) # 35% if including ocd
# d=d|>mutate(ratio=EARcvdk/EAR) # 25% if just cvdk
# gh1=geom_hline(yintercept=c(0.25,0.35),col="gray")
# d|>ggplot(aes(x=age,y=ratio))+gh0+gh1+geom_line()+geom_point(size=0.5)+ gyE+tc(12)+
#   scale_y_continuous(breaks=c(0,seq(0.05,0.35,0.1)))
# # average over times >2 years, TKI account for roughly 25% to 35% of total EAR.  
# # focus on times >10 years in 7B_time.R to get average ratio at steady state

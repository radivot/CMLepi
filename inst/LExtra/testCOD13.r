# testCOD13.R   
graphics.off();rm(list=ls())#clear plots and environment 
library(tidyverse)  
load("~/data/CMLepi/D13e.RData") 
names(D) 
#  [1] "age"   "year"  "t"     "sex"   "PY"    "ASH"   "BEN"   "CA"    "COPD"  "CV"    "DK"    "ILL"   "IN"    "LC"    "LIV"   "MCA"  
# [17] "MSPC"  "OCD"   "UNK"   "Eash"  "Eben"  "Eca"   "Ecopd" "Ecv"   "Edk"   "Eill"  "Ein"   "Elc"   "Eliv"  "Emca"  "Emspc" "Eocd" 
D=D|>group_by(t)|>summarize(Oash=sum(ASH),Eash=sum(Eash),Oben=sum(BEN),Eben=sum(Eben),
                             Oca=sum(CA),Eca=sum(Eca),Ocopd=sum(COPD),Ecopd=sum(Ecopd),
                             Ocv=sum(CV),Ecv=sum(Ecv),Odk=sum(DK),Edk=sum(Edk),
                             Oin=sum(IN),Ein=sum(Ein),Oill=sum(ILL),Eill=sum(Eill),
                             Olc=sum(LC),Elc=sum(Elc),Oliv=sum(LIV),Eliv=sum(Eliv),
                             Omca=sum(MCA),Emca=sum(Emca),Omspc=sum(MSPC),Emspc=sum(Emspc),
                             # Oocd=sum(OCD),Eocd=sum(Eocd),PY=sum(PY))|>
                            Oocd=sum(OCD+UNK),Eocd=sum(Eocd),PY=sum(PY))|>  #bump up just O with UNK => fig ends in OCD+.pdf
 mutate(EARash=(Oash-Eash)/PY,LLash=EARash-1.96*sqrt(Oash)/PY,ULash=EARash+1.96*sqrt(Oash)/PY, 
        EARben=(Oben-Eben)/PY,LLben=EARben-1.96*sqrt(Oben)/PY,ULben=EARben+1.96*sqrt(Oben)/PY, 
        EARca=(Oca-Eca)/PY,LLca=EARca-1.96*sqrt(Oca)/PY,ULca=EARca+1.96*sqrt(Oca)/PY, 
        EARcopd=(Ocopd-Ecopd)/PY,LLcopd=EARcopd-1.96*sqrt(Ocopd)/PY,ULcopd=EARcopd+1.96*sqrt(Ocopd)/PY, 
        EARcv=(Ocv-Ecv)/PY,LLcv=EARcv-1.96*sqrt(Ocv)/PY,ULcv=EARcv+1.96*sqrt(Ocv)/PY, 
        EARdk=(Odk-Edk)/PY,LLdk=EARdk-1.96*sqrt(Odk)/PY,ULdk=EARdk+1.96*sqrt(Odk)/PY, 
        EARill=(Oill-Eill)/PY,LLill=EARill-1.96*sqrt(Oill)/PY,ULill=EARill+1.96*sqrt(Oill)/PY, 
        EARin=(Oin-Ein)/PY,LLin=EARin-1.96*sqrt(Oin)/PY,ULin=EARin+1.96*sqrt(Oin)/PY, 
        EARlc=(Olc-Elc)/PY,LLlc=EARlc-1.96*sqrt(Olc)/PY,ULlc=EARlc+1.96*sqrt(Olc)/PY,
        EARliv=(Oliv-Eliv)/PY,LLliv=EARliv-1.96*sqrt(Oliv)/PY,ULliv=EARliv+1.96*sqrt(Oliv)/PY,
        EARmca=(Omca-Emca)/PY,LLmca=EARmca-1.96*sqrt(Omca)/PY,ULmca=EARmca+1.96*sqrt(Omca)/PY, 
        EARmspc=(Omspc-Emspc)/PY,LLmspc=EARmspc-1.96*sqrt(Omspc)/PY,ULmspc=EARmspc+1.96*sqrt(Omspc)/PY, 
        EARocd=(Oocd-Eocd)/PY,LLocd=EARocd-1.96*sqrt(Oocd)/PY,ULocd=EARocd+1.96*sqrt(Oocd)/PY) 


(Dash=D|>select(t,EARash,LLash,ULash)|>mutate(cause="ASH")|>rename(EAR=EARash,LL=LLash,UL=ULash))
(Dben=D|>select(t,EARben,LLben,ULben)|>mutate(cause="BEN")|>rename(EAR=EARben,LL=LLben,UL=ULben))
(Dca=D|>select(t,EARca,LLca,ULca)|>mutate(cause="CA")|>rename(EAR=EARca,LL=LLca,UL=ULca))
(Dcopd=D|>select(t,EARcopd,LLcopd,ULcopd)|>mutate(cause="COPD")|>rename(EAR=EARcopd,LL=LLcopd,UL=ULcopd))
(Dcv=D|>select(t,EARcv,LLcv,ULcv)|>mutate(cause="CV")|>rename(EAR=EARcv,LL=LLcv,UL=ULcv))
(Ddk=D|>select(t,EARdk,LLdk,ULdk)|>mutate(cause="DK")|>rename(EAR=EARdk,LL=LLdk,UL=ULdk))
(Dill=D|>select(t,EARill,LLill,ULill)|>mutate(cause="ILL")|>rename(EAR=EARill,LL=LLill,UL=ULill))
(Din=D|>select(t,EARin,LLin,ULin)|>mutate(cause="IN")|>rename(EAR=EARin,LL=LLin,UL=ULin))
(Dlc=D|>select(t,EARlc,LLlc,ULlc)|>mutate(cause="LC")|>rename(EAR=EARlc,LL=LLlc,UL=ULlc))
(Dliv=D|>select(t,EARliv,LLliv,ULliv)|>mutate(cause="LIV")|>rename(EAR=EARliv,LL=LLliv,UL=ULliv))
(Dmca=D|>select(t,EARmca,LLmca,ULmca)|>mutate(cause="MCA")|>rename(EAR=EARmca,LL=LLmca,UL=ULmca))
(Dmspc=D|>select(t,EARmspc,LLmspc,ULmspc)|>mutate(cause="MSPC")|>rename(EAR=EARmspc,LL=LLmspc,UL=ULmspc))
(Docd=D|>select(t,EARocd,LLocd,ULocd)|>mutate(cause="OCD")|>rename(EAR=EARocd,LL=LLocd,UL=ULocd))
Dt=bind_rows(Dash,Dben,Dca,Dcopd,Dcv,Ddk,Dlc,Dliv,Dill,Din,Dmca,Dmspc,Docd)|>mutate(t=t+0.5)

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
ccE=coord_cartesian(ylim=c(-0.003,0.02))
Dt|>filter(cause%in%c("CV","DK","ASH"))|>mutate(cause=factor(cause,levels=c("CV","DK","ASH")))|>
  ggplot(aes(x=t,y=EAR,col=cause,shape=cause))+ghE+geEAR+geom_point(size=2)+geom_line()+
  gyE+scale_y_continuous(minor_breaks=NULL,breaks=EARbrks)+ jco+tc(13)+gh0 + leg + gx +ccE
ggsave(file=paste0("LExtra/outs/7C_tki.pdf"),height=3,width=3) # TKI driven excess risks rise with aging
# ASH is TKI cost driven, so its excess risk should fall with calendar time and thus aging, as seen here 

Dt|>filter(cause%in%c("LIV","MSPC"))|>
  ggplot(aes(x=t,y=EAR,col=cause,shape=cause))+ghE+geEAR+geom_point(size=2)+geom_line()+
  gyE+scale_y_continuous(minor_breaks=NULL,breaks=EARbrks)+ jco+tc(13)+gh0 + leg + gx # +ccE
ggsave(file=paste0("LExtra/outs/7D_better.pdf"),height=3,width=3)
# surviving CML selects for those who drink less alcohol, so LIV down with time. 
#   91   2.58   #MSPC    Alzheimers            so MSPC is mostly alzheimers ... not clear why it would trend down
#  100   0.238  #MSPC   stomach ulcers
#  103   0.025  #MSPC   Complications of birth
#  104   0.568  #MSPC   congenital conditions
#  105   0.77   #MSPC   perinatl condiditions   

Dt|>filter(cause%in%c("CA","IN"))|>
  ggplot(aes(x=t,y=EAR,col=cause,shape=cause))+ghE+geEAR+geom_point(size=2)+geom_line()+
  gyE+scale_y_continuous(minor_breaks=NULL,breaks=EARbrks)+ jco+tc(13)+gh0 + leg + gx # +ccE
ggsave(file=paste0("LExtra/outs/7D_CA_IN.pdf"),height=3,width=3)
## excess in first 5 years is via coincidental. IN staying up longer could be CML weakening immunity against infections

Dt|>filter(cause%in%c("MCA","ILL"))|>mutate(cause=as_factor(cause))|>
  ggplot(aes(x=t,y=EAR,col=cause,shape=cause))+ghE+geEAR+geom_point(size=2)+geom_line()+
  gyE+scale_y_continuous(minor_breaks=NULL,breaks=EARbrks)+ jco+tc(13)+gh0 + leg + gx
ggsave(file=paste0("LExtra/outs/7_zoomOut.pdf"),height=3,width=3)
## mspc = miscellaneous specific cancers, so not expected to include CML, and does not
## mca = miscellaneous malignant cancers is a grab bag, so include zoomed out CMLs
## ILL = ill-defined cancers may also include zoomed out CMLs

ccE=coord_cartesian(ylim=c(-0.003,0.01))
Dt|>filter(cause%in%c("BEN","OCD"))|>
  ggplot(aes(x=t,y=EAR,col=cause,shape=cause))+ghE+geEAR+geom_point(size=2)+geom_line()+
  gyE+scale_y_continuous(minor_breaks=NULL,breaks=EARbrks)+ jco+tc(13)+gh0 + leg + gx +ccE
ggsave(file=paste0("LExtra/outs/7_BEN-OCD+.pdf"),height=3,width=3) # obs includes 252 = no DC
# ggsave(file=paste0("LExtra/outs/7_BEN-OCD.pdf"),height=3,width=3) # obs is strictly OCD
## BEN looks like coincidental CML brought about via benign tumors (e.g. meningioma) that can kill
## OCD increasing with time/aging could be CML deaths caused by restricting TKI to IM at old ages

Dt|>filter(cause%in%c("LC"))|>
  ggplot(aes(x=t,y=EAR,col=cause,shape=cause))+ghE+geEAR+geom_point(size=2)+geom_line()+
  gyE+scale_y_continuous(minor_breaks=NULL,breaks=EARbrks)+ jco+tc(13)+gh0 + leg + gx
ggsave(file=paste0("LExtra/outs/7_LC.pdf"),height=3,width=3)


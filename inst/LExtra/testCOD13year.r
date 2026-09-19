# testCOD13year.R   
graphics.off();rm(list=ls())#clear plots and environment 
library(tidyverse)  

load("~/data/CMLepi/D13e.RData") 
names(D) 
#  [1] "age"   "year"  "t"     "sex"   "PY"    "ASH"   "BEN"   "CA"    "COPD"  "CV"    "DK"    "ILL"   "IN"    "LC"    "LIV"   "MCA"  
# [17] "MSPC"  "OCD"   "UNK"   "Eash"  "Eben"  "Eca"   "Ecopd" "Ecv"   "Edk"   "Eill"  "Ein"   "Elc"   "Eliv"  "Emca"  "Emspc" "Eocd" 
D=D|>filter(age>=80,t>=5,year>2005)
D=D|>group_by(year)|>summarize(Oash=sum(ASH),Eash=sum(Eash),Oben=sum(BEN),Eben=sum(Eben),
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


(Dash=D|>select(year,EARash,LLash,ULash)|>mutate(cause="ASH")|>rename(EAR=EARash,LL=LLash,UL=ULash))
(Dben=D|>select(year,EARben,LLben,ULben)|>mutate(cause="BEN")|>rename(EAR=EARben,LL=LLben,UL=ULben))
(Dca=D|>select(year,EARca,LLca,ULca)|>mutate(cause="CA")|>rename(EAR=EARca,LL=LLca,UL=ULca))
(Dcopd=D|>select(year,EARcopd,LLcopd,ULcopd)|>mutate(cause="COPD")|>rename(EAR=EARcopd,LL=LLcopd,UL=ULcopd))
(Dcv=D|>select(year,EARcv,LLcv,ULcv)|>mutate(cause="CV")|>rename(EAR=EARcv,LL=LLcv,UL=ULcv))
(Ddk=D|>select(year,EARdk,LLdk,ULdk)|>mutate(cause="DK")|>rename(EAR=EARdk,LL=LLdk,UL=ULdk))
(Dill=D|>select(year,EARill,LLill,ULill)|>mutate(cause="ILL")|>rename(EAR=EARill,LL=LLill,UL=ULill))
(Din=D|>select(year,EARin,LLin,ULin)|>mutate(cause="IN")|>rename(EAR=EARin,LL=LLin,UL=ULin))
(Dlc=D|>select(year,EARlc,LLlc,ULlc)|>mutate(cause="LC")|>rename(EAR=EARlc,LL=LLlc,UL=ULlc))
(Dliv=D|>select(year,EARliv,LLliv,ULliv)|>mutate(cause="LIV")|>rename(EAR=EARliv,LL=LLliv,UL=ULliv))
(Dmca=D|>select(year,EARmca,LLmca,ULmca)|>mutate(cause="MCA")|>rename(EAR=EARmca,LL=LLmca,UL=ULmca))
(Dmspc=D|>select(year,EARmspc,LLmspc,ULmspc)|>mutate(cause="MSPC")|>rename(EAR=EARmspc,LL=LLmspc,UL=ULmspc))
(Docd=D|>select(year,EARocd,LLocd,ULocd)|>mutate(cause="OCD")|>rename(EAR=EARocd,LL=LLocd,UL=ULocd))
Dy=bind_rows(Dash,Dben,Dca,Dcopd,Dcv,Ddk,Dlc,Dliv,Dill,Din,Dmca,Dmspc,Docd)|>mutate(year=year+0.5)

myt=theme(legend.text=element_text(size=12),strip.text=element_text(size=12))
geEAR=geom_errorbar(aes(ymin=LL,ymax=UL),width=.2,col="gray")
EARbrks=c(0,0.01,0.05,0.1)
ghE=geom_hline(yintercept=c(0.01,.1),col="gray")
leg=theme(legend.margin=margin(0,0,0,0),legend.title=element_blank(),
          legend.position="top") #redfine leg without size increase
gyE=ylab("Excess Absolute Risk of Death")
geE=geom_errorbar(aes(ymin=LL,ymax=UL),width=0.2)#for absolute risks
jco=ggsci::scale_color_jco()
shps=scale_shape_manual(values = c("circle", "triangle", "square","diamond","square cross"))
tc=function(sz) theme_classic(base_size=sz)
gh0=geom_hline(yintercept=0)
gv90=geom_vline(xintercept=90,col="gray50")
gx=xlab("Year")
ccE=coord_cartesian(ylim=c(-0.02,0.11))

codA=c("LC","CV","OCD")
Dy|>filter(cause%in%codA)|>mutate(cause=factor(cause,levels=codA))|>
  ggplot(aes(x=year,y=EAR,col=cause,shape=cause))+gh0 +ghE+geEAR+geom_line()+geom_point(size=2)+shps+
  # scale_x_continuous(minor_breaks=NULL,breaks=seq(20,90,10))+ 
  gyE+scale_y_continuous(minor_breaks=NULL,breaks=EARbrks)+ jco+tc(13)+ leg + gx + ccE

ggsave(file=paste0("LExtra/outs/7E_LCtkiOCDyear.pdf"),height=3,width=5) # TKI driven excess risks rise with aging
# ASH is TKI cost driven, so its excess risk should fall with calendar time and thus aging, as seen here 

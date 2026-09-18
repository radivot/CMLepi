# mkCOD13.R   
graphics.off();rm(list=ls())#clear plots and environment 
library(survival)
library(tidyverse)  
options(pillar.sigfig = 5) # Shows 5 significant digits to get decimal after year in tibble prints
load("~/data/CMLepi/cml.RData") #made in 1_SEERdata.R  53254
(d=d|>filter(yrdx>=2000,agedx<90,surv<80,surv>0)) #43932
d=d|>mutate(COD13=as.factor(COD13)) 
levels(d$COD13) # "alive" "ASH"   "BEN"   "CA"    "COPD"  "CV"    "DK"    "ILL"   "IN"    "LC"    "LIV"   "MCA"   "MSPC"  "OCD"   "UNK" 
d=d|>mutate(agedx=agedx+0.5,yrdx=yrdx+0.5) 
(d=d|>mutate(adx=agedx,astart=adx,astop=adx+surv,.before=COD))
(Da=survSplit(Surv(astart,astop,COD13)~.,d,cut = 1:107)|>tibble()|>relocate(astart:COD13,.before=COD)) 
(Da=Da|>mutate(ystart = yrdx + astart - adx, ystop  = yrdx + astop - adx,.before=COD))
(Day=survSplit(Surv(ystart,ystop,COD13)~.,Da,cut=2001:2022)|>tibble()|>relocate(ystart:COD13,.before=COD)) 
(Day=Day|>mutate(astart = adx + ystart - yrdx, astop  = adx + ystop - yrdx)) # fix problems
(Day=Day|>mutate(tstart = astart-adx, tstop  = astop-adx,.before=status)) ## and bring in tstart and tstop
(Dayt=survSplit(Surv(tstart,tstop,COD13)~.,Day,cut=1:22,episode="Time")|>tibble()|>relocate(tstart:COD13,.before=COD)) 
(Dayt=Dayt|>mutate(t=Time-1,age=floor(astart),year=floor(ystart),PY=tstop-tstart,.before=COD)|>select(-Time))
(Dayt=Dayt|>mutate(ASH=ifelse(COD13=="ASH",1,0),.before=COD))
(Dayt=Dayt|>mutate(BEN=ifelse(COD13=="BEN",1,0),.before=COD))
(Dayt=Dayt|>mutate(CA=ifelse(COD13=="CA",1,0),.before=COD))
(Dayt=Dayt|>mutate(COPD=ifelse(COD13=="COPD",1,0),.before=COD))
(Dayt=Dayt|>mutate(CV=ifelse(COD13=="CV",1,0),.before=COD))
(Dayt=Dayt|>mutate(DK=ifelse(COD13=="DK",1,0),.before=COD))
(Dayt=Dayt|>mutate(ILL=ifelse(COD13=="ILL",1,0),.before=COD))
(Dayt=Dayt|>mutate(IN=ifelse(COD13=="IN",1,0),.before=COD))
(Dayt=Dayt|>mutate(LC=ifelse(COD13=="LC",1,0),.before=COD))
(Dayt=Dayt|>mutate(LIV=ifelse(COD13=="LIV",1,0),.before=COD))
(Dayt=Dayt|>mutate(MCA=ifelse(COD13=="MCA",1,0),.before=COD))
(Dayt=Dayt|>mutate(MSPC=ifelse(COD13=="MSPC",1,0),.before=COD))
(Dayt=Dayt|>mutate(OCD=ifelse(COD13=="OCD",1,0),.before=COD))
(Dayt=Dayt|>mutate(UNK=ifelse(COD13=="UNK",1,0),.before=COD))
system.time(dASH<-biostat3::survRate(Surv(PY,ASH)~age+year+t+sex, data=Dayt)|>tibble())#27s   
system.time(dBEN<-biostat3::survRate(Surv(PY,BEN)~age+year+t+sex, data=Dayt)|>tibble())#27s   
system.time(dCA<-biostat3::survRate(Surv(PY,CA)~age+year+t+sex, data=Dayt)|>tibble())#25s   
system.time(dCOPD<-biostat3::survRate(Surv(PY,COPD)~age+year+t+sex, data=Dayt)|>tibble())#25s   
system.time(dCV<-biostat3::survRate(Surv(PY,CV)~age+year+t+sex, data=Dayt)|>tibble())#25s   
system.time(dDK<-biostat3::survRate(Surv(PY,DK)~age+year+t+sex, data=Dayt)|>tibble())#25s   
system.time(dILL<-biostat3::survRate(Surv(PY,ILL)~age+year+t+sex, data=Dayt)|>tibble())#25s   
system.time(dIN<-biostat3::survRate(Surv(PY,IN)~age+year+t+sex, data=Dayt)|>tibble())#25s   
system.time(dLC<-biostat3::survRate(Surv(PY,LC)~age+year+t+sex, data=Dayt)|>tibble())#25s   
system.time(dLIV<-biostat3::survRate(Surv(PY,LIV)~age+year+t+sex, data=Dayt)|>tibble())#25s   
system.time(dMCA<-biostat3::survRate(Surv(PY,MCA)~age+year+t+sex, data=Dayt)|>tibble())#24s   
system.time(dMSPC<-biostat3::survRate(Surv(PY,MSPC)~age+year+t+sex, data=Dayt)|>tibble())#24s   
system.time(dOCD<-biostat3::survRate(Surv(PY,OCD)~age+year+t+sex, data=Dayt)|>tibble())#24s   
system.time(dUNK<-biostat3::survRate(Surv(PY,UNK)~age+year+t+sex, data=Dayt)|>tibble())#24s   

(dASH=dASH|>rename(PY=tstop,ASH=event)|>select(age:ASH))
(dBEN=dBEN|>rename(PY=tstop,BEN=event)|>select(age:BEN))
(dCA=dCA|>rename(PY=tstop,CA=event)|>select(age:CA))
(dCOPD=dCOPD|>rename(PY=tstop,COPD=event)|>select(age:COPD))
(dCV=dCV|>rename(PY=tstop,CV=event)|>select(age:CV))
(dDK=dDK|>rename(PY=tstop,DK=event)|>select(age:DK))
(dILL=dILL|>rename(PY=tstop,ILL=event)|>select(age:ILL))
(dIN=dIN|>rename(PY=tstop,IN=event)|>select(age:IN))
(dLC=dLC|>rename(PY=tstop,LC=event)|>select(age:LC))
(dLIV=dLIV|>rename(PY=tstop,LIV=event)|>select(age:LIV))
(dMCA=dMCA|>rename(PY=tstop,MCA=event)|>select(age:MCA))
(dMSPC=dMSPC|>rename(PY=tstop,MSPC=event)|>select(age:MSPC))
(dOCD=dOCD|>rename(PY=tstop,OCD=event)|>select(age:OCD))
(dUNK=dUNK|>rename(PY=tstop,UNK=event)|>select(age:UNK))

(D13=left_join(dASH,dBEN))
(D13=left_join(D13,dCA))
(D13=left_join(D13,dCOPD))
(D13=left_join(D13,dCV))
(D13=left_join(D13,dDK))
(D13=left_join(D13,dILL))
(D13=left_join(D13,dIN))
(D13=left_join(D13,dLC))
(D13=left_join(D13,dLIV))
(D13=left_join(D13,dMCA))
(D13=left_join(D13,dMSPC))
(D13=left_join(D13,dOCD))
(D13=left_join(D13,dUNK))
save(D13,file="~/data/CMLepi/D13.RData") #179 kb file 

(D=D13|>filter(age>20)) # G fits in mkG12 are only for ages over 20
# CMLepi::mkG13() # makes file read in below
load("~/data/CMLepi/G13.RData") 
(nms=names(G)) # "ASH"   "BEN"   "CA"    "COPD"  "CV"    "DK"    "ILL"   "IN"    "LC"    "MMC"   "MSC"   "OCD"  
D=D|>mutate(num=ASH,denom=PY)
D$Eash=as.numeric(exp(predict(G[["ASH"]],D)))
D=D|>mutate(num=BEN)
D$Eben=as.numeric(exp(predict(G[["BEN"]],D)))
D=D|>mutate(num=CA)
D$Eca=as.numeric(exp(predict(G[["CA"]],D)))
D=D|>mutate(num=COPD)
D$Ecopd=as.numeric(exp(predict(G[["COPD"]],D)))
D=D|>mutate(num=CV)
D$Ecv=as.numeric(exp(predict(G[["CV"]],D)))
D=D|>mutate(num=DK)
D$Edk=as.numeric(exp(predict(G[["DK"]],D)))
D=D|>mutate(num=ILL)
D$Eill=as.numeric(exp(predict(G[["ILL"]],D)))
D=D|>mutate(num=IN)
D$Ein=as.numeric(exp(predict(G[["IN"]],D)))
D=D|>mutate(num=LC)
D$Elc=as.numeric(exp(predict(G[["LC"]],D)))
D=D|>mutate(num=LIV)
D$Eliv=as.numeric(exp(predict(G[["LIV"]],D)))
D=D|>mutate(num=MCA)
D$Emca=as.numeric(exp(predict(G[["MCA"]],D)))
D=D|>mutate(num=MSPC)
D$Emspc=as.numeric(exp(predict(G[["MSPC"]],D)))
D=D|>mutate(num=OCD)
D$Eocd=as.numeric(exp(predict(G[["OCD"]],D)))
# note that there are no expected values for UNK since it isn't in the US Mortality CODs 
D=D|>select(-denom,-num)
save(D,file="~/data/CMLepi/D13e.RData") #3.6 MB file 


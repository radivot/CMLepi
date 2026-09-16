# mkCOD7.R   
graphics.off();rm(list=ls())#clear plots and environment 
library(survival)
library(tidyverse)  
options(pillar.sigfig = 5) # Shows 5 significant digits to get decimal after year in tibble prints
load("~/data/CMLepi/cml.RData") #made in mkSEER.R  53254
(d=d|>filter(yrdx>=2000,agedx<90,surv<80,surv>0)) #43932
d=d|>mutate(COD7=as.factor(COD7)) 
levels(d$COD7) #"alive" "ASH"   "CA"    "CV"    "DK"    "IN"    "LC"    "YOC" 
d=d|>mutate(agedx=agedx+0.5,yrdx=yrdx+0.5) 
(d=d|>mutate(adx=agedx,astart=adx,astop=adx+surv,.before=COD))
(Da=survSplit(Surv(astart,astop,COD7)~.,d,cut = 1:107)|>tibble()|>relocate(astart:COD7,.before=COD)) 
(Da=Da|>mutate(ystart = yrdx + astart - adx, ystop  = yrdx + astop - adx,.before=COD))
(Day=survSplit(Surv(ystart,ystop,COD7)~.,Da,cut=2000:2023)|>tibble()|>relocate(ystart:COD7,.before=COD)) 
(Day=Day|>mutate(astart = adx + ystart - yrdx, astop  = adx + ystop - yrdx)) # fix problems
(Day=Day|>mutate(tstart = astart-adx, tstop  = astop-adx,.before=status)) ## and bring in tstart and tstop
(Dayt=survSplit(Surv(tstart,tstop,COD7)~.,Day,cut=1:45,episode="Time")|>tibble()|>relocate(tstart:COD7,.before=COD)) 
(Dayt=Dayt|>mutate(t=Time-1,age=floor(astart),year=floor(ystart),PY=tstop-tstart,.before=COD)|>select(-Time))
(Dayt=Dayt|>mutate(ASH=ifelse(COD7=="ASH",1,0),.before=COD))
(Dayt=Dayt|>mutate(CA=ifelse(COD7=="CA",1,0),.before=COD))
(Dayt=Dayt|>mutate(CV=ifelse(COD7=="CV",1,0),.before=COD))
(Dayt=Dayt|>mutate(DK=ifelse(COD7=="DK",1,0),.before=COD))
(Dayt=Dayt|>mutate(IN=ifelse(COD7=="IN",1,0),.before=COD))
(Dayt=Dayt|>mutate(LC=ifelse(COD7=="LC",1,0),.before=COD))
(Dayt=Dayt|>mutate(YOC=ifelse(COD7=="YOC",1,0),.before=COD))
system.time(dASH<-biostat3::survRate(Surv(PY,ASH)~age+year+t+sex, data=Dayt)|>tibble())#27s   
system.time(dCA<-biostat3::survRate(Surv(PY,CA)~age+year+t+sex, data=Dayt)|>tibble())#25s   
system.time(dCV<-biostat3::survRate(Surv(PY,CV)~age+year+t+sex, data=Dayt)|>tibble())#25s   
system.time(dDK<-biostat3::survRate(Surv(PY,DK)~age+year+t+sex, data=Dayt)|>tibble())#25s   
system.time(dIN<-biostat3::survRate(Surv(PY,IN)~age+year+t+sex, data=Dayt)|>tibble())#25s   
system.time(dLC<-biostat3::survRate(Surv(PY,LC)~age+year+t+sex, data=Dayt)|>tibble())#25s   
system.time(dYOC<-biostat3::survRate(Surv(PY,YOC)~age+year+t+sex, data=Dayt)|>tibble())#24s   

(dASH=dASH|>rename(PY=tstop,ASH=event)|>select(age:ASH))
(dCA=dCA|>rename(PY=tstop,CA=event)|>select(age:CA))
(dCV=dCV|>rename(PY=tstop,CV=event)|>select(age:CV))
(dDK=dDK|>rename(PY=tstop,DK=event)|>select(age:DK))
(dIN=dIN|>rename(PY=tstop,IN=event)|>select(age:IN))
(dLC=dLC|>rename(PY=tstop,LC=event)|>select(age:LC))
(dYOC=dYOC|>rename(PY=tstop,YOC=event)|>select(age:YOC))

(D7=left_join(dYOC,dLC))
(D7=left_join(D7,dIN))
(D7=left_join(D7,dDK))
(D7=left_join(D7,dCV))
(D7=left_join(D7,dCA))
(D7=left_join(D7,dASH))
save(D7,file="~/data/CMLepi/D7.RData") #179 kb file 

(D=D7|>filter(age>20)) # G fits in mkG6 are only for ages over 20
load("~/data/CMLepi/G6.RData") # made by mkMorts.R
(nms=names(G)) # "LC"  "CA"  "IN"  "CV"  "DK"  "ASH" "YOC"
D=D|>mutate(num=YOC,denom=PY)
D$Eyoc=as.numeric(exp(predict(G[["YOC"]],D)))
D=D|>mutate(num=ASH)
D$Eash=as.numeric(exp(predict(G[["ASH"]],D)))
D=D|>mutate(num=DK)
D$Edk=as.numeric(exp(predict(G[["DK"]],D)))
D=D|>mutate(num=CV)
D$Ecv=as.numeric(exp(predict(G[["CV"]],D)))
D=D|>mutate(num=IN)
D$Ein=as.numeric(exp(predict(G[["IN"]],D)))
D=D|>mutate(num=CA)
D$Eca=as.numeric(exp(predict(G[["CA"]],D)))
D=D|>mutate(num=LC)
D$Elc=as.numeric(exp(predict(G[["LC"]],D)))
D=D|>select(-denom,-num)
save(D,file="~/data/CMLepi/D7e.RData") #2.2 MB file 
D



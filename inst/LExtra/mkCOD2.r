# mkCOD2.R   #extend Dt to LC and OC causes of deaths
graphics.off();rm(list=ls())#clear plots and environment 
library(survival)
library(tidyverse)  
options(pillar.sigfig = 5) # Shows 5 significant digits to get decimal after year in tibble prints
load("~/data/CMLepi/cml.RData") #made in mkSEER.R  53254
(d=d|>filter(yrdx>=2000,agedx<90,surv<80,surv>0)) #43932
d=d|>mutate(COD2=as.factor(COD2)) 
str(d$COD2[1:10])
d=d|>mutate(agedx=agedx+0.5,yrdx=yrdx+0.5) 
(d=d|>mutate(adx=agedx,astart=adx,astop=adx+surv,.before=COD))
(Da=survSplit(Surv(astart,astop,COD2)~.,d,cut = 1:107)|>tibble()|>relocate(astart:COD2,.before=COD)) 
(Da=Da|>mutate(ystart = yrdx + astart - adx, ystop  = yrdx + astop - adx,.before=COD))
(Day=survSplit(Surv(ystart,ystop,COD2)~.,Da,cut=2000:2023)|>tibble()|>relocate(ystart:COD2,.before=COD)) 
(Day=Day|>mutate(astart = adx + ystart - yrdx, astop  = adx + ystop - yrdx)) # fix problems
(Day=Day|>mutate(tstart = astart-adx, tstop  = astop-adx,.before=status)) ## and bring in tstart and tstop
(Dayt=survSplit(Surv(tstart,tstop,COD2)~.,Day,cut=1:45,episode="Time")|>tibble()|>relocate(tstart:COD2,.before=COD)) 
(Dayt=Dayt|>mutate(t=Time-1,age=floor(astart),year=floor(ystart),PY=tstop-tstart,.before=COD)|>select(-Time))
(Dayt=Dayt|>mutate(LC=ifelse(COD2=="LC",1,0),.before=COD))
(Dayt=Dayt|>mutate(OC=ifelse(COD2=="OC",1,0),.before=COD))
system.time(dLC<-biostat3::survRate(Surv(PY,LC)~age+year+t+sex, data=Dayt)|>tibble())#25s   
system.time(dOC<-biostat3::survRate(Surv(PY,OC)~age+year+t+sex, data=Dayt)|>tibble())#24s   
load("~/data/mrt/us_mort.RData")
(m<-us_mort |>filter(Sex != "Total", Year > 1999)|>select(year=Year,age=Age,sex=Sex,m=Mortality) ) 
(dLC=dLC|>rename(PY=tstop,LC=event)|>select(age:LC))
(dOC=dOC|>rename(PY=tstop,OC=event)|>select(age:OC))
(D2=left_join(dLC,dOC))
(D2=left_join(D2,m))
(D2=D2|>mutate(Elc=0,Eoc=m*PY)) #approximate LC mort as 0 and OC mort as AC
save(D2,file="~/data/CMLepi/D2.RData") #201 kb file 


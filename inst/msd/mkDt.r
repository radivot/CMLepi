# mkDt.R   #formula version of mkD0.R (exclude 90+, S=0 and S=NA and focus on >=2000) that includes time since Dx
graphics.off();rm(list=ls())#clear plots and environment 
library(survival)
library(tidyverse)  
options(pillar.sigfig = 5) # Shows 5 significant digits to get decimal after year in tibble prints
load("~/data/CMLepi/cml.RData") #made in mkSEER.R  53254
(d=d|>filter(yrdx>=2000,agedx<90,surv<80,surv>0)) #43932
d=d|>mutate(agedx=agedx+0.5,yrdx=yrdx+0.5) 
(d=d|>mutate(adx=agedx,astart=adx,astop=adx+surv,.before=COD))
# (Da=d|>survSplit(cut = 1:107, event = "status",start = "astart", end = "astop")|>tibble()|>relocate(astart:status,.before=COD)) 
(Da=survSplit(Surv(astart,astop,status)~.,d,cut = 1:107)|>tibble()|>relocate(astart:status,.before=COD)) 
(Da=Da|>mutate(ystart = yrdx + astart - adx, ystop  = yrdx + astop - adx,.before=COD))
# (Day=Da|>survSplit(cut=2000:2023,event="status",start="ystart",end="ystop")|>tibble()|>relocate(ystart:status,.before=COD)) 
(Day=survSplit(Surv(ystart,ystop,status)~.,Da,cut=2000:2023)|>tibble()|>relocate(ystart:status,.before=COD)) 
(Day=Day|>mutate(astart = adx + ystart - yrdx, astop  = adx + ystop - yrdx)) # fix problems
(Day=Day|>mutate(tstart = astart-adx, tstop  = astop-adx,.before=status)) ## and bring in tstart and tstop
# (Dayt=Day|>survSplit(cut=1:45,event="status",start="tstart",end="tstop",episode="Time")|>tibble()|>relocate(tstart:status,.before=COD)) 
(Dayt=survSplit(Surv(tstart,tstop,status)~.,Day,cut=1:45,episode="Time")|>tibble()|>relocate(tstart:status,.before=COD)) 
(Dayt=Dayt|>mutate(t=Time-1,age=floor(astart),year=floor(ystart),PY=tstop-tstart,.before=COD)|>select(-Time))
system.time(dayts <- biostat3::survRate(Surv(PY,status)~age+year+t+sex, data=Dayt)|>tibble())#24s   
load("~/data/mrt/us_mort.RData")
(m<-us_mort |>filter(Sex != "Total", Year > 1999)|>select(year=Year,age=Age,sex=Sex,m=Mortality) ) 
(Dt=dayts|>rename(PY=tstop,O=event)|>select(age:O))
(Dt=left_join(Dt,m))
save(Dt,file="~/data/CMLepi/Dt.RData") #201 kb file 


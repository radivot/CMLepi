# mkD0.R   #keep things as simple as possible by excluding 90+, S=0 and S=NA and focusing on >=2000
graphics.off();rm(list=ls())#clear plots and environment 
library(tidyverse)  
options(pillar.sigfig = 5) # Shows 5 significant digits to get decimal after year in tibble prints
load("~/data/CMLepi/cml.RData") #made in mkSEER.R  53254
(d=d|>filter(yrdx>=2000,agedx<90,surv<80,surv>0)) #43932
d=d|>mutate(agedx=agedx+0.5,yrdx=yrdx+0.5) #need to shift both to match preshift of just age in SEERaBomb matrix approach which ties them naturally 
(d=d|>mutate(adx=agedx,astart=adx,astop=adx+surv,.before=COD))
(Da=d|>survival::survSplit(cut = 1:107, event = "status",start = "astart", end = "astop")|>tibble()|>relocate(astart:status,.before=COD)) 
(Da=Da|>mutate(ystart = yrdx + astart - adx, ystop  = yrdx + astop - adx,.before=COD))
(Day=Da|>survival::survSplit(cut=2000:2023,event="status",start="ystart",end="ystop")|>tibble()|>relocate(ystart:status,.before=COD)) 
(Day=Day|>mutate(astart = adx + ystart - yrdx, astop  = adx + ystop - yrdx)) # fix problems
(Day=Day|>mutate(tstart = astart-adx, tstop  = astop-adx,.before=status)) ## and bring in tstart and tstop
(Dayt=Day|>survival::survSplit(cut=1:45,event="status",start="tstart",end="tstop",episode="Time")|>tibble()|>relocate(tstart:status,.before=COD)) 
(Dayt=Dayt|>mutate(t=Time-1,age=floor(astart),year=floor(ystart),PY=tstop-tstart,.before=COD)|>select(-Time))
system.time(days <- biostat3::survRate(survival::Surv(PY,status)~age+year+sex, data=Dayt)|>tibble())#2.5 secs   

load("~/data/mrt/us_mort.RData")
(m<-us_mort |>filter(Sex != "Total", Year > 1999)|>select(year=Year,age=Age,sex=Sex,m=Mortality) ) 
(D0=days|>rename(PY=tstop,O=event)|>select(age:O))
(D0=left_join(D0,m))
save(D0,file="~/data/CMLepi/D0.RData") #55 kb file 


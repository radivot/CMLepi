# mkNoShift.R does not add 0.5 years to agedx
graphics.off();rm(list=ls())#clear plots and environment 
library(biostat3)  # for survRate, and survSplit via survival 
library(tidyverse)  
options(pillar.sigfig = 5) # Shows 5 significant digits to get decimal after year in tibble prints
load("~/data/CMLepi/cml.RData") #made in mkSEER.R  53.2k
d=d|>mutate(agedx=ifelse(agedx==90,92.3,agedx)) # set over 90 to 92.3
d=d|>mutate(surv=ifelse(surv>80,0.01,surv)) #set NA surv to 0.01 (3.65 days)
d=d|>mutate(surv=ifelse(surv==0,0.01,surv)) #set  0 surv to 0.01 (else survSplit throws 'zero' parameter must be less than any observed times)
(d=d|>mutate(adx=agedx,astart=adx,astop=adx+surv,.before=COD)) #### was adx=agedx+0.5 before  ##############
(Da=d|>survSplit(cut = 1:107, event = "status",start = "astart", end = "astop")|>tibble()|>relocate(astart:status,.before=COD)) 
(Da=Da|>mutate(ystart = yrdx + astart - adx, ystop  = yrdx + astop - adx,.before=COD))
(Day=Da|>survSplit(cut=1975:2023,event="status",start="ystart",end="ystop")|>tibble()|>relocate(ystart:status,.before=COD)) 
(Day=Day|>mutate(astart = adx + ystart - yrdx, astop  = adx + ystop - yrdx)) # fix problems
(Day=Day|>mutate(tstart = astart-adx, tstop  = astop-adx,.before=status)) ## and bring in tstart and tstop
(Dayt=Day|>survSplit(cut=1:45,event="status",start="tstart",end="tstop",episode="Time")|>tibble()|>relocate(tstart:status,.before=COD)) 
(Dayt=Dayt|>mutate(t=Time-1,age=floor(astart),year=floor(ystart),PY=tstop-tstart,.before=COD)|>select(-Time))
(Dayt=Dayt|>filter(PY>1e-5))  #noise removed => back to size before, so we really didn't need the last split other than to create the time column t
system.time(dayts <- survRate(Surv(PY,status)~age+year+t+sex, data=Dayt)|>tibble()) #47 secs  89k
(day=day|>rename(PY=tstop,O=event)|>select(age:O))
load("~/data/mrt/us_mort.RData")
(m<-us_mort |>filter(Sex != "Total", Year > 1974)|>select(year=Year,age=Age,sex=Sex,m=Mortality) ) 
(dayts=dayts|>rename(PY=tstop,O=event)|>select(age:O))
(dayts=left_join(dayts,m))
save(dayts,file="~/data/CMLepi/daytsNoShift.RData") #311 kb file

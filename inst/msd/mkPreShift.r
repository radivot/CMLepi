# mkPreShift.R adds 0.5 years to agedx upfront
graphics.off();rm(list=ls())#clear plots and environment 
library(survival)  # for survRate, and survSplit via survival 
library(tidyverse)  
options(pillar.sigfig = 5) # Shows 5 significant digits to get decimal after year in tibble prints
load("~/data/CMLepi/cml.RData") #made in mkSEER.R  53.2k
d=d|>mutate(agedx=agedx+0.5,yrdx=yrdx+0.5) #need to shift both to match preshift of just age in SEERaBomb matrix approach which ties them naturally 
d=d|>mutate(agedx=ifelse(agedx>90,92.3,agedx)) # set over 90 to 92.3
d=d|>mutate(surv=ifelse(surv>80,0.01,surv)) #set NA surv to 0.01 (3.65 days)
d=d|>mutate(surv=ifelse(surv==0,0.01,surv)) #set  0 surv to 0.01 (else survSplit throws 'zero' parameter must be less than any observed times)
(d=d|>mutate(adx=agedx,astart=adx,astop=adx+surv,.before=COD)) #### was adx=agedx+0.5 before, now shit at top  ##############
(Da=d|>survSplit(cut = 1:110, event = "status",start = "astart", end = "astop")|>tibble()|>relocate(astart:status,.before=COD)) 
(Da=Da|>mutate(ystart = yrdx + astart - adx, ystop  = yrdx + astop - adx,.before=COD))
(Day=Da|>survSplit(cut=1975:2023,event="status",start="ystart",end="ystop")|>tibble()|>relocate(ystart:status,.before=COD)) 
(Day=Day|>mutate(astart = adx + ystart - yrdx, astop  = adx + ystop - yrdx)) # fix problems
load("~/data/mrt/us_mort.RData")
(m<-us_mort |>filter(Sex != "Total", Year > 1974)|>select(year=Year,age=Age,sex=Sex,m=Mortality) ) 
(Day=Day|>mutate(age=floor(astart),year=floor(ystart),PY=ystop-ystart,.before=COD))
system.time(days <- biostat3::survRate(Surv(PY,status)~age+year+sex, data=Day)|>tibble()) #4 secs 9k
(days=days|>rename(PY=tstop,O=event)|>select(age:O))
(days=left_join(days,m))
save(days,file="~/data/CMLepi/daysPreShift.RData") 



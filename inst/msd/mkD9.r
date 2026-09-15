# mk9.R   
graphics.off();rm(list=ls())#clear plots and environment 
library(tidyverse)  
options(pillar.sigfig = 5) # Shows 5 significant digits to get decimal after year in tibble prints
load("~/data/CMLepi/cml.RData") #53,254   made in mkSEER.R  
(d=d|>filter(yrdx>=2000)) #45,636 
d=d|>mutate(agedx=agedx+0.5,yrdx=yrdx+0.5) #need to shift both to match preshift of just age in SEERaBomb matrix approach which ties them naturally 
# d=d|>mutate(agedx=ifelse(agedx>90,92.3,agedx)) # set over 90 to 92.3  => ### No, this already creates diffs ###
d=d|>mutate(agedx=ifelse(agedx>90,92.5,agedx)) # try keeping it on the age grid  =>  #### yes, this eliminates the diffs!!! ######
(d=d|>filter(yrdx>=2000,surv<80,surv>0)) #exclude all others     
(d=d|>mutate(adx=agedx,astart=adx,astop=adx+surv,.before=COD))
Da=d|>survival::survSplit(cut = 1:107, event = "status",start = "astart", end = "astop")|>tibble()|>relocate(astart:status,.before=COD)
(Da=Da|>mutate(ystart = yrdx + astart - adx, ystop  = yrdx + astop - adx,.before=COD)) # 351,606 
(Day=Da|>survival::survSplit(cut=1975:2023,event="status",start="ystart",end="ystop")|>tibble()|>relocate(ystart:status,.before=COD)) 
(Day=Day|>mutate(astart = adx + ystart - yrdx, astop  = adx + ystop - yrdx)) # fix problems
(Day=Day|>mutate(tstart = astart-adx, tstop  = astop-adx,.before=status)) ## and bring in tstart and tstop
(Dayt=Day|>survival::survSplit(cut=1:45,event="status",start="tstart",end="tstop",episode="Time")|>tibble()|>relocate(tstart:status,.before=COD)) 
(Dayt=Dayt|>mutate(t=Time-1,age=floor(astart),year=floor(ystart),PY=tstop-tstart,.before=COD)|>select(-Time)) #637,263
system.time(dayts <- biostat3::survRate(survival::Surv(PY,status)~age+year+t+sex, data=Dayt)|>tibble())#25s   
load("~/data/mrt/us_mort.RData")
(m<-us_mort |>filter(Sex != "Total", Year > 1974)|>select(year=Year,age=Age,sex=Sex,m=Mortality) ) 
(D9=dayts|>rename(PY=tstop,O=event)|>select(age:O))
(D9=left_join(D9,m))
save(D9,file="~/data/CMLepi/D92.5.RData") 
save(D9,file="~/data/CMLepi/D92.3.RData")


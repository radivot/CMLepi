# mk5.R   
graphics.off();rm(list=ls())#clear plots and environment 
library(tidyverse)  
options(pillar.sigfig = 5) # Shows 5 significant digits to get decimal after year in tibble prints
load("~/data/CMLepi/cml.RData") #53,254   made in mkSEER.R  
(d=d|>filter(yrdx>=2000)) #45,636 
d=d|>mutate(agedx=agedx+0.5,yrdx=yrdx+0.5) #need to shift both to match preshift of just age in SEERaBomb matrix approach which ties them naturally 
d=d|>mutate(agedx=ifelse(agedx>90,92.5,agedx)) # based on mkD9.R, keep this on the matrix grid inside SEERaBomb::msd()
#if surv = NA (i.e. is >80), set it to 5 (no healthcare received), but first push back yrdx and agedx by 5 years 
d=d|>mutate(yrdx=ifelse(surv>80,yrdx-5,yrdx)) 
d=d|>mutate(agedx=ifelse(surv>80,agedx-5,agedx)) 
d|>filter(agedx<1) # one patient was 3.5 at Dx (now -1.5), so set this to 0.5 with a surv of 3 years
#        id sex    agedx   yrdx histo3 cancer      surv status   COD COD2  COD7  CODS                                        
# 1 21278145 Female   0.5 2008.5   9875 CML     0.046543      1   130 OC    CA    In situ, benign or unknown behavior neoplasm
# 2 52361446 Female  -1.5 1997.5   9863 CML    89.706         1    78 LC    LC    Chronic Myeloid Leukemia                    
# 3 52927570 Female   0.5 2009.5   9863 CML    13.435         0     0 alive alive Alive                                       
d=d|>mutate(surv=ifelse(agedx<0,3,surv)) #set this patient's survival time to 3 years  
d=d|>mutate(agedx=ifelse(agedx<0,0.5,agedx)) #and set negative age to 0.5 
d=d|>mutate(surv=ifelse(surv>80,5,surv)) # now set the rest of the NA survs to 5 years   
d=d|>mutate(surv=ifelse(surv==0,0.01,surv)) #set  0 surv to 0.01 (else survSplit throws 'zero' parameter must be less than any observed times)
(d=d|>mutate(adx=agedx,astart=adx,astop=adx+surv,.before=COD))  # A tibble: 45,636 × 15
Da=d|>survival::survSplit(cut = 1:107, event = "status",start = "astart", end = "astop")|>tibble()|>relocate(astart:status,.before=COD)
(Da=Da|>mutate(ystart = yrdx + astart - adx, ystop  = yrdx + astop - adx,.before=COD))
(Day=Da|>survival::survSplit(cut=1975:2023,event="status",start="ystart",end="ystop")|>tibble()|>relocate(ystart:status,.before=COD)) 
(Day=Day|>mutate(astart = adx + ystart - yrdx, astop  = adx + ystop - yrdx)) # fix problems
(Day=Day|>mutate(tstart = astart-adx, tstop  = astop-adx,.before=status)) ## and bring in tstart and tstop
(Dayt=Day|>survival::survSplit(cut=1:45,event="status",start="tstart",end="tstop",episode="Time")|>tibble()|>relocate(tstart:status,.before=COD)) 
(Dayt=Dayt|>mutate(t=Time-1,age=floor(astart),year=floor(ystart),PY=tstop-tstart,.before=COD)|>select(-Time))
system.time(dayts <- biostat3::survRate(survival::Surv(PY,status)~age+year+t+sex, data=Dayt)|>tibble())#25s   
load("~/data/mrt/us_mort.RData")
(m<-us_mort |>filter(Sex != "Total", Year > 1974)|>select(year=Year,age=Age,sex=Sex,m=Mortality) ) 
(D5=dayts|>rename(PY=tstop,O=event)|>select(age:O))
(D5=left_join(D5,m))
save(D5,file="~/data/CMLepi/D5.RData") #213 kb file 


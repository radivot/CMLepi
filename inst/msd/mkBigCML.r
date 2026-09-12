# mkBigCML   # extends biostat3/AgeYearTime.R to make a data file for all ages at Dx
graphics.off();rm(list=ls())#clear plots and environment 
library(tidyverse)  
options(pillar.sigfig = 5) # Shows 5 significant digits to get decimal after year in tibble prints
load("~/data/CMLepi/cml.RData") #made in mkSEER.R  53.2k
d=d|>mutate(agedx=agedx+0.5,yrdx=yrdx+0.5) #need to shift both to match preshift of just age in SEERaBomb matrix approach which ties them naturally 
d=d|>mutate(agedx=ifelse(agedx>90,92.3,agedx)) # set over 90 to 92.3
d=d|>mutate(surv=ifelse(surv>80,0.01,surv)) #set NA surv to 0.01 (3.65 days)
d=d|>mutate(surv=ifelse(surv==0,0.01,surv)) #set  0 surv to 0.01 (else survSplit throws 'zero' parameter must be less than any observed times)
min(d$surv) #0.00274 = 1 day
(d=d|>mutate(adx=agedx,astart=adx,astop=adx+surv,.before=COD))
# # A tibble: 53,254 × 15
#        id sex    agedx  yrdx histo3 cancer      surv status   adx astart  astop   COD COD2  COD7  CODS                                                
#     <int> <fct>  <dbl> <dbl>  <dbl> <chr>      <dbl>  <dbl> <dbl>  <dbl>  <dbl> <int> <chr> <chr> <fct>                                               
#  1   5253 Female    71  2002   9863 CML    13.996         1  71.5   71.5 85.496   208 OC    YOC   Other Cause of Death                                
#  2   9874 Male      33  2002   9875 CML    21.744         0  33.5   33.5 55.244     0 alive alive Alive                                               
max(d$astop) #107.6 years
(Da=d|>survival::survSplit(cut = 1:107, event = "status",start = "astart", end = "astop")|>tibble()|>relocate(astart:status,.before=COD)) 
(Da=Da|>mutate(ystart = yrdx + astart - adx, ystop  = yrdx + astop - adx,.before=COD))
# # A tibble: 400,612 × 17
#       id sex    agedx   yrdx histo3 cancer   surv   adx astart astop status ystart ystop   COD COD2  COD7  CODS                
#    <int> <fct>  <dbl>  <dbl>  <dbl> <chr>   <dbl> <dbl>  <dbl> <dbl>  <dbl>  <dbl> <dbl> <int> <chr> <chr> <fct>               
#  1  5253 Female  71.5 2002.5   9863 CML    13.996  71.5   71.5    72      0 2002.5  2003   208 OC    YOC   Other Cause of Death
#  2  5253 Female  71.5 2002.5   9863 CML    13.996  71.5   72      73      0 2003    2004   208 OC    YOC   Other Cause of Death
#  3  5253 Female  71.5 2002.5   9863 CML    13.996  71.5   73      74      0 2004    2005   208 OC    YOC   Other Cause of Death
(Day=Da|>survival::survSplit(cut=1975:2023,event="status",start="ystart",end="ystop")|>tibble()|>relocate(ystart:status,.before=COD)) 
# # A tibble: 401,966 × 17
#       id sex    agedx   yrdx histo3 cancer   surv   adx astart astop ystart ystop status   COD COD2  COD7  CODS                
#    <int> <fct>  <dbl>  <dbl>  <dbl> <chr>   <dbl> <dbl>  <dbl> <dbl>  <dbl> <dbl>  <dbl> <int> <chr> <chr> <fct>               
#  1  5253 Female  71.5 2002.5   9863 CML    13.996  71.5   71.5    72 2002.5  2003      0   208 OC    YOC   Other Cause of Death
#  2  5253 Female  71.5 2002.5   9863 CML    13.996  71.5   72      73 2003    2004      0   208 OC    YOC   Other Cause of Death
#  3  5253 Female  71.5 2002.5   9863 CML    13.996  71.5   73      74 2004    2005      0   208 OC    YOC   Other Cause of Death
(Day=Day|>mutate(astart = adx + ystart - yrdx, astop  = adx + ystop - yrdx)) # fix problems
(Day=Day|>mutate(tstart = astart-adx, tstop  = astop-adx,.before=status)) ## and bring in tstart and tstop
max(d$surv) #45.8 years
(Dayt=Day|>survival::survSplit(cut=1:45,event="status",start="tstart",end="tstop",episode="Time")|>tibble()|>relocate(tstart:status,.before=COD)) 
(Dayt=Dayt|>mutate(t=Time-1,age=floor(astart),year=floor(ystart),PY=tstop-tstart,.before=COD)|>select(-Time))
# # A tibble: 725,104 × 23
#       id sex    agedx   yrdx histo3 cancer   surv   adx astart astop ystart ystop tstart tstop status     t   age  year    PY   COD COD2  COD7  CODS  
#    <int> <fct>  <dbl>  <dbl>  <dbl> <chr>   <dbl> <dbl>  <dbl> <dbl>  <dbl> <dbl>  <dbl> <dbl>  <dbl> <dbl> <dbl> <dbl> <dbl> <int> <chr> <chr> <fct> 
#  1  5253 Female  71.5 2002.5   9863 CML    13.996  71.5   71.5    72 2002.5  2003    0     0.5      0     0    71  2002   0.5   208 OC    YOC   Other…
#  2  5253 Female  71.5 2002.5   9863 CML    13.996  71.5   72      73 2003    2004    0.5   1        0     0    72  2003   0.5   208 OC    YOC   Other…
#  3  5253 Female  71.5 2002.5   9863 CML    13.996  71.5   72      73 2003    2004    1     1.5      0     1    72  2003   0.5   208 OC    YOC   Other…
#  4  5253 Female  71.5 2002.5   9863 CML    13.996  71.5   73      74 2004    2005    1.5   2        0     1    73  2004   0.5   208 OC    YOC   Other…

# now take individual data to grouped data using biostat3 function survRate
system.time(day <- biostat3::survRate(survival::Surv(PY,status)~age+year, data=Dayt)|>tibble()) #2 secs  4.7k
system.time(days <- biostat3::survRate(survival::Surv(PY,status)~age+year+sex, data=Dayt)|>tibble()) #4 secs  8.8k
system.time(dayt <- biostat3::survRate(survival::Surv(PY,status)~age+year+t, data=Dayt)|>tibble()) #21 secs  55.2k
system.time(dayts <- biostat3::survRate(survival::Surv(PY,status)~age+year+t+sex, data=Dayt)|>tibble()) #49 secs  91.2k
system.time(daytsh <- biostat3::survRate(survival::Surv(PY,status)~age+year+t+sex+histo3, data=Dayt)|>tibble())#83 secs  129.9k   
# There is no reason to put agedx or yrdx on RHS of formula since they can be computed as age - t and year - t.
# WARNING: Using "~age+year+t+sex+histo3+yrdx+aged" with mem.maxVSize(vsize=Inf) after hitting 64GB lim swaps disk space for ever ... don't do it! 

#Now add m cols and save binaries
load("~/data/mrt/us_mort.RData")
(m<-us_mort |>filter(Sex != "Total", Year > 1974)|>select(year=Year,age=Age,sex=Sex,m=Mortality) ) 
(mt<-us_mort|>filter(Sex == "Total", Year > 1974)|>select(year=Year,age=Age,m=Mortality) ) 

(Day=day|>rename(PY=tstop,O=event)|>select(age:O))
(Day=left_join(Day,mt)) 
save(Day,file="~/data/CMLepi/Day.RData")  # 58kb file

(Days=days|>rename(PY=tstop,O=event)|>select(age:O)) 
(Days=left_join(Days,m))
save(Days,file="~/data/CMLepi/Days.RData")  # 79kb file

(Dayt=dayt|>rename(PY=tstop,O=event)|>select(age:O))
(Dayt=left_join(Dayt,mt))
save(Dayt,file="~/data/CMLepi/Dayt.RData") #220kb file

(Dayts=dayts|>rename(PY=tstop,O=event)|>select(age:O))
(Dayts=left_join(Dayts,m))
save(Dayts,file="~/data/CMLepi/Dayts.RData") #311 kb file

(Daytsh=daytsh|>rename(PY=tstop,O=event)|>select(age:O))
(Daytsh=left_join(Daytsh,m))
save(Daytsh,file="~/data/CMLepi/Daytsh.RData") #463 kb file 


# mkBigCML   # extends biostat3/AgeYearTime.R to make a data file for all ages at Dx
graphics.off();rm(list=ls())#clear plots and environment 
library(biostat3)  # for survRate, and survSplit via survival 
library(tidyverse)  
options(pillar.sigfig = 5) # Shows 5 significant digits to get decimal after year in tibble prints
load("~/data/CMLepi/cml.RData") #made in mkSEER.R  53.2k
d=d|>mutate(agedx=ifelse(agedx==90,92.3,agedx)) # set over 90 to 92.3
d=d|>mutate(surv=ifelse(surv>80,0.01,surv)) #set NA surv to 0.01 (3.65 days)
d=d|>mutate(surv=ifelse(surv==0,0.01,surv)) #set  0 surv to 0.01 (else survSplit throws 'zero' parameter must be less than any observed times)
min(d$surv) #0.00274 = 1 day
(d=d|>mutate(adx=agedx+0.5,astart=adx,astop=adx+surv,.before=COD))
# # A tibble: 53,254 × 15
#        id sex    agedx  yrdx histo3 cancer      surv status   adx astart  astop   COD COD2  COD7  CODS                                                
#     <int> <fct>  <dbl> <dbl>  <dbl> <chr>      <dbl>  <dbl> <dbl>  <dbl>  <dbl> <int> <chr> <chr> <fct>                                               
#  1   5253 Female    71  2002   9863 CML    13.996         1  71.5   71.5 85.496   208 OC    YOC   Other Cause of Death                                
#  2   9874 Male      33  2002   9875 CML    21.744         0  33.5   33.5 55.244     0 alive alive Alive                                               
max(d$astop) #108.1 years
(Da=d|>survSplit(cut = 1:107, event = "status",start = "astart", end = "astop")|>tibble()|>relocate(astart:status,.before=COD)) 
(Da=Da|>mutate(ystart = yrdx + astart - adx, ystop  = yrdx + astop - adx,.before=COD))
# # A tibble: 401,076 × 17
#       id sex    agedx  yrdx histo3 cancer   surv   adx astart astop status ystart  ystop   COD COD2  COD7  CODS                
#    <int> <fct>  <dbl> <dbl>  <dbl> <chr>   <dbl> <dbl>  <dbl> <dbl>  <dbl>  <dbl>  <dbl> <int> <chr> <chr> <fct>               
#  1  5253 Female    71  2002   9863 CML    13.996  71.5   71.5    72      0 2002   2002.5   208 OC    YOC   Other Cause of Death
#  2  5253 Female    71  2002   9863 CML    13.996  71.5   72      73      0 2002.5 2003.5   208 OC    YOC   Other Cause of Death
#  3  5253 Female    71  2002   9863 CML    13.996  71.5   73      74      0 2003.5 2004.5   208 OC    YOC   Other Cause of Death
#  4  5253 Female    71  2002   9863 CML    13.996  71.5   74      75      0 2004.5 2005.5   208 OC    YOC   Other Cause of Death
(Day=Da|>survSplit(cut=1975:2023,event="status",start="ystart",end="ystop")|>tibble()|>relocate(ystart:status,.before=COD)) 
# # A tibble: 724,214 × 17
#       id sex    agedx  yrdx histo3 cancer   surv   adx astart astop ystart  ystop status   COD COD2  COD7  CODS                
#    <int> <fct>  <dbl> <dbl>  <dbl> <chr>   <dbl> <dbl>  <dbl> <dbl>  <dbl>  <dbl>  <dbl> <int> <chr> <chr> <fct>               
#  1  5253 Female    71  2002   9863 CML    13.996  71.5   71.5    72 2002   2002.5      0   208 OC    YOC   Other Cause of Death
#  2  5253 Female    71  2002   9863 CML    13.996  71.5   72      73 2002.5 2003        0   208 OC    YOC   Other Cause of Death
#  3  5253 Female    71  2002   9863 CML    13.996  71.5   72      73 2003   2003.5      0   208 OC    YOC   Other Cause of Death
#  4  5253 Female    71  2002   9863 CML    13.996  71.5   73      74 2003.5 2004        0   208 OC    YOC   Other Cause of Death
#  5  5253 Female    71  2002   9863 CML    13.996  71.5   73      74 2004   2004.5      0   208 OC    YOC   Other Cause of Death
(Day=Day|>mutate(astart = adx + ystart - yrdx, astop  = adx + ystop - yrdx)) # fix problems
(Day=Day|>mutate(tstart = astart-adx, tstop  = astop-adx,.before=status)) ## and bring in tstart and tstop
max(d$surv) #45.8 years
(Dayt=Day|>survSplit(cut=1:45,event="status",start="tstart",end="tstop",episode="Time")|>tibble()|>relocate(tstart:status,.before=COD)) 
(Dayt=Dayt|>mutate(t=Time-1,age=floor(astart),year=floor(ystart),PY=tstop-tstart,.before=COD)|>select(-Time))
# # A tibble: 725,258 × 23
#       id sex    agedx  yrdx histo3 cancer   surv   adx astart astop ystart  ystop tstart tstop status     t   age  year    PY   COD COD2  COD7  CODS  
#    <int> <fct>  <dbl> <dbl>  <dbl> <chr>   <dbl> <dbl>  <dbl> <dbl>  <dbl>  <dbl>  <dbl> <dbl>  <dbl> <dbl> <dbl> <dbl> <dbl> <int> <chr> <chr> <fct> 
#  1  5253 Female    71  2002   9863 CML    13.996  71.5   71.5  72   2002   2002.5    0     0.5      0     0    71  2002   0.5   208 OC    YOC   Other…
#  2  5253 Female    71  2002   9863 CML    13.996  71.5   72    72.5 2002.5 2003      0.5   1        0     0    72  2002   0.5   208 OC    YOC   Other…
#  3  5253 Female    71  2002   9863 CML    13.996  71.5   72.5  73   2003   2003.5    1     1.5      0     1    72  2003   0.5   208 OC    YOC   Other…
#  4  5253 Female    71  2002   9863 CML    13.996  71.5   73    73.5 2003.5 2004      1.5   2        0     1    73  2003   0.5   208 OC    YOC   Other…
(Dayt=Dayt|>filter(PY>1e-5))  #noise removed => back to size before, so we really didn't need the last split other than to create the time column t
# # A tibble: 724,214 × 23
#       id sex    agedx  yrdx histo3 cancer   surv   adx astart astop ystart  ystop tstart tstop status     t   age  year    PY   COD COD2  COD7  CODS  
#    <int> <fct>  <dbl> <dbl>  <dbl> <chr>   <dbl> <dbl>  <dbl> <dbl>  <dbl>  <dbl>  <dbl> <dbl>  <dbl> <dbl> <dbl> <dbl> <dbl> <int> <chr> <chr> <fct> 
#  1  5253 Female    71  2002   9863 CML    13.996  71.5   71.5  72   2002   2002.5    0     0.5      0     0    71  2002   0.5   208 OC    YOC   Other…
#  2  5253 Female    71  2002   9863 CML    13.996  71.5   72    72.5 2002.5 2003      0.5   1        0     0    72  2002   0.5   208 OC    YOC   Other…
#  3  5253 Female    71  2002   9863 CML    13.996  71.5   72.5  73   2003   2003.5    1     1.5      0     1    72  2003   0.5   208 OC    YOC   Other…
# now take individual data to grouped data using biostat3 function survRate
system.time(day <- survRate(Surv(PY,status)~age+year, data=Dayt)|>tibble()) #2 secs  4.7k
system.time(dayt <- survRate(Surv(PY,status)~age+year+t, data=Dayt)|>tibble()) #21 secs  54k
system.time(dayts <- survRate(Surv(PY,status)~age+year+t+sex, data=Dayt)|>tibble()) #47 secs  89k
system.time(daytsh <- survRate(Surv(PY,status)~age+year+t+sex+histo3, data=Dayt)|>tibble()) #79 secs  127k   
#next two lines are not needed, since yrdx is year - t and agedx is age - t. So output should be no bigger, and compute time shouldn't be in hours
# system.time(daytshy <- survRate(Surv(PY,status)~age+year+t+sex+histo3+yrdx, data=Dayt)|>tibble()) #111 secs 127k  
# system.time(daytshya <- survRate(Surv(PY,status)~age+year+t+sex+histo3+yrdx+agedx, data=Dayt)|>tibble())#hit vec memory limit of 64 Gb,  
# mem.maxVSize() #65GB
# mem.maxVSize(vsize=Inf) # bad idea, took an hour and 50GB of swap space with no end in sight, so killed the job

(day=day|>rename(PY=tstop,O=event)|>select(age:O))
load("~/data/mrt/us_mort.RData")
(m<-us_mort |>filter(Sex == "Total", Year > 1974)|>select(age=Age,year=Year,m=Mortality) ) 
(day=left_join(day,m))
save(day,file="~/data/CMLepi/day.RData")  # 30kb file

(dayt=dayt|>rename(PY=tstop,O=event)|>select(age:O))
(dayt=left_join(dayt,m))
save(dayt,file="~/data/CMLepi/dayt.RData") #219kb file

(m<-us_mort |>filter(Sex != "Total", Year > 1974)|>select(year=Year,age=Age,sex=Sex,m=Mortality) ) 
(dayts=dayts|>rename(PY=tstop,O=event)|>select(age:O))
(dayts=left_join(dayts,m))
save(dayts,file="~/data/CMLepi/dayts.RData") #311 kb file

(daytsh=daytsh|>rename(PY=tstop,O=event)|>select(age:O))
(daytsh=left_join(daytsh,m))
save(daytsh,file="~/data/CMLepi/daytsh.RData") #463 kb file 

# (daytsh=daytsh|>select(age:O)) #keep bak mrt out file
# save(daytsh,file="~/data/CMLepi/daytsh0.RData") #370 kb file 
# (daytsh=daytsh|>mutate(sex=factor(sex),histo3=factor(histo3))) 
# save(daytsh,file="~/data/CMLepi/daytshF.RData") #361 kb not really worth it 
# (daytsh=daytsh|>mutate(age=as.integer(age),year=as.integer(year),t=as.integer(t),O=as.integer(O))) 
# save(daytsh,file="~/data/CMLepi/daytshI.RData") #345 kb, again, not really worth it 


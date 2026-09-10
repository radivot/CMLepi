# AgeYearTime.R   Inspired by biostat3 labs  
graphics.off();rm(list=ls())#clear plots and environment 
library(biostat3)  # loads survival (for Surv and survSplit) 
library(tidyverse)  
options(pillar.sigfig = 5) # Shows 5 significant digits to get decimal after year in tibble prints
load("~/data/CMLepi/cml.RData") #made in mkSEER.R  53.2k
d=d|>mutate(agedx=ifelse(agedx==90,92.3,agedx)) # set over 90 to 92.3
d=d|>mutate(surv=ifelse(surv>80,0.01,surv)) #set NA surv to 0.01 (3.65 days)
d=d|>mutate(surv=ifelse(surv==0,0.01,surv)) #set  0 surv to 0.01 (else survSplit throws 'zero' parameter must be less than any observed times)
min(d$surv)*365# = 1 day = 0.00274 years
d=d|>filter(agedx>=80)|>mutate(histo3=as_factor(histo3))|>select(-cancer,-(COD:CODS)) 
(d=d|>mutate(adx=agedx+0.5,astart=adx,astop=adx+surv))
# A tibble: 7,412 × 10
#      id sex    agedx  yrdx histo3      surv status   adx astart  astop
#   <int> <fct>  <dbl> <dbl> <fct>      <dbl>  <dbl> <dbl>  <dbl>  <dbl>
# 1  98271 Male      87  2005 9863   0.33949        1  87.5   87.5 87.839
# 2 125911 Female    87  2009 9863   0.33128        1  87.5   87.5 87.831
# 3 134825 Female    89  2000 9863   7.3949         1  89.5   89.5 96.895
max(d$astop) #108.1 years
(Da=d|>survSplit(cut = 80:110, event = "status",start = "astart", end = "astop")|>tibble())
# # A tibble: 23,659 × 10
#        id sex    agedx  yrdx histo3    surv   adx astart  astop status
#     <int> <fct>  <dbl> <dbl> <fct>    <dbl> <dbl>  <dbl>  <dbl>  <dbl>
#  1  98271 Male      87  2005 9863   0.33949  87.5   87.5 87.839      1
#  2 125911 Female    87  2009 9863   0.33128  87.5   87.5 87.831      1
#  3 134825 Female    89  2000 9863   7.3949   89.5   89.5 90          0
#  4 134825 Female    89  2000 9863   7.3949   89.5   90   91          0
#  5 134825 Female    89  2000 9863   7.3949   89.5   91   92          0
#  6 134825 Female    89  2000 9863   7.3949   89.5   92   93          0
#  7 134825 Female    89  2000 9863   7.3949   89.5   93   94          0
#  8 134825 Female    89  2000 9863   7.3949   89.5   94   95          0
#  9 134825 Female    89  2000 9863   7.3949   89.5   95   96          0
# 10 134825 Female    89  2000 9863   7.3949   89.5   96   96.895      1

# For each age time band from (a), we calculate the start and stop in calendar time 
# We calculate the time since diagnosis as difference between age at start/stop and 
# age at diagnosis, and add that interval to year at diagnosis
(Da=Da|>mutate(ystart = yrdx + astart - adx, ystop  = yrdx + astop - adx))
# # A tibble: 23,659 × 12
#        id sex    agedx  yrdx histo3    surv   adx astart  astop status ystart  ystop
#     <int> <fct>  <dbl> <dbl> <fct>    <dbl> <dbl>  <dbl>  <dbl>  <dbl>  <dbl>  <dbl>
#  1  98271 Male      87  2005 9863   0.33949  87.5   87.5 87.839      1 2005   2005.3
#  2 125911 Female    87  2009 9863   0.33128  87.5   87.5 87.831      1 2009   2009.3
#  3 134825 Female    89  2000 9863   7.3949   89.5   89.5 90          0 2000   2000.5
#  4 134825 Female    89  2000 9863   7.3949   89.5   90   91          0 2000.5 2001.5
#  5 134825 Female    89  2000 9863   7.3949   89.5   91   92          0 2001.5 2002.5
#  6 134825 Female    89  2000 9863   7.3949   89.5   92   93          0 2002.5 2003.5
#  7 134825 Female    89  2000 9863   7.3949   89.5   93   94          0 2003.5 2004.5
#  8 134825 Female    89  2000 9863   7.3949   89.5   94   95          0 2004.5 2005.5
#  9 134825 Female    89  2000 9863   7.3949   89.5   95   96          0 2005.5 2006.5
# 10 134825 Female    89  2000 9863   7.3949   89.5   96   96.895      1 2006.5 2007.4
# Prioritized +0.5 accuracy in age over year since things tend to change more by age than by calendar time
# ... I recall challenges arising in SEERaBomb when I tried to add +0.5 to both time scales
## Now we can split along the calendar time. 
(Day=Da|>survSplit(cut=1975:2023,event="status",start="ystart",end="ystop")|>tibble()) 
#  A tibble: 36,979 × 12
#        id sex    agedx  yrdx histo3    surv   adx astart  astop ystart  ystop status
#     <int> <fct>  <dbl> <dbl> <fct>    <dbl> <dbl>  <dbl>  <dbl>  <dbl>  <dbl>  <dbl>
#  1  98271 Male      87  2005 9863   0.33949  87.5   87.5 87.839 2005   2005.3      1
#  2 125911 Female    87  2009 9863   0.33128  87.5   87.5 87.831 2009   2009.3      1
#  3 134825 Female    89  2000 9863   7.3949   89.5   89.5 90     2000   2000.5      0
#  4 134825 Female    89  2000 9863   7.3949   89.5   90   91     2000.5 2001        0
#  5 134825 Female    89  2000 9863   7.3949   89.5   90   91     2001   2001.5      0  ## here astart and astop need to be fixed
#  6 134825 Female    89  2000 9863   7.3949   89.5   91   92     2001.5 2002        0
#  7 134825 Female    89  2000 9863   7.3949   89.5   91   92     2002   2002.5      0
#  8 134825 Female    89  2000 9863   7.3949   89.5   92   93     2002.5 2003        0
#  9 134825 Female    89  2000 9863   7.3949   89.5   92   93     2003   2003.5      0
# 10 134825 Female    89  2000 9863   7.3949   89.5   93   94     2003.5 2004        0
(Day=Day|>mutate(astart = adx + ystart - yrdx, astop  = adx + ystop - yrdx)) ## so fix those problems
(Day=Day|>mutate(tstart = astart-adx, tstop  = astop-adx)) ## and bring in tstart and tstop
(Day=Day|>mutate(age=floor(astart),year=floor(ystart),PY=ystop-ystart)|>tibble())#set up getting PY totals in each age-year bin
Day|>print(n=17)
# # A tibble: 36,979 × 17
#        id sex    agedx  yrdx histo3    surv   adx astart  astop ystart  ystop status tstart   tstop   age  year      PY
#     <int> <fct>  <dbl> <dbl> <fct>    <dbl> <dbl>  <dbl>  <dbl>  <dbl>  <dbl>  <dbl>  <dbl>   <dbl> <dbl> <dbl>   <dbl>
#  1  98271 Male      87  2005 9863   0.33949  87.5   87.5 87.839 2005   2005.3      1    0   0.33949    87  2005 0.33949
#  2 125911 Female    87  2009 9863   0.33128  87.5   87.5 87.831 2009   2009.3      1    0   0.33128    87  2009 0.33128
#  3 134825 Female    89  2000 9863   7.3949   89.5   89.5 90     2000   2000.5      0    0   0.5        89  2000 0.5    
#  4 134825 Female    89  2000 9863   7.3949   89.5   90   90.5   2000.5 2001        0    0.5 1          90  2000 0.5    
#  5 134825 Female    89  2000 9863   7.3949   89.5   90.5 91     2001   2001.5      0    1   1.5        90  2001 0.5    
#  6 134825 Female    89  2000 9863   7.3949   89.5   91   91.5   2001.5 2002        0    1.5 2          91  2001 0.5    
#  7 134825 Female    89  2000 9863   7.3949   89.5   91.5 92     2002   2002.5      0    2   2.5        91  2002 0.5    
#  8 134825 Female    89  2000 9863   7.3949   89.5   92   92.5   2002.5 2003        0    2.5 3          92  2002 0.5    
#  9 134825 Female    89  2000 9863   7.3949   89.5   92.5 93     2003   2003.5      0    3   3.5        92  2003 0.5    
# 10 134825 Female    89  2000 9863   7.3949   89.5   93   93.5   2003.5 2004        0    3.5 4          93  2003 0.5    
# 11 134825 Female    89  2000 9863   7.3949   89.5   93.5 94     2004   2004.5      0    4   4.5        93  2004 0.5    
# 12 134825 Female    89  2000 9863   7.3949   89.5   94   94.5   2004.5 2005        0    4.5 5          94  2004 0.5    
# 13 134825 Female    89  2000 9863   7.3949   89.5   94.5 95     2005   2005.5      0    5   5.5        94  2005 0.5    
# 14 134825 Female    89  2000 9863   7.3949   89.5   95   95.5   2005.5 2006        0    5.5 6          95  2005 0.5    
# 15 134825 Female    89  2000 9863   7.3949   89.5   95.5 96     2006   2006.5      0    6   6.5        95  2006 0.5    
# 16 134825 Female    89  2000 9863   7.3949   89.5   96   96.5   2006.5 2007        0    6.5 7          96  2006 0.5    
# 17 134825 Female    89  2000 9863   7.3949   89.5   96.5 96.895 2007   2007.4      1    7   7.3949     96  2007 0.39493
max(d$surv) #21.4 years
(Dayt=Day|>survSplit(cut=1:21,event="status",start="tstart",end="tstop",episode="Time")|>tibble()|>mutate(time=Time-1)) 
Dayt|>print(n=17)
# # A tibble: 38,023 × 19
#        id sex    agedx  yrdx histo3    surv   adx astart  astop ystart  ystop   age  year      PY tstart   tstop status  Time  time
#     <int> <fct>  <dbl> <dbl> <fct>    <dbl> <dbl>  <dbl>  <dbl>  <dbl>  <dbl> <dbl> <dbl>   <dbl>  <dbl>   <dbl>  <dbl> <dbl> <dbl>
#  1  98271 Male      87  2005 9863   0.33949  87.5   87.5 87.839 2005   2005.3    87  2005 0.33949    0   0.33949      1     1     0
#  2 125911 Female    87  2009 9863   0.33128  87.5   87.5 87.831 2009   2009.3    87  2009 0.33128    0   0.33128      1     1     0
#  3 134825 Female    89  2000 9863   7.3949   89.5   89.5 90     2000   2000.5    89  2000 0.5        0   0.5          0     1     0
#  4 134825 Female    89  2000 9863   7.3949   89.5   90   90.5   2000.5 2001      90  2000 0.5        0.5 1            0     1     0
#  5 134825 Female    89  2000 9863   7.3949   89.5   90.5 91     2001   2001.5    90  2001 0.5        1   1.5          0     2     1
#  6 134825 Female    89  2000 9863   7.3949   89.5   91   91.5   2001.5 2002      91  2001 0.5        1.5 2            0     2     1
#  7 134825 Female    89  2000 9863   7.3949   89.5   91.5 92     2002   2002.5    91  2002 0.5        2   2.5          0     3     2
#  8 134825 Female    89  2000 9863   7.3949   89.5   92   92.5   2002.5 2003      92  2002 0.5        2.5 3            0     3     2
#  9 134825 Female    89  2000 9863   7.3949   89.5   92.5 93     2003   2003.5    92  2003 0.5        3   3.5          0     4     3
# 10 134825 Female    89  2000 9863   7.3949   89.5   93   93.5   2003.5 2004      93  2003 0.5        3.5 4            0     4     3
# 11 134825 Female    89  2000 9863   7.3949   89.5   93.5 94     2004   2004.5    93  2004 0.5        4   4.5          0     5     4
# 12 134825 Female    89  2000 9863   7.3949   89.5   94   94.5   2004.5 2005      94  2004 0.5        4.5 5            0     5     4
# 13 134825 Female    89  2000 9863   7.3949   89.5   94.5 95     2005   2005.5    94  2005 0.5        5   5.5          0     6     5
# 14 134825 Female    89  2000 9863   7.3949   89.5   95   95.5   2005.5 2006      95  2005 0.5        5.5 6            0     6     5
# 15 134825 Female    89  2000 9863   7.3949   89.5   95.5 96     2006   2006.5    95  2006 0.5        6   6.5          0     7     6
# 16 134825 Female    89  2000 9863   7.3949   89.5   96   96.5   2006.5 2007      96  2006 0.5        6.5 7            0     7     6
# 17 134825 Female    89  2000 9863   7.3949   89.5   96.5 96.895 2007   2007.4    96  2007 0.39493    7   7.3949       1     8     7
Day|>filter(id==805804) # found example of problem by scrolling
#       id sex    agedx  yrdx histo3   surv   adx astart  astop ystart  ystop status     tstart   tstop   age  year      PY
# 1 805804 Female  92.3  2000 9863   1.5880  92.8 92.800 93     2000   2000.2      0 1.8474e-13 0.20000    92  2000 0.20000
# 2 805804 Female  92.3  2000 9863   1.5880  92.8 93     93.800 2000.2 2001        0 2.0000e- 1 1.0000     93  2000 0.80000
# 3 805804 Female  92.3  2000 9863   1.5880  92.8 93.800 94     2001   2001.2      0 1.0000e+ 0 1.2000     93  2001 0.20000
# 4 805804 Female  92.3  2000 9863   1.5880  92.8 94     94.388 2001.2 2001.6      1 1.2000e+ 0 1.5880     94  2001 0.38795
Dayt|>filter(id==805804)
#       id sex    agedx  yrdx histo3   surv   adx astart  astop ystart  ystop   age  year      PY     tstart   tstop status  Time  time
# 1 805804 Female  92.3  2000 9863   1.5880  92.8 92.800 93     2000   2000.2    92  2000 0.20000 1.8474e-13 0.20000      0     1     0
# 2 805804 Female  92.3  2000 9863   1.5880  92.8 93     93.800 2000.2 2001      93  2000 0.80000 2.0000e- 1 1            0     1     0
# 3 805804 Female  92.3  2000 9863   1.5880  92.8 93     93.800 2000.2 2001      93  2000 0.80000 1     e+ 0 1.0000       0     2     1 row has ~0 PY => remove it
# 4 805804 Female  92.3  2000 9863   1.5880  92.8 93.800 94     2001   2001.2    93  2001 0.20000 1.0000e+ 0 1.2000       0     2     1
# 5 805804 Female  92.3  2000 9863   1.5880  92.8 94     94.388 2001.2 2001.6    94  2001 0.38795 1.2000e+ 0 1.5880       1     2     1
Dayt=Dayt|>mutate(PY=tstop-tstart)  #recompute PY to find ones near zero
(Dayt=Dayt|>filter(PY>1e-5))  #noise removed => back to size before, so we really didn't need the last split other than to create the time column
# # A tibble: 36,979 × 19
#        id sex    agedx  yrdx histo3    surv   adx astart  astop ystart  ystop   age  year      PY tstart   tstop status  Time  time
#     <int> <fct>  <dbl> <dbl> <fct>    <dbl> <dbl>  <dbl>  <dbl>  <dbl>  <dbl> <dbl> <dbl>   <dbl>  <dbl>   <dbl>  <dbl> <dbl> <dbl>
#  1  98271 Male      87  2005 9863   0.33949  87.5   87.5 87.839 2005   2005.3    87  2005 0.33949    0   0.33949      1     1     0
#  2 125911 Female    87  2009 9863   0.33128  87.5   87.5 87.831 2009   2009.3    87  2009 0.33128    0   0.33128      1     1     0
#  3 134825 Female    89  2000 9863   7.3949   89.5   89.5 90     2000   2000.5    89  2000 0.5        0   0.5          0     1     0
#  4 134825 Female    89  2000 9863   7.3949   89.5   90   90.5   2000.5 2001      90  2000 0.5        0.5 1            0     1     0
#  5 134825 Female    89  2000 9863   7.3949   89.5   90.5 91     2001   2001.5    90  2001 0.5        1   1.5          0     2     1
#  6 134825 Female    89  2000 9863   7.3949   89.5   91   91.5   2001.5 2002      91  2001 0.5        1.5 2            0     2     1
#  7 134825 Female    89  2000 9863   7.3949   89.5   91.5 92     2002   2002.5    91  2002 0.5        2   2.5          0     3     2
#  8 134825 Female    89  2000 9863   7.3949   89.5   92   92.5   2002.5 2003      92  2002 0.5        2.5 3            0     3     2
system.time(D <- survRate(Surv(PY,status)~time+age+year+sex, data=Dayt)|>tibble()) #3 secs
(D=D|>rename(PY=tstop,O=event))
# # A tibble: 6,969 × 9
#     time   Age  Year Sex          PY     O     rate     lower    upper
#    <dbl> <dbl> <dbl> <fct>     <dbl> <dbl>    <dbl>     <dbl>    <dbl>
#  1     0    80  1975 Female 1.9117       3  1.5693   0.32362    4.5861
#  2     0    80  1975 Male   0.60951      2  3.2813   0.39738   11.853 
#  3     0    80  1976 Female 0.021903     2 91.313   11.058    329.85  
#  4     0    80  1976 Male   1.4682       2  1.3622   0.16497    4.9209
#  5     0    80  1977 Female 0.66256      2  3.0186   0.36557   10.904 
#  6     0    80  1977 Male   1            0  0        0          3.6889
#  7     0    80  1978 Male   0.63142      2  3.1675   0.38360   11.442 
#  8     0    80  1979 Female 1.8532       1  0.53961  0.013662   3.0065
#  9     0    80  1979 Male   0.56297      1  1.7763   0.044972   9.8969
# 10     0    80  1980 Female 0.51000      1  1.9608   0.049643  10.925 
range(D$Age) #80 to 108
range(D$Year)#1975 to 2023
range(D$time)#0 to 21. These are the left ends of time intervals.  21 include all bigger times
(d80=D|>select(time:O))
load("~/data/mrt/us_mort.RData")
(m<-us_mort |>filter(Sex != "Total", Year > 1974,Age>=80)|>select(year=Year,age=Age,sex=Sex,mrt=Mortality) ) 
(d80=left_join(d80,m))
save(d80,file="~/data/CMLepi/cml80.RData")

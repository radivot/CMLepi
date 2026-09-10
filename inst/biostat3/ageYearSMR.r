# AgeYearSMR.R   Inspired by lab 23  
graphics.off();rm(list=ls())#clear plots and environment 
library(biostat3)  # loads survival (for Surv and survSplit) 
library(tidyverse)   
load("~/data/CMLepi/cml.RData") #made in mkSEER.R  53.2k
d=d|>mutate(agedx=ifelse(agedx==90,92.3,agedx)) # set over 90 to 92.3
d=d|>mutate(surv=ifelse(surv>80,0.1,surv)) #set NA surv to 0.1
d=d|>mutate(surv=ifelse(surv==0,0.1,surv)) #set  0 surv to 0.1 (else survSplit throws 'zero' parameter must be less than any observed times)
d=d|>filter(agedx>=80)|>mutate(histo3=as_factor(histo3))|>select(-cancer,-(COD:CODS)) 
(d=d|>mutate(adx=agedx+0.5,astart=adx,astop=adx+surv))
max(d$astop) #108.1 years
(Da=d|>survSplit(cut = 80:110, event = "status",start = "astart", end = "astop")|>tibble())
# # A tibble: 23,659 × 10
#        id sex    agedx  yrdx histo3  surv   adx astart astop status
#     <int> <fct>  <dbl> <dbl> <fct>  <dbl> <dbl>  <dbl> <dbl>  <dbl>
#  1  98271 Male      87  2005 9863   0.339  87.5   87.5  87.8      1
#  2 125911 Female    87  2009 9863   0.331  87.5   87.5  87.8      1
#  3 134825 Female    89  2000 9863   7.39   89.5   89.5  90        0
#  4 134825 Female    89  2000 9863   7.39   89.5   90    91        0
#  5 134825 Female    89  2000 9863   7.39   89.5   91    92        0
#  6 134825 Female    89  2000 9863   7.39   89.5   92    93        0
#  7 134825 Female    89  2000 9863   7.39   89.5   93    94        0
#  8 134825 Female    89  2000 9863   7.39   89.5   94    95        0
#  9 134825 Female    89  2000 9863   7.39   89.5   95    96        0
# 10 134825 Female    89  2000 9863   7.39   89.5   96    96.9      1

# For each age time band from (a), we calculate the start and stop in calendar time 
# We calculate the time since diagnosis as difference between age at start/stop and 
# age at diagnosis, and add that interval to year at diagnosis
(Da=Da|>mutate(ystart = num(yrdx + astart - adx,digits=1), ystop  = num(yrdx + astop - adx,digits=1)))
#  A tibble: 23,659 × 12
#        id sex    agedx  yrdx histo3  surv   adx astart astop status    ystart     ystop
#     <int> <fct>  <dbl> <dbl> <fct>  <dbl> <dbl>  <dbl> <dbl>  <dbl> <num:.1!> <num:.1!>
#  1  98271 Male      87  2005 9863   0.339  87.5   87.5  87.8      1    2005.0    2005.3
#  2 125911 Female    87  2009 9863   0.331  87.5   87.5  87.8      1    2009.0    2009.3
#  3 134825 Female    89  2000 9863   7.39   89.5   89.5  90        0    2000.0    2000.5
#  4 134825 Female    89  2000 9863   7.39   89.5   90    91        0    2000.5    2001.5
#  5 134825 Female    89  2000 9863   7.39   89.5   91    92        0    2001.5    2002.5
#  6 134825 Female    89  2000 9863   7.39   89.5   92    93        0    2002.5    2003.5
#  7 134825 Female    89  2000 9863   7.39   89.5   93    94        0    2003.5    2004.5
#  8 134825 Female    89  2000 9863   7.39   89.5   94    95        0    2004.5    2005.5
#  9 134825 Female    89  2000 9863   7.39   89.5   95    96        0    2005.5    2006.5
# 10 134825 Female    89  2000 9863   7.39   89.5   96    96.9      1    2006.5    2007.4
# Prioritized +0.5 accuracy in age over year since things tend to change more by age than by calendar time
# ... I recall challenges arising in SEERaBomb when I tried to add +0.5 to both time scales
## Now we can split along the calendar time. 
(Day=Da|>survSplit(cut=1975:2023,event="status",start="ystart",end="ystop")) #tibble now back to data.frame
#        id    sex agedx yrdx histo3        surv  adx astart    astop ystart    ystop status
# 1   98271   Male  87.0 2005   9863 0.339493498 87.5   87.5 87.83949 2005.0 2005.339      1
# 2  125911 Female  87.0 2009   9863 0.331279945 87.5   87.5 87.83128 2009.0 2009.331      1
# 3  134825 Female  89.0 2000   9863 7.394934976 89.5   89.5 90.00000 2000.0 2000.500      0
# 4  134825 Female  89.0 2000   9863 7.394934976 89.5   90.0 91.00000 2000.5 2001.000      0
# 5  134825 Female  89.0 2000   9863 7.394934976 89.5   90.0 91.00000 2001.0 2001.500      0 ## here astart and astop need to be fixed
# 6  134825 Female  89.0 2000   9863 7.394934976 89.5   91.0 92.00000 2001.5 2002.000      0
# 7  134825 Female  89.0 2000   9863 7.394934976 89.5   91.0 92.00000 2002.0 2002.500      0 ## same here
# 8  134825 Female  89.0 2000   9863 7.394934976 89.5   92.0 93.00000 2002.5 2003.000      0
# 9  134825 Female  89.0 2000   9863 7.394934976 89.5   92.0 93.00000 2003.0 2003.500      0
# 10 134825 Female  89.0 2000   9863 7.394934976 89.5   93.0 94.00000 2003.5 2004.000      0
# 11 134825 Female  89.0 2000   9863 7.394934976 89.5   93.0 94.00000 2004.0 2004.500      0
# 12 134825 Female  89.0 2000   9863 7.394934976 89.5   94.0 95.00000 2004.5 2005.000      0
# 13 134825 Female  89.0 2000   9863 7.394934976 89.5   94.0 95.00000 2005.0 2005.500      0
# 14 134825 Female  89.0 2000   9863 7.394934976 89.5   95.0 96.00000 2005.5 2006.000      0
# 15 134825 Female  89.0 2000   9863 7.394934976 89.5   95.0 96.00000 2006.0 2006.500      0
# 16 134825 Female  89.0 2000   9863 7.394934976 89.5   96.0 96.89493 2006.5 2007.000      0
# 17 134825 Female  89.0 2000   9863 7.394934976 89.5   96.0 96.89493 2007.0 2007.395      1
(Day=Day|>mutate(astart = adx + ystart - yrdx, astop  = adx + ystop - yrdx)) ## so fix those problems here
(Day=Day|>mutate(age=floor(astart),year=floor(ystart),PY=ystop-ystart))#set up getting PY totals in each age-year bin
Day=Day|>tibble()|>mutate(ystart=num(ystart,digits=1),ystop=num(ystop,digits=1)) # and get back to tibbles
Day|>print(n=17)
# # A tibble: 36,979 × 15
#        id sex    agedx  yrdx histo3  surv   adx astart astop    ystart     ystop status   age  year    PY
#     <int> <fct>  <dbl> <dbl> <fct>  <dbl> <dbl>  <dbl> <dbl> <num:.1!> <num:.1!>  <dbl> <dbl> <dbl> <dbl>
#  1  98271 Male      87  2005 9863   0.339  87.5   87.5  87.8    2005.0    2005.3      1    87  2005 0.339
#  2 125911 Female    87  2009 9863   0.331  87.5   87.5  87.8    2009.0    2009.3      1    87  2009 0.331
#  3 134825 Female    89  2000 9863   7.39   89.5   89.5  90      2000.0    2000.5      0    89  2000 0.5  
#  4 134825 Female    89  2000 9863   7.39   89.5   90    90.5    2000.5    2001.0      0    90  2000 0.5  
#  5 134825 Female    89  2000 9863   7.39   89.5   90.5  91      2001.0    2001.5      0    90  2001 0.5  
#  6 134825 Female    89  2000 9863   7.39   89.5   91    91.5    2001.5    2002.0      0    91  2001 0.5  
#  7 134825 Female    89  2000 9863   7.39   89.5   91.5  92      2002.0    2002.5      0    91  2002 0.5  
#  8 134825 Female    89  2000 9863   7.39   89.5   92    92.5    2002.5    2003.0      0    92  2002 0.5  
#  9 134825 Female    89  2000 9863   7.39   89.5   92.5  93      2003.0    2003.5      0    92  2003 0.5  
# 10 134825 Female    89  2000 9863   7.39   89.5   93    93.5    2003.5    2004.0      0    93  2003 0.5  
# 11 134825 Female    89  2000 9863   7.39   89.5   93.5  94      2004.0    2004.5      0    93  2004 0.5  
# 12 134825 Female    89  2000 9863   7.39   89.5   94    94.5    2004.5    2005.0      0    94  2004 0.5  
# 13 134825 Female    89  2000 9863   7.39   89.5   94.5  95      2005.0    2005.5      0    94  2005 0.5  
# 14 134825 Female    89  2000 9863   7.39   89.5   95    95.5    2005.5    2006.0      0    95  2005 0.5  
# 15 134825 Female    89  2000 9863   7.39   89.5   95.5  96      2006.0    2006.5      0    95  2006 0.5  
# 16 134825 Female    89  2000 9863   7.39   89.5   96    96.5    2006.5    2007.0      0    96  2006 0.5  
# 17 134825 Female    89  2000 9863   7.39   89.5   96.5  96.9    2007.0    2007.4      1    96  2007 0.395

xtabs(PY ~ age + year, data=Day, subset = age<100 & year>=2014)
#     year
# age        2014       2015       2016       2017       2018       2019       2020       2021       2022       2023
#   80 13.7844627 14.4487337 12.3586585 17.8050650 12.7065024 19.3217659 11.5000000 14.0284052 18.1225188 13.5049281
#   81 38.2629706 36.7857632 37.7344969 37.8204654 36.9688569 40.2156057 38.6170431 38.6260096 43.8724162 24.5510609
#   82 46.5975359 55.4948665 58.3505133 51.6722108 53.2617385 57.9312115 51.8227926 56.6074606 54.7438056 30.7046543
#   83 56.7427105 61.8700205 73.2501711 64.0798084 56.4673511 66.8947981 69.3052704 60.3720055 65.2852841 37.4765229
#   84 63.1438741 64.8895277 73.7920602 74.8813142 68.0130732 67.3111567 75.7412731 76.6131417 72.5322382 38.8727584
#   85 57.7772758 65.8963723 70.0863107 72.4151951 76.7937029 73.4873374 70.6208077 74.1496235 82.5948665 42.2095825
#   86 71.7924709 63.7267625 69.4281999 68.3221081 73.8158795 77.4618070 74.7229295 66.7989049 74.8507187 41.2005476
#   87 78.5164956 74.2759069 70.8973306 65.4473648 67.9112252 77.5777550 70.7228611 71.4367556 66.8013005 35.7422998
#   88 59.9373032 70.9377139 73.6577002 62.6718686 58.2238877 67.5852841 68.8459959 63.7583847 65.5451745 33.2369610
#   89 54.0708419 54.3646133 67.9104038 70.0842574 56.4453799 55.1283368 62.0650924 62.4568789 58.9418207 32.5308693
#   90 55.9425051 45.7665982 45.7847365 59.7522245 58.1208077 48.1143053 46.4377139 52.6235455 49.0465435 21.6057495
#   91 33.4493498 40.0766598 34.2724162 35.0677618 44.1088296 44.5383299 36.5749487 33.6526352 38.5273785 19.8559206
#   92 25.7553730 31.0455852 40.5167009 35.4301848 35.6512663 40.9876112 37.6581793 37.2718001 33.3377823 22.9553730
#   93 36.6882957 36.2538672 45.3084189 49.9292950 44.0845311 43.7585216 42.3332649 52.6045859 52.7481862 23.8880219
#   94 26.6106092 23.0147159 22.7466804 31.1828884 30.3590007 26.0566051 31.6394935 30.6855578 37.3711841 18.7430527
#   95 19.7719370 19.1045175 16.8273785 16.7734428 21.7654346 22.6086927 17.7462697 19.2306639 22.1575633 11.4205339
#   96 10.4394251 15.4103354 12.9574949 10.7104723 10.7737166 15.2405886 13.9756331 10.0738535 11.4913758  6.5633128
#   97  3.5067077  7.0727584 10.4236140  7.7761123  8.1865161  8.1780287 11.1650924  8.1602327  6.8027379  3.8464750
#   98  3.0180014  1.5000000  4.7941821  7.9825462  5.7716632  4.8198494  6.0803559  6.4765229  5.1190965  2.4865161
#   99  0.9454483  1.2000000  1.5000000  3.1000000  6.3021218  3.2117728  1.1146475  2.1775496  4.0390828  1.9808350
xtabs(status ~ age + year, data=Day, subset = age<100 & year>=2014) ## count deaths 
#     year
# age  2014 2015 2016 2017 2018 2019 2020 2021 2022 2023
#   80    5    5    4    5    5    4    0    5    5    6
#   81   12   10   12   12    5   13    9    6    7    9
#   82   12   12   12   18   16   11   12    7   19    6
#   83   12   15   16   21   14   14   13   13   17    9
#   84   19   17   12   15   13   13   16   13   12   10
#   85   15   13   19   23   20   16   18   19   14   15
#   86   26   14   19   18   23   22   22   18   27   12
#   87   15   24   15   21   10   24   19   18   10   13
#   88   26   21   22   21   28   16   25   23   18   16
#   89   21   14   18   13   16   16   21   11   19    6
#   90   15   19    7   23   20    8   14   14   15    6
#   91   10   14    8    6   12   10   12   14    9    7
#   92   24   27   26   24   21   27   31   26   21   21
#   93   21   19   19   22   28   24   10   22   15   13
#   94   11   13    7    9   13   11   10   13    9    6
#   95    6    2    8    6   11    7    9   10    7    2
#   96    5    7    4    6    2    4   11    4    4    4
#   97    3    1    5    2    2    2    4    2    1    0
#   98    2    0    4    2    3    3    6    5    2    1
#   99    1    0    0    0    1    1    1    2    2    1
(Day=Day|>mutate(age10=cut(age,c(80,90,110), right=FALSE),year10 = cut(year, c(1975,1985,2000,2010,2024), right=FALSE)))
survRate(Surv(PY, status) ~ year10+age10, data=Day)|>tibble() 
#   year10      age10     tstop event  rate lower upper
#   <fct>       <fct>     <dbl> <dbl> <dbl> <dbl> <dbl>
# 1 [1975,1985) [80,90)   403.    368 0.914 0.823 1.01 
# 2 [1975,1985) [90,110)   74.8    76 1.02  0.800 1.27 
# 3 [1985,2000) [80,90)   969.    621 0.641 0.591 0.693
# 4 [1985,2000) [90,110)  211.    196 0.930 0.804 1.07 
# 5 [2000,2010) [80,90)  3208.   1433 0.447 0.424 0.470
# 6 [2000,2010) [90,110)  745.    521 0.699 0.640 0.762
# 7 [2010,2024) [80,90)  7529.   2019 0.268 0.257 0.280
# 8 [2010,2024) [90,110) 3054.   1266 0.415 0.392 0.438

##### rest is SMR 
load("~/data/mrt/us_mort.RData")
head(m<-us_mort |>filter(Sex != "Total", Year > 1974) )
#    Year   Age OpenInterval Sex    Deaths Exposures Population Mortality
#   <int> <int> <lgl>        <chr>   <dbl>     <dbl>      <dbl>     <dbl>
# 1  1975     0 FALSE        Female 21717.  1533580.   1517559.  0.0142    21717./1533580.
# 2  1975     1 FALSE        Female  1474.  1476741.   1486948.  0.000998
# 3  1975     2 FALSE        Female   989.  1507341.   1535882.  0.000656
# 4  1975     3 FALSE        Female   776.  1609269.   1675307.  0.000482
tail(m)
#    Year   Age OpenInterval Sex   Deaths Exposures Population Mortality
#   <int> <int> <lgl>        <chr>  <dbl>     <dbl>      <dbl>     <dbl>
# 1  2024   105 FALSE        Male    163.     325.       302.      0.501
# 2  2024   106 FALSE        Male     95      164.       147.      0.580
# 3  2024   107 FALSE        Male     45       75.2       69.0     0.599
# 4  2024   108 FALSE        Male     21       35.6       32.9     0.590
# 5  2024   109 FALSE        Male      6       18.3       16.1     0.327
# 6  2024   110 TRUE         Male     11       30.0       28.5     0.367
# Population: The total count of people alive at the start of the interval 
# Exposure: PY from count, modified to account for the exact timing of events.
163/325 #  0.5015385
(popmort=m|>mutate(sex=ifelse(Sex=="Male",1,2))|>select(sex,age=Age,year=Year,rate=Mortality))
#     sex   age  year     rate
#   <dbl> <int> <int>    <dbl>
# 1     2     0  1975 0.0142  
# 2     2     1  1975 0.000998
# 3     2     2  1975 0.000656
# 4     2     3  1975 0.000482
# 5     2     4  1975 0.000419

(pt=Day|>mutate(sex=unclass(sex))|>    # make sex integer to be in line with popmort 
  group_by(sex, age, year)|>summarise(PY=sum(PY),O=sum(status),.groups="keep")|>ungroup())  

(joint=left_join(pt, popmort))
# # A tibble: 1,745 × 6
#      sex   age  year     PY     O   rate
#    <dbl> <dbl> <dbl>  <dbl> <dbl>  <dbl>
#  1     1    80  1975 1.91       3 0.101 
#  2     1    80  1976 0.0219     2 0.102 
#  3     1    80  1977 0.663      2 0.0991
#  4     1    80  1979 1.85       1 0.0952
#  5     1    80  1980 0.600      1 0.0987
(joint <- mutate(joint, E = PY*rate)) 

calculate_smr = function(data)
  summarise(data,
            Observed = sum(O),
            Expected = sum(E)) |>
  mutate(SMR=Observed/Expected,
         poisson.ci(Observed,Expected))

calculate_smr(joint)
# # A tibble: 1 × 5
#   Observed Expected   SMR `2.5 %` `97.5 %`
#      <dbl>    <dbl> <dbl>   <dbl>    <dbl>
# 1     6500    2133.  3.05    2.97     3.12
SMR_byYear=joint|>group_by(year)|>calculate_smr()
library(tinyplot)
with(SMR_byYear,plt(SMR~year, type = "ribbon", ymin=`2.5 %`, ymax=`97.5 %`,ylim=c(0,20)))

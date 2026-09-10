# splitAgeTime.R   Inspired by lab 22  
graphics.off();rm(list=ls())#clear plots and environment 
library(biostat3)  # loads survival (for Surv and survSplit) 
library(tidyverse)   
load("~/data/CMLepi/cml20.RData") #made in mkSEER.R
(d=d20|>filter(histo3%in%c(9863,9875),agedx<90,surv<80,surv>0,yrdx>=2014)) # 21.3k CML cases
(d=d|>filter(agedx>=80)|>mutate(histo3=as_factor(histo3))|>select(-cancer,-(COD:CODS))) #2.1k
# # A tibble: 2,082 × 7
#         id sex    agedx  yrdx histo3   surv status
#      <int> <fct>  <int> <dbl> <fct>   <dbl>  <dbl>
#  1  233367 Female    81  2021 9875   2.09        0
#  2  320365 Male      87  2015 9875   0.0986      1
#  3  371865 Male      89  2017 9875   3.38        1
#  4  482907 Male      88  2014 9875   7.46        1
#  5  733354 Female    82  2021 9863   0.575       1
#  6  871256 Female    83  2016 9875   7.15        0
#  7  978678 Male      86  2017 9875   0.192       1
#  8 1013229 Male      80  2019 9863   3.32        1
#  9 1269332 Female    82  2014 9863   5.66        1
# 10 1296321 Male      88  2014 9875   0.271       1
max(d$surv) #9.99 years
time.cuts <- c(0,1,5,10)
(Dt=survSplit(Surv(surv,status)~.,data=d,cut=time.cuts,end="tstop",start="tstart",event="status")|>tibble())
# # A tibble: 3,732 × 8
#        id sex    agedx  yrdx histo3 tstart  tstop status
#     <int> <fct>  <int> <dbl> <fct>   <dbl>  <dbl>  <dbl>
#  1 233367 Female    81  2021 9875        0 1           0
#  2 233367 Female    81  2021 9875        1 2.09        0
#  3 320365 Male      87  2015 9875        0 0.0986      1
#  4 371865 Male      89  2017 9875        0 1           0
#  5 371865 Male      89  2017 9875        1 3.38        1
#  6 482907 Male      88  2014 9875        0 1           0
#  7 482907 Male      88  2014 9875        1 5           0
#  8 482907 Male      88  2014 9875        5 7.46        1
#  9 733354 Female    82  2021 9863        0 0.575       1
# 10 871256 Female    83  2016 9875        0 1           0
(Dt=Dt|>mutate(age_start=agedx+tstart,age_end=agedx+tstop))
# # A tibble: 3,732 × 10
#        id sex    agedx  yrdx histo3 tstart  tstop status age_start age_end
#     <int> <fct>  <int> <dbl> <fct>   <dbl>  <dbl>  <dbl>     <dbl>   <dbl>
#  1 233367 Female    81  2021 9875        0 1           0        81    82  
#  2 233367 Female    81  2021 9875        1 2.09        0        82    83.1
#  3 320365 Male      87  2015 9875        0 0.0986      1        87    87.1
#  4 371865 Male      89  2017 9875        0 1           0        89    90  
#  5 371865 Male      89  2017 9875        1 3.38        1        90    92.4
#  6 482907 Male      88  2014 9875        0 1           0        88    89  
#  7 482907 Male      88  2014 9875        1 5           0        89    93  
#  8 482907 Male      88  2014 9875        5 7.46        1        93    95.5
#  9 733354 Female    82  2021 9863        0 0.575       1        82    82.6
# 10 871256 Female    83  2016 9875        0 1           0        83    84  

age_cat <- c(80,85,90,100) # Split at these ages
(Dta=survSplit(Dt,cut=age_cat,start="age_start",end="age_end",event="status")|>tibble())
# # A tibble: 4,190 × 10
#        id sex    agedx  yrdx histo3 tstart  tstop age_start age_end status
#     <int> <fct>  <int> <dbl> <fct>   <dbl>  <dbl>     <dbl>   <dbl>  <dbl>
#  1 233367 Female    81  2021 9875        0 1             81    82        0
#  2 233367 Female    81  2021 9875        1 2.09          82    83.1      0
#  3 320365 Male      87  2015 9875        0 0.0986        87    87.1      1
#  4 371865 Male      89  2017 9875        0 1             89    90        0
#  5 371865 Male      89  2017 9875        1 3.38          90    92.4      1
#  6 482907 Male      88  2014 9875        0 1             88    89        0
#  7 482907 Male      88  2014 9875        1 5             89    90        0
#  8 482907 Male      88  2014 9875        1 5             90    93        0
#  9 482907 Male      88  2014 9875        5 7.46          93    95.5      1
# 10 733354 Female    82  2021 9863        0 0.575         82    82.6      1
(Dta=mutate(Dta,PY = age_end- age_start, # Creating new time at risk
               age = cut(age_end, age_cat)))   # Creating age band category
# # A tibble: 4,190 × 12
#        id sex    agedx  yrdx histo3 tstart  tstop age_start age_end status     PY age     
#     <int> <fct>  <int> <dbl> <fct>   <dbl>  <dbl>     <dbl>   <dbl>  <dbl>  <dbl> <fct>   
#  1 233367 Female    81  2021 9875        0 1             81    82        0 1      (80,85] 
#  2 233367 Female    81  2021 9875        1 2.09          82    83.1      0 1.09   (80,85] 
#  3 320365 Male      87  2015 9875        0 0.0986        87    87.1      1 0.0986 (85,90] 
#  4 371865 Male      89  2017 9875        0 1             89    90        0 1      (85,90] 
#  5 371865 Male      89  2017 9875        1 3.38          90    92.4      1 2.38   (90,100]
#  6 482907 Male      88  2014 9875        0 1             88    89        0 1      (85,90] 
#  7 482907 Male      88  2014 9875        1 5             89    90        0 1      (85,90] 
#  8 482907 Male      88  2014 9875        1 5             90    93        0 3      (90,100]
#  9 482907 Male      88  2014 9875        5 7.46          93    95.5      1 2.46   (90,100]
# 10 733354 Female    82  2021 9863        0 0.575         82    82.6      1 0.575  (80,85] 

## Calculate crude rates
survRate(Surv(PY, status) ~ age+sex+histo3, data=Dta)|>tibble()
# # A tibble: 12 × 8
#    age      sex    histo3 tstop event  rate lower upper
#    <fct>    <fct>  <fct>  <dbl> <dbl> <dbl> <dbl> <dbl>
#  1 (80,85]  Female 9863   582.    136 0.234 0.196 0.277
#  2 (80,85]  Female 9875   550.     73 0.133 0.104 0.167
#  3 (80,85]  Male   9863   615.    188 0.306 0.264 0.353 in all pairs
#  4 (80,85]  Male   9875   611.    121 0.198 0.164 0.236 male rates are higher
#  5 (85,90]  Female 9863   605.    180 0.297 0.256 0.344
#  6 (85,90]  Female 9875   545.    118 0.217 0.179 0.259
#  7 (85,90]  Male   9863   577.    218 0.378 0.329 0.432
#  8 (85,90]  Male   9875   500.    130 0.260 0.217 0.309
#  9 (90,100] Female 9863   152.     46 0.302 0.221 0.403
# 10 (90,100] Female 9875   133.     39 0.294 0.209 0.401
# 11 (90,100] Male   9863    88.3    34 0.385 0.267 0.538
# 12 (90,100] Male   9875   120.     39 0.325 0.231 0.444   in all pairs 9875 has a lower rate

survRate(Surv(PY, status) ~ age, data=Dta)|>tibble()
#   age      tstop event  rate lower upper
#   <fct>    <dbl> <dbl> <dbl> <dbl> <dbl>
# 1 (80,85]  2359.   518 0.220 0.201 0.239
# 2 (85,90]  2226.   646 0.290 0.268 0.313
# 3 (90,100]  494.   158 0.320 0.272 0.374  death rates go up with age

survRate(Surv(PY, status) ~ sex, data=Dta)|>tibble()
#   sex    tstop event  rate lower upper
#   <fct>  <dbl> <dbl> <dbl> <dbl> <dbl>
# 1 Female 2567.   592 0.231 0.212 0.250
# 2 Male   2511.   730 0.291 0.270 0.313 #higher in males, perhaps mostly via background being higher

summary(poi <- glm(status ~ age + sex + histo3 + offset(log(PY)), family=poisson,data=Dta))
eform(poi)
#             exp(beta)     2.5 %    97.5 %
# (Intercept) 0.2254754 0.2009132 0.2530405
# age(85,90]  1.3212625 1.1769238 1.4833030
# age(90,100] 1.5079998 1.2614428 1.8027480
# sexMale     1.2884622 1.1558110 1.4363376
# histo39875  0.6899669 0.6178334 0.7705221
summary(coxph(Surv(age_start, age_end, status) ~ sex +  histo3, data = Dta))
#                coef exp(coef) se(coef)      z Pr(>|z|)    
# sexMale     0.25083   1.28509  0.05552  4.518 6.24e-06 ***
# histo39875 -0.36777   0.69227  0.05645 -6.515 7.25e-11 ***

survRate(Surv(PY, status) ~ tstart, data=Dta)
#          tstart     tstop event      rate     lower     upper
# tstart=0      0 1617.8645   565 0.3492258 0.3210198 0.3792461
# tstart=1      1 2926.9391   622 0.2125087 0.1961341 0.2298855
# tstart=5      5  533.8248   135 0.2528920 0.2120334 0.2993283 rise back up via aging


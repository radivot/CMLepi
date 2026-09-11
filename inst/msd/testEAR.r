# testEAR.R   
graphics.off();rm(list=ls())#clear plots and environment 
library(tidyverse)  
load("~/data/CMLepi/cml.RData") #made in mkSEER.R  53.2k
d=d|>mutate(agedx=ifelse(agedx==90,92.3,agedx)) # set over 90 to 92.3
d=d|>mutate(surv=ifelse(surv>80,0.01,surv)) #set NA surv to 0.01 (3.65 days)
d=d|>mutate(surv=ifelse(surv==0,0.01,surv)) #set  0 surv to 0.01 (else survSplit throws 'zero' parameter must be less than any observed times)
load("~/data/mrt/mrtUSA.RData")#mrt is list of 3 matrices
d=d|>select(yrdx,agedx,sex,surv,status) # 43,932 CML cases
(D=SEERaBomb::msd(d,mrt,brkst=c(0)))
# # A tibble: 2 × 12
#   int         t     O     E      PY    EAR     LL     UL    RR   rrL   rrU sex   
#   <fct>   <dbl> <dbl> <dbl>   <dbl>  <dbl>  <dbl>  <dbl> <dbl> <dbl> <dbl> <chr> 
# 1 (0,100]  3.23 14913 4055. 194493. 0.0558 0.0546 0.0571  3.68  3.62  3.74 Male  
# 2 (0,100]  3.32 10918 2702. 153663. 0.0535 0.0521 0.0548  4.04  3.96  4.12 Female
SEERaBomb::foldD(D,keep=c("int"))
#  int         O     E      PY     t    EAR     LL     UL    RR   rrL   rrU
# 1 (0,100] 25831 6757. 348155.  3.28 0.0548 0.0539 0.0557  3.82  3.78  3.87

load("~/data/CMLepi/day.RData")#made by mkBigCML.R
day=day|>mutate(E=PY*m)
d0=day|>summarize(O=sum(O),E=sum(E),PY=sum(PY))
d0|>mutate(EAR=(O-E)/PY,LL=EAR-1.96*sqrt(O)/PY,UL=EAR+1.96*sqrt(O)/PY,
             RR=O/E,rrL=qchisq(.025,2*O)/(2*E),rrU=qchisq(.975,2*O+2)/(2*E)) 
#       O      E      PY      EAR       LL       UL     RR    rrL    rrU
#   <dbl>  <dbl>   <dbl>    <dbl>    <dbl>    <dbl>  <dbl>  <dbl>  <dbl>
# 1 25831 6837.0 348155. 0.054556 0.053651 0.055461 3.7781 3.7322 3.8245
# so O and PY are at the same, but E is off via m 

load("~/data/CMLepi/dayts.RData") #try separate sexes
dayts=dayts|>mutate(E=PY*m)
d0=dayts|>group_by(sex)|>summarize(O=sum(O),E=sum(E),PY=sum(PY))
d0|>mutate(EAR=(O-E)/PY,LL=EAR-1.96*sqrt(O)/PY,UL=EAR+1.96*sqrt(O)/PY,
             RR=O/E,rrL=qchisq(.025,2*O)/(2*E),rrU=qchisq(.975,2*O+2)/(2*E)) 
#   sex        O      E      PY      EAR       LL       UL     RR    rrL    rrU
#   <chr>  <dbl>  <dbl>   <dbl>    <dbl>    <dbl>    <dbl>  <dbl>  <dbl>  <dbl>
# 1 Female 10918 2840.8 153663. 0.052565 0.051232 0.053897 3.8433 3.7715 3.9161
# 2 Male   14913 4239.0 194493. 0.054881 0.053651 0.056112 3.5181 3.4618 3.5750
# so O and PY are still fine, but E is still off via m 

(f23=mrt[["Female"]][,"2023"])
load("~/data/mrt/us_mort.RData")
(d23<-us_mort|>filter(Sex == "Female", Year == 2023)|>select(age=Age,year=Year,m=Mortality) ) 
d23$f23=f23
View(d23) #identical in 2023

(f75=mrt[["Female"]][,"1975"])
(d75<-us_mort|>filter(Sex == "Female", Year == 1975)|>select(age=Age,year=Year,m=Mortality) ) 
d75$f75=f75
View(d75) #identical in 1975 also ... so what gives?




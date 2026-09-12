# testPreShift.R   #takehome:  shifting agedx upfront for msd equates to shifting both agedx and yrdx upfront for survSplit
graphics.off();rm(list=ls())#clear plots and environment 
library(tidyverse)  
load("~/data/CMLepi/cml.RData") #made in mkSEER.R  53.2k
d=d|>mutate(agedx=agedx+0.5)
d=d|>mutate(agedx=ifelse(agedx>90,92.3,agedx)) # set over 90 to 92.3
d=d|>mutate(surv=ifelse(surv>80,0.01,surv)) #set NA surv to 0.01 (3.65 days)
d=d|>mutate(surv=ifelse(surv==0,0.01,surv)) #set  0 surv to 0.01 (else survSplit throws 'zero' parameter must be less than any observed times)
d=d|>select(yrdx,agedx,sex,surv,status) # 43,932 CML cases
load("~/data/mrt/mrtUSA.RData")#mrt is list of 3 matrices
(D=SEERaBomb::msd(d,mrt,brkst=c(0)))
#  A tibble: 2 × 12
#   int         t     O     E      PY    EAR     LL     UL    RR   rrL   rrU sex   
#   <fct>   <dbl> <dbl> <dbl>   <dbl>  <dbl>  <dbl>  <dbl> <dbl> <dbl> <dbl> <chr> 
# 1 (0,100]  3.23 14913 4212. 194493. 0.0550 0.0538 0.0563  3.54  3.48  3.60 Male  
# 2 (0,100]  3.32 10918 2822. 153663. 0.0527 0.0514 0.0540  3.87  3.80  3.94 Female
load("~/data/CMLepi/daysPreShift.RData")
days=days|>mutate(E=PY*m)
d0=days|>group_by(sex)|>summarize(O=sum(O),E=sum(E),PY=sum(PY))
d0|>mutate(EAR=(O-E)/PY,LL=EAR-1.96*sqrt(O)/PY,UL=EAR+1.96*sqrt(O)/PY,
             RR=O/E,rrL=qchisq(.025,2*O)/(2*E),rrU=qchisq(.975,2*O+2)/(2*E)) 
# # A tibble: 2 × 10
#   sex        O     E      PY    EAR     LL     UL    RR   rrL   rrU
#   <chr>  <dbl> <dbl>   <dbl>  <dbl>  <dbl>  <dbl> <dbl> <dbl> <dbl>
# 1 Female 10918 2831. 153663. 0.0526 0.0513 0.0540  3.86  3.78  3.93
# 2 Male   14913 4231. 194493. 0.0549 0.0537 0.0562  3.52  3.47  3.58
# So, while O and PY are still the same,  E is now bigger: 4212=>4231 and 2822=>2831. 
range(dayts$year)

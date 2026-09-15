# testD0.R   #goal here is to see how close we can get them
graphics.off();rm(list=ls())#clear plots and environment 
library(tidyverse)  
options(pillar.sigfig = 6) # Show 6 significant digits 
load("~/data/CMLepi/cml.RData") #made in mkSEER.R  53254
(d=d|>filter(yrdx>=2000,agedx<90,surv<80,surv>0)) #43932
d=d|>mutate(agedx=agedx+0.5)
load("~/data/mrt/mrtUSA.RData")#mrt is list of 3 matrices
SEERaBomb::msd(d,mrt,brkst=0)|>select(O:sex)|>relocate(sex,.before=O)|>arrange(sex)
#   sex        O       E      PY       EAR        LL        UL      RR     rrL     rrU
# 1 Female  7089 2316.13 135286. 0.0352799 0.0340601 0.0364998 3.06071 2.98987 3.13280
# 2 Male   10132 3514.47 170732. 0.0387596 0.0376041 0.0399152 2.88294 2.82707 2.93963
load("~/data/CMLepi/D0.RData") 
dd=D0|>mutate(E=PY*m)
dd=dd|>group_by(sex)|>summarize(O=sum(O),E=sum(E),PY=sum(PY))
dd|>mutate(EAR=(O-E)/PY,LL=EAR-1.96*sqrt(O)/PY,UL=EAR+1.96*sqrt(O)/PY,
            RR=O/E,rrL=qchisq(.025,2*O)/(2*E),rrU=qchisq(.975,2*O+2)/(2*E)) 
#   sex        O       E      PY       EAR        LL        UL      RR     rrL     rrU
#   <chr>  <dbl>   <dbl>   <dbl>     <dbl>     <dbl>     <dbl>   <dbl>   <dbl>   <dbl>
# 1 Female  7089 2316.13 135286. 0.0352799 0.0340601 0.0364998 3.06071 2.98987 3.13280
# 2 Male   10132 3514.47 170732. 0.0387596 0.0376041 0.0399152 2.88294 2.82707 2.93963
# There are NO round off errors!
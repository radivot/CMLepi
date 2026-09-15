# testPreShift.R  Shows that adding 0.5 to agedx upfront with msd equals adding it to both agedx and yrdx upfront for survSplit
graphics.off();rm(list=ls())#clear plots and environment 
library(tidyverse)  
options(pillar.sigfig = 6) # Show 6 significant digits 
load("~/data/CMLepi/cml.RData") #made in mkSEER.R  53.2k
d=d|>mutate(agedx=agedx+0.5)
d=d|>mutate(agedx=ifelse(agedx>90,92.3,agedx)) # set over 90 to 92.3
d=d|>mutate(surv=ifelse(surv>80,0.01,surv)) #set NA surv to 0.01 (3.65 days)
d=d|>mutate(surv=ifelse(surv==0,0.01,surv)) #set  0 surv to 0.01 (else survSplit throws 'zero' parameter must be less than any observed times)
d=d|>select(yrdx,agedx,sex,surv,status) # 43,932 CML cases
load("~/data/mrt/mrtUSA.RData")#mrt is list of 3 matrices
SEERaBomb::msd(d,mrt,brkst=0)|>select(O:sex)|>relocate(sex,.before=O)|>arrange(sex)
#   sex        O       E      PY       EAR        LL        UL      RR     rrL     rrU
# 1 Female 10918 2822.12 153663. 0.0526861 0.0513533 0.0540188 3.86872 3.79649 3.94198
# 2 Male   14913 4211.69 194493. 0.0550217 0.0537911 0.0562524 3.54086 3.48426 3.59815
load("~/data/CMLepi/daysPreShift.RData")
days=days|>mutate(E=PY*m)
d0=days|>group_by(sex)|>summarize(O=sum(O),E=sum(E),PY=sum(PY))
d0|>mutate(EAR=(O-E)/PY,LL=EAR-1.96*sqrt(O)/PY,UL=EAR+1.96*sqrt(O)/PY,
             RR=O/E,rrL=qchisq(.025,2*O)/(2*E),rrU=qchisq(.975,2*O+2)/(2*E)) 
#   sex        O       E      PY       EAR        LL        UL      RR     rrL     rrU
# 1 Female 10918 2821.94 153663. 0.0526872 0.0513544 0.0540200 3.86897 3.79673 3.94223
# 2 Male   14913 4211.44 194493. 0.0550230 0.0537923 0.0562536 3.54107 3.48446 3.59836
### there are tiny diffs in E
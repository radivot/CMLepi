# testNoShift.R   Shows that if there are no agedx shifts, SEERaBomb's matrix approach matches the biostat3 approach  
graphics.off();rm(list=ls())#clear plots and environment 
library(tidyverse)  
options(pillar.sigfig = 6) # Show 6 significant digits 
load("~/data/CMLepi/cml.RData") #made in mkSEER.R  53.2k
d=d|>mutate(agedx=ifelse(agedx==90,92.3,agedx)) # set over 90 to 92.3
d=d|>mutate(surv=ifelse(surv>80,0.01,surv)) #set NA surv to 0.01 (3.65 days)
d=d|>mutate(surv=ifelse(surv==0,0.01,surv)) #set  0 surv to 0.01 (else survSplit throws 'zero' parameter must be less than any observed times)
load("~/data/mrt/mrtUSA.RData")#mrt is list of 3 matrices
d=d|>select(yrdx,agedx,sex,surv,status) # 43,932 CML cases
SEERaBomb::msd(d,mrt,brkst=0)|>select(O:sex)|>relocate(sex,.before=O)|>arrange(sex)
#  sex        O       E      PY       EAR        LL        UL      RR     rrL     rrU
# 1 Female 10918 2702.49 153663. 0.0534646 0.0521318 0.0547974 4.03998 3.96455 4.11648
# 2 Male   14913 4054.66 194493. 0.0558291 0.0545984 0.0570597 3.67799 3.61920 3.73751
load("~/data/CMLepi/daytsNoShift.RData") #adx=agedx, not adx=agedx+0.5 
dayts=dayts|>mutate(E=PY*m)
d0=dayts|>group_by(sex)|>summarize(O=sum(O),E=sum(E),PY=sum(PY))
d0|>mutate(EAR=(O-E)/PY,LL=EAR-1.96*sqrt(O)/PY,UL=EAR+1.96*sqrt(O)/PY,
             RR=O/E,rrL=qchisq(.025,2*O)/(2*E),rrU=qchisq(.975,2*O+2)/(2*E)) 
#  sex        O       E      PY       EAR        LL        UL      RR     rrL     rrU
# 1 Female 10918 2702.55 153663. 0.0534642 0.0521314 0.0547970 4.03989 3.96446 4.11639
# 2 Male   14913 4055.09 194493. 0.0558269 0.0545962 0.0570575 3.67760 3.61881 3.73711
# small diffs in E values could be round off errors, i.e. it seems msd() has been working pretty fine 


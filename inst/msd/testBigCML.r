# testBigCML.R   Using preShifts to make things match
graphics.off();rm(list=ls())#clear plots and environment 
library(tidyverse)  
load("~/data/CMLepi/cml.RData") #made in mkSEER.R  53.2k
d=d|>mutate(agedx=agedx+0.5)
d=d|>mutate(agedx=ifelse(agedx>90,92.3,agedx)) # set over 90 to 92.3
d=d|>mutate(surv=ifelse(surv>80,0.01,surv)) #set NA surv to 0.01 (3.65 days)
d=d|>mutate(surv=ifelse(surv==0,0.01,surv)) #set  0 surv to 0.01 (else survSplit throws 'zero' parameter must be less than any observed times)
d=d|>select(yrdx,agedx,sex,surv,status) # 43,932 CML cases
load("~/data/mrt/mrtUSA.RData")#mrt is list of 3 matrices
SEERaBomb::msd(d,mrt,brkst=0)|>select(O:sex)|>relocate(sex,.before=O)|>arrange(sex)
#   sex        O     E      PY    EAR     LL     UL    RR   rrL   rrU
#   <chr>  <dbl> <dbl>   <dbl>  <dbl>  <dbl>  <dbl> <dbl> <dbl> <dbl>
# 1 Female 10918 2822. 153663. 0.0527 0.0514 0.0540  3.87  3.80  3.94
# 2 Male   14913 4212. 194493. 0.0550 0.0538 0.0563  3.54  3.48  3.60

### now go through files with separate sexes made in mkBigCML.R
load("~/data/CMLepi/Days.RData") # 79kb file
load("~/data/CMLepi/Dayts.RData") #311 kb file
load("~/data/CMLepi/Daytsh.RData")#463 kb file 
getEARs=function(Ds) {
  dd=Ds|>mutate(E=PY*m)
  dd=dd|>group_by(sex)|>summarize(O=sum(O),E=sum(E),PY=sum(PY))
  dd|>mutate(EAR=(O-E)/PY,LL=EAR-1.96*sqrt(O)/PY,UL=EAR+1.96*sqrt(O)/PY,
              RR=O/E,rrL=qchisq(.025,2*O)/(2*E),rrU=qchisq(.975,2*O+2)/(2*E)) 
}
getEARs(Days)
#   sex        O     E      PY    EAR     LL     UL    RR   rrL   rrU
#   <chr>  <dbl> <dbl>   <dbl>  <dbl>  <dbl>  <dbl> <dbl> <dbl> <dbl>
# 1 Female 10918 2822. 153663. 0.0527 0.0514 0.0540  3.87  3.80  3.94
# 2 Male   14913 4211. 194493. 0.0550 0.0538 0.0563  3.54  3.48  3.60
getEARs(Dayts) #same
getEARs(Daytsh) #same => good so far

# now do the other two files
load("~/data/CMLepi/Day.RData")  # 58kb file
load("~/data/CMLepi/Dayt.RData") #220kb file
getEAR=function(Ds) {
  dd=Ds|>mutate(E=PY*m)
  dd=dd|>summarize(O=sum(O),E=sum(E),PY=sum(PY))
  dd|>mutate(EAR=(O-E)/PY,LL=EAR-1.96*sqrt(O)/PY,UL=EAR+1.96*sqrt(O)/PY,
              RR=O/E,rrL=qchisq(.025,2*O)/(2*E),rrU=qchisq(.975,2*O+2)/(2*E)) 
}
getEAR(Day)
#       O     E      PY    EAR     LL     UL    RR   rrL   rrU
#   <dbl> <dbl>   <dbl>  <dbl>  <dbl>  <dbl> <dbl> <dbl> <dbl>
# 1 25831 6796. 348155. 0.0547 0.0538 0.0556  3.80  3.75  3.85
getEAR(Dayt) #same
# and SEERabomb gives us
D=SEERaBomb::msd(d,mrt,brkst=0)
SEERaBomb::foldD(D)|>select(-int,-t)
#       O     E      PY    EAR     LL     UL    RR   rrL   rrU
#   <dbl> <dbl>   <dbl>  <dbl>  <dbl>  <dbl> <dbl> <dbl> <dbl>
# 1 25831 7034. 348155. 0.0540 0.0531 0.0549  3.67  3.63  3.72
# E is the correct here, not above.  Why? Sex-pooled morts used by survSplit()
# above weigh F > M at old ages due to women living longer, but CML incidence is
# ~60% higher for males, so this pop has a different sex ratio. 
# Lesson: Use sex-specific morts and add up Es, Os and PYs later.


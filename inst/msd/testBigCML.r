# testBigCML.R   Using preShifts to make things match
graphics.off();rm(list=ls())#clear plots and environment 
library(tidyverse)  
options(pillar.sigfig = 5) # Show 5 significant digits 
load("~/data/CMLepi/cml.RData") #made in mkSEER.R  53.2k
d=d|>mutate(agedx=agedx+0.5)
d=d|>mutate(agedx=ifelse(agedx>90,92.3,agedx)) # set over 90 to 92.3
d=d|>mutate(surv=ifelse(surv>80,0.01,surv)) #set NA surv to 0.01 (3.65 days)
d=d|>mutate(surv=ifelse(surv==0,0.01,surv)) #set  0 surv to 0.01 (else survSplit throws 'zero' parameter must be less than any observed times)
d=d|>select(yrdx,agedx,sex,surv,status) # 43,932 CML cases
load("~/data/mrt/mrtUSA.RData")#mrt is list of 3 matrices
SEERaBomb::msd(d,mrt,brkst=0)|>select(O:sex)|>relocate(sex,.before=O)|>arrange(sex)
#   sex        O      E      PY      EAR       LL       UL     RR    rrL    rrU
# 1 Female 10918 2822.1 153663. 0.052686 0.051353 0.054019 3.8687 3.7965 3.9420
# 2 Male   14913 4211.7 194493. 0.055022 0.053791 0.056252 3.5409 3.4843 3.5982
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
#  A tibble: 2 × 10
#   sex        O      E      PY      EAR       LL       UL     RR    rrL    rrU
# 1 Female 10918 2821.9 153663. 0.052687 0.051354 0.054020 3.8690 3.7967 3.9422 
# 2 Male   14913 4211.4 194493. 0.055023 0.053792 0.056254 3.5411 3.4845 3.5984 
# E diffs of 0.2 and 0.3 could be round off errors
getEARs(Dayts)  #same
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
#       O      E      PY      EAR       LL       UL     RR    rrL    rrU
# 1 25831 6795.7 348155. 0.054675 0.053770 0.055580 3.8011 3.7549 3.8477
#         6795.7 is not equal to 2821.9 + 4211.4 = 7033.3, so we really need to use sex specific mortalities
getEAR(Dayt) #same  #Note that sex-pooled mortalities give more weight to women at old ages, so E is low

D=SEERaBomb::msd(d,mrt,brkst=0) # SEERabomb always uses sex specific morts  
SEERaBomb::foldD(D)|>select(-int,-t)  # and then folds them (adds Es) afterwards
#       O      E      PY      EAR       LL       UL     RR    rrL    rrU
# 1 25831 7033.8 348155. 0.053991 0.053086 0.054896 3.6724 3.6278 3.7175
#         7033.8 = 2822.1 + 4211.7, which is what we want
# E is the correct here, not above.  Why? Sex-pooled morts used by survSplit()
# above weigh F > M at old ages due to women living longer, but CML incidence is
# ~60% higher for males, so the CML patient population has a different sex ratio. 
# Lesson: Use sex-specific morts and add up Es, Os and PYs later.


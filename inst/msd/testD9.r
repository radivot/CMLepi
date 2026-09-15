# testD9.R  see mkD5.R 
graphics.off();rm(list=ls())#clear plots and environment 
library(tidyverse)  
options(pillar.sigfig = 8) 
load("~/data/CMLepi/cml.RData") #53,254   made in mkSEER.R  
(d=d|>filter(yrdx>=2000)) #45,636 
d=d|>mutate(agedx=agedx+0.5,yrdx=yrdx+0.5) #need to shift both to match preshift of just age in SEERaBomb matrix approach which ties them naturally 
d=d|>mutate(agedx=ifelse(agedx>90,92.3,agedx)) # set over 90 to 92.3 and exclude all others
(d=d|>filter(yrdx>=2000,surv<80,surv>0)) #44,879    
d=d|>select(yrdx,agedx,sex,surv,status) # 43,932 CML cases
load("~/data/mrt/mrtUSA.RData")#mrt is list of 3 matrices
SEERaBomb::msd(d,mrt,brkst=0)|>select(O:sex)|>relocate(sex,.before=O)|>arrange(sex)
#  sex        O         E        PY         EAR          LL          UL        RR       rrL       rrU
# 1 Female  7579 2478.5186 136011.94 0.037500247 0.036245705 0.038754788 3.0578750 2.9894146 3.1275077
# 2 Male   10499 3657.7468 171268.56 0.039944594 0.038771987 0.041117201 2.8703463 2.8157011 2.9257853
load("~/data/CMLepi/D92.3.RData")
dd=D9|>mutate(E=PY*m)
dd=dd|>group_by(sex)|>summarize(O=sum(O),E=sum(E),PY=sum(PY))
dd|>mutate(EAR=(O-E)/PY,LL=EAR-1.96*sqrt(O)/PY,UL=EAR+1.96*sqrt(O)/PY,
            RR=O/E,rrL=qchisq(.025,2*O)/(2*E),rrU=qchisq(.975,2*O+2)/(2*E)) 
#   sex        O         E        PY         EAR          LL          UL        RR       rrL       rrU
# 1 Female  7579 2478.3248 136011.94 0.037501671 0.036247130 0.038756213 3.0581141 2.9896483 3.1277522
# 2 Male   10499 3657.5054 171268.56 0.039946004 0.038773397 0.041118610 2.8705357 2.8158869 2.9259784
###### so slightly off already

load("~/data/CMLepi/cml.RData") 
(d=d|>filter(yrdx>=2000)) 
d=d|>mutate(agedx=agedx+0.5,yrdx=yrdx+0.5)
d=d|>mutate(agedx=ifelse(agedx>90,92.5,agedx)) ############ try keeping it on the age grid at 92.5 ####
(d=d|>filter(yrdx>=2000,surv<80,surv>0)) 
d=d|>select(yrdx,agedx,sex,surv,status) 
load("~/data/mrt/mrtUSA.RData")
SEERaBomb::msd(d,mrt,brkst=0)|>select(O:sex)|>relocate(sex,.before=O)|>arrange(sex)
#   sex        O         E        PY         EAR          LL          UL        RR       rrL       rrU
# 1 Female  7579 2481.5342 136011.94 0.037478075 0.036223534 0.038732617 3.0541590 2.9857818 3.1237071
# 2 Male   10499 3660.1241 171268.56 0.039930714 0.038758107 0.041103320 2.8684820 2.8138723 2.9238850
load("~/data/CMLepi/D92.5.RData")
dd=D9|>mutate(E=PY*m)
dd=dd|>group_by(sex)|>summarize(O=sum(O),E=sum(E),PY=sum(PY))
dd|>mutate(EAR=(O-E)/PY,LL=EAR-1.96*sqrt(O)/PY,UL=EAR+1.96*sqrt(O)/PY,
            RR=O/E,rrL=qchisq(.025,2*O)/(2*E),rrU=qchisq(.975,2*O+2)/(2*E)) 
#   sex        O         E        PY         EAR          LL          UL        RR       rrL       rrU
# 1 Female  7579 2481.5342 136011.94 0.037478075 0.036223534 0.038732617 3.0541590 2.9857818 3.1237071
# 2 Male   10499 3660.1241 171268.56 0.039930714 0.038758107 0.041103320 2.8684820 2.8138723 2.9238850
## bingo, staying on the grid made the difference

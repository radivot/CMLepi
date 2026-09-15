# testD5.R  see mkD5.R 
graphics.off();rm(list=ls())#clear plots and environment 
library(tidyverse)  
options(pillar.sigfig = 8) 
load("~/data/CMLepi/cml.RData") #53,254   made in mkSEER.R  
(d=d|>filter(yrdx>=2000)) #45,636 
d=d|>mutate(agedx=agedx+0.5,yrdx=yrdx+0.5) #need to shift both to match preshift of just age in SEERaBomb matrix approach which ties them naturally 
d=d|>mutate(agedx=ifelse(agedx>90,92.5,agedx)) # set over 90 to 92.5 
#if surv = NA (i.e. is >80), set it to 5 (no healthcare received), but first push back yrdx and agedx by 5 years 
d=d|>mutate(yrdx=ifelse(surv>80,yrdx-5,yrdx)) 
d=d|>mutate(agedx=ifelse(surv>80,agedx-5,agedx)) 
d|>filter(agedx<1) # one patient was 3.5 at Dx (now -1.5), so set this to 0.5 with a surv of 3 years
d=d|>mutate(surv=ifelse(agedx<0,3,surv)) #set this patient's survival time to 3 years  
d=d|>mutate(agedx=ifelse(agedx<0,0.5,agedx)) #and set negative age to 0.5 
d=d|>mutate(surv=ifelse(surv>80,5,surv)) # now set the rest of the NA survs to 5 years   
d=d|>mutate(surv=ifelse(surv==0,0.01,surv)) #set  0 surv to 0.01 (else survSplit throws 'zero' parameter must be less than any observed times)
d=d|>select(yrdx,agedx,sex,surv,status) # 43,932 CML cases
load("~/data/mrt/mrtUSA.RData")#mrt is list of 3 matrices
SEERaBomb::msd(d,mrt,brkst=0)|>select(O:sex)|>relocate(sex,.before=O)|>arrange(sex)
#   sex        O         E        PY         EAR          LL          UL        RR       rrL       rrU
# 1 Female  7893 2575.2368 137255.92 0.038743415 0.037474753 0.040012077 3.0649608 2.9977130 3.1333369
# 2 Male   10855 3757.0375 172689.82 0.041102379 0.039919871 0.042284888 2.8892446 2.8351448 2.9441171
load("~/data/CMLepi/D5.RData")
dd=D5|>mutate(E=PY*m)
dd=dd|>group_by(sex)|>summarize(O=sum(O),E=sum(E),PY=sum(PY))
dd|>mutate(EAR=(O-E)/PY,LL=EAR-1.96*sqrt(O)/PY,UL=EAR+1.96*sqrt(O)/PY,
            RR=O/E,rrL=qchisq(.025,2*O)/(2*E),rrU=qchisq(.975,2*O+2)/(2*E)) 
#   sex        O         E        PY         EAR          LL          UL        RR       rrL       rrU
# 1 Female  7893 2575.2368 137255.92 0.038743415 0.037474753 0.040012077 3.0649608 2.9977130 3.1333369
# 2 Male   10855 3757.0375 172689.82 0.041102379 0.039919871 0.042284888 2.8892446 2.8351448 2.9441171

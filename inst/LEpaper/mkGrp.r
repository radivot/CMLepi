# mkGrp.R  EE = exclude both;  KK = keep both (of S=NA and S=0)  
graphics.off();rm(list=ls())#clear plots and environment 
library(survival)
library(tidyverse)  
library(mgcv)  # if not here, predict uses lm, which throws "lm object does not have a proper 'qr' component"
options(pillar.sigfig = 5) # Shows 5 significant digits to get decimal after year in tibble prints
load("~/data/CMLepi/G12.RData") # made in mkMorts.R
load("~/data/CMLepi/cml.RData") #made in mkSEER.R  53.2k
d=d|>mutate(COD12=as.factor(COD12)) 
levels(d$COD12) # alive" "ASH"   "BEN"   "CA"    "COPD"  "CV"    "DK"    "ILL"   "IN"    "LC"    "MCA"   "MSPC"  "OCD"  
(d0=d|>filter(yrdx>=2000)) #45636  
d0|>filter(agedx==90) #1,085 × 12
d0|>filter(surv>80) #533 × 12
d0|>filter(surv==0) #224 × 12
(dKK=d0|>filter(agedx<90)) #44551 
dKK|>filter(surv>80) #416 × 12
dKK|>filter(surv==0) #203 × 12
dKK=dKK|>mutate(agedx=agedx+0.5,yrdx=yrdx+0.5) # keep in mind that doing this pushes some people's surv times into 2024.  
sum(dKK$status) # 17754 deaths
i="KKold"  # keep everything
i="EE"
i="DIFF"
i="OLD"
i="KK"
if (i=="EE") # exlude both
  (d=dKK|>filter(surv<80,surv>0)) #43932
if (i=="KK") {
  d=dKK|>mutate(surv=ifelse(surv>80,0.001,surv)) 
  d=d|>mutate(surv=ifelse(surv==0,0.001,surv)) 
} 
if (i=="DIFF") {
  d=dKK|>filter(surv>80 | surv==0) # 619 gets all since can't have both in same person
  d=d|>mutate(surv=0.001) 
} 
if (i=="OLD") {
  d=d0|>filter(agedx==90) # 1085
  d|>filter(surv>80) #117 × 12
  d|>filter(surv==0) #21 × 12  
  d=d|>mutate(surv=ifelse(surv>80,0.001,surv)) # go with KK as making the most sense
  d=d|>mutate(surv=ifelse(surv==0,0.001,surv)) 
  d=d|>mutate(agedx=92.6,yrdx=yrdx+0.5) # keep in mind that doing this pushes some people's surv times beyond 2024 (will pool all as if 2024)  
} 

if (i=="KKold") {
  d=d0|>mutate(agedx=agedx+0.5,yrdx=yrdx+0.5) 
  d=d|>mutate(agedx=ifelse(agedx>90,92.3,agedx)) # set over 90 (those at 90.5) to 92.3
  d=d|>mutate(surv=ifelse(surv>80,0.001,surv)) # go with KK as making the most sense
  d=d|>mutate(surv=ifelse(surv==0,0.001,surv)) 
} 

(d=d|>mutate(adx=agedx,astart=adx,astop=adx+surv,.before=COD))
(Da=survSplit(Surv(astart,astop,COD12)~.,d,cut = 20:110)|>tibble()|>relocate(astart:COD12,.before=COD)) 
(Da=Da|>mutate(ystart = yrdx + astart - adx, ystop  = yrdx + astop - adx,.before=COD))
(Day=survSplit(Surv(ystart,ystop,COD12)~.,Da,cut=2000:2024)|>tibble()|>relocate(ystart:COD12,.before=COD)) 
(Day=Day|>mutate(astart = adx + ystart - yrdx, astop  = adx + ystop - yrdx)) # fix problems
(Day=Day|>mutate(tstart = astart-adx, tstop  = astop-adx,.before=status)) ## and bring in tstart and tstop
(Dayt=survSplit(Surv(tstart,tstop,COD12)~.,Day,cut=1:24,episode="Time")|>tibble()|>relocate(tstart:COD12,.before=COD)) 
(Dayt=Dayt|>mutate(t=Time-1,age=floor(astart),year=floor(ystart),PY=tstop-tstart,.before=COD)|>select(-Time))

(Dayt=Dayt|>mutate(ASH=ifelse(COD12=="ASH",1,0),.before=COD))
(Dayt=Dayt|>mutate(BEN=ifelse(COD12=="BEN",1,0),.before=COD))
(Dayt=Dayt|>mutate(CA=ifelse(COD12=="CA",1,0),.before=COD))
(Dayt=Dayt|>mutate(COPD=ifelse(COD12=="COPD",1,0),.before=COD))
(Dayt=Dayt|>mutate(CV=ifelse(COD12=="CV",1,0),.before=COD))
(Dayt=Dayt|>mutate(DK=ifelse(COD12=="DK",1,0),.before=COD))
(Dayt=Dayt|>mutate(ILL=ifelse(COD12=="ILL",1,0),.before=COD))
(Dayt=Dayt|>mutate(IN=ifelse(COD12=="IN",1,0),.before=COD))
(Dayt=Dayt|>mutate(LC=ifelse(COD12=="LC",1,0),.before=COD))
(Dayt=Dayt|>mutate(MCA=ifelse(COD12=="MCA",1,0),.before=COD))
(Dayt=Dayt|>mutate(MSPC=ifelse(COD12=="MSPC",1,0),.before=COD))
(Dayt=Dayt|>mutate(OCD=ifelse(COD12=="OCD",1,0),.before=COD))
system.time(dASH<-biostat3::survRate(Surv(PY,ASH)~age+year+t+sex, data=Dayt)|>tibble())   
system.time(dBEN<-biostat3::survRate(Surv(PY,BEN)~age+year+t+sex, data=Dayt)|>tibble())   
system.time(dCA<-biostat3::survRate(Surv(PY,CA)~age+year+t+sex, data=Dayt)|>tibble())   
system.time(dCOPD<-biostat3::survRate(Surv(PY,COPD)~age+year+t+sex, data=Dayt)|>tibble())   
system.time(dCV<-biostat3::survRate(Surv(PY,CV)~age+year+t+sex, data=Dayt)|>tibble())   
system.time(dDK<-biostat3::survRate(Surv(PY,DK)~age+year+t+sex, data=Dayt)|>tibble())   
system.time(dILL<-biostat3::survRate(Surv(PY,ILL)~age+year+t+sex, data=Dayt)|>tibble())   
system.time(dIN<-biostat3::survRate(Surv(PY,IN)~age+year+t+sex, data=Dayt)|>tibble())   
system.time(dLC<-biostat3::survRate(Surv(PY,LC)~age+year+t+sex, data=Dayt)|>tibble())   
system.time(dMCA<-biostat3::survRate(Surv(PY,MCA)~age+year+t+sex, data=Dayt)|>tibble())   
system.time(dMSPC<-biostat3::survRate(Surv(PY,MSPC)~age+year+t+sex, data=Dayt)|>tibble())   
system.time(dOCD<-biostat3::survRate(Surv(PY,OCD)~age+year+t+sex, data=Dayt)|>tibble())   
(dASH=dASH|>rename(PY=tstop,ASH=event)|>select(age:ASH))
(dBEN=dBEN|>rename(PY=tstop,BEN=event)|>select(age:BEN))
(dCA=dCA|>rename(PY=tstop,CA=event)|>select(age:CA))
(dCOPD=dCOPD|>rename(PY=tstop,COPD=event)|>select(age:COPD))
(dCV=dCV|>rename(PY=tstop,CV=event)|>select(age:CV))
(dDK=dDK|>rename(PY=tstop,DK=event)|>select(age:DK))
(dILL=dILL|>rename(PY=tstop,ILL=event)|>select(age:ILL))
(dIN=dIN|>rename(PY=tstop,IN=event)|>select(age:IN))
(dLC=dLC|>rename(PY=tstop,LC=event)|>select(age:LC))
(dMCA=dMCA|>rename(PY=tstop,MCA=event)|>select(age:MCA))
(dMSPC=dMSPC|>rename(PY=tstop,MSPC=event)|>select(age:MSPC))
(dOCD=dOCD|>rename(PY=tstop,OCD=event)|>select(age:OCD))
(D12=left_join(dASH,dBEN))
(D12=left_join(D12,dCA))
(D12=left_join(D12,dCOPD))
(D12=left_join(D12,dCV))
(D12=left_join(D12,dDK))
(D12=left_join(D12,dILL))
(D12=left_join(D12,dIN))
(D12=left_join(D12,dLC))
(D12=left_join(D12,dMCA))
(D12=left_join(D12,dMSPC))
(D12=left_join(D12,dOCD))
(D=D12|>filter(age>20)) # G fits in mkG13 are only for ages over 20
D$is_2020 <- as.numeric(D$year == 2020)
D$is_2021 <- as.numeric(D$year == 2021)
D$is_2022 <- as.numeric(D$year == 2022)
save(D,file=paste0("~/data/CMLepi/D_",i,".RData")) #save intermediate since above takes time
D=D|>mutate(num=ASH,denom=PY)
D$Eash=as.numeric(exp(predict(G[["ASH"]],D,qr=TRUE)))
D=D|>mutate(num=BEN)
D$Eben=as.numeric(exp(predict(G[["BEN"]],D)))
D=D|>mutate(num=CA)
D$Eca=as.numeric(exp(predict(G[["CA"]],D)))
D=D|>mutate(num=COPD)
D$Ecopd=as.numeric(exp(predict(G[["COPD"]],D)))
D=D|>mutate(num=CV)
D$Ecv=as.numeric(exp(predict(G[["CV"]],D)))
D=D|>mutate(num=DK)
D$Edk=as.numeric(exp(predict(G[["DK"]],D)))
D=D|>mutate(num=ILL)
D$Eill=as.numeric(exp(predict(G[["ILL"]],D)))
D=D|>mutate(num=IN)
D$Ein=as.numeric(exp(predict(G[["IN"]],D)))
D=D|>mutate(num=LC)
D$Elc=as.numeric(exp(predict(G[["LC"]],D)))
D=D|>mutate(num=MCA)
D$Emca=as.numeric(exp(predict(G[["MCA"]],D)))
D=D|>mutate(num=MSPC)
D$Emspc=as.numeric(exp(predict(G[["MSPC"]],D)))
D=D|>mutate(num=OCD)
D$Eocd=as.numeric(exp(predict(G[["OCD"]],D)))
D=D|>select(-denom,-num)
D=D|>select(-is_2020,-is_2021,-is_2022)
load("~/data/mrt/us_mort.RData") 
(m<-us_mort |>filter(Sex != "Total", Year > 1974)|>select(year=Year,age=Age,sex=Sex,m=Mortality) ) 
(D=left_join(D,m))
D=D|>mutate(O=ASH+BEN+CA+COPD+CV+DK+ILL+IN+LC+MCA+MSPC+OCD,E=m*PY, # O and E are for all cause (AC) mortality
            Echk=Eash+Eben+Eca+Ecopd+Ecv+Edk+Eill+Ein+Elc+Emca+Emspc+Eocd) # check that fits sum into alighment with HMD totals
save(D,file=paste0("~/data/CMLepi/grp_",i,".RData"))
D|>select(ASH:OCD,O)|>colSums()
  # ASH   BEN    CA  COPD    CV    DK   ILL    IN    LC   MCA  MSPC   OCD     O 
  # 472   432  1359   502  3390   634   105   585  7610   302   296  1982 17669 
# 17669  is close to 17754 deaths at top. 17754-17669= 85 deaths missing
D12|>filter(age<=20)|>select(ASH:OCD)|>colSums()|>sum() # yes, 85 young deaths cut out
D|>select(E,Echk,O)|>colSums()
 #        E      Echk         O 
 # 5828.062  5799.980 17669.000   
# expectations 30 off is not bad since fits were not constrained to sum to the total, and totals are SEER vs HMD  


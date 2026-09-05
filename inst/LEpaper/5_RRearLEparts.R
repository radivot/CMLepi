## 5_RRearLEparts.R
graphics.off();rm(list=ls())#clear plots and environment 
library(tidyverse)
library(vital)  
library(tsibble)
load("~/data/mrt/us_mort.RData") #us_mort is a class vital object (single tibble-like)  
load("~/data/mrt/mrtUSA.RData")#mrt is list of 3 matrices
load("~/data/CMLepi/cml20.RData") #made in mkSEER.R
(d20=d20|>filter(histo3%in%c(9863,9875),agedx<90,surv<80,surv>0)) #43932
d20=d20|>select(yrdx,agedx,sex,surv,status) # 43,932 CML cases

####  A and B use all ages, so d20 is used in the line below 
(D=SEERaBomb::msd(d20,mrt,brkst=c(0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21)))
D=D|>rename(Group="sex")|>select(Group,int,everything())
(D=SEERaBomb::foldD(D,keep=c("int")))

#### plotting acronyms used in 4 plots below (i.e. 5ABCD)
leg=theme(legend.margin=margin(0,0,0,0),legend.position="top",legend.text=element_text(size=12))
gx=xlab("Years Since Dx")
gyE=ylab("Excess Absolute Risk of Death")
cc1=coord_cartesian(xlim=c(0.0,40))
gh0=geom_hline(yintercept=0)
gh1=geom_hline(yintercept=1)
gh1p8=geom_hline(yintercept=1.8,col="gray")
ghp017=geom_hline(yintercept=0.017,col="gray")
geRR=geom_errorbar(aes(ymin=rrL,ymax=rrU),width=.2)
geE=geom_errorbar(aes(ymin=LL,ymax=UL),width=0.2)#for absolute risks
gl=geom_line()
gp=geom_point()
tc=function(sz) theme_classic(base_size=sz)

D|>ggplot(aes(x=t,y=RR)) +gp+gl+gx+gh1+gh1p8+tc(13)+leg+scale_y_continuous(breaks=c(1,1.8,5,10,15))+geRR+
  scale_color_manual(values = c("gray","black"),name = "Year of Dx")+ylab("Relative Risk of Death")
ggsave("LE/outs/5A_RRallAges.pdf",width=2.5,height=3)

D|>ggplot(aes(x=t,y=EAR))+gp+gl+gx+gyE+gh0+ghp017+tc(13)+leg+geE+
  scale_x_continuous(breaks=c(0,5,10,15,20,25))+
  scale_y_continuous(breaks=c(0,0.017,0.05,0.1))
ggsave("LE/outs/5B_EARallAges.pdf",width=2.5,height=3)

####### losses below are added to the figure at powerpoint level
vital_vars(us_mort)   #it knows from here to use Exposures to form Mortality rate
us_mort|>filter(Sex == "Total", Year == 2024)|>life_table()|>filter(Age==59) # 24.5 years
us_mort|>filter(Sex == "Total", Year == 2024)|>mutate(Mortality=1.8*Mortality)|>life_table()|>filter(Age==59) # 19.5 years
us_mort|>filter(Sex == "Total", Year == 2024)|>mutate(Mortality=0.017+Mortality)|>life_table()|>filter(Age==59) # 19.4 years
# thus the steady state losses are 5 and 5.1 years for A and B, respectively.

### now get the losses of the transient components
(Vit=us_mort |>filter(Sex == "Total", Year == 2024))
Vit=Vit|>filter(Age>58) # save with capital V for both panels
vit=Vit
(lenV=length(vit$Mortality)) #52 long
(lenD=length(D$EAR)) #22 long, so add 30 zeros
vit$Mortality=vit$Mortality*c(D$RR,rep(1.8,30))
vit|>life_table()|>filter(Age==59) # 17.5 years so 6.9 are missing! 5 from ss and 1.9 from transients
vit=Vit
vit$Mortality=vit$Mortality+c(D$EAR,rep(0.017,30))
vit|>life_table()|>filter(Age==59) # 16.3 years so 8.2 are missing! 5.1 from ss and 3.1 from transients

###########  C and D
d59=d20|>filter(agedx>=54,agedx<=64) #9k
(D=SEERaBomb::msd(d59,mrt,brkst=c(0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21)))
D=D|>rename(Group="sex")|>select(Group,int,everything())
(D=SEERaBomb::foldD(D,keep=c("int")))

D|>ggplot(aes(x=t,y=RR)) +gp+gl+gx+gh1+gh1p8+tc(13)+leg+scale_y_continuous(breaks=c(1,1.8,5,10,15))+geRR+
  scale_color_manual(values = c("gray","black"),name = "Year of Dx")+ylab("Relative Risk of Death")
ggsave("LE/outs/5C_RR54to64.pdf",width=2.5,height=3)

D|>ggplot(aes(x=t,y=EAR))+gp+gl+gx+gyE+gh0+ghp017+tc(13)+leg+geE+
  scale_x_continuous(breaks=c(0,5,10,15,20,25))+
  scale_y_continuous(breaks=c(0,0.017,0.05,0.1))
ggsave("LE/outs/5D_EAR54to64.pdf",width=2.5,height=3)

####### losses below are added to the figure at powerpoint level

(Vit=us_mort |>filter(Sex == "Total", Year == 2024))
Vit=Vit|>filter(Age>58)
vit=Vit  ## 5C
vit=vit|>mutate(Mortality=1.8*Mortality)
vit|>life_table()|>filter(Age==59) # 19.5 years   5 lost to steady state, as in 5A
(lenV=length(vit$Mortality)) #52 long
(lenD=length(D$EAR)) #22 long, so add 30 zero
vit=Vit
vit$Mortality=vit$Mortality*c(D$RR,rep(1.8,30))
vit|>life_table()|>filter(Age==59) # 17.0 years so 7.5 are missing! 5 from ss and 2.5 from transients
vit=Vit  #5D
vit=vit|>mutate(Mortality=0.017+Mortality)
vit|>life_table()|>filter(Age==59) # 19.4 years  5.1 lost to steady state, as in 5B
vit=Vit  #reset vit to Vit
vit$Mortality=vit$Mortality+c(D$EAR,0.05+0.03/5*(1:30))
vit|>life_table()|>filter(Age==59) # 16.3 years so 8.2 are missing! 5.1 from ss and 3.1 from transients
vit=Vit
vit$Mortality=vit$Mortality+c(D$EAR,c(0.05+(0.05/9)*1:9,rep(0.10,21)))
vit|>life_table()|>filter(Age==59) # same answer if cap at 0.10 at ages over 90


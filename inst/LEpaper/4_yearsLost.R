# 4_yearsLost.R   (This script produces Figures 4A, 4B and 4C of an upcoming CML LE paper)
graphics.off();rm(list=ls()) 
library(tidyverse)
tc=function(sz) theme_classic(base_size=sz)
load("~/data/CMLepi/cml12.RData") 
head(d<-d12%>%filter(histo3%in%c(9863,9875))) #15325 
table(d$histo3,d$yrdx) #see new code really only starts in 2001 
d=d|>filter(histo3==9863|(histo3=9875)&(yrdx>2000)) 
(d=d|>mutate(histo3=as_factor(histo3)))
D=d|>group_by(yrdx,histo3)|>summarize(mage=mean(agedx))#both 
gh=geom_hline(yintercept=c(50,59),col="gray")
gv=geom_vline(xintercept=c(1992,2007),col="gray")
myt=theme(
  legend.key.size = unit(0.4, "cm"),
  legend.position = c(0.20,0.4),
  legend.text = element_text(size = 8),
  legend.title = element_text(size = 8, margin = margin(b = 1, unit = "pt")),
  legend.key.spacing.y = unit(0.0, "cm"),
  legend.margin = margin(t = 0, r = 0, b = 0, l = 0, unit = "pt")
)
D|>ggplot(aes(x=yrdx,y=mage,group=histo3,col=histo3)) + gh + gv + geom_line()+tc(13)+
  labs(y="Mean Age at Dx",x="Year of Dx",col="ICD-O-3")+
  scale_x_continuous(breaks=c(1992,2000,2007,2023))+
  scale_y_continuous(breaks=c(50,53,56,59))+myt
ggsave("LE/outs/4A_dxAges.pdf",width=2.5,height=2) 
ggsave("LE/outs/4A_dxAges.png",width=2.5,height=2) 

# get life expectancy at 59 as function of years
# https://github.com/robjhyndman/vital  updated for US HMD data
library(vital)  
library(tsibble)
files=c("Deaths_1x1.txt", "Exposures_1x1.txt", "Population.txt", "Mx_1x1.txt")
files=paste0("~/data/hmd_countries/USA/stats/",files) #put files above from HMD into here
us_mort=read_hmd_files(files) # function in R package vital
# save(us_mort,file="~/data/mrt/us_mort.RData")
# load("~/data/mrt/us_mort.RData")
(D=us_mort |>filter(Sex == "Total", Year >= 1992)|>life_table()|>filter(Age==59))
gh=geom_hline(yintercept=24.5,col="gray")
gv=geom_vline(xintercept=c(2019,2024),col="gray")
D|>ggplot(aes(x=Year,y=ex))+gh+gv+ 
  geom_line()+tc(13)+labs(y="LE at Age 59",x="Year") + 
  scale_x_continuous(breaks=c(1992,2000,2010,2019,2024))+
  theme(axis.text.x = element_text(size = 9) )
ggsave("LE/outs/4B_normalLEs.pdf",width=3,height=2)
ggsave("LE/outs/4B_normalLEs.png",width=3,height=2)

head(d<-d12%>%filter(histo3%in%c(9863,9875))) #15325 
d=d|>mutate(surv=ifelse(surv>50,0,surv)) 
table(d$histo3,d$yrdx)
d=d|>filter(histo3==9863|(histo3=9875)&(yrdx>2000)) # >50 deaths by 9875 after year 2000
(Dt=d|>mutate(ageDth=surv+agedx)) 
(Dt=Dt|>mutate(yrDth=floor(surv+yrdx),histo3=as_factor(histo3)))
(Dt=Dt|>filter(status==1))
(DT=Dt|>group_by(yrDth,histo3)|>summarize(mage=mean(ageDth)))
gh=geom_hline(yintercept=c(65,75),col="gray")
DT|>ggplot(aes(x=yrDth,y=mage,group=histo3,col=histo3)) + gh + geom_line()+tc(13)+
  labs(y="Mean Age at Death",x="Year of Death",col="ICD-O-3") + 
  scale_y_continuous(breaks=seq(65,75,5))+geom_smooth()+
  scale_x_continuous(breaks=c(1975,1986,1992,2000,2010,2023))+ myt +
theme(legend.position = c(0.85,0.18))
ggsave("LE/outs/4C_deathsAges.pdf",width=2.5,height=2)
ggsave("LE/outs/4C_deathsAges.png",width=2.5,height=2)


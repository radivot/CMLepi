## 6_dxDecades.R
graphics.off();rm(list=ls())#clear plots and environment 
library(tidyverse)
library(vital)  
library(tsibble)
load("~/data/mrt/us_mort.RData") #us_mort is a class vital object (single tibble-like)  
load("~/data/mrt/mrtUSA.RData")#mrt is list of 3 matrices
load("~/data/CMLepi/cml20.RData") #made in mkSEER.R
(d20=d20|>filter(histo3%in%c(9863,9875),agedx<90,surv<80,surv>0)) #43932
d20=d20|>select(yrdx,agedx,sex,surv,status) # 43,932 CML cases

dt=d20|>filter(agedx>=80,agedx<90) #4.7k
(D=SEERaBomb::msd(dt,mrt,brkst=c(0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15)))
D=D|>rename(Group="sex")|>select(Group,int,everything())
(D=SEERaBomb::foldD(D,keep=c("int")))
(dL=tibble(t=5:25,EAR=c(rep(0.1,21))))

gx=xlab("Years Since Dx")
gl=geom_line()
gp=geom_point()
gyE=ylab("Excess Absolute Risk of Death")
tc=function(sz) theme_classic(base_size=sz)
geE=geom_errorbar(aes(ymin=LL,ymax=UL),width=0.2)#for absolute risks
gh0=geom_hline(yintercept=0)
ccEAR=coord_cartesian(ylim=c(0,0.4))
txt="Age at Dx in 80-89"
ghp1=geom_hline(yintercept=0.1,col="gray")

D|>ggplot(aes(x=t,y=EAR))+gp+gl+gx+gyE+gh0+ghp1+tc(13)+geE+ccEAR+
  geom_line(data=dL,col="red",linewidth=1)+
  ggtitle(txt)+  theme(plot.title = element_text(size = 10))
ggsave("LE/outs/6A_85.pdf",width=2.5,height=3)
(Vit=us_mort|>filter(Sex == "Total", Year == 2024))
Vit=Vit|>filter(Age>84)
(n85=Vit|>life_table()|>filter(Age==85)) # 6.91 years
vit=Vit|>mutate(Mortality=0.1+Mortality)
vit|>life_table()|>filter(Age==85) #4.5 years 
dD=D|>filter(t<5)
(lenV=length(vit$Mortality)) #26 long
(lenD=length(dD$EAR)) #5 long, so add 21
vit=Vit  #reset Vit
vit$Mortality=vit$Mortality+c(dD$EAR,rep(0.1,21)) 
(v85=vit|>life_table()|>filter(Age==85)) # 3.32 years so 6.91-3.32 = 3.6 are missing, 2.4 ss, 1.2 transients

txt="Age at Dx in 70-79"
ccEAR=coord_cartesian(ylim=c(0,0.2))
ghp06=geom_hline(yintercept=0.06,col="gray")

dt=d20|>filter(agedx>=70,agedx<80) #7.7k
(D=SEERaBomb::msd(dt,mrt,brkst=c(0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21)))
D=D|>rename(Group="sex")|>select(Group,int,everything())
(D=SEERaBomb::foldD(D,keep=c("int")))
(dL=tibble(t=15:35,EAR=c(rep(0.1,21))))
D|>ggplot(aes(x=t,y=EAR))+gp+gl+gx+gyE+gh0+ghp06+tc(13)+ccEAR+geE+
  geom_line(data=dL,col="red",linewidth=1)+
  ggtitle(txt)+  theme(plot.title = element_text(size = 10))
ggsave("LE/outs/6A_75.pdf",width=3.5,height=3)
(Vit=us_mort|>filter(Sex == "Total", Year == 2024))
Vit=Vit|>filter(Age>74)
(n75=Vit|>life_table()|>filter(Age==75)) # 12.6 years
vit=Vit|>mutate(Mortality=0.06+Mortality)
vit|>life_table()|>filter(Age==75) #8.2 years =>lost 4.2
dD=D|>filter(t<15)
(lenV=length(vit$Mortality)) #36 long
(lenD=length(dD$EAR)) #15 long, so add 21
vit=Vit  #reset Vit
vit$Mortality=vit$Mortality+c(dD$EAR,rep(0.1,21)) 
(v75=vit|>life_table()|>filter(Age==75)) # 7 years so 12.6-7 = 5.6 are missing 

txt="Age at Dx in 60-69"
ccEAR=coord_cartesian(ylim=c(0,0.1))
ghp025=geom_hline(yintercept=0.025,col="gray")
dt=d20|>filter(agedx>=60,agedx<70) #7.7k
(D=SEERaBomb::msd(dt,mrt,brkst=c(0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21)))
D=D|>rename(Group="sex")|>select(Group,int,everything())
(D=SEERaBomb::foldD(D,keep=c("int")))
(dL=tibble(t=15:45,EAR=c(0.05+(0.05/10)*(0:10),rep(0.1,20))))
D%>%ggplot(aes(x=t,y=EAR))+gp+gl+gx+gyE+gh0+ghp025+tc(13)+geE+ccEAR+
  geom_line(data=dL,col="red",linewidth=1)+
  ggtitle(txt)+  theme(plot.title = element_text(size = 10))
ggsave("LE/outs/6A_65.pdf",width=4.5,height=3)
(Vit=us_mort|>filter(Sex == "Total", Year == 2024))
Vit=Vit|>filter(Age>64)
(n65=Vit|>life_table()|>filter(Age==65)) # 19.8 years
vit=Vit|>mutate(Mortality=0.025+Mortality)
vit|>life_table()|>filter(Age==65) #15 years =>lost 4.8 by min
dD=D|>filter(t<15)
(lenV=length(vit$Mortality)) #46 long
(lenD=length(dD$EAR)) #15 long, so add 31
vit=Vit  #reset Vit
vit$Mortality=vit$Mortality+c(dD$EAR,dL$EAR) 
(v65=vit|>life_table()|>filter(Age==65)) # 12.6 years so 19.8-12.6 = 7.2 are missing 

txt="Age at Dx in 50-59"
dt=d20|>filter(agedx>=50,agedx<60) #8.2k
ghp012=geom_hline(yintercept=0.012,col="gray")
(D=SEERaBomb::msd(dt,mrt,brkst=c(0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21)))
D=D|>rename(Group="sex")|>select(Group,int,everything())
(D=SEERaBomb::foldD(D,keep=c("int")))
(dL=tibble(t=15:55,EAR=c(0.015+(0.035/10)*(0:10),0.05+(0.05/10)*(1:10),rep(0.1,20))))
D|>ggplot(aes(x=t,y=EAR))+gp+gl+gx+gyE+gh0+ghp012+tc(13)+geE+ccEAR+
  geom_line(data=dL,col="red",linewidth=1)+
  scale_y_continuous(breaks=c(0,0.012,0.025,0.05,0.075,0.1)) +
  ggtitle(txt)+  theme(plot.title = element_text(size = 10))
ggsave("LE/outs/6A_55.pdf",width=5.5,height=3)
(Vit=us_mort|>filter(Sex == "Total", Year == 2024))
Vit=Vit|>filter(Age>54)
(n55=Vit|>life_table()|>filter(Age==55)) # 27.8 years
vit=Vit|>mutate(Mortality=0.015+Mortality)
vit|>life_table()|>filter(Age==55) #22.1 years =>lost 5.7 by min
dD=D|>filter(t<15)
length(vit$Mortality) #46 long
length(dD$EAR) #15 long, so add 31
vit=Vit  #reset Vit
vit$Mortality=vit$Mortality+c(dD$EAR,dL$EAR) 
(LT55=vit|>life_table())
(v55=LT55|>filter(Age==55)) # 19.3 years so 27.8-19.3 = 8.5 are missing 

txt="Age at Dx in 40-49"
dt=d20|>filter(agedx>=40,agedx<50) #8.2k
ghp01=geom_hline(yintercept=0.01,col="gray")
(D=SEERaBomb::msd(dt,mrt,brkst=c(0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21)))
D=D|>rename(Group="sex")|>select(Group,int,everything())
(D=SEERaBomb::foldD(D,keep=c("int")))
(dL=tibble(t=15:65,EAR=c(0.01+(0.005/10)*(0:10),0.015+(0.035/10)*(1:10),0.05+(0.05/10)*(1:10),rep(0.1,20))))
D|>ggplot(aes(x=t,y=EAR))+gp+gl+gx+gyE+gh0+ghp01+tc(13)+geE+ccEAR+
  geom_line(data=dL,col="red",linewidth=1)+
  scale_y_continuous(breaks=c(0,0.01,0.025,0.05,0.075,0.1)) +
  ggtitle(txt)+  theme(plot.title = element_text(size = 10))
ggsave("LE/outs/6A_45.pdf",width=6.5,height=3)
(Vit=us_mort|>filter(Sex == "Total", Year == 2024))
Vit=Vit|>filter(Age>44)
(n45=Vit|>life_table()|>filter(Age==45)) # 36.6 years
vit=Vit|>mutate(Mortality=0.01+Mortality)
vit|>life_table()|>filter(Age==45) #30.1 years =>lost 5.5 if nadir held flat
dD=D|>filter(t<15)
length(vit$Mortality) #66 long
length(dD$EAR) #15 long, so add 31
vit=Vit  #reset Vit
vit$Mortality=vit$Mortality+c(dD$EAR,dL$EAR) 
(v45=vit|>life_table()|>filter(Age==45)) # 25.9 years so 36.6-25.9 = 10.7 are missing 
paste0("Normal LE=",n45["ex"],", CML LE=",v45["ex"],", Loss=",n45["ex"]-v45["ex"])
(v=bind_rows(v45,v55,v65,v75,v85))
(n=bind_rows(n45,n55,n65,n75,n85))
(dLE=tibble(LE0=n$ex,LE=v$ex,LLE=LE0-LE))
dLE$Agedx=c(45,55,65,75,85)
dLE$name="LLE"
dLE$Sex="Both"
dLE$paper="this study"
(dLE=dLE|>mutate(value=LLE,PLLE=LLE/LE0)) #our LLE and PLLE in B and C 
#     LE0    LE   LLE Agedx name  Sex   paper      value  PLLE
# 1 36.6  25.9  10.7     45 LLE   Both  this study 10.7  0.292
# 2 27.8  19.3   8.54    55 LLE   Both  this study  8.54 0.307
# 3 19.8  12.6   7.20    65 LLE   Both  this study  7.20 0.364
# 4 12.6   6.98  5.67    75 LLE   Both  this study  5.67 0.448
# 5  6.91  3.32  3.59    85 LLE   Both  this study  3.59 0.519

# next get values of BJH paper 
# install.packages("docxtractr")
library(docxtractr)
doc <- read_docx("~/pdfs/lifeExp/cml/cmlLifeExpSEER_BJH22sup.docx") #from Shivarov & Grigorova online
class(doc)
tbls <- docx_extract_all_tbls(doc)
(b=tbls[[3]])
(b=b|>filter(Period=="2010s")|>select(Agedx=Specific.age,Race,Sex,bLE=LE))
(b=b|>filter(Race=="White")|>select(-Race,Sex)|>mutate(Agedx=as.numeric(Agedx),bLE=as.numeric(bLE)))
# in BJH paper LE are OK, but losses relative to a fixed old age need to be relative to background LE 
(Vit=us_mort |>filter(Sex%in%c("Male","Female"), Year == 2024)) #which we take here from 2024
(LT=Vit|>life_table()|>filter(Age%in%c(55,65,75,85)))
(bD=LT|>select(Sex,Agedx=Age,LE0=ex)|>arrange(Agedx)|>as_tibble())
(bjh=left_join(bD,b)|>select(-Year))
(bjh=bjh|>mutate(bLLE=LE0-bLE))
(bjh=bjh|>select(Sex,Agedx,LE0,LE=bLE,LLE=bLLE)|>mutate(PLLE=LLE/LE0))
bjh$paper="Shivarov-Grigorova"
bjh
#   Sex    Agedx   LE0    LE   LLE  PLLE paper             
#   <chr>  <dbl> <dbl> <dbl> <dbl> <dbl> <chr>             
# 1 Female    55 29.3  17.6  11.7  0.400 Shivarov-Grigorova
# 2 Male      55 26.2  15.0  11.2  0.426 Shivarov-Grigorova
# 3 Female    65 20.9  11.9   9.05 0.432 Shivarov-Grigorova
# 4 Male      65 18.5   9.67  8.84 0.478 Shivarov-Grigorova
# 5 Female    75 13.4   6.47  6.89 0.516 Shivarov-Grigorova
# 6 Male      75 11.8   6.00  5.78 0.491 Shivarov-Grigorova
# 7 Female    85  7.24  3.70  3.54 0.489 Shivarov-Grigorova
# 8 Male      85  6.41  2.75  3.66 0.571 Shivarov-Grigorova

(en13=tibble(Agedx=rep(seq(45,85,10),2),Sex=rep(c("Male","Female"),each=5),name="LLE", paper="Bower to 2013",
             LE0=c(40.5,29.9,20.1,11.7,5.7,42.9,32.6,22.8,14.0,6.9),
             LLE=c(3,2.6,2.5,2.2,1.6,3.3,2.9,2.9,2.6,2.0),  
             LLE_L=c(0.8,1.0,1.2,1.2,0.9,0.9,1.2,1.4,1.4,1.2),  
             LLE_U=c(5.2,4.1,3.8,3.2,2.3,5.7,4.6,4.4,3.8,2.8),  
             PLLE=c(.07,.09,.13,.18,.28,.08,.09,.13,.19,.28),  
             PLLE_L=c(.02,.04,.06,.10,.16,.02,.04,.06,.10,.17),  
             PLLE_U=c(.13,.14,.20,.27,.40,.13,.14,.19,.27,.40)  
))

(en20=tibble(Agedx=rep(seq(46,86,10),2),Sex=rep(c("Male","Female"),each=5),name="LLE", paper="Bower to 2020",
             LE0=c(40.6,29.9,20.1,11.6,5.6,43.0,32.4,22.5,13.7,6.7),
             LLE=c(4.7,4.5,4.3,3.3,2.1,4.6,4.5,4.5,3.8,2.4),  
             LLE_L=c(2.3,3.2,3.4,2.8,1.7,2.1,3.2,3.6,3.1,2.0),  
             LLE_U=c(7.1,5.7,5.2,3.9,2.5,7.1,5.8,5.5,4.4,2.9),  
             PLLE=c(.11,.15,.21,.29,.38,.11,.14,.20,.27,.36),  
             PLLE_L=c(.06,.11,.17,.24,.31,.05,.10,.16,.22,.29),  
             PLLE_U=c(.17,.19,.26,.34,.44,.16,.18,.24,.32,.43)  
))


gv=geom_vline(xintercept=59,col="gray")
myt=theme(legend.key.size = unit(0.4, "cm"),
  legend.position = c(0.84,0.82),
  legend.text = element_text(size = 7),
  legend.title=element_blank(),
  legend.key.spacing.y = unit(0.0, "cm"),
  legend.margin = margin(t = 0, r = 0, b = 0, l = 0, unit = "pt") )

bjh|>ggplot(aes(x=Agedx,y=LLE,col=Sex,pch=paper)) + gv+ geom_line()+ geom_point()+
  tc(13)+xlab("Age at Dx")+ylab("Years of Life Lost")+gh0+
  scale_y_continuous(breaks=seq(0,12,2))+
  geom_point(data=en13,size=2)+
  geom_point(data=dLE,size=2)+
  geom_line(data=dLE)+ 
  geom_line(data=en13)+
  geom_errorbar(aes(ymin=LLE_L,ymax=LLE_U),data=en13,width=0.2)+
  geom_point(data=en20,size=2)+
  geom_line(data=en20)+
  geom_errorbar(aes(ymin=LLE_L,ymax=LLE_U),data=en20,width=0.2)+myt
ggsave("LE/outs/6B_LLE.pdf",width=4,height=4)

cc1=coord_cartesian(ylim=c(0.0,1.05))
dLE|>ggplot(aes(x=Agedx,y=PLLE,col=Sex,pch=paper))+gv +geom_line()+geom_point()+gh0+
  geom_point(data=en13,size=2)+cc1+
  geom_point(data=bjh,size=1.5)+
  geom_line(data=bjh)+
  geom_line(data=en13)+
  geom_errorbar(aes(ymin=PLLE_L,ymax=PLLE_U),data=en13,width=0.2)+
  geom_point(data=en20,size=2)+
  geom_line(data=en20)+
  geom_errorbar(aes(ymin=PLLE_L,ymax=PLLE_U),data=en20,width=0.2)+
  tc(13)+xlab("Age at Dx") +   ylab("Proportion of Life Lost") + myt
ggsave("LE/outs/6C_PLLE.pdf",width=4,height=4)


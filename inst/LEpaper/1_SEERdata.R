# 1_SEERdata.R  (makes Figure 1A and 1B of a CML LE paper in the works)
graphics.off();rm(list=ls()) 
library(tidyverse)
load("~/data/CMLepi/cml20.RData") 
head(d20<-d20%>%filter(histo3%in%c(9863,9875))) #45636
d20$db="SEER20"
D20=d20|>summarize(n=n(),.by = c(yrdx, histo3,db))|>mutate(histo3=factor(histo3))

load("~/data/CMLepi/cml12.RData") 
head(d12<-d12%>%filter(histo3%in%c(9863,9875))) #15325 
d12$db="SEER12"
D12=d12|>summarize(n=n(),.by = c(yrdx, histo3,db))|>mutate(histo3=factor(histo3))

load("~/data/CMLepi/cml8.RData") 
head(d8<-d8%>%filter(histo3%in%c(9863,9875))) #14919 
d8$db="SEER8"
D8=d8|>summarize(n=n(),.by = c(yrdx, histo3,db))|>mutate(histo3=factor(histo3))
D=bind_rows(D20,D12,D8)|>mutate(db=as_factor(db)) 
Dtot=D|>summarize(n=sum(n),.by = c(yrdx,db))
myt=theme(
  legend.key.size = unit(0.5, "cm"),
  legend.position = c(0.2,0.65),
  legend.text = element_text(size = 9),
  legend.title = element_text(size = 9, margin = margin(b = 1, unit = "pt")),
  legend.key.spacing.y = unit(0.0, "cm"),
  legend.margin = margin(t = 0, r = 0, b = 0, l = 0, unit = "pt")) 
Dtot|>ggplot(aes(x=yrdx,y=n,col=db))+geom_point(size=1)+
  labs(y="Number of Cases",x="Year of Diagnosis",col="Database")+
  ylim(c(0,NA))+theme_classic(base_size=14)+ myt
ggsave("LE/outs/1A_caseCounts.pdf",width=4,height=2.5) 
ggsave("LE/outs/1A_caseCounts.png",width=4,height=2.5) 

D|>ggplot(aes(x=yrdx,y=n,col=db,shape=histo3))+geom_point(size=1)+
  labs(y="Number of Cases",x="Year of Diagnosis",col="Database",shape="ICD-O-3")+
  scale_shape_manual(values = c(6,2))+
  ylim(c(0,NA))+theme_classic(base_size=14)+myt
ggsave("LE/outs/1B_countsByCodes.pdf",width=4,height=2.5) 
ggsave("LE/outs/1B_countsByCodes.png",width=4,height=2.5) 

d=bind_rows(d20,d12,d8) #75880
d=d|>distinct(pick(-db)) # don't count it as different via db being different
d # 53254 unique cases = number in figure legend; others there were read off by eye

hist(d$surv) #outlier cluster past 80 are unknown survival times
(du=d|>filter(surv>80)) # 668 with unknown survival times = 89.7
table(du$status) # 666 of 668 are dead
(d0=d|>filter(surv==0)) #261 likely left registry area right after Dx
table(d$agedx)

dOld=d|>filter(agedx>=90) #1295
table(dOld$status)# 91 alive, 1204 dead
table(dOld$COD2) #635 by LC, 569 by OC
dOldOC=dOld|>filter(COD2=="OC")
sort(table(dOldOC$CODS)) # heart=223, cerebroVasc=25, athero=11 => 223+25+11=259 in text
table(dOld$COD7) #263 by CVD (3 more via hypertension + 1 via aortic aneurism)
table(dOld$status,dOld$yrdx) # tells us to take it up thru 2017 to have most cases dead
dOld=dOld|>filter(surv<80) # 1159 => lost 136 to survival unknown
dOld=dOld|>filter(surv>0) # 1135  => lost 24 more to survival = 0
dOld|>group_by(yrdx)|>summarize(mn=mean(surv))|>t() # and 2017 is where LE peaks before censoring brings it back down
dOld=dOld|>filter(yrdx<=2017)
table(dOld$status) # 870 dead, 8 still alive
(dOld=dOld|>mutate(yrG=cut(yrdx,breaks=c(1975,1990,2000,2005,2011,2017),include.lowest=T,dig.lab=4)))
dOld|>group_by(yrG)|>summarize(mn=mean(surv),sd=sd(surv)) 
summary(lm(surv~yrdx,data=dOld))
summary(lmG<-lm(surv~0+yrG,data=dOld))
(ci=round(cbind(coef(lmG),confint(lmG)),2))
paste0(ci[,1]," (",ci[,2],", ",ci[,3],")",collapse=", ")
# "1.02 (0.67, 1.38), 0.92 (0.6, 1.24), 1 (0.73, 1.27), 1.32 (1.08, 1.57), 1.81 (1.58, 2.03)"
1/3.8 #26%
1.81/4.0 #45%
table(d20$COD2) # 8194 by leukemia
table(d20$CODS) 
table(d20$CODS=="Chronic Myeloid Leukemia")       # 5210 by CML
table(d20$CODS=="Aleukemic, Subleukemic and NOS") # 1117 by sub-Leu
table(d20$CODS=="Acute Myeloid Leukemia")          # 772 by AML
table(d20$CODS=="Acute Monocytic Leukemia")        # 8 AMoL => AML=780
table(d20$CODS=="Chronic Lymphocytic Leukemia")   # 220 by CLL
table(d20$CODS=="Acute Lymphocytic Leukemia")     # 134 by ALL
table(d20$CODS=="Other Myeloid/Monocytic Leukemia") # 488 by OML
table(d20$CODS=="Other Acute Leukemia")          # 226 by OAL
table(d20$CODS=="Other Lymphocytic Leukemia") # 19 OLL
1117+772+8+220+134+488+226+19 # 2984
488+226+19 #733 by OL
1117+780+220+134+733#2984 = sum of subleu=1117,AML=780,CLL=220,ALL=134,OL=733
8194*2.4# 19665.6 US CML patient deaths by a leukemia in 2000-2023 among those Diagnosed in 2000-2023
# to put this in perspective, 2984 deaths by other leukemias is way too high relative to  
table(d20$CODS=="Lung and Bronchus")# 351 lung cancer deaths

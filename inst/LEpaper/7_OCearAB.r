# 7_OCearAB.R
graphics.off();rm(list=ls())
library(tidyverse)
load("~/data/mrt/mrtUSA.RData")#mrt is list of 3 matrices
F=reshape2:::melt(mrt[[1]],value.name="rate")
M=reshape2:::melt(mrt[[2]],value.name="rate")
F$sex="Female"
M$sex="Male"
bh<-bind_rows(F,M)|>tibble()
names(bh)[1:2]=c("age","year")
bh=bh|>filter(age!="110+")
bh=bh|>mutate(age=as.numeric(as.character(age))+0.5)
bh  #background hazard as a long tibble

load("~/data/CMLepi/cml20.RData") #made in mkSEER.R
(d20=d20|>filter(histo3%in%c(9863,9875),agedx<90,surv<80,surv>0)) #43932

yrs=1975:2023; ages=0.5:125.5
Z=matrix(0,ncol=length(yrs),nrow=length(ages))
colnames(Z)=yrs
rownames(Z)=ages
head(Z)# all zeros initially
brks=c(0,1,2,3,4,5,6,7,8,10,12,15,18)
(binS<-levels(cut(brks+0.1,breaks=c(brks,100))))
L=NULL  # initiate big list L
L$brks=brks
L$binS=binS
L$d=d20
L

for (i in c("LC","OC"))
  for (bin in binS) 
    for (j in c("Male","Female")) { 
      # j="Male"; bin="(0,0.5]"; i="OC"
      (binIndx=SEERaBomb::getBinInfo(bin,binS)["index"])
      (bin=binS[binIndx])
      (LL=SEERaBomb::getBinInfo(bin,binS)["LL"])
      head(d<-L$d%>%filter(sex==j))
      d=d|>mutate(survC=cut(surv,breaks=c(-1,brks,100),include.lowest = TRUE)) 
      d=d|>filter(surv>LL) 
      d=d|>mutate(py=SEERaBomb::getPY(surv,bin,binS,brks)) # getpy leaves zeros when surv end is left of LL
      d=d|>filter(py>0)  #get rid of such rows upfront
      d=d|>mutate(ageL=agedx+LL) 
      d$year=floor(d$yrdx+LL)
      PYin=d|>select(py,ageL,year)
      if(dim(PYin)[1]==0) binMidPnt=LL else binMidPnt=LL+sum(PYin$py)/dim(PYin)[1]/2
      PYin=as.matrix(PYin)
      PY=Z+0 # add 0 to make sure not a lazy copy of just the pointer 
      SEERaBomb::fillPYM(PYin, PY)
      dO=d|>filter(survC==bin) 
      (nDead=sum(dO$COD2==i)) 
      N=1e4 ## shrink by N big then multiply by it later, so each PY of 1 lands in one bin
      head(Od<-dO%>%mutate(py=(COD2==i)/N)|>select(py,ageL,year))
      head(Oin<-as.matrix(Od))
      O=Z+0
      SEERaBomb::fillPYM(Oin, O)
      O=O*N
      head(dO<-reshape2:::melt(O,value.name="Obs"))
      head(dP<-reshape2:::melt(PY,value.name="PY"))
      head(dd<-left_join(dO,dP))
      names(dd)[1:2]<-c("age","year")
      (D=dd|>filter(PY>0))
      D$sex=j
      D$cause=i
      (D=D|>tibble())
      D=left_join(D,bh)
      if (i=="LC") D$Eus=0  
      if (i=="OC")  D=D|>mutate(Eus=PY*rate)
      D=D|>rename(num=Obs,denom=PY)
      D=D|>filter(age<90,age>20) 
      L[[i]][[bin]][[j]]$mid=binMidPnt
      L[[i]][[bin]][[j]]$D=D
    }
d

str(L[["OC"]])
str(L[["LC"]])
getMid=function(x) mean(x$Male$mid,x$Female$mid)
sapply(L[["LC"]],getMid)
(mids=sapply(L[["LC"]],getMid)) # same PY at risk for each death type
poolMF=function(x) bind_rows(x$Male$D,x$Female$D)
(Doc=lapply(L[["OC"]],poolMF)) 
(Dlc=lapply(L[["LC"]],poolMF)) 
sumM=function(x) x|>summarize(O=sum(num),PY=sum(denom),E=sum(Eus),cause=toupper(cause[1]))
(D=tibble(t=mids,int=binS,lc=Dlc,oc=Doc))
D=D|>mutate(lc1=map(lc,sumM))
D=D|>mutate(oc1=map(oc,sumM))
D=D|>select(-lc,-oc)
Dl=D|>select(-oc1)
Do=D|>select(-lc1)
(Dl=Dl|>unnest(lc1))
(Do=Do|>unnest(oc1))
D=bind_rows(Do,Dl)
(D=D|>mutate(EAR=(O-E)/PY,LL=EAR-1.96*sqrt(O)/PY,UL=EAR+1.96*sqrt(O)/PY))
D=D|>mutate(cause=if_else(cause=="OC","Other Cause","Leukemic Cause"))
D=D|>mutate(cause=as.factor(cause)) 

leg=theme(legend.margin=margin(0,0,0,0),legend.title=element_blank(),
          legend.position="top") #redfine leg without size increase
EARbrks=c(0,0.01,0.05,0.1)
ccE=coord_cartesian(ylim=c(0,0.08))
tc=function(sz) theme_classic(base_size=sz)
ghp01=geom_hline(yintercept=c(0.01),col="gray")
gh0=geom_hline(yintercept=0)
gx=xlab("Years Since CML Diagnosis")
gyE=ylab("Excess Absolute Risk of Death")
geE=geom_errorbar(aes(ymin=LL,ymax=UL),width=0.2)#for absolute risks
gl=geom_line()

D|>ggplot(aes(x=t,y=EAR,col = cause))+ghp01+geE+geom_point(size=0.7)+gyE +gl+ccE+
  scale_y_continuous(minor_breaks=NULL,breaks=EARbrks) +tc(12)+gh0+leg+gx+
  scale_color_manual(values = c("gray","black"),name = "cause")
ggsave(file=paste0("LE/outs/7A_OCtime.pdf"),height=3,width=3.2) 

##### now do B
d=d20|>mutate(ageC=cut(agedx,c(0,seq(20,90,10)),include.lowest = TRUE, right = FALSE))
(L=split.data.frame(d,d$ageC))
Lf=lapply(L,function(x) SEERaBomb::msd(x,mrt,brkst=c(0,1))|>SEERaBomb::foldD())
(Db=bind_rows(Lf,.id="ageGroup")|>filter(int=="(1,100]"))
(D1=lapply(L,function(x) x|>filter(surv>1)))
Db$OC=sapply(D1,function(x) sum(x$COD2=="OC"))
(Db$age=c(10,seq(25,85,by=10)))
(Db1=Db|>mutate(EARoc = (OC - E)/PY, LLoc = EARoc - 1.96 * sqrt(OC)/PY, ULoc = EARoc + 1.96 * sqrt(OC)/PY))
EARbrks=c(0,0.003,0.01,0.04)
ccE=coord_cartesian(ylim=c(0,0.085))
ghp01=geom_hline(yintercept=c(0,0.003,0.01,0.04),col="gray")
geEoc=geom_errorbar(aes(ymin=LLoc,ymax=ULoc),width=0.2)#for absolute risks
Db1|>ggplot(aes(x=age,y=EARoc))+ghp01+geEoc+geom_point(size=0.7)+gyE+gl+#ccE+
  scale_y_continuous(minor_breaks=NULL,breaks=EARbrks) +tc(12)+
  gh0+leg+xlab("Decade of Age at Dx")+
  ggtitle("Excess OC Risks at >1y after Dx")+
  theme(plot.title = element_text(size = 10, # Set the font size
                                  hjust = 0.5)) # Center the title
ggsave(file="LE/outs/7B_OCage.pdf",height=3,width=3.2) 

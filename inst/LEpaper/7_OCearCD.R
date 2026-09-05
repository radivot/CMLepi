## 7_OCeearCD.R
graphics.off();rm(list=ls()) 
library(tidyverse)
library(mgcv)
load("~/data/CMLepi/cml20.RData") #made in mkSEER.R
(d20=d20|>filter(histo3%in%c(9863,9875),agedx<90,surv<80,surv>0)) #43932
yrs=1975:2023
ages=0.5:125.5
Z=matrix(0,ncol=length(yrs),nrow=length(ages))
colnames(Z)=yrs
rownames(Z)=ages
head(Z)# all zeros initially
survCut=0
brks=c(0,0.5,1,2,3,4,5,6,8,10,12,14) 
print(binS<-levels(cut(brks+0.1,breaks=c(brks,100))))  
load("~/data/CMLepi/G6.RData") # made by mkMorts.R
L=NULL
L$G=G
L$brks=brks
L$binS=binS
L$d=d20
L
(nms=names(L$G))

for (i in nms)
  for (bin in binS) 
    for (j in c("Male","Female")) 
    {
      # j="Male"; bin="(0,0.5]"; i="LC"
      (binIndx=SEERaBomb::getBinInfo(bin,binS)["index"])
      (bin=binS[binIndx])
      (LL=SEERaBomb::getBinInfo(bin,binS)["LL"])
      head(d<-L$d%>%filter(sex==j),2)
      d=d|>mutate(survC=cut(surv,breaks=c(-1,brks,100),include.lowest = TRUE)) 
      d=d|>filter(surv>LL) 
      d=d|>mutate(py=SEERaBomb::getPY(surv,bin,binS,brks)) # getpy leaves zeros when surv end is left of LL
      d=d|>filter(py>0)  #get rid of such rows upfront
      d=d|>mutate(ageL=agedx+LL) 
      d$year=floor(d$yrdx+LL)
      PYin=d%>%select(py,ageL,year)
      if(dim(PYin)[1]==0) binMidPnt=LL else binMidPnt=LL+sum(PYin$py)/dim(PYin)[1]/2
      PYin=as.matrix(PYin)
      PY=Z+0 # add 0 to make sure not a lazy copy of just the pointer 
      SEERaBomb::fillPYM(PYin, PY)
      dO=d|>filter(survC==bin) 
      N=1e4 ## shrink by N big then multiply by it later, so each PY of 1 lands in one bin
      head(Od<-dO%>%mutate(py=(COD7==i)/N)|>select(py,ageL,year),2)
      head(Oin<-as.matrix(Od),2)
      O=Z+0
      SEERaBomb::fillPYM(Oin, O)
      O=O*N
      head(dO<-reshape2:::melt(O,value.name="Obs"),2)
      head(dP<-reshape2:::melt(PY,value.name="PY"),2)
      head(dd<-left_join(dO,dP),2)
      names(dd)[1:2]<-c("age","year")
      D=dd|>filter(PY>0)
      D$sex=j
      D=D|>rename(num=Obs,denom=PY)
      D$Eus=as.numeric(exp(predict(L$G[[i]],D)))
      D$cause=i
      tibble(D)
      sum(D$num)
      D=D|>filter(age<90,age>20)
      L[[i]][[bin]][[j]]$mid=binMidPnt
      L[[i]][[bin]][[j]]$D=D
    }
names(L)

getMid=function(x) mean(x$Male$mid,x$Female$mid)
(mids=sapply(L[["YOC"]],getMid)) # same PY at risk for each death type
poolMF=function(x) bind_rows(x$Male$D,x$Female$D)
D=vector(length=length(nms),mode="list")
names(D)<-nms
for (i in nms) D[[i]]=lapply(L[[i]],poolMF)
sumM=function(x) x|>summarize(O=sum(num),PY=sum(denom),E=sum(Eus),cause=toupper(cause[1]))
(DD=tibble(t=mids,int=binS))
for (i in nms) DD[[i]]=lapply(D[[i]],sumM)  
Di=vector(length=length(nms),mode="list")
names(Di)<-nms
for (i in nms) Di[[i]]=(DD[,c("t","int",i)])|>unnest(all_of(i))
DD=bind_rows(Di)
(DD=DD|>mutate(EAR=(O-E)/PY,LL=EAR-1.96*sqrt(O)/PY,UL=EAR+1.96*sqrt(O)/PY))
(DD=DD|>filter(cause!="0"))

myt=theme(legend.text=element_text(size=12),strip.text=element_text(size=12))
geEAR=geom_errorbar(aes(ymin=LL,ymax=UL),width=.2,col="gray")
EARbrks=c(0,0.01,0.05,0.1)
ghE=geom_hline(yintercept=c(0.01),col="gray")
leg=theme(legend.margin=margin(0,0,0,0),legend.title=element_blank(),
          legend.position="top") #redfine leg without size increase
gyE=ylab("Excess Absolute Risk of Death")
geE=geom_errorbar(aes(ymin=LL,ymax=UL),width=0.2)#for absolute risks
jco=ggsci::scale_color_jco()
tc=function(sz) theme_classic(base_size=sz)
gh0=geom_hline(yintercept=0)
gx=xlab("Years Since CML Diagnosis")

DD|>filter(cause%in%c("CV","DK","ASH"))|>mutate(cause=as_factor(cause))|>
  ggplot(aes(x=t,y=EAR,col=cause,shape=cause))+ghE+geEAR+geom_point(size=2)+geom_line()+
  gyE+scale_y_continuous(minor_breaks=NULL,breaks=EARbrks)+ jco+tc(13)+gh0 + leg + gx
ggsave(file=paste0("LE/outs/7_OCearC_cvd.pdf"),height=3,width=3)

DD|>filter(cause%in%c("CA","IN","YOC"))|>
  ggplot(aes(x=t,y=EAR,col=cause,shape=cause))+ghE+geEAR+geom_point(size=2)+geom_line()+
  gyE+scale_y_continuous(minor_breaks=NULL,breaks=EARbrks)+ jco+tc(13)+gh0 + leg + gx
ggsave(file=paste0("LE/outs/7_OCearD_yoc.pdf"),height=3,width=3)


graphics.off();rm(list=ls())#clear plots and environment 
library(tidyverse)
tc=function(sz) theme_classic(base_size=sz)
load("~/data/CMLepi/cml20.RData") 
(d=d20|>filter(histo3%in%c(9863,9875),agedx<90,surv<80,surv>0)) #43932
agedxUL=64  
agedxLL=54  
agedxGap=10  
(d=d|>filter(agedx>=agedxLL,agedx<=agedxUL))
(d=d|>group_by(agedx,sex)|>summarize(n=n())|>ungroup())
(d=d|>mutate(p=n/sum(n)))
load("~/data/mrt/mrtUSA.RData")#loads US mortality data mrt thru 2024
myf=function(agedx,sex)  mrt[[sex]][(agedx+1):(agedx+45),"2024"]
(M=matrix(nrow=45,ncol=20))
for (i in 1:20) {
  M[,i]=d$p[i]*myf(d$agedx[i],d$sex[i])
}
M
(h59=apply(M,1,sum))
hb=data.frame(time=0:44,hazard=h59)
save(hb,file="~/data/mrt/hBack26.RData")

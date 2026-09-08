# bshazard.R  Inspired by Lab 3 in biostat3, this script shows how the shape of the hazard has been evolving for CML 
library(biostat3)
library(tidyverse) #get select() from here, not MASS above
library(bshazard) # smooth hazards
load("~/data/CMLepi/cml8.RData") #made in mkSEER.R  # 19250 cases
(d=d8|>filter(histo3%in%c(9863,9875),agedx<90,agedx>=80,surv<80,surv>0)) #1.7k
(d=d|>mutate(yrdxG=cut(yrdx,breaks=c(1975,1985,2000,2010,2023),include.lowest = T)))
(sf=survfit(Surv(surv, status==1) ~ yrdxG, data = d))
str(sf)
rainbow=c("red","orange","green","blue","violet")
par(mfrow=c(1,1))
sf|>plot(col = rainbow, xlab = "Years since Dx",ylab = "Survival",main = "K-M estimates")
legend("topright", levels(d$yrdxG), col=rainbow, lty = 1,bty="n")
par(mfrow=c(2,2))
for(level in levels(d$yrdxG)) 
  plot(bshazard(Surv(surv,status)~1,d|>filter(yrdxG==level),verbose=T), main=level, xlim=c(0,10),xlab="Years since Dx")
group_by(d,yrdxG)|>do(as.data.frame(bshazard(Surv(surv, status)~1, data=.,verbose=FALSE)))|>ungroup()|>
  ggplot(aes(x=time,y=hazard,group=yrdxG)) + geom_line(aes(col=yrdxG)) +
  geom_ribbon(aes(ymin=conf.low, ymax=conf.high, fill=yrdxG), alpha=0.3) +
  xlab('Years since Dx') + ylab('Hazard')
biostat3::survRate(Surv(surv,status)~yrdxG, data=d)
#                         yrdxG     tstop event      rate     lower     upper
# yrdxG=[1975,1985] [1975,1985]  564.2628   450 0.7975007 0.7255094 0.8747026
# yrdxG=(1985,2000] (1985,2000]  900.0219   521 0.5788748 0.5302288 0.6307840
# yrdxG=(2000,2010] (2000,2010] 1154.5243   373 0.3230768 0.2911179 0.3575862
# yrdxG=(2010,2023] (2010,2023] 1276.0027   303 0.2374603 0.2114728 0.2657598

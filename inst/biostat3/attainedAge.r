#attainedAge.r  Use attained age in survSplit. Inspired by Lab 8 in biostat3
# Also here and not in survSplit.R is use of formulas in survSplit.
graphics.off();rm(list=ls())#clear plots and environment 
library(biostat3)  # loads survival (for Surv and survSplit) 
library(bshazard)  # for bshazard
library(broom)     # tidy
library(tidyverse)     # tidy
load("~/data/CMLepi/cml.RData") #53.2k;  made in 1_mkSEERdata.R (unique cases of CML, with CMML cases removed)
d=d|>filter(agedx<90,agedx>=20) #50.7k  cases
d=d|>filter(surv>0) # 50.46k# clear these to avoid log(0) errors
(d=d|>mutate(yrdxG=cut(yrdx,breaks=c(1975,1985,2000,2010,2023),include.lowest = T),.after=yrdx))
levels(d$yrdxG) #"[1975,1985]" "(1985,2000]" "(2000,2010]" "(2010,2023]"
if (0) {#set to 1 to see interesting plots of hazards vs age in different eras, else code is slightly slow
  load("~/data/mrt/us_mort.RData")
  head(m<-us_mort |>filter(Sex == "Total", Year > 2000, Age>=20) )
  plot(bshazard(Surv(agedx, agedx+surv, status) ~ 1, verbose=FALSE,
                data=subset(d,yrdxG=="[1975,1985]")),ylim=c(0,0.6), xlim=c(20,85),conf.int=TRUE, xlab="Attained age (years)")
  lines(bshazard(Surv(agedx, agedx+surv, status) ~ 1,verbose=FALSE,
                 data=subset(d,yrdxG=="(1985,2000]")), col="red", conf.int=TRUE)
  lines(bshazard(Surv(agedx, agedx+surv, status) ~ 1,verbose=FALSE,
                 data=subset(d,yrdxG=="(2000,2010]")), col="orange", conf.int=FALSE)
  lines(bshazard(Surv(agedx, agedx+surv, status) ~ 1,verbose=FALSE,
                 data=subset(d,yrdxG=="(2010,2023]")), col="green", conf.int=FALSE)
  with(subset(m,Year==2005),lines(Age,Mortality, col="blue"))
  with(subset(m,Year==2017),lines(Age,Mortality, col="violet")) #similar recent progress with all causes of death
  legend("topleft",legend=c(levels(d$yrdxG),'Gen Pop 2005','Gen Pop 2017'),y.intersp=2,cex=1,col=c("black","red","orange","green","blue","violet"),lty=1,bty="n")
  dev.copy2pdf(file = "Rpacks/biostat3/outs/Q8attainedAge.pdf",width=5,height=7)#see very little recent progress and that 60 is a good cut-off age
}

summary(poi<-glm(status ~ yrdxG + offset(log(surv)), family=poisson, data=d))
# broom::tidy(poi, conf.int=TRUE, exponentiate=TRUE) # tidyverse (profile-based CIs)
biostat3::eform(poi) # Wald-based CIs

(d=d|>select(-(COD:CODS))) #trim off COD cols not used in this script
(d=d|>mutate(agedxG=cut(agedx, breaks=c(20,60,90),right = FALSE),.after=yrdx)|>select(-cancer))
# # A tibble: 50,462 × 9      # so this is the input into survSplit()
#        id sex    agedx  yrdx agedxG  yrdxG       histo3    surv status
#     <int> <fct>  <int> <dbl> <fct>   <fct>        <dbl>   <dbl>  <dbl>
#  1   5253 Female    71  2002 [60,90) (2000,2010]   9863 14.0         1
#  2   9874 Male      33  2002 [20,60) (2000,2010]   9875 21.7         0
#  3  98271 Male      87  2005 [60,90) (2000,2010]   9863  0.339       1
#  4 125911 Female    87  2009 [60,90) (2000,2010]   9863  0.331       1
#  5 134825 Female    89  2000 [60,90) (1985,2000]   9863  7.39        1
table(d$agedxG)#pretty even, consistent with median age at dx of ~59
d=d|>mutate(histo3=factor(histo3)) #for plots 
age.cuts <- c(20,60,90) 
(D=survSplit(Surv(agedx,agedx+surv,status)~.,data=d,cut=age.cuts,start="astart",end="astop")|>tibble())
 # A tibble: 59,544 × 9 #output here replaces surv and status above with tstart, tstop and a new status
 #       id sex     yrdx agedxG  yrdxG       histo3 astart astop status
 #    <int> <fct>  <dbl> <fct>   <fct>        <fct>  <dbl> <dbl>  <dbl>
 # 1   5253 Female  2002 [60,90) (2000,2010]   9863     71  85.0      1
 # 2   9874 Male    2002 [20,60) (2000,2010]   9875     33  54.7      0
 # 3  98271 Male    2005 [60,90) (2000,2010]   9863     87  87.3      1
 # 4 125911 Female  2009 [60,90) (2000,2010]   9863     87  87.3      1
 # 5 134825 Female  2000 [60,90) (1985,2000]   9863     89  90        0
 # 6 134825 Female  2000 [60,90) (1985,2000]   9863     90  96.4      1 #first case where survSplit() did some work
 # 7 158608 Female  2005 [60,90) (2000,2010]   9863     81  86.4      1
dups=D[duplicated(D$id),]
D[D$id%in%dups$id[1:4],] # show first 4 cases that were split
#       id sex     yrdx agedxG  yrdxG       histo3 astart astop status
#    <int> <fct>  <dbl> <fct>   <fct>        <fct>  <dbl> <dbl>  <dbl>
# 1 134825 Female  2000 [60,90) (1985,2000]   9863     89  90        0
# 2 134825 Female  2000 [60,90) (1985,2000]   9863     90  96.4      1
# 3 232693 Female  2003 [20,60) (2000,2010]   9875     49  60        0
# 4 232693 Female  2003 [20,60) (2000,2010]   9875     60  63.8      1
# 5 313894 Female  2004 [60,90) (2000,2010]   9863     77  90        0
# 6 313894 Female  2004 [60,90) (2000,2010]   9863     90  96.2      0
# 7 371865 Male    2017 [60,90) (2010,2023]   9875     89  90        0
# 8 371865 Male    2017 [60,90) (2010,2023]   9875     90  92.4      1
D=D|>mutate(PY = astop - astart, .after=status )
(D=D|>mutate(age = as.factor(astart),.after=status))
(rates=survRate(Surv(PY,status)~age, data=D)|>tibble()|>
  mutate(start=as.numeric(levels(age))[age],midA=start+0.5))
#  A tibble: 71 × 8
#    age   tstop event   rate  lower  upper start  midA
#    <fct> <dbl> <dbl>  <dbl>  <dbl>  <dbl> <dbl> <dbl>
#  1 20    1689.    30 0.0178 0.0120 0.0254    20  20.5   30/1689 is 0.01776199
#  2 21    1885.    44 0.0233 0.0170 0.0313    21  21.5
#  3 22    2298.    60 0.0261 0.0199 0.0336    22  22.5
#  4 23    2554.    69 0.0270 0.0210 0.0342    23  23.5
#  5 24    2816.    61 0.0217 0.0166 0.0278    24  24.5
#  6 25    2866.    83 0.0290 0.0231 0.0359    25  25.5
#  7 26    3204.    67 0.0209 0.0162 0.0266    26  26.5
#  8 27    3330.    81 0.0243 0.0193 0.0302    27  27.5
#  9 28    3287.    86 0.0262 0.0209 0.0323    28  28.5
# 10 29    3540.    91 0.0257 0.0207 0.0316    29  29.5
# ℹ 61 more rows
(rates=survRate(Surv(PY,status)~age+histo3, data=D)|>tibble()|>
  mutate(start=as.numeric(levels(age))[age],midA=start+0.5))
# # A tibble: 142 × 9
#    age   histo3 tstop event    rate   lower  upper start  midA
#    <fct>  <dbl> <dbl> <dbl>   <dbl>   <dbl>  <dbl> <dbl> <dbl>
#  1 20      9863 1053.    21 0.0199  0.0123  0.0305    20  20.5
#  2 20      9875  636.     9 0.0142  0.00647 0.0269    20  20.5
#  3 21      9863 1296.    31 0.0239  0.0163  0.0340    21  21.5
#  4 21      9875  589.    13 0.0221  0.0118  0.0378    21  21.5
(rates=survRate(Surv(PY,status)~age+histo3, data=D|>filter(yrdx>2014))|>tibble()|>
  mutate(start=as.numeric(levels(age))[age],midA=start+0.5))

rates|>ggplot(aes(x=midA,y=rate,group=histo3)) + geom_line(aes(col=histo3)) +
  geom_ribbon(aes(ymin=lower, ymax=upper, fill=histo3), alpha=0.3) +
  xlab('Years since Dx') + ylab('Hazard')  #9863 is worse at ages >55

##  To finish this script, we go back to time since Dx (as in survSplit.R, but now using formulas) 
time.cuts <- c(0, 5, 10, 20, 100)
(Dt=survSplit(Surv(surv,status)~.,data=d,cut=time.cuts,end="tstop",start="tstart",event="status")|>tibble())
# # A tibble: 91,737 × 10
#        id sex    agedx  yrdx agedxG  yrdxG       histo3 tstart  tstop status
#     <int> <fct>  <int> <dbl> <fct>   <fct>       <fct>   <dbl>  <dbl>  <dbl>
#  1   5253 Female    71  2002 [60,90) (2000,2010] 9863        0  5          0
#  2   5253 Female    71  2002 [60,90) (2000,2010] 9863        5 10          0
#  3   5253 Female    71  2002 [60,90) (2000,2010] 9863       10 14.0        1
#  4   9874 Male      33  2002 [20,60) (2000,2010] 9875        0  5          0
#  5   9874 Male      33  2002 [20,60) (2000,2010] 9875        5 10          0
#  6   9874 Male      33  2002 [20,60) (2000,2010] 9875       10 20          0
#  7   9874 Male      33  2002 [20,60) (2000,2010] 9875       20 21.7        0
Dt=Dt|>mutate(PY = tstop - tstart, .after=status )
(Dt=Dt|>mutate(fu = as.factor(tstart),.before=tstart))
levels(Dt$fu)
(ratesFU <- survRate(Surv(PY,status)~fu, data=Dt) |>
  mutate(start=as.numeric(levels(fu))[fu],
         mid=ifelse(start==0,2.5,ifelse(start==5,7.5,ifelse(start==10,15,30)))))
#       fu     tstop event       rate      lower      upper start  mid
# fu=0   0 174210.29 16102 0.09242852 0.09100634 0.09386736     0  2.5
# fu=5   5  92817.48  4579 0.04933338 0.04791470 0.05078339     5  7.5   #rates are high since using all years
# fu=10 10  70195.04  2678 0.03815085 0.03671946 0.03962374    10 15.0
# fu=20 20  44397.65   845 0.01903254 0.01777074 0.02036029    20 30.0

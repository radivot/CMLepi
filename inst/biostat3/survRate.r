# survRate.R    This biostat3 function is like SEERaBomb::fold() but for rates O/PY, instead of EAR=(O-E)/PY and RR=O/E
library(biostat3)         
library(tidyverse)   
load("~/data/CMLepi/cml20.RData") #made in mkSEER.R
(d=d20|>filter(histo3%in%c(9863,9875),agedx<90,surv<80,surv>0,yrdx>=2014)) # 21.3k CML cases
(d=d|>filter(agedx>=80)|>mutate(histo3=as_factor(histo3))|>select(-cancer)) #2.1k
(D <- survSplit(d, cut=0:9, end="surv", start="start", event="status")|>tibble())# 2.1 => 6.3
D=D|>mutate(PY = surv - start, .after=surv )
(D=D|>mutate(fu = as.factor(start),.before=start))
(rates=D|>group_by(fu)|>summarize(PY=sum(PY),O=sum(status),rate=O/PY)) #getting rate is easy enough, but
(rates|>mutate(LL=rate-1.96*sqrt(O)/PY,UL=rate+1.96*sqrt(O)/PY)) ##negative lower CI at bottom here is no good
# and doing it right is a bit uglier than calling survRate()
rates=rates|>mutate(test=map2(O,PY,function(O,PY) poisson.test(O,T=PY))) 
rates=rates|>mutate(LL=map_dbl(test,function(x) x$conf.int[[1]]))
rates=rates|>mutate(UL=map_dbl(test,function(x) x$conf.int[[2]]))
rates=rates|>select(-test)
rates
(ratesFU=survRate(Surv(PY,status)~fu, data=D)|>tibble()) # ... using survRate() will keep scripts cleaner  

#note that you can use survRate directly on the SEER data in d to sum up rates over all times and ages
load("~/data/CMLepi/cml.RData") #made in 1_mkSEERdata.R (unique cases of CML, with CMML cases removed)
d=d|>filter(agedx<90,agedx>=20,surv>0) # 50.46k
(d=d|>mutate(yrdxG=cut(yrdx,breaks=c(1975,1985,2000,2010,2023),include.lowest = T),.after=yrdx))
(d=d|>mutate(status=ifelse(status==1&surv<=13,1,0),surv=pmin(surv,13))) #truncate time to 13 years
survRate(Surv(agedx, agedx+surv, status) ~ sex+yrdxG, data=d) # i.e. to sum up progress over the decades
#                                  sex       yrdxG     tstop event       rate      lower      upper
# sex=Female, yrdxG=[1975,1985] Female [1975,1985]  3869.795  1059 0.27365793 0.25742210 0.29064943
# sex=Female, yrdxG=(1985,2000] Female (1985,2000] 13514.541  1938 0.14340109 0.13708701 0.14993098
# sex=Female, yrdxG=(2000,2010] Female (2000,2010] 55932.003  3269 0.05844597 0.05645943 0.06048456
# sex=Female, yrdxG=(2010,2023] Female (2010,2023] 58434.916  2845 0.04868664 0.04691388 0.05050925
# sex=Male, yrdxG=[1975,1985]     Male [1975,1985]  4521.066  1436 0.31762418 0.30140668 0.33448756
# sex=Male, yrdxG=(1985,2000]     Male (1985,2000] 18242.055  2656 0.14559763 0.14011255 0.15124240
# sex=Male, yrdxG=(2000,2010]     Male (2000,2010] 70273.264  4420 0.06289732 0.06105658 0.06477945
# sex=Male, yrdxG=(2010,2023]     Male (2010,2023] 75055.627  4364 0.05814354 0.05643112 0.05989472

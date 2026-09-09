# survSplit.R   Inspired by lab 7 in biostat3. 
graphics.off();rm(list=ls())#clear plots and environment 
library(biostat3) #loads survival, which contains survSplit()
library(tidyverse) # get this select()
library(broom)    # tidy()
library(bshazard)
load("~/data/CMLepi/cml20.RData") #made in mkSEER.R
(d=d20|>filter(histo3%in%c(9863,9875),agedx<90,surv<80,surv>0,yrdx>=2013)) # 23k CML cases
(d=d|>filter(agedx>=80)|>mutate(histo3=as_factor(histo3))|>select(-cancer)) #2.3k
(D <- survSplit(d, cut=0:9, end="surv", start="start", event="status")|>tibble())# 2.3k rows => 7.2k
# # A tibble: 7,191 × 12
#        id sex    agedx  yrdx histo3   COD COD2  COD7  CODS                           start   surv status
#     <int> <fct>  <int> <dbl> <fct>  <int> <chr> <chr> <fct>                          <dbl>  <dbl>  <dbl>
#  1 233367 Female    81  2021 9875       0 alive alive Alive                              0 1           0
#  2 233367 Female    81  2021 9875       0 alive alive Alive                              1 2           0
#  3 233367 Female    81  2021 9875       0 alive alive Alive                              2 2.09        0
#  4 320365 Male      87  2015 9875      78 LC    LC    Chronic Myeloid Leukemia           0 0.0986      1
#  5 371865 Male      89  2017 9875     154 OC    CV    Diseases of Heart                  0 1           0
#  6 371865 Male      89  2017 9875     154 OC    CV    Diseases of Heart                  1 2           0
#  7 371865 Male      89  2017 9875     154 OC    CV    Diseases of Heart                  2 3           0
#  8 371865 Male      89  2017 9875     154 OC    CV    Diseases of Heart                  3 3.38        1
#  9 482907 Male      88  2014 9875      86 OC    CA    Miscellaneous Malignant Cancer     0 1           0
# 10 482907 Male      88  2014 9875      86 OC    CA    Miscellaneous Malignant Cancer     1 2           0
#  Each patient is expanded into its contributions into each time-since-Dx bin 
#  Calculate person-years (PY) in each bin and recode start time as a factor
D=D|>mutate(PY = surv - start, .after=surv )
D=D|>mutate(fu = as.factor(start),.before=start)
## Calculate the incidence rate by observation year
(ratesFU <- survRate(Surv(PY,status)~fu, data=D) |>
  mutate(start=as.numeric(levels(fu))[fu],mid=start+0.5))

par(mfrow=c(1,1))
library(tinyplot) # lightweight base graphics extension
par(mfrow=c(1,1))
ratesFU|> # df with 10 rows, one for each fu interval
  with({
    plt(rate~mid, ymin=lower, ymax=upper,type="ribbon",
        ylab="Mortality rate per person-year",xlab="Years since Dx")
  })  #as in lifetab2.r, picture is consistent with background going from 0.1 to 0.2 between 85 and 95

summary(poi <- glm(status ~ fu + offset(log(PY)),family = poisson,data = D))
eform(poi) #from biostat3
#             exp(beta)     2.5 %    97.5 %
# (Intercept) 0.3444526 0.3183588 0.3726852  starts at 0.34 via bad cases
# fu1         0.5862868 0.5070283 0.6779348  
# fu2         0.5532571 0.4686466 0.6531435  drops to min of just over half that, i.e ~0.20
# fu3         0.6636325 0.5570133 0.7906599
# fu4         0.6838575 0.5574430 0.8389397
# fu5         0.7098441 0.5601483 0.8995451
# fu6         0.7602462 0.5756699 1.0040030
# fu7         0.7696617 0.5422833 1.0923793
# fu8         0.6487741 0.3887472 1.0827290
# fu9         1.0553044 0.6092994 1.8277835  and rises back up 1st year risks

summary(poi <- glm(status ~ fu  + histo3 + sex + offset(log(PY)),family = poisson,data = D))
eform(poi) 
#            exp(beta)     2.5 %    97.5 %
# (Intercept) 0.3477225 0.3131398 0.3861244
# fu1         0.5906960 0.5108329 0.6830449
# fu2         0.5629183 0.4768008 0.6645898
# fu3         0.6777298 0.5687928 0.8075307
# fu4         0.6970238 0.5681236 0.8551699
# fu5         0.7219793 0.5696940 0.9149721
# fu6         0.7689066 0.5821878 1.0155096
# fu7         0.7809994 0.5502045 1.1086061
# fu8         0.6619999 0.3965977 1.1050085
# fu9         1.0707455 0.6180697 1.8549623
# histo39875  0.7183570 0.6474249 0.7970605  #better to be 9875
# sexMale     1.2577225 1.1364114 1.3919834  #worse to be male, perhaps via background

D=D|>mutate(agedxG=as_factor(agedx),yrdxG=as_factor(yrdx))
summary(poi<-glm(status~fu+histo3+sex+agedxG+yrdxG+offset(log(PY)),family=poisson,data=D))
eform(poi) 
#             exp(beta)     2.5 %    97.5 %
# (Intercept) 0.1804174 0.1436996 0.2265171
# fu1         0.6225277 0.5371146 0.7215232
# fu2         0.5953831 0.5023866 0.7055941
# fu3         0.7156843 0.5977469 0.8568911
# fu4         0.7486725 0.6066881 0.9238858
# fu5         0.7886945 0.6182511 1.0061267
# fu6         0.8629990 0.6486163 1.1482401
# fu7         0.8933984 0.6246903 1.2776902
# fu8         0.7680309 0.4566865 1.2916333
# fu9         1.3647718 0.7775429 2.3954974
# histo39875  0.7218173 0.6494576 0.8022390  #holds
# sexMale     1.3122122 1.1847300 1.4534121  #holds
# agedxG81    1.2165330 0.9913838 1.4928148
# agedxG82    1.3569920 1.1016817 1.6714695
# agedxG83    1.4243932 1.1561220 1.7549151
# agedxG84    1.4247714 1.1461424 1.7711356
# agedxG85    1.7180642 1.3772095 2.1432793
# agedxG86    1.8255820 1.4786769 2.2538727
# agedxG87    1.6923997 1.3448901 2.1297032
# agedxG88    2.6318097 2.0938995 3.3079057
# agedxG89    2.3471808 1.8561645 2.9680870  #worse with age
# yrdxG2014   1.2340890 1.0074005 1.5117876
# yrdxG2015   1.3741773 1.1199282 1.6861468
# yrdxG2016   1.0993783 0.8951773 1.3501600
# yrdxG2017   1.4112368 1.1391366 1.7483322
# yrdxG2018   1.2715114 1.0188692 1.5867996
# yrdxG2019   1.2685227 1.0161087 1.5836396
# yrdxG2020   1.1709951 0.9127915 1.5022372
# yrdxG2021   0.9581247 0.7253343 1.2656275
# yrdxG2022   1.2906623 0.9777975 1.7036342
# yrdxG2023   1.8031372 1.2919186 2.5166477

## Test if the effect of yrdxG is significant using a likelihood ratio test
drop1(poi, ~yrdxG, test="Chisq")  # P = 0.002877 **  (all worse than 2013)
summary(poi<-glm(status~fu+histo3+sex+agedxG+yrdx+offset(log(PY)),family=poisson,data=D))
# Coefficients:
#               Estimate Std. Error z value Pr(>|z|)    
# (Intercept) -28.709411  19.752006  -1.453 0.146087    
# fu1          -0.499249   0.074335  -6.716 1.87e-11 ***
# fu2          -0.532263   0.085437  -6.230 4.67e-10 ***
# fu3          -0.330344   0.090873  -3.635 0.000278 ***
# fu4          -0.281836   0.106664  -2.642 0.008235 ** 
# fu5          -0.230692   0.123771  -1.864 0.062341 .  
# fu6          -0.155129   0.145357  -1.067 0.285870    
# fu7          -0.131376   0.182249  -0.721 0.470995    
# fu8          -0.285215   0.264692  -1.078 0.281240    
# fu9           0.207161   0.284615   0.728 0.466697    
# histo39875   -0.312227   0.053490  -5.837 5.31e-09 ***
# sexMale       0.262024   0.051997   5.039 4.68e-07 ***
# agedxG81      0.179002   0.104042   1.720 0.085346 .  
# agedxG82      0.301244   0.106130   2.838 0.004533 ** 
# agedxG83      0.346534   0.106089   3.266 0.001089 ** 
# agedxG84      0.339422   0.110650   3.068 0.002158 ** 
# agedxG85      0.524273   0.112508   4.660 3.16e-06 ***
# agedxG86      0.592748   0.107167   5.531 3.18e-08 ***
# agedxG87      0.507540   0.116531   4.355 1.33e-05 ***
# agedxG88      0.934279   0.116326   8.032 9.62e-16 ***
# agedxG89      0.835741   0.119526   6.992 2.71e-12 ***
# yrdx          0.013491   0.009788   1.378 0.168132      slope is not sigificant

summary(poi<-glm(status~histo3+sex+fu*agedxG+offset(log(PY)),family=poisson,data=D))
# Coefficients:
#                Estimate Std. Error z value Pr(>|z|)    
# (Intercept)    -1.59839    0.13237 -12.076  < 2e-16 ***
# histo39875     -0.31438    0.05352  -5.874 4.24e-09 ***
# sexMale         0.25566    0.05229   4.889 1.01e-06 ***
# ...  with tons of interaction nuissance params

library(splines)
D=D|>mutate(mid=start+0.5)
summary(poi<-glm(status~histo3+sex+ns(mid,df=3)+agedxG+yrdxG+offset(log(PY)),family=poisson,data=D))
(df0=tibble(agedxG="85", yrdxG="2019",mid=1:10, PY=1,histo3="9875",sex="Male")) #baseline hazard predictions
pred0=predict(poi,newdata=df0, se.fit=TRUE)
(df0=with(pred0, bind_cols(df0,fit=exp(fit),conf.low=exp(fit-1.96*se.fit),conf.high=exp(fit+1.96*se.fit))))
# # A tibble: 10 × 9
#    agedxG yrdxG   mid    PY histo3 sex     fit conf.low conf.high
#    <chr>  <chr> <int> <dbl> <chr>  <chr> <dbl>    <dbl>     <dbl>
#  1 85     2019      1     1 9875   Male  0.281    0.219     0.360
#  2 85     2019      2     1 9875   Male  0.221    0.171     0.285
#  3 85     2019      3     1 9875   Male  0.243    0.189     0.313
#  4 85     2019      4     1 9875   Male  0.268    0.206     0.348
#  5 85     2019      5     1 9875   Male  0.290    0.222     0.379
#  6 85     2019      6     1 9875   Male  0.311    0.237     0.409
#  7 85     2019      7     1 9875   Male  0.331    0.248     0.442
#  8 85     2019      8     1 9875   Male  0.349    0.251     0.486
#  9 85     2019      9     1 9875   Male  0.368    0.248     0.546
# 10 85     2019     10     1 9875   Male  0.386    0.240     0.622
library(tinyplot)
with(df0,plt(fit~mid,ymin=conf.low,ymax=conf.high,type="ribbon",ylab="Rate",xlab="Years since Dx"))


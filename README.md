
# CMLepi

<!-- badges: start -->
[![R-CMD-check](https://github.com/radivot/CMLepi/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/radivot/CMLepi/actions/workflows/R-CMD-check.yaml)
<!-- badges: end -->

The immediate goal of CMLepi is to help people estimate chronic myeloid leukemia (CML) patient life expectancies (LEs) from 
Surveillance Epidemiology and End Results (SEER) data. 

## Installation

You can install the development version of CMLepi like so:

``` r
remotes::install_github("radivot/CMLepi")
```

To use it you must first gain access to the SEER data via a Windows program called SEER\*stat. 
This involves requesting access (no need for the Plus version) and waiting a day to get it. 
Once you have it,  work through the SEER\*Stat Case Listing tutorial and create one using 
Incidence - SEER Research Data, 8 Registries, Nov 2025 Sub (1975-2023), i.e. SEER8. Select 
cases by Site and Morphology.Site recode ICD-O-3/WHO 2008 = Chronic Myeloid Leukemia and choose as 
column variables: Patient ID, Sex, Age recode with single ages and 90+, 
Year of diagnosis, ICD-O-3 Hist/behav, Survival Days, and COD to site recode. Execute the session. When it
completes, select all, right-click on the header, pick display unformatted, and
export to the files `cml8.txt` and `cml8.dic`. Repeat this changing the database to SEER12 and 
SEER20 (i.e. SEER21 excluding IL) and save the results to `cml12.txt/cml12.dic` and `cml20.txt/cml20.dic`. Put these six
files in the folder `~/data/CMLepi`. 



## Introduction

Bringing dic/txt file pairs into single R tibbles is done using `seer2r()`. This function is a wrapper around 
R package **SEER2R**'s function `read.SeerStat()`. The following work is done by `seer2r()`. 
``` r
# pak::pak("cran/SEER2R") #installs fine from  GitHub even though no longer on CRAN
library(SEER2R)
n8 = read.SeerStat("~/data/CMLepi/cml8.dic",UseVarLabelsInData=FALSE) #get numbers(n)
head(n8<-attr(n8,"assignColNames")(n8,c("id","sex","Age","Year","ICDO3","surv","COD")))
library(tidyverse)
(n8=n8|>rename(yrdx=Year,agedx=Age)|>as_tibble())
n8=n8|>mutate(yrdx=yrdx+1800,surv=surv/365.25)
n8=n8|>mutate(status=as.numeric(COD>0),.after=surv)
n8=n8|>mutate(sex=ifelse(sex==1,"Male","Female"))
(n8=n8|>mutate(sex=as_factor(sex)))

c8 = read.SeerStat("~/data/CMLepi/cml8.dic",UseVarLabelsInData=TRUE) 
head(c8<-attr(c8,"getSubDataByVarName")(c8,c("ICDO3","site")))
(c8=c8|>rename(histS=ICDO3,CODS=site)|>as_tibble())
(d8=bind_cols(n8,c8))
(d8=d8|>mutate(histo3=as.numeric(str_sub(histS,end=4)),.after=ICDO3))
table(d8$histo3)
#  9863  9875  9876  9945  9946 
# 10856  4065   116  4331    63
d8=d8|>filter(histo3%in%c(9863,9875,9945)) #leave out rare jCMML=9946 and atypical CML=9876
table(d8$histo3,d8$yrdx)
#      1975 1976 1977 1978 1979 1980 1981 1982 1983 1984 1985 1986 1987 1988 1989 1990 1991
# 9863  244  264  222  237  270  250  242  239  254  251  304  241  263  279  230  246  263
# 9875    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0
# 9945    0    0    0    0    0    0    0    0    0    1    0   34   34   33   42   45   63
#      1992 1993 1994 1995 1996 1997 1998 1999 2000 2001 2002 2003 2004 2005 2006 2007 2008
# 9863  273  273  291  296  261  275  266  275  275  254  212  220  210  222  239  241  240
# 9875    0    0    1    2    1    1    2    4    2   38   35   39   58   48   68   74  102
# 9945   50   70   72   76   90   95   72   98   91   89  111  109  121   95  108  111  103
#      2009 2010 2011 2012 2013 2014 2015 2016 2017 2018 2019 2020 2021 2022 2023
# 9863  209  190  165  142  154  127  137  157  134  138  170  107  137  131  136
# 9875   82  164  211  236  223  258  241  291  276  230  284  274  307  259  254
# 9945  143  121  132  122  144  150  145  163  182  181  184  185  221  198  247
(d8=d8|>mutate(cancer=ifelse(histo3%in%c(9863,9875),"CML","CMML"),.after=histo3))
(d8=d8|>select(-histS))
(d8=d8|>mutate(COD2=ifelse(COD==0,"alive",ifelse((COD>=74)&(COD<=85)|(COD==89),"LC","OC")),.after=COD))
(d8=d8|>mutate(CODS=as_factor(CODS))) #to save a little memory

mapCOD7=function(D){
  COD=D$COD #start with vec of integers. Map to a vec of Strings
  CODS=rep("UNK",dim(D)[1]) #set default to "unknown" type of death
  CODS[COD==0]="alive"
  CODS[(COD>=1)&(COD<=73)|(COD==86)|(COD==90)]="CA"
  CODS[(COD>=74)&(COD<=85)|(COD==89)]="LC"
  CODS[COD==130]="CA" #in situ (benign)
  CODS[(COD>=133)&(COD<=145)]="IN" # infection
  CODS[COD==148]="DK"  # diabetes
  CODS[COD==151]="YOC"  # alzheimers -> yet other causes
  CODS[COD==154]="CV"  # heart disease
  CODS[COD==157]="CV" # hypertension without HD
  CODS[COD==160]="CV"  #cerebroVasc"
  CODS[COD==163]= "CV" #"athero"
  CODS[COD==166]= "CV"     #"aoritic aneurysm"
  CODS[COD==169]= "CV"  #"other disease of Vasc"
  CODS[COD==172]= "IN" #"pneumonia"
  CODS[COD==175]="YOC" # COPD, no signal so smoking makes both. chronic obstructive pulminary disease
  CODS[COD==178]="YOC" # ulcer
  CODS[COD==181]="YOC" # liver disease
  CODS[COD==184]="DK" # kidney disease
  CODS[COD==199]="ASH" #"accidents"
  CODS[COD==202]="ASH"  #"suicide"
  CODS[COD==205]="ASH" # homocide"
  CODS[COD%in%c(187,190,193)]="YOC" # perinatal conditions
  CODS[COD%in%c(196,208,252)]="YOC" # other causes, including ill-defined and unknown
  # CODS[COD==252]="UNK" # same if no comment => all accounted for
  D$COD7=as.factor(CODS)
  D|>relocate(COD7, .after = COD2)
}
(d8=mapCOD7(d8))
# # A tibble: 19,252 × 13
#      id sex    agedx  yrdx ICDO3 histo3 cancer  surv status   COD COD2  COD7  CODS                          
#   <int> <fct>  <int> <dbl> <int>  <dbl> <chr>  <dbl>  <dbl> <int> <chr> <fct> <fct>                         
# 1  2075 Female    80  1990  7783   9945 CMML    1.79      1    78 LC    LC    Chronic Myeloid Leukemia      
# 2  2226 Male      86  1988  7455   9863 CML     1.11      1   154 OC    CV    Diseases of Heart             
# 3  4093 Female    50  1989  7455   9863 CML     8.10      1   154 OC    CV    Diseases of Heart             
# 4  5253 Female    71  2002  7455   9863 CML    14.0       1   208 OC    YOC   Other Cause of Death          
# 5  6614 Female    81  2014  7783   9945 CMML    1.68      1    85 LC    LC    Aleukemic, Subleukemic and NOS
# 6  8674 Male      63  1998  7503   9875 CML     3.72      1    78 LC    LC    Chronic Myeloid Leukemia      
# 7  8686 Female    77  1998  7455   9863 CML     4.39      1    78 LC    LC    Chronic Myeloid Leukemia      
# 8  8767 Male      51  1995  7455   9863 CML     3.67      1    78 LC    LC    Chronic Myeloid Leukemia      
# 9  8938 Male      58  1997  7455   9863 CML     6.16      1    78 LC    LC    Chronic Myeloid Leukemia      
#10  8958 Male      52  1997  7455   9863 CML     4.36      1    78 LC    LC    Chronic Myeloid Leukemia      
```

We use `seer2r()` to make SEER CML incidence binary files that load much faster 
downstream as follows. 


``` r
# mkSEERincid.R  (name of this R script)
library(tidyverse)
library(CMLepi)
d8=seer2r("cml8")
system.time(save(d8,file="~/data/CMLepi/cml8.RData")) 
load("~/data/CMLepi/cml8.RData")
d12=seer2r("cml12") # few secs
system.time(save(d12,file="~/data/CMLepi/cml12.RData")) 
load("~/data/CMLepi/cml12.RData") 
d20=seer2r("cml20") # 5 secs
system.time(save(d20,file="~/data/CMLepi/cml20.RData")) 
load("~/data/CMLepi/cml20.RData") 

######## using prevalence sessions in SEER*stat I get on jan1 2023, pop average of 2022 and 2023 is 140,038,579.0
# total is (336.75+334)/2 = 335   based on https://www.multpl.com/united-states-population/table/by-year
335/140.04 #2.392 = multiplier for SEER20 pop of 140,038,579.0
2.392*25451 #61k (histo3 sum)   23-year limited duration prevalence
335/40.355 #8.3 is multiplier  SEER12 pop of 40,355,231.5
8.3*7419.4 #61.5 k     31-year limited duration prevalence
335/27.69 #12.1 is multiplier for SEER8 pop of 27,692,133.5
12.1*5205.2 #63k  48-year limited duration prevalence. Why not 66k? why not use 5500 alive at end of 2023?
# maybe they think loss of follow up doesn't mean confirmed alive, just not confirmed dead yet ... tut 1 supports this  

dc=d8|>filter(yrdx<2000,cancer=="CML")
(dc=dc|>mutate(year=yrdx+surv,status=ifelse(year>2000,0,1)))
table(dc$status) # 1206 alive; 12*1206 = 14472 is prevalence of cases alive entering 2000 
dc=d8|>filter(cancer=="CML")
table(dc$status) # 5500 alive, 9419 dead;12*5500 = 66000 is prevalence on dec31 2023

# now getnumbers of new cases each year in each database
round(table(d8$cancer,d$yrdx)*12.1)
#      1975 1976 1977 1978 1979 1980 1981 1982 1983 1984 1985 1986 1987 1988 1989 1990 1991 1992 1993 1994 1995 1996 1997 1998 1999 2000
# CML  2952 3194 2686 2868 3267 3025 2928 2892 3073 3037 3678 2916 3182 3376 2783 2977 3182 3303 3303 3533 3606 3170 3340 3243 3376 3352
# CMML    0    0    0    0    0    0    0    0    0   12    0  411  411  399  508  544  762  605  847  871  920 1089 1150  871 1186 1101
#      2001 2002 2003 2004 2005 2006 2007 2008 2009 2010 2011 2012 2013 2014 2015 2016 2017 2018 2019 2020 2021 2022 2023
# CML  3533 2989 3134 3243 3267 3715 3812 4138 3521 4283 4550 4574 4562 4658 4574 5421 4961 4453 5493 4610 5372 4719 4719
# CMML 1077 1343 1319 1464 1150 1307 1343 1246 1730 1464 1597 1476 1742 1815 1754 1972 2202 2190 2226 2238 2674 2396 2989
round(table(d12$cancer,d12$yrdx)*8.3)
#      1992 1993 1994 1995 1996 1997 1998 1999 2000 2001 2002 2003 2004 2005 2006 2007 2008 2009 2010 2011 2012 2013 2014 2015 2016 2017
# CML  3420 3312 3494 3777 3279 3420 3378 3486 3378 3362 2980 3196 3204 3486 3619 3760 4075 3544 4274 4241 4457 4532 4532 4308 5063 5005
# CMML  622  805  797  930 1021 1087  872 1112 1029 1071 1170 1237 1270 1154 1320 1212 1237 1519 1386 1411 1328 1519 1685 1594 1751 2025
#      2018 2019 2020 2021 2022 2023
# CML  4424 5013 4698 5171 4606 4714
# CMML 2108 2224 2050 2523 2258 2872
round(table(d20$cancer,d20$yrdx)*2.4) 
#      2000 2001 2002 2003 2004 2005 2006 2007 2008 2009 2010 2011 2012 2013 2014 2015 2016 2017 2018 2019 2020 2021 2022 2023
# CML  3698 3650 3206 3605 3715 3746 3739 3888 4063 4320 4438 4754 4687 4987 5210 5218 5297 5167 5218 5398 5093 5623 5261 5549
# CMML  883  926  934 1068 1044 1013 1099 1116 1154 1363 1354 1433 1313 1606 1651 1606 1819 1819 1939 2134 2076 2318 2278 2645
```

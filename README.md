
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

To use it you must first gain access to the SEER data via a Windows program called SEER*stat. 
This involves requesting access (no need for the Plus version) and waiting a day to get it. 



## Introduction

There is overhead. First work through SEER*Stat's Case Listing tutorial. 
Next, create a Case Listing session using 
Incidence - SEER Research Data, 8 Registries, Nov 2025 Sub (1975-2023) as the database (i.e. SEER8), selecting 
cases with Site and Morphology.Site recode ICD-O-3/WHO 2008} = '      Chronic Myeloid Leukemia'.
Next, choose the following variables as columns: Patient ID, Sex, Age recode with single ages and 90+, 
Year of diagnosis, ICD-O-3 Hist/behav, Survival Days, and COD to site recode. Then execute (under actions)
to create the listing. Select all and right-clicking on the header, display as unformatted raw numbers. Finally,
export, changing the file names from export.txt and export.dic to cml8.txt and cml8.dic. Repeat for SEER12 
and SEER20 (i.e. SEER21 excluding IL), calling those files cml12.txt and cml12.dic, and cml20.txt and cml20.dic.

To bring these files into R use the R package SEER2R. This package is no longer on CRAN
but it is still useful and can be installed from source via 15-year old CRAN read-only files on GitHub.
``` r
# pak::pak("cran/SEER2R") installs fine even though not on CRAN anymore due to 6 help page notes on a check
library(SEER2R)
n8 = read.SeerStat("Rpacks/SEER2R/cml8.dic",UseVarLabelsInData=FALSE) #get numbers(n)
head(n8<-attr(n8,"assignColNames")(n8,c("id","sex","Age","Year","ICDO3","surv","COD")))
library(tidyverse)
(n8=n8|>rename(yrdx=Year,agedx=Age)|>as_tibble())
n8=n8|>mutate(yrdx=yrdx+1800,status=as.numeric(COD>0),surv=surv/365.25)
n8=n8|>mutate(sex=ifelse(sex==1,"Male","Female"))
(n8=n8|>mutate(sex=as_factor(sex)))

c8 = read.SeerStat("Rpacks/SEER2R/cml8.dic",UseVarLabelsInData=TRUE) 
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

mapCOD6=function(D){
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
  D$COD6=as.factor(CODS)
  D
}
(d8=mapCOD6(d8))
# # A tibble: 19,252 × 13
#      id sex    agedx  yrdx ICDO3 histo3 cancer  surv   COD COD2  status CODS                           COD6 
#   <int> <fct>  <int> <dbl> <int>  <dbl> <chr>  <dbl> <int> <chr>  <dbl> <chr>                          <fct>
# 1  2075 Female    80  1990  7783   9945 CMML    1.79    78 LC         1 Chronic Myeloid Leukemia       LC   
# 2  2226 Male      86  1988  7455   9863 CML     1.11   154 OC         1 Diseases of Heart              CV   
# 3  4093 Female    50  1989  7455   9863 CML     8.10   154 OC         1 Diseases of Heart              CV   
# 4  5253 Female    71  2002  7455   9863 CML    14.0    208 OC         1 Other Cause of Death           YOC  
# 5  6614 Female    81  2014  7783   9945 CMML    1.68    85 LC         1 Aleukemic, Subleukemic and NOS LC   
# 6  8674 Male      63  1998  7503   9875 CML     3.72    78 LC         1 Chronic Myeloid Leukemia       LC   
# 7  8686 Female    77  1998  7455   9863 CML     4.39    78 LC         1 Chronic Myeloid Leukemia       LC   
# 8  8767 Male      51  1995  7455   9863 CML     3.67    78 LC         1 Chronic Myeloid Leukemia       LC   
# 9  8938 Male      58  1997  7455   9863 CML     6.16    78 LC         1 Chronic Myeloid Leukemia       LC   
#10  8958 Male      52  1997  7455   9863 CML     4.36    78 LC         1 Chronic Myeloid Leukemia       LC 

```


# mkSEERincid.R  (name of this R script)
library(tidyverse)
library(CMLepi)
# system.time({
#   d8=seer2r("cml8")
#   save(d8,file="~/data/CMLepi/cml8.RData")
#   d12=seer2r("cml12")
#   save(d12,file="~/data/CMLepi/cml12.RData")
#   d20=seer2r("cml20") 
#   save(d20,file="~/data/CMLepi/cml20.RData") 
# })  # 10 secs, so comment out
load("~/data/CMLepi/cml8.RData")
load("~/data/CMLepi/cml12.RData") 
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
round(table(d8$cancer,d8$yrdx)*12.1)
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

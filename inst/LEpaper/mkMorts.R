# mkMorts.R  
# Go into SEER*stat freq session with crude rate on mort data and let cod be pages, sex and age be rows, and year be a column.
# Rate, count and pop will then be a triplet within each matrix data cell. Export to ratesSex.txt and ratesSex.dic
# and put them in ~/data/CMLepi
search()
library(tidyverse)
search()
# vars=c("COD","sex","a20","year","rate","num","denom")
# d=readr::read_csv("~/data/CMLepi/ratesSex.txt", col_names=vars,skip=0) # less robust => switch to SEER2R
head(n<-SEER2R::read.SeerStat("~/data/CMLepi/ratesSex.dic",UseVarLabelsInData=FALSE),3) #get numbers(n)
head(n<-attr(n,"assignColNames")(n,c("cause","sex","age","rate")),3)
head(n<-n|>rename(COD=cause,year=Year_of_death,a20=age,num=Count,denom=Population)|>as_tibble(),3)
#     COD   sex   a20  year  rate     num     denom
#   <int> <int> <int> <int> <dbl>   <int>     <int>
# 1     0     0     0     0  937. 1959092 209143907
# 2     0     0     0     1 2148.   75073   3494217
# 3     0     0     0     2 2128.   74667   3508906
# head(d,3) #same answer. Using dic is more robust to how export was done in SEER*stat (comma vs tab delim)
d=n|>mutate(year=year+1968)
d=d|>filter(year!=1968)  #1968 is the sum over all years 
d=d|>filter(year>=1975)  
d=d|>filter(a20!=20)  #only known ages 
d=d|>mutate(ageMid=ifelse(a20==19,92.3,ifelse(a20==0,0.5,ifelse(a20==1,3,(a20-1)*5+2.5))))
(d=d|>mutate(rate2=num/denom,sex=ifelse(sex==0,"Both",ifelse(sex==1,"Male","Female"))))
tail(d)
SEERaBomb::mkMrtLocal()  #makes ~/data/mrt/mrtUSA.RData from Mx_1x1.txt in ~/data/hmd_countries/USA" 4/6/26 => thru 2024
load("~/data/mrt/mrtUSA.RData")# mrt is a list of 3 matrices. Used by SEERaBomb::msd() in Figs 5, 6, and 7
mrt[["Both"]][1:5,1:5]
SEERaBomb::mkAges() #makes ~/data/mrt/Ages.RData from ~/data/mrt/mrtUSA.RData 
load("~/data/mrt/Ages.RData") # Ages = list of PY-weighted 5-year bin age-group midpoints 
getAge=function(a20,year,sex) Ages[[sex]][[as.character(a20),as.character(year)]]
(d=d|>mutate(age=pmap_dbl(list(a20,year,sex),getAge))) #update ages 
save(d,file="~/data/CMLepi/seerMrt.RData") # 336,000 × 10

#the following lines make Generalized Additive Models (GAM): G2 for LC, OC and AC, G6 for CV, IN, CA, DK, ASH and YOC instead of OC
if(0) { # switch to 1 if you need to run these to make G*.RData files ... they take 80 secs
  system.time(SEERaBomb::mkG2(seerHome="~/data/CMLepi"))#30s so comment after run. Puts .RData files Glc, Goc and Gac in ~/data/CMLepi
  system.time(SEERaBomb::mkG6(seerHome="~/data/CMLepi"))#50s. Puts G6.RData in ~/data/CMLepi 
}
# G6.RData is used by Figure 7CD. Gac.RData is used by Figures 3. 

# also do this block once a year for ageDx specific LT estimates of LEs in 
files=c("Deaths_1x1.txt", "Exposures_1x1.txt", "Population.txt", "Mx_1x1.txt")
files=paste0("~/data/hmd_countries/USA/stats/",files)
us_mort=vital::read_hmd_files(files)
save(us_mort,file="~/data/mrt/us_mort.RData") #us_mort is one big tibble used by life table methods in R package vital



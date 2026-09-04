# 2_MortalityDataAB.R
library(tidyverse)
# In SEER*stat, use a frequency session to create the mortality table.
# Note that case listing makes no sense here since there is no individual level mortality data 
nms=c("a20","year","count")
system.time(d <-read_tsv("~/data/CMLepi/CMLdeaths.txt", col_names=nms)) #
d=d|>mutate(age=ifelse(a20==19,92.5,ifelse(a20==0,0.5,ifelse(a20==1,3,(a20-1)*5+2.5))))
d=d|>mutate(year=year+1968)
d=d|>filter(year!=1968,a20!=20)  #1968 is the sum over all years 
head(d)
#     a20  year count   age
#   <dbl> <dbl> <dbl> <dbl>
# 1     0  1969     1   0.5
# 2     0  1970     4   0.5
sum(d$count) # 92610
d|>group_by(year)|>summarize(n=sum(count))|>
  ggplot(aes(x=year,y=n))+geom_point()+labs(y="Number of Deaths by CML",x="Year of Death")+
  ylim(c(0,NA))+theme_classic(base_size=14)+
  geom_hline(yintercept=c(1200),col="gray")+
  annotate('text',x=1980,y=1300,label='1200 deaths/year',col="gray50",size=4.5)+
  annotate('text',x=2004.7,y=20,label='2007',col="gray50",size=3)+
  annotate('text',x=2014.7,y=20,label='2017',col="gray50",size=3)+
  annotate('text',x=2026.3,y=20,label='2024',col="gray50",size=3)+
  annotate('text',x=1971.3,y=20,label='1969',col="gray50",size=3)+
  geom_vline(xintercept=c(1969,2000,2007,2017,2020,2024),col="gray")
ggsave("LE/outs/2A_deathCounts.pdf",width=5,height=3)
ggsave("LE/outs/2A_deathCounts.png",width=5,height=3)
d|>group_by(year)|>summarize(mage=weighted.mean(age,count))|> filter(year>1998)|>
  ggplot(aes(x=year,y=mage))+geom_point()+labs(y="Mean age at death",x="Year of Deaths")+
  scale_y_continuous(breaks=c(65,70,75))+
  scale_x_continuous(breaks=c(1999,2010,2024))+
  theme_classic(base_size=14)
ggsave("LE/outs/2B_deathAges.pdf",width=3,height=3)  
ggsave("LE/outs/2B_deathAges.png",width=3,height=3)  
## Compare death counts in 2000-2023 nationwide to those estimated via SEER incidence COD info.
d|>filter(year%in%c(2000:2023))|>summarize(n=sum(count)) #total of 28085 deaths by CML in 2000-2023
5210*2.4# 12.5k by CML with Dx in 2000-2023 + 14.5k alive in 2000 implies an upper limit of 27k dead by CML. 
# 19.5k + 14.5k = 34k is at least greater than 28k, so Mortality data likely fixed the problem of missclassifications  
# of deaths as by other leukemias in the Incidence Data. It likely also does a better job of picking up
# deaths by intense therapy of intense disease, as 8.5k deaths by CML out of 14.5k alive in 2000 still seems high.  
# Unclear is if the 28k deaths include any via higher rates of CVD deaths caused by chronic use of tyrosine kinase inhibitors. 


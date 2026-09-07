# 63v75.R diffs in survival via diffs in care quality??
library(tidyverse)  
library(survival)  
load("~/data/CMLepi/cml20.RData") #made in mkSEER.R  60,882 cases including CMML (9845)
table(d20$histo3)
#  9863  9875  9945 
# 27632 18004 15246 
(d20=d20|>filter(histo3%in%c(9863,9875))|>mutate(histo3=as_factor(histo3))) # 45,636 CML cases
table(d20$histo3,d20$yrdx)
#      2000 2001 2002 2003 2004 2005 2006 2007 2008 2009 2010 2011 2012 2013 2014 2015 2016 2017 2018 2019 2020 2021 2022 2023
# 9863 1505 1262 1112 1214 1249 1251 1243 1282 1298 1319 1190 1128 1099 1142 1126 1145 1065 1040 1040 1067  940 1007  953  955
# 9875   36  259  224  287  299  310  315  338  395  481  659  853  853  936 1045 1029 1142 1113 1134 1182 1182 1336 1239 1357
#### shows SEER20 data in Figure 1B, i.e. rising use of 9875 and falling use of 9863, with changes now being very small 

dNA=d20|>filter(surv>80) # 89.7 => surv not available (NA). This mostly amounts to Dx via a death certificate
table(dNA$histo3,dNA$yrdx)
#      2000 2001 2002 2003 2004 2005 2006 2007 2008 2009 2010 2011 2012 2013 2014 2015 2016 2017 2018 2019 2020 2021 2022 2023
# 9863   25   20   32   18   22   18   27   14   19   16   19   21   17   14   34   20   26   27   25   23   23   23   27   15
# 9875    0    1    0    0    0    0    1    0    0    0    1    2    1    0    0    0    0    0    0    1    0    0    1    0
# So 9875 correlates with getting the Dx before death  
d0=d20|>filter(surv==0) # 0 => Dx while alive (i.e. via heath care reporting), but no second visit 
table(d0$histo3,d0$yrdx) # here it seems more random, as 9875 numbers with S=0 increase with their totals 
# 2000 2001 2002 2003 2004 2005 2006 2007 2008 2009 2010 2011 2012 2013 2014 2015 2016 2017 2018 2019 2020 2021 2022 2023
# 9863   12   24   10   11    8    5    4    8   10    5    7    8    9    8    5    5    2    7    6    5    0    9    9   15
# 9875    0    1    0    0    0    0    1    0    0    0    1    2    1    0    0    0    1    1    5    0    4    2    2   11
d90=d20|>filter(agedx==90) #1085 only known to be over 90
table(d90$histo3,d90$yrdx) # from 2016 on, 2X more 9863 in years when 9875 alread surpassed 9863 => significant bias to use 9863 in elderly 
#      2000 2001 2002 2003 2004 2005 2006 2007 2008 2009 2010 2011 2012 2013 2014 2015 2016 2017 2018 2019 2020 2021 2022 2023
# 9863   41   40   46   41   31   40   38   33   41   31   43   39   32   33   38   41   38   39   32   32   30   44   31   30
# 9875    0    3    1    0    0    0    2    6    1    2    6   12   13    9    6   11   20    9   16   17   13   16   22   16
(d=d20|>filter(agedx<90,surv<80,surv>0)) # 43,932 CML cases
(d=d|>filter(agedx>=80,yrdx>2015)) #1.6k
#        id sex    agedx  yrdx histo3 cancer    surv status   COD COD2  COD7  CODS                         
#     <int> <fct>  <int> <dbl> <fct>  <chr>    <dbl>  <dbl> <int> <chr> <fct> <fct>                        
# 1  233367 Female    81  2021 9875   CML    2.09         0     0 alive alive Alive                        
# 2  371865 Male      89  2017 9875   CML    3.38         1   154 OC    CV    Diseases of Heart            
# 3  733354 Female    82  2021 9863   CML    0.575        1   199 OC    ASH   Accidents and Adverse Effects
# 4  871256 Female    83  2016 9875   CML    7.15         0     0 alive alive Alive                        
# 5  978678 Male      86  2017 9875   CML    0.192        1   154 OC    CV    Diseases of Heart            
# 6 1013229 Male      80  2019 9863   CML    3.32         1   208 OC    YOC   Other Cause of Death         
# 7 1371723 Male      81  2018 9863   CML    5.05         0     0 alive alive Alive                        
# 8 1420914 Female    84  2017 9875   CML    1.02         1   154 OC    CV    Diseases of Heart            
# 9 1427312 Male      86  2016 9863   CML    6.46         1    78 LC    LC    Chronic Myeloid Leukemia     
#10 1431561 Male      81  2017 9863   CML    0.00274      1   172 OC    IN    Pneumonia and Influenza      
table(d$histo3)
# 9863 9875 
#  914  736 
table(d$histo3)/sum(table(d$histo3))
#      9863      9875 
# 0.5539394 0.4460606 
(sf=survfit(Surv(surv, status) ~ histo3, data = d)) 
par(mfrow=c(1, 1))
sf|>plot(col=1:2,xlab = "Years Since Dx", ylab = "Survival")
legend("topright", levels(d$histo3), col=1:2, lty = 1) # 9875 has better survival
# also for specific ages
(sf=survfit(Surv(surv, status) ~ histo3, data = d|>filter(agedx==80))) 
sf|>plot(col=1:2,xlab = "Years Since Dx", ylab = "Survival")
legend("topright", levels(d$histo3), col=1:2, lty = 1) 
(sf=survfit(Surv(surv, status) ~ histo3, data = d|>filter(agedx==85))) 
sf|>plot(col=1:2,xlab = "Years Since Dx", ylab = "Survival")
legend("topright", levels(d$histo3), col=1:2, lty = 1) 
# but by 89 at Dx, background dominates, so gap tightens
(sf=survfit(Surv(surv, status) ~ histo3, data = d|>filter(agedx==89))) 
sf|>plot(col=1:2,xlab = "Years Since Dx", ylab = "Survival")
legend("topright", levels(d$histo3), col=1:2, lty = 1) 
#### but hazard is still seen to be bigger
library(bshazard) # smooth hazards
as.data.frame.bshazard <- function(x, ...)  #code from labs in R package biostat3
  with(x, data.frame(Time=time,Hazard=hazard,lower.ci,upper.ci)) 
(bs=bshazard(Surv(surv, status)~histo3, d|>filter(agedx==89))) 
bs|>plot(overall=FALSE, col=1:2, lty=1, ylim=c(0,1))
legend("topright", levels(d$histo3), col=1:2, lty = 1) 
#plot above assumed proportional hazards between 63 and 75, next plot does not
library(tinyplot)
(H=d|>group_by(histo3)|>do(as.data.frame(bshazard(Surv(surv,status)~1,data=.,verbose=FALSE)))|>ungroup())
with(H, plt(Hazard~Time|histo3, ymin=lower.ci, ymax=upper.ci,type="ribbon",lty=1,  ylim=0:1))  #so gap is 
ggplot(H,aes(x=Time,y=Hazard,group=histo3)) + geom_line(aes(col=histo3)) +
  geom_ribbon(aes(ymin=lower.ci, ymax=upper.ci, fill=histo3), alpha=0.3) +ylim(0,1) 

biostat3::survRate(Surv(surv, status) ~ histo3, data=d) 
#             histo3    tstop event      rate     lower     upper
# histo3=9863   9863 1894.439   562 0.2966577 0.2726349 0.3222299
# histo3=9875   9875 1739.132   376 0.2161998 0.1948965 0.2391962
d|>group_by(histo3) |>
  summarise(D = sum(status), M = sum(surv), Rate = D/M,
            CI_low = stats::poisson.test(D,M)$conf.int[1],
            CI_high = stats::poisson.test(D,M)$conf.int[2]) 
#   histo3     D     M  Rate CI_low CI_high
# 1 9863     562 1894. 0.297  0.273   0.322
# 2 9875     376 1739. 0.216  0.195   0.239
biostat3::survRate(Surv(surv, status) ~ histo3+sex, data=d) 
#                        histo3    sex    tstop event      rate     lower     upper
# histo3=9863, sex=Female   9863 Female 926.0424   246 0.2656466 0.2334850 0.3010005
# histo3=9863, sex=Male     9863   Male 968.3970   316 0.3263125 0.2913222 0.3643479
# histo3=9875, sex=Female   9875 Female 872.8706   171 0.1959053 0.1676423 0.2275691
# histo3=9875, sex=Male     9875   Male 866.2615   205 0.2366491 0.2053611 0.2713564
(Freq <- xtabs(~histo3+sex,d))
chisq.test(Freq) #indeed, no sex diff in 63 vs 75

par(mfrow=c(1, 2))
survfit(Surv(surv, COD2=="LC") ~ histo3, data = d) |>
  plot(col=1:2,xlab = "Years since Dx",ylab = "Survival",main = "K-M of LC death")
survfit(Surv(surv, COD2=="OC") ~ histo3, data = d) |>   
  plot(col=1:2,xlab = "Years since Dx",ylab = "Survival",main = "K-M of OC death")
legend("topright", levels(d$histo3), col=1:2, lty = 1,bty="n")
# diff in LC more than OC
table(d$COD7)  # go to COD7 to zoom in on OCs
# alive   ASH    CA    CV    DK    IN    LC   YOC 
#   712    25    84   217    35    30   341   206 
par(mfrow=c(1, 2))
survfit(Surv(surv, COD7=="CV") ~ histo3, data = d) |>
  plot(col=1:2,xlab = "Years since Dx",ylab = "Survival",main = "K-M of CVD death") 
legend("topright", levels(d$histo3), col=1:2, lty = 1,bty="n")
survfit(Surv(surv, COD7=="YOC") ~ histo3, data = d) |>   
  plot(col=1:2,xlab = "Years since Dx",ylab = "Survival",main = "K-M of YOC death")
# diff in YOC more than CVD
par(mfrow=c(1, 1))
survfit(Surv(surv, COD7=="CA") ~ histo3, data = d) |>
  plot(col=1:2,xlab = "Years since Dx",ylab = "Survival",main = "K-M of CA death")
legend("topright", levels(d$histo3), col=1:2, lty = 1,bty="n")
# better care of other cancer that killed the patient correlates with molecular testing that picks up CML
dca=d|>filter(COD7=="CA")
sort(table(dca$CODS),decreasing=T)  
# Must believe the top 2 here are really deaths by CML

# In situ, benign or unknown behavior neoplasm 
#                                           18 
#               Miscellaneous Malignant Cancer 
#                                           17 
#                            Lung and Bronchus 
#                                           10 
#                         Non-Hodgkin Lymphoma 
#                                            7 
#                       Colon excluding Rectum 
#                                            5 
#                                     Pancreas 
#                                            5 
#                                       Breast 
#                                            4 
#                                    Esophagus 
#                                            3 
#                      Kidney and Renal Pelvis 
#                                            3 
#                                        Liver 
#                                            3 
#                                     Prostate 
#                                            2 
#                              Urinary Bladder 
#                                            2 
#                            Non-Melanoma Skin 
#                                            2 
#                                      Myeloma 
#                                            1 
#                                      Stomach 
#                                            1 
#                Other Oral Cavity and Pharynx 
#                                            1 
#                         Other Cause of Death 
#                                            0 
#               Aleukemic, Subleukemic and NOS 
#                                            0 

# Also suspect CT monitoring of NHL is causing some CML 
# and that lung radiation therapy is causing CML, or radiation is cause both CML and lung cancer in the same person 

dca
#          id sex    agedx  yrdx histo3 cancer  surv status   COD COD2  COD7  CODS                                        
#       <int> <fct>  <int> <dbl> <fct>  <chr>  <dbl>  <dbl> <int> <chr> <fct> <fct>                                       
#  1  5276857 Male      86  2019 9875   CML    0.611      1    11 OC    CA    Esophagus                                   
#  2  6160396 Female    81  2020 9863   CML    2.04       1    29 OC    CA    Liver                                       
#  3 10926214 Female    89  2020 9875   CML    0.214      1   130 OC    CA    In situ, benign or unknown behavior neoplasm
#  4 11430215 Female    89  2016 9863   CML    1.18       1    14 OC    CA    Colon excluding Rectum                      
#  5 11605287 Male      85  2017 9863   CML    1.35       1   130 OC    CA    In situ, benign or unknown behavior neoplasm
#  6 14688847 Male      88  2017 9875   CML    1.06       1    86 OC    CA    Miscellaneous Malignant Cancer              
#  7 16629184 Male      86  2016 9863   CML    2.43       1    86 OC    CA    Miscellaneous Malignant Cancer              
#  8 21213220 Male      83  2016 9875   CML    0.342      1    86 OC    CA    Miscellaneous Malignant Cancer              
#  9 23739529 Male      81  2022 9875   CML    0.159      1    39 OC    CA    Lung and Bronchus                           
# 10 26019051 Male      87  2016 9875   CML    2.85       1    54 OC    CA    Prostate  

d=d|>mutate(COD2=ifelse(COD==130,"LC",ifelse(COD==86,"LC",COD2)))
d=d|>mutate(COD7=ifelse(COD==130,"LC",ifelse(COD==86,"LC",COD7)))

par(mfrow=c(1, 2))
survfit(Surv(surv, COD2=="LC") ~ histo3, data = d) |>
  plot(col=1:2,xlab = "Years since Dx",ylab = "Survival",main = "K-M of LC death")
survfit(Surv(surv, COD2=="OC") ~ histo3, data = d) |>   
  plot(col=1:2,xlab = "Years since Dx",ylab = "Survival",main = "K-M of OC death")
legend("topright", levels(d$histo3), col=1:2, lty = 1,bty="n")
# diff in LC more than OC is now even more exaggerated

par(mfrow=c(1, 1))
survfit(Surv(surv, COD7=="CA") ~ histo3, data = d) |>
  plot(col=1:2,xlab = "Years since Dx",ylab = "Survival",main = "K-M of CA death")
legend("topright", levels(d$histo3), col=1:2, lty = 1,bty="n")
# difference is now much smaller
dca=d|>filter(COD7=="CA")
sort(table(dca$CODS),decreasing=T)  
                                #       Lung and Bronchus                                    Non-Hodgkin Lymphoma 
                                #                      10                                                       7 
                                #  Colon excluding Rectum                                                Pancreas 
                                #                       5                                                       5 
                                #                  Breast                                               Esophagus 
                                #                       4                                                       3 
                                # Kidney and Renal Pelvis                                                   Liver 
                                #                       3                                                       3 
                                #                Prostate                                         Urinary Bladder 
                                #                       2                                                       2 
                                #       Non-Melanoma Skin                                                 Myeloma 
                                #                       2                                                       1 
                                #                 Stomach                           Other Oral Cavity and Pharynx 
                                #                       1                                                       1 
                                #    Other Cause of Death                          Aleukemic, Subleukemic and NOS 
                                #                       0                                                       0 

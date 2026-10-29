# mkSEER_h8.R  (name of this R script)
graphics.off();rm(list=ls())#clear plots and environment
library(tidyverse)
library(CMLepi)
system.time({
  d8=seer2r("h8")
  save(d8,file="~/data/CMLepi/h8.RData")
})  # 10 secs, so comment out
load("~/data/CMLepi/h8.RData")
head(d8<-d8|>filter(histo3%in%c(9863,9875))) #14919 


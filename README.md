
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
export, changing the file names from export.txt and export.dic to cml8.txt and cml8.dic. Repeat  
for SEER12 and SEER20 (i.e. SEER21 excluding IL), calling those files cml12.txt and cml12.dic, and cml20.txt and cml20.dic.

To bring these files into R use the R package SEER2R. This package is no longer on CRAN
but it is still useful and can be installed from source via 15-year old CRAN read-only files on GitHub.
``` r
pak::pak("cran/SEER2R") #Installs fine: check yields 6 help page notes (no errors or warnings). 
library(SEER2R)
c8 = read.SeerStat("Rpacks/SEER2R/cml8.dic",UseVarLabelsInData=FALSE) 
(DICInfo1 = attr(c8, "DICInfo"))
str(c8)

```


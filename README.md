
# CMLepi

<!-- badges: start -->
[![R-CMD-check](https://github.com/radivot/CMLepi/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/radivot/CMLepi/actions/workflows/R-CMD-check.yaml)
<!-- badges: end -->

The immediate goal of CMLepi is to help people estimate chronic myeloid leukemia patient life expectancies from 
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
The create a case listing  across all cancers including: Patient_ID, Sex, Race_recode_White_Black_Other,
Year_of_diagnosis, Agerecodewith1_year_olds_and_90, Agerecodewithsingle_ages_and_90, Site_recode_ICD_O_3_WHO_2008,
Histologic_Type_ICD_O_3, Survival_days  (32765 days = 89.7 years means "Unknown"), and COD_to_site_recode. 
Convert the listing to raw code before saving as a CSV file, including a SAS file as well to document the mapping
to raw code. Repeat this for the SEER8, SEER12 and SEER20 (i.e. SEER21 excluding IL) databases, call the files s8.txt, s12.txt
and s20.txt, and place them into the folder ~/data/seer26/csvs.


``` r
library(CMLepi)
## basic example code
```


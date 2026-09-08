# lifetab2.R Application of Lab 1 in biostat3 to CML. Shows how to generate actuarial survival life tables with hazards 
library(biostat3)   #loads survival (for surv and survfit) and MASS (includes boxcox, datasets and unfortunately, select)
library(tidyverse)   #masks MASS::select, want tidy version, so first load biostat3
load("~/data/CMLepi/cml20.RData") #made in mkSEER.R
(d=d20|>filter(histo3%in%c(9863,9875),agedx<90,surv<80,surv>0)) # 43,932 CML cases
(d=d|>filter(agedx>=80)|>select(yrdx,agedx,sex,surv,status)) #4.7k
(lt=biostat3::lifetab2(Surv(surv,status)~1,data=d, breaks=0:12))
#        tstart tstop nsubs nlost  nrisk nevent       surv        pdf    hazard     se.surv      se.pdf   se.hazard
# 0-1         0     1  4726   170 4641.0   1612 1.00000000 0.34733894 0.4203390 0.000000000 0.006988996 0.010235456
# 1-2         1     2  2944   143 2872.5    602 0.65266106 0.13678049 0.2341046 0.006988996 0.005168189 0.009475802 #haz at low 0.2's
# 2-3         2     3  2199   116 2141.0    428 0.51588057 0.10312792 0.2221069 0.007421762 0.004699230 0.010669534
# 3-4         3     4  1655    91 1609.5    346 0.41275265 0.08873092 0.2408632 0.007425805 0.004517908 0.012854640
# 4-5         4     5  1218    77 1179.5    241 0.32402173 0.06620537 0.2275732 0.007200394 0.004078689 0.014564068
# 5-6         5     6   900    47  876.5    187 0.25781636 0.05500474 0.2388250 0.006877119 0.003857484 0.017339652
# 6-7         6     7   666    35  648.5    147 0.20281161 0.04597272 0.2556522 0.006480309 0.003643654 0.020912860
# 7-8         7     8   484    43  462.5    108 0.15683889 0.03662400 0.2643819 0.006019327 0.003390451 0.025216903
# 8-9         8     9   333    26  320.0     72 0.12021489 0.02704835 0.2535211 0.005550315 0.003071572 0.029636738
# 9-10        9    10   235    18  226.0     60 0.09316654 0.02473448 0.3061224 0.005135937 0.003057567 0.039054560
# 10-11      10    11   157    15  149.5     36 0.06843206 0.01647862 0.2737643 0.004660540 0.002643114 0.045197902 #haz now up to ~0.3
# 11-12      11    12   106    13   99.5     30 0.05195344 0.01566435 0.3550296 0.004271523 0.002715086 0.063789787 # so gained ~0.1 in ~10 years
# 12-Inf     12   Inf    63    20   53.0     43 0.03628909         NA        NA 0.003822963          NA          NA # likely all via aging (see below)
(ft=survfit(Surv(surv,status)~1,data=d)) 
plot(ft) #plots KM curve made above
lines(lt$tstart,lt$surv,col="red") #actuarial survival plot is right on top KM
(ft=survfit(Surv(floor(surv),status)~1,data=d)) #floor of surv makes it fall below the actuarial
plot(ft)
lines(lt$tstart,lt$surv,col="red") #with both floored actuarial curve starts out hitting tops of steps but eventually goes through them
load("~/data/mrt/mrtUSA.RData")#mrt is list of 3 matrices
D=SEERaBomb::msd(d,mrt,brkst=c(0,1,2,3,4,5,6,7,8,9,10))
D=D|>rename(Group="sex")|>select(Group,int,everything())
(D=SEERaBomb::foldD(D,keep=c("int")))
#    int          O     E    PY      t    EAR     LL    UL    RR   rrL   rrU
#    <fct>    <dbl> <dbl> <dbl>  <dbl>  <dbl>  <dbl> <dbl> <dbl> <dbl> <dbl>
#  1 (0,1]     1612 302.  3566.  0.377 0.367  0.345  0.389  5.34  5.08  5.61
#  2 (1,2]      602 232.  2535.  1.43  0.146  0.127  0.165  2.59  2.39  2.81
#  3 (2,3]      428 191.  1907.  2.43  0.124  0.103  0.145  2.24  2.03  2.46
#  4 (3,4]      347 158.  1442.  3.44  0.131  0.106  0.157  2.20  1.97  2.44
#  5 (4,5]      240 125.  1047.  4.43  0.109  0.0804 0.138  1.91  1.68  2.17
#  6 (5,6]      187 102.   783.  5.43  0.108  0.0738 0.142  1.83  1.57  2.11
#  7 (6,7]      147  83.4  574.  6.43  0.111  0.0694 0.152  1.76  1.49  2.07
#  8 (7,8]      108  65.1  404.  7.42  0.106  0.0556 0.156  1.66  1.36  2.00
#  9 (8,9]       72  49.9  281.  8.42  0.0786 0.0194 0.138  1.44  1.13  1.82
# 10 (9,10]      60  38.0  195.  9.41  0.113  0.0348 0.191  1.58  1.20  2.03
# 11 (10,100]   109  75.6  327. 11.0   0.102  0.0396 0.165  1.44  1.18  1.74
#shows EAR leveling at 0.1 by 5 years, so background haz above goes from ~0.1 to ~0.2
### let's see if we can get this out of us_mort
load("~/data/mrt/us_mort.RData") #us_mort is a class vital object (single tibble-like)  
library(vital)  
Vit=us_mort|>filter(Sex == "Total", Year == 2024)
(n85=Vit|>filter(Age>84)|>life_table())  # assume if 80-90 at Dx, ave Dx at 85
#    Year   Age Sex       mx     qx    lx     dx    Lx    Tx    ex    rx    nx    ax
#   <int> <int> <chr>  <dbl>  <dbl> <dbl>  <dbl> <dbl> <dbl> <dbl> <dbl> <dbl> <dbl>
# 1  2024    85 Total 0.0814 0.0782 1     0.0782 0.961  6.91  6.91 0.961     1   0.5
# 2  2024    86 Total 0.0913 0.0873 0.922 0.0805 0.882  5.95  6.45 0.917     1   0.5
# 3  2024    87 Total 0.102  0.0969 0.841 0.0815 0.801  5.06  6.02 0.908     1   0.5
# 4  2024    88 Total 0.115  0.109  0.760 0.0829 0.718  4.26  5.61 0.897     1   0.5
# 5  2024    89 Total 0.124  0.117  0.677 0.0791 0.637  3.55  5.24 0.887     1   0.5
# 6  2024    90 Total 0.136  0.127  0.598 0.0759 0.560  2.91  4.87 0.878     1   0.5
# 7  2024    91 Total 0.153  0.142  0.522 0.0743 0.485  2.35  4.50 0.866     1   0.5
# 8  2024    92 Total 0.172  0.159  0.448 0.0711 0.412  1.86  4.16 0.850     1   0.5
# 9  2024    93 Total 0.193  0.176  0.377 0.0664 0.343  1.45  3.86 0.833     1   0.5
#10  2024    94 Total 0.215  0.194  0.310 0.0602 0.280  1.11  3.57 0.816     1   0.5
# Indeed, background h = mx or ~qx, goes up from ~0.1 to ~0.2 over these 10 years 

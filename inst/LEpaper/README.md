These are the scripts that we used to create the figures (e.g. N_XXX.R makes Figure N) in
**Diagnosis of Chronic Myeloid Leukemia at a Mean Age of 59 Implies a Mean Loss of Roughly 8 of 24 Years of Life in the United States**.
Before running thesse scripts, first run the ones on the main [page](https://github.com/radivot/CMLepi). Then run `mkMorts.R`. 
It should then be possible to run the other scripts. Supplementary figure script names begin with an S.  Note that
`S7_mkhBack.r` must be run before `S7_RP.R`. 

Figure 7 scripts are difficult to grasp because they contain inner workings of  
`SEERaBomb::msd()` (Mortality Since Diagnosis) never meant to  be exposed. These scripts are first attempts at extending `msd()` to
handle specific causes of death.  Work toward simplifying `msd()` and extending it to specific causes are in the adjacent folder `msd`.



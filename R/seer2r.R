#' Convert db.txt db.dic files generated in a SEER*stat Case Listing session into an R tibble
#'
#'
#'
#' @param db database name.
#' @param path folder name.
#' @returns A tibble.
#' @importFrom SEER2R read.SeerStat
#' @importFrom dplyr mutate rename filter select bind_cols as_tibble relocate
#' @importFrom stringr str_sub
#' @importFrom forcats as_factor
#' @export
seer2r=function(db,path="~/data/CMLepi") {
  # pak::pak("cran/SEER2R") #installs fine from  GitHub even though no longer on CRAN
  # library(SEER2R)

  # gimmick to get rid of unwanted notes in R CMD check
  id=sex=agedx=histo3=surv=yrdx=cancer=COD=CODS=COD2=COD7=histS=site=Year=Age=ICDO3=mapCOD7=NULL

  (inF=file.path(path,paste0(db,".dic")))
  # n = SEER2R::read.SeerStat(inF,UseVarLabelsInData=FALSE) #get numbers(n)
  n = read.SeerStat(inF,UseVarLabelsInData=FALSE) #get numbers(n)
  n<-attr(n,"assignColNames")(n,c("id","sex","Age","Year","ICDO3","surv","COD"))
  (n=n|>rename(yrdx=Year,agedx=Age)|>as_tibble())
  n=n|>mutate(yrdx=yrdx+1800,surv=surv/365.25)
  n=n|>mutate(status=as.numeric(COD>0),.after=surv)
  n=n|>mutate(sex=ifelse(sex==1,"Male","Female"))
  (n=n|>mutate(sex=as_factor(sex)))

  # c = SEER2R::read.SeerStat(inF,UseVarLabelsInData=TRUE)
  c = read.SeerStat(inF,UseVarLabelsInData=TRUE)
  c<-attr(c,"getSubDataByVarName")(c,c("ICDO3","site"))
  (c=c|>rename(histS=ICDO3,CODS=site)|>as_tibble())
  (d=bind_cols(n,c))
  (d=d|>mutate(histo3=as.numeric(str_sub(histS,end=4)),.after=ICDO3))
  d=d|>filter(histo3%in%c(9863,9875,9945)) #leave out rare jCMML=9946 and atypical CML=9876
  (d=d|>mutate(cancer=ifelse(histo3%in%c(9863,9875),"CML","CMML"),.after=histo3))
  d=d|>select(-histS)
  d=d|>mutate(COD2=ifelse(COD==0,"alive",ifelse((COD>=74)&(COD<=85)|(COD==89),"LC","OC")),.after=COD)
  d=d|>mutate(CODS=as_factor(CODS)) #to save a little memory

  mapCOD7=function(D){
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
    D$COD7=as.factor(CODS)
    D|>relocate(COD7, .after = COD2)
  }
  mapCOD7(d)
}

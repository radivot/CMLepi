#' Makes generalized additive models for causes of death that were aggregated into 13 groups
#'
#'
#'
#' @param seerHome folder name.
#' @param inFile input file.
#' @param outFile output file.
#' @returns Nothing is returned. Run for the side effect of creating a list of gam fits in the output file.
#' @importFrom dplyr mutate rename filter select relocate last_col summarize group_by
#' @importFrom forcats as_factor
#' @importFrom mgcv gam
#' @importFrom stats poisson
#' @export
mkG13<-function(seerHome="~/data/CMLepi",
               inFile="seerMrt.RData",
               outFile="G13.RData"){
  # # Use mkSEERmrt() to make inFile
  # seerHome="~/data/CMLepi"
  # inFile="seerMrt.RData"
  # outFile="G12.RData"
  # require(dplyr)
  # require(forcats)
  # require(mgcv)
  d=COD=num=denom=year=sex=Ages=age=ageMid=a20=rate=rate2=W=NULL
  load(file.path(seerHome,inFile)) # get d with dims 336k × 7
  (d=d|>select(-rate,-rate2,-a20,-ageMid)) #drop rates since pois reg as in mkG2
  (d=d|>relocate(c(denom,num), .after = last_col())) # get num ready to sum over CODs grouped together
  # show how total of 120M deaths breaks down, and how 1M leukemic deaths distribute
  # d|>filter(sex=="Both")|>group_by(COD)%>%summarise(num=sum(num),.groups="drop")|>mutate(num=round(num/1e6,3))|>print(n=112)
#    2   0.437  #CA oral
#   13   6.66   #CA digestive
#   30   7.15   #CA resp [i.e. lung]  # numbers of deaths here are in millions
#   36   0.07   #CA bones
#   37   0.187  #CA soft tissue
#   38   0.491  #CA skin
#   41   2.05   #CA breast
#   42   1.34   #CA Female
#   51   1.49   #CA Male
#   56   1.23   #CA urinary
#   61   0.015  #CA eye
#   62   0.649  #CA brain
#   63   0.114  #CA endocrine
#   66   1.01   #CA lymphoma
#   69   0.491  #CA multiple myeloma
  L=NULL
  tits=NULL
  dCA=d|>filter(COD%in%c(2,13,30,36:38,41:42,51,56,61:63,66,69))|>mutate(COD="CA")
  (L[["CA"]]=dCA|>group_by(COD,year,age,sex,denom)|>summarize(num=sum(num),.groups="drop"))
  tits["CA"]="Cancer"
#   70   1.02  #LC  leukemic causes (includes CML zoom outs to other leukemias, BP => AML, & getting only chronic right => CLL
  dLC=d|>filter(COD==70)|>mutate(COD="LC")
  (dLC=dLC|>group_by(COD,year,age,sex,denom)|>summarize(num=sum(num),.groups="drop")) # need this to order it same as below
  (L[["LC"]]=dLC)
  tits["LC"]="Leukemic Cause"
#   83   1.94 #MCA  Miscellaneous malignant CAncer (could include zoom outs of CML to some cancer)
  dMCA=d|>filter(COD==83)|>mutate(COD="MCA")
  (L[["MCA"]]=dMCA|>group_by(COD,year,age,sex,denom)|>summarize(num=sum(num),.groups="drop"))
  tits["MCA"]="Misc CAncer"
#   84   0.552 #BEN  Benign (also could include zoom outs of CML, if someone thinks it is pretty benign)
  dBEN=d|>filter(COD==84)|>mutate(COD="BEN")
  (L[["BEN"]]=dBEN|>group_by(COD,year,age,sex,denom)|>summarize(num=sum(num),.groups="drop"))
  tits["BEN"]="BEN"
#   85   0.06   #IN  TB
#   86   0.004  #IN  syphillus
#   87   0.554  #IN  HIV
#   88   1.34   #IN  sepsis
#   89   0.703  #IN  parasites
  (dIN=d|>filter(COD%in%c(85:89,98))|>mutate(COD="IN")) #INfections TB, syph, hiv, septicemia, parasites 98=flu
  (L[["IN"]]=dIN|>group_by(COD,year,age,sex,denom)|>summarize(num=sum(num),.groups="drop")) # 21k => 3k
  tits["IN"]="Infection"
#   90   3.12 #DK  Diabetes Mellitus
  (dDK=d|>filter(COD%in%c(90,102))|>mutate(COD="DK"))  #90=DM,  # 102 = kidney disease
  (L[["DK"]]=dDK|>group_by(COD,year,age,sex,denom)|>summarize(num=sum(num),.groups="drop"))
  tits["DK"]="DM & Kidney Disease"
#   91   2.58      #MSPC    Alzheimers
#   92  34.8       #CV  Heart Disease
#   93   0.99      #CV  Hypertension without Heart Disease
#   94   7.67      #CV  Cerebrovascular Diseases
#   95   0.749     #CV  Atherosclerosis
#   96   0.668     #CV  Aortic Aneurysm
#   97   0.464     #CV  Other Diseases of Arteries, Arterioles, Capillaries
  (dCV=d|>filter(COD%in%c(92:97))|>mutate(COD="CV"))
  (L[["CV"]]=dCV|>group_by(COD,year,age,sex,denom)|>summarize(num=sum(num),.groups="drop")) # 18k to 3k
  tits["CV"]="Cardio and Vascular"
#   98   3.09   #IN     flu
#   99   5.49   #COPD   Chronic Obstructive Pulmonary Disease and Allied Conditions
  dCOPD=d|>filter(COD==99)|>mutate(COD="COPD")
  (L[["COPD"]]=dCOPD|>group_by(COD,year,age,sex,denom)|>summarize(num=sum(num),.groups="drop"))
  tits["COPD"]="COPD"
#  100   0.238  #MSPC   stomach ulcers
#  101   1.61   #LIV   liver disease
  dLIV=d|>filter(COD==101)|>mutate(COD="LIV")
  (L[["LIV"]]=dLIV|>group_by(COD,year,age,sex,denom)|>summarize(num=sum(num),.groups="drop"))
  tits["LIV"]="LIVer disease"
#  102   1.72   #DK    kidney disease
#  103   0.025  #MSPC   Complications of birth
#  104   0.568  #MSPC   congenital conditions
#  105   0.77   #MSPC   perinatl condiditions
  dMSPC=d|>filter(COD%in%c(91,100,103:105))|>mutate(COD="MSPC")
  (L[["MSPC"]]=dMSPC|>group_by(COD,year,age,sex,denom)|>summarize(num=sum(num),.groups="drop"))
  tits["MSPC"]="Miscellaneous SPecific Causes"
#  106   1.56   #ILL   ill defined (i.e. could include zoom out of CML)
  dILL=d|>filter(COD==106)|>mutate(COD="ILL")
  (L[["ILL"]]=dILL|>group_by(COD,year,age,sex,denom)|>summarize(num=sum(num),.groups="drop"))
  tits["ILL"]="Ill-defined COD"
#  107   6.07   #ASH   accidents
#  108   1.74   #ASH   suicides
#  109   1.03   #ASH   homocides
  dASH=d|>filter(COD%in%c(107:109))|>mutate(COD="ASH") #accidents, suicides and homocides
  (L[["ASH"]]=dASH|>group_by(COD,year,age,sex,denom)|>summarize(num=sum(num),.groups="drop"))
  tits["ASH"]="Accidents, Suicides and Homocides"
#  110   1.04  #OCD  covid
#  111  14.7   #OCD Other Cause of Death (i.e. could include zoom out of CML)
  dOCD=d|>filter(COD%in%c(100,111))|>mutate(COD="OCD")
  (L[["OCD"]]=dOCD|>group_by(COD,year,age,sex,denom)|>summarize(num=sum(num),.groups="drop"))
  tits["OCD"]="Other Cause of Death"
  (nms=names(L))
  G=vector(mode="list",length=length(L))
  names(G)<-nms
  for (i in nms) {
    # i="LC"
    print(i)
    D=L[[i]]
    (D=D|>filter(age>20))
    D$sex=as_factor(D$sex)
    (Df=D|>filter(year<2020)) # data for fitting
    print(summary(G[[i]]<-mgcv::gam(num ~ sex+s(age,year,by=sex)+ti(age,year)+offset(log(denom)),family=poisson(),data=Df)))
  }
  save(G,L,tits,file=file.path(seerHome,outFile))
}

# WARNING: the 13 Mortality Data based COD groups defined above need to be synced up with 13 incidence data COD13 defs

# US mort defs
#  0 = "All Causes of Death"
#  1 = "  All Malignant Cancers"   # all cancers
#  2 = "    Oral Cavity and Pharynx"
#  3 = "      Lip"
#  4 = "      Tongue"
#  5 = "      Salivary Gland"
#  6 = "      Floor of Mouth"
#  7 = "      Gum and Other Mouth"
#  8 = "      Nasopharynx"
#  9 = "      Tonsil"
# 10 = "      Oropharynx"
# 11 = "      Hypopharynx"
# 12 = "      Other Oral Cavity and Pharynx"
# 13 = "    Digestive System"
# 14 = "      Esophagus"
# 15 = "      Stomach"
# 16 = "      Small Intestine"
# 17 = "      Colon and Rectum"
# 18 = "        Colon excluding Rectum"
# 19 = "        Rectum and Rectosigmoid Junction"
# 20 = "      Anus, Anal Canal and Anorectum"
# 21 = "      Liver and Intrahepatic Bile Duct"
# 22 = "        Liver"
# 23 = "        Intrahepatic Bile Duct"
# 24 = "      Gallbladder"
# 25 = "      Other Biliary"
# 26 = "      Pancreas"
# 27 = "      Retroperitoneum"
# 28 = "      Peritoneum, Omentum and Mesentery"
# 29 = "      Other Digestive Organs"
# 30 = "    Respiratory System"
# 31 = "      Nose, Nasal Cavity and Middle Ear"
# 32 = "      Larynx"
# 33 = "      Lung and Bronchus"
# 34 = "      Pleura"
# 35 = "      Trachea, Mediastinum and Other Respiratory Organs"
# 36 = "    Bones and Joints"
# 37 = "    Soft Tissue including Heart"
# 38 = "    Skin"
# 39 = "      Melanoma of the Skin"
# 40 = "      Non-Melanoma Skin"
# 41 = "    Breast"
# 42 = "    Female Genital System"
# 43 = "      Cervix Uteri"
# 44 = "      Corpus and Uterus, NOS"
# 45 = "        Corpus Uteri"
# 46 = "        Uterus, NOS"
# 47 = "      Ovary"
# 48 = "      Vagina"
# 49 = "      Vulva"
# 50 = "      Other Female Genital Organs"
# 51 = "    Male Genital System"
# 52 = "      Prostate"
# 53 = "      Testis"
# 54 = "      Penis"
# 55 = "      Other Male Genital Organs"
# 56 = "    Urinary System"
# 57 = "      Urinary Bladder"
# 58 = "      Kidney and Renal Pelvis"
# 59 = "      Ureter"
# 60 = "      Other Urinary Organs"
# 61 = "    Eye and Orbit"
# 62 = "    Brain and Other Nervous System"
# 63 = "    Endocrine System"
# 64 = "      Thyroid"
# 65 = "      Other Endocrine including Thymus"
# 66 = "    Lymphoma"
# 67 = "      Hodgkin Lymphoma"
# 68 = "      Non-Hodgkin Lymphoma"
# 69 = "    Myeloma"
# 70 = "    Leukemia"
# 71 = "      Lymphocytic Leukemia"
# 72 = "        Acute Lymphocytic Leukemia"
# 73 = "        Chronic Lymphocytic Leukemia"
# 74 = "        Other Lymphocytic Leukemia"
# 75 = "      Myeloid and Monocytic Leukemia"
# 76 = "        Acute Myeloid Leukemia"
# 77 = "        Acute Monocytic Leukemia"
# 78 = "        Chronic Myeloid Leukemia"
# 79 = "        Other Myeloid/Monocytic Leukemia"
# 80 = "      Other Leukemia"
# 81 = "        Other Acute Leukemia"
# 82 = "        Aleukemic, Subleukemic and NOS"
# 83 = "    Miscellaneous Malignant Cancer"
# 84 = "  In situ, benign or unknown behavior neoplasm"
# 85 = "  Tuberculosis"
# 86 = "  Syphilis"
# 87 = "  Human Immunodeficiency Virus (HIV) (1987+)"
# 88 = "  Septicemia"
# 89 = "  Other Infectious and Parasitic Diseases"
# 90 = "  Diabetes Mellitus"
# 91 = "  Alzheimers (ICD-9 and 10 only)"
# 92 = "  Diseases of Heart"
# 93 = "  Hypertension without Heart Disease"
# 94 = "  Cerebrovascular Diseases"
# 95 = "  Atherosclerosis"
# 96 = "  Aortic Aneurysm and Dissection"
# 97 = "  Other Diseases of Arteries, Arterioles, Capillaries"
# 98 = "  Pneumonia and Influenza"
# 99 = "  Chronic Obstructive Pulmonary Disease and Allied Cond"
#100 = "  Stomach and Duodenal Ulcers"
#101 = "  Chronic Liver Disease and Cirrhosis"
#102 = "  Nephritis, Nephrotic Syndrome and Nephrosis"
#103 = "  Complications of Pregnancy, Childbirth, Puerperium"
#104 = "  Congenital Anomalies"
#105 = "  Certain Conditions Originating in Perinatal Period"
#106 = "  Symptoms, Signs and Ill-Defined Conditions"
#107 = "  Accidents and Adverse Effects"
#108 = "  Suicide and Self-Inflicted Injury"
#109 = "  Homicide and Legal Intervention"
#110 = "  COVID-19 (2020+)"
#111 = "  Other Cause of Death"

##### so pairing  86 below with  83 above
#####     and    130 below with  84 above
#####     and    208 below with 111 above

### but where does 110=COVID go below, guess 145
### and where does 252 below (no DC) go above? guess 111

# [Format=COD to site recode]
# 0=Alive
# 1=Lip
# 2=Tongue
# 3=Salivary Gland
# 4=Floor of Mouth
# 5=Gum and Other Mouth
# 6=Nasopharynx
# 7=Tonsil
# 8=Oropharynx
# 9=Hypopharynx
# 10=Other Oral Cavity and Pharynx
# 11=Esophagus
# 12=Stomach
# 13=Small Intestine
# 14=Colon excluding Rectum
# 24=Rectum and Rectosigmoid Junction
# 27=Anus, Anal Canal and Anorectum
# 29=Liver
# 30=Intrahepatic Bile Duct
# 31=Gallbladder
# 32=Other Biliary
# 33=Pancreas
# 34=Retroperitoneum
# 35=Peritoneum, Omentum and Mesentery
# 36=Other Digestive Organs
# 37=Nose, Nasal Cavity and Middle Ear
# 38=Larynx
# 39=Lung and Bronchus
# 40=Pleura
# 41=Trachea, Mediastinum and Other Respiratory Organs
# 42=Bones and Joints
# 43=Soft Tissue including Heart
# 44=Melanoma of the Skin
# 45=Non-Melanoma Skin
# 46=Breast
# 47=Cervix Uteri
# 48=Corpus Uteri
# 49=Uterus, NOS
# 50=Ovary
# 51=Vagina
# 52=Vulva
# 53=Other Female Genital Organs
# 54=Prostate
# 55=Testis
# 56=Penis
# 57=Other Male Genital Organs
# 58=Urinary Bladder
# 59=Kidney and Renal Pelvis
# 60=Ureter
# 61=Other Urinary Organs
# 62=Eye and Orbit
# 90=Brain and Other Nervous System
# 65=Thyroid
# 66=Other Endocrine including Thymus
# 67=Hodgkin Lymphoma
# 70=Non-Hodgkin Lymphoma
# 73=Myeloma
# 74=Acute Lymphocytic Leukemia
# 75=Chronic Lymphocytic Leukemia
# 76=Other Lymphocytic Leukemia
# 77=Acute Myeloid Leukemia
# 80=Acute Monocytic Leukemia
# 78=Chronic Myeloid Leukemia
# 89=Other Myeloid/Monocytic Leukemia
# 83=Other Acute Leukemia
# 85=Aleukemic, Subleukemic and NOS
# 86=Miscellaneous Malignant Cancer
# 130=In situ, benign or unknown behavior neoplasm
# 133=Tuberculosis
# 136=Syphilis
# 142=Septicemia
# 145=Other Infectious and Parasitic Diseases including HIV
# 148=Diabetes Mellitus
# 151=Alzheimers (ICD-9 and 10 only)
# 154=Diseases of Heart
# 157=Hypertension without Heart Disease
# 160=Cerebrovascular Diseases
# 163=Atherosclerosis
# 166=Aortic Aneurysm and Dissection
# 169=Other Diseases of Arteries, Arterioles, Capillaries
# 172=Pneumonia and Influenza
# 175=Chronic Obstructive Pulmonary Disease and Allied Cond
# 178=Stomach and Duodenal Ulcers
# 181=Chronic Liver Disease and Cirrhosis
# 184=Nephritis, Nephrotic Syndrome and Nephrosis
# 187=Complications of Pregnancy, Childbirth, Puerperium
# 190=Congenital Anomalies
# 193=Certain Conditions Originating in Perinatal Period
# 196=Symptoms, Signs and Ill-Defined Conditions
# 199=Accidents and Adverse Effects
# 202=Suicide and Self-Inflicted Injury
# 205=Homicide and Legal Intervention
# 208=Other Cause of Death
# 252=State DC not available or state DC available but no COD


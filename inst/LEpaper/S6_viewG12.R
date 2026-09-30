seerHome="~/data/CMLepi"
inFile="seerMrt.RData"
library(dplyr)
library(forcats)
library(mgcv)
library(rgl)  #or call funcs with ::
options(rgl.useNULL = TRUE)
options(rgl.printRglwidget = TRUE)
# load(file.path(seerHome,inFile)) # get d with dims 336k × 7
# (d=d|>select(-rate,-rate2,-a20,-ageMid)) #drop rates since pois reg as in mkG2
# (d=d|>relocate(c(denom,num), .after = last_col())) # get num ready to sum over CODs grouped together
# L=NULL
# tits=NULL
# dCA=d|>filter(COD%in%c(2,13,30,36:38,41:42,51,56,61:63,66,69))|>mutate(COD="CA")
# (L[["CA"]]=dCA|>group_by(COD,year,age,sex,denom)|>summarize(num=sum(num),.groups="drop"))
# tits["CA"]="Cancer"
# #   70   1.02  #LC  leukemic causes (includes CML zoom outs to other leukemias, BP => AML, & getting only chronic right => CLL
# dLC=d|>filter(COD==70)|>mutate(COD="LC")
# (dLC=dLC|>group_by(COD,year,age,sex,denom)|>summarize(num=sum(num),.groups="drop")) # need this to order it same as below
# (L[["LC"]]=dLC)
# tits["LC"]="Leukemic Cause"
# #   83   1.94 #MCA  Miscellaneous malignant CAncer (could include zoom outs of CML to some cancer)
# dMCA=d|>filter(COD==83)|>mutate(COD="MCA")
# (L[["MCA"]]=dMCA|>group_by(COD,year,age,sex,denom)|>summarize(num=sum(num),.groups="drop"))
# tits["MCA"]="Misc CAncer"
# #   84   0.552 #BEN  Benign (also could include zoom outs of CML, if someone thinks it is pretty benign)
# dBEN=d|>filter(COD==84)|>mutate(COD="BEN")
# (L[["BEN"]]=dBEN|>group_by(COD,year,age,sex,denom)|>summarize(num=sum(num),.groups="drop"))
# tits["BEN"]="BEN"
# #   85   0.06   #IN  TB
# #   86   0.004  #IN  syphillus
# #   87   0.554  #IN  HIV
# #   88   1.34   #IN  sepsis
# #   89   0.703  #IN  parasites
# (dIN=d|>filter(COD%in%c(85:89,98))|>mutate(COD="IN")) #INfections TB, syph, hiv, septicemia, parasites 98=flu
# (L[["IN"]]=dIN|>group_by(COD,year,age,sex,denom)|>summarize(num=sum(num),.groups="drop")) # 21k => 3k
# tits["IN"]="Infection"
# #   90   3.12 #DK  Diabetes Mellitus
# (dDK=d|>filter(COD%in%c(90,102))|>mutate(COD="DK"))  #90=DM,  # 102 = kidney disease
# (L[["DK"]]=dDK|>group_by(COD,year,age,sex,denom)|>summarize(num=sum(num),.groups="drop"))
# tits["DK"]="DM & Kidney Disease"
# #   91   2.58      #MSPC    Alzheimers
# #   92  34.8       #CV  Heart Disease
# #   93   0.99      #CV  Hypertension without Heart Disease
# #   94   7.67      #CV  Cerebrovascular Diseases
# #   95   0.749     #CV  Atherosclerosis
# #   96   0.668     #CV  Aortic Aneurysm
# #   97   0.464     #CV  Other Diseases of Arteries, Arterioles, Capillaries
# (dCV=d|>filter(COD%in%c(92:97))|>mutate(COD="CV"))
# (L[["CV"]]=dCV|>group_by(COD,year,age,sex,denom)|>summarize(num=sum(num),.groups="drop")) # 18k to 3k
# tits["CV"]="Cardio and Vascular"
# #   98   3.09   #IN     flu
# #   99   5.49   #COPD   Chronic Obstructive Pulmonary Disease and Allied Conditions
# dCOPD=d|>filter(COD==99)|>mutate(COD="COPD")
# (L[["COPD"]]=dCOPD|>group_by(COD,year,age,sex,denom)|>summarize(num=sum(num),.groups="drop"))
# tits["COPD"]="COPD"
# #  100   0.238  #MSPC   stomach ulcers
# #  101   1.61   #MSPC   liver disease
# #  102   1.72   #DK    kidney disease
# #  103   0.025  #MSPC   Complications of birth
# #  104   0.568  #MSPC   congenital conditions
# #  105   0.77   #MSPC   perinatl condiditions
# dMSPC=d|>filter(COD%in%c(91,100:101,103:105))|>mutate(COD="MSPC")
# (L[["MSPC"]]=dMSPC|>group_by(COD,year,age,sex,denom)|>summarize(num=sum(num),.groups="drop"))
# tits["MSPC"]="Miscellaneous SPecific Causes"
# #  106   1.56   #ILL   ill defined (i.e. could include zoom out of CML)
# dILL=d|>filter(COD==106)|>mutate(COD="ILL")
# (L[["ILL"]]=dILL|>group_by(COD,year,age,sex,denom)|>summarize(num=sum(num),.groups="drop"))
# tits["ILL"]="Ill-defined COD"
# #  107   6.07   #ASH   accidents
# #  108   1.74   #ASH   suicides
# #  109   1.03   #ASH   homocides
# dASH=d|>filter(COD%in%c(107:109))|>mutate(COD="ASH") #accidents, suicides and homocides
# (L[["ASH"]]=dASH|>group_by(COD,year,age,sex,denom)|>summarize(num=sum(num),.groups="drop"))
# tits["ASH"]="Accidents, Suicides and Homocides"
# #  110   1.04  #OCD  covid
# #  111  14.7   #OCD Other Cause of Death (i.e. could include zoom out of CML)
# dOCD=d|>filter(COD%in%c(110,111))|>mutate(COD="OCD")
# (L[["OCD"]]=dOCD|>group_by(COD,year,age,sex,denom)|>summarize(num=sum(num),.groups="drop"))
# tits["OCD"]="Other Cause of Death"
load("~/data/CMLepi/G12.RData") # made in mkMorts.R via CMLepi::mkG12()
(nms=names(L))  #bin above has all three of these, i.e. L, G and tits
(nms=names(G)) # "CA"   "LC"   "MCA"  "BEN"  "IN"   "DK"   "CV"   "COPD" "MSPC" "ILL"  "ASH"  "OCD" 
tits
######## new rgl problem on mac doesn't let this work in a for loop anymore
i="LC"  # very slight covid impact
i="CA"  
i="MCA"  
i="BEN"  
i="IN"  
i="DK"  
i="CV"  
i="COPD"  
i="MSPC"  
i="ASH"  # also big ridge
i="OCD"  # big ridge, so lots of covid deaths in here
i="ILL"  
print(i)
D=L[[i]]
D  # scroll through to see that data goes into low 90's with ages pegged by mkAges
# print(summary(G[[i]]<-mgcv::gam(num ~ sex+s(age,year,by=sex)+ti(age,year)+
#                 s(age, by = is_2020, bs = "cr") +
#                 s(age, by = is_2021, bs = "cr") + # age-specific extra
#                 s(age, by = is_2022, bs = "cr") + # for each covid year
#                 offset(log(denom)),family=poisson(),data=D)))
# D$E=exp(predict(G[[i]],D))  # with rgl not working right, code below is no longer a useful check
# head(D<-D%>%mutate(Eincid=E/denom))
## plot the data used in the fit
head(D<-D%>%mutate(incid=num/denom))
head(D<-D%>%mutate(incid=if_else(incid==0,0.0001,incid)))
with(D,plot3d(year,age,log10(incid),col=ifelse(sex=="Female","red","blue"),
              xlab="",ylab="Age",zlab="log10(M)",alpha=1,size=4))
#### set up fine mesh to plot the fit
(Ages=seq(min(D$age),max(D$age)))
(Ages=20:95)
(Years=seq(min(D$year),max(D$year))) # fine at 1975 to 2024
head(nD<-expand.grid(Ages,Years))
names(nD)<-c("age","year")
nD$denom=1
nDf=nD|>mutate(sex="Female")
nDm=nD|>mutate(sex="Male")
nD=bind_rows(nDf,nDm)
nD$sex=as_factor(nD$sex)
head(nD)
dim(nD)  #7600 = 2*76*50
nD$is_2020 <- as.numeric(nD$year == 2020)
nD$is_2021 <- as.numeric(nD$year == 2021)
nD$is_2022 <- as.numeric(nD$year == 2022)
nD$E=exp(predict(G[[i]],nD))
M=reshape2::acast(nD%>%filter(sex=="Female")%>%select(year,age,E), year~age, value.var="E")
surface3d(Years,Ages,log10(M),col="red",alpha=0.4) # M for Matrix
M=reshape2::acast(nD%>%filter(sex=="Male")%>%select(year,age,E), year~age, value.var="E")
surface3d(Years,Ages,log10(M),col="blue",alpha=0.4)
light3d(theta = 0, phi = 75)
bgplot3d({
  plot.new()
  title(main = paste0(i), line = 0,cex.main=1)
})

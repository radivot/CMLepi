# 2_MortalityDataC.R  #### shows death rate drops in 80's being small 
library(tidyverse)
nms=c("COD","a20","year","rate","num","denom")
d=read_tsv("~/data/CMLepi/rates.txt", col_names=nms) 
d=d|>mutate(age=ifelse(a20==19,92.5,ifelse(a20==0,0.5,ifelse(a20==1,3,(a20-1)*5+2.5))))
d=d|>mutate(year=year+1968)
d=d|>filter(year!=1968,a20!=20)  #1968 is the sum over all years 
d=d|>filter(COD==78) # 78="        Chronic Myeloid Leukemia" (see rates.dic file)
(d=d|>filter(age>=30))
(d=d|>filter(year>=1999))
(D=d|>mutate(Age=cut(age,c(30,60,80,100)))|>group_by(year,Age)|>summarize(n=sum(num),d=sum(denom)))
dput(levels(D$Age))
(ord=c("(30,60]", "(60,80]", "(80,100]")[3:1])
labs=c("30-59", "60-79", "80-99")[3:1]
D=D|>mutate(incid=n/d,Age=factor(Age,levels=ord,labels=labs))
dput(levels(D$Age))
lh=theme(legend.direction="horizontal")
D|>ggplot(aes(x=year,y=incid,col=Age))+geom_line()+
  labs(y="CML Mortality Rate",x="Year",col="Age at Death")+scale_y_log10() +
  scale_x_continuous(breaks=c(1999,2010,2024))+theme_classic(base_size=14)+
  theme(
    legend.key.size = unit(0.5, "cm"),
    legend.position = c(0.70,0.3),
    legend.text = element_text(size = 9),
    legend.title = element_text(size = 9, margin = margin(b = 3, unit = "pt")),
    legend.key.spacing.y = unit(0.0, "cm"),
    legend.margin = margin(t = 0, r = 0, b = 0, l = 0, unit = "pt")
  )
ggsave("LE/outs/2C_noHelpAt80.pdf",width=3,height=3)
ggsave("LE/outs/2C_noHelpAt80.png",width=3,height=3)



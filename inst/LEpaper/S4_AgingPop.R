# S4_AgingPop.R
library(dplyr)
library(ggplot2)
load("~/data/mrt/us_mort.RData")
(v=us_mort|>filter(Sex == "Total",Year%in%c(2007,2024)))
library(scales)
gv=geom_vline(xintercept=c(43,60),col=c(alpha("#00BFC4",0.4),alpha("#F8766D",0.4)))
v|>tibble()|>mutate(n=Population/1e6)|>mutate(Year=factor(Year,levels=c("2024","2007")))|>
  ggplot(aes(x=Age,y=n,col=Year,group=Year))+gv+geom_line()+labs(y="Number of people (Millions)",x="Age")+
  scale_x_continuous(breaks=c(0,20,43,60,80,100))+
  theme_classic(base_size=14)+
  theme(
    legend.key.size = unit(0.4, "cm"),
    legend.position = c(0.8,0.8),
    legend.text = element_text(size = 10),
    legend.title = element_text(size = 10, margin = margin(b = 2, unit = "pt")),
    legend.key.spacing.y = unit(0.0, "cm"),
    legend.margin = margin(t = 0, r = 0, b = 0, l = 0, unit = "pt")
  )
ggsave("LE/outs/supp/S4_agingPop.pdf",width=4,height=3)


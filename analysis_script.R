####Analysis of Clio channel and knight inlet data collected summer 2024#
#Last updated August 15, 2024
library("tidyverse")

data <- read.csv("~/Documents/GitHub/clio/clio.csv", skip=42)

data$Quantity <- as.numeric(data$Quantity) #Make quantities numeric
data<-data %>%
  mutate(site=case_when(str_detect(Sample.Name, "CL")~"clio",
                        str_detect(Sample.Name, "KN")~"knight")) 

boxplot <- ggplot(data=data, aes(x=site,y=Quantity, fill=Target.Name))+
  geom_boxplot()
boxplot

scatter <- ggplot(data=data, aes(x=site, y=Quantity, color=Target.Name))+
  geom_point()+
  geom_jitter()
scatter

#Dirty anova on Atlantic salmon
sasa <- subset(data, data$Target.Name=="sasa")
test.sasa <- aov(data=sasa, Quantity~site) 
#Significant difference between sites

#Dirty anova on Tmar
tmar <- subset(data, data$Target.Name=="Te_mar")
test.tmar <- aov(data=tmar, Quantity~site) 
#no significant difference

#Dirty anova on Tfin
tfin <- subset(data, data$Target.Name=="Te_fin")
test.tfin <- aov(data=tfin, Quantity~site) 
#no significant difference

#Data preparation

library("tidyverse")

data1 <- read.csv("~/Documents/GitHub/clio/clio.csv", skip=42)
data2 <- read.csv("~/Documents/GitHub/clio/clio2.csv", skip=42)
data3 <- read.csv("~/Documents/GitHub/clio/clio3.csv", skip=42)

data1$Quantity <- as.numeric(data1$Quantity) #Make quantities numeric
data2$Quantity <- as.numeric(data2$Quantity) #Make quantities numeric
data3$Quantity <- as.numeric(data3$Quantity) #Make quantities numeric
data1$Quantity.SD <- as.numeric(data1$Quantity.SD) #Make quantity SDs numeric
data2$Quantity.SD <- as.numeric(data2$Quantity.SD) #Make quantity SDs numeric
data3$Quantity.SD <- as.numeric(data3$Quantity.SD) #Make quantity SDs numeric

data <- bind_rows(data1, data2, data3)


data<-data %>%
  mutate(site=case_when(str_detect(Sample.Name, "CL")~"clio",
                        str_detect(Sample.Name, "KN")~"knight",
                        TRUE~'control')) %>%
  filter(data$Task != "STANDARD") #remove standards
NAs <- which(is.na(data$Quantity))
data$Quantity[NAs] <- 0
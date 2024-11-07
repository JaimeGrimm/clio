####Analysis of Clio channel and knight inlet data collected summer 2024#
#Last updated November 6, 2024
library("tidyverse")
library("MASS")
library("vcdExtra")  
library("boot")


####Data preparation####
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

####Visualization####
boxplot <- ggplot(data=data, aes(x=site,y=log(Quantity+1), fill=Target.Name))+
  geom_boxplot()+
  theme_bw()+
  theme(axis.title.x = element_text(size=18),
  axis.title.y = element_text(size=18),
  legend.text = element_text(size=18),
  axis.text = element_text(size=16))+
  ylab("eDNA concentration (copies/μL)")+
  scale_x_discrete(limits=c("clio", "knight", "control"))+
  scale_fill_discrete(name="Species")#, labels=c("Atlantic salmon calicivirus", "Cutthroat trout virus-2",
                                                #"Erythrocytic necrosis virus","P. salmonis", "Atlantic salmon calcivirus","Atlantic salmon", 
                                               #"Candidatus S. salmonis", "T. finnmarkense", "T. maritimum"))

boxplot

scatter <- ggplot(data=data, aes(x=site, y=Quantity, color=Target.Name))+
  geom_point()+
  geom_jitter()+
  theme(axis.title.x = element_text(size=18),
        axis.title.y = element_text(size=18),
        legend.text = element_text(size=18),
        axis.text = element_text(size=20))+
  ylab("eDNA concentration (copies/μL)")+
  scale_x_discrete(limits=c("clio", "knight", "control"))+
  scale_color_discrete(name="Species", labels=c("Atlantic salmon calicivirus", "Cutthroat trout virus-2",
                                                "Erythrocytic necrosis virus","P. salmonis", "Atlantic salmon calcivirus","Atlantic salmon", 
                                                "Candidatus S. salmonis", "T. finnmarkense", "T. maritimum"))+
  theme_bw()
scatter

####Hurdle model####
#subset data for only Atlantic salmon
sasa <- subset(data, data$Target.Name=="sasa" & data$site!="control")
sasa$logquant <- log(sasa$Quantity+1)
sasa$non_zero <- ifelse(sasa$Quantity>0, 1, 0)
sasa.plot <- ggplot(data=sasa, aes(x=site, y=Quantity, colour=site))+
  geom_boxplot()

sasa.clio <- subset(sasa, site == "clio")
sasa.knight <- subset(sasa, site == "knight")

#Model non-zero data using a binomial:
m.bin.clio <- glm(non_zero ~ 1, data = sasa.clio, family=binomial(link="logit"))
m.bin.knight <- glm(non_zero ~ 1, data = sasa.knight, family=binomial(link="logit"))

#Model zero data using a gamma:
m.gam.clio <- glm(Quantity ~ 1, data = subset(sasa.clio, non_zero == 1), family=Gamma(link="log"))
m.gam.knight <- glm(Quantity ~ 1, data = subset(sasa.knight, non_zero == 1), family=Gamma(link="log"))

#Extract model coefficients
m.bin.clio.coef <- plogis(coef(m.bin.clio)) #binomial coefficients are on logit scale so we inverse them
m.bin.knight.coef <- plogis(coef(m.bin.knight)) #binomial coefficients are on logit scale so we inverse them
m.gam.clio.coef <- exp(coef(m.gam.clio)) #gamma coefficients are on a log scale so we exponentiate them
m.gam.knight.coef <- exp(coef(m.gam.knight)) #gamma coefficients are on a log scale so we exponentiate them

#We combine binomial and gamma models by adding means on log scale and re-exponentiating them
pred.clio <- exp(log(m.bin.clio.coef) + log(m.gam.clio.coef))
pred.knight <- exp(log(m.bin.knight.coef) + log(m.gam.knight.coef))

sasa.plot + geom_hline(yintercept=pred.clio) + geom_hline(yintercept=pred.knight)

####Calculate confidence intervals on our model prediction####
# We are going to use a parametric bootstrapping approach using the boot function in the boot package
#First we create a function specifying how to generate random variables to hand the boot function
hurdle.fun <- function(data, i) {
  dat.boot <- data[i, ]
  m1 <- glm(non_zero ~ 1, data = dat.boot, family=binomial(link="logit"))
  m2 <- glm(Quantity ~ 1, data = subset(dat.boot, non_zero == 1), family=Gamma(link="log"))
  bin_coef <- plogis(coef(m1))
  gamma_coef <- exp(coef(m2))
  exp(log(bin_coef) + log(gamma_coef))
}

#Generate bootstraps and bootstrapped confidence intervals
set.seed(1)
b.clio <- boot(sasa.clio, hurdle.fun, R = 1000)
b.knight <- boot(sasa.knight, hurdle.fun, R = 1000)
b.ci.clio <- boot.ci(b.clio, type = "bca")
b.ci.knight <- boot.ci(b.knight, type = "bca")

preds <- data.frame(site = c("clio", "knight"),
                      meanpred = c(pred.clio[[1]], pred.knight[[1]]),
                      lower = c(b.ci.clio$bca[[4]], b.ci.knight$bca[[4]]),
                      upper = c(b.ci.clio$bca[[5]], b.ci.knight$bca[[5]]))

ggplot() + 
  geom_jitter(data=sasa, aes(x=site, y=Quantity, colour=site))+ 
  geom_point(data=preds, aes(x=site, y=meanpred))+
  geom_errorbar(data=preds, aes(x=site, ymin=lower, ymax=upper))


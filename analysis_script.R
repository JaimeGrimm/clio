####Analysis of Clio channel and knight inlet data collected summer 2024#
#Last updated October 8, 2024
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
                        TRUE~'control'))%>%
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
                                                "Erythrocytic necrosis virus","P. salmonis", "Atlantic salmon calcivirus","Atlantic salmon", 
                                               "Candidatus S. salmonis", "T. finnmarkense", "T. maritimum"))
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
sasa$non_zero <- ifelse(sasa$Quantity>0, 1, 0)
sasa.plot <- ggplot(data=sasa, aes(x=site, y=Quantity))+
  geom_boxplot()

####Fit models to each treatment separately####

#Model non-zero data using a binomial:
m1 <- glm(non_zero ~ site, data = sasa, family=binomial(link="logit"))

#Model zero data using a gamma:
m2 <- glm(Quantity ~ site, data = subset(sasa, non_zero == 1), family=Gamma(link="log"))

#Extract model coefficients
m1.coef <- plogis(coef(m1)) #binomial coefficients are on logit scale so we inverse them
m2.coef <- exp(coef(m2)) #gamma coefficients are on a log scale so we exponentiate them
#Extract model predictions
pred1 <- predict(m1, se=TRUE, type="link")
pred2 <- predict(m2, se=TRUE, type="link")

#We combine binomial and gamma models by adding means on log scale and re-exponentiating them
pred <- exp(log(m1.coef) + log(m2.coef))

sasa.plot + geom_hline(yintercept=pred)

####Calculate confidence intervals on our model prediction####
# We are going to use a parametric bootstrapping approach using the boot function in the boot package

#First we create a function specifying how to generate random variables to hand the boot function
hurdle.fun <- function(data, i) {
  dat.boot <- data[i, ]
  m1 <- glm(non_zero ~ site, data = dat.boot, family=binomial(link="logit"))
  m2 <- glm(Quantity ~ site, data = subset(dat.boot, non_zero == 1), family=Gamma(link="log"))
  bin_coef <- plogis(coef(m1)[[1]])
  gamma_coef <- exp(coef(m2)[[1]])
  exp(log(bin_coef) + log(gamma_coef))
}

#Generate bootstraps and bootstrapped confidence intervals
b <- boot(sasa, hurdle.fun, R = 500)
b.ci <- boot.ci(b, type = "bca")

sasa.plot +
  geom_hline(yintercept = pred) +
  geom_hline(yintercept = b.ci$bca[c(4:5)],
             colour = "darkgrey")

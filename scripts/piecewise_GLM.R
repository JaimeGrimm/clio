####Analysis of Clio channel and knight inlet data collected summer 2024#
#Last updated February 3, 2025
library("tidyverse")
library("MASS")
library("vcdExtra")  
library("boot")

#Load clean data
source("~/Documents/GitHub/clio/scripts/data_preparation.R")

####Visualization####
boxplot <- data %>% filter(Site != "Control") %>% 
  ggplot(aes(x=Site, y=log(Conc+1), fill=Target))+
  geom_boxplot() +
  theme_bw()+
  theme(axis.title.x = element_text(size=18),
  axis.title.y = element_text(size=18),
  legend.text = element_text(size=18),
  axis.text = element_text(size=16))+
  ylab("log(eDNA concentration + 1) (copies/μL)") +
  scale_fill_discrete(name="Species", labels = c(env = "Erythrocytic necrosis virus", pa_ther = "Paranucleospora theridion",
                                                 pisck_sal = "Piscirickettsia salmonis", sasa = "Atlantic salmon", 
                                      sch = "Candidatus Syngnamydia salmonis", Te_fin = "Tenacibaculum finnmarkense",
                                      Te_mar = "Tenacibaculum maritimum"))

boxplot

scatter <- ggplot(data=data, aes(x=Site, y=Conc, color=Target))+
  geom_point()+
  geom_jitter()+
  theme(axis.title.x = element_text(size=18),
        axis.title.y = element_text(size=18),
        legend.text = element_text(size=18),
        axis.text = element_text(size=20))+
  ylab("eDNA concentration (copies/μL)")+
  scale_color_discrete(name="Species", labels=c(env = "Erythrocytic necrosis virus", pa_ther = "Paranucleospora theridion",
                                                pisck_sal = "Piscirickettsia salmonis", sasa = "Atlantic salmon", 
                                                  sch = "Candidatus Syngnamydia salmonis", Te_fin = "Tenacibaculum finnmarkense",
                                                Te_mar = "Tenacibaculum maritimum"))+
  theme_bw()
scatter

##Model fits ----
envdata$logquant <- log(envdata$Conc+1)
envdata$non_zero <- ifelse(envdata$Conc > 0, 1, 0)

targets <- unique(envdata$Target)
clio.datasets <- list()
knight.datasets <- list()

for (i in 1:length(targets)){  
  clio.datasets[[i]] <- envdata %>% filter(Target == targets[i] & Site == "Clio")
  knight.datasets[[i]] <- envdata %>% filter(Target == targets[i] & Site == "Knight")
}

names(clio.datasets) <- targets
names(knight.datasets) <- targets

##Hurdle model ----
#Model non-zero data using a binomial and zero data using a gamma:
bin.models.clio <- list()
bin.models.knight <- list()
gam.models.clio <- list()
gam.models.knight <- list()
coef.bin.clio <- c()
coef.bin.knight <- c()
coef.gam.clio <- c()
coef.gam.knight <- c()

#Run the model with each dataset (one for each target)
#There isn't enough non-zero data for Te fin. The tryCatch function gives the error but continues the loop for the remaining iterations
for (j in 1:length(targets)){
  tryCatch({
  bin.models.clio[[j]] <- glm(non_zero ~ 1 + temp + salinity + turbidity, data = clio.datasets[[j]], family=binomial(link="logit"))
  bin.models.knight[[j]] <- glm(non_zero ~ 1 + temp + salinity + turbidity, data = knight.datasets[[j]], family=binomial(link="logit"))
  gam.models.clio[[j]] <- glm(Conc ~ 1 + temp + salinity + turbidity, data = subset(clio.datasets[[j]], non_zero == 1), family=Gamma(link="log"))
  gam.models.knight[[j]] <- glm(Conc ~ 1 + temp + salinity + turbidity, data = subset(knight.datasets[[j]], non_zero == 1), family=Gamma(link="log"))
  coef.bin.clio[j] <- plogis(coef(bin.models.clio[[j]]))
  coef.bin.knight[j] <- plogis(coef(bin.models.knight[[j]]))
  coef.gam.clio[j] <- exp(coef(gam.models.clio[[j]]))
  coef.gam.knight[j] <- exp(coef(gam.models.knight[[j]]))
  #extract coefficienct and backtransform according to link functions
  }, error = function(e){cat("ERROR :",conditionMessage(e), "\n")})
}

names(bin.models.clio) <- targets
names(bin.models.knight) <- targets
names(gam.models.clio) <- targets
names(gam.models.knight) <- targets

#Put coefficients into a dataframe
coefs <- data.frame(Target = targets, 
                    knight.binom = coef.bin.knight,
                    clio.binom = coef.bin.clio,
                    knight.gam = coef.gam.knight,
                    clio.gam = coef.gam.clio)

#We combine binomial and gamma models by adding means on log scale and re-exponentiating them
coefs$mean.clio <- exp(log(coefs$clio.binom) + log(coefs$clio.gam))
coefs$mean.knight <- exp(log(coefs$knight.binom) + log(coefs$knight.gam))

#For models that include environmental data this is trickier - need to predict with set values for env data
newdata <- expand.grid(Site = c("Clio", "Knight"),
                         Target = sort(unique(envdata$Target)),
                         temp = mean(envdata$temp),
                         salinity = mean(envdata$salinity),
                         turbidity = mean(envdata$turbidity, na.rm = TRUE),
                         non_zero = c(0, 1))


#Use models and new data to predict mean at average environmental variables
mean.clio <- c()
mean.knight <- c()
for (z in 1:length(targets)){
  tryCatch({
  mean.clio[z] <- predict(bin.models.clio[[z]], newdata, type = "response") * predict(gam.models.clio[[z]], newdata, type = "response")
  mean.knight[z] <-  predict(bin.models.knight[[z]], newdata, type = "response") * predict(gam.models.knight[[z]], newdata, type = "response")
  }, error = function(e){cat("ERROR :",conditionMessage(e), "\n")})
}

coefs$mean.clio <- mean.clio
coefs$mean.knight <- mean.knight

#Calculate confidence intervals on our model prediction----
# We are going to use a parametric bootstrapping approach using the boot function in the boot package
#First we create a function specifying how to generate random variables to hand the boot function
hurdle.fun <- function(data, i) {
  dat.boot <- envdata[i, ]
  m1 <- glm(non_zero ~ 1 + temp + salinity + turbidity, data = dat.boot, family=binomial(link="logit"))
  m2 <- glm(Conc ~ 1 + temp + salinity + turbidity, data = subset(dat.boot, non_zero == 1), family=Gamma(link="log"))
  bin_coef <- plogis(coef(m1))
  gamma_coef <- exp(coef(m2))
  exp(log(bin_coef) + log(gamma_coef))
}

#Generate bootstraps and bootstrapped confidence intervals
b.clio <- list()
b.knight <- list()
b.ci.knight <- list()
b.ci.clio <- list()

for (k in 1:length(targets)){
  tryCatch({
  set.seed(1)
  b.clio[[k]] <- boot(clio.datasets[[k]], hurdle.fun, R = 1000)
  b.knight[[k]] <- boot(knight.datasets[[k]], hurdle.fun, R = 1000)
  b.ci.clio[[k]] <- boot.ci(b.clio[[k]], type = "bca")
  b.ci.knight[[k]] <- boot.ci(b.knight[[k]], type = "bca")
  }, error = function(e){cat("ERROR :",conditionMessage(e), "\n")})
}

lower.knight <- c()
upper.knight <- c()
lower.clio <- c()
upper.clio <- c()
for (m in 1:length(targets)){
  tryCatch({
  lower.knight[m] <- b.ci.knight[[m]]$bca[[4]]
  upper.knight[m] <- b.ci.knight[[m]]$bca[[5]]
  lower.clio[m] <- b.ci.clio[[m]]$bca[[4]]
  upper.clio[m] <- b.ci.clio[[m]]$bca[[5]]
  }, error = function(e){cat("ERROR :",conditionMessage(e), "\n")})
}


##Following Sean bootstrapping----
#Simulate new data from models
sims.gamma <- simulate(gam.models.clio[[1]], 1000)
sims.binom <- simulate(bin.models.clio[[1]], 1000)

#Parametric bootstrap
refit <- function(x1, x2){
  data.binom.clio <- bin.models.clio[[1]]$data 
 # data.binom.knight <- bin.models.knight[[1]]$data
  clio.datasets[[1]]$non_zero <- x1
  newmod.binom.clio <- update(bin.models.clio[[1]], data = data.binom.clio)
  #newmod.binom.knight <- update(bin.models.knight[[1]], data = data.binom.knight)
  
  
  data.gamma.clio <- gam.models.clio[[1]]$data
  #data.gamma.knight <- gam.models.knight[[1]]$data
  clio.datasets[[1]]$Conc <- x2
  newmod.gamma.clio <- update(gam.models.clio[[1]], data = data.gamma.clio)
  #newmod.gamma.knight <- update(gam.models.knight[[1]], data = data.gamma.knight)
  
  predict(newmod.binom.clio, newdata, type = response) * predict(newmod.gamma.clio, newdata, type = response)
  #predict(newmod.binom.knight, newdata, type = response) * predict(newmod.gamma.knight, newdata, type = response)
}

boot <- mapply(refit, sims.binom, sims.gamma)


#Plots----

coefs <- coefs %>% 
  mutate(lower.knight = lower.knight, upper.knight = upper.knight,
         lower.clio = lower.clio, upper.clio = upper.clio) %>% 
  pivot_longer(cols = !Target, names_to = "quantity", values_to = "value") %>% 
  mutate(Site = case_when(str_detect(quantity, "clio")~"Clio",
    str_detect(quantity, "knight")~"Knight")) %>% 
  pivot_wider(names_from = quantity, values_from = value) %>% 
  mutate(binomial.coefficient = ifelse(is.na(knight.binom), clio.binom, knight.binom),
         gamma.coefficient = ifelse(is.na(knight.gam), clio.gam, knight.gam),
         mean = ifelse(is.na(mean.knight), mean.clio, mean.knight),
         upper = ifelse(is.na(upper.knight), upper.clio, upper.knight),
         lower = ifelse(is.na(lower.knight), lower.clio, lower.knight)) %>% 
  dplyr::select(Target, Site, binomial.coefficient, gamma.coefficient, mean, upper, lower)

#Main plot----
target.names <- c(env = "Erythrocytic necrosis virus", pa_ther = "Paranucleospora theridion",
                  pisck_sal = "Piscirickettsia salmonis", sasa = "Atlantic salmon", 
                  sch = "Candidatus Syngnamydia \n salmonis", Te_fin = "Tenacibaculum finnmarkense",
                  Te_mar = "Tenacibaculum maritimum")

plot.data <- data %>% filter(Site != "Control")
plot.coefs <- coefs %>% filter(Target !="Te_fin" & Target != "Te_mar")
ggplot() + 
  geom_jitter(data = plot.data, aes(x=Site, y=Conc, colour=Site))+ 
  geom_point(data = plot.coefs, aes(x=Site, y=mean))+
  geom_errorbar(data = plot.coefs, aes(x=Site, ymin=lower, ymax=upper)) +
  facet_wrap(~Target, scales = "free", labeller = as_labeller(target.names)) +
  ylab("eDNA Concentration (copies/μL)") +
  xlab(NULL)+
  theme(legend.position = "none")




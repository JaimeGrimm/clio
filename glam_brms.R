####Analysis of Clio channel and knight inlet data collected summer 2024#
#Bayesian regression models using stan. The goal is to model species distributions with species ID as a random effect
#Created November 18, 2024
#Last updated November 18, 2024

library(tidyr)
library(brms)
library(broom)
library(broom.mixed) 
library(emmeans)
library(tseries)

####Load and clean data----
source("~/Documents/GitHub/clio/scripts/data_preparation.R")

#Count data for Poisson models
data.ct <- data %>% mutate(across(Quantity, round)) %>% #round quantity into count data
                    dplyr::select(Sample.Name, Target.Name, Quantity, site, seq)

#Continuous data for gamma models
gamdat <- data %>%  dplyr::select(Sample.Name, Target.Name, Quantity, site, seq) %>% 
  group_by(Target.Name, Sample.Name, site, seq) %>% 
  dplyr::summarise(mean = mean(Quantity)) %>% 
  ungroup()

#Explore Atlantic salmon data
sasa <- data.ct %>% filter(Target.Name == "sasa")

sasafit <- brm(Quantity ~ site, data.ct = sasa, family = zero_inflated_poisson("log"), chain = 2,
           cores = 2)

####Plot data ----
#Histograms by species
ggplot(data = data, aes(x=Quantity))+
  geom_histogram()+
  facet_wrap(~ Target.Name)
  
#Boxplots by species
ggplot(data = data, aes(x=site, y=Quantity))+
  geom_boxplot()+
  facet_wrap(~ Target.Name)

####Set sampling regime and priors----
#Priors


#### Model 1 ----
#Count data
#Random effect by species at intercept
fit <- brm(Quantity ~ site + (1|Target.Name), data = data.ct, family = zero_inflated_poisson("log"), chain = 2,
           cores = 2)

#### Model 2 ----
#Count data
#Random effect by species at slope and intercept
m2 <- brm(Quantity ~ site + (site|Target.Name), data = data.ct, family = zero_inflated_poisson("log"), chain = 2,
          cores = 2)

#Put output into a tidy dataframe:
m2df <- tidy(m2)
#Posterior predictive check:
pp_check(m2) #Not great

#### Model 3 ----
#Gamma hurdle without defining the binary hurdle part (to build intuition)
m3 <- brm(bf(Quantity ~ site + (site|Target.Name), hu ~ 1), data = data, family = hurdle_gamma(),
          chain = 2, cores = 2)
m3df <- tidy(m3) 

#What proportion of our data are zeros?
hu <- m3df %>% filter(term == "hu_(Intercept)") %>% pull(estimate)
plogis(hu) #39%

conditional_effects(m3)
pp_check(m3)

#### Autocorrelation----
#Explore autocorrelation and identify correct ARIMA model
#Plot mean versus sample order
ggplot(data=gamdat)+
  geom_point(aes(x=seq, y=mean, col=Target.Name, group=Target.Name))+
  geom_line(aes(x=seq, y=mean, col=Target.Name, group=Target.Name))+
  facet_wrap(~ Target.Name)

# Remove trend by calculating difference between consecutive points for each species
gamdat <- gamdat %>% arrange(by_group=Target.Name, seq) %>% 
  group_by(Target.Name) %>% 
  mutate(diff = mean-lag(mean)) %>% 
  mutate(diff2 = diff-lag(diff)) %>% 
  ungroup() 
gamdat$seq <- as.numeric(gamdat$seq)

uniqueTarget <- unique(gamdat$Target.Name)

#Check if time series is stationary
for (i in 1:length(uniqueTarget)){  
  temp <- gamdat %>% filter(Target.Name==uniqueTarget[i])
  print(adf.test(gamdat$mean))
}

#What about with species pooled?
adf.test(gamdat$mean)
#not stationary

overalldiff <- c(diff(gamdat$mean, differences=1))
adf.test(overalldiff)
#stationary with first order differences

#All are non-stationary so we have to take differences and test again
#Non-stationary data indicates there is some time-dependent structure that does not have constant variance over time

#Plot differences versus sample order
ggplot(data=gamdat)+
  geom_point(aes(x=seq, y=diff2, col=Target.Name, group=Target.Name))+
  geom_line(aes(x=seq, y=diff2, col=Target.Name, group=Target.Name))+
  facet_wrap(~ Target.Name)

#Check if differences are stationary. 
for (i in 1:length(uniqueTarget)){
  temp <- gamdat %>% filter(Target.Name==uniqueTarget[i]) %>% 
    drop_na()
  print(adf.test(temp$diff2))
}

#First order differences were not stationary so took the difference again

par(mfrow=c(2,4))
for (i in 1:length(uniqueTarget)){
  m <- gamdat %>% filter(Target.Name==uniqueTarget[i])
  auto <- pacf(m$mean, na.action = na.pass, plot= FALSE)
  plot(auto, main=uniqueTarget[i])
}

#acf and pacf with species pooled:
acf(overalldiff, plot=TRUE)
pacf(overalldiff, plot=TRUE)

#### Model 4 ----
#Gamma hurdle with site being a predictor of zero/not zeros process and autocorrelation
gamdat$seq <- as.numeric(gamdat$seq)

m4 <- brm(bf(mean ~ site + (1 + site|Target.Name), hu ~ site, autocor = ~ar(time = seq, 
            gr = Target.Name:Sample.Name, p = 1, cov=TRUE)), 
          data = gamdat, family = hurdle_gamma(), chain = 2, cores = 2)
m4df <- tidy(m4) 

#What proportion of our data are zeros?
hu_knight <- m4df %>% filter(term == "hu_(Intercept)") %>% pull(estimate)
hu_clio <- m4df %>% filter(term == "hu_siteknight") %>% pull(estimate)
pp_check(m4)
plot(conditional_effects(m4), points=TRUE) 

####Model diagnostics ----
mcmc_plot(m4, type="trace")
mcmc_plot(m4, type="pairs")
mcmc_plot(m4, type="dens")

ggplot()+
  geom_jitter(data=gamdat, aes(x=site, y=mean, colour=site)) +
  geom_point(aes(x=gamdat$site, y=m4df$estimate))#+
  geom_errorbar(data=preds, aes(x=site, ymin=lower, ymax=upper))

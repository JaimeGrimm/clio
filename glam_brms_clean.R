####Analysis of Clio channel and knight inlet data collected summer 2024#
#Bayesian regression models using stan. The goal is to model species distributions with species ID as a random effect
#Created November 18, 2024
#Last updated November 18, 2024

library(tidyr)
library(tidyverse)
library(brms)
library(broom)
library(broom.mixed) 
library(emmeans)
library(tseries)
library(jtools)
library(gghalves)
library(ggthemes)
library(tidybayes)
library(egg)


####Load and clean data----
source("~/Documents/GitHub/clio/scripts/data_preparation.R")

#Mean data for gamma models
data <- data %>%  dplyr::select(Sample.Name, Target.Name, Quantity, site, seq) %>% 
  group_by(Target.Name, Sample.Name, site, seq) %>% 
  dplyr::summarise(mean = mean(Quantity)) %>% 
  ungroup()

#Add column to dataframe specifying previous sample
gamdat <- data %>% group_by(Target.Name) %>% 
  arrange(by_group = seq) %>% 
  mutate(prevmean = lag(mean, n=1)) %>%
  mutate(prevmean = ifelse((seq == "01" | seq == "8.5"), 0, prevmean)) %>% 
  ungroup() %>% 
  mutate(prevmean = replace_na(prevmean, 0)) %>%
  mutate(clio = ifelse(site == "clio", 1, 0)) %>%
  mutate(knight = ifelse(site == "knight", 1, 0))
gamdat$seq <- as.numeric(gamdat$seq)

####Plot data ----
mycolours <- c("#855C75","#D9AF6B","#AF6458","#736F4C","#526A83","#625377","#68855C","#9C9C5E",
               "#A06177","#8C785D","#467378","#7C7C7C")

#Histograms by species
ggplot(data = data, aes(x=Quantity))+
  geom_histogram()+
  facet_wrap(~ Target.Name)

#Boxplots by species
ggplot(data = data, aes(x=site, y=Quantity))+
  geom_half_point(aes(color=Target.Name), side="l")+
  geom_half_boxplot(aes(fill=Target.Name), side="r")+
  facet_wrap(~ Target.Name, scales = "free")

#### Full model ----
chains <- 4
iter <- 4000
warmup <- 1000
seed <- 1234

#**Fit model----
#Gamma hurdle with site being a predictor of zero/not zeros process and autocorrelation
m4 <- brm(
  bf(mean ~ clio + knight + prevmean + (1 + site|Target.Name),
     hu ~ clio + knight + prevmean + (1 + site|Target.Name)),
  data = gamdat,
  family = hurdle_gamma(),
  chains = chains,
  iter = iter,
  warmup = warmup,
  seed = seed,
  core  = 4,
  control = list(max_treedepth = 13),
  prior = set_prior("normal(0.50, 10)", class = "b", lb = 0))

m4df <- tidy(m4) 

#**Extract draws----
#Random (species) effects
drawsr <- spread_draws(m4,
                       r_Target.Name[condition, term],
                       r_Target.Name__hu[condition, term])
drawsr <- drawsr %>% mutate(mean.adj = exp(r_Target.Name)) %>% 
  mutate(hu.adj = 1 - plogis(r_Target.Name__hu)) %>% 
  mutate(mh = mean.adj * hu.adj)  %>% 
  select(c(condition, term, mh, .chain, .iteration, .draw)) %>% 
  pivot_wider(names_from = term, values_from = mh)

#Atlantic salmon posteriors
sasa <- drawsr %>% filter(condition == "sasa")
sasa$adj_clio <- sasa$Intercept
sasa$adj_control <- sasa$Intercept + sasa$sitecontrol
sasa$adj_knight <- sasa$Intercept + sasa$siteknight

sasa.longr <- sasa %>% pivot_longer(cols = c(adj_clio, adj_control, adj_knight))

#Atlantic salmon posterior plot
sasa.longr %>% 
  ggplot(aes(x=value, y=0, fill=name, alpha = 0.5)) +
  stat_halfeye(aes(color = name)) +
  theme_bw()+
  labs(x = "Relative DNA concentration", 
       title = "Atlantic salmon", y = NULL)+
  guides(alpha = "none")+
  xlim(0, 5)

#Calculate the probability that clio > knight
sasa.diff <- sasa$adj_clio - sasa$adj_knight
sum(sasa.diff < 0)


diffsasa <- compare_levels(drawsr, variable = r_Target.Name, by = term) %>% 
  filter(condition == "sasa" & term == "siteknight - Intercept")

hypothesis(m4, "Target.Name[sasa,Intercept] > Target.Name[sasa,siteknight]", class = "r")
hypothesis(sasa.longr, "Intercept > siteknight")

#Fixed effects   
drawsf <- spread_draws(m4, b_Intercept, b_hu_Intercept,
                       b_clio, b_knight, b_prevmean, b_hu_clio,
                       b_hu_knight, b_hu_prevmean)

#Full model posteriors
drawsf$B_control <- exp(drawsf$b_Intercept) * (1 - plogis(drawsf$b_hu_Intercept))
drawsf$B_clio <- exp(drawsf$b_Intercept + drawsf$b_clio) * (1 - plogis(drawsf$b_hu_Intercept + drawsf$b_hu_clio))
drawsf$B_knight <- exp(drawsf$b_Intercept + drawsf$b_knight) * (1 - plogis(drawsf$b_hu_Intercept + drawsf$b_hu_knight))

#Hurdle posteriors
drawsf$HU_control <- (plogis(drawsf$b_hu_Intercept))
drawsf$HU_clio <- (plogis(drawsf$b_hu_Intercept + drawsf$b_hu_clio))
drawsf$HU_knight <- (plogis(drawsf$b_hu_Intercept + drawsf$b_hu_knight))

draws.longf <- drawsf %>%  pivot_longer(everything())

#**Posterior plots----
#Full model
draws.longf %>% 
  filter(name == "B_clio"| name == "B_knight" | name == "B_control") %>% 
  ggplot(aes(x=value, y=0, fill=name, alpha = 0.5)) +
  stat_halfeye(aes(color = name)) +
  theme_bw()+
  # scale_fill_manual(values = mycolours)+
  xlim(-5,200) +
  labs(x = "Relative DNA concentration", 
       title = "mu and hu parts of the model", y = NULL)+
  guides(alpha = "none")

#Difference between Clio and Knight
draws.longf %>% 
  filter(name == "diff_clio"| name == "diff_knight") %>% 
  ggplot(aes(x=value, y=0, fill=name, alpha = 0.5)) +
  stat_halfeye(aes(color = name)) +
  theme_bw()+
  # scale_fill_manual(values = mycolours)+
  xlim(-5,200) +
  labs(x = "Relative DNA concentration", 
       title = "mu and hu parts of the model", y = NULL)+
  guides(alpha = "none")


#Hurdle component only
draws.longf %>% 
  filter(name == "HU_clio"| name == "HU_knight" | name == "HU_control") %>% 
  ggplot(aes(x=value, y=0, fill=name, alpha = 0.5)) +
  stat_halfeye(aes(color = name)) +
  theme_bw()+
  labs(x = "Proportion of zeros in data", 
       title = "hu part of the model", y = NULL)+
  guides(alpha = "none")

#What proportion of our data are zeros?
hu_clio <- m4df %>% filter(term == "hu_(Intercept)") %>% pull(estimate)
hu_knight <- m4df %>% filter(term == "hu_siteknight") %>% pull(estimate)
hu_control <- m4df %>% filter(term=="hu_sitecontrol") %>% pull(estimate)
pp_check(m4)

#**Plots of conditional effects----
#This is for both model components (non-zero and zeros)
cond.data <- plot(conditional_effects(m4, re_formula = NULL), plot=FALSE)[[1]]$data 

ggplot()+
  geom_jitter(data=gamdat, aes(x=site, y=mean, color=Target.Name), width=0.2, alpha=0.4)+
  geom_point(data = cond.data, aes(x=effect1__, y=estimate__), size=2, color="black")+
  geom_errorbar(data= cond.data, aes(x= effect1__, ymin=lower__, ymax=upper__),color="black", width=0.5)+
  labs(x="Site", y="Predicted mean DNA concentration", subtitle = "Combined parts of the model (\"mu\" and \"hu\")")+
  theme_bw()

#This is for just the hurdle components
hu.data <- plot(conditional_effects(m4, dpar="hu", re_formula = NULL), plot=FALSE)[[1]]$data
ggplot()+
  #  geom_jitter(data=gamdat, aes(x=site, y=mean, color=Target.Name), width=0.2, alpha=0.4)+
  geom_point(data = hu.data, aes(x=effect1__, y=estimate__), size=2, color="black")+
  geom_errorbar(data= hu.data, aes(x= effect1__, ymin=lower__, ymax=upper__),color="black", width=0.5)+
  labs(x="Site", y="Predicted mean DNA concentration", subtitle = "Hurdle part of the model (\"hu\") only")+
  theme_bw()

#**Model diagnostics ----
mcmc_plot(m4, type="trace")
mcmc_plot(m4, type="pairs", variable=variables(m4)[1:8])
mcmc_plot(m4, type="pairs", variable=variables(m4)[9:20])
mcmc_plot(m4, type="dens")

ggplot()+
  geom_jitter(data=gamdat, aes(x=site, y=mean, colour=site)) +
  geom_point(aes(x=gamdat$site, y=m4df$estimate))#+
geom_errorbar(data=preds, aes(x=site, ymin=lower, ymax=upper))


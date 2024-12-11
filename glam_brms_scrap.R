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

#### Some exploratory models ----
# **Model 1 
#Count data
#Random effect by species at intercept
fit <- brm(Quantity ~ site + (1|Target.Name), data = data.ct, family = zero_inflated_poisson("log"), chain = 2,
           cores = 2)

# ** Model 2 
#Count data
#Random effect by species at slope and intercept
m2 <- brm(Quantity ~ site + (site|Target.Name), data = data.ct, family = zero_inflated_poisson("log"), chain = 2,
          cores = 2)

#Put output into a tidy dataframe:
m2df <- tidy(m2)
#Posterior predictive check:
pp_check(m2) #Not great

#** Model 3 
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
  facet_wrap(~ Target.Name, scales = "free")

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

#### Full model ----
chains <- 4
iter <- 4000
warmup <- 1000
seed <- 1234

#Add column to dataframe specifying previous sample
gamdat <- gamdat %>% group_by(Target.Name) %>% 
  arrange(by_group = seq) %>% 
  mutate(prevmean = lag(mean, n=1)) %>%
  mutate(prevmean = ifelse((seq == "01" | seq == "8.5"), 0, prevmean)) %>% 
  mutate(prevmean = replace_na(prevmean, 0)) %>%
  mutate(clio = ifelse(site == "clio", 1, 0)) %>%
  mutate(knight = ifelse(site == "knight", 1, 0))
  
#Gamma hurdle with site being a predictor of zero/not zeros process and autocorrelation
gamdat$seq <- as.numeric(gamdat$seq)
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

#Put draws into a dataframe
#Random (species) effects
drawsr <- spread_draws(m4,
                      r_Target.Name[condition, term],
                      r_Target.Name__hu[condition, term])
drawsr <- drawsr %>% mutate(mean.adj = exp(r_Target.Name)) %>% 
                     mutate(hu.adj = 1 - plogis(r_Target.Name__hu)) %>% 
                     mutate(mh = mean.adj * hu.adj) # %>% 
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
  
#**Pulling out posteriors from random effects----

  
  
#**Posterior predictions----
#Posterior predictions includes all sources of variance (residuals, fixed effects, grouping factors)
#Posterior estimates includes uncertainty in model parameters but not observation-level variance
new.sasa <- expand.grid(site = c("clio", "knight"),
                        Target.Name = c("sasa"),
                        Sample.Name = c("KN0925"))
new.sasa$seq = c(11.0, 10.0)

#Expected values of the posterior predictive distribution for sasa:
pred.sasa <- m4 %>% predicted_draws(newdata = new.sasa)

epred.sasa <- m4 %>% epred_draws(newdata = new.sasa)

#Plot predicted draws and expectation of predicted draws:
plot.sasa <- bind_rows(
  "Predicted draws" = pred.sasa,
  "Expectation of predicted draws" = rename(epred.sasa, .prediction = .epred),
  .id = "draw_type") %>% 
  mutate(draw_type = fct_inorder(draw_type))

plot.sasa %>% group_by(draw_type, site) %>% 
  median_hdi(.prediction)

ggplot(plot.sasa, aes(x = .prediction, fill = site, alpha = 0.3)) +
  stat_halfeye()+
  labs(x = "Predicted DNA density", y = "Density", fill = "Site") +
  facet_wrap(vars(draw_type), scales="free_x", ncol=1) +
  guides(size="legend", alpha= "none")+
  theme_bw() +
  theme(panel.grid.major = element_blank(),
           panel.grid.minor = element_blank())

#Calculate differences between expected values of predicted draws:
sasa_effect_draws <- m4 %>% 
  emmeans(~site, at=list(Target.Name = "sasa"), epred=TRUE) %>% 
  contrast() %>% 
  gather_emmeans_draws()

#Plot differences
ggplot(sasa_effect_draws, aes(x = .value)) +
  stat_halfeye() +
  labs(x ="Average marginal effect of site on Atlantic salmon DNA denisty", y= "Density")+
  theme_bw()

#Calculate global average marginal effect of site
newdata <- expand.grid(site = c("knight", "clio", "control"),
                      Target.Name = c("env", "pa_ther", "pisck_sal", "sasa", "sch",
                                      "Te_fin", "Te_mar"),
                      Sample.Name = c("KN0925"))
newdata$seq = c(1:21)

grand_mean <- m4 %>% 
  epred_draws(newdata = newdata,
              re_formula = NULL)

plot_grand_mean <- ggplot(grand_mean, aes(x=.epred, y= Target.Name, fill = site)) +
  stat_halfeye() +
  labs(x = "Predicted DNA density", y = NULL, fill = "Species", subtitle = "Posterior predictions")+
  theme_bw()+
  theme(legend.position = "bottom")+
  xlim(0, 1000)

#No-autocor model ----
m5 <- brm(
  bf(mean ~ site + (1 + site|Target.Name),
     hu ~ site),
  data = gamdat,
  family = hurdle_gamma(),
  chains = chains,
  iter = iter,
  warmup = warmup,
  seed = seed)
m5df <- tidy(m5) 

m5df

#Plot conditional effects, comparing model 4 and 5:
cond.m5 <- plot(conditional_effects(m5, re_formula = NULL), plot=FALSE)[[1]]$data 

ggplot()+
  geom_jitter(data=gamdat, aes(x=site, y=mean, color=Target.Name), width=0.2, alpha=0.4)+
  geom_point(data = cond.m5, aes(x=effect1__, y=estimate__), color = "black", size=2,
             position = position_nudge(x=-0.1))+
  geom_errorbar(data= cond.m5, aes(x= effect1__, ymin=lower__, ymax=upper__, color="no autocorrelation"), width=0.5,
                position = position_nudge(x=-0.1))+
  geom_point(data = cond.data, aes(x=effect1__, y=estimate__), size=2, 
             position = position_nudge(x=0.1))+
  geom_errorbar(data= cond.data, aes(x= effect1__, ymin=lower__, ymax=upper__, color = "autocorrelation"), 
                position = position_nudge(x=0.1), width=0.5)+
  labs(x="Site", y="Predicted mean DNA concentration", subtitle = "Combined parts of the model (\"mu\" and \"hu\")")+
  theme_bw() +
  scale_color_manual(name = "Model", breaks = c("no autocorrelation", "autocorrelation"),
                     values = c("no autocorrelation" = "black", "autocorrelation" = "red"))

#*Model diagnostics ----
mcmc_plot(m5, type="trace")
mcmc_plot(m5, type="pairs", variable=variables(m5)[1:6])
mcmc_plot(m5, type = "pairs", variable =variables(m5)[7:12])
mcmc_plot(m5, type="dens")

#*Plots of conditional and marginal effects ----

#Conditional effects
ggplot()+
  geom_jitter(data=gamdat, aes(x=site, y=mean, color=Target.Name), width=0.2, alpha=0.4)+
  geom_point(data = cond.m5, aes(x=effect1__, y=estimate__), color = "black", size=2)+
  geom_errorbar(data= cond.m5, aes(x= effect1__, ymin=lower__, ymax=upper__), color = "black", width=0.5)+
  labs(x="Site", y="Predicted mean DNA concentration", subtitle = "Combined parts of the model (\"mu\" and \"hu\")")+
  theme_bw() 

#Posterior estimates includes uncertainty in model parameters but not observation-level variance
new.sasa <- expand.grid(site = c("clio", "knight"),
                        Target.Name = c("sasa"),
                        Sample.Name = c("KN0925"))

#Expected values of the posterior predictive distribution for sasa:
pred.sasa <- m5 %>% predicted_draws(newdata = new.sasa)

epred.sasa <- m5 %>% epred_draws(newdata = new.sasa)

#Grand mean (no species information)

m5.grandmean <- m5 %>% epred_draws(newdata = expand.grid(site = c("knight", "clio"),
                                                         Sample.Name = c("KN0925")),
                                   re_formula = NA)
  
grandmean.plot <- ggplot(m5.grandmean, aes(x = .epred, y= "Grand mean", fill= site, alpha = 0.3)) +
  stat_halfeye() +
  labs(x = "Predicted DNA density", y = NULL, fill = "Site", subtitle = "Posterior predictions") +
  theme_bw() +
  guides(alpha = "none")+
  theme(legend.position = "bottom")
  
#Marginal effect (clio - knight)
m5.grand.ME <- m5 %>% 
  emmeans(~site, 
          at = list(Sample.Name = "KN0925"),
          epred = TRUE, re_formula = NA) %>% 
  contrast() %>% 
  gather_emmeans_draws()

grand.mean.me.plot <- ggplot(m5.grand.ME, aes(x= .value, y="Grand AME")) +
  stat_halfeye() +
  labs(x = "Average marginal effect of site", y=NULL,
       subtitle = "Marginal effect (clio-knight)")+
  theme_bw()+
  theme(legend.position = "bottom")

#Combine onto one plot
ggarrange(grandmean.plot, grand.mean.me.plot, nrow= 1)

#Conditional effects for species
all.species.dist <- m5 %>% 
  epred_draws(newdata = expand_grid(site = c("clio", "knight"),
                                    Sample.Name = "KN0925",
                                    Target.Name = c("env", "pa_ther", "pisck_sal", "sasa", "sch", "Te_mar")),
              re_formula=NULL)

plot.all.species <- ggplot(all.species.dist, aes(x = .epred, y=Target.Name, fill = site, alpha = 0.3)) +
  stat_halfeye(point_size = 1.5)+
  labs(x = "Predicted DNA density", y=NULL, fill = "Site", subtitle = "Posterior predictions") +
  theme_bw()+
  guides(alpha = "none")+
  theme(legend.position = "bottom")
plot.all.species

#All species average marginal effects
all.species.ame <- m5 %>% 
  emmeans(~ site + Target.Name,
          at = list(site = c("clio", "knight"),
                    Sample.Name = "KN0925",
                    Target.Name = c("env", "pa_ther", "pisck_sal", "sasa", "sch", "Te_mar")),
          epred=TRUE, re_formula=NULL) %>% 
  contrast(by = "Target.Name") %>% 
  gather_emmeans_draws()
         
plot.all.species.ame <- ggplot(all.species.ame, aes(x = .value, y=Target.Name)) +
  stat_halfeye(point_size = 1.5) +
  labs(x = "Average marginal effect of site", y = NULL, subtitle = "Marginal effect (clio-knight)") +
  theme_bw()+
  theme(legend.position = "bottom")
plot.all.species.ame

#Combine plots
ggarrange(plot.all.species, plot.all.species.ame, nrow=1)

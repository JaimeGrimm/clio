####Analysis of Clio channel and knight inlet data collected summer 2024#
#Bayesian regression models using stan. The goal is to model species distributions with species ID as a random effect
#Created November 18, 2024
#Last updated January 6, 2025

#Housekeeping ----
library(tidyr)
library(tidyverse)
library(brms)
library(emmeans)
library(gghalves)
library(ggthemes)
library(tidybayes)
library(scales)
library(ggtext)
library(repr)
library(marginaleffects)
options(repr.plot.width=12, repr.plot.height=8, repr.plot.dpi=300)

#Load data ----
source("~/Documents/GitHub/clio/scripts/data_preparation.R")

#Plot data ----
#Histograms by species
envdata_nc %>% ggplot(aes(x=Conc*100, fill = Target))+
  geom_histogram(bins = 10) +
  facet_grid(Site ~ Target, scales = "free") +
  theme_bw() + theme(legend.position = "none")

#Boxplots by species
envdata_nc %>% ggplot(aes(x=Target, y=Conc*100))+
  geom_half_point(aes(color=Target), side="l", alpha = 0.2,
                  transformation = PositionIdentity)+
  geom_half_boxplot(aes(fill=Target), side="r")+
  facet_wrap(~ Site) +
  ylab("log(1 + DNA concentration)")+
  xlab("Species")+
  scale_y_continuous(trans="log1p")+
  theme_bw() + theme(legend.position = "none")+
  scale_x_discrete(labels = c(Te_mar = "T. maritimum", Te_fin = "T. finnmarkense", sch = "Cand. S. salmonis",
                              sasa = "Atlantic salmon", pisck_sal = "P. salmonis", pa_ther = "P. theridion",
                              env = "Erythrocytic necrosis virus"), guide = guide_axis(angle = 45))


envdata %>% ggplot(aes(x=Site, y=temp)) +
  geom_boxplot()

envdata %>% ggplot(aes(x=Site, y=salinity))+
  geom_boxplot()

envdata %>% ggplot(aes(x=Site, y=turbidity))+
  geom_boxplot()

#Model fitting ----
env.formula <- bf(
  Conc ~ Site + PrevConc + temp + salinity + turbidity + (1 + Site | Target) + (1|Sample),
  hu ~ Site + PrevConc + temp + salinity + turbidity + (1 + Site | Target) + (1|Sample)
)


family <- hurdle_gamma()

# what's the default prior look like?
get_prior(env.formula, data = envdata_nc, family = family)
# -

fit <- brm(
  env.formula, data = envdata_nc, family = family,
  chains = 6, iter = 4000, warmup = 2000, seed = 1,
  core = 6, control = list(max_treedepth = 10, adapt_delta = 0.9),
  prior = c(set_prior("normal(0, 100)", class = "b")),  # add a wide prior to the coefficients
  silent = 0
)
summary(fit)

tidyfit <- tidy(fit)

#Conditional effects ----
library(ggpubr)
#Temp
ce_temp <- plot(conditional_effects(fit, "temp"), plot=FALSE)[[1]]$data %>% 
  ggplot(aes(x = effect1__, y = estimate__*1000)) +
  geom_ribbon(aes(ymin = lower__*1000, ymax = upper__*1000), alpha = 0.3, colour = NA) +
  geom_line(size = 1) +
  xlab("Seawater temperature °C") +
  ylab("Expected eDNA \n concentration (copies/μL)") +
  theme_bw(base_size = 18)

#What is the effect of temp on zeros #As salinity increases, so do zeros
ce_hu_temp <- plot(conditional_effects(fit, "temp", dpar = "hu"), plot = FALSE)[[1]]$data %>% 
  ggplot(aes(x = effect1__, y = estimate__)) +
  geom_ribbon(aes(ymin = lower__, ymax = upper__), alpha = 0.3, colour = NA) +
  geom_line(size = 1) +
  xlab("Seawater temperature °C") +
  ylab("Predicted probability of \n zero eDNA detection") +
  theme_bw(base_size = 18)
hypothesis(fit, "temp > 0")

ggarrange(ce_temp, ce_hu_temp,
          labels = c("a", "b"),
          ncol = 1, nrow = 2)

#Salinity
ce_sal <- plot(conditional_effects(fit, "salinity"), plot=FALSE)[[1]]$data %>% 
  ggplot(aes(x = effect1__, y = estimate__*1000)) +
  geom_ribbon(aes(ymin = lower__*1000, ymax = upper__*1000), alpha = 0.3, colour = NA) +
  geom_line(size = 1) +
  xlab("Salinity (ppt)") +
  ylab("Expected eDNA \n concentration (copies/μL)") +
  theme_bw(base_size = 18)
#What is the effect of salinity on zeros #No effect!
ce_hu_sal <- plot(conditional_effects(fit, "salinity", dpar = "hu"), plot = FALSE)[[1]]$data %>% 
  ggplot(aes(x = effect1__, y = estimate__)) +
  geom_ribbon(aes(ymin = lower__, ymax = upper__), alpha = 0.3, colour = NA) +
  geom_line(size = 1) +
  xlab("Salinity (ppt)") +
  ylab("Predicted probability of \n zero eDNA detection") +
  theme_bw(base_size = 18)

ggarrange(ce_sal, ce_hu_sal,
          labels = c("a", "b"),
          ncol = 1, nrow = 2)
hypothesis(fit, "salinity > 0")

#Turbidity
ce_turb <- plot(conditional_effects(fit, "turbidity"), plot=FALSE)[[1]]$data %>% 
  ggplot(aes(x = effect1__, y = estimate__*1000)) +
  geom_ribbon(aes(ymin = lower__*1000, ymax = upper__*1000), alpha = 0.3, colour = NA) +
  geom_line(size = 1) +
  xlab("Secchi depth (m)") +
  ylab("Expected eDNA \n concentration (copies/μL)") +
  theme_bw(base_size = 18)

ce_hu_turb <- plot(conditional_effects(fit, "turbidity", dpar = "hu"), plot = FALSE)[[1]]$data %>% 
  ggplot(aes(x = effect1__, y = estimate__)) +
  geom_ribbon(aes(ymin = lower__, ymax = upper__), alpha = 0.3, colour = NA) +
  geom_line(size = 1) +
  xlab("Secchi depth (m)") +
  ylab("Predicted probability of \n zero eDNA detection") +
  theme_bw(base_size = 18)

ggarrange(ce_turb, ce_hu_turb,
          labels = c("a", "b"),
          ncol = 1, nrow = 2)

hypothesis(fit, "turbidity > 0")

#
ce_pc <- plot(conditional_effects(fit, "PrevConc"), plot=FALSE)[[1]]$data %>% 
  ggplot(aes(x = effect1__*1000, y = estimate__*1000)) +
  geom_ribbon(aes(ymin = lower__*1000, ymax = upper__*1000), alpha = 0.3, colour = NA) +
  geom_line(size = 1) +
  xlab("eDNA concentraction \n in previous sample (copies/μL)") +
  ylab("Expected eDNA \n concentration (copies/μL)") +
  theme_bw(base_size = 18)

ce_hu_pc <- plot(conditional_effects(fit, "PrevConc", dpar = "hu"), plot = FALSE)[[1]]$data %>% 
  ggplot(aes(x = effect1__*1000, y = estimate__)) +
  geom_ribbon(aes(ymin = lower__, ymax = upper__), alpha = 0.3, colour = NA) +
  geom_line(size = 1) +
  xlab("eDNA concentraction \n in previous sample (copies/μL)") +
  ylab("Predicted probability of \n zero eDNA detection") +
  theme_bw(base_size = 18)

ggarrange(ce_pc, ce_hu_pc,
          labels = c("a", "b"),
          ncol = 1, nrow = 2)
hypothesis(fit, "PrevConc > 0")


#Pairs plots ----
pairs(fit, variable = c("b_Intercept", "b_SiteKnight", "b_temp", "b_salinity", "b_turbidity"))

pairs(fit,
      variable = c("shape", "b_Intercept", "b_SiteKnight", "b_hu_Intercept"),
      # variable = c("shape", "b_Intercept", "b_SiteKnight", "b_PrevConc", "b_hu_Intercept"),
      diag_fun="dens",
      off_diag_fun="hex",
)

pairs(fit,
      variable = c("b_hu_Intercept", "b_hu_SiteKnight", "r_Target__hu\\[.*,Intercept\\]"),
      regex = TRUE,
      diag_fun="dens",
      off_diag_fun="hex",
)

pairs(fit,
      variable = c("b_hu_Intercept", "b_hu_SiteKnight", "r_Target__hu\\[.*,SiteKnight\\]"),
      regex = TRUE,
      diag_fun="dens",
      off_diag_fun="hex",
)

# The fixed effect and shape parameter pair plots look fine.
# The random effects appear totally degenerate with the fixed effect (not unexpected, we're effectively fitting for the effect of site in four different places in this model!) and thus with eachother.

# a quick posterior predictive check, to make sure we're near the right track
pp_check(fit, ndraws = 20) + scale_x_continuous(trans = "log1p") + theme_bw()

#Extract draws----
drawsf <- spread_draws(fit, 
                       b_Intercept,
                       b_SiteKnight,
                       b_PrevConc,
                       b_hu_Intercept,
                       b_hu_SiteKnight,
                       b_hu_PrevConc,
                       b_hu_temp,
                       b_hu_salinity,
                       b_hu_turbidity,
                       b_temp)


#Species posterior plots  ----
#Just mu
fit_draws %>% ggplot(aes(x=mu, color=Target, linetype=Site)) +
  geom_density() +
  scale_x_log10() +
  labs(x = "Expected DNA conc.", y = "Density") + theme_bw() + theme(legend.position = "bottom")

# hurdle parameter has way too much freedom, which leads to these silly posteriors:
fit_draws %>% ggplot(aes(x=hu, color=Target, linetype=Site)) +
  geom_density() +
  #scale_x_log10() +
  scale_y_log10() +
  coord_cartesian(ylim=c(1e-1, 1e+3)) +
  labs(x = "Hurdle (zero-inflation) probability", y = "Density") + theme_bw() + theme(legend.position = "bottom")

# mu and hu species effects
newdata <- expand_grid(
  Site = c("Clio", "Knight"),
  Target = sort(unique(envdata$Target)),
  PrevConc = mean(envdata$PrevConc),
  temp = mean(envdata$temp),
  salinity = mean(envdata$salinity),
  turbidity = mean(envdata$turbidity, na.rm = TRUE),
  Sample = (unique(envdata$Sample))
)
epred <- epred_draws(fit, newdata=newdata, allow_new_levels = TRUE)


epred.long <- epred %>% pivot_wider(id_cols = c(Target, Sample, .draw), names_from = Site, values_from = .epred) %>% 
  mutate(ratio = Clio/Knight) %>% 
  mutate(fullname = case_when(Target == "Te_mar" ~ "T. maritimum", 
                              Target == "Te_fin" ~ "T. finnmarkense", 
                              Target == "sch" ~ "Cand. S. salmonis",
                              Target == "sasa" ~ "Atlantic salmon", 
                              Target == "pisck_sal" ~ "P. salmonis", 
                              Target == "pa_ther" ~ "P. theridion",
                              Target == "env" ~ "ENV")
  )

ggplot(data = epred.long) +
  geom_density(aes(x= ratio, color=fullname), size = 0.75, trim = FALSE) +
  #geom_density(aes(x=ratio), size = 0.5, linetype = 2, trim = FALSE) +
  xlim(-1, 15) +
  labs(x = "Ratio of expected DNA concentration at \n Clio vs. Knight", y = "Density") + 
  theme_bw() + 
  theme(legend.position = "bottom", legend.title = element_text(size=14),
        legend.text = element_text(size=14),
        axis.text=element_text(size=12),
        axis.title=element_text(size=14))+
  labs(color = "Target species")+
  geom_vline(xintercept=1, linetype = 2)





####Analysis of Clio channel and knight inlet data collected summer 2024#
#Bayesian regression models using stan. The goal is to model species distributions with species ID as a random effect
#Created November 18, 2024
#Last updated January 6, 2025

#Housekeeping ----
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
library(scales)
library(ggtext)
library(repr)
library(marginaleffects)
options(repr.plot.width=12, repr.plot.height=8, repr.plot.dpi=300)

#Load and prepare data ----
source("~/Documents/GitHub/clio/scripts/data_preparation.R")

data_nc <- data %>% filter(Site != "Control")
#Plot data ----
#Histograms by species
data_nc %>% ggplot(aes(x=Conc*100, fill = Target))+
  geom_histogram(bins = 10) +
  facet_grid(Site ~ Target, scales = "free") +
  theme_bw() + theme(legend.position = "none")

#Boxplots by species
data_nc %>% ggplot(aes(x=Target, y=Conc*100))+
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


#Model fitting ----
formula <- bf(
    Conc ~ Site + PrevConc + (1 + Site | Target) + (1|Sample),
    hu ~ Site + PrevConc + (1 + Site | Target) + (1|Sample)
)

family <- hurdle_gamma()

# what's the default prior look like?
default_prior(formula, data = data_nc, family = family)
# -

fit1 <- brm(
    formula, data = data_nc, family = family,
    chains = 6, iter = 4000, warmup = 2000, seed = 1,
    core = 6, control = list(max_treedepth = 10, adapt_delta = 0.9),
    prior = c(set_prior("normal(0, 100)", class = "b")),  # add a wide prior to the coefficients
    silent = 0
)
summary(fit1)

tidyfit <- tidy(fit1)

#Conditional effects ----

conditional_effects(fit, "Site")
hypothesis(envfit, "Site > 0")

conditional_effects(fit, "Site", method = "posterior_linpred")

fit %>% emmeans(~ Site, epred = FALSE) %>% contrast("pairwise")

hypothesis(fit, "SiteKnight > 0")

# For the effect of previous concentration
#Conditional effects method
conditional_effects(fit, "PrevConc")

#Marginal effects predictions method
conditional_pc <- predictions(
  fit, 
  newdata = datagrid(PrevConc = data_nc$PrevConc), 
  by = "PrevConc", 
  re_formula = NA
) %>% 
  posterior_draws() 
conditional_pc <- as.data.frame(conditional_pc)
conditional_pc %>% ggplot(aes(x=temp, y = estimate)) +
  geom_line()+
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = 0.3)+
  theme_bw() +
  ylab("Effect size") +
  xlab("eDNA concentration of previous sample (copies/μL")

hypothesis(fit, "PrevConc > 0")

#Pairs plots ----
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
drawsf <- spread_draws(fit1, 
  b_Intercept,
  b_SiteKnight,
  b_PrevConc,
  b_hu_Intercept,
  b_hu_SiteKnight,
  b_hu_PrevConc)

#Keep it simple for global means/hu
drawsf$B_clio <- exp(drawsf$b_Intercept) * (1 - plogis(drawsf$b_hu_Intercept))
drawsf$B_knight <- exp(drawsf$b_Intercept + drawsf$b_SiteKnight) * (1 - plogis(drawsf$b_hu_Intercept + drawsf$b_hu_SiteKnight))

drawsf$ratio <- drawsf$B_clio/drawsf$B_knight

#Global prev conc
drawsf$PC <- exp(drawsf$b_PrevConc + drawsf$b_Intercept) * (1-plogis(drawsf$b_hu_Intercept + drawsf$b_hu_PrevConc))

#Hurdle posteriors
#drawsf$HU_control <- (plogis(drawsf$b_hu_Intercept))
drawsf$HU_clio <- (plogis(drawsf$b_hu_Intercept))
drawsf$HU_knight <- (plogis(drawsf$b_hu_Intercept + drawsf$b_hu_SiteKnight))
drawsf$HU_ratio <- drawsf$HU_clio/drawsf$HU_knight

draws.longf <- drawsf %>% pivot_longer(everything())

#Posterior plots----
#Clio:Knight
draws.longf %>% 
  filter(name == "ratio") %>% 
  ggplot(aes(x=value, y=0, alpha = 0.5, color = "grey40")) +
  stat_halfeye(color = "grey20") +
  theme_bw()+
  coord_cartesian(xlim=c(0,10))+
  labs(x = "Proportion of zeros in DNA concentration, Clio:Knight", y = NULL)+
  guides(fill = "none", alpha = "none")+
  geom_vline(xintercept=1, linetype =2)

hypothesis(fit, "Intercept > SiteKnight")
#Clio is greater than knight at 86% probability

#Full model
draws.longf %>% 
  filter(name == "B_clio"| name == "B_knight") %>% # | name == "B_control") %>% 
  ggplot(aes(x=value, alpha = 0.5, fill=name)) +
  #stat_halfeye(color = "grey20", scale = 0.9) +
  stat_slabinterval(show_interval = FALSE) +
  theme_bw()+
  labs(x = "eDNA concentration", y= NULL, fill = "Site")+
  guides(alpha = "none") +
  coord_cartesian(xlim=c(0,15)) +
  scale_fill_discrete(labels = c(B_clio = "Clio channel", B_knight = "Knight inlet", B_control ="Control"), 
                      type = c("#1b9e77", "#d95f02", "#7570b3"))+
  theme(legend.position = c(0.8, 0.8), legend.background = element_rect(fill = alpha("white", 0)))

#Hurdle component only
draws.longf %>% 
  filter(name == "HU_clio"| name == "HU_knight") %>% 
  ggplot(aes(x=value, y=name,  alpha = 0.5, color = "grey40")) +
  stat_slabinterval(color="grey20") +
  theme_bw()+
  labs(x = "Proportion of zeros in data", y = NULL)+
  guides(alpha = "none") +
  scale_y_discrete(labels = c(HU_clio = "Clio channel", HU_knight = "Knight inlet", HU_control ="Control"))

#Species effects ----
fit_draws <- fit %>% spread_draws(
  b_Intercept,
  b_SiteKnight,
  b_PrevConc,
  b_hu_Intercept,
  b_hu_SiteKnight,
  b_hu_PrevConc,
  r_Target[Target,Site],
  r_Target__hu[Target,Site],
)

fit_draws_clio <- fit_draws %>%
  filter(Site == "Intercept") %>%
  mutate(
    r_Total = r_Target,
    r_Total__hu = r_Target__hu,
  )

fit_draws_knight <- fit_draws %>%
  filter(Site == "SiteKnight") 

fit_draws_knight$r_Total = fit_draws_clio$r_Target + fit_draws_knight$r_Target
fit_draws_knight$r_Total__hu = fit_draws_clio$r_Target__hu + fit_draws_knight$r_Target__hu

fit_draws <- bind_rows(fit_draws_knight, fit_draws_clio)

fit_draws <- fit_draws %>%
  mutate(
    Knight = as.numeric(Site == "SiteKnight"),
    logmu = b_Intercept + (b_SiteKnight * Knight) + r_Total,
    mu = exp(logmu),
    logithu = b_hu_Intercept + (b_hu_SiteKnight * Knight) + r_Total__hu,
    hu = plogis(logithu),
    su = 1 - hu
  )

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

#Just mu - percent change species effects
# mu_cl = exp(b_intercept + r_target[intercept])
# mu_kn = exp(mu_cl + b_siteknight + r_target[knight])
# 100 * (mu_cl / mu_kn - 1) = 100/exp(b_siteknight + r_target) - 100

ggplot(data = fit_draws_knight) +
  geom_density(aes(x=1 / exp(b_SiteKnight + r_Target), color=Target),
               size = 0.4) +
  geom_density(aes(x=1 / exp(b_SiteKnight + r_Target)), size = 1) +
  xlim(0, 5) +
  #coord_cartesian(xlim=c(0, 500)) +
  labs(x = "Ratio of expected DNA concentration at Clio vs. Knight", y = "Density") + theme_bw() + theme(legend.position = "bottom")

# mu and hu species effects
newdata1 <- expand_grid(
  Site = c("Clio", "Knight"),
  Target = sort(unique(data_nc$Target)),
  PrevConc = 0,
  Sample = unique(data_nc$Sample)
)
epred1 <- epred_draws(fit1, newdata=newdata)


epred1 <- epred1 %>% pivot_wider(id_cols = c(Target, .draw, Sample), names_from = Site, values_from = .epred) %>% 
        mutate(ratio = Clio/Knight) %>% 
        mutate(fullname = case_when(Target == "Te_mar" ~ "T. maritimum", 
                                    Target == "Te_fin" ~ "T. finnmarkense", 
                                    Target == "sch" ~ "Cand. S. salmonis",
                                    Target == "sasa" ~ "Atlantic salmon", 
                                    Target == "pisck_sal" ~ "P. salmonis", 
                                    Target == "pa_ther" ~ "P. theridion",
                                    Target == "env" ~ "ENV")
        )

ggplot(data = epred1) +
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
  
  
  

               
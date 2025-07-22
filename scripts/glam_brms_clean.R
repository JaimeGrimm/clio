####Analysis of Clio channel and knight inlet data collected summer 2024
#Written by Jaime Grimm
#Bayesian regression models using stan. 
#The goal is to model species distributions as a hierarchical model with species 
  #IDs as random effects
#Created November 18, 2024
#Last updated July 22, 2025

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
library(bayestestR)
options(repr.plot.width=12, repr.plot.height=8, repr.plot.dpi=300)

#Load and prepare data ----
source("./scripts/data_preparation.R")

data_nc <- data %>% filter(Site != "Control") #remove controls

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
data_nc$Site <- relevel(factor(data_nc$Site), "Knight") #Make Knight the baseline

formula <- bf(
    Conc ~ Site + PrevConc + (1 + Site | Target) + (1|Sample),
    hu ~ Site + PrevConc + (1 + Site | Target) + (1|Sample)
)

family <- hurdle_gamma()

# what's the default prior look like?
#default_prior(formula, data = data_nc, family = family)
# -

fit1 <- brm(
    formula, data = data_nc, family = family,
    chains = 6, iter = 4000, warmup = 2000, seed = 1,
    core = 6, control = list(max_treedepth = 10, adapt_delta = 0.9),
    prior = c(set_prior("normal(0, 100)", class = "b")),  # add a wide prior to the coefficients
    silent = 0
)
summary(fit1)

#Conditional effects ----

conditional_effects(fit1, "Site")
hypothesis(fit1, "Intercept > SiteKnight")

conditional_effects(fit, "Site", method = "posterior_linpred")

fit1 %>% emmeans(~ Site, epred = FALSE) %>% contrast("pairwise")

hypothesis(fit1, "SiteKnight > 0")

# For the effect of previous concentration
#Conditional effects method
conditional_effects(fit1, "PrevConc")

#Marginal effects predictions method - checking that this produces the same result
  #as the conditional effect function
conditional_pc <- predictions(
  fit1, 
  newdata = datagrid(PrevConc = data_nc$PrevConc), 
  by = "PrevConc", 
  re_formula = NULL
) %>% 
  posterior_draws() 
conditional_pc <- as.data.frame(conditional_pc)
conditional_pc %>% ggplot(aes(x=PrevConc, y = estimate)) +
  geom_line()+
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = 0.3)+
  theme_bw(base_size = 16) +
  ylab("Estimated eDNA concentration (copies/μL)") +
  xlab("eDNA concentration of previous sample (copies/μL)")

#If we were to deviate from Bayesian norms, is the effect of the PrevConc on the
  #subsequent sample concentration "significantly" greater than 0?
hypothesis(fit1, "PrevConc > 0")

#Pairs plots to check for correlation and divergent transitions----
pairs(fit,
      variable = c("shape", "b_Intercept", "b_SiteKnight", "b_PrevConc", "b_hu_Intercept"),
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
# The random effects appear totally degenerate with the fixed effect 
  #(not unexpected, we're effectively fitting for the effect of site in four 
  #different places in this model!) and thus with eachother.

# a quick posterior predictive check, to make sure we're near the right track
pp_check(fit1, ndraws = 20)  + scale_x_continuous(trans = "log1p") + theme_bw()

#Extract draws----
drawsf <- spread_draws(fit1, 
  b_Intercept,
  b_SiteClio,
  b_PrevConc,
  b_hu_Intercept,
  b_hu_SiteClio,
  b_hu_PrevConc)

#Keep it simple for global means/hu.
drawsf$B_clio <- exp(drawsf$b_Intercept) * (1 - plogis(drawsf$b_hu_Intercept))
drawsf$B_knight <- exp(drawsf$b_Intercept + drawsf$b_SiteKnight) * (1 - plogis(drawsf$b_hu_Intercept + drawsf$b_hu_SiteKnight))

#Ratio clio/knight full model
drawsf$ratio <- drawsf$B_clio/drawsf$B_knight

#Hurdle posteriors
#drawsf$HU_control <- (plogis(drawsf$b_hu_Intercept))
drawsf$HU_clio <- (plogis(drawsf$b_hu_Intercept))
drawsf$HU_knight <- (plogis(drawsf$b_hu_Intercept + drawsf$b_hu_SiteKnight))
drawsf$HU_ratio <- drawsf$HU_clio/drawsf$HU_knight

#drawsf <- drawsf %>% filter(B_knight != 0)

draws.longf <- drawsf %>% pivot_longer(everything())

#Posterior plots----

#Clio:Knight
draws.longf %>% 
  filter(name == "ratio") %>% 
  ggplot(aes(x=value, y=0, alpha = 0.5, color = "grey40")) +
  stat_halfeye(color = "grey20", point_interval = "mean_qi",  .width = c(0.66, 0.95)) +
  theme_bw(base_size = 16)+
  coord_cartesian(xlim=c(0,10))+
  labs(x = "Clio:Knight", y = NULL)+
  guides(fill = "none", alpha = "none")+
  geom_vline(xintercept=1, linetype =2)

#Calculate how much of the posterior ratio falls above 1
drawsf %>% mutate(pos = ratio > 1) %>% summarize(100*mean(pos)) 
#79.6%

#Find mean and credible intervals to from ratio distribution:
mean_qi(drawsf$ratio)

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
  stat_slabinterval(color="grey20", point_interval = "mean_qi") +
  theme_bw(base_size = 16)+
  labs(x = "Proportion of zeros in data", y = NULL)+
  guides(alpha = "none") +
  scale_y_discrete(labels = c(HU_clio = "Clio channel", HU_knight = "Knight inlet", HU_control ="Control"))

# mu and hu species effects
newdata1 <- expand_grid(
  Site = c("Clio", "Knight"),
  Target = sort(unique(data_nc$Target)),
  PrevConc = 0,
  Sample = unique(data_nc$Sample)
)
epred1 <- epred_draws(fit1, newdata=newdata1)

epred1.long <- epred1 %>% pivot_wider(id_cols = c(Target, .draw, Sample), names_from = Site, values_from = .epred) %>% 
        mutate(ratio = Clio/Knight) %>% 
        mutate(fullname = case_when(Target == "Te_mar" ~ "T. maritimum", 
                                    Target == "Te_fin" ~ "T. finnmarkense", 
                                    Target == "sch" ~ "Cand. S. salmonis",
                                    Target == "sasa" ~ "Atlantic salmon", 
                                    Target == "pisck_sal" ~ "P. salmonis", 
                                    Target == "pa_ther" ~ "P. theridion",
                                    Target == "env" ~ "ENV")
        )
#Calculate how much more Atlantic salmon DNA there is in Clio than Knight
#We have to deal with long tails on the posteriors making it difficult to visualize
#Lets plot the 95% CI of the posterior rather than the full posterior
ci(epred1.long$ratio, method = "HDI")
#The upper CI limit when considering all species together is 11.71.
#To be conservative, what's the max upper CI amongst all species?
targets <- unique(epred1.long$Target)
cis <- list()
for (i in 1:length(targets)){
cis[[i]] <- ci(epred1.long$ratio[epred1.long$Target == targets[i]], method = "HDI")
}

#The largest upper CI is for Atlantic salmon, with 95% CIs between 1.35 and 30.2.
# We will plot the 95% CI distribution for each species

foo <- epred1.long %>% filter(case_when(Target == "env" ~ between(ratio, 0.48, 2.09),
                                        Target == "pa_ther" ~ between(ratio, 0.54, 2.16),
                                        Target == "pisck_sal" ~ between(ratio, 0.63, 1.85),
                                        Target == "sasa" ~ between(ratio, 1.35, 30.23),
                                        Target == "sch" ~ between(ratio, 0.58, 1.72),
                                        Target == "Te_fin" ~ between(ratio, 0, 13.04),
                                        Target == "Te_mar" ~ between(ratio, 0, 3.69)
                                          ))

ggplot(data = foo) +
  #stat_halfeye(aes(x = ratio, color = fullname)) +
  geom_density(aes(x= ratio, color=fullname), size = 0.75, trim = FALSE) +
  #geom_density(aes(x=ratio), size = 0.5, linetype = 2, trim = FALSE) +
  #xlim(-1, 15) +
  #ylim(0, 1.6) +
  #labs(x = NULL, y = NULL, title = "Without environmental variables") +
  labs(x = "Ratio of expected DNA concentration at \n Clio vs. Knight", y = "Density") + 
  theme_bw() + 
  theme(legend.position = "bottom") +
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14))+
  labs(color = "Target species")+
  geom_vline(xintercept=1, linetype = 2)
  

# # Analysis of Clio channel and Knight inlet data collected summer 2024
#
# Bayesian regression models using stan. The goal is to model species distributions with species ID as a random effect.

# ## preamble

# +
library(tidyr)
library(tidyverse)
library(bayesplot)
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
options(repr.plot.width=12, repr.plot.height=8, repr.plot.dpi=300)
# -


# ## loading & cleaning data

# +
# using readr's read_csv for better typing
data1 <- read_csv("./clio.csv", skip = 42, show_col_types = FALSE)
data2 <- read_csv("./clio2.csv", skip = 42, show_col_types = FALSE)
data3 <- read_csv("./clio3.csv", skip = 42, show_col_types = FALSE)

data <- bind_rows(data1, data2, data3) %>%
  # simplify column names:
  rename(Target = `Target Name`, Sample = `Sample Name`, Amp = `Amp Status`) %>%
  # add `Site` column:
  mutate(Site = case_when(
      str_detect(Sample, "CL")~"Clio",
      str_detect(Sample, "KN")~"Knight",
      TRUE~'Control')) %>%
  # remove standards and extraneous samples:
  filter(
      Task != "STANDARD",
      Target != "ascv",
      Target != "ctv-2",
      Sample != "EB0108",
      Sample != "EB0108.1",
      Sample != "CT0804",
      Sample != "qPCRBlank",
      Sample != "qPCRblank",
      Sample != "Messedup",
      Sample != "water") %>%
  # grab sequence numbers:
  mutate(Sequence = as.numeric(substr(str_extract(Sample, "[[:digit:]]+"), 1, 2))) %>%
  # controls were run at the beginning & end of each day; increment second day's sequence number by 2:
  mutate(Sequence = if_else(Sequence < 9, Sequence, Sequence + 2)) %>%
  # grab only the columns we need:
  select(Sample, Sequence, Site, Target, Quantity, Amp)

# Manually assign sequence numbers for controls:
data$Sequence[data$Sample=="CT0884"] <- 0 #beginning of day 1
data$Sequence[data$Sample=="CT0837"] <- 9 #end of day 1
data$Sequence[data$Sample=="CT0931"] <- 10 #beginning of day 2
data$Sequence[data$Sample=="CT0917"] <- 27 #end of day 2
# -

# ### adding a `PrevConc` covariate
#
# There's some fear that previously-taken samples may 'contaminate' the subsequent samples (i.e., autocorrelated noise in sample sequence).
# To accomodate this possible effect in the model, we add a covariate corresponding to the (nan-)mean of the previous samples' observed concentrations.
# (We need to average over the triplicated samples in the previous sequence.)
#
# We also drop some samples at this point:
#
# + `Control` samples, whose fitted parameters are of no interest to us
# + the `KN1922` sample targeting "Te_mar", which is a significant outlier
## JG note - Jan 6. Rather than just removing the outlier at this point, I'm removing all failed amplifications that still have a non-zero quantity associated with them
# +
data$Conc <- if_else(is.na(data$Quantity), 0, data$Quantity)

# check out avg. DNA concentration against sequence number
data %>% group_by(Target, Sequence) %>%
    summarise(MeanConc=mean(Conc, na.rm = TRUE), .groups = "keep") %>%
    ggplot(aes(x=Sequence, y=MeanConc, color=Target)) + geom_line() + geom_point() + scale_y_continuous(trans = "log1p")

# +
data$PrevConc = 0

for (t in unique(data$Target)) {
    for (s in 1:max(data$Sequence)) {
        pq <- mean(
            na.omit(filter(data, Sequence == s - 1, Target == t)$Quantity)
        )
        data$PrevConc[(data$Sequence == s) & (data$Target == t)] <- if_else(is.na(pq), 0, pq)
    }
}

data$PrevConc[data$Sequence == 10] <- 0  #decontaminated before start of day 2

data_nc <- data %>%
    filter(Site != "Control") %>% # filter out control trials, for fitting
    #filter(!(Sample == "KN1922" & Target == "Te_mar")) %>% # drop the outlier sequence
    filter(!(Amp == "No Amp" & Conc > 0)) %>% 
    mutate(Conc = Conc / 100, PrevConc = PrevConc / 100) # scale for easier sampling
# -

# ## plotting data

#Histograms by species
data_nc %>% ggplot(aes(x=Conc, fill = Target))+
  geom_histogram(bins = 10) +
  facet_grid(Site ~ Target, scales = "free") +
  theme_bw() + theme(legend.position = "none")

#Boxplots by species
data_nc %>% ggplot(aes(x=Target, y=Conc))+
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


# ## model fitting

# +
formula <- bf(
    Conc ~ Site + PrevConc + (1 + Site | Target),
    hu ~ Site + PrevConc + (1 + Site | Target)
)

family <- hurdle_gamma()

# what's the default prior look like?
default_prior(formula, data = data_nc, family = family)
# -

fit <- brm(
    formula, data = data_nc, family = family,
    chains = 6, iter = 4000, warmup = 2000, seed = 1,
    core = 6, control = list(max_treedepth = 10, adapt_delta = 0.9),
    prior = c(set_prior("normal(0, 100)", class = "b")),  # add a wide prior to the coefficients
    silent = 0
)
#threads = threading(2),
summary(fit)

# ## conditional effects

conditional_effects(fit, "Site", method = "posterior_linpred")

fit %>% emmeans(~ Site, epred = FALSE) %>% contrast("pairwise")

hypothesis(fit, "SiteKnight > 0")

conditional_effects(fit, "PrevConc", method = "posterior_linpred")

hypothesis(fit, "PrevConc > 0")

# ### notes on conditional effects
#
# + The `Site` marginal effect posterior's median is about exp(0.2) ~= 1.2 times higher at Clio than Knight (i.e., the "best fit" `b_SiteKnight ~= -0.24`), so we expect ~30% more DNA (in expectation, not observation!) at Clio vs. Knight. (This difference is driven almost entirely by the `sasa` and `Te_fin` assays, the others don't see the same effect size.) The 95% CI for the `b_SiteKnight` coefficient is still consistent with 0, though.
# + The `PrevConc` marginal effect is statistically (at 99% credibility) distinct from zero, although the effect *size* is consistant with 0 for the vast majority of data samples (where PrevConc < 10; see the conditional effect plot above). This is likely an artifact of the outlying "KN1922" sample, which we've dropped from the fit but left in the `PrevConc` calculation; we should do something more clever there

# ## pairs plots

pairs(fit,
      variable = c("shape", "b_Intercept", "b_SiteKnight", "b_hu_Intercept"),
      # variable = c("shape", "b_Intercept", "b_SiteKnight", "b_PrevConc", "b_hu_Intercept"),
      diag_fun="dens",
      off_diag_fun="hex",
)

pairs(fit,
      variable = c("b_Intercept", "b_SiteKnight", "r_Target\\[.*,Intercept\\]"),
      regex = TRUE,
      diag_fun="dens",
      off_diag_fun="hex",
)

pairs(fit,
      variable = c("b_Intercept", "b_SiteKnight", "r_Target\\[.*,SiteKnight\\]"),
      regex = TRUE,
      diag_fun="dens",
      off_diag_fun="hex",
)

# The fixed effect and shape parameter pair plots look fine.
# The random effects appear totally degenerate with the fixed effect (not unexpected, we're effectively fitting for the effect of site in four different places in this model!) and thus with eachother.
# Should probably revise the model to have clearer dependence on both target and site.

# a quick posterior predictive check, to make sure we're near the right track
pp_check(fit, ndraws = 20) + scale_x_continuous(trans = "log1p") + theme_bw()

# ## other posterior plots

# +
# consider both sites and all target assays, assume no contamination (PrevConc = 0)
newdata <- expand_grid(
    Site = c("Clio", "Knight"),
    Target = sort(unique(data$Target)),
    PrevConc = 0,
)

linpred <- linpred_draws(fit, newdata = newdata)
postpred <- predicted_draws(fit, newdata=newdata)
epred <- epred_draws(fit, newdata=newdata)

ggplot() +
    geom_histogram(data = data_nc, aes(x = Conc, y = after_stat(density)), fill = "grey80", bins = 10) + 
    geom_density(data = data_nc, aes(x = Conc,  linetype = "Empirical DNA density")) +
    geom_density(data = postpred, aes(x = .prediction, color = Target, linetype = "Posterior prediction")) +
    #labs(x = "Expected DNA concentration", fill = "Site") +
    theme_bw() + theme(legend.position = "bottom") +
    facet_grid(rows = vars(Site), cols = vars(Target)) +
    coord_cartesian(xlim = c(0,10), ylim = c(0, 10)) +
    labs(linetype="Line type")
    
# -
# ## marginal effect of site?
#
# We're mainly interested in the marginal effect of "site" (Clio vs. Knight).
# The model makes this a little tough to answer, since the site variable enters in four places:
# + the mean site coefficient (e.g., `b_SiteKnight`)
# + the mean site random effect (which is conditioned on target assay: `r_Target[Site,Target]`)
# + the hurdle coefficient (e.g.,`b_hu_SiteKnight`)
# + the hurdle site random effect (also conditioned on target: `r_Target__hu[Site,Target]`)

# +
newdata_clio <- newdata %>% filter(Site == "Clio")
newdata_knight <- newdata %>% filter(Site == "Knight")

lpd <- fit %>%
    linpred_draws(newdata_clio) %>%
    rename(logmu_clio=.linpred)#

lpd$logmu_knight <- linpred_draws(fit, newdata_knight)$.linpred
lpd <- lpd %>% mutate(
    logmu_diff = logmu_clio - logmu_knight,
    mu_ratio = exp(logmu_diff)
    )
# -

lpd %>%
    ggplot(aes(x = 100 * mu_ratio - 100)) +
    geom_density(aes(color = Target)) + 
    geom_density(color="black") + 
    labs(x = "percent difference in expected DNA concentration at Clio vs. Knight", y = "posterior density", fill = "Site") +
    xlim(-100,400) +
    theme_bw() + theme(legend.position = "bottom")

fit %>% emmeans(~ Site, epred = TRUE) %>% contrast("pairwise") %>% 
    gather_emmeans_draws() %>%
    ggplot(aes(x = 100 * exp(.value) - 100)) +
    stat_halfeye() +
    labs(x = "Percent increase in expected DNA concentration at Clio vs. Knight", y = "posterior density") +
    xlim(-100, 400) +
    theme_bw() + theme(legend.position = "bottom")

# ### can we get the same thing directly from the MCMC draws?
#
# The above plots rely on the `linpred_draws` function.
# To make sure I know what it does, I'll try to recover the same results directly from the MCMC posterior draws.

# +
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
    mutate(r_Total = r_Target, r_Total__hu = r_Target__hu)

fit_draws_knight <- fit_draws %>%
    filter(Site == "SiteKnight") 

fit_draws_knight$r_Total = fit_draws_knight$r_Target + fit_draws_clio$r_Total
fit_draws_knight$r_Total__hu = fit_draws_knight$r_Target__hu + fit_draws_clio$r_Total__hu

fit_draws <- bind_rows(fit_draws_knight, fit_draws_clio)

fit_draws <- fit_draws %>%
    mutate(
        Knight = as.numeric(Site == "SiteKnight"),
        logmu = b_Intercept + (b_SiteKnight * Knight) + r_Total,
        mu = exp(logmu),
        logithu = b_hu_Intercept + (b_hu_SiteKnight * Knight) + r_Total__hu,
        hu = plogis(logithu)
    )
# -

fit_draws %>% ggplot(aes(x=mu, color=Target, linetype=Site)) +
    geom_density() +
    scale_x_log10() +
    labs(x = "Expected DNA conc.", y = "Density") + theme_bw() + theme(legend.position = "bottom")

# hurdle parameter has way too much freedom, which leads to these silly posteriors:
fit_draws %>% ggplot(aes(x=hu, color=Target, linetype=Site)) +
    geom_density() +
    # scale_x_log10() +
    scale_y_log10() +
    coord_cartesian(ylim=c(1e-1, 1e+3)) +
    labs(x = "Hurdle (zero-inflation) probability", y = "Density") + theme_bw() + theme(legend.position = "bottom")

ggplot(data = fit_draws_knight) +
    geom_density(aes(x=100 / exp(b_SiteKnight + r_Target) - 100, color=Target),
                 size = 0.4) +
    geom_density(aes(x=100 / exp(b_SiteKnight + r_Target) - 100), size = 1) +
    xlim(-100, 500) +
    labs(x = "Percent increase in expected DNA concentration at Clio vs. Knight", y = "Density") + theme_bw() + theme(legend.position = "bottom")




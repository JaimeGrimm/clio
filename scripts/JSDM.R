#Joint species distribution modeling in HMSC
#Clio and knight inlet data collected in summer 2024
# Last updated November 18, 2024

library(tidyr)
library(lme4)
library(brms)
#Load clean data
source("~/Documents/GitHub/clio/scripts/data_preparation.R")

####CLEAN DATA####
rawdata <- filter(data, Target.Name != "ascv" & Target.Name != "ctv-2") %>% #remove assays that didn't amplify
        mutate(dummy = recode(site, "control" = 0,
                                    "knight" = 1,
                                    "clio" = 2))
nocontrol <- filter(data, Target.Name != "ascv" & Target.Name != "ctv-2",
                    site != "control") %>%
            dplyr::select(site, seq, Quantity, Target.Name, Sample.Name) %>% 
             mutate(dummy = recode(site, "knight" = 0,
                        "clio" = 1)) %>% 
              mutate_if(is.numeric, round)
  
####MODEL 1 - species as random effects, no temporal autocorrelation----
#*Model 1 - Hmsc----
###Organize data into right format for Hmsc###
XData1 <- data.frame(site=nocontrol$dummy)
YData1 <- nocontrol$Quantity

studyDesign1 <- data.frame(sample = as.factor(nocontrol$Sample.Name), species = as.factor(nocontrol$Target.Name))
rl <- HmscRandomLevel(units = studyDesign1$species)

nChains <- 2
thin <- 5
samples <- 1000
transient <- 500*thin
verbose <- 500*thin

#Define and fit model

m1 <- Hmsc(Y=YData1, XData = XData1, XFormula = ~site, distr = "lognormal poisson",
          studyDesign = studyDesign1, ranLevels = list(sample = rl))
m1 <- sampleMcmc(m1, thin = thin, samples = samples, transient = transient,
                nChains = nChains, verbose = verbose)

#Pull out posterior and model fit
M1F <- evaluateModelFit(hM = m1, predY = preds)
M1F$SR2
postBeta1 <- getPostEstimate(m1, parName = "Beta")
#extract posterior estimates
mpost1 <- convertToCodaObject(m1) 

#compute posterior distribution of predicted values
preds1 <- computePredictedValues(m1)

#Checking convergence
plot(mpost1$Beta)
effectiveSize(mpost1$Beta)
gelman.diag(mpost1$Beta, multivariate = FALSE)$psrf

#Checking assumptions of linear model
preds.mean1 <- apply(preds1, FUN=mean, MARGIN=1)
nres1 <- scale(nocontrol$Quantity-preds.mean1)
par(mfrow = c(1,2))
hist(nres1)
plot(preds.mean1, nres1)

#Visualize
plotBeta(m1, post=postBeta, param = "Mean")


#MODEL 2----
##Pivot wider to separate species for fixed effects
data.wide <- pivot_wider(nocontrol, names_from = Target.Name, #Make wide
                         values_from = Quantity)
select.data <- dplyr::select(data.wide, env, pa_ther, pisck_sal, sasa, sch, Te_fin, Te_mar,
                            Sample.Name, dummy, seq) 
              
            
cdata <- select.data %>% group_by(Sample.Name, dummy, seq) %>%
            summarise(across(env:Te_mar, mean, na.rm = TRUE)) %>%
            ungroup() %>%
            mutate_if(is.numeric, round) #round data



###Organize data into right format for Hmsc###
XData <- data.frame(site=cdata$dummy)
YData <- dplyr::select(cdata, env:Te_mar)

studyDesign <- data.frame(as.factor(cdata$Sample.Name))
colnames(studyDesign) <- "sample"
rl <- HmscRandomLevel(units = as.factor(unique(cdata$seq)))

#Fit a basic model
nChains <- 2
thin <- 5
samples <- 1000
transient <- 500*thin
verbose <- 500*thin

m <- Hmsc(Y=YData, XData = XData, XFormula = ~site, distr = "lognormal poisson",
          studyDesign = studyDesign, ranLevels = list(sample = rl))
m <- sampleMcmc(m, thin = thin, samples = samples, transient = transient,
                nChains = nChains, verbose = verbose)

#extract posterior estimates
mpost <- convertToCodaObject(m1) 

#compute posterior distribution of predicted values
preds <- computePredictedValues(m1)


#Evalulate model fit
MF <- evaluateModelFit(hM = m1, predY = preds)
MF$SR2

postBeta <- getPostEstimate(m1, parName = "Beta")
plotBeta(m1, post=postBeta, param = "Mean")



data.long <- pivot_longer(cdata, cols= c(env, pa_ther, pisck_sal, sasa, sch, Te_fin, Te_mar), 
                          names_to = "assay") %>%
             mutate(site = recode(dummy, `0` = "knight",
                          `1` = "clio")) %>% 
            arrange(assay)

#Extract quantiles from posterior
qmat <- unname(summary(mpost$Beta)$quantiles) 
post <- data.frame(assay = unique(data.long$assay),
                     knight2.5 = c(qmat[1,1], qmat[3,1], qmat[5,1], qmat[7,1],
                                 qmat[9,1], qmat[11,1], qmat[13,1]),
                     knight97.5 = c(qmat[1,5], qmat[3,5], qmat[5,5], qmat[7,5],
                                  qmat[9,5], qmat[11,5], qmat[13,5]),
                     clio2.5 = c(qmat[2,1], qmat[4,1], qmat[6,1], qmat[8,1],
                                   qmat[10,1], qmat[12,1], qmat[14,1]),
                     clio97.5 = c(qmat[2,5], qmat[4,5], qmat[6,5], qmat[8,5],
                                    qmat[10,5], qmat[12,5], qmat[14,5]),
                    knight1 = postBeta$mean[1,],
                    clio1 = postBeta$mean[2,]
                                 )
rownames(post) <- NULL

post.long <- post %>% gather(key = "name", value = "value", 2:7 ) %>% 
             mutate(site = case_when(str_detect(name, "clio") ~ "clio",
                                        str_detect(name, "knight") ~ "knight")) %>% 
                mutate(quant = case_when(str_detect(name, "2.5")~ "lower",
                                        str_detect(name, "97.5")~ "upper",
                                        str_detect(name, "1")~ "mean",
                                        )) %>% 
                select(-c(name)) %>% 
                pivot_wider(names_from = quant, values_from = value)

                  


ggplot() +
  geom_jitter(data = data.long, aes(x=assay, y = log(value), col=site))+
  geom_point(data = foo1, aes(x=assay, y= mean, col=site))+
  geom_errorbar(data=foo1, aes(x=assay, ymin=lower, ymax=upper, col=site))

#Model constrution
#Three models exploring joint species distributions across two sites
#Model 1 is a lognormal Poisson
#Models 2 and 3 together are a hurdle model


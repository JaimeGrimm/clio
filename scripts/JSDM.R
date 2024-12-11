#Joint species distribution modeling in HMSC
#Clio and knight inlet data collected in summer 2024
# Last updated November 18, 2024

library(tidyr)
library(Hmsc)
library(reshape2)
#Load clean data
source("~/Documents/GitHub/clio/scripts/data_preparation.R")

####CLEAN DATA####
rawdata <- data %>% mutate(dummy = recode(site, "control" = 0,
                                    "knight" = 1,
                                    "clio" = 2))
nocontrol <- filter(data, Target.Name != "ascv" & Target.Name != "ctv-2",
                    site != "control") %>%
            dplyr::select(site, seq, Quantity, Target.Name, Sample.Name) %>% 
             mutate(dummy = recode(site, "knight" = 0,
                        "clio" = 1)) %>% 
              mutate_if(is.numeric, round)

####Set up model estimation regime----
nChains <- 2
thin <- 5
samples <- 2000
transient <- 1000
verbose <- 500*thin

####MODEL 1 - species as random effects, no temporal autocorrelation----
###Organize data into right format for Hmsc###
XData1 <- data.frame(site=nocontrol$dummy)
YData1 <- nocontrol$Quantity

studyDesign1 <- data.frame(sample = as.factor(nocontrol$Sample.Name), species = as.factor(nocontrol$Target.Name))
rl <- HmscRandomLevel(units = studyDesign1$species)

#Define and fit model

m1 <- Hmsc(Y=YData1, XData = XData1, XFormula = ~site, distr = "lognormal poisson",
          studyDesign = studyDesign1, ranLevels = list(sample = rl))
m1 <- sampleMcmc(m1, thin = thin, samples = samples, transient = transient,
                nChains = nChains, verbose = verbose)

#compute posterior distribution of predicted values
preds1 <- computePredictedValues(m1)

#Pull out posterior and model fit
M1F <- evaluateModelFit(hM = m1, predY = preds1)
M1F$SR2
postBeta1 <- getPostEstimate(m1, parName = "Beta")
#extract posterior estimates
mpost1 <- convertToCodaObject(m1) 

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
plotBeta(m1, post=postBeta1, param = "Mean")

#*Predictions----
gradient <- constructGradient(m1, focalVariable = "site")
predY1 <- predict(m1, Gradient = gradient, expected=TRUE)
plotGradient(m1, gradient, pred=predY1, measure="Y", index=1, showData = TRUE)

####MODEL 2----
#Species as fixed effects
##Pivot wider to separate species for fixed effects
data.wide <- nocontrol %>% dplyr::select(Sample.Name, Target.Name, Quantity, site, seq, dummy) %>% 
  group_by(Target.Name, Sample.Name, site, seq, dummy) %>% 
  dplyr::summarise(mean = mean(Quantity)) %>% ungroup() %>% 
  pivot_wider(names_from = Target.Name, values_from = mean) %>% 
  mutate_if(is.numeric, round) #round data
              
###Organize data into right format for Hmsc###
XData <- data.frame(site=data.wide$dummy)
YData <- dplyr::select(data.wide, Te_fin:sch)
data.wide$seq <- as.numeric(data.wide$seq)

xy = as.matrix(data.wide$seq)

studyDesign <- data.frame(sample = as.factor(data.wide$Sample.Name))
rownames(xy) = studyDesign[,1]
rl <- HmscRandomLevel(sData = studyDesign$time)

m <- Hmsc(Y=YData, XData = XData, XFormula = ~site, distr = "lognormal poisson",
          studyDesign = studyDesign, ranLevels = list(sample = rl))
m <- sampleMcmc(m, thin = thin, samples = samples, transient = transient,
                nChains = nChains, verbose = verbose)

#extract posterior estimates
mpost <- convertToCodaObject(m) 

#compute posterior distribution of predicted values
preds <- computePredictedValues(m)


#Evalulate model fit
MF <- evaluateModelFit(hM = m, predY = preds)
MF$SR2

postBeta <- getPostEstimate(m, parName = "Beta")
plotBeta(m, post=postBeta, param = "Mean")

data.long <- pivot_longer(data.wide, cols= c(env, pa_ther, pisck_sal, sasa, sch, Te_fin, Te_mar), 
                          names_to = "assay") %>%
             mutate(site = recode(dummy, `0` = "knight",
                          `1` = "clio")) %>% 
            arrange(assay)

#*Model diagnostics ----
effectiveSize(mpost$Beta)
hist(effectiveSize(mpost$Beta))
gelman.diag(mpost$Beta, multivariate = FALSE)$psrf

#*Plotting----
#Plot post estimates

#Get species and covariate names from model
spNames <- m$spNames
covNames <- recode(m$covNames, "(Intercept)" = "Knight", "site" = "Clio")

#Pull data from postBeta
mbeta <- postBeta$mean
betaP <- postBeta$support
supportlevel <- 0.95

#Fotmat the data as a sign and narrow to the covariates to desired support level
toPlot<- sign(mbeta)
toPlot <- toPlot * ((betaP > supportlevel) + (betaP < (1-supportlevel)) > 0)

#Format the data as a matrix and add column and row names
betaMat <- matrix(toPlot, nrow = m$nc, ncol= ncol(m$Y))
colnames(betaMat) <- spNames
rownames(betaMat) <- covNames

#remove intercept 
betaMat <- as.data.frame(betaMat) %>% 
  dplyr::slice(-1)

#reformat for ggplot
betaMatmelt <- as.data.frame(melt(as.matrix(betaMat)))

ggplot(betaMatmelt, aes(x = Var1, y = Var2, fill = factor(value))) +
  labs(x = "Site", y = "Species", fill = "Sign") +
  geom_tile(color = 'gray60')+
  theme(plot.background = element_rect(fill = "white"),
        legend.background = element_rect(fill = "white"),
        panel.border = element_rect(fill = NA, color = NA),
        legend.margin = margin(l = 1, unit = 'cm'),
        legend.title = element_text(hjust = 0.1),
        legend.key.width = unit(1, 'cm'),
        legend.key.height = unit(4, 'cm'),
        legend.text = element_text(size = 12),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))+
  scale_x_discrete(expand = c(0, 0)) + #because his x and y axis are numbers
  scale_y_discrete(expand = c(0, 0))

#Extract quantiles from posterior
qmat <- unname(summary(mpost$Beta)$quantiles) 
knight2.5 = c(qmat[1,1], qmat[3,1], qmat[5,1], qmat[7,1],
              qmat[9,1], qmat[11,1], qmat[13,1])
knight97.5 = c(qmat[1,5], qmat[3,5], qmat[5,5], qmat[7,5],
               qmat[9,5], qmat[11,5], qmat[13,5])
post <- data.frame(assay = unique(data.long$assay),
                    knight2.5 = knight2.5,
                    knight97.5 = knight97.5,
                    clio2.5 = knight2.5 + c(qmat[2,1], qmat[4,1], qmat[6,1], qmat[8,1],
                                   qmat[10,1], qmat[12,1], qmat[14,1]),
                    clio97.5 = knight97.5 + c(qmat[2,5], qmat[4,5], qmat[6,5], qmat[8,5],
                                    qmat[10,5], qmat[12,5], qmat[14,5]),
                    knight1 = postBeta$mean[1,],
                    clio1 = postBeta$mean[1,]+postBeta$mean[2,]
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
  geom_point(data = post.long, aes(x=assay, y= mean, col=site))+
  geom_errorbar(data=post.long, aes(x=assay, ymin=lower, ymax=upper, col=site))+
  labs(ylab = )

#Variance partitioning
VP <- computeVariancePartitioning(m)
plotVariancePartitioning(m, VP)

#Plot variance partitioning manually so it's less ugly:
VP.long <- as.data.frame(VP$vals) %>% rownames_to_column(var = "effect") %>% 
  pivot_longer(cols=!effect)

ggplot(VP.long, (aes(x= name, y=value, fill = effect))) +
  geom_bar(stat = "identity")+
  labs(y = "Variance proportion", x = "Species")
  theme_bw()

####Analysis of Clio channel and knight inlet data collected summer 2024
#Written by Jaime Grimm
#Bayesian regression models using stan. 
#The goal is to model eDNA concentration between the sites for each species individually 
#Last updated July 22, 2025


library(tidyr)
library(tidyverse)
library(brms)
library(tidybayes)
library(scales)
library(ggtext)
library(repr)
library(marginaleffects)
library(bayesplot)
options(repr.plot.width=12, repr.plot.height=8, repr.plot.dpi=300)

#Load data ----
source("./data_preparation.R")

#Model fitting ----
#Specify the model formulation
env.formula <- bf(
  Conc ~ Site + (1|Sample),
  hu ~ Site + (1|Sample)
)
family <- hurdle_gamma()

#Separate datasets by species
targets <- unique(envdata$Target) #target = unique species ID
datasets <- list()

for (j in 1:length(targets)){
  datasets[[j]] <- envdata_nc %>% filter(Target == targets[j])
}

#Fit for each species' dataset
fits <- list()
for (i in 1:length(datasets)){
  fits[[i]] <- fit <- brm(
    env.formula, data = datasets[[i]], family = family,
    chains = 6, iter = 4000, warmup = 2000, seed = 1,
    core = 6, control = list(max_treedepth = 10, adapt_delta = 0.9),
    prior = c(set_prior("normal(0, 100)", class = "b")),  # add a wide prior to the coefficients
    silent = 0
  )
}
names(fits) <- targets

#Check summary and posterior predictive one by one. Generally they look good, except for the Tenacibaculum spp 
#which is unsurprising given the data sparsity
summary(fits[[3]])

pp_check(fits[[7]], ndraws = 20) + scale_x_continuous(trans = "log1p") + theme_bw()


#Site effect plots ----
library(ggpubr)
library(PNWColors)
pal <- pnw_palette("Starfish")
targets
fullnames <- c("Atlantic salmon", "Tenacibaculum maritimum", "Paranucleospora theridion", "Tenacibaculum finnmarkense",
               "Candidatus Syngnamydia salmonis", "Piscirickettsia salmonis", "Erythrocytic necrosis virus")
label <- c("n = 72", "n = 72", "n = 69", "n = 72", "n = 72", "n = 72", "n = 72")

plotdata <- list()
plot <- list()

for (k in 1:length(targets)){
plotdata[[k]] <- 
  plot(conditional_effects(fits[[k]], "Site"), plot=FALSE, re_formula=NULL)[[1]]$data

plot[[k]] <- ggplot() +
  geom_jitter(data = datasets[[k]], aes(x = Site, y = Conc*100, color = Site), alpha = 0.5) +
  geom_point(data = plotdata[[k]], aes(x = effect1__, y = estimate__*100)) +
  geom_errorbar(data = plotdata[[k]], aes(x = effect1__, ymin = lower__*100, ymax = upper__*100)) +
  xlab(NULL) +
  ylab(NULL) +
  labs(title = fullnames[k]) +
  theme_bw() +
  theme(legend.position = "none")+
  scale_color_manual(values=c(pal[3], pal[7])) +
  annotate("text", x = 2.3, y = max(datasets[[k]]$Conc*100)*0.98, label = label[k])
  }

library(grid)
foo <- ggpubr::ggarrange(plot[[1]], plot[[2]], plot[[3]], plot[[4]], plot[[5]], plot[[6]], plot[[7]], align = "hv")# Temar_site, Tefin_site, align = "hv")
annotate_figure(foo, left = textGrob("eDNA concentration (copies/μL)", rot = 90, gp = gpar(fontsize =16)), 
                bottom = textGrob("Area", gp = gpar(fontsize =16)))


#Checking out environmental variables' effects----------------------------------
#Temperaure effects
#tempplotdata <- list()
#tempplot <- list()
#tempplotdatahu <- list()
#tempplothu <- list()
#
#for (k in 1:length(targets)){
#  tempplotdata[[k]] <- 
#    plot(conditional_effects(fits[[k]], "temp"), plot=FALSE, re_formula=NULL)[[1]]$data
#  tempplotdatahu[[k]] <- plot(conditional_effects(fits[[k]], "temp", dpar = "hu"), plot=FALSE, re_formula=NULL)[[1]]$data
#  
#  tempplot[[k]] <- ggplot() +
#   #geom_jitter(data = datasets[[k]], aes(x = temp, y = Conc*100, color = Site)) +
#    geom_line(data = tempplotdata[[k]], aes(x = effect1__, y = estimate__*100)) +
#    geom_ribbon(data = tempplotdata[[k]], aes(x = effect1__, ymin = lower__*100, ymax = upper__*100), alpha = 0.3) +
#    xlab(NULL) +
#    ylab(NULL) +
#    labs(title = fullnames[k]) +
#    theme_bw() +
#    theme(legend.position = "none")
#  
#  tempplothu[[k]] <- ggplot() +
#    geom_line(data = tempplotdatahu[[k]], aes(x = effect1__, y = estimate__)) +
#    geom_ribbon(data = tempplotdatahu[[k]], aes(x = effect1__, ymin = lower__, ymax = upper__), alpha = 0.3) +
#    xlab(NULL) +
#    ylab(NULL) +
#    #labs(title = fullnames[k]) +
#    theme_bw() +
#    theme(legend.position = "none")
#}
#
#temp_plot <- ggarrange(tempplot[[1]], tempplothu[[1]], tempplot[[3]], tempplothu[[3]], 
#                       tempplot[[4]], tempplothu[[4]], tempplot[[6]], tempplothu[[6]], 
#                       tempplot[[7]], tempplothu[[7]], ncol =2, nrow = 5, align = "hv")
#annotate_figure(temp_plot, left = textGrob("Expected eDNA concentration (copies/μL)", rot = 90, gp = gpar(fontsize =16)), 
#                bottom = textGrob("Seawater temperature °C", gp = gpar(fontsize =16)), 
#                right = textGrob("Predicted probability of zero eDNA detection", rot = 90, gp = gpar(fontsize =16)))
#
#Salinity effects
#salplotdata <- list()
#salplot <- list()
#salplotdatahu <- list()
#salplothu <- list()
#
#for (k in 1:length(targets)){
#  salplotdata[[k]] <- 
#    plot(conditional_effects(fits[[k]], "salinity"), plot=FALSE, re_formula=NULL)[[1]]$data
#  salplotdatahu[[k]] <- plot(conditional_effects(fits[[k]], "salinity", dpar = "hu"), plot=FALSE, re_formula=NULL)[[1]]$data
#  
#  salplot[[k]] <- ggplot() +
    #geom_jitter(data = datasets[[k]], aes(x = temp, y = Conc*100, color = Site)) +
#    geom_line(data = salplotdata[[k]], aes(x = effect1__, y = estimate__*100)) +
#    geom_ribbon(data = salplotdata[[k]], aes(x = effect1__, ymin = lower__*100, ymax = upper__*100), alpha = 0.3) +
#    xlab(NULL) +
#    ylab(NULL) +
#    labs(title = fullnames[k]) +
#    theme_bw() +
#    theme(legend.position = "none")
#  
#  salplothu[[k]] <- ggplot() +
#    geom_line(data = salplotdatahu[[k]], aes(x = effect1__, y = estimate__)) +
#    geom_ribbon(data = salplotdatahu[[k]], aes(x = effect1__, ymin = lower__, ymax = upper__), alpha = 0.3) +
#    xlab(NULL) +
#    ylab(NULL) +
#    #labs(title = fullnames[k]) +
#    theme_bw() +
#    theme(legend.position = "none")
#}
#
#sal_plot <- ggarrange(salplot[[1]], salplothu[[1]], salplot[[3]], salplothu[[3]], 
#                      salplot[[4]], salplothu[[4]], salplot[[6]], salplothu[[6]], 
#                      salplot[[7]], salplothu[[7]], ncol =2, nrow = 5, align = "hv")
#annotate_figure(sal_plot, left = textGrob("Expected eDNA concentration (copies/μL)", rot = 90, gp = gpar(fontsize =16)), 
#                bottom = textGrob("Salinity (ppt)", gp = gpar(fontsize =16)), 
#                right = textGrob("Predicted probability of zero eDNA detection", rot = 90, gp = gpar(fontsize =16)))
#
##Turbidity effects
#turplotdata <- list()
#turplot <- list()
#turplotdatahu <- list()
#turplothu <- list()

#for (k in 1:length(targets)){
#  turplotdata[[k]] <- 
#    plot(conditional_effects(fits[[k]], "turbidity"), plot=FALSE, re_formula=NULL)[[1]]$data
#  turplotdatahu[[k]] <- plot(conditional_effects(fits[[k]], "turbidity", dpar = "hu"), plot=FALSE, re_formula=NULL)[[1]]$data
#  
#  turplot[[k]] <- ggplot() +
#    #geom_jitter(data = datasets[[k]], aes(x = temp, y = Conc*100, color = Site)) +
#    geom_line(data = turplotdata[[k]], aes(x = effect1__, y = estimate__*100)) +
#    geom_ribbon(data = turplotdata[[k]], aes(x = effect1__, ymin = lower__*100, ymax = upper__*100), alpha = 0.3) +
#    xlab(NULL) +
#    ylab(NULL) +
#    labs(title = fullnames[k]) +
#    theme_bw() +
#    theme(legend.position = "none")
#  
#  turplothu[[k]] <- ggplot() +
#    geom_line(data = turplotdatahu[[k]], aes(x = effect1__, y = estimate__)) +
#    geom_ribbon(data = turplotdatahu[[k]], aes(x = effect1__, ymin = lower__, ymax = upper__), alpha = 0.3) +
#    xlab(NULL) +
#    ylab(NULL) +
#    #labs(title = fullnames[k]) +
#    theme_bw() +
#    theme(legend.position = "none")
#}
#
#tur_plot <- ggarrange(turplot[[1]], turplothu[[1]], turplot[[3]], turplothu[[3]], 
#                      turplot[[4]], turplothu[[4]], turplot[[6]], turplothu[[6]], 
#                      turplot[[7]], turplothu[[7]], ncol =2, nrow = 5, align = "hv")
#annotate_figure(tur_plot, left = textGrob("Expected eDNA concentration (copies/μL)", rot = 90, gp = gpar(fontsize =16)), 
#                bottom = textGrob("Secchi depth (m)", gp = gpar(fontsize =16)), 
#                right = textGrob("Predicted probability of zero eDNA detection", rot = 90, gp = gpar(fontsize =16)))

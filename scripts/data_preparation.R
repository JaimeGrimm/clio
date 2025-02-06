#Load and prepare data ----
# using readr's read_csv for better typing
data1 <- read_csv("~/Documents/GitHub/clio/data/clio.csv", skip = 42, show_col_types = FALSE)
data2 <- read_csv("~/Documents/GitHub/clio/data/clio2.csv", skip = 42, show_col_types = FALSE)
data3 <- read_csv("~/Documents/GitHub/clio/data/clio3.csv", skip = 42, show_col_types = FALSE)
envdata <- read_csv("~/Documents/GitHub/clio/data/rawsamplingdata.csv")

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
  dplyr::select(Sample, Sequence, Site, Target, Quantity, Amp) %>% 
  mutate(fullname = case_when(Target == "Te_mar" ~ "T. maritimum", 
                              Target == "Te_fin" ~ "T. finnmarkense", 
                              Target == "sch" ~ "Cand. S. salmonis",
                              Target == "sasa" ~ "Atlantic salmon", 
                              Target == "pisck_sal" ~ "P. salmonis", 
                              Target == "pa_ther" ~ "P. theridion",
                              Target == "env" ~ "ENV"))

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
# failed amplifications that still have a non-zero quantity associated with them
# +
data$Conc <- if_else(is.na(data$Quantity), 0, data$Quantity)
data <- data %>% filter(!(Amp == "No Amp" & Conc > 0)) 

# check out avg. DNA concentration against sequence number
data %>% group_by(Target, Sequence) %>%
  mutate(MeanConc= mean(Conc, na.rm = TRUE), .groups = "keep") %>%
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

#Combine with environmental data
envdata <- envdata %>% rename(Site = region, Sequence = Seq) %>% 
  filter(Site != "NA") %>% 
  dplyr::select(Site, Sequence, temp, salinity, turbidity)

envdata <- merge(data, envdata)

envdata_nc <- envdata %>%
  filter(Site != "Control") %>% # filter out control trials, for fitting
  #filter(!(Sample == "KN1922" & Target == "Te_mar")) %>% # drop the outlier sequence
  mutate(Conc = Conc / 100, PrevConc = PrevConc / 100) # scale for easier sampling

#Look at contamination in controls ----
controls <- bind_rows(data1, data2, data3) %>% 
  rename(Target = `Target Name`, Sample = `Sample Name`, Amp = `Amp Status`) %>%
  mutate(Site = case_when(
    str_detect(Sample, "CL")~"Clio",
    str_detect(Sample, "KN")~"Knight",
    TRUE~'Control')) %>% 
  mutate(fullname = case_when(Target == "Te_mar" ~ "T. maritimum", 
                              Target == "Te_fin" ~ "T. finnmarkense", 
                              Target == "sch" ~ "Cand. S. salmonis",
                              Target == "sasa" ~ "Atlantic salmon", 
                              Target == "pisck_sal" ~ "P. salmonis", 
                              Target == "pa_ther" ~ "P. theridion",
                              Target == "env" ~ "ENV"))
controls$Conc <- if_else(is.na(controls$Quantity), 0, controls$Quantity)
    
cont.only<- controls %>% filter(Site != "Clio", Site != "Knight", Conc > 0, Amp == "Amp",
                               Task != "STANDARD") 
temp <- data %>% filter(fullname == "Atlantic salmon" | fullname == "Cand. S. salmonis" | fullname == "P. salmonis" | 
                        fullname == "P. theridion")
ggplot()+
  geom_violin(data = temp, aes(y = Conc, x = fullname)) +
  geom_jitter(data = cont.only, aes(y=Conc, x = fullname, color = Sample), width = 0.1, height = 0, size = 2) +
  ylab("log(1 + eDNA contamination (copies/μL))") +
  xlab("Target") +
  scale_y_continuous(trans="log1p") +
  theme_bw(base_size = 18)+
  theme(axis.text.x = element_text(angle = 45, vjust =0.6)) +
  scale_color_discrete(labels = c("CT0884" = "Start of day 1", "CT0837" = "End of day 1", "CT0917" = "End of day 2", 
                                  "EB0108" = "Extraction blank 1", "EB0108.1" = "Extraction blank 2", 
                                 "qPCRblank" = "qPRC blank"))
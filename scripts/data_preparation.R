#Load and prepare data ----
library(tidyverse)
# using readr's read_csv for better typing
data1 <- read_csv("./data/clio.csv", skip = 42, show_col_types = FALSE) %>% 
  filter(`Sample Name` != "KN1922", `Sample Name` != "CL2037") #These two filters were outliers. Reran in data4 to confirm. Will use data4 amplifications
data2 <- read_csv("./data/clio2.csv", skip = 42, show_col_types = FALSE)
data3 <- read_csv("./data/clio3.csv", skip = 42, show_col_types = FALSE)
data4 <- read_csv("./data/clio4.csv", skip = 42, show_col_types = FALSE) 

envdata <- read_csv("~/Documents/GitHub/clio/data/rawsamplingdata.csv")

data <- bind_rows(data1, data2, data3, data4) %>%
  # simplify column names:
  rename(Target = `Target Name`, Sample = `Sample Name`, Amp = `Amp Status`, well = 'Well Position') %>%
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
    Sample != "qPCRBlank",
    Sample != "qPCRblank",
    Sample != "qPRCBlank",
    Sample != "Messedup",
    Sample != "water",
    mAMP != "ATYPICAL") %>%
  # grab sequence numbers:
  mutate(Sequence = as.numeric(substr(str_extract(Sample, "[[:digit:]]+"), 1, 2))) %>%
  # controls were run at the beginning & end of each day; increment second day's sequence number by 2:
  mutate(Sequence = if_else(Sequence < 9, Sequence, Sequence + 2)) %>%
  # grab only the columns we need:
  dplyr::select(Sample, Sequence, Site, Target, Quantity, Amp, well, CT, mAMP) %>% 
  mutate(fullname = case_when(Target == "Te_mar" ~ "T. maritimum", 
                              Target == "Te_fin" ~ "T. finnmarkense", 
                              Target == "sch" ~ "Cand. S. salmonis",
                              Target == "sasa" ~ "Atlantic salmon", 
                              Target == "pisck_sal" ~ "P. salmonis", 
                              Target == "pa_ther" ~ "P. theridion",
                              Target == "env" ~ "ENV"))
data$Sample[data$Sample =="CT0804"] <- "CT0884"

# Manually assign sequence numbers for controls:
data$Sequence[data$Sample=="CT0884"] <- 0 #beginning of day 1
data$Sequence[data$Sample=="CT0837"] <- 9 #end of day 1
data$Sequence[data$Sample=="CT0931"] <- 10 #beginning of day 2
data$Sequence[data$Sample=="CT0917"] <- 27 #end of day 2
# -
data <- data %>% filter(!grepl("CT", Sample))

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
foo <- data %>% filter((mAMP == "NOAMP" & Conc > 0)) 

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
library(PNWColors)
pal <- pnw_palette("Starfish")

controls <- bind_rows(data1, data2, data3, data4) %>% 
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
                              Target == "env" ~ "ENV")) %>% 
  filter(Target != "ascv", Target != "ctv-2", mAMP != "ATYPICAL")
controls$Conc <- if_else(is.na(controls$Quantity), 0, controls$Quantity)
controls$Sample[controls$Sample=="qPCRblank"] <- "qPCRBlank"
controls$Sample[controls$Sample=="qPRCBlank"] <- "qPCRBlank"
controls$Sample[controls$Sample=="water"] <- "qPCRBlank"
controls$Sample[controls$Sample =="CT0804"] <- "CT0884"
    
cont.only<- controls %>% filter(Site != "Clio", Site != "Knight", #Amp == "Amp",
                               Task != "STANDARD", Sample != "Messedup") 
temp <- data %>% filter(fullname == "Atlantic salmon" | fullname == "Cand. S. salmonis" | fullname == "P. salmonis" | 
                        fullname == "P. theridion")
ggplot()+
  geom_jitter(data = data, aes(y = Conc, x = fullname), color = "grey70", alpha=0.3, width = 0.3, size = 2.5) +
  geom_point(data = cont.only, aes(y=Conc, x = fullname, color = Sample), position = position_jitter(seed = 2, width = 0.2, height = 0), size = 3.5) +
  ylab("eDNA concentration (copies/μL)") +
  xlab("Target") +
  scale_y_continuous(trans="sqrt") +
  theme_bw(base_size = 18)+
  theme(axis.text.x = element_text(angle = 45, vjust =0.6),
        axis.line.y.right = element_blank()) +
  scale_color_discrete(labels = c("CT0884" = "Start-of-day 1", "CT0837" = "End-of-day 1", "CT0931" = "Start-of-day 2",
                                  "CT0917" = "End-of-day 2", "EB0108" = "Extraction blank 1", "EB0108.1" = "Extraction blank 2", 
                                 "qPCRBlank" = "qPRC blank")) +
  scale_y_break(breaks = c(2000, 30000), scales = 0.35, expand = TRUE)
  #scale_color_manual(values = pal)


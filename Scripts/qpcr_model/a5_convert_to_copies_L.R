## NGN Running stan model for QPCR data
# Author: Eily Allan but really Ryan Kelly
# Person running: Eily
# Last modified: 10/19/22 by Eily
# Date of run: 03/21/2023 by Eily 

# Overview 
# This script 

# Inputs: 
# 1) 

# Outputs: 
# 1) 

####################################################################
# Set up
####################################################################

# Load packages
library(tidyverse)
library(rstan)
library(here)
library(readxl)
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

# Call functions for Stan model
source(here("Scripts","functions", "calibrate_qPCR.R")) 

# Find metadata and all data for cutthroat and coho qpcr
#### Metadata file is for both cutthroat and coho 
# qPCRmeta <- here("Input","qpcr","qPCR_samples_ALL.xlsx")

wsdotnames <- read.csv(here("Input","qpcr","wsdot_names_wrangle.csv"))
# redometa <- read_xlsx(qPCRmeta, sheet = "info_samples")
# 
# qPCRmeta2 <- redometa %>% 
#   rename(new_name = sample) %>% 
#   left_join(wsdotnames %>% filter(str_detect(plate,"CUT")), by="new_name") %>% 
#   select(-plate) %>% 
#   mutate(qPCR_name = case_when(is.na(qPCR_name) ~ new_name, 
#                                !is.na(qPCR_name) ~ qPCR_name)) %>% 
#   select(-new_name) %>% 
#   rename(sample=qPCR_name)
# 
# write.csv(qPCRmeta2, here("Input","qpcr","qPCR_samples_ALL_NEWNAMES.csv"), row.names=FALSE)

qPCRmeta <- here("Input","qpcr","qPCR_samples_ALL_NEWNAMES.csv")

adj_vol_filtered <- read.csv(here("Output","qpcr","adj_vol_filtered.csv"))

#### Cutthroat and coho file directory 
cut_files <- here("Input","qpcr","Results","CUT")
coho_files <- here("Input","qpcr","Results","COHO")
n_cut <- length(list.files(path = cut_files, pattern = "*data", recursive = T, full.names = T))
n_coho <- length(list.files(path = coho_files, pattern = "*data", recursive = T, full.names = T))


####################################################################
# Use functions to read in files and format for Stan model
####################################################################

## cutthroat 
cutcleandata <- list()
for (i in 1:n_cut) {
  plate_num <- i
  cutcleandata[[i]] <- input_qPCR_data(cut_files, qPCRmeta, plate_num)
}
cut_data_for_stan <- do.call(rbind, cutcleandata)  
cut_data_for_stan <- cut_data_for_stan %>% 
  select(-DateStamp)

# remove plate 98 standard 7 (because looks wrong) 
cut_data_for_stan <- cut_data_for_stan %>% 
  filter(! (Plate==98 & Sample == "St7"))

cut_data_for_stan2 <- cut_data_for_stan %>% 
  rename(qPCR_name = Sample) %>%
  left_join(wsdotnames %>% filter(str_detect(plate,"CUT")) , by="qPCR_name") %>%
  select(-plate) %>%
  mutate(new_name = case_when(is.na(new_name) ~ qPCR_name,
                             !is.na(new_name) ~ new_name)) %>%
  select(-qPCR_name) %>%
  rename(Sample=new_name) %>%
  mutate(dilution = case_when(dilution == NA ~ 1,
                              dilution == "NA" ~ 1,
                             dilution != "NA" ~ dilution)) %>% 
  mutate(dilution = as.numeric(dilution))
  

check <- cut_data_for_stan2 %>% 
  separate(Sample, into=c("time","creek","stn","bio")) %>% 
  unite(creekstnbio, c(creek,stn,bio)) 

check2 <- check %>% 
  group_by(creekstnbio) %>% 
  summarize(n=n())

check3 <- cut_data_for_stan2 %>% 
  filter(!str_detect(Sample, "St")) %>% 
  group_by(Sample) %>% 
  mutate(sdCt = sd(Ct, na.rm=TRUE))

check4 <- check3 %>% 
  select(c(dilution,Sample,sdCt)) %>% 
  distinct()

keeplowsds <- check3 %>% 
  filter(sdCt < 1.5) %>% 
  distinct()

keepstds <- cut_data_for_stan2 %>% 
  filter(str_detect(Sample, "St"))

cut_data_for_stan2 <- cut_data_for_stan2 %>% 
  filter(Sample %in% keeplowsds$Sample) %>% 
  mutate(Sample = case_when(Sample == "822.5Sqm.Up.3" ~ "0822.5Sqm.Up.3",
                            Sample == "822.5Sqm.Up.2" ~ "0822.5Sqm.Up.2",
                            Sample == "822.4Pad.Up5.3" ~ "0822.4Pad.Up5.3",
                            Sample == "822.4Pad.Up5.1" ~ "0822.4Pad.Up5.1", 
                            Sample == "722.5Sqm.Up.3" ~ "0722.5Sqm.Up.3",
                            Sample == "722.1Prt.Up.3" ~ "0722.1Prt.Up.3",
                            Sample == "722.1Prt.Up.2" ~ "0722.1Prt.Up.2",
                           TRUE ~ Sample))

cut_data_for_stan3 <- rbind(cut_data_for_stan2, keepstds)

# write.csv(cut_data_for_stan3, here("Output","qpcr","cut_data_for_stan.csv"), row.names=FALSE)

####################################################################
# Run Stan model for cutthroat 
####################################################################

#cut_data_for_stan <- read_csv(here("Output","qpcr","cut_data_for_stan.csv"))

cut_qMod_out <- run_qPCR_model(here("Output","qpcr","cut_data_for_stan.csv"),
                               here("Scripts", "functions", "qPCR_calibration_enchilada.stan"))

#save whole model output
write_rds(cut_qMod_out, here("Output","qpcr","20230505_cut_qPCR_modelFit.RDS"))

# save just concentrations -- units here are copies / L water (not yet flow corrected) 
cut_modeled_conc <- cut_qMod_out$results_qPCR
write_rds(cut_modeled_conc, here("Output","qpcr","20230505_cut_modeled_conc.RDS"))

####################################################################
# Check output to make sure it looks reasonable  
####################################################################

#cut_modeled_conc <- readRDS(here("Output","qpcr","20230426_cut_modeled_conc.RDS"))

checkcut <- cut_modeled_conc %>% 
  unite(Sample, c(time,creek,station,biorep), sep =".") %>% 
  select(-c(Plate, dilution, Adj_Vol)) %>% 
  left_join(cut_data_for_stan3, by = "Sample") %>% 
  separate(Sample, into=c("time","creek","station","biorep"))

checkcut <- cut_data_for_stan3 %>% 
  left_join( (cut_modeled_conc %>% 
                unite(Sample, c(time,creek,station,biorep), sep =".") %>% 
                select(-c(Plate, dilution, Adj_Vol))), by="Sample" ) %>% 
  separate(Sample, into=c("time","creek","station","biorep"))

checkcut %>% 
  filter(Plate < 10) %>% 
  filter(!str_detect(time, "St")) %>% 
ggplot(aes(y=Ct, x=mean_concentration_est, color=Plate)) +
  geom_point() +
  scale_x_log10() +
  facet_wrap(~dilution~Plate) +
  scale_color_continuous(type = "viridis") 

checkcut %>% 
  filter(Plate > 10 & Plate < 90) %>% 
  filter(!str_detect(time, "St")) %>% 
  ggplot(aes(y=Ct, x=mean_concentration_est, color=Plate)) +
  geom_point() +
  scale_x_log10() +
  facet_wrap(~dilution~Plate) +
  scale_color_continuous(type = "viridis") 

checkcut %>% 
  filter(Plate > 90) %>% 
  filter(!str_detect(time, "St")) %>% 
  ggplot(aes(y=Ct, x=mean_concentration_est, color=Plate)) +
  geom_point() +
  scale_x_log10() +
  facet_grid(~dilution~Plate) +
  scale_color_continuous(type = "viridis") 

stds_only <- checkcut %>% 
  filter(str_detect(time, "St")) %>% 
  select(c(Plate, Ct, time)) %>% 
  rename(Sample = time) %>% 
  mutate(Quantity = case_when(Sample == "St1" ~ 100000,
                              Sample == "St2" ~ 10000,
                              Sample == "St3" ~ 1000,
                              Sample == "St4" ~ 100,
                              Sample == "St5" ~ 10,
                              Sample == "St6" ~ 5,
                              Sample == "St7" ~ 3,
                              Sample == "St8" ~ 1))

ggplot(stds_only, aes(y=Ct, x=Quantity)) +
  geom_point() +
  scale_x_log10() +
  facet_wrap(~Plate) +
  scale_color_continuous(type = "viridis") 

fitted_models <- stds_only %>% group_by(Plate) %>% do(model = lm(Ct ~ log10(Quantity), data = .))

####################################################################
# Plot cutthroat environmental samples 
####################################################################

cut_modeled_conc %>% 
  ggplot(aes(x = time, y = log(mean_concentration_est), col = station)) +
  geom_point() +
  geom_segment(aes(x = time, xend = time, y = log(ci25_concentration_est), yend = log(ci75_concentration_est))) +
  facet_grid(rows=vars(creek)) +
  labs(y= "Log copies/L water", title = "Cutthroat Trout") + 
  theme_bw()



simplecompare <- cut_modeled_conc %>% 
  unite(Sample, c("time","creek","station","biorep"), sep=".") %>% 
  left_join(cut_data_for_stan2 %>% select(c(Sample,Ct)), by = "Sample")


simplecompare %>% 
  filter(Plate<8) %>% 
  #filter(dilution==1) %>% 
  #filter(Adj_Vol==2) %>% 
  filter(mean_concentration_est > 1) %>% 
ggplot(aes(y=log(mean_concentration_est), x=Ct, color=Plate)) +
  geom_point() +
  facet_wrap(~dilution~Plate)


padup5compare1 <- cut_data_for_stan %>% 
  select(c(Sample,Ct)) %>% 
  filter(str_detect(Sample, "Up5")) %>% 
  rename(qPCR_name = Sample) %>% 
  left_join(wsdotnames, by="qPCR_name") %>%
  mutate(new_name = case_when(is.na(new_name) ~ qPCR_name,
                              !is.na(new_name) ~ new_name)) %>%
  select(-c(plate,qPCR_name)) %>%
  rename(Sample=new_name) %>% 
  drop_na()

padup5compare2 <- cut_data_for_stan3 %>% 
  select(c(Sample,Ct)) %>% 
  filter(str_detect(Sample, "Up5"))

padup5_highsds <- check3 %>% 
  filter(str_detect(Sample,"Up5"))
  

padup5compare1 %>% 
  #filter(mean_concentration_est > 1) %>% 
  ggplot(aes(y=log(mean_concentration_est), x=Ct, color=Plate)) +
  geom_point() +
  facet_wrap(~dilution~Plate)

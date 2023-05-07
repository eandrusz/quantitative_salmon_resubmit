## NGN qPCR Inhibition Checking
# Author: Eily Allan 
# Person running: Eily
# Last modified: 10/19/22 by Eily
# Date of run: 03/21/23 by Eily 

# Overview 
# This script takes 

# Inputs: 
# 1) 

# Outputs: 
# 1)

####################################################################
# Set up
####################################################################

# Load packages
library(tidyverse)
library(here)
library(readxl)
source(here("Scripts","functions","inhibition_functions.R"))

## Point to files for cutthroat and inhibition (IPC)
qPCRfiles <- list.files(path = here("Input","qpcr","Results","CUT"), pattern = "*data", recursive = T, full.names = T)
qPCRmeta <- read.csv(file=here("Input","qpcr","qPCR_samples_ALL_NEWNAMES.csv"))

####################################################################
# Run inhibition testing functions for all plates
####################################################################

redoinhib <- list()
thresh <- vector()
for (i in 1:length(qPCRfiles)) {
  plate_num <- i
  plate <- simple_IPC_plate(plate_num)
  thresh[i] <- determine_thresh(plate)
  inhibited <- inhib_samp(plate, thresh[i]) 
  redoinhib[[i]] <- still_inhib_samp(plate_num,inhibited)
}

# plate_num = 18 --> plate 23 does not have IPC (only redoing standard deviations)
# plate_num = 24 --> plate 25 does not have IPC (only redoing standard deviations)
# plates 30 and 31 --> IPC was run before so these are all not inhibited

redoinhib <- do.call(rbind, redoinhib)  

# write.csv(redoinhib, here("Output","qpcr","cut_fourth_inhibition.csv"), row.names=FALSE)
# This *SHOULD* have nothing left 

####################################################################
# Run high std dev check for all plates
####################################################################

badsd <- list()
for (i in 1:length(qPCRfiles)) {
  plate_num <- i
  plate <- simple_IPC_plate(plate_num)
  thresh[i] <- determine_thresh(plate)
  badsd[[i]] <- bad_sd_samps(plate_num,plate, thresh[i])
}
badsd <- do.call(rbind, badsd)  


### Manually check mayas because they won't go through because there is no IPC data 
plate98 <- read_excel(qPCRfiles[31], sheet = "Results", skip=6)
colnames(plate98)[7] <- "Ct"  
colnames(plate98)[2] <- "Sample"
colnames(plate98)[3] <- "Assay"
plate98_sd <- plate98 %>% 
  #filter(!str_detect(Sample,"St")) %>% 
  filter(Ct != "Undetermined") %>% 
  select(c(Sample, Ct)) %>% 
  mutate(Ct = as.numeric(Ct)) %>% 
  group_by(Sample) %>% 
  mutate(meanCt = mean(Ct)) %>% 
  mutate(sdCt = sd(Ct)) %>% 
  select(c(Sample,meanCt,sdCt)) %>% 
  distinct() 

plate99 <- read_excel(qPCRfiles[32], sheet = "Results", skip=6)
colnames(plate99)[7] <- "Ct"  
colnames(plate99)[2] <- "Sample"
colnames(plate99)[3] <- "Assay"
plate99_sd <- plate99 %>% 
  #filter(!str_detect(Sample,"St")) %>% 
  filter(Ct != "Undetermined") %>% 
  select(c(Sample, Ct)) %>% 
  mutate(Ct = as.numeric(Ct)) %>% 
  group_by(Sample) %>% 
  mutate(meanCt = mean(Ct)) %>% 
  mutate(sdCt = sd(Ct)) %>% 
  select(c(Sample,meanCt,sdCt)) %>% 
  distinct()

# write.csv(redoinhib, here("Output","qpcr","cut_fourth_inhibition.csv"), row.names=FALSE)
# This *SHOULD* have nothing left 




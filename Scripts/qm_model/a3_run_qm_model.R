## NGN Run QM Model
# Author: Eily Allan - but really all from Ryan Kelly and Ole Shelton 
# Person running: Eily
# Last modified: 10/17/22 by Eily
# Date of run: 03/21/2023 by Eily 

# Overview 
# 

# Inputs: 
# 1)  

# Outputs: 
# 1) 

# Run separately for each creek. 
# Remove O. nerka for Portage - otherwise it won't converge. 

####################################################################
# Set up
####################################################################

## Load packages
library(tidyverse)
library(rstan)
library(MCMCpack) #for rdirichelet function
library(here)
library(gridExtra)
library(unikn)
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

source(here("Scripts","functions", "calibrate_metabarcoding.R")) 

# read in taxa table 
taxa_table <- read.csv(here("Output","metabarcoding", "taxa_table.csv"))
mock <- readRDS(here("Output","metabarcoding","20230321_mockdatatocalibrate.RDS"))
fixnames <- read.csv("/Users/elizabethandruszkiewicz/GoogleDrive/UW/GitHub/wsdot_salmonids/samples_rename.csv")

wsdot <- taxa_table %>% 
  filter(! str_detect(Sample_name, "MC")) %>%
  filter(!str_detect(Sample_name, "MiFish")) %>% 
  left_join(fixnames, by="Sample_name") %>% 
  dplyr::select(-Sample_name) %>%
  separate(New_name, into=c("creek", "station","time","bio","tech")) %>% 
  mutate(creek=case_when(creek == "Pad" ~ "4Pad",
                         creek == "Prt" ~ "1Prt",
                         creek == "Sqm" ~ "5Sqm")) %>% 
  mutate(marker="MiFish") %>% 
  unite(Sample_name, c("marker", "time", "creek", "station", "bio", "tech"), sep=".") %>% 
  filter(totalReads > 0) %>% 
  dplyr::select(c(Sample_name, species, totalReads))

ngn <- taxa_table %>% 
  filter(! str_detect(Sample_name, "MC")) %>%
  filter(str_detect(Sample_name, "MiFish")) 

enviro <- rbind(ngn,wsdot)

enviro <- enviro %>% 
  separate(Sample_name, into=c("marker", "time", "creek", "station", "biol", "tech"), remove=TRUE) %>% 
  #filter(station != "Up5") %>% 
  # mutate(station = case_when(station == "Up11" ~ "Up", 
  #                             TRUE ~ station)) %>% 
  dplyr::rename(Nreads = totalReads) %>% 
  dplyr::select(-marker) %>% 
  filter(creek == "1Prt") #### HERE CHANGE TO RUN EACH CREEK INDIVIDUALLY

# only focus on four salmonids
mock <- mock %>%
  #filter(species %in% c("Oncorhynchus clarkii","Oncorhynchus kisutch","Oncorhynchus mykiss","Oncorhynchus nerka"))
  filter(species %in% c("Oncorhynchus clarkii","Oncorhynchus kisutch","Oncorhynchus mykiss"))

# Prepare for stan model 
qmdata <- format_metabarcoding_data(enviro, mock)
stan_metabarcoding_data <- makeDesign(qmdata, N_pcr_cycles = 43)

# Set color palette so don't change when plot
o_i_colors <- c(#rgb(230, 159,   0, maxColorValue = 255),  # orange
                rgb( 86, 180, 233, maxColorValue = 255),  # skyblue
                rgb(  0, 158, 115, maxColorValue = 255),  # green
                #rgb(240, 228,  66, maxColorValue = 255),  # yellow
                rgb(  0, 114, 178, maxColorValue = 255))  # blue
                #rgb(204, 121, 167, maxColorValue = 255)   # purple
)
#o_i_colors = scales::hue_pal()(length(unique(mock$species)))
pal_okabe_ito <- newpal(col = o_i_colors,
                        names = unique(mock$species))

##################################################################
# Run Bayesian QM model
####################################################################

bayes_out <- QM_bayes(here("Scripts", "functions", "quant_metabar_rosetta_noSampleEta.stan"), stan_metabarcoding_data)
#write_rds(bayes_out, here("Output","metabarcoding", "qm_model", "20230406_bayesout_3spp_1Prt.RDS"))

summaryout <- summary(bayes_out$Bayes_modelfit)$summary
#write.csv(summaryout, here("Output","metabarcoding", "qm_model", "20230406_bayesout_3spp_summary_1Prt.csv"))


bayes_out$Bayes_estimates %>%
  rownames_to_column("sample") %>%
  pivot_longer(-sample, names_to = "species") %>%
  separate(col = sample, into = c("time", "creek", "station", "biol"), remove = FALSE) %>% 
  mutate(station = case_when(station == 1 ~ "Down",
                             station == 2 ~ "Up")) %>% 
  separate(time, into = c("month","year"), sep = 2, remove=FALSE) %>% 
  unite(newtime, c(year,month), sep="-", remove=FALSE) %>% 
  filter(value > 0.001) %>%
  mutate(creek = case_when(creek == "1Prt" ~ "Portage",
                           creek == "2Brn" ~ "Barnes",
                           creek == "3Chk" ~ "Chuckanut",
                           creek == "4Pad" ~ "Padden",
                           creek == "5Sqm" ~ "Squalicum",
                           TRUE ~ creek))  %>% 
  mutate(facetorder = factor(creek, levels=c('Padden','Portage','Chuckanut','Squalicum', 'Barnes'))) %>% 
  unite(creekstntime, c(creek,station,time), remove=FALSE) %>% 
  group_by(creekstntime) %>% 
  mutate(sumprop = sum(value)) %>% 
  mutate(avgprop = value/sumprop) %>%
  ggplot(aes(x = newtime, fill = species, y = avgprop)) +
  geom_col() +
  #theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))  +
  facet_grid(~facetorder ~station) +
  scale_fill_manual(values = pal_okabe_ito) + 
  scale_color_manual(values = pal_okabe_ito) +
  ylab("Proportion of DNA after QM correction") +
  labs(x="Date (YY-MM)", y="Proportion of DNA after QM correction", fill="Species", color="") %>% 
  scale_x_discrete(guide = guide_axis(angle = -45)) +
  theme_bw()

#ggsave(here("Output","SupplementalFigures","20221121_proportions_after_qm.png"))

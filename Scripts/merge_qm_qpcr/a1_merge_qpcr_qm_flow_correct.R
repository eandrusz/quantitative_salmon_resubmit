## NGN Merging QM results with QPCR results -- FLOW CORRECTED 
# Author: Eily Allan and Ryan Kelly
# Person running: Eily
# Last modified: 10/20/22 by Eily
# Date of run: 03/23/2023 by Eily 

# Overview 
# This script 

# Inputs: 
# 1) 

# Outputs: 
# 1) 

####################################################################
# Set up
####################################################################

## Load packages
library(tidyverse)
library(here)
library(gtools)

source(here("Scripts", "functions","sample_posteriors.R"))

## Pull posteriors of Bayesian model for QM
bayes_out_prt <- readRDS(here("Output","metabarcoding", "qm_model", "20230406_bayesout_3spp_1Prt.RDS"))
bayes_out_brn <- readRDS(here("Output","metabarcoding", "qm_model", "20230406_bayesout_4spp_2Brn.RDS"))
bayes_out_chk <- readRDS(here("Output","metabarcoding", "qm_model", "20230406_bayesout_4spp_3Chk.RDS"))
bayes_out_pad <- readRDS(here("Output","metabarcoding", "qm_model", "20230406_bayesout_4spp_4Pad.RDS"))
bayes_out_sqm <- readRDS(here("Output","metabarcoding", "qm_model", "20230406_bayesout_4spp_5Sqm.RDS"))

## Pull posteriors of Bayesian model for qPCR
qMod_out <- readRDS(here("Output","qpcr","20230505_cut_qPCR_modelFit.RDS"))
cut_names <- readRDS(here("Output","qpcr","20230505_cut_modeled_conc.RDS"))

## Flow data to correct
flow <- read.csv(here("Output","qpcr","flowrates_touse_extended.csv"))

####################################################################
# Deal with QM model output
####################################################################

b_prt <- qm_pull_posteriors(bayes_out_prt)
b_brn <- qm_pull_posteriors(bayes_out_brn)
b_chk <- qm_pull_posteriors(bayes_out_chk)
b_pad <- qm_pull_posteriors(bayes_out_pad)
b_sqm <- qm_pull_posteriors(bayes_out_sqm)

b <- rbind(b_prt, b_brn, b_chk, b_pad, b_sqm)

####################################################################
# Deal with qPCR model output 
####################################################################

cut_names <- cut_names %>%
  # mutate(creek = case_when(creek == "4Pad" & station == "Dn" ~ "4Pad11",
  #                          creek == "4Pad" & station == "Up11" ~ "4Pad11",
  #                          TRUE ~ creek)) %>% 
  # mutate(creek = case_when(creek == "4Pad" & station == "Up5" ~ "4Pad5",
  #                          TRUE ~ creek)) %>% 
  # mutate(station = case_when(creek == "4Pad11" & station == "Up11" ~ "Up",
  #                            TRUE ~ station)) %>% 
  # mutate(station = case_when(creek == "4Pad5" & station == "Up5" ~ "Up",
  #                            TRUE ~ station)) %>% 
  mutate(station = case_when(station == "Dn" ~ 1,
                             station == "Up" ~ 2,
                             station == "Up11" ~ 2,
                             station == "Up5" ~ 3)) %>%
  unite(Sample, c("time","creek","station","biorep"), sep="_") %>% 
  dplyr::select(c(Sample,dilution,Adj_Vol))

qgrabthese <- grep("envir_concentration", names(qMod_out$qMod@sim$samples[[1]]))
qpostSamples <- qMod_out$qMod@sim$samples[[1]][qgrabthese]


####################################################################
# Deal with flow
####################################################################

flowtojoin <- flow %>% 
  pivot_longer(!daysampletoavg, names_to="creek", values_to="flow_m3s") %>% 
  mutate(creek = case_when(creek == "pad_flowm3s" ~ "4Pad",
                           creek == "flow_chk_by_padden" ~ "3Chk",
                           creek == "flow_sqm_by_padden" ~ "5Sqm",
                           creek == "flow_prt_by_padden" ~ "1Prt",
                           creek == "flow_brn_by_padden" ~ "2Brn",
                           TRUE ~ creek)) %>% 
  separate(daysampletoavg, c("year","month","day")) %>% 
  mutate(year= str_sub(year, 3, 4))  %>%
  mutate(month = case_when(month == "3" ~ "03",
                           month == "4" ~ "04",
                           month == "5" ~ "05",
                           month == "6" ~ "06",
                           month == "7" ~ "07",
                           month == "8" ~ "08",
                           month == "9" ~ "09",
                           month == "1" ~ "01",
                           month == "2" ~ "02",
                           TRUE ~ month)) %>% 
  unite(time, c(month, year), sep="") %>% 
  unite(timecreek, c(time, creek)) %>%
  mutate(flow_m3s = ifelse(is.na(flow_m3s), 0.00141584235, flow_m3s))

qa <- as.data.frame(qpostSamples) 
qb <- t(qa) %>% 
  as_tibble() %>% 
  add_column(cut_names) %>% 
  #add_column(coho_names) %>% 
  group_by(Sample, dilution, Adj_Vol) %>% 
  nest() %>% 
  rename(cutqpcr = data) %>% 
  rename(bottle = Sample) %>% 
  # for flow correction
  separate(bottle, into=c("time","creek","station","bio"), remove=FALSE) %>% 
  mutate(time = case_when(time == "822" ~ "0822",
                           time == "722" ~ "0722",
                           TRUE ~ time)) %>% 
  unite(timecreek, c(time,creek), remove=FALSE) %>%  
  unite(bottle, c("time","creek","station","bio"), remove=FALSE) %>% 
  left_join(flowtojoin, by="timecreek") %>% 
  select(-c(timecreek,station,bio, bottle)) %>% 
  group_by(bottle) %>% 
  arrange(desc(dilution)) %>%
  slice(1)

####################################################################
# Check for high and low percentages and remove, and low qpcr
####################################################################

mergedb <- b %>% 
  left_join(qb, by="bottle") %>% 
  drop_na() %>% 
  mutate(cut_dnacopy_uL = map2(.x = cutqpcr, .y = dilution, .f = function(x, y) y*(10^x))) %>% 
  mutate(cut_dnacopy_L = map2(.x = cut_dnacopy_uL, .y = Adj_Vol, .f = function(x,y) x*(50/y))) %>% 
  ungroup()

mergedb <- mergedb %>%
  group_by(bottle, species) %>% 
  mutate(meanpropreads = mean(unlist(data))) %>%
  mutate(meancutqpcr = mean(unlist(cut_dnacopy_uL))) %>% 
  mutate(meancutdna = mean(unlist(cut_dnacopy_L))) %>% 
  separate(bottle, into=c("time","creek","station", "bio"), remove=FALSE)

mergedb %>%
  #filter(meanpropreads < 0.05) %>% 
ggplot(aes(x=meanpropreads*100)) +
  geom_histogram() +
  facet_wrap(~species) +
  labs(x="Percentage of Sequencing Reads", y="Frequency", title="Percentage of Reads < 5%") +
  theme_bw()
ggsave(here("Output","SupplementalFigures","histograms_percent_reads_by_species.png"))
#ggsave(here("Output","SupplementalFigures","histograms_percent_reads_by_species_less5.png"))

mergedb_1_99pct <- mergedb %>% 
  filter(meanpropreads > 0.01 & meanpropreads < 0.99)

mergedb %>%
  filter(species=="Oncorhynchus clarkii") %>% 
  ggplot(aes(x=log10(meancutqpcr))) +
  geom_histogram() +
  #facet_wrap(~creek) +
  labs(x="Log10 Cutthroat trout concentration (copies/uL)", y="Frequency") +
  theme_bw()
ggsave(here("Output","SupplementalFigures","histogram_cut_conc_copy_ul.png"))

mergedb %>%
  filter(species=="Oncorhynchus clarkii") %>%
  filter(meancutqpcr < 50) %>% 
  ggplot(aes(x=meancutqpcr)) +
  geom_histogram() +
  #facet_wrap(~creek) +
  labs(x="Log10 Cutthroat trout concentration (copies/uL)", y="Frequency") +
  theme_bw()
ggsave(here("Output","SupplementalFigures","histogram_cut_conc_copy_ul_50.png"))

mergedb_1_99pct_2copies <- mergedb_1_99pct %>% 
  filter(meancutqpcr > 2)

####################################################################
# Actually merge
####################################################################

#calculate mean total DNA, only using cutthroat
meantotdna <- mergedb_1_99pct_2copies %>% 
  filter(species == "Oncorhynchus clarkii") %>% 
  mutate(meantotdna = meancutdna/meanpropreads) %>% 
  dplyr::select(bottle, meantotdna)

mergedc <- mergedb_1_99pct_2copies
mergedc$meantotdna <- meantotdna$meantotdna[match(mergedc$bottle, meantotdna$bottle)]

#check; are all values for total dna the same?
mergedc %>% 
  dplyr::select(bottle, species, meanpropreads, meancutdna, meantotdna) %>% 
  filter(bottle == "0321_1Prt_2_1")

mergedc %>%
  ggplot(aes(x=log10(meantotdna))) +
  geom_histogram() +
  labs(x="Log10 total DNA (copies/L)", y="Frequency") +
  theme_bw()

mergedc <- mergedc %>% 
  drop_na() %>% 
  ungroup() %>% 
  mutate(meandnaconc = meanpropreads*meantotdna)

#check; do concentration values make sense?
mergedc %>% 
  dplyr::select(bottle, species, meanpropreads, meancutdna, meantotdna, meandnaconc) %>% 
  filter(bottle == "0421_1Prt_1_1")


quants_to_plot <- mergedc %>% 
  mutate(station = case_when(station == 1 ~ "Down",
                             station == 2 ~ "Up",
                             station == 3 ~ "Up5")) %>%
  separate(time, into = c("month","year"), sep = 2) %>% 
  unite(newtime, c(year,month), sep="-") %>% 
  mutate(creek = case_when(creek == "1Prt" ~ "Portage",
                           creek == "2Brn" ~ "Barnes",
                           creek == "3Chk" ~ "Chuckanut",
                           creek == "4Pad" ~ "Padden",
                           creek == "5Sqm" ~ "Squalicum",
                           TRUE ~ creek)) %>% 
  # add on flow multiplier
  mutate(meandnaconcflow = meandnaconc*flow_m3s*1000) # convert 1000 L = m3

ggplot(quants_to_plot, aes(x=newtime, y=log10(meandnaconcflow), color=station)) +
  geom_point() +
  facet_grid(~creek~species) +
  # scale_color_manual(values = pal_okabe_ito) +
  guides(color= 'none') + 
  labs(x="Date (YY-MM)", y= "Log10 copies/s") + 
  theme_bw() + 
  scale_x_discrete(guide = guide_axis(angle = -45))

ggplot(quants_to_plot, aes(x=newtime, y=log10(meantotdna))) +
  geom_point() +
  facet_grid(~creek ~station) +
  labs(x="Date (YY-MM)", y= "Log10 copies total salmonid DNA/s") + 
  theme_bw() + 
  scale_x_discrete(guide = guide_axis(angle = -45))

#ggsave(here("Output","SupplementalFigures","totaldna_ts_after_qm.png"))


simple <- quants_to_plot %>% 
  select(c(bottle, newtime, creek, station, bio, species, meanpropreads, meancutdna, meantotdna, meandnaconc, meandnaconcflow))

#write_rds(simple, here("Output","20230505_abundance_flowcorrected.RDS"))


## NGN Visualizing qPCR results before / after Stan model
# Author: Eily Allan 
# Person running: Eily
# Last modified: 10/19/22 by Eily
# Date of run: 10/19/22 by Eily 

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
library(here)
library(tidyverse)

cut_data_for_stan <- read_csv(here("Output","qpcr","cut_data_for_stan.csv"))
cut_modeled_conc <- readRDS(here("Output","qpcr","20230505_cut_modeled_conc.RDS"))
monthlyflow <- read.csv(here("Output","qpcr","monthly_flow.csv"))
closestflow <- read.csv(here("Output","qpcr","closest_flow.csv"))
flowtouse <- read.csv(here("Output","qpcr","flowrates_touse_extended.csv"))

####################################################################
# qPCR: standard curves before Stan model
####################################################################

cutplot <- cut_data_for_stan %>% 
  filter(str_detect(Sample, "St")) %>% 
  mutate(Quantity = 0) %>% 
  mutate(Quantity = case_when(Sample == "St1" ~ 100000,
                                 Sample == "St2" ~ 10000,
                                 Sample == "St3" ~ 1000,
                                 Sample == "St4" ~ 100,
                                 Sample == "St5" ~ 10,
                                 Sample == "St6" ~ 5,
                                 Sample == "St7" ~ 3,
                                 Sample == "St8" ~ 1,
                              TRUE ~ Quantity)) 
  #filter(Plate == 2) %>% 
  #filter(z == 1) %>% 
  ggplot(cutplot, aes(x= log(Quantity), y = Ct, color = as.factor(Plate))) +
  geom_point() +
  geom_smooth(method = "lm") +
  theme_bw()


####################################################################
# qPCR: cutthroat modeled concentration -- up/down on same plot
####################################################################
  
cut_modeled_conc_plot <- cut_modeled_conc %>% 
  mutate(time = case_when(time == "822" ~ "0822",
                          time == "722" ~ "0722",
                             TRUE ~ time)) %>% 
  separate(time, into = c("month","year"), sep = 2) %>% 
  unite(newtime, c(year,month), sep="-") %>% 
  mutate(creek = case_when(creek == "1Prt" ~ "Portage",
                           creek == "2Brn" ~ "Barnes",
                           creek == "3Chk" ~ "Chuckanut",
                           creek == "4Pad" ~ "Padden",
                           creek == "5Sqm" ~ "Squalicum",
                           TRUE ~ creek)) %>% 
  mutate(station = case_when(station == "Dn" ~ "Down",
                             TRUE ~ station)) %>% 
  mutate(facetorder = factor(creek, levels=c('Padden','Portage','Chuckanut','Squalicum', 'Barnes')))

brnonly <- cut_modeled_conc_plot %>% 
    filter(creek=="Barnes")
chkonly <- cut_modeled_conc_plot %>% 
  filter(creek=="Chuckanut") %>% 
  filter(station=="Down")

ggplot(cut_modeled_conc_plot, aes(x = newtime, y = log(mean_concentration_est), color=station)) +
  geom_point() +
  stat_smooth(method = "loess", span=0.5, aes(group=station)) +
  geom_segment(aes(x = newtime, xend = newtime, y = log(ci25_concentration_est), yend = log(ci75_concentration_est))) +
  geom_rect(data=brnonly, aes(xmin = 12, xmax = 20, ymin = -Inf, ymax = Inf), alpha = .01, color=NA, fill = "grey") +
  geom_rect(data=chkonly, aes(xmin = 12, xmax = 20, ymin = -Inf, ymax = Inf), alpha = .01, color=NA, fill = "grey") +
  facet_grid(rows=vars(facetorder)) +
  labs(x="Date (YY-MM)", y= "Log copies/L water", title = "Cutthroat Trout", color="Station") + 
  theme_bw() + 
  scale_x_discrete(guide = guide_axis(angle = -45)) +
  scale_color_manual(values=c("darkgrey", "black", "black", "blue"))

# ggsave(here("Output","SupplementalFigures", "modeled_cut_qpcr_updown-line.png"))


####################################################################
# qPCR: cutthroat modeled concentration MULTIPLIED BY FLOW
####################################################################

## for now, use padden for all of them
  flowmerge <- flowtouse %>% 
    pivot_longer(!daysampletoavg, names_to="creek", values_to="flow_m3s") %>% 
    mutate(creek = case_when(creek == "pad_flowm3s" ~ "4Pad",
                             creek == "flow_chk_by_padden" ~ "3Chk",
                             creek == "flow_sqm_by_padden" ~ "5Sqm",
                             creek == "flow_prt_by_padden" ~ "1Prt",
                             creek == "flow_brn_by_padden" ~ "2Brn",
                             TRUE ~ creek)) %>% 
    separate(daysampletoavg, c("year","month","day")) %>% 
    mutate(year= str_sub(year, 3, 4))  %>%
    mutate(month = as.numeric(month)) %>%
    unite(time, c(month, year), sep="", remove=FALSE) %>%
    mutate(creek = case_when(creek == "1Prt" ~ "Portage",
                             creek == "2Brn" ~ "Barnes",
                             creek == "3Chk" ~ "Chuckanut",
                             creek == "4Pad" ~ "Padden",
                             creek == "5Sqm" ~ "Squalicum",
                             TRUE ~ creek)) %>% 
    unite(timecreek, c(time, creek), remove=FALSE) %>%
  unite(creekyrmo, c(creek,year,month)) %>%
  mutate(flow_m3s = ifelse(is.na(flow_m3s), 0.00141584235, flow_m3s)) # when no metered flow, use half the lowest verified flow


cut_flow_corrected <- cut_modeled_conc %>%
  mutate(time = case_when(time == "822" ~ "0822",
                          time == "722" ~ "0722",
                          TRUE ~ time)) %>% 
  separate(time, into = c("month","year"), sep = 2) %>%
  unite(newtime, c(year,month), sep="-", remove=FALSE) %>%
  mutate(year = as.numeric(year)) %>%
  mutate(month = as.numeric(month)) %>%
  mutate(creek = case_when(creek == "1Prt" ~ "Portage",
                           creek == "2Brn" ~ "Barnes",
                           creek == "3Chk" ~ "Chuckanut",
                           creek == "4Pad" ~ "Padden",
                           creek == "5Sqm" ~ "Squalicum",
                           TRUE ~ creek)) %>%
  mutate(station = case_when(station == "Dn" ~ "Down",
                             TRUE ~ station)) %>%
  unite(creekyrmo, c(creek,year,month), remove=FALSE) %>%
  left_join(flowmerge, by="creekyrmo") %>% 
  mutate(mean_conc_flow_correct = mean_concentration_est*flow_m3s*1000) %>%
  mutate(ci25_flow_correct = ci25_concentration_est*flow_m3s*1000) %>%
  mutate(ci75_conc_flow_correct = ci25_concentration_est*flow_m3s*1000) %>% 
  mutate(LOQ = (100/(Adj_Vol/2))*1000*flow_m3s)  %>% # 
  mutate(belowLOQ = mean_conc_flow_correct < LOQ) %>% 
  mutate(facetorder = factor(creek, levels=c('Padden','Portage','Chuckanut','Squalicum', 'Barnes')))

#write_rds(cut_flow_corrected, here("Output","qPCR", "modeled_cut_qpcr_flowcorrected_monthlyavg.RDS"))
#write_rds(cut_flow_corrected, here("Output","qPCR", "modeled_cut_qpcr_flowcorrected_allpadden.RDS"))

brnonly_flow <- cut_flow_corrected %>% 
  filter(creek=="Barnes")
chkonly_flow <- cut_flow_corrected %>% 
  filter(creek=="Chuckanut") %>% 
  filter(station=="Down")

cut_flow_corrected  %>% 
  ggplot(aes(x = newtime, y = log(mean_conc_flow_correct), color=station)) +
  geom_point() +
  stat_smooth(method = "loess", span=0.5, aes(group=station)) +
  geom_segment(aes(x = newtime, xend = newtime, y = log(ci25_flow_correct), yend = log(ci75_conc_flow_correct))) +
  geom_point(aes(x=newtime, y=log(LOQ)), color="red", shape=3) +
  geom_rect(data=brnonly_flow, aes(xmin = 12, xmax = 20, ymin = -Inf, ymax = Inf), alpha = .01, color=NA, fill = "grey") +
  geom_rect(data=chkonly_flow, aes(xmin = 12, xmax = 20, ymin = -Inf, ymax = Inf), alpha = .01, color=NA, fill = "grey") +
  facet_grid(rows=vars(facetorder)) +
  labs(x="Date (YY-MM)", y= "Log copies/s (flow corrected)", title = "Cutthroat Trout", color="Station") + 
  theme_bw() + 
  scale_x_discrete(guide = guide_axis(angle = -45)) +
  scale_color_manual(values=c("darkgrey", "black", "black", "blue"))

ggsave(here("Output","Figures", "modeled_cut_qpcr_updown_flowcorrected.png"), units="in", width=5, height=7)


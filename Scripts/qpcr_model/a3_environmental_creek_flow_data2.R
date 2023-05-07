## NGN Flow data
# Author: Eily Allan 
# Person running: Eily
# Last modified: 03/21/23 by Eily
# Date of run: 03/21/23 by Eily 

# Overview 
# This script does more advanced work with flow data from three gauges (obtained by request from the City of Bellingham).
# Here, we make correction factors for the ohter creeks using historical data to fill in missing values. 

# Inputs: 
# 1-3) flow data for Padden Creek, Squalicum Creek, Chuckanut Creek 

# Outputs: 
# 1-3) only relevant dates and converted to m3/s


####################################################################
# Set up
####################################################################

# Load packages
library(here)
library(tidyverse)
library(lubridate) # package that makes handling dates much easier 
library(ggplot2)
library(gridExtra)

source(here("Scripts","functions","flow_functions.R"))

# Read in files
padflow1 <- read.csv(here("Input","qpcr","flow","old", "Padden_Creek_Discharge.Time_Series_Data.2021110822572368.csv"), skip=29)[,1:2]
sqmflow1 <- read.csv(here("Input","qpcr","flow","old","Squalicum_Creek_Discharge.Time_Series_Data.2021110822504744.csv"), skip=27)[,1:2]
chkflow1 <- read.csv(here("Input","qpcr","flow","old","Chuckanut_Creek.Time_Series_Data.2021110823201271.csv"), skip=30)[,1:2]

padflow2 <- read.csv(here("Input","qpcr","flow","Padden_Creek_Computed_Discharge_(cfs).Time_Series_Data.2023041017032436.csv"), skip=30)[,1:2]
sqmflow2 <- read.csv(here("Input","qpcr","flow","Squalicum_Creek_Computed_Discharge_(cfs).Time_Series_Data.2022102422132690.csv"), skip=27)[,1:2]
chkflow2 <- read.csv(here("Input","qpcr","flow","Chuckanut_Creek_Computed_Discharge_(cfs).Time_Series_Data.2022102422194411.csv"), skip=30)[,1:2]

# read in when filtering occurred 
#filtermeta <- readRDS(here("Output","qpcr","adj_vol_filtered.RDS"))
filtermeta <- read.csv(here("Output","qpcr", "adj_vol_filtered.csv"))
wsdotnames <- read.csv(here("Input","qpcr","wsdot_names_wrangle.csv"))

filtermeta2 <- filtermeta %>% 
  dplyr::rename(qPCR_name = Sample) %>%
  left_join(wsdotnames %>% filter(!str_detect(plate,"COHO")), by="qPCR_name") %>%
  dplyr::select(-plate) %>%
  mutate(new_name = case_when(is.na(new_name) ~ qPCR_name,
                              !is.na(new_name) ~ new_name)) %>%
  dplyr::select(-qPCR_name) %>%
  dplyr::rename(Sample=new_name)

####################################################################
# Clean up gauge data
####################################################################
padflow2 <- padflow2 %>% 
  dplyr::rename(Discharge.Fairhaven.Park..cfs. = Discharge.Fairhaven.Park)

padflow <- rbind(padflow1, padflow2)
chkflow <- rbind(chkflow1, chkflow2)
sqmflow <- rbind(sqmflow1, sqmflow2)

pad.flow <- padflow %>% 
  distinct() %>% 
  add_column(creek = "Padden") %>% 
  dplyr::rename(flow_cfs = Discharge.Fairhaven.Park..cfs.) 

chk.flow <- chkflow %>% 
  distinct() %>% 
  add_column(creek = "Chuckanut") %>% 
  dplyr::rename(flow_cfs = Discharge.Arroyo.Park..cfs.) 

sqm.flow <- sqmflow %>% 
  distinct() %>% 
  add_column(creek = "Squalicum") %>% 
  dplyr::rename(flow_cfs = Discharge.West.Street..cfs.) 

####################################################################
# Clean up sampling data
####################################################################

filtermetaplot <- filtermeta2 %>% 
  mutate(DateStamp = gsub("UTC","", DateStamp)) %>% 
  mutate(timeplot = parse_date_time(DateStamp, "ymd_HMS")) %>% 
  dplyr::select(-DateStamp) %>% 
  filter(!str_detect("Up5", Sample)) %>% 
  separate(Sample, into=c("time","creek","station","biol"), remove=FALSE) %>% 
  mutate(creek = case_when(creek == "3Chk" ~ "Chuckanut",
                           creek == "4Pad" ~ "Padden",
                           creek == "5Sqm" ~ "Squalicum",
                           TRUE ~ creek)) %>% 
  dplyr::select(-c(biol,station,time, Adj_Vol)) %>% 
  mutate(timeplot = round_date(timeplot, unit="15 minutes")) %>% 
  mutate(daysampletoavg = floor_date(timeplot, unit="day"))


####################################################################
# Use functions to manipulate and plot 
####################################################################

pad.past.flow <- format_flow(pad.flow)
pad.sampling.flow <- format_flow2(pad.flow)
pad.year.avg <- yearly_avg_flow(pad.past.flow)
pad.discrete <- flow_discrete(filtermetaplot, "Padden", pad.sampling.flow)
pad.ts.plot <- plot_ts(pad.past.flow, "Padden")
pad.ts.facet.plot <- plot_facet_year(pad.past.flow, "Padden")
pad.ts.year.plot <- plot_year_avg(pad.year.avg, "Padden")
pad.ts.year.discrete.plot <- plot_year_avg_discrete(pad.year.avg, pad.discrete %>% filter(str_detect(Sample, "Dn")), "Padden")
pad.monthly.avg <- monthly_avg_flow(pad.past.flow)
pad.ts.monthly.discrete.plot <- plot_year_avg_discrete(pad.monthly.avg, pad.discrete %>% filter(str_detect(Sample, "Dn")), "Padden")
pad.ts.sampling.discrete.plot <-  plot_sampling_discrete(pad.sampling.flow, pad.discrete %>% filter(str_detect(Sample, "Dn")), "Padden")


sqm.past.flow <- format_flow(sqm.flow)
sqm.sampling.flow <- format_flow2(sqm.flow)
sqm.year.avg <- yearly_avg_flow(sqm.past.flow)
sqm.discrete <- flow_discrete(filtermetaplot, "Squalicum", sqm.sampling.flow)
sqm.ts.plot <- plot_ts(sqm.past.flow, "Squalicum")
sqm.ts.facet.plot <- plot_facet_year(sqm.past.flow, "Squalicum")
sqm.ts.year.plot <- plot_year_avg(sqm.year.avg, "Squalicum")
sqm.ts.year.discrete.plot <- plot_year_avg_discrete(sqm.year.avg, sqm.discrete %>% filter(str_detect(Sample, "Dn")), "Squalicum")
sqm.monthly.avg <- monthly_avg_flow(sqm.past.flow)
sqm.ts.monthly.discrete.plot <- plot_year_avg_discrete(sqm.monthly.avg, sqm.discrete %>% filter(str_detect(Sample, "Dn")), "Squalicum")
sqm.ts.sampling.discrete.plot <-  plot_sampling_discrete(sqm.sampling.flow, sqm.discrete %>% filter(str_detect(Sample, "Dn")), "Squalicum")

chk.past.flow <- format_flow(chk.flow)
chk.sampling.flow <- format_flow2(chk.flow)
chk.year.avg <- yearly_avg_flow(chk.past.flow)
chk.discrete <- flow_discrete(filtermetaplot, "Chuckanut", chk.sampling.flow)
chk.ts.plot <- plot_ts(chk.past.flow, "Chuckanut")
chk.ts.facet.plot <- plot_facet_year(chk.past.flow, "Chuckanut")
chk.ts.year.plot <- plot_year_avg(chk.year.avg, "Chuckanut")
chk.ts.year.discrete.plot <- plot_year_avg_discrete(chk.year.avg, chk.discrete %>% filter(str_detect(Sample, "Dn")), "Chuckanut")
chk.monthly.avg <- monthly_avg_flow(chk.past.flow)
chk.ts.monthly.discrete.plot <- plot_year_avg_discrete(chk.monthly.avg, chk.discrete %>% filter(str_detect(Sample, "Dn")), "Chuckanut")
chk.ts.sampling.discrete.plot <-  plot_sampling_discrete(chk.sampling.flow, chk.discrete %>% filter(str_detect(Sample, "Dn")), "Chuckanut")

flow_CFs <- pad.monthly.avg %>%
  filter(datetime > '2022-01-01') %>% 
  dplyr::rename(padflow = val) %>% 
  left_join(chk.monthly.avg, by="datetime") %>% 
  dplyr::rename(chkflow = val) %>% 
  left_join(sqm.monthly.avg, by="datetime") %>% 
  dplyr::rename(sqmflow = val) %>% 
  mutate(chk_cf = chkflow/padflow)%>% 
  mutate(sqm_cf = sqmflow/padflow)

use_padden_and_cfs <- pad.discrete %>% 
  filter(str_detect(Sample, "Dn")) %>% 
  mutate(monthplot = month(timeplot)) %>% 
  mutate(datetime = lubridate::make_datetime(yearplot, monthplot, 15)) %>%
  left_join(flow_CFs %>% dplyr::select(c(monthplot, chk_cf, sqm_cf)), by= "monthplot") %>% 
  mutate(flow_chk_by_padden = flow_m3s*chk_cf) %>% 
  mutate(flow_sqm_by_padden = flow_m3s*sqm_cf) 

compare_chk_sqm <- use_padden_and_cfs %>% 
  dplyr::select(c(daysampletoavg, flow_m3s,flow_chk_by_padden,flow_sqm_by_padden)) %>%
  dplyr::rename(pad_flowm3s = flow_m3s) %>% 
  left_join(chk.discrete %>% dplyr::select(c(daysampletoavg, flow_m3s)), by= "daysampletoavg") %>% 
  dplyr::rename(chk_flowm3s = flow_m3s) %>% 
  left_join(sqm.discrete %>% dplyr::select(c(daysampletoavg, flow_m3s)), by= "daysampletoavg") %>% 
  dplyr::rename(sqm_flowm3s = flow_m3s) %>% 
  distinct()

ggplot(compare_chk_sqm) +
  geom_line(aes(x = daysampletoavg, y = sqm_flowm3s)) +
  geom_point(aes(x = daysampletoavg, y = flow_sqm_by_padden)) +
  #labs(y=bquote('Discharge '(m^3/s)), x="Date", title = creek, colour = "Year") + 
  #facet_wrap(~yearplot, nrow=7) +
  #guides(color='none') +
  theme_bw()

flowratestouse <- compare_chk_sqm %>% 
  dplyr::select(c(daysampletoavg, pad_flowm3s, flow_chk_by_padden,flow_sqm_by_padden)) %>% 
  mutate(flow_brn_by_padden = pad_flowm3s*1) %>%  #for now, just make something up.... 
  mutate(flow_prt_by_padden = pad_flowm3s*1) %>%  #for now, just make something up.... 
  distinct()

cf_plot <- flowratestouse %>% 
  rename(flow_pad = pad_flowm3s) %>% 
  pivot_longer(cols = starts_with("flow"),
               values_to = "flow_m3s") %>% 
  mutate(creek = case_when( name == "flow_pad" ~ "Padden",
                           name == "flow_prt_by_padden" ~ "Portage",
                           name == "flow_brn_by_padden" ~ "Barnes",
                           name == "flow_chk_by_padden" ~ "Chuckanut",
                           name == "flow_sqm_by_padden" ~ "Squalicum")) %>% 
  separate(daysampletoavg, into=c("yearplot", "monthplot", "dayplot")) %>% 
  mutate(yearplot = as.numeric(yearplot), monthplot = as.numeric(monthplot), dayplot = as.numeric(dayplot)) %>% 
  mutate(datetime = lubridate::make_datetime(yearplot, monthplot, dayplot)) 
  

# write.csv(flowratestouse, here("Output","qpcr","flowrates_touse_extended.csv"), row.names=FALSE)

pad.cf <- ggplot(pad.sampling.flow, aes(x = datetime, y = flow_m3s)) + 
  geom_line() + 
  labs(y=bquote('Discharge '(m^3/s)), x="Month", title = "Padden Creek") + 
  theme_bw() +
  geom_point(data=cf_plot %>% filter(creek=="Padden"), aes(x=datetime, y=flow_m3s), size=3, bg="black", pch=21) 

sqm.cf <- ggplot(sqm.sampling.flow, aes(x = timeplot, y = flow_m3s)) + 
  geom_line() + 
  labs(y=bquote('Discharge '(m^3/s)), x="Month", title = "Squalicum Creek") + 
  theme_bw() +
  geom_point(data=cf_plot %>% filter(creek=="Squalicum"), aes(x=datetime, y=flow_m3s), size=3, bg="blue", pch=21) 

chk.cf <- ggplot(chk.sampling.flow, aes(x = timeplot, y = flow_m3s)) + 
  geom_line() + 
  labs(y=bquote('Discharge '(m^3/s)), x="Month", title = "Chuckanut Creek") + 
  theme_bw() +
  geom_point(data=cf_plot %>% filter(creek=="Chuckanut"), aes(x=datetime, y=flow_m3s), size=3, bg="red", pch=21) 

g <- arrangeGrob(pad.cf, chk.cf, sqm.cf, nrow=3) #generates g
ggsave(file=here("Output","Figures","flow_gauges.png"),g, units= "in", height=6, width=8)


####################################################################
# Plot and save  
####################################################################

g <- arrangeGrob(pad.ts.plot, chk.ts.plot, sqm.ts.plot, nrow=3) #generates g
ggsave(file=here("Output","SupplementalFigures","historical_flow.png"),g)

g <- arrangeGrob(pad.ts.year.discrete.plot, chk.ts.year.discrete.plot, sqm.ts.year.discrete.plot, nrow=3) #generates g
ggsave(file=here("Output","SupplementalFigures","historical_flow_year_avg_sampling.png"),g, units="in", height=7, width=9)

# g <- arrangeGrob(pad_cf_plot, chk_cf_plot, sqm_cf_plot, nrow=3) #generates g
# ggsave(file=here("Output","SupplementalFigures","flow_by_cfs.png"),g)


####################################################################
# Compare closest measurement, daily average, monthly average, historical average (?)  
####################################################################

plot_monthly_avg <- read.csv(here("Output", "qpcr", "monthly_flow.csv"))
plot_daily_avg <-  read.csv(here("Output", "qpcr", "daily_flow.csv"))
plot_closest <- read.csv(here("Output", "qpcr", "closest_flow.csv"))

plot_monthly_avg <- plot_monthly_avg %>%  
  dplyr::rename(flow=meanmonthflow) %>% 
  select(-timestring) %>% 
  mutate(dayplot = 15) %>% 
  mutate(type="Monthly Average")

plot_daily_avg <- plot_daily_avg %>%  
  dplyr::rename(flow=meandayflow, creek = creek.x) %>% 
  select(-c(Sample, creek.y, daysampletoavg)) %>% 
  mutate(yearplot = year(timeplot)) %>% 
  mutate(monthplot = month(timeplot)) %>% 
  mutate(dayplot = 15) %>% 
  mutate(type="Daily Average")

plot_closest <- plot_closest %>% 
  filter(!timeplot=="2021-04-08 22:15:00") %>% 
  select(-Sample) %>% 
  distinct() %>% 
  dplyr::rename(flow=flow_m3s) %>% 
  select(-daysampletoavg) %>% 
  mutate(dayplot = 15) %>% 
  mutate(type="Closest Value")

plot_cf_values <- flowratestouse %>% 
  pivot_longer(!daysampletoavg, names_to="creek", values_to="flow") %>% 
  mutate(creek = case_when(creek == "flow_chk_by_padden" ~ "Chuckanut",
                           creek == "flow_sqm_by_padden" ~ "Squalicum",
                           creek == "pad_flowm3s" ~ "Padden")) %>% 
  drop_na() %>% 
  mutate(type="Correction Factor") %>% 
  dplyr::rename(timeplot = daysampletoavg) %>% 
  mutate(yearplot = year(timeplot)) %>% 
  mutate(monthplot = month(timeplot)) %>% 
  mutate(dayplot = 15)

plot_all <- rbind(plot_monthly_avg, plot_daily_avg, plot_closest, plot_cf_values) 


plot_all %>% 
  mutate(datetime = make_date(yearplot, monthplot,dayplot)) %>% 
  filter(datetime > "2021-03-10") %>% 
  distinct() %>% 
  ggplot(aes(x=datetime, y=flow, color=type), alpha=0.5) + 
  geom_point(alpha=0.8, size=2) +
  scale_x_date(date_labels = "%b %y") +
  scale_color_manual(values=c("#999999","#E69F00","#56B4E9","#009E73")) +
  theme_bw() +
  facet_wrap(~creek, nrow=1, ncol=3) +
  labs(y=bquote('Discharge '(m^3/s)), x="Date of Sampling", color="Source of Flow Data")

ggsave(here("Output","SupplementalFigures","flow_by_source.png"), units="in", width=10, height=5)

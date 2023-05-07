#spline model for NGN resubmit
library(tidyverse)
library(here)
library(mgcv)
library(brms)
library(gratia)
library(tidybayes)
library(rstan)
library(rstanarm)
options(mc.cores = parallel::detectCores())


d <- #readRDS(here("Output","20230323_abundance_flowcorrected.RDS")) %>% 
  readRDS(here("Output/20230418_abundance_flowcorrected.RDS")) %>% 
  filter(creek != "Barnes") %>%   #omit Barnes because we have few datapoints for that creek in terms of abs concentration
  mutate(station_idx = case_when(station == "Up" ~ 2, 
                                 station == "Down" ~ 1,
                                 station == "Up5" ~ 3),
         bio = as.integer(bio)) %>% 
  dplyr::select(-station) %>% 
  arrange(creek) %>% 
  arrange(newtime) %>% 
  mutate(station = as.factor(station_idx))
d$creek_idx <- match(d$creek, unique(d$creek))
d$time_idx <- match(d$newtime, unique(d$newtime))
d$species_idx <- match(d$species, unique(d$species))

q <- d %>% 
  filter(species == "Oncorhynchus clarkii") %>% 
  filter(creek == "Padden") 
  
  
#single-species version w intercept varying by station
m1 <- stan_glmer(log(meandnaconcflow) ~ newtime + (1|station),
                 data = q)
  
  q %>%
  add_predicted_draws(m1) %>%
  ggplot(aes(y = log(meandnaconcflow), x= time_idx)) +
  stat_lineribbon(aes(y = .prediction)) +
  geom_point(data = q) +
  scale_fill_brewer(palette = "Greys") +
  facet_grid(~station_idx)

#single-species hierarchical version, effect of time varying by station
m2 <- stan_glmer(log(meandnaconcflow) ~ (1 + newtime|station),
                 data = q)

q %>%
  add_predicted_draws(m2) %>%
  ggplot(aes(y = log(meandnaconcflow), x= time_idx)) +
  stat_lineribbon(aes(y = .prediction)) +
  geom_point(data = q) +
  scale_fill_brewer(palette = "Greys") +
  facet_grid(~station)

loo_compare(loo(m1), loo(m2))

#####multicreek

f <- d %>% 
  filter(species == "Oncorhynchus clarkii") %>% 
  mutate(station = as.factor(station_idx))

m3 <- stan_glmer(log(meandnaconcflow) ~ (newtime|creek) + (1|station:creek),
                 data = f)

f %>%
  add_predicted_draws(m3) %>%
  ggplot(aes(y = log(meandnaconcflow), x= time_idx, color = station)) +
  stat_lineribbon(aes(y = .prediction)) +
  geom_point(data = f) +
  scale_fill_brewer(palette = "Greys") +
  facet_grid(~creek)


##multicreek, multispecies -- effect of month varies by creek and species, plus intercept for each station/creek/species

flowdata <- read.csv(here("Output","qpcr","flowrates_touse_extended.csv"))

LOQ_creek_month <- flowdata %>% 
  pivot_longer(cols=2:6, names_to = "creek", values_to = "flowrate") %>% 
  mutate(creek = case_when(creek == "pad_flowm3s" ~ "Padden",
                           creek == "flow_chk_by_padden" ~ "Chuckanut",
                           creek == "flow_sqm_by_padden" ~ "Squalicum", 
                           creek == "flow_brn_by_padden" ~ "Barnes",
                           creek == "flow_prt_by_padden" ~ "Portage")) %>% 
  separate(daysampletoavg, c("year","month","day")) %>% 
  mutate(year= str_sub(year, 3, 4))  %>%
  unite(newtime, c(year,month), sep="-") %>% 
  select(-day) %>% 
  mutate(LOQ = log(flowrate*100*1000)) %>% 
  select(-flowrate)

aboveloq <- d %>%
  left_join(LOQ_creek_month, by = c("newtime", "creek")) %>% 
  mutate(logconc = log(meandnaconcflow)) %>% 
  mutate(tomodel = logconc > LOQ) %>% 
  filter(tomodel == TRUE) 

m4 <- stan_glmer(log(meandnaconcflow) ~ (1 + newtime|creek:species) + (1|station:creek:species:newtime),
                 data = d)



summary(m4, pars="b[(Intercept) station:creek:species:newtime:1:Chuckanut:Oncorhynchus_clarkii:21-03]")

test <- as.data.frame(m4)
testnames <- colnames(test)
 
test3 <- m4$stan_summary
names <- rownames(test3)
rownames(test3) <- NULL
test3 <- cbind(names,test3)
test3 <- as.tibble(test3)

etas <- test3 %>% 
  filter(str_detect(names, "Intercept")) %>% 
  filter(!str_detect(names, "Sigma")) %>% 
  filter(names != "(Intercept)") %>% 
  filter(!str_detect(names,"NEW")) %>% 
  separate(names, into=c("a","b","c","d","e","f","g","h","i","j","k","l","m")) %>% 
  select(-c("a","b","c","d","e","f","i","m")) %>% 
  rename(station = g, creek = h, species = j, year = k, month = l) %>% 
  unite(newtime, c("year","month"), sep="-") %>% 
  drop_na() %>% 
  mutate(mean = as.numeric(mean)) %>% 
  mutate(p025 = as.numeric(`2.5%`)) %>% 
  mutate(p975 = as.numeric(`97.5%`))


etas %>% 
  filter(creek=="Padden") %>% 
ggplot(aes(x=newtime, y=mean)) +
  geom_point() + 
  geom_segment(aes(x = newtime, xend = newtime, y = p025, yend = p975)) +
  facet_grid(~species~station)

deltas <- etas %>% 
  unite(creeksptime, c(creek,species,newtime)) %>% 
  select(c(creeksptime,station,mean)) %>% 
  pivot_wider(names_from=station, values_from=mean) %>% 
  mutate(diff21 = `2` - `1`) %>% 
  mutate(diff32 = `3` - `2`) %>% 
  mutate(percent21 = exp(diff21)) %>% 
  mutate(percent32 = exp(diff32)) %>% 
  separate(creeksptime, into=c("creek","species","year", "month")) %>% 
  unite(newtime, c("year","month"), sep="-")

flowfordeltas <- flowtojoin %>% 
  select(-day) %>% 
  separate(timecreek, c("time","creek")) %>% 
  separate(time, c("month", "year"), sep=2) %>% 
  mutate(creek = case_when(creek == "4Pad" ~ "Padden",
                           creek == "3Chk" ~ "Chuckanut",
                           creek == "5Sqm" ~ "Squalicum", 
                           creek == "2Brn" ~ "Barnes",
                           creek == "1Prt" ~ "Portage")) %>% 
  unite(newtime, c("year","month"), sep="-")

deltas <- deltas %>% 
  left_join(flowfordeltas, by=c("newtime","creek"))

deltas %>% 
  ggplot(aes(x=newtime, y=percent21, size=flow_m3s)) +
  geom_point() + 
  geom_point(aes(x = newtime,y=percent32), color="red") +
  geom_hline(yintercept=1, linetype=2) +
  facet_grid(~species~creek) + 
  theme_bw() +
  scale_x_discrete(guide = guide_axis(angle = -45)) 

deltas %>% 
  #filter(species !="nerka") %>% 
  ggplot(aes(x=log(flow_m3s), y=percent21, color=species)) +
  geom_point() + 
  facet_wrap(~creek, scales="free") +
  theme_bw()


m4loq_noup5 <- stan_glmer(log(meandnaconcflow) ~ (newtime|creek:species) + (1|station:creek:species),
                 data = aboveloq)

m4loq_nodn <- stan_glmer(log(meandnaconcflow) ~ (newtime|creek:species) + (1|station:creek:species),
                          data = aboveloq)


d %>%
  left_join(LOQ_creek_month, by = c("newtime", "creek")) %>% 
  add_predicted_draws(m4) %>%
  ggplot(aes(y = log(meandnaconcflow), x= time_idx, color = station)) +
  stat_lineribbon(aes(y = .prediction)) +
  geom_line(aes(x=time_idx, y=LOQ)) +
  #geom_line(data=LOQ_creek_month, aes(x=time_idx, y=LOQ)) +
  #geom_point(data = d) +
  scale_fill_brewer(palette = "Greys") +
  facet_grid(species~creek)


d %>%
  left_join(LOQ_creek_month, by = c("newtime", "creek")) %>% 
  #add_predicted_draws(m4) %>%
  ggplot(aes(y = log(meandnaconcflow), x= time_idx, color = station)) +
  #stat_lineribbon(aes(y = .prediction)) +
  #geom_line(aes(x=time_idx, y=LOQ)) +
  geom_point(data = d) +
  #scale_fill_brewer(palette = "Greys") +
  facet_grid(species~creek)


d %>%
  left_join(LOQ_creek_month, by = c("newtime", "creek")) %>% 
  #add_predicted_draws(m4loq) %>%
  ggplot(aes(y = log(meandnaconcflow), x= time_idx, color = station)) +
  #stat_lineribbon(aes(y = .prediction)) +
  #geom_line(aes(x=time_idx, y=LOQ)) +
  geom_point(data = aboveloq) +
  scale_fill_brewer(palette = "Greys") +
  facet_grid(species~creek)

d %>%
  left_join(LOQ_creek_month, by = c("newtime", "creek")) %>% 
  #add_predicted_draws(m4loq) %>%
  ggplot(aes(y = log(meandnaconc), x= time_idx, color = station)) +
  #stat_lineribbon(aes(y = exp(.prediction))) +
  geom_line(aes(x=time_idx, y=LOQ)) +
  geom_point(data = aboveloq) +
  scale_fill_brewer(palette = "Greys") +
  facet_grid(species~creek) 


summary(m4loq)



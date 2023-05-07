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

d <- readRDS(here("Output/20230505_abundance_flowcorrected.RDS")) %>% 
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

### SET UP MATRICES FOR TIMES NOT SAMPLED OR SAMPLED BUT NOT DOWN / UP  

alltimes <- c("21-03", "21-04", "21-05", "21-06", "21-07", "21-08", "21-09", "21-10", "21-11", "21-12", 
              "22-01", "22-02", "22-03", "22-04", "22-05", "22-06", "22-07", "22-08", "22-12")
allspecies <- unique(d$species)
allcreeks <- unique(d$creek)
notsampledtimes <- c("22-09", "22-10", "22-11")
notsampledtimes2 <- c("22-03", "22-04", "22-05", "22-06", "22-07", "22-08", "22-12")

notsampled <- expand.grid(notsampledtimes, allspecies, allcreeks)  
notsampled <- notsampled %>% 
  mutate(value=-5) %>% 
  rename(newtime=Var1, species=Var2, creek=Var3) %>% 
  mutate(species = case_when(species == "Oncorhynchus clarkii" ~ "clarkii",
                           species == "Oncorhynchus mykiss" ~ "mykiss",
                           species == "Oncorhynchus kisutch" ~ "kisutch",
                           species == "Oncorhynchus nerka" ~ "nerka",
                           TRUE ~ species))

notsampled2 <- expand.grid(notsampledtimes2, allspecies, c("Chuckanut","Barnes"))
notsampled2 <- notsampled2 %>% 
  mutate(value=-5) %>% 
  rename(newtime=Var1, species=Var2, creek=Var3) %>% 
  mutate(species = case_when(species == "Oncorhynchus clarkii" ~ "clarkii",
                             species == "Oncorhynchus mykiss" ~ "mykiss",
                             species == "Oncorhynchus kisutch" ~ "kisutch",
                             species == "Oncorhynchus nerka" ~ "nerka",
                             TRUE ~ species))

### LOOK AT FLOW DATA TO DETERMINE LOQ FOR EACH CREEK / TIME POINT

flowdata <- read.csv(here("Output","qpcr","flowrates_touse_extended.csv"))
flowtojoin <- flowdata %>% 
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
  mutate(flowrate = ifelse(is.na(flowrate), 0.00141584235, flowrate))%>% 
  mutate(LOQ = log(flowrate*100*1000)) %>% 
  select(-flowrate)

aboveloq <- d %>%
  left_join(LOQ_creek_month, by = c("newtime", "creek")) %>% 
  mutate(logconc = log(meandnaconcflow)) %>% 
  mutate(tomodel = logconc > LOQ) %>% 
  filter(tomodel == TRUE) 

### RUN LINEAR MODEL 

m4 <- stan_glmer(log(meandnaconcflow) ~ (1 + newtime|creek:species) + (1|station:creek:species:newtime),
                 data = aboveloq)

b <- as.data.frame(m4$stanfit)
dd <- posterior_predict(m4)
e <- aboveloq %>% 
  bind_cols(summarise_draws(dd))

## only padden for shading construction
e_pad <- expand.grid(alltimes, unique(e$species), unique(e$station)) %>% 
  mutate(creek="Padden") %>% 
  rename(newtime=Var1, species=Var2, station=Var3)

## only chuckanut for times not sampled
e_chk <- expand.grid(alltimes, unique(e$species), unique(e$station)) %>% 
  mutate(creek="Chuckanut") %>% 
  rename(newtime=Var1, species=Var2, station=Var3)

## only squalicum for times not sampled 
e_brn <- expand.grid(alltimes, unique(e$species), unique(e$station)) %>% 
  mutate(creek="Barnes") %>% 
  rename(newtime=Var1, species=Var2, station=Var3)


multispeciesTrendsplot_simple <- ggplot() +  #treating meandnaconc as observations
  geom_rect(data=e_pad, aes(xmin = 6, xmax = 8, ymin = -Inf, ymax = Inf), alpha = .01, color=NA, fill = "wheat") +
  geom_rect(data=e_pad, aes(xmin = 16, xmax = 17, ymin = -Inf, ymax = Inf), alpha = .01, color=NA, fill = "wheat") +
  geom_rect(data=e_chk, aes(xmin = 12, xmax = 20, ymin = -Inf, ymax = Inf), alpha = .01, color=NA, fill = "grey") +
  geom_rect(data=e_brn, aes(xmin = 12, xmax = 20, ymin = -Inf, ymax = Inf), alpha = .01, color=NA, fill = "grey") +
  geom_point(data=e, aes(x = newtime, y = mean, color = station), size = 1.5, alpha = .7) + #mean estimates of mu
  #geom_smooth(data=e, aes(x = newtime, y = mean, color = station), se = F, size = 1, alpha = .5, span = .2) +
  geom_segment(data=e, aes(x = newtime, xend = newtime, y = q5, yend = q95, color = station), size = .7, alpha = .4) +
  facet_grid(species~creek) +
  xlab("Date (YY-MM)") + ylab("Log DNA Concentration (copies/s)") +
  labs(color='Station') +
  scale_color_manual(values = c("grey", "black", "blue"), labels = c('Downstream', 'Upstream (SR-11)', 'Upstream (I-5)')) + 
  # ggtitle(paste0(unique(f$species)[i], "\n(predicted mean +/- 95 CI)")) +
  theme_bw() +
  theme(panel.grid.minor = element_blank(),
        panel.grid.major.x = element_blank(),
        legend.position = "bottom",
        legend.box.spacing = unit(0, "mm")) + 
  scale_x_discrete(guide = guide_axis(angle = -45))  

ggsave(multispeciesTrendsplot_simple, file = here("Output/Figures/multispeciesTrends_flowcorrected.png"), width=12, height=7)


### NOW LOOK AT CULVERT EFFECT 

etas <- b[,grep("station:creek:species:newtime:", colnames(b))]

etas_post <- as_tibble(cbind(name = names(etas), t(etas)))  %>% 
  separate(name, into=c("a","b","c","d","e","f","g","h","i","j","k","l","m")) %>% 
  select(-c("a","b","c","d","e","f","i","m")) %>% 
  rename(station = g, creek = h, species = j, year = k, month = l) %>% 
  unite(newtime, c("year","month"), sep="-") %>% 
  filter(station != "NEW") %>% 
  filter(station !="Intercept")

etas_post_nest <- etas_post %>% 
  group_by(newtime,creek,species,station) %>% 
  nest()

downs <- etas_post_nest %>% 
  filter(station==1) %>% 
  ungroup() %>% 
  select(-station) %>% 
  dplyr::rename(downpost = data) %>% 
  rowwise() %>% 
  mutate(downpost = list(as.numeric(downpost, use.names=FALSE)))

ups <- etas_post_nest %>% 
  filter(station==2) %>% 
  ungroup() %>% 
  select(-station) %>% 
  dplyr::rename(uppost = data) %>% 
  rowwise() %>% 
  mutate(uppost = list(as.numeric(uppost, use.names=FALSE)))

up5s <- etas_post_nest %>% 
  filter(station==3) %>% 
  ungroup() %>% 
  select(-station) %>% 
  dplyr::rename(up5post = data) %>% 
  rowwise() %>% 
  mutate(up5post = list(as.numeric(up5post, use.names=FALSE)))

dud <- downs %>% 
  left_join(ups, by=c("newtime","creek","species")) %>% 
  filter(!is.null(uppost)) %>% 
  filter(!is.null(downpost)) %>% 
  rowwise() %>% 
  mutate(up_down_diff = list(as_tibble(uppost - downpost))) %>% 
  select(-c(uppost,downpost)) %>% 
  unnest(cols = c(up_down_diff)) %>% 
  group_by(creek,species,newtime) %>% 
  mutate(meandiff = mean(value)) %>% 
  unite(creektimespecies, c("creek","newtime","species"), remove=FALSE) %>% 
  left_join(flowfordeltas, by=c("newtime","creek"))

du511 <- ups %>% 
  left_join(up5s, by=c("newtime","creek","species")) %>% 
  filter(!is.null(uppost)) %>% 
  filter(!is.null(up5post)) %>% 
  rowwise() %>% 
  mutate(up5_up11_diff = list(as_tibble(up5post - uppost))) %>% 
  select(-c(uppost,up5post)) %>% 
  unnest(cols = c(up5_up11_diff)) %>% 
  group_by(creek,species,newtime) %>% 
  mutate(meandiff = mean(value)) %>% 
  unite(creektimespecies, c("creek","newtime","species"), remove=FALSE) %>% 
  left_join(flowfordeltas, by=c("newtime","creek"))

du_zeros <- expand.grid(alltimes, allspecies, allcreeks)  
du_zeros <- du_zeros %>% 
  mutate(value=-5) %>% 
  rename(newtime=Var1, species=Var2, creek=Var3) %>% 
  unite(creektimespecies, c("creek","newtime","species"), remove=FALSE) %>% 
  filter(! creektimespecies %in% dud$creektimespecies) %>% 
  mutate(species = case_when(species == "Oncorhynchus clarkii" ~ "clarkii",
                             species == "Oncorhynchus mykiss" ~ "mykiss",
                             species == "Oncorhynchus kisutch" ~ "kisutch",
                             species == "Oncorhynchus nerka" ~ "nerka",
                             TRUE ~ species))

du511_zeros <- du_zeros %>% 
  filter(! creektimespecies %in% du511$creektimespecies) %>% 
  filter(creek == "Padden") %>% 
  mutate(species = case_when(species == "Oncorhynchus clarkii" ~ "clarkii",
                             species == "Oncorhynchus mykiss" ~ "mykiss",
                             species == "Oncorhynchus kisutch" ~ "kisutch",
                             species == "Oncorhynchus nerka" ~ "nerka",
                             TRUE ~ species))


#### PLOTS 

## Difference between upstream and downstream for all creeks and species 

# Box plot with all smashed down --- only passable culverts (Chk, Sqm, Prt)
dud %>% 
  filter(creek %in% c("Chuckanut", "Squalicum","Portage")) %>% 
  ggplot(aes(x = newtime, y = value)) +
  geom_boxplot(fill = "grey") +
  #facet_grid(~creek~species) +
  geom_hline(yintercept=0, color="black", linetype=2) +
  geom_point(data=du_zeros %>% filter(creek %in% c("Chuckanut", "Squalicum","Portage")), aes(x=newtime, y=value), pch=8) + 
  geom_point(data=notsampled %>% filter(creek %in% c("Chuckanut", "Squalicum","Portage")), aes(x=newtime, y=value), pch=8, color="grey") + 
  geom_vline(xintercept=12, linetype=3) + 
  annotate("label", x = 6, y = 6, label = "n = 3 creeks") +
  annotate("label", x = 17, y = 6, label = "n = 2 creeks") +
  #geom_point(data=notsampled2, aes(x=newtime, y=value), color="grey") + 
  xlab("Date (YY-MM)") + ylab("Log-Fold Change \n Upstream to Downstream") +
  theme_bw() + 
  scale_x_discrete(guide = guide_axis(angle = -45))
ggsave(filename=here("Output","Figures","culvert_boxplot_passable.png"), units="in", width=7, height=5)

# Box plot with separated (for supplemental) 
xlabels <- sort(unique(c(alltimes, notsampledtimes)))
xlabels[seq(2, length(xlabels), 2)] <- ""
pad_only <- dud %>% filter(creek=="Padden")

dud %>% 
  ggplot(aes(x = newtime, y = value)) +
  geom_boxplot(fill = "grey") +
  facet_grid(~species~creek) +
  geom_rect(data=pad_only, aes(xmin = 6, xmax = 8, ymin = -Inf, ymax = Inf), alpha = .01, color=NA, fill = "wheat") +
  geom_rect(data=pad_only, aes(xmin = 17, xmax = 20, ymin = -Inf, ymax = Inf), alpha = .01, color=NA, fill = "wheat") +
  geom_hline(yintercept=0, color="black", linetype=2) +
  geom_point(data=du_zeros, aes(x=newtime, y=value)) + 
  geom_point(data=notsampled, aes(x=newtime, y=value), color="grey") + 
  geom_point(data=notsampled2, aes(x=newtime, y=value), color="grey") + 
  xlab("Date (YY-MM)") + ylab("Log-Fold Change \n Upstream to Downstream") +
  theme_bw() + 
  scale_x_discrete(labels = xlabels, guide = guide_axis(angle = -45))
ggsave(filename=here("Output","SupplementalFigures","culvert_boxplot_separated.png"), units="in", width=13, height=4)

# Plotting means to show relationship to flow 
dud %>% 
  ggplot() +
  geom_point(aes(x=newtime, y=meandiff, size=flow_m3s)) + 
  #geom_point(aes(x = newtime,y=percent32), color="red") +
  geom_hline(yintercept=0, linetype=2) +
  geom_point(data=du_zeros, aes(x=newtime, y=value)) + 
  geom_point(data=notsampled, aes(x=newtime, y=value), color="grey") + 
  geom_point(data=notsampled2, aes(x=newtime, y=value), color="grey") + 
  geom_rect(data=pad_only, aes(xmin = 6, xmax = 8, ymin = -Inf, ymax = Inf), alpha = .01, color=NA, fill = "wheat") +
  geom_rect(data=pad_only, aes(xmin = 17, xmax = 20, ymin = -Inf, ymax = Inf), alpha = .01, color=NA, fill = "wheat") +
  facet_grid(~species~creek) + 
  theme_bw() +
  scale_x_discrete(labels = xlabels, guide = guide_axis(angle = -45)) 
ggsave(filename=here("Output","SupplementalFigures","culvert_means_flow_separated.png"), units="in", width=12, height=5)

# kind of busy so remove the times when not sampled and collapse by species / creek
dud %>% 
  ggplot() +
  geom_point(aes(x=newtime, y=meandiff, size=flow_m3s, color=creek)) + 
  geom_hline(yintercept=0, linetype=2) +
  #geom_point(data=dud_zeros, aes(x=newtime, y=value)) + 
  #geom_point(data=notsampled, aes(x=newtime, y=value), color="grey") + 
  #geom_point(data=notsampled2, aes(x=newtime, y=value), color="grey") + 
  geom_rect(data=pad_only, aes(xmin = 6, xmax = 8, ymin = -Inf, ymax = Inf), alpha = .01, color=NA, fill = "wheat") +
  geom_rect(data=pad_only, aes(xmin = 17, xmax = 20, ymin = -Inf, ymax = Inf), alpha = .01, color=NA, fill = "wheat") +
  #facet_grid(~species~creek) + 
  labs(x="Date (YY-MM)", y="Log-Fold Change \n Upstream to Downstream", size="Flow Rate (m3/s)", color="Creek") + 
  theme_bw() +
  scale_x_discrete(labels = xlabels, guide = guide_axis(angle = -45)) 
ggsave(filename=here("Output","Figures","culvert_means_flow.png"), units="in", width=7, height=5)

  
### CONTROL CREEKS VERSUS PADDEN 

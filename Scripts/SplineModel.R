#spline model for NGN resubmit
library(tidyverse)
library(here)
library(mgcv)
library(brms)
library(gratia)
library(tidybayes)
library(ggdist)
library(rstan)
options(mc.cores = parallel::detectCores())


d <- readRDS(here("Output","20230323_abundance_flowcorrected.RDS")) %>% 
  #readRDS("/Users/rpk/Desktop/NGNresubmit 2/Output/20230323_abundance_flowcorrected.RDS") %>% 
  filter(creek != "Barnes") %>%   #omit Barnes because we have few datapoints for that creek in terms of abs concentration
  mutate(station_idx = ifelse(station == "Up", 2, 1),
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
m1 <- brm(bf(log(meandnaconcflow) ~ 
               station + 
               s(time_idx, bs="cc", k =8)),
          data = q, family = gaussian(), cores = 4, seed = 17,
          iter = 4000, warmup = 1000, thin = 10, refresh = 0,
          control = list(adapt_delta = 0.99),
          silent = 0,
          save_model = "mod_m1.stan")
q %>%
  add_predicted_draws(m1) %>%
  ggplot(aes(y = log(meandnaconcflow), x= time_idx)) +
  stat_lineribbon(aes(y = .prediction)) +
  geom_point(data = q) +
  scale_fill_brewer(palette = "Greys") +
  facet_grid(~station_idx)



#single-species hierarchical version, on stations
m2 <- brm(bf(log(meandnaconcflow) ~ 
               s(time_idx, bs="cc") + 
                s(time_idx, by=station, m=1, bs="cc") +  #main effect, station
                s(station, bs="re")), #random effect, station
          data = q, family = gaussian(), cores = 4, seed = 17,
          iter = 4000, warmup = 1000, thin = 10, refresh = 0,
          control = list(adapt_delta = 0.99),
          silent = 0,
          save_model = "mod_m2.stan")


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


m3 <- brm(bf(log(meandnaconcflow) ~ 
               s(time_idx, bs="cc") + 
               s(time_idx, by=station, m=1, bs="cc")+  #main effect, station
               s(station, bs="re") +  #random effect, station
               s(time_idx, by=creek, m=1, bs="cc")+  #main effect, creek
               s(creek, bs="re")), #random effect, creek
          data = f, family = gaussian(), cores = 4, seed = 17,
          iter = 4000, warmup = 1000, thin = 10, refresh = 100,
          control = list(adapt_delta = 0.99),
          silent = 0,
          save_model = "mod_m3.stan")

f %>%
  add_predicted_draws(m3) %>%
  ggplot(aes(y = log(meandnaconcflow), x= time_idx)) +
  stat_lineribbon(aes(y = .prediction)) +
  geom_point(data = f) +
  scale_fill_brewer(palette = "Greys") +
  facet_grid(station~creek)


##multicreek, multispecies
m4 <- brm(bf(log(meandnaconcflow) ~ 
               s(time_idx, bs="cc") + 
               s(time_idx, by=station, m=1, bs="cc")+  #main effect, station
               s(station, bs="re") +  #random effect, station
               s(time_idx, by=creek, m=1, bs="cc")+  #main effect, creek
               s(creek, bs="re") + #random effect, creek
               s(time_idx, by=species, m=1, bs="cc")+  #main effect, species
               s(species, bs="re")), #random effect, species
          data = d, family = gaussian(), cores = 4, seed = 17,
          iter = 1000, warmup = 500, thin = 10, refresh = 100,
          control = list(adapt_delta = 0.99),
          silent = 0,
          save_model = "mod_m4.stan")

d %>%
  add_predicted_draws(m4) %>%
  ggplot(aes(y = log(meandnaconcflow), x= time_idx, color = station)) +
  stat_lineribbon(aes(y = .prediction)) +
  geom_point(data = d) +
  scale_fill_brewer(palette = "Greys") +
  facet_grid(species~creek)


##multicreek, indiv species
m5_clarkii <- brm(bf(log(meandnaconcflow) ~ 
               s(time_idx, bs="cc") + 
               s(time_idx, by=station, m=1, bs="cc")+  #main effect, station
               s(station, bs="re") +  #random effect, station
               s(time_idx, by=creek, m=1, bs="cc")+  #main effect, creek
               s(creek, bs="re")),  #random effect, creek
          data = d %>% filter(species == "Oncorhynchus clarkii"), family = gaussian(), cores = 4, seed = 17,
          iter = 1000, warmup = 500, thin = 10, refresh = 100,
          control = list(adapt_delta = 0.99),
          silent = 0,
          save_model = "m5_clarkii.stan")

d %>%
  filter(species == "Oncorhynchus clarkii") %>% 
  add_predicted_draws(m5_clarkii) %>%
  ggplot(aes(y = log(meandnaconc), x= time_idx, color = station)) +
  stat_lineribbon(aes(y = .prediction)) +
  geom_point(data = d %>% filter(species == "Oncorhynchus clarkii")) +
  scale_fill_brewer(palette = "Greys") +
  facet_grid(.~creek)

m5_mykiss <- brm(bf(log(meandnaconc) ~ 
                       s(time_idx, bs="cc") + 
                       s(time_idx, by=station, m=1, bs="cc")+  #main effect, station
                       s(station, bs="re") +  #random effect, station
                       s(time_idx, by=creek, m=1, bs="cc")+  #main effect, creek
                       s(creek, bs="re")),  #random effect, creek
                  data = d %>% filter(species == "Oncorhynchus mykiss"), family = gaussian(), cores = 4, seed = 17,
                  iter = 1000, warmup = 500, thin = 10, refresh = 100,
                  control = list(adapt_delta = 0.99),
                  silent = 0,
                  save_model = "m5_mykiss.stan")

df <- m5_mykiss$fit@sim$samples %>% as.data.frame()
culvert_effect(m5_mykiss, time_idx = 1)

d %>%
  filter(species == "Oncorhynchus mykiss") %>% 
  add_predicted_draws(m5_mykiss) %>%
  ggplot(aes(y = log(meandnaconc), x= time_idx, color = station)) +
  stat_lineribbon(aes(y = .prediction)) +
  geom_point(data = d %>% filter(species == "Oncorhynchus mykiss")) +
  scale_fill_brewer(palette = "Greys") +
  facet_grid(.~creek)

m5_kisutch <- brm(bf(log(meandnaconc) ~ 
                      s(time_idx, bs="cc") + 
                      s(time_idx, by=station, m=1, bs="cc")+  #main effect, station
                      s(station, bs="re") +  #random effect, station
                      s(time_idx, by=creek, m=1, bs="cc")+  #main effect, creek
                      s(creek, bs="re")),  #random effect, creek
                 data = d %>% filter(species == "Oncorhynchus kisutch"), family = gaussian(), cores = 4, seed = 17,
                 iter = 1000, warmup = 500, thin = 10, refresh = 100,
                 control = list(adapt_delta = 0.99),
                 silent = 0,
                 save_model = "m5_kisutch.stan")

d %>%
  filter(species == "Oncorhynchus kisutch") %>% 
  add_predicted_draws(m5_kisutch) %>%
  ggplot(aes(y = log(meandnaconc), x= time_idx, color = station)) +
  stat_lineribbon(aes(y = .prediction)) +
  geom_point(data = d %>% filter(species == "Oncorhynchus kisutch")) +
  scale_fill_brewer(palette = "Greys") +
  facet_grid(.~creek)


m5_nerka <- brm(bf(log(meandnaconc) ~ 
                       s(time_idx, bs="cc") + 
                       s(time_idx, by=station, m=1, bs="cc")+  #main effect, station
                       s(station, bs="re")),   #random effect, station
                       # s(time_idx, by=creek, m=1, bs="cc")+  #main effect, creek
                       # s(creek, bs="re")),  #random effect, creek
                  data = d %>% filter(species == "Oncorhynchus nerka"), family = gaussian(), cores = 4, seed = 17,
                  iter = 1000, warmup = 500, thin = 10, refresh = 100,
                  control = list(adapt_delta = 0.99),
                  silent = 0,
                  save_model = "m5_nerka.stan")

d %>%
  filter(species == "Oncorhynchus nerka") %>% 
  add_predicted_draws(m5_nerka) %>%
  ggplot(aes(y = log(meandnaconc), x= time_idx, color = station)) +
  stat_lineribbon(aes(y = .prediction)) +
  geom_point(data = d %>% filter(species == "Oncorhynchus nerka")) +
  scale_fill_brewer(palette = "Greys") +
  facet_grid(.~creek)


culvert_effect <- function(model, time_idx, species, creek){
  temp <- brms::posterior_predict(model, 
                                  newdata = expand_grid("time_idx" = time_idx,
                                                        "species" = species,
                                                        "creek" = creek,
                                                        "station" = as.factor(c(1,2))))
  return((temp[,2]-temp[,1])/temp[,2]*100) #difference upstream - downstream
  
}


all_culverts <- matrix(nrow=192, ncol=4)

for (i in 1:12) {
  for (j in 1:4) {
    species = case_when(j==1 ~ "Oncorhynchus clarkii",
                        j==2 ~ "Oncorhynchus kisutch",
                        j==3 ~ "Oncorhynchus mykiss", 
                        j==4 ~ "Oncorhynchus nerka")
    for (k in 1:4) {
      creek = case_when(k==1 ~ "Chuckanut",
                         k==2 ~ "Padden",
                          k==3 ~ "Portage", 
                          k==4 ~ "Squalicum")
      counter = 1
      this_culvert <- culvert_effect(six,i,species,creek)
      all_culverts[counter,1] <- i
      all_culverts[counter,2] <- species
      all_culverts[counter,3] <- creek
      all_culverts[counter,4] <- mean(this_culvert)
      counter= counter+1
    }
  }
}

culvert_effect(m5_clarkii, 12, "Padden") %>% hist()

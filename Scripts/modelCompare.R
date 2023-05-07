one <- brm(bf(log(meandnaconcflow) ~ 
                       s(time_idx, bs="cc") + 
                       #s(time_idx, by=station, m=1, bs="cc")+  #main effect, station
                       s(station, bs="re") +  #random effect, station
                       #s(time_idx, by=creek, m=1, bs="cc")+  #main effect, creek
                       s(creek, bs="re")),  #random effect, creek
                  data = d %>% filter(species == "Oncorhynchus kisutch"), family = gaussian(), cores = 4, seed = 17,
                  iter = 1000, warmup = 500, thin = 10, refresh = 100,
                  control = list(adapt_delta = 0.99),
                  silent = 0
                  )
two <- brm(bf(log(meandnaconcflow) ~ 
                s(time_idx, bs="cc") + 
                s(time_idx, by=station, m=1, bs="cc")+  #main effect, station
                s(station, bs="re") +  #random effect, station
                s(time_idx, by=creek, m=1, bs="cc")+  #main effect, creek
                s(creek, bs="re")),  #random effect, creek
           data = d %>% filter(species == "Oncorhynchus kisutch"), family = gaussian(), cores = 4, seed = 17,
           iter = 1000, warmup = 500, thin = 10, refresh = 100,
           control = list(adapt_delta = 0.99),
           silent = 0
           )
three <- brm(bf(log(meandnaconcflow) ~ 
                s(time_idx, bs="cc") + 
                s(time_idx, by=station, m=1, bs="cc")+  #main effect, station
                s(station, bs="re") +  #random effect, station
                s(time_idx, by=creek, m=1, bs="cc")),  #main effect, creek
                # s(creek, bs="re")),  #random effect, creek
           data = d %>% filter(species == "Oncorhynchus kisutch"), family = gaussian(), cores = 4, seed = 17,
           iter = 1000, warmup = 500, thin = 10, refresh = 100,
           control = list(adapt_delta = 0.99),
           silent = 0
)
four <- brm(bf(log(meandnaconcflow) ~ 
                  station +
                  s(time_idx, bs="cc") + 
                  #s(time_idx, by=station, m=1, bs="cc")+  #main effect, station
                  #s(station, bs="re") +  #random effect, station
                  s(time_idx, by=creek, m=1, bs="cc") +  #main effect, creek
                  s(creek, bs="re")),  #random effect, creek
             data = d %>% filter(species == "Oncorhynchus kisutch"), family = gaussian(), cores = 4, seed = 17,
             iter = 1000, warmup = 500, thin = 10, refresh = 100,
             control = list(adapt_delta = 0.99),
             silent = 0)

five <- brm(bf(log(meandnaconcflow) ~ 
                 s(time_idx, bs="cc") + 
                 #s(time_idx, by=station, m=1, bs="cc")+  #main effect, station
                 s(station, bs="re") +  #random effect, station
                 s(time_idx, by=creek, m=1, bs="cc") +  #main effect, creek
                 s(creek, bs="re")),  #random effect, creek
            data = d %>% filter(species == "Oncorhynchus kisutch"), family = gaussian(), cores = 4, seed = 17,
            iter = 1000, warmup = 500, thin = 10, refresh = 100,
            control = list(adapt_delta = 0.99),
            silent = 0)

six <- brm(bf(log(meandnaconcflow) ~ 
                 s(time_idx, station, bs="fs") + 
                 #s(time_idx, by=station, m=1, bs="cc")+  #main effect, station
                 #s(station, bs="re") +  #random effect, station
                 s(time_idx, by=creek, m=1, bs="cc") +  #main effect, creek
                 s(creek, bs="re")),  #random effect, creek
            data = d %>% filter(species == "Oncorhynchus kisutch"), family = gaussian(), cores = 4, seed = 17,
            iter = 1000, warmup = 500, thin = 10, refresh = 100,
            control = list(adapt_delta = 0.99),
            silent = 0)

seven <- brm(bf(log(meandnaconcflow) ~ 
                s(time_idx, station, bs="fs") + 
                s(creek, bs="fs")),
                #s(time_idx, by=station, m=1, bs="cc")+  #main effect, station
                #s(station, bs="re") +  #random effect, station
                #s(time_idx, by=creek, m=1, bs="cc") +  #main effect, creek
                #s(creek, bs="re")),  #random effect, creek
           data = d %>% filter(species == "Oncorhynchus kisutch"), family = gaussian(), cores = 4, seed = 17,
           iter = 1000, warmup = 500, thin = 10, refresh = 100,
           control = list(adapt_delta = 0.99),
           silent = 0)

loo_compare(loo(one), loo(two), loo(three), loo(four), loo(five), loo(six), loo(seven))

plot(five, variable = "station")
plot(four)

d %>%
  filter(species == "Oncorhynchus kisutch") %>% 
  add_predicted_draws(six) %>%
  ggplot(aes(y = log(meandnaconcflow), x= time_idx, color = station)) +
  stat_lineribbon(aes(y = .prediction)) +
  geom_point(data = d %>% filter(species == "Oncorhynchus kisutch")) +
  scale_fill_brewer(palette = "Greys") +
  facet_grid(.~creek)

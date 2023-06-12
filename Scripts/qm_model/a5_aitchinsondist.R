# Load packages
library(here)
library(tidyverse)
library(unikn)
library(vegan)

bayes_out_prt <- readRDS(here("Output","metabarcoding", "qm_model", "20230406_bayesout_3spp_1Prt.RDS"))
bayes_out_brn <- readRDS(here("Output","metabarcoding", "qm_model", "20230406_bayesout_4spp_2Brn.RDS"))
bayes_out_chk <- readRDS(here("Output","metabarcoding", "qm_model", "20230406_bayesout_4spp_3Chk.RDS"))
bayes_out_pad <- readRDS(here("Output","metabarcoding", "qm_model", "20230406_bayesout_4spp_4Pad.RDS"))
bayes_out_sqm <- readRDS(here("Output","metabarcoding", "qm_model", "20230406_bayesout_4spp_5Sqm.RDS"))

####################################################################
# Set color palette
####################################################################

# Set color palette so don't change when plot
o_i_colors <- c(#rgb(230, 159,   0, maxColorValue = 255),  # orange
  rgb( 86, 180, 233, maxColorValue = 255),  # skyblue
  rgb(  0, 158, 115, maxColorValue = 255),  # green
  #rgb(240, 228,  66, maxColorValue = 255),  # yellow
  rgb(  0, 114, 178, maxColorValue = 255),  # blue
  rgb(204, 121, 167, maxColorValue = 255), # purple
  rgb(211,211,211, maxColorValue = 255) # GREY FOR NOT SAMPLED
)
pal_okabe_ito <- newpal(col = o_i_colors,
                        names = c("Oncorhynchus clarkii","Oncorhynchus kisutch","Oncorhynchus mykiss","Oncorhynchus nerka", "Not Sampled"))

####################################################################
# Calculate proportions of salmonids after QM correction
####################################################################

prtpost <- bayes_out_prt$Bayes_estimates %>%
  rownames_to_column("sample") %>%
  pivot_longer(-sample, names_to = "species") 

brnpost <- bayes_out_brn$Bayes_estimates %>%
  rownames_to_column("sample") %>%
  pivot_longer(-sample, names_to = "species") 

chkpost <- bayes_out_chk$Bayes_estimates %>%
  rownames_to_column("sample") %>%
  pivot_longer(-sample, names_to = "species") 

padpost <- bayes_out_pad$Bayes_estimates %>%
  rownames_to_column("sample") %>%
  pivot_longer(-sample, names_to = "species") 

sqmpost <- bayes_out_sqm$Bayes_estimates %>%
  rownames_to_column("sample") %>%
  pivot_longer(-sample, names_to = "species") 

allpost <- rbind(prtpost,brnpost,chkpost,padpost,sqmpost)

postplot <- allpost %>%
  separate(col = sample, into = c("time", "creek", "station", "biol"), remove = FALSE) %>% 
  mutate(station = case_when(station == 1 ~ "Down",
                             station == 2 ~ "Up", 
                             station == 3 ~ "Up5")) %>%
  mutate(creek = case_when(creek == "1Prt" ~ "Portage",
                           creek == "2Brn" ~ "Barnes",
                           creek == "3Chk" ~ "Chuckanut",
                           creek == "4Pad" ~ "Padden",
                           creek == "5Sqm" ~ "Squalicum",
                           TRUE ~ creek)) %>% 
  unite(creekstntimesp, c(creek,station,time,species), remove=FALSE) %>% 
  group_by(creekstntimesp) %>% 
  mutate(avgprop = mean(value)) %>%
  #mutate(sumprop = sum(value)) %>% 
  #mutate(avgprop = value/sumprop) %>%
  mutate(facetorder = factor(creek, levels=c('Padden','Portage','Chuckanut','Squalicum', 'Barnes'))) %>% 
  filter(value > 0.001) %>%
  separate(time, into = c("month","year"), sep = 2, remove=FALSE) %>% 
  unite(newtime, c(year,month), sep="-", remove=FALSE) %>% 
  unite(creekstntime, c(creek,station,time), remove=FALSE) %>%  
  ungroup() %>% 
  select(c(creekstntime,species,avgprop)) %>% 
  distinct()


####################################################################
# Calculate Aitchison distance
####################################################################

paddenfordist <- postplot %>% 
  pivot_wider(names_from=species, values_from = avgprop) %>%
  filter(str_detect(creekstntime, "Padden")) %>% 
  column_to_rownames(var="creekstntime") %>% 
  mutate(`Oncorhynchus mykiss` = replace_na(`Oncorhynchus mykiss`, 0)) %>% 
  mutate(`Oncorhynchus clarkii` = replace_na(`Oncorhynchus clarkii`, 0)) %>% 
  mutate(`Oncorhynchus nerka` = replace_na(`Oncorhynchus nerka`, 0)) %>% 
  mutate(`Oncorhynchus kisutch` = replace_na(`Oncorhynchus kisutch`, 0)) 

aitch.dist <- vegdist(paddenfordist, method="aitchison", na.rm=TRUE, pseudocount = 1e-5)
aitch.dist.mat <- as.matrix(aitch.dist)
aitchtoplot <- as_tibble(aitch.dist.mat) 
aitchtoplot$sample2 <- row.names(aitch.dist.mat)
  
aitchlong <- aitchtoplot %>% 
  pivot_longer(!sample2, names_to="sample1", values_to="aitchdist")

dn_up11_dist <- aitchlong %>% 
  filter(str_detect(sample1, "Down")) %>% 
  filter(str_detect(sample2, "Up")) %>% 
  filter(!str_detect(sample2, "Up5")) %>% 
  separate(sample1, into=c("creek1","stn1","time1")) %>% 
  separate(sample2, into=c("creek2","stn2","time2")) %>% 
  select(-c(creek1,creek2)) %>% 
  filter(time1==time2) %>% 
  select(c(time1,aitchdist)) %>% 
  rename(dnup11dist = aitchdist)

up11_up5_dist <- aitchlong %>% 
  filter(str_detect(sample1, "Up5")) %>% 
  filter(str_detect(sample2, "Up")) %>% 
  filter(!str_detect(sample2, "Up5")) %>% 
  separate(sample1, into=c("creek1","stn1","time1")) %>% 
  separate(sample2, into=c("creek2","stn2","time2")) %>% 
  select(-c(creek1,creek2)) %>% 
  filter(time1==time2) %>% 
  select(c(time1,aitchdist)) %>% 
  rename(up11up5dist = aitchdist)

dn_up5_dist <- aitchlong %>% 
  filter(str_detect(sample1, "Down")) %>% 
  filter(str_detect(sample2, "Up5")) %>% 
  separate(sample1, into=c("creek1","stn1","time1")) %>% 
  separate(sample2, into=c("creek2","stn2","time2")) %>% 
  select(-c(creek1,creek2)) %>% 
  filter(time1==time2) %>% 
  select(c(time1,aitchdist)) %>% 
  rename(dnup5dist = aitchdist)

padaitchdistplot <- dn_up11_dist %>% 
  left_join(up11_up5_dist, by="time1") %>% 
  left_join(dn_up5_dist, by="time1") %>% 
  mutate(time1 = case_when(time1 == "0321" ~ "21-03",
                          time1 == "0421" ~ "21-04",
                          time1 == "0521" ~ "21-05",
                          time1 == "0621" ~ "21-06",
                          time1 == "0721" ~ "21-07",
                          time1 == "0821" ~ "21-08",
                          time1 == "0921" ~ "21-09",
                          time1 == "1021" ~ "21-10",
                          time1 == "1121" ~ "21-11",
                          time1 == "1221" ~ "21-12",
                          time1 == "0122" ~ "22-01",
                          time1 == "0222" ~ "22-02",
                          time1 == "0322" ~ "22-03",
                          time1 == "0422" ~ "22-04",
                          time1 == "0522" ~ "22-05",
                          time1 == "0622" ~ "22-06",
                          time1 == "0722" ~ "22-07",
                          time1 == "0822" ~ "22-08",
                          time1 == "1222" ~ "22-12",
                          TRUE ~ time1)) %>% 
  add_row(time1="22-09",dnup11dist=NA,up11up5dist=NA,dnup5dist=NA) %>% 
  add_row(time1="22-10",dnup11dist=NA,up11up5dist=NA,dnup5dist=NA) %>% 
  add_row(time1="22-11",dnup11dist=NA,up11up5dist=NA,dnup5dist=NA)

padaitchdistplotlong <- padaitchdistplot %>% 
  arrange(time1) %>% 
  pivot_longer(!time1, names_to="pairwise_sites") %>% 
  rename(distvalue = value) %>% 
  mutate(consttype=c(rep(1,times=3*5),rep(2,times=3*3), #1 is pre-SR11, 2 is during SR-11
                     rep(3,times=3*8), rep(4,times=3*3), rep(5,times=3*3))) # 3 is pre-I5, 4 is during I5, 5 is post I5

ggplot(padaitchdistplotlong %>% filter(!str_detect(pairwise_sites,"dnup5")), aes(x=time1, y=distvalue, color=pairwise_sites)) +
  geom_point(size=2, pch=16) + 
  geom_line() +
  #geom_smooth(method = "loess", na.rm=TRUE, level=0.75, aes(group=pairwise_sites, fill=pairwise_sites)) + 
  #geom_point(aes(x=time1, y=up11up5dist, color="Up11/Up5"), size=3, pch=17) +
  #geom_smooth(aes(x=time1, y=up11up5dist, color="Up11/Up5"), method = "loess") + 
  #geom_point(aes(x=time1, y=dnup5dist, color="Dn/Up5"), size=3, pch=15) +
  geom_rect(aes(xmin = 6, xmax = 8, ymin = -Inf, ymax = Inf), alpha = .01, color=NA, fill = "wheat") +
  geom_rect(aes(xmin = 17, xmax = 19, ymin = -Inf, ymax = Inf), alpha = .01, color=NA, fill = "wheat") +
  geom_rect(aes(xmin = 19, xmax = 21.5, ymin = -Inf, ymax = Inf), alpha = .01, color=NA, fill = "grey") +
  theme_bw() +
  labs(x="Date (YY-MM)", y="Aitchinson Distance") + 
  scale_x_discrete(guide = guide_axis(angle = -45)) +
  scale_fill_manual(values=c("black","blue")) +
  guides(fill="none") +
  scale_color_manual(name="Pairwise Sites", values=c("black","blue"), labels=c('Down/Up11', 'Up11/Up5'))
ggsave(here("Output","Figures","padden_aitchinson_dist.png"), units="in", width=8,height=6)


padaitchdistplotlong %>% 
  filter(str_detect(pairwise_sites,"dnup11")) %>% 
  mutate(consttype = case_when(consttype == 1 ~ 1,
                               consttype > 1 ~ 2)) %>% 
  ggplot() +
  geom_histogram(aes(x=distvalue, color=factor(consttype)), fill="white", binwidth=.5) + 
  theme_bw() +
  labs(title= "SR-11 Culvert Replacement", x="Aitchinson Distance", y="Frequency") +
  scale_color_manual(name="Pre/Post Construction", values=c("red","black"), labels=c('Pre-Construction', 'Post-Construction'))


padaitchdistplotlong %>% 
  filter(str_detect(pairwise_sites,"up11up5")) %>% 
  mutate(consttype = case_when(consttype == 5 ~ 2,
                               consttype < 5 ~ 1)) %>% 
  ggplot() +
  geom_histogram(aes(x=distvalue, color=factor(consttype)), fill="white", binwidth=.5) + 
  theme_bw() +
  labs(title= "I-5 Culvert Replacement", x="Aitchinson Distance", y="Frequency") +
  scale_color_manual(name="Pre/Post Construction", values=c("red","black"), labels=c('Pre-Construction', 'Post-Construction'))


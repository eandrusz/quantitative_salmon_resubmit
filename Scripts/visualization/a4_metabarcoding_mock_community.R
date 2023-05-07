## NGN Visualizing metabarcoding results before / after Stan model
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
library(unikn)

metabeforeqm <- read.csv(here("Output","metabarcoding", "taxa_table.csv"))
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
# Plot proportions of salmonids before QM correction
####################################################################

fixnames <- read.csv(here("Input", "samples_rename.csv"))

wsdot <- metabeforeqm %>% 
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

ngn <- metabeforeqm %>% 
  filter(! str_detect(Sample_name, "MC")) %>%
  filter(str_detect(Sample_name, "MiFish")) 

metabeforeqm <- rbind(ngn,wsdot)
  
metabeforeqmplot <- metabeforeqm %>% 
  filter(! str_detect(Sample_name, "MC")) %>% 
  filter(! str_detect(Sample_name, "Kangaroo")) %>% 
  #filter(! str_detect(Sample_name, "Up5")) %>% 
  filter(species %in% c("Oncorhynchus clarkii","Oncorhynchus kisutch","Oncorhynchus mykiss","Oncorhynchus nerka")) %>%  ### ONLY FOUR
  # group_by(Sample_name) %>% 
  # mutate(ReadDepth = sum(totalReads)) %>% 
  # mutate(propReads = totalReads/ReadDepth) %>% 
  separate(Sample_name, c("marker","time","creek","station","bio", "tech")) %>% 
  dplyr::select(-marker) %>% 
  # mutate(station = case_when(station == "Up11" ~ "Up",
  #                            TRUE ~ station)) %>%
  mutate(station = case_when(station == "Dn" ~ "Down",
                             TRUE ~ station)) %>%
  mutate(creek = case_when(creek == "1Prt" ~ "Portage",
                           creek == "2Brn" ~ "Barnes",
                           creek == "3Chk" ~ "Chuckanut",
                           creek == "4Pad" ~ "Padden",
                           creek == "5Sqm" ~ "Squalicum",
                           TRUE ~ creek))  %>% 
  separate(time, into = c("month","year"), sep = 2, remove=FALSE) %>% 
  unite(newtime, c(year,month), sep="-", remove=FALSE) %>% 
  filter(tech==1) %>%
  unite(creekstntime, c(creek,station,time), remove=FALSE) %>% 
  group_by(creekstntime) %>% 
  mutate(ReadDepth = sum(totalReads)) %>% 
  mutate(propReads = totalReads/ReadDepth) %>%
  mutate(facetorder = factor(creek, levels=c('Padden','Portage','Chuckanut','Squalicum', 'Barnes'))) %>% 
  unite(creektimestation, c("creek","newtime","station"), remove=FALSE)

### ADD when we didn't sample - to distinguish from when there were no salmonid reads (really for Barnes)
alltimes <- c("21-03", "21-04", "21-05", "21-06", "21-07", "21-08", "21-09", "21-10", "21-11", "21-12", 
              "22-01", "22-02", "22-03", "22-04", "22-05", "22-06", "22-07", "22-08", "22-12")
allcreeks <- unique(metabeforeqmplot$creek)
allstations <- unique(metabeforeqmplot$station)
notsampledtimes <- c("22-09", "22-10", "22-11")
notsampledtimes2 <- c("22-03", "22-04", "22-05", "22-06", "22-07", "22-08", "22-12")

notsamp <- expand.grid(notsampledtimes2, allcreeks, allstations)  
notsamp <- notsamp %>% 
  mutate(species="Not Sampled") %>% 
  rename(newtime=Var1, creek=Var2, station=Var3) %>% 
  unite(creektimestation, c("creek","newtime","station"), remove=FALSE) %>% 
  filter(! creektimestation %in% metabeforeqmplot$creektimestation) %>% 
  mutate(propReads=1) %>%
  mutate(facetorder = factor(creek, levels=c('Padden','Portage','Chuckanut','Squalicum', 'Barnes'))) %>% 
  add_row(creektimestation="Squalicum_21-09_Up", newtime = "21-09", creek="Squalicum", station="Up", species="Not Sampled", propReads=1, facetorder="Squalicum")

metabeforeqmplot_nopad <- metabeforeqmplot %>% 
  filter(creek != "Padden") %>%
  filter(station != "Up5") %>%
  select(c(newtime,propReads,species, station, facetorder)) 

metabeforeqmplot_nopad_nosample <- rbind(metabeforeqmplot_nopad, (notsamp %>% filter(creek != "Padden") %>% filter(station != "Up11") %>% filter(station != "Up5")))

metabeforeqmplot_nopad_nosample %>% 
ggplot(aes(x=newtime, y=propReads, fill=species, color=species)) +
  geom_col() +
  facet_grid(~station~facetorder) +
  scale_fill_manual(values = pal_okabe_ito) +
  scale_color_manual(values = pal_okabe_ito) +
  guides(color= 'none') + 
  labs(x="Date (YY-MM)", y= "Proportion of DNA before QM correction", title="Control Creeks") + 
  theme_bw() + 
  #geom_vline(data=filter(postplot, facetorder=="Padden"), aes(xintercept=6.5), color="black", linetype=2) + 
  scale_x_discrete(guide = guide_axis(angle = -45))

ggsave(here("Output","SupplementalFigures","proportions_before_qm_nopadden.png"), units="in", width=14, height=5)

metabeforeqmplot %>% 
  filter(creek == "Padden") %>% 
  ggplot(aes(x=newtime, y=propReads, fill=species, color=species)) +
  geom_col() +
  facet_grid(rows="station") +
  scale_fill_manual(values = pal_okabe_ito) +
  scale_color_manual(values = pal_okabe_ito) +
  guides(color= 'none') + 
  labs(x="Date (YY-MM)", y= "Proportion of DNA before QM correction", title="Padden Creek") + 
  theme_bw() + 
  #geom_vline(data=filter(postplot, facetorder=="Padden"), aes(xintercept=6.5), color="black", linetype=2) + 
  scale_x_discrete(guide = guide_axis(angle = -45))

ggsave(here("Output","SupplementalFigures","proportions_before_qm_padden.png"), units="in", width=8, height=5)


####################################################################
# Plot proportions of salmonids after QM correction
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
  unite(creekstntime, c(creek,station,time), remove=FALSE) %>% 
  group_by(creekstntime) %>% 
  mutate(sumprop = sum(value)) %>% 
  mutate(avgprop = value/sumprop) %>%
  mutate(facetorder = factor(creek, levels=c('Padden','Portage','Chuckanut','Squalicum', 'Barnes'))) %>% 
  filter(value > 0.001) %>%
  separate(time, into = c("month","year"), sep = 2, remove=FALSE) %>% 
  unite(newtime, c(year,month), sep="-", remove=FALSE)

postplot_nopad <- postplot %>% 
  filter(creek != "Padden") %>%
  select(c(newtime, avgprop, species, station, facetorder)) 

postplot_nopad_nosample <- rbind(postplot_nopad, (notsamp %>% rename(avgprop = propReads) %>% filter(creek != "Padden") %>% filter(station != "Up11") %>% filter(station != "Up5")))


postplot_nopad_nosample %>% 
  ggplot(aes(x = newtime, fill = species, color = species, y = avgprop)) +
  geom_col() +
  facet_grid(~station ~facetorder) +
  scale_fill_manual(values = pal_okabe_ito) + 
  scale_color_manual(values = pal_okabe_ito) +
  guides(color= 'none') +
  ylab("Proportion of DNA after QM correction") +
  labs(x="Date (YY-MM)", fill="Species") +
  ggtitle("Control Creeks") +
  scale_x_discrete(guide = guide_axis(angle = -45)) +
  #geom_vline(data=filter(postplot, facetorder=="Padden"), aes(xintercept=6.5), color="black", linetype=2) + 
  theme_bw() +
  theme(legend.position="bottom")

ggsave(here("Output","Figures","proportions_after_qm_nopadden.png"), units="in", width=14, height=5)

paddenpostplot <- postplot %>% 
  filter(creek == "Padden")  

paddenpostplot %>% 
  ggplot(aes(x = newtime, fill = species, color = species, y = avgprop)) +
  geom_col() +
  facet_grid(rows="station") +
  scale_fill_manual(values = pal_okabe_ito) + 
  scale_color_manual(values = pal_okabe_ito) +
  guides(color= 'none') +
  ylab("Proportion of DNA after QM correction") +
  labs(x="Date (YY-MM)", fill="Species") +
  ggtitle("Padden Creek") +
  scale_x_discrete(guide = guide_axis(angle = -45)) +
  geom_vline( aes(xintercept=6.5), color="black", linetype=2) + 
  geom_vline( aes(xintercept=7.5), color="black", linetype=2) + 
  geom_vline(aes(xintercept=16.5), color="black", linetype=2) + 
  geom_vline( aes(xintercept=18.5), color="black", linetype=2) + 
  theme_bw()

ggsave(here("Output","Figures","proportions_after_qm_padden.png"), units="in", width=8, height=5)


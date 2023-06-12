

####################################################################
# Set up
####################################################################

## Load packages
library(tidyverse)
library(here)
library(cowplot)

# read in hash key and asv table
hash_key <- read.csv(here("Input","metabarcoding","wsdot_combined_hash.csv"))
ASVs <- read.csv(here("Input","metabarcoding","wsdot_combined_asv.csv"))
taxonomy <- read.csv(here("Output","metabarcoding","taxonomy_hash_key.csv"))
meta_metadata <- read.csv(here("Input","metabarcoding","all_mifish_metadata.csv"))
fixnames <- read.csv(here("Input","samples_rename.csv"))

wsdot <- ASVs %>% 
  filter(! str_detect(Sample_name, "MC")) %>%
  filter(! str_detect(Sample_name, "Kangaroo")) %>%
  filter(!str_detect(Sample_name, "MiFish")) %>% 
  left_join(fixnames, by="Sample_name") %>% 
  dplyr::select(-Sample_name) %>%
  separate(New_name, into=c("creek", "station","time","bio","tech")) %>% 
  mutate(creek=case_when(creek == "Pad" ~ "4Pad",
                         creek == "Prt" ~ "1Prt",
                         creek == "Sqm" ~ "5Sqm")) %>% 
  mutate(marker="MiFish") %>% 
  unite(Sample_name, c("marker", "time", "creek", "station", "bio", "tech"), sep=".") %>% 
  filter(nReads > 0) %>% 
  dplyr::select(c(Sample_name, Hash, nReads))

wsdotkangaroo <- ASVs  %>% 
  filter(!str_detect(Sample_name, "MiFish")) %>% 
  filter(str_detect(Sample_name, "Kangaroo")) %>%
  filter(nReads > 0) %>% 
  mutate(Sample_name=case_when(Sample_name == "Kangaroo" ~ "MiFish.Kangaroo.WSDOTRun23",
                               Sample_name == "Kangaroo-0522" ~ "MiFish.Kangaroo.WSDOTRun1")) %>% 
  dplyr::select(c(Sample_name, Hash, nReads))


ngn <- ASVs %>% 
  filter(! str_detect(Sample_name, "MC")) %>%
  filter(str_detect(Sample_name, "MiFish")) %>% 
  dplyr::select(c(Sample_name, Hash, nReads))

MC <- ASVs %>% 
  filter(str_detect(Sample_name, "MC")) %>%
  dplyr::select(c(Sample_name, Hash, nReads))

ASVs <- rbind(ngn, wsdot, wsdotkangaroo, MC)


finalasvtotaxa <- ASVs %>% 
  left_join(taxonomy, by="Hash") %>% 
  group_by(Sample_name) %>% 
  mutate(ReadDepth = sum(nReads)) %>% 
  mutate(TotalASVs = length(unique(Hash)))

envsamples <- finalasvtotaxa %>% 
  filter(! str_detect(Sample_name, "MC")) %>% 
  filter(! str_detect(Sample_name, "Kangaroo")) 

## plot read depths
medianRD = median(unique(envsamples$ReadDepth))
meanRD = mean(unique(envsamples$ReadDepth))
envsamples %>% 
  select(c(Sample_name, ReadDepth)) %>% 
  distinct() %>% 
  separate(Sample_name, into=c("marker","time","creek","station","bio", "tech")) %>% 
  mutate(creek = case_when(creek == "1Prt" ~ "Portage",
                              creek == "2Brn" ~ "Barnes",
                              creek == "3Chk" ~ "Chuckanut",
                              creek == "4Pad" ~ "Padden",
                              creek == "5Sqm" ~ "Squalicum",
                              TRUE ~ creek)) %>% 
  ggplot(aes(x=ReadDepth, color=creek, fill=creek)) +
  geom_histogram(position="stack", alpha=0.5, bins=40) +
  geom_vline(xintercept=medianRD, linetype="dashed") +
  labs(y="Number of Samples", x="Read Depth", fill="Creek", title="Read Depth of Environmental Samples") +
  guides(color="none") +
  theme_bw()
ggsave(filename = here("Output", "SupplementalFigures","read_depths.png"))

# annotation
envsamplesannotated <- envsamples %>% 
  filter(! is.na(species)) 

mcsamples <- finalasvtotaxa %>% 
  filter(str_detect(Sample_name, "MC")) 

mcsamplesannotated <- mcsamples %>% 
  filter(! is.na(species)) 

persample_annotation <- finalasvtotaxa %>% 
  filter(is.na(species)) %>% 
  group_by(Sample_name) %>% 
  mutate(NAreads = sum(nReads)) %>% 
  mutate(NAasvs = length(unique(Hash))) %>% 
  select(c(Sample_name, ReadDepth, TotalASVs, NAreads, NAasvs)) %>%
  distinct() %>% 
  mutate(Areads = ReadDepth-NAreads) %>% 
  mutate(Aasvs = TotalASVs-NAasvs) %>% 
  mutate(percAreads = Areads/ReadDepth*100) %>% 
  mutate(percAasvs = Aasvs/TotalASVs*100)
   
 
env_annotation <- persample_annotation %>% 
  filter(! str_detect(Sample_name, "MC")) %>% 
  filter(! str_detect(Sample_name, "Kangaroo")) %>% 
  #filter(! str_detect(Sample_name, "Up5")) %>% 
  separate(Sample_name, into=c("marker","time","creek","station","bio","tech"), remove=FALSE) %>% 
  # mutate(., station = case_when(station == "Up11" ~ "Up",
  #                               TRUE ~ station)) %>% 
  mutate(., creek = case_when(creek == "1Prt" ~ "Portage Creek",
                              creek == "2Brn" ~ "Barnes Creek",
                              creek == "3Chk" ~ "Chuckanut Creek",
                              creek == "4Pad" ~ "Padden Creek",
                              creek == "5Sqm" ~ "Squalicum Creek",
                              TRUE ~ creek)) 

env_annotation_stats <- env_annotation %>% 
  ungroup() %>% 
  select(c(ReadDepth, TotalASVs, percAreads, percAasvs)) %>% 
  summarize(meanpercAreads = mean(percAreads), meanpercAasvs = mean(percAasvs),
            medianpercAreads = median(percAreads), medianpercAasvs = median(percAasvs),
            minpercreads = min(percAreads), minpercasvs = min(percAasvs),
            maxpercreads = max(percAreads), maxpercasvs = max(percAasvs),
            meanTOTALreads = mean(ReadDepth), meanTOTALasvs = mean(TotalASVs),
            minTOTALreads = min(ReadDepth), minTOTALasvs = min(TotalASVs),
            maxTOTALreads = max(ReadDepth), maxTOTALasvs = max(TotalASVs))
  
env_annotation  %>% 
  ggplot(aes(x=ReadDepth, y=percAreads, color=time)) +
  geom_point() +
  scale_x_log10(limits = c(1E3, 1E6), guide = guide_axis(angle = -45)) +
  facet_grid(~station~creek) +
  labs(x="Total number of reads", y= "Percent of reads annotated") + 
  theme_bw() 

ggsave(here("Output","SupplementalFigures","percentannotatedbytime.png"), units="in", width=10, height=4)

# meta_to_join <- meta_metadata %>% 
#   select(-Sample_name) %>% 
#   dplyr::rename(Sample_name = Sample_ID) %>% 
#   separate(Sample_name, into=c("marker","time","creek","station","bio")) %>% 
#   mutate(tech=1) %>% 
#   unite(Sample_name, c("marker","time","creek","station","bio", "tech"), sep=".")

  
salmonidsused <- envsamples %>% 
  filter(species %in% c("Oncorhynchus kisutch", "Oncorhynchus nerka", "Oncorhynchus mykiss", "Oncorhynchus clarkii")) 

enviroreads <- sum(envsamples$nReads)
enviroasvs <- length(unique(envsamples$Hash))

enviroannotatedreads <- sum(envsamplesannotated$nReads)
enviroannotatedasvs <- length(unique(envsamplesannotated$Hash))

enviropercreadsannotated <- enviroannotatedreads/enviroreads*100
enviropercasvsannotated<- enviroannotatedasvs/enviroasvs*100

totalsalmonidreads <- sum(salmonidsused$nReads)
percentsalmonidreadstotal <- totalsalmonidreads/enviroreads*100
percentsalmonidreadsannotated <- totalsalmonidreads/enviroannotatedreads*100

mcreads <- sum(mcsamples$nReads)
mcannotatedreads <- sum(mcsamplesannotated$nReads)
mcpercannotated <- mcannotatedreads/mcreads*100


salmonidspersample <- salmonidsused %>% 
  select(c(Sample_name,species, nReads)) %>% 
  group_by(Sample_name,species) %>% 
  summarize(salmonreads = sum(nReads)) %>% 
  distinct() 

salmonidreaddepth <- salmonidspersample %>% 
  group_by(Sample_name) %>% 
  summarize(totalsalmonreads = sum(salmonreads))
  
salmonstats <- salmonidspersample %>% 
  group_by(species) %>% 
  summarize(numsamples = n()) %>% 
  mutate(percentsamples = numsamples/length(unique(envsamples$Sample_name))*100)

oclarkiistats <- salmonidspersample %>% 
  group_by(Sample_name) %>% 
  mutate(foursalmonreads=sum(salmonreads)) %>% 
  mutate(percentofeach = salmonreads/foursalmonreads*100) %>% 
  filter(species == "Oncorhynchus clarkii")

overhalfclarkii <- oclarkiistats %>% filter(percentofeach >50)
percentoverhalfclarkii = nrow(overhalfclarkii)/nrow(env_annotation)*100

## plot other taxa found to demonstrate more than just salmonids were in the samples 
taxatoplot <- envsamplesannotated %>% 
  select(c(Sample_name, nReads, species)) %>% 
  separate(Sample_name, into=c("marker","time","creek","station","bio","tech")) %>% 
  select(-marker) %>% 
  unite(newsample, c(creek,station), remove=FALSE) %>% 
  # mutate(., station = case_when(station == "Up11" ~ "Up",
  #                               TRUE ~ station)) %>% 
  mutate(., creek = case_when(creek == "1Prt" ~ "Portage",
                              creek == "2Brn" ~ "Barnes",
                              creek == "3Chk" ~ "Chuckanut",
                              creek == "4Pad" ~ "Padden",
                              creek == "5Sqm" ~ "Squalicum",
                              TRUE ~ creek)) %>%
  mutate(time = case_when(time == "0321" ~ "21-03",
                          time == "0421" ~ "21-04",
                          time == "0521" ~ "21-05",
                          time == "0621" ~ "21-06",
                          time == "0721" ~ "21-07",
                          time == "0821" ~ "21-08",
                          time == "0921" ~ "21-09",
                          time == "1021" ~ "21-10",
                          time == "1121" ~ "21-11",
                          time == "1221" ~ "21-12",
                          time == "0122" ~ "22-01",
                          time == "0222" ~ "22-02",
                          time == "0322" ~ "22-03",
                          time == "0422" ~ "22-04",
                          time == "0522" ~ "22-05",
                          time == "0622" ~ "22-06",
                          time == "0722" ~ "22-07",
                          time == "0822" ~ "22-08",
                          time == "1222" ~ "22-12",
                             TRUE ~ time)) %>% 
  # mutate(station = case_when(station == "Dn" ~ "Downstream",
  #                            station == "Up" ~ "Upstream")) %>%
  mutate(facetorder = factor(creek, levels=c('Padden','Portage','Chuckanut','Squalicum', 'Barnes'))) 

classlist <- read.csv(here("Input","metabarcoding","specieslist.csv"))

taxatoplot <- taxatoplot %>% 
  #filter(!str_detect(newsample, "Up11")) %>%
  left_join(classlist, by="species") %>% 
  filter(nReads > 2 ) 

nshowsup <- taxatoplot %>% 
  group_by(species) %>% 
  summarize(n=n())

moreten <- nshowsup %>% 
  filter(n>9)

fishcon <- taxatoplot %>% 
  filter(creek != "Padden") %>% 
  filter(class == "fish") %>% 
  filter(species %in% moreten$species) %>% 
  ggplot(aes(x = time, y = species, fill = nReads)) +
  geom_tile() +
  facet_grid(~facetorder~station) +
  geom_vline(data = taxatoplot %>% filter(creek == "Chuckanut"), aes(xintercept = 12.5)) +
  geom_vline(data = taxatoplot %>% filter(creek == "Barnes"), aes(xintercept = 12.5)) +
  theme_bw() + 
  labs(x="Date (YY-MM)", y="Species",fill="Number of Reads") +
  ggtitle("Control Creeks: Fish (in more than 10 samples)") + 
  scale_x_discrete(guide = guide_axis(angle = -45))
ggsave(here("Output", "SupplementalFigures","fish_controlcreeks.png"), units="in", width=10, height=8)

fishpad <- taxatoplot %>% 
  filter(creek == "Padden") %>% 
  filter(class == "fish") %>% 
  filter(species %in% moreten$species) %>% 
  ggplot(aes(x = time, y = species, fill = nReads)) +
  geom_tile() +
  facet_grid(~facetorder~station) +
  theme_bw() + 
  labs(x="Date (YY-MM)", y="Species",fill="Number of Reads") +
  ggtitle("Padden Creek: Fish (in more than 10 samples)") + 
  scale_x_discrete(guide = guide_axis(angle = -45))
ggsave(here("Output", "SupplementalFigures","fish_padden.png"), units="in", width=12, height=5)


mamcon <- taxatoplot %>% 
  filter(creek != "Padden") %>% 
  filter(class == "mammal") %>% 
  filter(species %in% moreten$species) %>% 
  ggplot(aes(x = time, y = species, fill = nReads)) +
  geom_tile() +
  facet_grid(~facetorder~station) +
  geom_vline(data = taxatoplot %>% filter(creek == "Chuckanut"), aes(xintercept = 12.5)) +
  geom_vline(data = taxatoplot %>% filter(creek == "Barnes"), aes(xintercept = 12.5)) +
  theme_bw() + 
  labs(x="Date (YY-MM)", y="Species",fill="Number of Reads") +
  ggtitle("Control Creeks: Mammals (in more than 10 samples)") + 
  scale_x_discrete(guide = guide_axis(angle = -45))
ggsave(here("Output", "SupplementalFigures","mammals_control.png"), units="in", width=10, height=8)

mampad <- taxatoplot %>% 
  filter(creek == "Padden") %>% 
  filter(class == "mammal") %>% 
  filter(species %in% moreten$species) %>% 
  ggplot(aes(x = time, y = species, fill = nReads)) +
  geom_tile() +
  facet_grid(~facetorder~station) +
  theme_bw() + 
  labs(x="Date (YY-MM)", y="Species",fill="Number of Reads") +
  ggtitle("Padden Creek: Mammals (in more than 10 samples)") + 
  scale_x_discrete(guide = guide_axis(angle = -45))
ggsave(here("Output", "SupplementalFigures","mammals_padden.png"), units="in", width=12, height=5)

birdcon <- taxatoplot %>% 
  filter(creek != "Padden") %>% 
  filter(class == "bird") %>% 
  filter(species %in% moreten$species) %>% 
  ggplot(aes(x = time, y = species, fill = nReads)) +
  geom_tile() +
  facet_grid(~facetorder~station) +
  geom_vline(data = taxatoplot %>% filter(creek == "Chuckanut"), aes(xintercept = 12.5)) +
  geom_vline(data = taxatoplot %>% filter(creek == "Barnes"), aes(xintercept = 12.5)) +
  theme_bw() + 
  labs(x="Date (YY-MM)", y="Species",fill="Number of Reads") +
  ggtitle("Control Creeks: Birds (in more than 10 samples)") + 
  scale_x_discrete(guide = guide_axis(angle = -45))
ggsave(here("Output", "SupplementalFigures","birds_controls.png"), units="in", width=10, height=8)

birdpad <- taxatoplot %>% 
  filter(creek == "Padden") %>% 
  filter(class == "bird") %>% 
  filter(species %in% moreten$species) %>% 
  ggplot(aes(x = time, y = species, fill = nReads)) +
  geom_tile() +
  facet_grid(~facetorder~station) +
  theme_bw() + 
  labs(x="Date (YY-MM)", y="Species",fill="Number of Reads") +
  ggtitle("Padden Creek: Birds (in more than 10 samples)") + 
  scale_x_discrete(guide = guide_axis(angle = -45))
ggsave(here("Output", "SupplementalFigures","birds_padden.png"), units="in", width=12, height=5)

ampcon <- taxatoplot %>% 
  filter(creek != "Padden") %>% 
  filter(class == "amphibian") %>% 
  filter(species %in% moreten$species) %>% 
  ggplot(aes(x = time, y = species, fill = nReads)) +
  geom_tile() +
  geom_vline(data = taxatoplot %>% filter(creek == "Chuckanut"), aes(xintercept = 12.5)) +
  geom_vline(data = taxatoplot %>% filter(creek == "Barnes"), aes(xintercept = 12.5)) +
  facet_grid(~facetorder~station) +
  theme_bw() + 
  labs(x="Date (YY-MM)", y="Species",fill="Number of Reads") +
  ggtitle("Control Creeks: Amphibians (in more than 10 samples)") + 
  scale_x_discrete(guide = guide_axis(angle = -45))
ggsave(here("Output", "SupplementalFigures","amphibs_control.png"), units="in", width=10, height=8)

amppad <- taxatoplot %>% 
  filter(creek == "Padden") %>% 
  filter(class == "amphibian") %>% 
  filter(species %in% moreten$species) %>% 
  ggplot(aes(x = time, y = species, fill = nReads)) +
  geom_tile() +
  facet_grid(~facetorder~station) +
  theme_bw() + 
  labs(x="Date (YY-MM)", y="Species",fill="Number of Reads") +
  ggtitle("Padden Creek: Amphibians (in more than 10 samples)") + 
  scale_x_discrete(guide = guide_axis(angle = -45))
ggsave(here("Output", "SupplementalFigures","amphibs_padden.png"), units="in", width=10, height=8)

pad <- cowplot::plot_grid(fishpad, mampad, amppad, birdpad, nrow=2, ncol=2)
ggsave(here("Output", "SupplementalFigures","ALLTAXA_padden.png"), units="in", width=20, height=8)
con <- cowplot::plot_grid(fishcon, mamcon, ampcon, birdcon, nrow=2, ncol=2)
ggsave(here("Output", "SupplementalFigures","ALLTAXA_controls.png"), units="in", width=20, height=12)


####################################################################
# Write functions to speed up
####################################################################

## Write function to simplify output to be just IPC data 
simple_IPC_plate <- function(plate_num) {
  # first just get file open - they are listed by date so plate number is not in order
  filename = qPCRfiles[plate_num]
  filename2 = unlist(strsplit(filename, "_"))[4]
  filename3 = unlist(strsplit(filename2,"-"))[3]
  filename4 = str_sub(filename3, start=1, end =3)
  plate <- read_excel(qPCRfiles[plate_num], sheet = "Results", skip=6)
  colnames(plate)[7] <- "Ct"  
  colnames(plate)[2] <- "Sample"
  colnames(plate)[3] <- "Assay"
  plate$Plate <- as.numeric(filename4) # now right the actual plate number that it is (i.e., what Megan called it)
  
  # now only look at IPC data and Ct values 
  IPCdata <- plate %>% 
    filter(Assay == "IPC") %>% 
    select(c("Plate", "Well","Sample", "Ct")) %>% 
    filter(Ct != "Undetermined") %>% 
    mutate(Ct = as.numeric(Ct)) 
  
  return(IPCdata)
}

## Write function to read in plate, only look at IPC data, look at the IPC value in the NTC and the standards (which should not be inhibited), determine threshold for what is deemed inhibited 
determine_thresh <- function(IPCdata) {
  
  # using ntcs, find mean and 2 sds, ct threshold = mean + 2 sds
  IPCntc <- IPCdata %>% 
    filter(str_detect(Sample,"NTC")) %>% 
    filter(! is.na(Ct)) 
  IPCntc_mean <- mean(IPCntc$Ct)
  IPCntc_sd <- sd(IPCntc$Ct)
  IPCntcthresh <- IPCntc_mean + 2*IPCntc_sd
  
  # standards are also basically NTCs for this purpose because they should have NO inhibitors (they are just gBlocks)
  IPCstds <- IPCdata %>% 
    filter(str_detect(Sample,"St")) %>% 
    filter(! is.na(Ct))
  IPCstds_mean <- mean(IPCstds$Ct)
  IPCstds_sd <- sd(IPCstds$Ct)
  IPCstdsthresh <- IPCstds_mean + 2*IPCstds_sd
  
  # take whichever is lower (ntc or stds)
  thresh <- min(IPCntcthresh, IPCstdsthresh)
  return(thresh)
}

## Write a function to select which samples are in fact inhibited
inhib_samp <- function(IPCdata,thresh) {
  inhibited <- IPCdata %>% 
    filter(!str_detect(Sample,"St")) %>% 
    filter(! is.na(Ct)) %>% 
    group_by(Sample) %>% 
    mutate(meanCt = mean(Ct)) %>% 
    mutate(inhib = meanCt > thresh) %>% 
    filter(inhib == TRUE) %>% 
    select(c(Plate,Sample,meanCt,inhib)) %>% 
    distinct()
  
  return(inhibited)
}

# now write a function to take only samples in metadata that we are actually using 
still_inhib_samp <- function(plate_num, inhibited) {
  filename = qPCRfiles[plate_num]
  filename2 = unlist(strsplit(filename, "_"))[2]
  filename3 = unlist(strsplit(filename2,"-"))[3]
  filename4 = str_sub(filename3, start=1, end =3)
  
  touse <- qPCRmeta %>% 
    filter(qPCR_no == filename2) %>% 
    filter(target == "cutt") %>%  # all the inhibition is on cutthroat plates
    filter(use == 1) %>% # sample by sample, plate by plate (only use good SC plates, and only use the samples we THOUGHT were not inhibited using the 3 ct threshold)
    filter(! str_detect(sample, "St")) %>% 
    filter(! str_detect(sample, "NTC")) %>% 
    filter(! str_detect(sample, "EB")) 
  
  stillinhibited <- inhibited %>% 
    filter(Sample %in% touse$sample)
  
  return(stillinhibited)
}

####################################################################
# Write function to check if samples have high std dev in tech reps
####################################################################

bad_sd_samps <- function(plate_num, IPCdata, thresh) {
  filename = qPCRfiles[plate_num]
  filename2 = unlist(strsplit(filename, "_"))[2]
  
  tousestd <- qPCRmeta %>% 
    filter(qPCR_no == filename2) %>% 
    filter(target == "cutt") %>%  # all the inhibition is on cutthroat plates
    filter(use == 1) %>% # sample by sample, plate by plate (only use good SC plates, and only use the samples we THOUGHT were not inhibited using the 3 ct threshold)
    filter(str_detect(sample, "St")) 
  
  tousesamp <- qPCRmeta %>% 
    filter(qPCR_no == filename2) %>% 
    filter(target == "cutt") %>%  # all the inhibition is on cutthroat plates
    filter(use == 1) %>% # sample by sample, plate by plate (only use good SC plates, and only use the samples we THOUGHT were not inhibited using the 3 ct threshold)
    filter(! str_detect(sample, "St")) %>% 
    filter(! str_detect(sample, "NTC")) %>% 
    filter(! str_detect(sample, "EB")) 
  
  notinhib <- IPCdata %>% 
    filter(!str_detect(Sample,"St")) %>% 
    filter(! is.na(Ct)) %>% 
    group_by(Sample) %>% 
    mutate(meanCt = mean(Ct)) %>% 
    mutate(sdCt = sd(Ct)) %>% 
    mutate(inhib = meanCt > thresh) %>% 
    filter(inhib == FALSE) %>% 
    select(c(Plate,Sample,meanCt,sdCt,inhib)) %>% 
    distinct() %>% 
    filter(Sample %in% tousesamp$sample)
  
  std_sd_thresh <- IPCdata %>% 
    filter(str_detect(Sample,"St")) %>% 
    filter(! is.na(Ct)) %>% 
    group_by(Sample) %>% 
    mutate(meanCt = mean(Ct)) %>% 
    mutate(sdCt = sd(Ct)) %>% 
    select(c(Plate,Sample,meanCt,sdCt)) %>% 
    distinct() 
  
  std_sd_cutoff <- mean(2*std_sd_thresh$sdCt)
  
  badsd <- notinhib %>% 
    mutate(std_sd_cutoff=std_sd_cutoff) %>% 
    mutate(good_sd = sdCt < 1.5) %>% 
    filter(good_sd == FALSE)
  
  return(badsd)
}

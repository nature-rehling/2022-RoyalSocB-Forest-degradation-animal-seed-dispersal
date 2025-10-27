# Get file paths
data_dir <- dirname(file.path("01_data", "filename"))
script_dir <- dirname(file.path("02_scripts",  "filename"))

# Load packages
pckgs<-list("tidyverse","hypervolume","doParallel", "alphahull","glue")

lapply(pckgs,require,character.only=T)

# Load data
load(paste(getwd(),data_dir,"data_dryad.rda",sep="/"))

bchd_all # Contains 4 datasets:
# bchd_all[[1]][[1]] is seed deposition sites in intact forest
# bchd_all[[1]][[2]] is forest microhabitat in intact forest
# bchd_all[[2]][[1]] is seed deposition sites in degraded forest
# bchd_all[[2]][[2]] is forest microhabitat in degraded forest
# Columns:
# disperser is one of four dispersers: Misc, Sylvia_atricapilla, Turdus_merula and Turdus philomelos
# disperser_full: This contains the real disperser, since all rare dispersers are grouped as 'Misc'
# gps: Identifier for transect segment
# cc and ccz: Canopy cover at deposition site and the z-transformed equivalent
# low_cov and low_covz: Ground vegetation at deposition site and the z-transformed equivalent
# sum_vol and sum_volz: Deadwood volume at deposition site and the z-transformed equivalent


######### List of abbrevations
# lo = degraded forest
# og = intact forest
# mi = Mixed disperser (all 'rare' species combined)
# sa = Sylvia atricapilla
# tm = Turdus merula
# tp = Turdus philomelos
# hd = Habitat
# com = Community
# hv = Hypervolume
# cc = canopy cover
# gc = ground cover
# dw = deadwood volume


####################
#### Bootstrap #####
####################


#####################################
##### Set number of subsamples ######
#####################################
nrep = 1 # For manuscript n = 200 was used but it took several hours to compile.

# Create subsamples by bootstrapping

# Intact
bchd_all[[1]][[1]] %>% count(disperser)

nbs_og = 55 # Number of scats per disperser in each subsample
#nrep= 5 # Number of subsamples

###### og #####
og_com <- bchd_all[[1]][[1]]
og_hd <- bchd_all[[1]][[2]]
# bootstrap
set.seed("13021948")
og_com_bs <- map(vector(mode="list",length=nrep), ~map(og_com %>% group_split(disperser),~slice_sample(.x,n = nbs_og)))
og_hd_bs <- map(vector(mode="list",length=nrep), ~map(list(og_hd), ~slice_sample(.x,n = nbs_og)))


# Degraded
bchd_all[[2]][[1]] %>% count(disperser)


nbs_lo = 150 # Number of scats per disperser in each subsample
#nrep= 5 # Number of scats per disperser in each subsample

lo_com <- bchd_all[[2]][[1]]
lo_hd <- bchd_all[[2]][[2]]

# bootstrap
set.seed("13021948")
lo_com_bs <- map(vector(mode="list",length=nrep), ~map(lo_com %>% group_split(disperser),~slice_sample(.x,n = nbs_lo)))
lo_hd_bs <- map(vector(mode="list",length=nrep), ~map(list(lo_hd), ~slice_sample(.x,n = nbs_lo)))





####### Compiling the hypervolumes  for volume and centroid analyses ######

# Get mean bandwidth for single frugivores in degraded and intact forests 
mean_bw_species_lo <- map_dfc(lo_com_bs,~map(.x,~estimate_bandwidth(cbind(.x$ccz,.x$low_covz,.x$sum_volz)))) %>% rowMeans()
mean_bw_species_og <- map_dfc(og_com_bs,~map(.x,~estimate_bandwidth(cbind(.x$ccz,.x$low_covz,.x$sum_volz)))) %>% rowMeans()
mean_bw_species <- rowMeans(cbind(mean_bw_species_lo,mean_bw_species_og))


###############################
##### Intact HVs###########
###############################

#nrep=5 # Number of hypervolumes to be calculated. Caution: Calculating the hypervolumes takes several hours (if nrep = 200) even when run in parallel.

set.seed(13021948)

cluster<-makeCluster(detectCores() - 1) #set parameters for calculation in parallel
registerDoParallel(cluster)

Sys.time()
full_ft_list_og <-foreach(i=1:nrep,.packages=c("hypervolume","tidyverse")) %dopar% {
  #1
  #Community HV
  com_hv <- hypervolume_gaussian(rbind(og_com_bs[[i]][[1]],og_com_bs[[i]][[2]],og_com_bs[[i]][[3]],og_com_bs[[i]][[4]])[,c("ccz","low_covz","sum_volz")],
                                 name="Full_com_og",verbose=FALSE,kde.bandwidth = mean_bw_species)
  #No Misc HV
  no_mi_hv <- hypervolume_gaussian(rbind(og_com_bs[[i]][[2]],og_com_bs[[i]][[3]],og_com_bs[[i]][[4]])[,c("ccz","low_covz","sum_volz")],
                                   name="nomi_com_og",verbose=FALSE,kde.bandwidth = mean_bw_species)
    #No SA HV
  no_sa_hv <- hypervolume_gaussian(rbind(og_com_bs[[i]][[1]],og_com_bs[[i]][[3]],og_com_bs[[i]][[4]])[,c("ccz","low_covz","sum_volz")],
                                   name="nosa_com_og",verbose=FALSE,kde.bandwidth = mean_bw_species)
  #No TM HV
  no_tm_hv <- hypervolume_gaussian(rbind(og_com_bs[[i]][[1]],og_com_bs[[i]][[2]],og_com_bs[[i]][[4]])[,c("ccz","low_covz","sum_volz")],
                                   name="notm_com_og",verbose=FALSE,kde.bandwidth = mean_bw_species)

  #No TP HV
  no_tp_hv <- hypervolume_gaussian(rbind(og_com_bs[[i]][[1]],og_com_bs[[i]][[2]],og_com_bs[[i]][[3]])[,c("ccz","low_covz","sum_volz")],
                                   name="notp_com_og",verbose=FALSE,kde.bandwidth = mean_bw_species)

  #Misc HV
  mi_hv <- hypervolume_gaussian(og_com_bs[[i]][[1]][,c("ccz","low_covz","sum_volz")], name = "mi_hv_og", verbose = FALSE, kde.bandwidth = mean_bw_species)
  #SA HV
  sa_hv <- hypervolume_gaussian(og_com_bs[[i]][[2]][,c("ccz","low_covz","sum_volz")], name = "sa_hv_og", verbose = FALSE, kde.bandwidth = mean_bw_species)
  #TM HV
  tm_hv <- hypervolume_gaussian(og_com_bs[[i]][[3]][,c("ccz","low_covz","sum_volz")], name = "tm_hv_og", verbose = FALSE, kde.bandwidth = mean_bw_species)
  #TP
  tp_hv <- hypervolume_gaussian(og_com_bs[[i]][[4]][,c("ccz","low_covz","sum_volz")], name = "tp_hv_og", verbose = FALSE, kde.bandwidth = mean_bw_species)
  
  # Hyperdisp
  hd_hv <- hypervolume_gaussian(og_hd_bs[[i]][[1]][,c("ccz","low_covz","sum_volz")], name = "hd_hv_og", verbose = FALSE, kde.bandwidth = mean_bw_species)
  
  
  list(com_hv, no_mi_hv, no_sa_hv,no_tm_hv,no_tp_hv,mi_hv,sa_hv,tm_hv,tp_hv,hd_hv)

}
Sys.time()

hv_vol_centr_og <- map_dfr(full_ft_list_og,~map_dfr(.x,~data.frame(disperser = .x@Name, volume = .x@Volume, ccz_centro = get_centroid(.x)[[1]], lcz_centro = get_centroid(.x)[[2]], dwz_centro = get_centroid(.x)[[3]])))



###########################
##### Degraded HVs###########
###########################

#nrep=5 # Number of hypervolumes to be calculated. Caution: Calculating the hypervolumes takes several hours (if nrep = 200) even when run in parallel.


set.seed(13021948)

cluster<-makeCluster(detectCores()-1) # Set parameters for calculation in parallel
registerDoParallel(cluster)

Sys.time()
full_ft_list_lo <-foreach(i=1:nrep,.packages=c("hypervolume","tidyverse")) %dopar% {
  #1
  #Community HV
  com_hv <- hypervolume_gaussian(rbind(lo_com_bs[[i]][[1]],lo_com_bs[[i]][[2]],lo_com_bs[[i]][[3]],lo_com_bs[[i]][[4]])[,c("ccz","low_covz","sum_volz")],
                                 name="Full_com_lo",verbose=FALSE,kde.bandwidth = mean_bw_species)
  #No Misc HV
  no_mi_hv <- hypervolume_gaussian(rbind(lo_com_bs[[i]][[2]],lo_com_bs[[i]][[3]],lo_com_bs[[i]][[4]])[,c("ccz","low_covz","sum_volz")],
                                   name="nomi_com_lo",verbose=FALSE,kde.bandwidth = mean_bw_species)

  #No SA HV
  no_sa_hv <- hypervolume_gaussian(rbind(lo_com_bs[[i]][[1]],lo_com_bs[[i]][[3]],lo_com_bs[[i]][[4]])[,c("ccz","low_covz","sum_volz")],
                                   name="nosa_com_lo",verbose=FALSE,kde.bandwidth = mean_bw_species)

  #No TM HV
  no_tm_hv <- hypervolume_gaussian(rbind(lo_com_bs[[i]][[1]],lo_com_bs[[i]][[2]],lo_com_bs[[i]][[4]])[,c("ccz","low_covz","sum_volz")],
                                   name="notm_com_lo",verbose=FALSE,kde.bandwidth = mean_bw_species)

  #No TP HV
  no_tp_hv <- hypervolume_gaussian(rbind(lo_com_bs[[i]][[1]],lo_com_bs[[i]][[2]],lo_com_bs[[i]][[3]])[,c("ccz","low_covz","sum_volz")],
                                   name="notp_com_lo",verbose=FALSE,kde.bandwidth = mean_bw_species)

  #Misc HV
  mi_hv <- hypervolume_gaussian(lo_com_bs[[i]][[1]][,c("ccz","low_covz","sum_volz")], name = "mi_hv_lo", verbose = FALSE, kde.bandwidth = mean_bw_species)
  #SA HV
  sa_hv <- hypervolume_gaussian(lo_com_bs[[i]][[2]][,c("ccz","low_covz","sum_volz")], name = "sa_hv_lo", verbose = FALSE, kde.bandwidth = mean_bw_species)
  #TM HV
  tm_hv <- hypervolume_gaussian(lo_com_bs[[i]][[3]][,c("ccz","low_covz","sum_volz")], name = "tm_hv_lo", verbose = FALSE, kde.bandwidth = mean_bw_species)
  #TP
  tp_hv <- hypervolume_gaussian(lo_com_bs[[i]][[4]][,c("ccz","low_covz","sum_volz")], name = "tp_hv_lo", verbose = FALSE, kde.bandwidth = mean_bw_species)
  
  # Hyperdisp
  hd_hv <- hypervolume_gaussian(lo_hd_bs[[i]][[1]][,c("ccz","low_covz","sum_volz")], name = "hd_hv_lo", verbose = FALSE, kde.bandwidth = mean_bw_species)
  
  
  list(com_hv, no_mi_hv, no_sa_hv,no_tm_hv,no_tp_hv,mi_hv,sa_hv,tm_hv,tp_hv,hd_hv)

}
Sys.time()

hv_vol_centr_lo <- map_dfr(full_ft_list_lo,~map_dfr(.x,~data.frame(disperser = .x@Name, volume = .x@Volume, ccz_centro = get_centroid(.x)[[1]], lcz_centro = get_centroid(.x)[[2]], dwz_centro = get_centroid(.x)[[3]])))



############# ANALYSES OF CENTROIDS AND THE VOLUME OF BOOTSTRAPPED HYPERVOLUMES ################
dispnames <- c("Misc","Sylvia_atricapilla","Turdus_merula","Turdus_philomelos")
frugmat <- t(combn(dispnames,m=2))
nummat <- t(combn(1:4,2))
# Function for pairwise comparisons
pFun <- function(x1, x2) {
  p <- ifelse(sum(x1 < x2)/length(x1) > 0.5, 1 - sum(x1 < x2)/length(x1), sum(x1 < x2)/length(x1)) * 2
  return(p)
}




# Table S4: Pairwise comparisons of the centroids of the forest microhabitat space in the intact and degraded forest of Bialowieza forest, Eastern Poland
hv_vol_centr_og %>% filter(str_detect(disperser,"hd")) %>% summarize_at(c("ccz_centro","lcz_centro","dwz_centro"),mean) #means in intact forest (z-transformed)
hv_vol_centr_lo %>% filter(str_detect(disperser,"hd")) %>% summarize_at(c("ccz_centro","lcz_centro","dwz_centro"),mean) #means in degraded forest (z-transformed)

# Canopy cover test
pFun(hv_vol_centr_og %>% filter(str_detect(disperser,"hd")) %>% select(ccz_centro) %>% pull(),
     hv_vol_centr_lo %>% filter(str_detect(disperser,"hd")) %>% select(ccz_centro) %>% pull())
# Ground cover test
pFun(hv_vol_centr_og %>% filter(str_detect(disperser,"hd")) %>% select(lcz_centro) %>% pull(),
     hv_vol_centr_lo %>% filter(str_detect(disperser,"hd")) %>% select(lcz_centro) %>% pull())
# Deadwood volume test
pFun(hv_vol_centr_og %>% filter(str_detect(disperser,"hd")) %>% select(dwz_centro) %>% pull(),
     hv_vol_centr_lo %>% filter(str_detect(disperser,"hd")) %>% select(dwz_centro) %>% pull())



##########################



# Table S5: Pairwise comparisons of the environmental centroids of the forest microhabitat space and the deposition microhabitat space of the frugivore community in the intact and degraded forest
# Intact
hv_vol_centr_og %>% filter(str_detect(disperser,"hd")) %>% summarize_at(c("ccz_centro","lcz_centro","dwz_centro"),mean) #means of forest microhabitat space (z-transformed)
hv_vol_centr_og %>% filter(str_detect(disperser,"Full_com")) %>% summarize_at(c("ccz_centro","lcz_centro","dwz_centro"),mean) #means of deposition mh space of frugivore community (z-transformed)

# Tests
pFun(hv_vol_centr_og %>% filter(str_detect(disperser,"hd")) %>% select(ccz_centro) %>% pull(),
     hv_vol_centr_og %>% filter(str_detect(disperser,"Full_com")) %>% select(ccz_centro) %>% pull()) # CC
pFun(hv_vol_centr_og %>% filter(str_detect(disperser,"hd")) %>% select(lcz_centro) %>% pull(),
     hv_vol_centr_og %>% filter(str_detect(disperser,"Full_com")) %>% select(lcz_centro) %>% pull()) # GC
pFun(hv_vol_centr_og %>% filter(str_detect(disperser,"hd")) %>% select(dwz_centro) %>% pull(),
     hv_vol_centr_og %>% filter(str_detect(disperser,"Full_com")) %>% select(dwz_centro) %>% pull()) # DW

# Degraded
hv_vol_centr_lo %>% filter(str_detect(disperser,"hd")) %>% summarize_at(c("ccz_centro","lcz_centro","dwz_centro"),mean) #means of forest microhabitat space (z-transformed)
hv_vol_centr_lo %>% filter(str_detect(disperser,"Full_com")) %>% summarize_at(c("ccz_centro","lcz_centro","dwz_centro"),mean) #means of deposition mh space of frugivore community (z-transformed)

# Tests
pFun(hv_vol_centr_lo %>% filter(str_detect(disperser,"hd")) %>% select(ccz_centro) %>% pull(),
     hv_vol_centr_lo %>% filter(str_detect(disperser,"Full_com")) %>% select(ccz_centro) %>% pull()) # CC
pFun(hv_vol_centr_lo %>% filter(str_detect(disperser,"hd")) %>% select(lcz_centro) %>% pull(),
     hv_vol_centr_lo %>% filter(str_detect(disperser,"Full_com")) %>% select(lcz_centro) %>% pull()) # GC
pFun(hv_vol_centr_lo %>% filter(str_detect(disperser,"hd")) %>% select(dwz_centro) %>% pull(),
     hv_vol_centr_lo %>% filter(str_detect(disperser,"Full_com")) %>% select(dwz_centro) %>% pull()) # DW




#######################




# Table S6 Pairwise comparison of the volume of the deposition microhabitat space [SD?] between different combinations of frugivores in the intact and degraded forest. 
# Intact
hv_vol_centr_og %>% filter(str_detect(disperser,"mi_hv")) %>% summarize_at(c("volume"),mean) # mean volume of of deposition mh space of mi (z-transformed)
hv_vol_centr_og %>% filter(str_detect(disperser,"sa_hv")) %>% summarize_at(c("volume"),mean) # mean volume of of deposition mh space of sa (z-transformed)
hv_vol_centr_og %>% filter(str_detect(disperser,"tm_hv")) %>% summarize_at(c("volume"),mean) # mean volume of of deposition mh space of tm (z-transformed)
hv_vol_centr_og %>% filter(str_detect(disperser,"tp_hv")) %>% summarize_at(c("volume"),mean) # mean volume of of deposition mh space of tp (z-transformed)


pFun(hv_vol_centr_og %>% filter(str_detect(disperser,"mi_hv")) %>% select(volume) %>% pull(),
     hv_vol_centr_og %>% filter(str_detect(disperser,"sa_hv")) %>% select(volume) %>% pull()) # mi vs sa
pFun(hv_vol_centr_og %>% filter(str_detect(disperser,"mi_hv")) %>% select(volume) %>% pull(),
     hv_vol_centr_og %>% filter(str_detect(disperser,"tm_hv")) %>% select(volume) %>% pull()) # mi vs tm
pFun(hv_vol_centr_og %>% filter(str_detect(disperser,"mi_hv")) %>% select(volume) %>% pull(),
     hv_vol_centr_og %>% filter(str_detect(disperser,"tp_hv")) %>% select(volume) %>% pull()) # mi vs tp
pFun(hv_vol_centr_og %>% filter(str_detect(disperser,"sa_hv")) %>% select(volume) %>% pull(),
     hv_vol_centr_og %>% filter(str_detect(disperser,"tm_hv")) %>% select(volume) %>% pull()) # sa vs tm
pFun(hv_vol_centr_og %>% filter(str_detect(disperser,"sa_hv")) %>% select(volume) %>% pull(),
     hv_vol_centr_og %>% filter(str_detect(disperser,"tp_hv")) %>% select(volume) %>% pull()) # sa vs tp
pFun(hv_vol_centr_og %>% filter(str_detect(disperser,"tm_hv")) %>% select(volume) %>% pull(),
     hv_vol_centr_og %>% filter(str_detect(disperser,"tp_hv")) %>% select(volume) %>% pull()) # tm vs tp

# Degraded
hv_vol_centr_lo %>% filter(str_detect(disperser,"mi_hv")) %>% summarize_at(c("volume"),mean) # mean volume of of deposition mh space of mi (z-transformed)
hv_vol_centr_lo %>% filter(str_detect(disperser,"sa_hv")) %>% summarize_at(c("volume"),mean) # mean volume of of deposition mh space of sa (z-transformed)
hv_vol_centr_lo %>% filter(str_detect(disperser,"tm_hv")) %>% summarize_at(c("volume"),mean) # mean volume of of deposition mh space of tm (z-transformed)
hv_vol_centr_lo %>% filter(str_detect(disperser,"tp_hv")) %>% summarize_at(c("volume"),mean) # mean volume of of deposition mh space of tp (z-transformed)


pFun(hv_vol_centr_lo %>% filter(str_detect(disperser,"mi_hv")) %>% select(volume) %>% pull(),
     hv_vol_centr_lo %>% filter(str_detect(disperser,"sa_hv")) %>% select(volume) %>% pull()) # mi vs sa
pFun(hv_vol_centr_lo %>% filter(str_detect(disperser,"mi_hv")) %>% select(volume) %>% pull(),
     hv_vol_centr_lo %>% filter(str_detect(disperser,"tm_hv")) %>% select(volume) %>% pull()) # mi vs tm
pFun(hv_vol_centr_lo %>% filter(str_detect(disperser,"mi_hv")) %>% select(volume) %>% pull(),
     hv_vol_centr_lo %>% filter(str_detect(disperser,"tp_hv")) %>% select(volume) %>% pull()) # mi vs tp
pFun(hv_vol_centr_lo %>% filter(str_detect(disperser,"sa_hv")) %>% select(volume) %>% pull(),
     hv_vol_centr_lo %>% filter(str_detect(disperser,"tm_hv")) %>% select(volume) %>% pull()) # sa vs tm
pFun(hv_vol_centr_lo %>% filter(str_detect(disperser,"sa_hv")) %>% select(volume) %>% pull(),
     hv_vol_centr_lo %>% filter(str_detect(disperser,"tp_hv")) %>% select(volume) %>% pull()) # sa vs tp
pFun(hv_vol_centr_lo %>% filter(str_detect(disperser,"tm_hv")) %>% select(volume) %>% pull(),
     hv_vol_centr_lo %>% filter(str_detect(disperser,"tp_hv")) %>% select(volume) %>% pull()) # tm vs tp




######################################################################




# Table S7 Pairwise comparisons of the volume reduction [SD?] due to forest degradation between different combinations of frugivores.
vol_slope <- hv_vol_centr_lo %>% filter(!(str_detect(disperser,"hd") | str_detect(disperser,"com"))) %>% pull(volume) -
  hv_vol_centr_og %>% filter(!(str_detect(disperser,"hd") | str_detect(disperser,"com"))) %>% pull(volume)

vol_slope <- data.frame(disperser = dispnames, sl = vol_slope) #sl is the difference in deposition microhabitat volume between the two forest types: Degraded - Intact

# Means
vol_slope %>% group_by(disperser) %>% summarize_at(c("sl"),mean)

pFun(vol_slope %>% filter(str_detect(disperser,"Misc")) %>% select(sl) %>% pull(),
     vol_slope %>% filter(str_detect(disperser,"Syl")) %>% select(sl) %>% pull()) # mi vs sa
pFun(vol_slope %>% filter(str_detect(disperser,"Misc")) %>% select(sl) %>% pull(),
     vol_slope %>% filter(str_detect(disperser,"merula")) %>% select(sl) %>% pull()) # mi vs tm
pFun(vol_slope %>% filter(str_detect(disperser,"Misc")) %>% select(sl) %>% pull(),
     vol_slope %>% filter(str_detect(disperser,"phil")) %>% select(sl) %>% pull()) # mi vs tp
pFun(vol_slope %>% filter(str_detect(disperser,"Syl")) %>% select(sl) %>% pull(),
     vol_slope %>% filter(str_detect(disperser,"merula")) %>% select(sl) %>% pull()) # sa vs tm
pFun(vol_slope %>% filter(str_detect(disperser,"Syl")) %>% select(sl) %>% pull(),
     vol_slope %>% filter(str_detect(disperser,"phil")) %>% select(sl) %>% pull()) # sa vs tp
pFun(vol_slope %>% filter(str_detect(disperser,"merula")) %>% select(sl) %>% pull(),
     vol_slope %>% filter(str_detect(disperser,"phil")) %>% select(sl) %>% pull()) # tm vs tp




######################################################




# Table S8 Pairwise comparisons of the environmental centroids of the deposition microhabitat space between different combinations of frugivore species in the intact forest. 
# Intact
hv_vol_centr_og %>% filter(str_detect(disperser,"mi_hv")) %>% summarize_at(c("ccz_centro","lcz_centro","dwz_centro"),mean) # mean centroids of of deposition mh space of mi (z-transformed)
hv_vol_centr_og %>% filter(str_detect(disperser,"sa_hv")) %>% summarize_at(c("ccz_centro","lcz_centro","dwz_centro"),mean) # mean centroids of of deposition mh space of sa (z-transformed)
hv_vol_centr_og %>% filter(str_detect(disperser,"tm_hv")) %>% summarize_at(c("ccz_centro","lcz_centro","dwz_centro"),mean) # mean centroids of of deposition mh space of tm (z-transformed)
hv_vol_centr_og %>% filter(str_detect(disperser,"tp_hv")) %>% summarize_at(c("ccz_centro","lcz_centro","dwz_centro"),mean) # mean centroids of of deposition mh space of tp (z-transformed)

# Canopy cover
pFun(hv_vol_centr_og %>% filter(str_detect(disperser,"mi_hv")) %>% select(ccz_centro) %>% pull(),
     hv_vol_centr_og %>% filter(str_detect(disperser,"sa_hv")) %>% select(ccz_centro) %>% pull()) # mi vs sa
pFun(hv_vol_centr_og %>% filter(str_detect(disperser,"mi_hv")) %>% select(ccz_centro) %>% pull(),
     hv_vol_centr_og %>% filter(str_detect(disperser,"tm_hv")) %>% select(ccz_centro) %>% pull()) # mi vs tm
pFun(hv_vol_centr_og %>% filter(str_detect(disperser,"mi_hv")) %>% select(ccz_centro) %>% pull(),
     hv_vol_centr_og %>% filter(str_detect(disperser,"tp_hv")) %>% select(ccz_centro) %>% pull()) # mi vs tp
pFun(hv_vol_centr_og %>% filter(str_detect(disperser,"sa_hv")) %>% select(ccz_centro) %>% pull(),
     hv_vol_centr_og %>% filter(str_detect(disperser,"tm_hv")) %>% select(ccz_centro) %>% pull()) # sa vs tm
pFun(hv_vol_centr_og %>% filter(str_detect(disperser,"sa_hv")) %>% select(ccz_centro) %>% pull(),
     hv_vol_centr_og %>% filter(str_detect(disperser,"tp_hv")) %>% select(ccz_centro) %>% pull()) # sa vs tp
pFun(hv_vol_centr_og %>% filter(str_detect(disperser,"tm_hv")) %>% select(ccz_centro) %>% pull(),
     hv_vol_centr_og %>% filter(str_detect(disperser,"tp_hv")) %>% select(ccz_centro) %>% pull()) # tm vs tp

# Ground cover
pFun(hv_vol_centr_og %>% filter(str_detect(disperser,"mi_hv")) %>% select(lcz_centro) %>% pull(),
     hv_vol_centr_og %>% filter(str_detect(disperser,"sa_hv")) %>% select(lcz_centro) %>% pull()) # mi vs sa
pFun(hv_vol_centr_og %>% filter(str_detect(disperser,"mi_hv")) %>% select(lcz_centro) %>% pull(),
     hv_vol_centr_og %>% filter(str_detect(disperser,"tm_hv")) %>% select(lcz_centro) %>% pull()) # mi vs tm
pFun(hv_vol_centr_og %>% filter(str_detect(disperser,"mi_hv")) %>% select(lcz_centro) %>% pull(),
     hv_vol_centr_og %>% filter(str_detect(disperser,"tp_hv")) %>% select(lcz_centro) %>% pull()) # mi vs tp
pFun(hv_vol_centr_og %>% filter(str_detect(disperser,"sa_hv")) %>% select(lcz_centro) %>% pull(),
     hv_vol_centr_og %>% filter(str_detect(disperser,"tm_hv")) %>% select(lcz_centro) %>% pull()) # sa vs tm
pFun(hv_vol_centr_og %>% filter(str_detect(disperser,"sa_hv")) %>% select(lcz_centro) %>% pull(),
     hv_vol_centr_og %>% filter(str_detect(disperser,"tp_hv")) %>% select(lcz_centro) %>% pull()) # sa vs tp
pFun(hv_vol_centr_og %>% filter(str_detect(disperser,"tm_hv")) %>% select(lcz_centro) %>% pull(),
     hv_vol_centr_og %>% filter(str_detect(disperser,"tp_hv")) %>% select(lcz_centro) %>% pull()) # tm vs tp

# Deadwood
pFun(hv_vol_centr_og %>% filter(str_detect(disperser,"mi_hv")) %>% select(dwz_centro) %>% pull(),
     hv_vol_centr_og %>% filter(str_detect(disperser,"sa_hv")) %>% select(dwz_centro) %>% pull()) # mi vs sa
pFun(hv_vol_centr_og %>% filter(str_detect(disperser,"mi_hv")) %>% select(dwz_centro) %>% pull(),
     hv_vol_centr_og %>% filter(str_detect(disperser,"tm_hv")) %>% select(dwz_centro) %>% pull()) # mi vs tm
pFun(hv_vol_centr_og %>% filter(str_detect(disperser,"mi_hv")) %>% select(dwz_centro) %>% pull(),
     hv_vol_centr_og %>% filter(str_detect(disperser,"tp_hv")) %>% select(dwz_centro) %>% pull()) # mi vs tp
pFun(hv_vol_centr_og %>% filter(str_detect(disperser,"sa_hv")) %>% select(dwz_centro) %>% pull(),
     hv_vol_centr_og %>% filter(str_detect(disperser,"tm_hv")) %>% select(dwz_centro) %>% pull()) # sa vs tm
pFun(hv_vol_centr_og %>% filter(str_detect(disperser,"sa_hv")) %>% select(dwz_centro) %>% pull(),
     hv_vol_centr_og %>% filter(str_detect(disperser,"tp_hv")) %>% select(dwz_centro) %>% pull()) # sa vs tp
pFun(hv_vol_centr_og %>% filter(str_detect(disperser,"tm_hv")) %>% select(dwz_centro) %>% pull(),
     hv_vol_centr_og %>% filter(str_detect(disperser,"tp_hv")) %>% select(dwz_centro) %>% pull()) # tm vs tp




##################################




# Table S9 Pairwise comparisons of the environmental centroids of the deposition microhabitat space between different combinations of frugivore species in the degraded forest
# Degraded
hv_vol_centr_lo %>% filter(str_detect(disperser,"mi_hv")) %>% summarize_at(c("ccz_centro","lcz_centro","dwz_centro"),mean) # mean centroids of of deposition mh space of mi (z-transformed)
hv_vol_centr_lo %>% filter(str_detect(disperser,"sa_hv")) %>% summarize_at(c("ccz_centro","lcz_centro","dwz_centro"),mean) # mean centroids of of deposition mh space of sa (z-transformed)
hv_vol_centr_lo %>% filter(str_detect(disperser,"tm_hv")) %>% summarize_at(c("ccz_centro","lcz_centro","dwz_centro"),mean) # mean centroids of of deposition mh space of tm (z-transformed)
hv_vol_centr_lo %>% filter(str_detect(disperser,"tp_hv")) %>% summarize_at(c("ccz_centro","lcz_centro","dwz_centro"),mean) # mean centroids of of deposition mh space of tp (z-transformed)

# Canopy cover
pFun(hv_vol_centr_lo %>% filter(str_detect(disperser,"mi_hv")) %>% select(ccz_centro) %>% pull(),
     hv_vol_centr_lo %>% filter(str_detect(disperser,"sa_hv")) %>% select(ccz_centro) %>% pull()) # mi vs sa
pFun(hv_vol_centr_lo %>% filter(str_detect(disperser,"mi_hv")) %>% select(ccz_centro) %>% pull(),
     hv_vol_centr_lo %>% filter(str_detect(disperser,"tm_hv")) %>% select(ccz_centro) %>% pull()) # mi vs tm
pFun(hv_vol_centr_lo %>% filter(str_detect(disperser,"mi_hv")) %>% select(ccz_centro) %>% pull(),
     hv_vol_centr_lo %>% filter(str_detect(disperser,"tp_hv")) %>% select(ccz_centro) %>% pull()) # mi vs tp
pFun(hv_vol_centr_lo %>% filter(str_detect(disperser,"sa_hv")) %>% select(ccz_centro) %>% pull(),
     hv_vol_centr_lo %>% filter(str_detect(disperser,"tm_hv")) %>% select(ccz_centro) %>% pull()) # sa vs tm
pFun(hv_vol_centr_lo %>% filter(str_detect(disperser,"sa_hv")) %>% select(ccz_centro) %>% pull(),
     hv_vol_centr_lo %>% filter(str_detect(disperser,"tp_hv")) %>% select(ccz_centro) %>% pull()) # sa vs tp
pFun(hv_vol_centr_lo %>% filter(str_detect(disperser,"tm_hv")) %>% select(ccz_centro) %>% pull(),
     hv_vol_centr_lo %>% filter(str_detect(disperser,"tp_hv")) %>% select(ccz_centro) %>% pull()) # tm vs tp

# Ground cover
pFun(hv_vol_centr_lo %>% filter(str_detect(disperser,"mi_hv")) %>% select(lcz_centro) %>% pull(),
     hv_vol_centr_lo %>% filter(str_detect(disperser,"sa_hv")) %>% select(lcz_centro) %>% pull()) # mi vs sa
pFun(hv_vol_centr_lo %>% filter(str_detect(disperser,"mi_hv")) %>% select(lcz_centro) %>% pull(),
     hv_vol_centr_lo %>% filter(str_detect(disperser,"tm_hv")) %>% select(lcz_centro) %>% pull()) # mi vs tm
pFun(hv_vol_centr_lo %>% filter(str_detect(disperser,"mi_hv")) %>% select(lcz_centro) %>% pull(),
     hv_vol_centr_lo %>% filter(str_detect(disperser,"tp_hv")) %>% select(lcz_centro) %>% pull()) # mi vs tp
pFun(hv_vol_centr_lo %>% filter(str_detect(disperser,"sa_hv")) %>% select(lcz_centro) %>% pull(),
     hv_vol_centr_lo %>% filter(str_detect(disperser,"tm_hv")) %>% select(lcz_centro) %>% pull()) # sa vs tm
pFun(hv_vol_centr_lo %>% filter(str_detect(disperser,"sa_hv")) %>% select(lcz_centro) %>% pull(),
     hv_vol_centr_lo %>% filter(str_detect(disperser,"tp_hv")) %>% select(lcz_centro) %>% pull()) # sa vs tp
pFun(hv_vol_centr_lo %>% filter(str_detect(disperser,"tm_hv")) %>% select(lcz_centro) %>% pull(),
     hv_vol_centr_lo %>% filter(str_detect(disperser,"tp_hv")) %>% select(lcz_centro) %>% pull()) # tm vs tp

# Deadwood
pFun(hv_vol_centr_lo %>% filter(str_detect(disperser,"mi_hv")) %>% select(dwz_centro) %>% pull(),
     hv_vol_centr_lo %>% filter(str_detect(disperser,"sa_hv")) %>% select(dwz_centro) %>% pull()) # mi vs sa
pFun(hv_vol_centr_lo %>% filter(str_detect(disperser,"mi_hv")) %>% select(dwz_centro) %>% pull(),
     hv_vol_centr_lo %>% filter(str_detect(disperser,"tm_hv")) %>% select(dwz_centro) %>% pull()) # mi vs tm
pFun(hv_vol_centr_lo %>% filter(str_detect(disperser,"mi_hv")) %>% select(dwz_centro) %>% pull(),
     hv_vol_centr_lo %>% filter(str_detect(disperser,"tp_hv")) %>% select(dwz_centro) %>% pull()) # mi vs tp
pFun(hv_vol_centr_lo %>% filter(str_detect(disperser,"sa_hv")) %>% select(dwz_centro) %>% pull(),
     hv_vol_centr_lo %>% filter(str_detect(disperser,"tm_hv")) %>% select(dwz_centro) %>% pull()) # sa vs tm
pFun(hv_vol_centr_lo %>% filter(str_detect(disperser,"sa_hv")) %>% select(dwz_centro) %>% pull(),
     hv_vol_centr_lo %>% filter(str_detect(disperser,"tp_hv")) %>% select(dwz_centro) %>% pull()) # sa vs tp
pFun(hv_vol_centr_lo %>% filter(str_detect(disperser,"tm_hv")) %>% select(dwz_centro) %>% pull(),
     hv_vol_centr_lo %>% filter(str_detect(disperser,"tp_hv")) %>% select(dwz_centro) %>% pull()) # tm vs tp




##################################





# Table S10: Pairwise comparisons of the environmental centroids of the deposition microhabitat space and the forest microhabitat space for different combinations of frugivores and forest type.
# See means above
# Intact
# Canopy cover
pFun(hv_vol_centr_og %>% filter(str_detect(disperser,"mi_hv")) %>% select(ccz_centro) %>% pull(),
     hv_vol_centr_og %>% filter(str_detect(disperser,"hd")) %>% select(ccz_centro) %>% pull()) # mi vs forest
pFun(hv_vol_centr_og %>% filter(str_detect(disperser,"sa_hv")) %>% select(ccz_centro) %>% pull(),
     hv_vol_centr_og %>% filter(str_detect(disperser,"hd")) %>% select(ccz_centro) %>% pull()) # sa vs forest
pFun(hv_vol_centr_og %>% filter(str_detect(disperser,"tm_hv")) %>% select(ccz_centro) %>% pull(),
     hv_vol_centr_og %>% filter(str_detect(disperser,"hd")) %>% select(ccz_centro) %>% pull()) # tm vs forest
pFun(hv_vol_centr_og %>% filter(str_detect(disperser,"tp_hv")) %>% select(ccz_centro) %>% pull(),
     hv_vol_centr_og %>% filter(str_detect(disperser,"hd")) %>% select(ccz_centro) %>% pull()) # tp vs forest

# Ground cover
pFun(hv_vol_centr_og %>% filter(str_detect(disperser,"mi_hv")) %>% select(lcz_centro) %>% pull(),
     hv_vol_centr_og %>% filter(str_detect(disperser,"hd")) %>% select(lcz_centro) %>% pull()) # mi vs forest
pFun(hv_vol_centr_og %>% filter(str_detect(disperser,"sa_hv")) %>% select(lcz_centro) %>% pull(),
     hv_vol_centr_og %>% filter(str_detect(disperser,"hd")) %>% select(lcz_centro) %>% pull()) # sa vs forest
pFun(hv_vol_centr_og %>% filter(str_detect(disperser,"tm_hv")) %>% select(lcz_centro) %>% pull(),
     hv_vol_centr_og %>% filter(str_detect(disperser,"hd")) %>% select(lcz_centro) %>% pull()) # tm vs forest
pFun(hv_vol_centr_og %>% filter(str_detect(disperser,"tp_hv")) %>% select(lcz_centro) %>% pull(),
     hv_vol_centr_og %>% filter(str_detect(disperser,"hd")) %>% select(lcz_centro) %>% pull()) # tp vs forest

# Deadwood
pFun(hv_vol_centr_og %>% filter(str_detect(disperser,"mi_hv")) %>% select(dwz_centro) %>% pull(),
     hv_vol_centr_og %>% filter(str_detect(disperser,"hd")) %>% select(dwz_centro) %>% pull()) # mi vs forest
pFun(hv_vol_centr_og %>% filter(str_detect(disperser,"sa_hv")) %>% select(dwz_centro) %>% pull(),
     hv_vol_centr_og %>% filter(str_detect(disperser,"hd")) %>% select(dwz_centro) %>% pull()) # sa vs forest
pFun(hv_vol_centr_og %>% filter(str_detect(disperser,"tm_hv")) %>% select(dwz_centro) %>% pull(),
     hv_vol_centr_og %>% filter(str_detect(disperser,"hd")) %>% select(dwz_centro) %>% pull()) # tm vs forest
pFun(hv_vol_centr_og %>% filter(str_detect(disperser,"tp_hv")) %>% select(dwz_centro) %>% pull(),
     hv_vol_centr_og %>% filter(str_detect(disperser,"hd")) %>% select(dwz_centro) %>% pull()) # tp vs forest


# Degraded
# Canopy cover
pFun(hv_vol_centr_lo %>% filter(str_detect(disperser,"mi_hv")) %>% select(ccz_centro) %>% pull(),
     hv_vol_centr_lo %>% filter(str_detect(disperser,"hd")) %>% select(ccz_centro) %>% pull()) # mi vs forest
pFun(hv_vol_centr_lo %>% filter(str_detect(disperser,"sa_hv")) %>% select(ccz_centro) %>% pull(),
     hv_vol_centr_lo %>% filter(str_detect(disperser,"hd")) %>% select(ccz_centro) %>% pull()) # sa vs forest
pFun(hv_vol_centr_lo %>% filter(str_detect(disperser,"tm_hv")) %>% select(ccz_centro) %>% pull(),
     hv_vol_centr_lo %>% filter(str_detect(disperser,"hd")) %>% select(ccz_centro) %>% pull()) # tm vs forest
pFun(hv_vol_centr_lo %>% filter(str_detect(disperser,"tp_hv")) %>% select(ccz_centro) %>% pull(),
     hv_vol_centr_lo %>% filter(str_detect(disperser,"hd")) %>% select(ccz_centro) %>% pull()) # tp vs forest

# Ground cover
pFun(hv_vol_centr_lo %>% filter(str_detect(disperser,"mi_hv")) %>% select(lcz_centro) %>% pull(),
     hv_vol_centr_lo %>% filter(str_detect(disperser,"hd")) %>% select(lcz_centro) %>% pull()) # mi vs forest
pFun(hv_vol_centr_lo %>% filter(str_detect(disperser,"sa_hv")) %>% select(lcz_centro) %>% pull(),
     hv_vol_centr_lo %>% filter(str_detect(disperser,"hd")) %>% select(lcz_centro) %>% pull()) # sa vs forest
pFun(hv_vol_centr_lo %>% filter(str_detect(disperser,"tm_hv")) %>% select(lcz_centro) %>% pull(),
     hv_vol_centr_lo %>% filter(str_detect(disperser,"hd")) %>% select(lcz_centro) %>% pull()) # tm vs forest
pFun(hv_vol_centr_lo %>% filter(str_detect(disperser,"tp_hv")) %>% select(lcz_centro) %>% pull(),
     hv_vol_centr_lo %>% filter(str_detect(disperser,"hd")) %>% select(lcz_centro) %>% pull()) # tp vs forest

# Deadwood
pFun(hv_vol_centr_lo %>% filter(str_detect(disperser,"mi_hv")) %>% select(dwz_centro) %>% pull(),
     hv_vol_centr_lo %>% filter(str_detect(disperser,"hd")) %>% select(dwz_centro) %>% pull()) # mi vs forest
pFun(hv_vol_centr_lo %>% filter(str_detect(disperser,"sa_hv")) %>% select(dwz_centro) %>% pull(),
     hv_vol_centr_lo %>% filter(str_detect(disperser,"hd")) %>% select(dwz_centro) %>% pull()) # sa vs forest
pFun(hv_vol_centr_lo %>% filter(str_detect(disperser,"tm_hv")) %>% select(dwz_centro) %>% pull(),
     hv_vol_centr_lo %>% filter(str_detect(disperser,"hd")) %>% select(dwz_centro) %>% pull()) # tm vs forest
pFun(hv_vol_centr_lo %>% filter(str_detect(disperser,"tp_hv")) %>% select(dwz_centro) %>% pull(),
     hv_vol_centr_lo %>% filter(str_detect(disperser,"hd")) %>% select(dwz_centro) %>% pull()) # tp vs forest



########################################




####################  Redundancy analysis ##########################
# Create new hypervolumes and hypervolume sets for intact forest:
#satu_list_og <- vector(mode="list",length=nrep)
satu_list_sh_og <- vector(mode="list",length=nrep)
satu_list_fullcom_vol_og <- vector(mode="list",length=nrep)

# loop will create 15 HVs: 
# 4 HVs with only one disperser
# 6 HVs with all possible double combinations of the four dispersers
# 4 HVs with all possible triple combinations of the four dispersers
# 1 HV  with all disperser (= full community hypervolume)

# Then hypervolume_set() is used to compute the intersection between the first 14 HVs and the full community hypervolume.

# Loop below will only save the intersection, because HVs and Sets take up quite a bit of memory (23 GB if nrep=200)
Sys.time()
for(j in 1:nrep)
{
  
  #1
  mi <- og_com_bs[[j]][[1]]
  sa <- og_com_bs[[j]][[2]]
  tm <- og_com_bs[[j]][[3]]
  tp <- og_com_bs[[j]][[4]]
  
  #2
  misa <- rbind(og_com_bs[[j]][[1]],og_com_bs[[j]][[2]])
  mitm <- rbind(og_com_bs[[j]][[1]],og_com_bs[[j]][[3]])
  mitp <- rbind(og_com_bs[[j]][[1]],og_com_bs[[j]][[4]])
  satm <- rbind(og_com_bs[[j]][[2]],og_com_bs[[j]][[3]])
  satp <- rbind(og_com_bs[[j]][[2]],og_com_bs[[j]][[4]])
  tmtp <- rbind(og_com_bs[[j]][[3]],og_com_bs[[j]][[4]])
  
  #3
  misatm <- rbind(og_com_bs[[j]][[1]],og_com_bs[[j]][[2]],og_com_bs[[j]][[3]])
  misatp <- rbind(og_com_bs[[j]][[1]],og_com_bs[[j]][[2]],og_com_bs[[j]][[4]])
  mitmtp <- rbind(og_com_bs[[j]][[1]],og_com_bs[[j]][[3]],og_com_bs[[j]][[4]])
  satmtp <- rbind(og_com_bs[[j]][[2]],og_com_bs[[j]][[3]],og_com_bs[[j]][[4]])
  
  # 4
  misatmtp <- rbind(og_com_bs[[j]][[1]],og_com_bs[[j]][[2]],og_com_bs[[j]][[3]],og_com_bs[[j]][[4]])
  
  disp_list <- list(mi,sa,tm,tp,misa,mitm,mitp,satm,satp,tmtp,misatm,misatp,mitmtp,satmtp,misatmtp)
  
  names<-c("mi","sa","tm","tp","misa","mitm","mitp","satm","satp","tmtp","misatm","misatp","mitmtp","satmtp","misatmtp")
  
  
  
  
  redundancy_hv_combs_og <-foreach(i=1:15,.packages=c("hypervolume","tidyverse")) %dopar% {
    hypervolume_gaussian(disp_list[[i]][,c("ccz","low_covz","sum_volz")],
                         name=names[i],verbose=FALSE,kde.bandwidth = mean_bw_species)
  }
  
  setlist <- foreach(i=1:14,.packages=c("hypervolume")) %dopar% {
    hypervolume_set(redundancy_hv_combs_og[[15]],redundancy_hv_combs_og[[i]],check.memory=F,verbose= F)
  }
  
  
  intdf <- data.frame(y=unlist(map(setlist,~.x[[3]]@Volume)),x=c(1,1,1,1,2,2,2,2,2,2,3,3,3,3),n=j)
  
  #satu_list_og[[j]] <- list(disp_list,redundancy_hv_combs_og,setlist,intdf)
  satu_list_sh_og[[j]] <- intdf
  satu_list_fullcom_vol_og[[j]] <- redundancy_hv_combs_og[[15]]@Volume
  print(j)
  
}
Sys.time()

intersec_og <- map2_df(satu_list_sh_og, satu_list_fullcom_vol_og, ~data.frame(disperser=names[1:14], intersect=.x$y, nrep=.x$n, full_com_vol = .y))


# Create new hypervolumes and hypervolume sets for degraded forest:
#satu_list_lo <- vector(mode="list",length=nrep)
satu_list_sh_lo <- vector(mode="list",length=nrep)
satu_list_fullcom_vol_lo <- vector(mode="list",length=nrep)


Sys.time()
for(j in 1:nrep)
{
  
  #1
  mi <- lo_com_bs[[j]][[1]]
  sa <- lo_com_bs[[j]][[2]]
  tm <- lo_com_bs[[j]][[3]]
  tp <- lo_com_bs[[j]][[4]]
  
  #2
  misa <- rbind(lo_com_bs[[j]][[1]],lo_com_bs[[j]][[2]])
  mitm <- rbind(lo_com_bs[[j]][[1]],lo_com_bs[[j]][[3]])
  mitp <- rbind(lo_com_bs[[j]][[1]],lo_com_bs[[j]][[4]])
  satm <- rbind(lo_com_bs[[j]][[2]],lo_com_bs[[j]][[3]])
  satp <- rbind(lo_com_bs[[j]][[2]],lo_com_bs[[j]][[4]])
  tmtp <- rbind(lo_com_bs[[j]][[3]],lo_com_bs[[j]][[4]])
  
  #3
  misatm <- rbind(lo_com_bs[[j]][[1]],lo_com_bs[[j]][[2]],lo_com_bs[[j]][[3]])
  misatp <- rbind(lo_com_bs[[j]][[1]],lo_com_bs[[j]][[2]],lo_com_bs[[j]][[4]])
  mitmtp <- rbind(lo_com_bs[[j]][[1]],lo_com_bs[[j]][[3]],lo_com_bs[[j]][[4]])
  satmtp <- rbind(lo_com_bs[[j]][[2]],lo_com_bs[[j]][[3]],lo_com_bs[[j]][[4]])
  
  # 4
  misatmtp <- rbind(lo_com_bs[[j]][[1]],lo_com_bs[[j]][[2]],lo_com_bs[[j]][[3]],lo_com_bs[[j]][[4]])
  
  disp_list <- list(mi,sa,tm,tp,misa,mitm,mitp,satm,satp,tmtp,misatm,misatp,mitmtp,satmtp,misatmtp)
  
  names<-c("mi","sa","tm","tp","misa","mitm","mitp","satm","satp","tmtp","misatm","misatp","mitmtp","satmtp","misatmtp")
  
  
  
  
  redundancy_hv_combs_lo <-foreach(i=1:15,.packages=c("hypervolume","tidyverse")) %dopar% {
    hypervolume_gaussian(disp_list[[i]][,c("ccz","low_covz","sum_volz")],
                         name=names[i],verbose=FALSE,kde.bandwidth = mean_bw_species)
  }
  
  setlist <- foreach(i=1:14,.packages=c("hypervolume")) %dopar% {
    hypervolume_set(redundancy_hv_combs_lo[[15]],redundancy_hv_combs_lo[[i]],check.memory=F,verbose= F)
  }
  
  
  intdf <- data.frame(y=unlist(map(setlist,~.x[[3]]@Volume)),x=c(1,1,1,1,2,2,2,2,2,2,3,3,3,3),n=j)
  
  #satu_list_lo[[j]] <- list(disp_list,redundancy_hv_combs_lo,setlist,intdf)
  satu_list_sh_lo[[j]] <- intdf
  satu_list_fullcom_vol_lo[[j]] <- redundancy_hv_combs_lo[[15]]@Volume
  print(j)
  
}
Sys.time()

intersec_lo <- map2_df(satu_list_sh_lo, satu_list_fullcom_vol_lo, ~data.frame(disperser=names[1:14], intersect=.x$y, nrep=.x$n, full_com_vol = .y))




############# Uniqueness analysis ####################
# To estimate the volume of the microhabitat space that is unique to - for example - Sylvia atricapilla we use hypervolume_set() to compare the following hypervolumes:
# HV1: Full community hypervolume.
# HV2: Hypervolume of community lacking Sylvia atricapilla.
# The unique microhabitat space of HV1 is then attributed to Sylvia atricapilla.


# Intact forest
# Create sets:

Sys.time()
full_nospe_sets_og <- foreach(i=1:nrep,.packages=c("hypervolume","tidyverse")) %dopar% {
  no_mi_set <- hypervolume_set(full_ft_list_og[[i]][[1]],full_ft_list_og[[i]][[2]],check.memory=F)
  no_sa_set <- hypervolume_set(full_ft_list_og[[i]][[1]],full_ft_list_og[[i]][[3]],check.memory=F)
  no_tm_set <- hypervolume_set(full_ft_list_og[[i]][[1]],full_ft_list_og[[i]][[4]],check.memory=F)
  no_tp_set <- hypervolume_set(full_ft_list_og[[i]][[1]],full_ft_list_og[[i]][[5]],check.memory=F)
  
  uniq_prop_og <- data.frame(disperser=dispnames,rel_uniq=c(
    no_mi_set[[5]]@Volume / no_mi_set[[1]]@Volume,
    no_sa_set[[5]]@Volume / no_sa_set[[1]]@Volume,
    no_tm_set[[5]]@Volume / no_tm_set[[1]]@Volume,
    no_tp_set[[5]]@Volume / no_tp_set[[1]]@Volume
  ))
  
  print(i)
  
  list(list(no_mi_set,no_sa_set,no_tm_set,no_tp_set),list=uniq_prop_og)
}
Sys.time()

uniq_prop_df_og <- map_df(full_nospe_sets_og,~.x[[2]])


  
# Degraded forest
# Create sets:

Sys.time()
full_nospe_sets_lo <- foreach(i=1:nrep,.packages=c("hypervolume","tidyverse")) %dopar% {
  no_mi_set <- hypervolume_set(full_ft_list_lo[[i]][[1]],full_ft_list_lo[[i]][[2]],check.memory=F)
  no_sa_set <- hypervolume_set(full_ft_list_lo[[i]][[1]],full_ft_list_lo[[i]][[3]],check.memory=F)
  no_tm_set <- hypervolume_set(full_ft_list_lo[[i]][[1]],full_ft_list_lo[[i]][[4]],check.memory=F)
  no_tp_set <- hypervolume_set(full_ft_list_lo[[i]][[1]],full_ft_list_lo[[i]][[5]],check.memory=F)
  
  uniq_prop_lo <- data.frame(disperser=dispnames,rel_uniq=c(
    no_mi_set[[5]]@Volume / no_mi_set[[1]]@Volume,
    no_sa_set[[5]]@Volume / no_sa_set[[1]]@Volume,
    no_tm_set[[5]]@Volume / no_tm_set[[1]]@Volume,
    no_tp_set[[5]]@Volume / no_tp_set[[1]]@Volume
  ))
  print(i)
  list(list(no_mi_set,no_sa_set,no_tm_set,no_tp_set),list=uniq_prop_lo)
}
Sys.time()

uniq_prop_df_lo <- map_df(full_nospe_sets_lo,~.x[[2]])


####### Table S11: Pairwise comparison of the proportion of uniqueness between dispersers in intact and degraded forest. 
# Pairwise comparisons intact
uniq_prop_og_mean <- uniq_prop_df_og %>% group_by(disperser) %>% summarize(mean_uniq = mean(rel_uniq))

pw_uniq_df_og <- data.frame(disperser1 = t(dispnames %>% combn(m = 2))[,1],
                             disperser2 = t(dispnames %>% combn(m = 2))[,2])

pw_uniq_df_og$disp1_uniq <- unlist(sapply(1:6, function(x) uniq_prop_og_mean[which(uniq_prop_og_mean$disperser == pw_uniq_df_og$disperser1[x]),"mean_uniq"]))
pw_uniq_df_og$disp2_uniq <- unlist(sapply(1:6, function(x) uniq_prop_og_mean[which(uniq_prop_og_mean$disperser == pw_uniq_df_og$disperser2[x]),"mean_uniq"]))

pw_uniq_df_og$pw_comparison <- sapply(1:6, function(x) pFun(uniq_prop_df_og[uniq_prop_df_og$disperser == pw_uniq_df_og$disperser1[x],"rel_uniq"],uniq_prop_df_og[uniq_prop_df_og$disperser == pw_uniq_df_og$disperser2[x],"rel_uniq"]))
                                      
             
#### Pairwise comparison degraded                    
uniq_prop_lo_mean <- uniq_prop_df_lo %>% group_by(disperser) %>% summarize(mean_uniq = mean(rel_uniq))

pw_uniq_df_lo <- data.frame(disperser1 = t(dispnames %>% combn(m = 2))[,1],
                            disperser2 = t(dispnames %>% combn(m = 2))[,2])

pw_uniq_df_lo$disp1_uniq <- unlist(sapply(1:6, function(x) uniq_prop_lo_mean[which(uniq_prop_lo_mean$disperser == pw_uniq_df_lo$disperser1[x]),"mean_uniq"]))
pw_uniq_df_lo$disp2_uniq <- unlist(sapply(1:6, function(x) uniq_prop_lo_mean[which(uniq_prop_lo_mean$disperser == pw_uniq_df_lo$disperser2[x]),"mean_uniq"]))

pw_uniq_df_lo$pw_comparison <- sapply(1:6, function(x) pFun(uniq_prop_df_lo[uniq_prop_df_lo$disperser == pw_uniq_df_lo$disperser1[x],"rel_uniq"],uniq_prop_df_lo[uniq_prop_df_lo$disperser == pw_uniq_df_lo$disperser2[x],"rel_uniq"]))



# FRAGEN: Brauche ich jemals bc_data, bc_data_bind, bc_hyperdisp oder nur bchd_all?

# MACHEN: nrep nur 1 mal definieren! und die anderen l?schen.
# Nomi etc -> uniqueness raus




# # Compile the full (i.e. non-bootstrapped) hypervolumes for the intact and degraded forest. These are briefly described in the results (Volume of 118 SD^3 vs 50 SD^3) and are depicted in Fig 1ABC :
# load(paste(getwd(),data_dir,"compiled_bc_data_and_hyperdisperser_loog.rda",sep="/"))
# #
# hd_df <- list(bchd_all[[1]][[2]] %>% mutate(forest_type = "og", used = ifelse(paste(ccz,low_covz,sum_volz) %in% paste(bchd_all[[1]][[1]]$ccz,bchd_all[[1]][[1]]$low_covz,bchd_all[[1]][[1]]$sum_volz),1,0)),
#               bchd_all[[2]][[2]] %>% mutate(forest_type = "lo", used = ifelse(paste(ccz,low_covz,sum_volz) %in% paste(bchd_all[[2]][[1]]$ccz,bchd_all[[2]][[1]]$low_covz,bchd_all[[2]][[1]]$sum_volz),1,0)))
# # #
# # # Create hypervolumes #
# bws_hd <- colMeans(rbind(estimate_bandwidth(bchd_all[[1]][[2]][c("ccz","low_covz","sum_volz")]),estimate_bandwidth(bchd_all[[2]][[2]][c("ccz","low_covz","sum_volz")])))
# # #
# set.seed(13021948)
# # #
# hdog_hv <- hypervolume_gaussian(bchd_all[[1]][[2]][c("ccz","low_covz","sum_volz")], name="hd_og",verbose=FALSE,kde.bandwidth = bws_hd)
# hdlo_hv <- hypervolume_gaussian(bchd_all[[2]][[2]][c("ccz","low_covz","sum_volz")], name="hd_lo",verbose=FALSE,kde.bandwidth = bws_hd)
# # #
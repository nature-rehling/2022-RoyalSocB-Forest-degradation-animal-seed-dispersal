
R script used to compile the analyses presented in the manuscript 
################################################################################
Rehling & Schlautmann et al. (2022):
Forest degradation limits the complementarity and quality of animal seed dispersal
doi: 10.1098/rspb.2022.0391
################################################################################

Parent folder contains the project
01_data contains an .rda file with the data
02_scripts contains an .r file which creates the subsamples and hypervolumes and the analysis.

For the script to run smoothly please open Hypervolume.Rproj.

The data is explained at the start of Hypervolume_analysis.R or below in this doc.

The number of subsamples can be set in line 50. For the paper we used nrep = 200 but it takes several hours to run.

data_dryad.rda will open the list 'bchd_all' in R.

bchd_all contains 4 datasets:
bchd_all[[1]][[1]] is seed deposition sites in intact forest
bchd_all[[1]][[2]] is forest microhabitat in intact forest
bchd_all[[2]][[1]] is seed deposition sites in degraded forest
bchd_all[[2]][[2]] is forest microhabitat in degraded forest
 Columns:
 disperser is one of four dispersers: Misc, Sylvia_atricapilla, Turdus_merula and Turdus philomelos
 disperser_full: This contains the real disperser, since all rare dispersers are grouped as 'Misc'
 gps: Identifier for transect segment
 cc and ccz: Canopy cover at deposition site and the z-transformed equivalent
 low_cov and low_covz: Ground vegetation at deposition site and the z-transformed equivalent
 sum_vol and sum_volz: Deadwood volume at deposition site and the z-transformed equivalent


List of abbrevations used in the script:
 lo = degraded forest
 og = intact forest
 mi = Mixed disperser (all 'rare' species combined)
 sa = Sylvia atricapilla
 tm = Turdus merula
 tp = Turdus philomelos
 hd = Habitat
 com = Community
 hv = Hypervolume
 cc = canopy cover
 gc = ground cover
 dw = deadwood volume

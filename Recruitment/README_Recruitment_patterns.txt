################################################################################
################################################################################
This package contains data from:

################################################################################
Rehling & Schlautmann et al. (2022):
Forest degradation limits the complementarity and quality of animal seed dispersal
doi: 10.1098/rspb.2022.0391

################################################################################
The package contains two data files (each in .csv format) that were used in the analyses of the above mentioned publication.

The structures of the data files are listed below with detailed descriptions of their contents.
################################################################################
Structure and description of data files that were used in the main analyses
################################################################################

1. Peak_recruitment_plants.csv
Recruitment pattern of the fleshy-fruited plant along a canopy cover gradient and a ground vegetation gradient in Bialowieza forest (Poland)
'data.frame':	332 obs. of  11 variables:
 $ X 	  : int, empty column
 $ plot   : factor, number for each recruitment plot
 $ cc     : num, forest canopy cover in [%], 
 $ z.cc   : num, z-transformed forest canopy cover
 $ bk1    : num, vegetation cover at ground level in [%]
 $ veg1   : num, z-transformed vegetation cover at ground level
 $ forest : chr, forest type; "l" = logged/degraded, "o" = old-growth/intact
 $ species: chr, plant species; "Eueu" = Euonymus europaeus, "Fral" = Frangula alnus, "Prpa" = Prunus padus, 
			"Rhca" = Rhamnus cathartica, "Rini" = Ribes nigrum, "Risp" = Ribes spicatum, 
			"Sani" = Sambucus nigra, "Soau" = Sorbus aucuparia, "Viop" = Viburnum opulus
 $ absrec2: int, number of recruited seedlings with cotyledons in year of peak recruitment
 $ denom2 : int, number of seeds remaining in the soil in the year of peak recruitment
 $ olre   : int, observation level random effect

2. Early_survival_plants.csv
'data.frame':	1751 obs. of  16 variables:
 $ gps     : chr, the lowest resolution of spatial information for the data, i.e. "plot_transect_section" -> "011_T3_50"
 $ sites   : int, study site
 $ forest  : chr, forest type; "l" = logged/degraded, "o" = old-growth/intact
 $ ind_id  : chr, ID for each individual
 $ tag_id  : chr, label ID for each individual
 $ species : chr, plant species; "Eueu" = Euonymus europaeus, "Fral" = Frangula alnus, "Prpa" = Prunus padus, 
			"Rhca" = Rhamnus cathartica, "Rini" = Ribes nigrum, "Risp" = Ribes spicatum, 
			"Sani" = Sambucus nigra, "Soau" = Sorbus aucuparia, "Viop" = Viburnum opulus
 $ h       : num, height of seedlings in [cm]
 $ year    : int, first study year
 $ hnext   : num, height of seedlings after 1 year of growth in [cm]
 $ yearnext: int, second study year
 $ surv    : int, 0 = no, 1 = yes
 $ olre    : int, observation level random effect
 $ cc      : num, forest canopy cover in [%], 
 $ z.cc    : num, z-transformed forest canopy cover
 $ bk1     : num, vegetation cover at ground level in [%]
 $ z.veg   : num, z-transformed vegetation cover at ground level
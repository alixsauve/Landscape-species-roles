# this script shows how to calculate the {z, c}-values of species observed 
# interacting within a landscape (a set of habitats)

# clear workspace
rm(list=ls(all=TRUE))

# load functions
source("Functions/mw_adjacencies.R")
source("Functions/network_role.R")
source("Functions/feature_scaling.R")

# load data
# interaction list for the mosaic of habitats in Hengistbury Head
download.file("https://ndownloader.figshare.com/files/15122768",
              "data/HH_interaction_list.csv")
HH_ListObs <- read.csv("data/HH_interaction_list.csv", stringsAsFactors = FALSE)

# normalise interaction frequencies
IntTypes <- unique(HH_ListObs$upper_guild)
HH_ListObs$freq_norm <- 0
for (it in IntTypes){
  freq_int <- HH_ListObs$freq[HH_ListObs$upper_guild == it]
  HH_ListObs$freq_norm[HH_ListObs$upper_guild == it] <- round(feature_scaling(freq_int, min(HH_ListObs$freq), max(HH_ListObs$freq)), 3) # feature scaling, 3 digits for fractional part
}

# species list for the mosaic of habitats in Hengistbury Head
# must be a data frame with at least one column with species names ('taxon'), it can be customised to your taste with more information
HH_SpList <- data.frame(taxon = unique(sort(c(HH_ListObs$lower_taxon, HH_ListObs$upper_taxon))))

# metaweb, all interaction types
a_obs <- mw_adjacencies(HH_ListObs, HH_SpList, FALSE, FALSE, TRUE) # 3D-array, dim = {Sp+Sa, Sp+Sa, H}, Sp = No. plant species, Sa = No. insect species, H = No. habitats
aw_obs <- mw_adjacencies(HH_ListObs, HH_SpList, TRUE, FALSE, TRUE) # weighted

# Among-habitat connectivity
C <- c_value(a_obs) # qualitative
Cw <- c_value(aw_obs) # quantitative

# Weighted mean of within-habitat degree/weight
Z <- z_value_landscape(a_obs) # qualitative
Zw <- z_value_landscape(aw_obs) # quantitative

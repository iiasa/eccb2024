## prepare Z data

## planning units
PU <- rast("data/other_layers/area_mask.tif")
PU <- as.numeric(PU)
# plot(PU)
writeRaster(PU, file = "data/PlanningUnits.tif")

## features
dir.create("data/SpeciesDistributions/")

sdms_raw <- list.files("data/sdms/future_rcp85_2065/", recursive = T, full.names = T)
## transform to 0 just for the harmonization of datasets, but for the prioritisation PU cells value must be =1 for the budget.

PU[PU>0] <- 0

for (i in 1:length(sdms_raw)){

  sppname <- gsub("ens-sdms_cur2005_prob_pot", "current", sdms_raw[i])
  sppname <- gsub("ens-sdms_rcp85_fut2065_prob_pot", "rcp85",sppname)
  sppname <- gsub("data/sdms/current/", "", sppname)
  sppname <- gsub("data/sdms/future_rcp85_2065/", "", sppname)

  spp1 <- rast(sdms_raw[i])

  spp <- harmonize_raster(spp1, PU, data_type = "conti")

  spp <- (spp - minmax(spp)[1]) / (minmax(spp)[2] - minmax(spp)[1])

  writeRaster(spp, file = paste0("data/SpeciesDistributions/", sppname))
  }



## protected areas

PA_1k <- rast("../../Dropbox/NaturaConnect/data/SpatialData/ProtectedAreas/Europe_1km2_PA_all.tif")

PA <- harmonize_raster(PA_1k, PU, data_type = "conti")

PA <- PA/minmax(PA)[2]
writeRaster(PA, file = "data/protectedareas.tif")


StPA_1k <- rast("../../Dropbox/NaturaConnect/data/SpatialData/ProtectedAreas/Europe_1km2_cat_I_II.tif")
StPA <- harmonize_raster(StPA_1k, PU, data_type = "conti")
StPA <- StPA/minmax(StPA)[2]
writeRaster(StPA, file = "data/protectedareas_I_II.tif")




## global human modification map
gHM_raw<- rast("DataThiago/other_layers/gHM_Europe.tif")
plot(gHM_raw)
gHM <- harmonize_raster(gHM_raw, PU, data_type = "conti")
# standardize between 0 and 1
gHM <- (gHM - minmax(gHM)[1]) / (minmax(gHM)[2] - minmax(gHM)[1])

writeRaster(gHM, file = "data/gHM.tif")


## land cover and land use

landsystems <- rast("../../Dropbox/PC/Documents/LECA/PhD/2022_crosswalk_habitats_land use/data/landSystem_withforesttype.tif")

landsystemsclasses <- read.csv("../../Dropbox/PC/Documents/LECA/PhD/2022_crosswalk_habitats_land use/FILTERING_DISTRIBUTIONS/_2022/INPUT/LandUseCodesCorrespondance.csv", sep = ";")



## reclassify urban
rcl.urban <- matrix(ncol = 3,
                    c(22, 23, 1), byrow = T)

urban.bin <- terra::classify(landsystems,rcl.urban,  include.lowest =T)

urban.bin[urban.bin>1] <- 0

urban.bin.10k <- aggregate(urban.bin, fact = 10, fun = "sum")
urban.bin.10k_rp <- harmonize_raster(urban.bin.10k, PU, data_type = "conti")

# urban.bin.10k_rp[urban.bin.10k_rp<50] <- 0

urban.bin.10k_rp <- 100*(urban.bin.10k_rp - minmax(urban.bin.10k_rp)[1]) / (minmax(urban.bin.10k_rp)[2] - minmax(urban.bin.10k_rp)[1])

writeRaster(urban.bin.10k_rp, file = "data/urban_prct.tif")

## forests
rcl.forest <- matrix(ncol = 3,
                     c(431, 432, 1), byrow= T
)
HI.forest.bin <- terra::classify(landsystems, rcl.forest, include.lowest = T)

HI.forest.bin[HI.forest.bin>1] <- 0

HI.forest.bin.10k <- aggregate(HI.forest.bin, fact = 10, fun = "sum")
HI.forest.bin.10k_rp <- harmonize_raster(HI.forest.bin.10k, PU, data_type = "conti")

# HI.forest.bin.10k_rp[HI.forest.bin.10k_rp<50] <- 0

HI.forest.bin.10k_rp <- 100*(HI.forest.bin.10k_rp - minmax(HI.forest.bin.10k_rp)[1]) / (minmax(HI.forest.bin.10k_rp)[2] - minmax(HI.forest.bin.10k_rp)[1])

writeRaster(HI.forest.bin.10k_rp, file = "data/HI_forest_prct.tif", overwrite = T)

##  NDVI
ndvi_raw<- rast("DataThiago/other_layers/ndvi_europe.tif")
plot(ndvi_raw)
ndvi <- harmonize_raster(ndvi_raw, PU, data_type = "conti")
# standardize between 0 and 1
ndvi <- (ndvi - minmax(ndvi)[1]) / (minmax(ndvi)[2] - minmax(ndvi)[1])
plot(ndvi)
writeRaster(ndvi, file = "data/ndvi.tif")



## red list for weights
# load red list data
redlist <- read.csv("../../Dropbox/NaturaConnect/data/Lists_species_habitats/EU_species_RedList_Art12_Art17.csv")
redlist$spp_name <- gsub(" ", "_", redlist$speciesname.EU.RedList)

global_redlist <- read.csv("../../Dropbox/NaturaConnect/data/Lists_species_habitats/spp-lists-fromLaetitia/IUCNRL_NatGlob/IUCN_GLOBALspeciesdata.csv")
global_redlist$spp_name <- gsub(" ", "_", global_redlist$ScientificName)

# retrieve species names from filenames
spp.list <- list.files(path = "data/SpeciesDistributions/", full.names = T, recursive = T, pattern = "tif$")
sppnames <- gsub("data/SpeciesDistributions/current/", "", spp.list)
sppnames <- gsub("_current.tif", "", sppnames)

redlist_trees_europe <- subset(redlist, redlist$spp_name %in% sppnames)

redlist_trees_global <- subset(global_redlist, global_redlist$spp_name %in% sppnames)

redlist_trees_global <- redlist_trees_global[, c("kingdom", "phylum","class","order", "family",  "spp_name", "Global", "Europe")]


## add missing species
t <- setdiff(sppnames, redlist_trees_global$spp_name)
spp_taxo <- data.frame(kingdom = rep(NA, length(t)),
                       phylum = rep(NA, length(t)),
                       class = rep(NA, length(t)),
                       order =  rep(NA, length(t)),
                       family =  rep(NA, length(t)),
                       spp_name = t)

for(i in 1:nrow(spp_taxo)){
  taxo_i <- classification(spp_taxo$spp_name[i], db = "gbif", rows= 1)
  if(!is.na(taxo_i)){
      spp_taxo$order[i] = taxo_i[[1]]$name[which(taxo_i[[1]]$rank == "order")]
      spp_taxo$class[i] = taxo_i[[1]]$name[which(taxo_i[[1]]$rank == "class")]
      spp_taxo$family[i] = taxo_i[[1]]$name[which(taxo_i[[1]]$rank == "family")]
      spp_taxo$phylum[i] = taxo_i[[1]]$name[which(taxo_i[[1]]$rank == "phylum")]
      spp_taxo$kingdom[i] = taxo_i[[1]]$name[which(taxo_i[[1]]$rank == "kingdom")]

    # spp_taxo$Genus[i] = taxo_i[[1]]$name[which(taxo_i[[1]]$rank == "genus")]
  }
}

spp_taxo$Global <- NA
spp_taxo$Europe <- NA


redlist_trees_global <- rbind(redlist_trees_global,
                              spp_taxo)

# assume NA = least concern..
redlist_trees_global$Europe[is.na(redlist_trees_global$Europe)] <- "Least Concern"
redlist_trees_global$Global[is.na(redlist_trees_global$Global)] <- "Least Concern"

unique(redlist_trees_global$Europe)
redlist_trees_global$weight <- ifelse(redlist_trees_global$Global == "Vulnerable" | redlist_trees_global$Europe == "Vulnerable", 4,
                                      ifelse(redlist_trees_global$Global == "Near Threatened" | redlist_trees_global$Europe == "Near Threatened", 2,
                                             ifelse(redlist_trees_global$Global == "Data Deficient" | redlist_trees_global$Europe == "Data Deficient", 2, 1)))


write.csv(redlist_trees_global, file = "data/species_red_list.csv", row.names = F)












#######  function to prepare spatial data given parameters ########

harmonize_raster <- function(x, PU, data_type = c("conti", "categ")){ ## continuous or categorical

  if(data_type == "conti"){
    method_rp <- "bilinear"
    if ( res(x)[1] < res(PU)[1] ){ ## aggregate to coarser resolution if needed before reprojection
      f<- (res(PU)[1])/(res(x)[1])
      x<- terra::aggregate(x, fact = f, fun = "sum", na.rm=T)
      x <- x/(f^2)
    }
  }
  ## x is the raw raster to convert
  ## PU is a raster with the right resolution, projection, extent, origin
  ## Study Area is a vector that has the outline of the spatial area included (e.g. member states)

  ## resolution
  else if (data_type == "categ") {
    method_rp <- "near"
    if ( res(x)[1] < res(PU)[1] ){
      f<- (res(PU)[1])/(res(x)[1])
      x<- terra::aggregate(x, fact = f, fun = "modal", na.rm=T)
    }
  }

  ## project to right grid
  x <- terra::project(x, PU, method = method_rp)

  ### mask out unnecessary countries based on shapefile of study area
  x <- mask(x, PU)

  if (!is.factor(x[cells(x)][1, 1])){
    ### fill in the gaps to have the exact same number of grid cells
    x <- sum(x, PU, na.rm = T)
  }
  return (x)
}



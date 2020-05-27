# data mgmt
library(data.table)
library(plyr)
library(lattice)
library(gridExtra)
library(viridis)
library(ggpubr)
library(foreign)
library(tidyr)
library(graphics)
library(tibble)
library(dplyr)
library(stargazer)
library(ggplot2)
library(ggthemes)
library(scales)
library(magrittr)
#devtools::install_github("cardiomoon/ztable")
library(ztable)
library(ggplot2)
library(extrafont)
library(grid)
library(gridExtra)
library(xtable)


# spatial
library(raster)
library(rgdal)
library(velox)
library(rgeos)
library(usdm)

# statistics
library(gstat)
library(car)
library(PerformanceAnalytics)
library(rcompanion)
library(olsrr)
library(fastDummies)
library(lme4)
library(MuMIn)
library(spdep)


#Load raster data for all dates
setwd("H:/My Drive/Research/Canopy_Morphology/Tobacco_Project/2019")

tobacco_area <- readOGR("H:/My Drive/Research/Canopy_Morphology/Tobacco_Project/2019/layers/boundary", "wilson19_boundary", stringsAsFactors = F)
tobacco_area <- spTransform(tobacco_area, CRS("+proj=lcc +lat_1=36.16666666666666 +lat_2=34.33333333333334 +lat_0=33.75 +lon_0=-79 +x_0=609601.22 +y_0=0 +datum=NAD83 +units=m +no_defs +ellps=GRS80 +towgs84=0,0,0"))

csm_619 <- crop(raster('H:/My Drive/Research/Canopy_Morphology/Tobacco_Project/2019/csm/wilson19_619_csm.tif'), tobacco_area)
csm_703 <- crop(raster('H:/My Drive/Research/Canopy_Morphology/Tobacco_Project/2019/csm/wilson19_703_csm.tif'), tobacco_area)
csm_717 <- crop(raster('H:/My Drive/Research/Canopy_Morphology/Tobacco_Project/2019/csm/wilson19_717_csm.tif'), tobacco_area)
dem <- crop(raster('H:/My Drive/Research/Canopy_Morphology/Tobacco_Project/2019/lidar/wilson19_dem.tif'), tobacco_area)

nosoil_619 <- crop(raster('H:/My Drive/Research/Canopy_Morphology/Tobacco_Project/2019/nosoil/csm_nosoil_619.tif'), tobacco_area)
nosoil_703 <- crop(raster('H:/My Drive/Research/Canopy_Morphology/Tobacco_Project/2019/nosoil/csm_nosoil_703.tif'), tobacco_area)
nosoil_717 <- crop(raster('H:/My Drive/Research/Canopy_Morphology/Tobacco_Project/2019/nosoil/csm_nosoil_717.tif'), tobacco_area)

csm_masked_619 <- mask(csm_619, nosoil_619)
csm_masked_703 <- mask(csm_619, nosoil_703)
csm_masked_717 <- mask(csm_619, nosoil_717)

csm_619_velox <- velox(nosoil_619)
csm_703_velox <- velox(nosoil_703)
csm_717_velox <- velox(nosoil_717)

plots <- readOGR("H:/My Drive/Research/Canopy_Morphology/Tobacco_Project/2019/layers/plots", "wilson19_plots", stringsAsFactors = F)
#plots$fertilizer <- as.factor(plots$fertilizer)
#plots$treatment <- as.factor(plots$treatment)
#plots$group <- as.factor(plots$group)
plots <- spTransform(plots, CRS("+proj=lcc +lat_1=36.16666666666666 +lat_2=34.33333333333334 +lat_0=33.75 +lon_0=-79 +x_0=609601.22 +y_0=0 +datum=NAD83 +units=m +no_defs +ellps=GRS80 +towgs84=0,0,0"))
plots@data <-  unite(plots@data, spad_id, block, treatment, sep="_", remove=F)
plots$soil <- as.factor(plots$soil)
plots_df <- dummy_cols(plots@data, select_columns = "soil")
plots@data <- plots_df

#fert_plot_ids <- plots$plotid[grep("-C", c(plots$plotid, "-C"))]
#fert_plots <- plots$Id %in% fert_plot_ids

#plots$nutr_id <- substr(plots$Id, 1, 3)
#plots$nutr_id[fert_plots] <- fert_plot_ids

nutrient <- fread("H:/My Drive/Research/Canopy_Morphology/Tobacco_Project/2019/wilson19_nutrients.csv", stringsAsFactors = F)
#names(nutrient) <- c("Plot","treatment", "Sample #",  "N","P","K","Mg", "Ca","S","Zn", "Mn", "Cu", "Fe", "B", "Al", "Na") 


spad_617 <- fread("H:/My Drive/Research/Canopy_Morphology/Tobacco_Project/2019/wilson19_617_spad.csv", stringsAsFactors = F)
spad_701 <- fread("H:/My Drive/Research/Canopy_Morphology/Tobacco_Project/2019/wilson19_617_spad.csv", stringsAsFactors = F)
spad_718 <- fread("H:/My Drive/Research/Canopy_Morphology/Tobacco_Project/2019/wilson19_617_spad.csv", stringsAsFactors = F)


red_619 <- crop(raster('H:/My Drive/Research/Canopy_Morphology/Tobacco_Project/2019/ortho/wilson/wilson19_619_red_georef.tif'), tobacco_area)
green_619 <- crop(raster('H:/My Drive/Research/Canopy_Morphology/Tobacco_Project/2019/ortho/wilson/wilson19_619_green_georef.tif'), tobacco_area)
blue_619 <- crop(raster('H:/My Drive/Research/Canopy_Morphology/Tobacco_Project/2019/ortho/wilson/wilson19_619_blue_georef.tif'), tobacco_area)
red_703 <- crop(raster('H:/My Drive/Research/Canopy_Morphology/Tobacco_Project/2019/ortho/wilson/wilson19_703_red_georef.tif'), tobacco_area)
green_703 <- crop(raster('H:/My Drive/Research/Canopy_Morphology/Tobacco_Project/2019/ortho/wilson/wilson19_703_green_georef.tif'), tobacco_area)
blue_703 <- crop(raster('H:/My Drive/Research/Canopy_Morphology/Tobacco_Project/2019/ortho/wilson/wilson19_703_blue_georef.tif'), tobacco_area)
red_717 <- crop(raster('H:/My Drive/Research/Canopy_Morphology/Tobacco_Project/2019/ortho/wilson/wilson19_717_red_georef.tif'), tobacco_area)
green_717 <- crop(raster('H:/My Drive/Research/Canopy_Morphology/Tobacco_Project/2019/ortho/wilson/wilson19_717_green_georef.tif'), tobacco_area)
blue_717 <- crop(raster('H:/My Drive/Research/Canopy_Morphology/Tobacco_Project/2019/ortho/wilson/wilson19_717_blue_georef.tif'), tobacco_area)

red_masked_619 <- mask(red_619, nosoil_619)
green_masked_619 <- mask(green_619, nosoil_619)
blue_masked_619 <- mask(blue_619, nosoil_619)

red_masked_703 <- mask(red_703, nosoil_703)
green_masked_703 <- mask(green_703, nosoil_703)
blue_masked_703 <- mask(blue_703, nosoil_703)

red_masked_717 <- mask(red_717, nosoil_717)
green_masked_717 <- mask(green_717, nosoil_717)
blue_masked_717 <- mask(blue_717, nosoil_717)

#----------------------------------------------------Spectral indicies------------------------------------------------------------

# Calculate VARI
varInd_619 <- (green_masked_619-red_masked_619)/(green_masked_619+red_masked_619-blue_masked_619)
varInd_619[varInd_619 < -1] <- -1
varInd_619[varInd_619 > 1] <- 1
varInd_703 <- (green_masked_703-red_masked_703)/(green_masked_703+red_masked_703-blue_masked_703)
varInd_703[varInd_703 < -1] <- -1
varInd_703[varInd_703 > 1] <- 1
varInd_717 <- (green_masked_717-red_masked_703)/(green_masked_717+red_masked_717-blue_masked_717)
varInd_717[varInd_717 < -1] <- -1
varInd_717[varInd_717 > 1] <- 1

# Calculate TGI
tgi_619 <- (green_masked_619-(0.39*red_masked_619)-(0.61*blue_masked_619))/max(c((cellStats(red_masked_619, stat="max")),(cellStats(green_masked_619, stat="max")),(cellStats(blue_masked_619, stat="max"))))
tgi_703 <- (green_masked_703-(0.39*red_masked_703)-(0.61*blue_masked_703))/max(c((cellStats(red_masked_703, stat="max")),(cellStats(green_masked_703, stat="max")),(cellStats(blue_masked_703, stat="max"))))
tgi_717 <- (green_masked_717-(0.39*red_masked_717)-(0.61*blue_masked_717))/max(c((cellStats(red_masked_717, stat="max")),(cellStats(green_masked_717, stat="max")),(cellStats(blue_masked_717, stat="max"))))

varIndRast_list <- list(varInd_619,varInd_703,varInd_717)
tgiRast_list <- list(tgi_619,tgi_703,tgi_717)

# RGB index zonal statistics calculation and exploration

varInd_metrics <- function(x){
  veloxRast <- velox(x)
  extract <- veloxRast$extract(sp=plots)
  names(extract) <- plots$plotid
  varInd_median <- sapply(extract, median, na.rm=T)
  varInd_mean <- sapply(extract, mean, na.rm=T)
  varInd_sd <- sapply(extract, sd, na.rm=T)
  varInd_iqr <- sapply(extract, IQR, na.rm=T)
  varIndDF <- data.frame(plots=plots$plotid, varInd_median= varInd_median, varInd_mean=varInd_mean, varInd_sd = varInd_sd, varInd_iqr = varInd_iqr)
  varIndDF
}

varInd_df_list <- lapply(varIndRast_list, varInd_metrics)

tgi_metrics <- function(x){
  veloxRast <- velox(x)
  extract <- veloxRast$extract(sp=plots)
  names(extract) <- plots$plotid
  tgi_median <- sapply(extract, median, na.rm=T)
  tgi_mean <- sapply(extract, mean, na.rm=T)
  tgi_sd <- sapply(extract, sd, na.rm=T)
  tgi_iqr <- sapply(extract, IQR, na.rm=T)
  tgiDF <- data.frame(plots=plots$plotid, tgi_median= tgi_median, tgi_mean=tgi_mean, tgi_sd = tgi_sd, tgi_iqr = tgi_iqr)
  tgiDF
}

tgi_df_list <- lapply(tgiRast_list, tgi_metrics)

#--------------------------------------------------------- SIMWE water depth model---------------------------------------------------------------
water <- crop(raster("water/wilson19_water.tif"), tobacco_area)

water_velox <- velox(water)
water_extract <- water_velox$extract(sp=plots)
names(water_extract) <- plots$plotid

watermedian <- sapply(water_extract, median, na.rm=TRUE)
watermean <- sapply(water_extract, mean, na.rm=TRUE)
watersd <- sapply(water_extract, sd, na.rm=TRUE)
watersum <- sapply(water_extract, sum, na.rm=TRUE)

waterDF <- data.frame(plots=plots$plotid, watersum = watersum, watersd=watersd, watermean=watermean)
# norm_water <- as.data.frame(scale(waterDF[,2:4]))
# norm_water$plots <- waterDF$plots

#---------------------------------------------------------- Canopy relief ratio----------------------------------------------------------------
crr <- function(x){
  (mean(x, na.rm=T)-min(x, na.rm=T))/(max(x, na.rm=T)-min(x, na.rm=T))
}

crrRast_619 <- focal(csm_619, w=matrix(1/25, nc=5, nr=5), crr)
crrRast_703 <- focal(csm_703, w=matrix(1/25, nc=5, nr=5), crr)
crrRast_717 <- focal(csm_717, w=matrix(1/25, nc=5, nr=5), crr)
crrRast_list <- list(crrRast_619,crrRast_703,crrRast_717)

writeRaster(crrRast_619, "H/My Drive/Research/Canopy_Morphology/Tobacco_Project/2019/results/crr_619.tif", format="GTiff", overwrite = T)
writeRaster(crrRast_703, "H/My Drive/Research/Canopy_Morphology/Tobacco_Project/2019/results/crr_703.tif", format="GTiff", overwrite = T)
writeRaster(crrRast_717, "H/My Drive/Research/Canopy_Morphology/Tobacco_Project/2019/results/crr_717.tif", format="GTiff", overwrite = T)


crr_metrics <- function(x){
  veloxRast <- velox(x)
  extract <- veloxRast$extract(sp=plots)
  names(extract) <- plots$plotid
  crrmean <- sapply(extract, mean, na.rm=T)
  crrmedian <- sapply(extract, median, na.rm=T)
  crrsd <- sapply(extract, sd, na.rm=T)
  crriqr <- sapply(extract, IQR, na.rm=T)
  crrDF <- data.frame(plots=plots$plotid, crrmedian= crrmedian, crrmean=crrmean, crrsd = crrsd, crriqr = crriqr)
  crrDF
}

crr_df_list <- lapply(crrRast_list, crr_metrics)

#---------------------------------------------------------- Crop height ----------------------------------------------------------
plot_extract_619 <- csm_619_velox$extract(sp=plots)
plot_extract_703 <- csm_703_velox$extract(sp=plots)
plot_extract_717 <- csm_717_velox$extract(sp=plots)
names(plot_extract_619) <- plots$plotid
names(plot_extract_703) <- plots$plotid
names(plot_extract_717) <- plots$plotid

extract_list <- list(plot_extract_619, plot_extract_703, plot_extract_717)

ch_metrics <- function(x){
  CHmedian <- sapply(x, median, na.rm=TRUE)
  CHmean <- sapply(x, mean, na.rm=TRUE)
  CHsd <- sapply(x, sd, na.rm=TRUE)
  CHiqr <- sapply(x, IQR, na.rm=TRUE)
  CHvar <- sapply(x, var, na.rm=TRUE)
  CHsum <- sapply(x, sum, na.rm=TRUE)
  CHskew <- sapply(x, skewness, na.rm=TRUE)
  CHkurt <- sapply(x, kurtosis, na.rm=TRUE)
  CHcrr <- sapply(x, crr)
  CHrange <- sapply(x, range, na.rm=TRUE)
  CHrange <- CHrange[2,]-CHrange[1,]
  cropHeightDF <- data.frame(plots=plots$plotid, median= CHmedian, mean=CHmean, sd = CHsd, iqr = CHiqr, var = CHvar, sum = CHsum, skew = CHskew, kurt= CHkurt, PlotCRR = CHcrr, chRange = CHrange)
  cropHeightDF
}

ch_df_list <- lapply(extract_list, ch_metrics)

histogram(plot_extract_703$`503_A`, nint=40) # High kurtosis
histogram(plot_extract_703$`504_G`, nint=40) # Low kurtosis

kurtosis(plot_extract_703$`503_A`, method = "excess")
kurtosis(plot_extract_703$`504_G`, method = "excess")

plotNormalHistogram(plot_extract_703$`503_A`, breaks = 40)
plotNormalHistogram(plot_extract_703$`306_D`, breaks = 40)
#--------------------------------------------------------- rumple index -----------------------------------------------------------------
library(lidR)

rumple <- function(x){
  plotNums <- plots$plotid
  rumple_all <- numeric(length(plotNums))
  names(rumple_all) <- "rumple"
  count <- 1
  for (num in plotNums){
    tmp <- plots@data
    position <- match(num, tmp$plotid)
    extract <- mask(x, plots[position,])
    extract_df <- rasterToPoints(extract)
    rumple <- rumple_index(extract_df[,1], extract_df[,2], extract_df[,3])
    rumple_all[count] <- rumple
    count <- count+1
  }
  rumple_df <- data.frame(plots=plots$plotid, rumple= rumple_all)
  rumple_df
}


csm_list <- list(csm_619, csm_703, csm_717)
rumple_df_list <- lapply(csm_list, rumple)


# ------------------------------------------------ spatial autocorrelation -----------------------------------------------------



autocor_metrics <- function(x, w){
  plotNums <- plots$plotid
  morans_all <- numeric(length(plotNums))
  geary_all <- numeric(length(plotNums))
  count <- 1
  
  for (num in plotNums){
    tmp <- plots@data
    position <- match(num, tmp$plotid)
    extract <- crop(x, plots[position,])
    moran <- Moran(extract, w)
    morans_all[count] <- moran
    geary <- Geary(extract, w)
    geary_all[count] <- geary
    count <- count+1
  }
  
  names(morans_all) <- plots$plotid
  names(geary_all) <- plots$plotid
  autocor_df <- data.frame(plots=plots$plotid, moran= morans_all, geary=geary_all)
  autocor_df
}

#Using approx 65 cm x 65 cm neighborhood matrix (approx in row plant spacing)
autocor_619 <- autocor_metrics(csm_619, matrix(c(rep.int(1,60), 0, rep.int(1,60)), nc=11, nr=11))
autocor_703 <- autocor_metrics(csm_703, matrix(c(rep.int(1,60), 0, rep.int(1,60)), nc=11, nr=11))
autocor_717 <- autocor_metrics(csm_717, matrix(c(rep.int(1,60), 0, rep.int(1,60)), nc=11, nr=11))

local_moran_11x11_619 <- MoranLocal(csm_619, matrix(c(rep.int(1,60), 0, rep.int(1,60)), nc=11, nr=11))
writeRaster(local_moran_11x11_619, "results/local_moran_11x11_619.tif")                  

local_moran_11x11_703 <- MoranLocal(csm_703, matrix(c(rep.int(1,60), 0, rep.int(1,60)), nc=11, nr=11))
writeRaster(local_moran_11x11_703, "results/local_moran_11x11_703.tif")                  

local_moran_11x11_717 <- MoranLocal(csm_717, matrix(c(rep.int(1,60), 0, rep.int(1,60)), nc=11, nr=11))
writeRaster(local_moran_11x11_717, "results/local_moran_11x11_717.tif")                  

# # Create rasters representing local Moran's I for visualization
# autocor_rasters <- function(x, w){
#   plotNums <- plots$plotid
#   moran_local_all <- as.list(x)
#   count <- 1
# 
#   for (num in plotNums){
#     tmp <- plots@data
#     position <- match(num, tmp$plotid)
#     cropped <- crop(x, plots[position,])
#     #masked <- mask(cropped, plots[position,])
#     moran_local <- MoranLocal(cropped, w)
#     #NAvalue(moran_local) <- 0
#     moran_local_all <- append(moran_local_all, moran_local)
#     count <- count+1
#   }
# 
#   moran_local_all <- moran_local_all[-1]
#   names(moran_local_all) <- plots$plotid
#   moran_local_all
# }
# 
# local_moran_plots_619 <- autocor_rasters(csm_619, matrix(c(rep.int(1,60), 0, rep.int(1,60)), nc=11, nr=11))
# names(local_moran_plots_619) <- NULL
# local_moran_plots_619$filename <- 'moranlocal_11x11_619.tif'
# local_moran_plots_619$overwrite <- TRUE
# local_moran_plots_619$fun <- sum
# 
# #merged_local_moran_plots_619 <- do.call(raster::merge, local_moran_plots_619)
# #writeRaster(merged_local_moran_plots_619, "results/local_moran_plots_11x11_619.tif", format="GTiff", overwrite = T)
# 
# local_moran_plots_703 <- autocor_rasters(csm_703, matrix(c(1,1,1,1,1,1,1,1,1,1,1,1,0,1,1,1,1,1,1,1,1,1,1,1,1), nc=5, nr=5))
# names(local_moran_plots_703) <- NULL
# local_moran_plots_703$filename <- 'moran_703.tif'
# local_moran_plots_703$overwrite <- TRUE
# local_moran_plots_703$fun <- sum
# merged_local_moran_plots_703 <- do.call(raster::merge, local_moran_plots_703)
# #writeRaster(merged_local_moran_plots_703, "H/My Drive/Research/Canopy_Morphology/Tobacco_Project/canopy_analysis/layers/moran/local_moran_plots_703.tif", format="GTiff", overwrite = T)
# 
# local_moran_plots_717 <- autocor_rasters(csm_717, matrix(c(1,1,1,1,1,1,1,1,1,1,1,1,0,1,1,1,1,1,1,1,1,1,1,1,1), nc=5, nr=5))
# names(local_moran_plots_717) <- NULL
# local_moran_plots_717$filename <- 'moran_717.tif'
# local_moran_plots_717$overwrite <- TRUE
# local_moran_plots_717$fun <- sum
# merged_local_moran_plots_717 <- do.call(raster::merge, local_moran_plots_717)
# #writeRaster(merged_local_moran_plots_717, "H/My Drive/Research/Canopy_Morphology/Tobacco_Project/canopy_analysis/layers/moran/local_moran_plots_717.tif", format="GTiff", overwrite = T)

# ---------------------------------------- Merge data to plots ----------------------------------------

plots_619 <- plots
plots_703 <- plots
plots_717 <- plots

plots_619 <- merge(plots_619, spad_617, by.x="spad_id", by.y='JOIN')
plots_703 <- merge(plots_703, spad_701, by.x="spad_id", by.y='JOIN')
plots_717 <- merge(plots_717, spad_718, by.x="spad_id", by.y='JOIN')


plots_619 <- merge(plots_619, varInd_df_list[[1]], by.x="plotid", by.y="plots")
plots_703 <- merge(plots_703, varInd_df_list[[2]], by.x="plotid", by.y="plots")
plots_717 <- merge(plots_717, varInd_df_list[[3]], by.x="plotid", by.y="plots")

plots_619 <- merge(plots_619, tgi_df_list[[1]], by.x="plotid", by.y="plots")
plots_703 <- merge(plots_703, tgi_df_list[[2]], by.x="plotid", by.y="plots")
plots_717 <- merge(plots_717, tgi_df_list[[3]], by.x="plotid", by.y="plots")


plots_619 <- merge(plots_619, waterDF, by.x="plotid", by.y="plots")
plots_703 <- merge(plots_703, waterDF, by.x="plotid", by.y="plots")
plots_717 <- merge(plots_717, waterDF, by.x="plotid", by.y="plots")


plots_619 <- merge(plots_619, crr_df_list[[1]], by.x="plotid", by.y="plots")
plots_703 <- merge(plots_703, crr_df_list[[2]], by.x="plotid", by.y="plots")
plots_717 <- merge(plots_717, crr_df_list[[3]], by.x="plotid", by.y="plots")


#merge all metrics into spatial polygons
plots_619 <- merge(plots_619, ch_df_list[[1]], by.x="plotid", by.y="plots")
plots_703 <- merge(plots_703, ch_df_list[[2]], by.x="plotid", by.y="plots")
plots_717 <- merge(plots_717, ch_df_list[[3]], by.x="plotid", by.y="plots")


#merge all metrics into spatial polygons
plots_619 <- merge(plots_619, rumple_df_list[[1]], by.x="plotid", by.y="plots")
plots_703 <- merge(plots_703, rumple_df_list[[2]], by.x="plotid", by.y="plots")
plots_717 <- merge(plots_717, rumple_df_list[[3]], by.x="plotid", by.y="plots")

plots_619 <- merge(plots_619, autocor_619, by.x="plotid", by.y="plots")
plots_703 <- merge(plots_703, autocor_703, by.x="plotid", by.y="plots")
plots_717 <- merge(plots_717, autocor_717, by.x="plotid", by.y="plots")

spplot(plots_619, zcol = 'kurt')
spplot(plots_619, zcol = 'mean')
spplot(plots_619, zcol = 'crrmedian')
spplot(plots_619, zcol = 'crriqr')
spplot(plots_619, zcol = 'rumple')
spplot(plots_619, zcol = 'moran')

spplot(plots_703, zcol = 'kurt')
spplot(plots_703, zcol = 'mean')
spplot(plots_703, zcol = 'crrmedian')
spplot(plots_703, zcol = 'crriqr')
spplot(plots_703, zcol = 'rumple')
spplot(plots_703, zcol = 'moran')

spplot(plots_717, zcol = 'kurt')
spplot(plots_717, zcol = 'mean')
spplot(plots_717, zcol = 'crrmedian')
spplot(plots_717, zcol = 'crriqr')
spplot(plots_717, zcol = 'rumple')
spplot(plots_717, zcol = 'moran')

write.csv(plots_619@data, "results/data_619.csv")
write.csv(plots_703@data, "results/data_703.csv")
write.csv(plots_717@data, "results/data_717.csv")

# ------------------------------------- Aggregate plots to match nutrient observations --------------------------------------
agg_619 <- aggregate(plots_619@data[,18:46], by = list(plots_619$join_id), FUN = mean)
agg_703 <- aggregate(plots_703@data[,18:46], by = list(plots_703$join_id), FUN = mean)
agg_717 <- aggregate(plots_717@data[,18:46], by = list(plots_717$join_id), FUN = mean)

agg_619_nut1 <- merge(agg_619, nutrient[Date==617], by.x="Group.1", by.y='Join_ID')
names(agg_619_nut1)[2] <- "spad"

agg_703_nut1 <- merge(agg_703, nutrient[Date==617], by.x="Group.1", by.y='Join_ID')
names(agg_703_nut1)[2] <- "spad"

agg_703_nut2 <- merge(agg_703, nutrient[Date==701], by.x="Group.1", by.y='Join_ID')
names(agg_703_nut2)[2] <- "spad"

agg_717_nut1 <- merge(agg_717, nutrient[Date==617], by.x="Group.1", by.y='Join_ID')
names(agg_717_nut1)[2] <- "spad"

agg_717_nut2 <- merge(agg_717, nutrient[Date==701], by.x="Group.1", by.y='Join_ID')
names(agg_717_nut2)[2] <- "spad"

agg_717_nut3 <- merge(agg_717, nutrient[Date==718], by.x="Group.1", by.y='Join_ID')
names(agg_717_nut3)[2] <- "spad"

# # ------------------------------------ Correlation tests ------------------------------------------
library(ggpubr)
# 
pairs(z_619[,c(20, 21, 26, 29, 31, 32, 34, 35,37, 38)], main="Simple Scatterplot Matrix")

cor(z_619[,c(19:35,37:39)])

pairs(agg_619_nut1[,c("N","P","K","B")])
pairs(agg_703_nut2[,c("N","P","K","B")])
pairs(agg_717_nut3[,c("N","P","K","B")])

cor.test(agg_619_nut1[,c("N")], agg_619_nut1[,c("P")], method = "pearson")
cor.test(agg_619_nut1[,c("N")], agg_619_nut1[,c("K")], method = "pearson")
cor.test(agg_619_nut1[,c("N")], agg_619_nut1[,c("B")], method = "pearson")
cor.test(agg_619_nut1[,c("P")], agg_619_nut1[,c("K")], method = "pearson")
cor.test(agg_619_nut1[,c("P")], agg_619_nut1[,c("B")], method = "pearson")
cor.test(agg_619_nut1[,c("K")], agg_619_nut1[,c("B")], method = "pearson")

cor.test(agg_703_nut2[,c("N")], agg_703_nut2[,c("P")], method = "pearson")
cor.test(agg_703_nut2[,c("N")], agg_703_nut2[,c("K")], method = "pearson")
cor.test(agg_703_nut2[,c("N")], agg_703_nut2[,c("B")], method = "pearson")
cor.test(agg_703_nut2[,c("P")], agg_703_nut2[,c("K")], method = "pearson")
cor.test(agg_703_nut2[,c("P")], agg_703_nut2[,c("B")], method = "pearson")
cor.test(agg_703_nut2[,c("K")], agg_703_nut2[,c("B")], method = "pearson")

cor.test(agg_717_nut3[,c("N")], agg_717_nut3[,c("P")], method = "pearson")
cor.test(agg_717_nut3[,c("N")], agg_717_nut3[,c("K")], method = "pearson")
cor.test(agg_717_nut3[,c("N")], agg_717_nut3[,c("B")], method = "pearson")
cor.test(agg_717_nut3[,c("P")], agg_717_nut3[,c("K")], method = "pearson")
cor.test(agg_717_nut3[,c("P")], agg_717_nut3[,c("B")], method = "pearson")
cor.test(agg_717_nut3[,c("K")], agg_717_nut3[,c("B")], method = "pearson")

correlation.matrix <- cor(agg_619_nut1[,c("N","P","K","B")])
stargazer(correlation.matrix, title="Correlation Matrix")

correlation.matrix <- cor(agg_703_nut2[,c("N","P","K","B")])
stargazer(correlation.matrix, title="Correlation Matrix")

correlation.matrix <- cor(agg_717_nut3[,c("N","P","K","B")])
stargazer(correlation.matrix, title="Correlation Matrix")


# #plot(plots_619$mean, plots_619$BuAc)
# plot(plots_619$N, plots_619$K, xlab = "Median Height (m)", ylab = "Yield (bu/ac)", cex = 1,cex.lab = 2, cex.axis = 2, pch=16)
# lm_median <- lm(BuAc ~ median, data = plots_1cm)
# text(par()$usr[1], par()$usr[4], bquote(R^2 ~ "=" ~ .(round(summary(lm_median)$r.squared, 2))), 
#      cex = 2, col = 2, adj = c(-2.1, 1.6))
# shapiro.test(plots_1cm$BuAc)
# 
# cor.test(plots_619$N, plots_619$K)
# ggscatter(plots_1cm@data, x = "median", y = "BuAc", 
#           add = "reg.line", conf.int = TRUE, 
#           cor.coef = TRUE, cor.method = "pearson",
#           xlab = "Median Height (m)", ylab = "Yield (bu/ac)")
# 
# #plot(plots_619$sum, plots_619$BuAc)
# plot(plots_619$iqr, plots_619$BuAc, xlab = "Interquartile Range (m)", ylab = "Yield (bu/ac)", cex = 1,cex.lab = 2, cex.axis = 2, pch=16)
# lm_iqr <- lm(BuAc ~ iqr, data = plots_1cm)
# text(par()$usr[1], par()$usr[4], bquote(R^2 ~ "=" ~ .(round(summary(lm_iqr)$r.squared, 2))), 
#      cex = 2, col = 2, adj = c(-2.1, 1.6))
# 
# plot(plots_619$rumple, plots_619$BuAc, xlab = "Rumple Index", ylab = "Yield (bu/ac)", cex = 1,cex.lab = 2, cex.axis = 2, pch=16)
# lm_rumple <- lm(BuAc ~ rumple, data = plots_1cm)
# text(par()$usr[1], par()$usr[4], bquote(R^2 ~ "=" ~ .(round(summary(lm_rumple)$r.squared, 2))), 
#      cex = 2, col = 2, adj = c(-2.1, 1.6))
# 
# plot(plots_619$skew, plots_619$BuAc, xlab = "Skewness", ylab = "Yield (bu/ac)", cex = 1,cex.lab = 2, cex.axis = 2, pch=16)
# lm_skew <- lm(BuAc ~ skew, data = plots_1cm)
# text(par()$usr[1], par()$usr[4], bquote(R^2 ~ "=" ~ .(round(summary(lm_skew)$r.squared, 2))), 
#      cex = 2, col = 2, adj = c(-2.1, 1.6))
# 
# plot(plots_619$kurt, plots_619$BuAc, xlab = "Kurtosis", ylab = "Yield (bu/ac)", cex = 1,cex.lab = 2, cex.axis = 2, pch=16)
# lm_kurt <- lm(BuAc ~ kurt, data = plots_1cm)
# text(par()$usr[1], par()$usr[4], bquote(R^2 ~ "=" ~ .(round(summary(lm_kurt)$r.squared, 2))), 
#      cex = 2, col = 2, adj = c(-2.1, 1.6))
# 
# #plot(plots_619$sd, plots_619$BuAc)
# plot(plots_619$crrmean, plots_619$BuAc, xlab = "CRR Mean", ylab = "Yield (bu/ac)", cex = 1,cex.lab = 2, cex.axis = 2, pch=16)
# lm_crrmean <- lm(BuAc ~ crrmean, data = plots_1cm)
# text(par()$usr[1], par()$usr[4], bquote(R^2 ~ "=" ~ .(round(summary(lm_crrmean)$r.squared, 2))), 
#      cex = 2, col = 2, adj = c(-2.1, 1.6))
# 
# plot(plots_619$crrsd, plots_619$BuAc, xlab = "CRR Standard Deviation", ylab = "Yield (bu/ac)", cex = 1,cex.lab = 2, cex.axis = 2, pch=16)
# lm_crrsd <- lm(BuAc ~ crrsd, data = plots_1cm)
# text(par()$usr[1], par()$usr[4], bquote(R^2 ~ "=" ~ .(round(summary(lm_crrsd)$r.squared, 2))), 
#      cex = 2, col = 2, adj = c(-3.4, 1.6))
# 
# #plot(plots_619$PlotCRR, plots_619$BuAc)
# # plot(plots_619$moran, plots_619$BuAc)
# # lm_moran <- lm(BuAc ~ moran, data = plots_1cm)
# # text(par()$usr[1], par()$usr[4], bquote(R^2 ~ "=" ~ .(round(summary(lm_moran)$r.squared, 2))), 
# #      cex = 2, col = 2, adj = c(-2.1, 1.6))
# 
# plot(plots_1cm$geary, plots_619$BuAc, xlab = "Geary's C", ylab = "Yield (bu/ac)", cex = 1,cex.lab = 2, cex.axis = 2, pch=16)
# lm_geary <- lm(BuAc ~ geary, data = plots_1cm)
# text(par()$usr[1], par()$usr[4], bquote(R^2 ~ "=" ~ .(round(summary(lm_geary)$r.squared, 2))), 
#      cex = 2, col = 2, adj = c(-3.4, 1.6))
# 
# #tmp <- lm(BuAc~median, data=plots_619)
# #summary(tmp)
# #plot(plots_619$sill, plots_619$BuAc)
# #plot(plots_619$range, plots_619$BuAc)
# #plot(plots_619$Trt, plots_619$BuAc)
# plot(plots_619$AUDPC, plots_619$BuAc, xlab = "AUDPC", ylab = "Yield (bu/ac)")
# lm_audpc <- lm(BuAc ~ AUDPC, data = plots_1cm)
# text(par()$usr[1], par()$usr[4], bquote(R^2 ~ "=" ~ .(round(summary(lm_audpc)$r.squared, 2))), 
#      cex = 1.5, col = 2, adj = c(-6.7, 1.4))
# #plot(plots_619$watersum, plots_619$BuAc)
# #plot(plots_619$watersd, plots_619$BuAc)
# 
# par(mfrow=c(1,1))
# 
# # Plot histogram distribution of each variable
# histogram(agg_619$N, nint=40)
# histogram(agg_619$P, nint=40)
# histogram(agg_619$K, nint=40)
# histogram(agg_619$Mg, nint=40)
# histogram(agg_619$Ca, nint=40)
# histogram(agg_619$S, nint=40)
# histogram(agg_619$Zn, nint=40)
# histogram(agg_619$Mn, nint=40)
# histogram(agg_619$Cu, nint=40)
# histogram(agg_619$Fe, nint=40)
# histogram(agg_619$B, nint=40)
# histogram(agg_619$Al, nint=40)
# 

# --------------------------------------- Stats data prep --------------------------------------------------------------------------------
# grouped_soils <- function(x)
# {
#   paste(as.vector(sort(plots_619$soil[plots_619$join_id == x]))[1],as.vector(sort(plots_619$soil[plots_619$join_id == x]))[2],as.vector(sort(plots_619$soil[plots_619$join_id == x]))[3], sep = "_")
#   
# }
# 
# groups <- unique(plots_619$join_id)
# soils <- sapply(groups, grouped_soils)
# groups <- as.data.frame(groups)
# groups$soils <- as.factor(soils)
# 
# agg_619 <- merge(agg_619, groups, by.x="Group.1", by.y="groups")
# agg_619$treatment <- tstrsplit(agg_619$Group.1, "_")[[2]]
# agg_619$fert <- tstrsplit(agg_619$Group.1, "_")[[1]]
# agg_703 <- merge(agg_703, groups, by.x="Group.1", by.y="groups")
# agg_703$treatment <- tstrsplit(agg_703$Group.1, "_")[[2]]
# agg_703$fert <- tstrsplit(agg_703$Group.1, "_")[[1]]
# agg_717 <- merge(agg_717, groups, by.x="Group.1", by.y="groups")
# agg_717$treatment <- tstrsplit(agg_717$Group.1, "_")[[2]]
# agg_717$fert <- tstrsplit(agg_717$Group.1, "_")[[1]]
# 
# write.csv(agg_619, "results/data_619.csv")
# write.csv(agg_703, "results/data_703.csv")
# write.csv(agg_717, "results/data_717.csv")

normalize <- function(x)
{
  return((x- min(x)) /(max(x)-min(x)))
}


#minmax normalization of data
norm_619_nut1 <- as.data.frame(normalize(agg_619_nut1[,c(2:30,39:51)]))
norm_619_nut1$Id <- agg_619_nut1$Group.1
norm_619_nut1$treatment <- agg_619_nut1$treatment
norm_619_nut1$fert <- agg_619_nut1$fert
#norm_619$soils <- agg_619$soils
#norm_619 <- dummy_cols(norm_619, select_columns = "soils")

norm_703_nut1 <- as.data.frame(normalize(agg_703_nut1[,c(2:30,39:51)]))
norm_703_nut1$Id <- agg_703_nut1$Group.1
norm_703_nut1$treatment <- agg_703_nut1$treatment
norm_703_nut1$fert <- agg_703_nut1$fert

norm_703_nut2 <- as.data.frame(normalize(agg_703_nut2[,c(2:30,39:51)]))
norm_703_nut2$Id <- agg_703_nut2$Group.1
norm_703_nut2$treatment <- agg_703_nut2$treatment
norm_703_nut2$fert <- agg_703_nut2$fert

norm_717_nut1 <- as.data.frame(normalize(agg_717_nut1[,c(2:30,39:51)]))
norm_717_nut1$Id <- agg_717_nut1$Group.1
norm_717_nut1$treatment <- agg_717_nut1$treatment
norm_717_nut1$fert <- agg_717_nut1$fert

norm_717_nut2 <- as.data.frame(normalize(agg_717_nut2[,c(2:30,39:51)]))
norm_717_nut2$Id <- agg_717_nut2$Group.1
norm_717_nut2$treatment <- agg_717_nut2$treatment
norm_717_nut2$fert <- agg_717_nut2$fert

norm_717_nut3 <- as.data.frame(normalize(agg_717_nut3[,c(2:30,39:51)]))
norm_717_nut3$Id <- agg_717_nut3$Group.1
norm_717_nut3$treatment <- agg_717_nut3$treatment
norm_717_nut3$fert <- agg_717_nut3$fert

#z score normalization of data
z_619_nut1 <- as.data.frame(scale(agg_619_nut1[,c(2:30,39:51)]))
z_619_nut1$Id <- agg_619_nut1$Group.1
z_619_nut1$treatment <- agg_619_nut1$treatment
z_619_nut1$fert <- agg_619_nut1$fert
#z_619$soils <- agg_619$soils
#z_619 <- dummy_cols(z_619, select_columns = "soils")

z_703_nut1 <- as.data.frame(scale(agg_703_nut1[,c(2:30,39:51)]))
z_703_nut1$Id <- agg_703_nut1$Group.1
z_703_nut1$treatment <- agg_703_nut1$treatment
z_703_nut1$fert <- agg_703_nut1$fert

z_703_nut2 <- as.data.frame(scale(agg_703_nut2[,c(2:30,39:51)]))
z_703_nut2$Id <- agg_703_nut2$Group.1
z_703_nut2$treatment <- agg_703_nut2$treatment
z_703_nut2$fert <- agg_703_nut2$fert

z_717_nut1 <- as.data.frame(scale(agg_717_nut1[,c(2:30,39:51)]))
z_717_nut1$Id <- agg_717_nut1$Group.1
z_717_nut1$treatment <- agg_717_nut1$treatment
z_717_nut1$fert <- agg_717_nut1$fert

z_717_nut2 <- as.data.frame(scale(agg_717_nut2[,c(2:30,39:51)]))
z_717_nut2$Id <- agg_717_nut2$Group.1
z_717_nut2$treatment <- agg_717_nut2$treatment
z_717_nut2$fert <- agg_717_nut2$fert

z_717_nut3 <- as.data.frame(scale(agg_717_nut3[,c(2:30,39:51)]))
z_717_nut3$Id <- agg_717_nut3$Group.1
z_717_nut3$treatment <- agg_717_nut3$treatment
z_717_nut3$fert <- agg_717_nut3$fert


#---------------------------------------- LM structural only (partial mask) ------------------------------------------
lmp <- function (modelobject) {
  if (class(modelobject) != "lm") stop("Not an object of class 'lm' ")
  f <- summary(modelobject)$fstatistic
  p <- pf(f[1],f[2],f[3],lower.tail=F)
  attributes(p) <- NULL
  return(p)
}


vif_func<-function(in_frame,thresh=10,trace=T,...){
  
  library(fmsb)
  
  if(any(!'data.frame' %in% class(in_frame))) in_frame<-data.frame(in_frame)
  
  #get initial vif value for all comparisons of variables
  vif_init<-NULL
  var_names <- names(in_frame)
  for(val in var_names){
    regressors <- var_names[-which(var_names == val)]
    form <- paste(regressors, collapse = '+')
    form_in <- formula(paste(val, '~', form))
    vif_init<-rbind(vif_init, c(val, VIF(lm(form_in, data = in_frame, ...))))
  }
  vif_max<-max(as.numeric(vif_init[,2]), na.rm = TRUE)
  
  if(vif_max < thresh){
    if(trace==T){ #print output of each iteration
      prmatrix(vif_init,collab=c('var','vif'),rowlab=rep('',nrow(vif_init)),quote=F)
      cat('\n')
      cat(paste('All variables have VIF < ', thresh,', max VIF ',round(vif_max,2), sep=''),'\n\n')
    }
    return(var_names)
  }
  else{
    
    in_dat<-in_frame
    
    #backwards selection of explanatory variables, stops when all VIF values are below 'thresh'
    while(vif_max >= thresh){
      
      vif_vals<-NULL
      var_names <- names(in_dat)
      
      for(val in var_names){
        regressors <- var_names[-which(var_names == val)]
        form <- paste(regressors, collapse = '+')
        form_in <- formula(paste(val, '~', form))
        vif_add<-VIF(lm(form_in, data = in_dat, ...))
        vif_vals<-rbind(vif_vals,c(val,vif_add))
      }
      max_row<-which(vif_vals[,2] == max(as.numeric(vif_vals[,2]), na.rm = TRUE))[1]
      
      vif_max<-as.numeric(vif_vals[max_row,2])
      
      if(vif_max<thresh) break
      
      if(trace==T){ #print output of each iteration
        prmatrix(vif_vals,collab=c('var','vif'),rowlab=rep('',nrow(vif_vals)),quote=F)
        cat('\n')
        cat('removed: ',vif_vals[max_row,1],vif_max,'\n\n')
        flush.console()
      }
      
      in_dat<-in_dat[,!names(in_dat) %in% vif_vals[max_row,1]]
      
    }
    
    return(names(in_dat))
    
  }
  
}

xvar_619 <- z_619[,c(19:35,37:39)] # omit PlotCRR and VARI
xvar_703 <- z_703[,c(19:35,37:39)]
xvar_717 <- z_717[,c(19:35,37:39)]

xkeep_all_619 <- vif_func(xvar_619, 4, F)
xkeep_all_703 <- vif_func(xvar_703, 4, F)
xkeep_all_717 <- vif_func(xvar_717, 4, F)

xkeep_struct_619 <- xkeep_all_619[3:9]
xkeep_spec_619 <- xkeep_all_619[1:2]

xkeep_struct_703 <- xkeep_all_703[3:9]
xkeep_spec_703 <- xkeep_all_703[1:2]

xkeep_struct_717 <- xkeep_all_717[3:9]
xkeep_spec_717 <- xkeep_all_717[1:2]

yvar <- c("N", "P","K", "B")


lm_compare <- function(y, xkeep, data){
  #control <- yvar[!(yvar %in% y)]
  xvars <- paste(c("0", xkeep),collapse='+')
  eq<-as.formula(paste(y, xvars, sep="~"))
  mod1<-lm(eq,data=data, na.action=na.fail)
  modsel <- dredge(mod1)
  topmod <- get.models(modsel, subset = 1)[[1]]
  res <- as.data.frame(modsel[1])
  res$adj_r2 <- summary(topmod)$adj.r.squared
  res$yvar <- y
  res
}

topmod_res <- function(y, xkeep, data){
  #control <- yvar[!(yvar %in% y)]
  xvars <- paste(c("0", xkeep),collapse='+')
  eq<-as.formula(paste(y, xvars, sep="~"))
  mod1<-lm(eq,data=data, na.action=na.fail)
  modsel <- dredge(mod1)
  topmod <- get.models(modsel, subset = 1)[[1]]
  topmod
}

xvars_all_fixed <- c("mean","kurt","rumple","crrmedian", "crriqr","moran","watersum","tgi_mean","tgi_sd")
xvars_spec_fixed <- xvars_all_fixed[c(8:9)]
xvars_struc_fixed <- xvars_all_fixed[1:7]

# ---------------------------------- Structural LMs ----------------------------------------------
struc_619_nut1 <- lapply(X=yvar, FUN=lm_compare, xkeep=xvars_struc_fixed, data=z_619_nut1)
struc_619_nut1 <- do.call(rbind, lapply(struc_619_nut1, data.frame))
struc_mods_619_nut1 <- lapply(yvar, topmod_res, xvars_struc_fixed, z_619_nut1)

struc_703_nut1 <- lapply(X=yvar, FUN=lm_compare, xkeep=xvars_struc_fixed, data=z_703_nut1)
struc_703_nut1 <- do.call(rbind, lapply(struc_703_nut1, data.frame))
struc_mods_703_nut1 <- lapply(yvar, topmod_res, xvars_struc_fixed, z_703_nut1)

struc_703_nut2 <- lapply(X=yvar, FUN=lm_compare, xkeep=xvars_struc_fixed, data=z_703_nut2)
struc_703_nut2 <- do.call(rbind, lapply(struc_703_nut2, data.frame))
struc_mods_703_nut2 <- lapply(yvar, topmod_res, xvars_struc_fixed, z_703_nut2)

struc_717_nut1 <- lapply(X=yvar, FUN=lm_compare, xkeep=xvars_struc_fixed, data=z_717_nut1)
struc_717_nut1 <- do.call(rbind, lapply(struc_717_nut1, data.frame))
struc_mods_717_nut1 <- lapply(yvar, topmod_res, xvars_struc_fixed, z_717_nut1)

struc_717_nut2 <- lapply(X=yvar, FUN=lm_compare, xkeep=xvars_struc_fixed, data=z_717_nut2)
struc_717_nut2 <- do.call(rbind, lapply(struc_717_nut2, data.frame))
struc_mods_717_nut2 <- lapply(yvar, topmod_res, xvars_struc_fixed, z_717_nut2)

struc_717_nut3 <- lapply(X=yvar, FUN=lm_compare, xkeep=xvars_struc_fixed, data=z_717_nut3)
struc_717_nut3 <- do.call(rbind, lapply(struc_717_nut3, data.frame))
struc_mods_717_nut3 <- lapply(yvar, topmod_res, xvars_struc_fixed, z_717_nut3)

# ---------------------------------- Spectral lm_compares ----------------------------------------------

spec_619_nut1 <- lapply(X=yvar, FUN=lm_compare, xkeep=xvars_spec_fixed, data=z_619_nut1)
spec_619_nut1 <- do.call(rbind, lapply(spec_619_nut1, data.frame))
spec_mods_619_nut1 <- lapply(yvar, topmod_res, xvars_spec_fixed, z_619_nut1)

spec_703_nut1 <- lapply(X=yvar, FUN=lm_compare, xkeep=xvars_spec_fixed, data=z_703_nut1)
spec_703_nut1 <- do.call(rbind, lapply(spec_703_nut1, data.frame))
spec_mods_703_nut1 <- lapply(yvar, topmod_res, xvars_spec_fixed, z_703_nut1)

spec_703_nut2 <- lapply(X=yvar, FUN=lm_compare, xkeep=xvars_spec_fixed, data=z_703_nut2)
spec_703_nut2 <- do.call(rbind, lapply(spec_703_nut2, data.frame))
spec_mods_703_nut2 <- lapply(yvar, topmod_res, xvars_spec_fixed, z_703_nut2)

spec_717_nut1 <- lapply(X=yvar, FUN=lm_compare, xkeep=xvars_spec_fixed, data=z_717_nut1)
spec_717_nut1 <- do.call(rbind, lapply(spec_717_nut1, data.frame))
spec_mods_717_nut1 <- lapply(yvar, topmod_res, xvars_spec_fixed, z_717_nut1)

spec_717_nut2 <- lapply(X=yvar, FUN=lm_compare, xkeep=xvars_spec_fixed, data=z_717_nut2)
spec_717_nut2 <- do.call(rbind, lapply(spec_717_nut2, data.frame))
spec_mods_717_nut2 <- lapply(yvar, topmod_res, xvars_spec_fixed, z_717_nut2)

spec_717_nut3 <- lapply(X=yvar, FUN=lm_compare, xkeep=xvars_spec_fixed, data=z_717_nut3)
spec_717_nut3 <- do.call(rbind, lapply(spec_717_nut3, data.frame))
spec_mods_717_nut3 <- lapply(yvar, topmod_res, xvars_spec_fixed, z_717_nut3)

# ---------------------------------- Structural and Spectral lm_compares ----------------------------------------------

spec_struc_619_nut1 <- lapply(X=yvar, FUN=lm_compare, xkeep=xvars_all_fixed, data=z_619_nut1)
spec_struc_619_nut1 <- do.call(rbind, lapply(spec_struc_619_nut1, data.frame))
spec_struc_mods_619_nut1 <- lapply(yvar, topmod_res, xvars_all_fixed, z_619_nut1)

spec_struc_703_nut1 <- lapply(X=yvar, FUN=lm_compare, xkeep=xvars_all_fixed, data=z_703_nut1)
spec_struc_703_nut1 <- do.call(rbind, lapply(spec_struc_703_nut1, data.frame))
spec_struc_mods_703_nut1 <- lapply(yvar, topmod_res, xvars_all_fixed, z_703_nut1)

spec_struc_703_nut2 <- lapply(X=yvar, FUN=lm_compare, xkeep=xvars_all_fixed, data=z_703_nut2)
spec_struc_703_nut2 <- do.call(rbind, lapply(spec_struc_703_nut2, data.frame))
spec_struc_mods_703_nut2 <- lapply(yvar, topmod_res, xvars_all_fixed, z_703_nut2)

spec_struc_717_nut1 <- lapply(X=yvar, FUN=lm_compare, xkeep=xvars_all_fixed, data=z_717_nut1)
spec_struc_717_nut1 <- do.call(rbind, lapply(spec_struc_717_nut1, data.frame))
spec_struc_mods_717_nut1 <- lapply(yvar, topmod_res, xvars_all_fixed, z_717_nut1)

spec_struc_717_nut2 <- lapply(X=yvar, FUN=lm_compare, xkeep=xvars_all_fixed, data=z_717_nut2)
spec_struc_717_nut2 <- do.call(rbind, lapply(spec_struc_717_nut2, data.frame))
spec_struc_mods_717_nut2 <- lapply(yvar, topmod_res, xvars_all_fixed, z_717_nut2)

spec_struc_717_nut3 <- lapply(X=yvar, FUN=lm_compare, xkeep=xvars_all_fixed, data=z_717_nut3)
spec_struc_717_nut3 <- do.call(rbind, lapply(spec_struc_717_nut3, data.frame))
spec_struc_mods_717_nut3 <- lapply(yvar, topmod_res, xvars_all_fixed, z_717_nut3)


# ------------------------------------------- Write results -----------------------------------------------
col.order.struc<- c("mean", "kurt", "rumple", "crrmedian", "crriqr", "moran", "watersum", "df", "logLik", "AICc", "delta", "weight", "adj_r2", "yvar")
col.order.spec <- c("tgi_mean", "tgi_sd", "df", "logLik", "AICc", "delta", "weight", "adj_r2", "yvar")
col.order.spec_struc <- c("mean", "kurt", "rumple", "crrmedian", "crriqr", "moran", "watersum", "tgi_mean", "tgi_sd","df", "logLik", "AICc", "delta", "weight", "adj_r2", "yvar")

struc_619_nut1 <- struc_619_nut1[, col.order.struc]
struc_703_nut1 <- struc_703_nut1[, col.order.struc]
struc_703_nut2 <- struc_703_nut2[, col.order.struc]
struc_717_nut1 <- struc_717_nut1[, col.order.struc]
struc_717_nut2 <- struc_717_nut2[, col.order.struc]
struc_717_nut3 <- struc_717_nut3[, col.order.struc]

spec_619_nut1 <- spec_619_nut1[, col.order.spec]
spec_703_nut1 <- spec_703_nut1[, col.order.spec]
spec_703_nut2 <- spec_703_nut2[, col.order.spec]
spec_717_nut1 <- spec_717_nut1[, col.order.spec]
spec_717_nut2 <- spec_717_nut2[, col.order.spec]
spec_717_nut3 <- spec_717_nut3[, col.order.spec]

spec_struc_619_nut1 <- spec_struc_619_nut1[, col.order.spec_struc]
spec_struc_703_nut1 <- spec_struc_703_nut1[, col.order.spec_struc]
spec_struc_703_nut2 <- spec_struc_703_nut2[, col.order.spec_struc]
spec_struc_717_nut1 <- spec_struc_717_nut1[, col.order.spec_struc]
spec_struc_717_nut2 <- spec_struc_717_nut2[, col.order.spec_struc]
spec_struc_717_nut3 <- spec_struc_717_nut3[, col.order.spec_struc]

write.csv(struc_619_nut1, "results/struc_619_nut1.csv")
write.csv(struc_703_nut1, "results/struc_703_nut1.csv")
write.csv(struc_703_nut2, "results/struc_703_nut2.csv")
write.csv(struc_717_nut1, "results/struc_717_nut1.csv")
write.csv(struc_717_nut2, "results/struc_717_nut2.csv")
write.csv(struc_717_nut3, "results/struc_717_nut3.csv")

write.csv(spec_619_nut1, "results/spec_619_nut1.csv")
write.csv(spec_703_nut1, "results/spec_703_nut1.csv")
write.csv(spec_703_nut2, "results/spec_703_nut2.csv")
write.csv(spec_717_nut1, "results/spec_717_nut1.csv")
write.csv(spec_717_nut2, "results/spec_717_nut2.csv")
write.csv(spec_717_nut3, "results/spec_717_nut3.csv")

write.csv(spec_struc_619_nut1, "results/spec_struc_619_nut1.csv")
write.csv(spec_struc_703_nut1, "results/spec_struc_703_nut1.csv")
write.csv(spec_struc_703_nut2, "results/spec_struc_703_nut2.csv")
write.csv(spec_struc_717_nut1, "results/spec_struc_717_nut1.csv")
write.csv(spec_struc_717_nut2, "results/spec_struc_717_nut2.csv")
write.csv(spec_struc_717_nut3, "results/spec_struc_717_nut3.csv")

get_coef <- function(mod){
  df <- as.data.frame(coef(summary(mod)))
  df$xvar <- names(mod$coefficients)
  df
}



struc_619_all_nut1 <- lapply(struc_mods_619_nut1, get_coef)
struc_619_all_nut1 = bind_rows(struc_619_all_nut1, .id = "column_label")
write.csv(struc_619_all_nut1, "results/struc_619_coefs_nut1.csv")

spec_619_all_nut1 <- lapply(spec_mods_619_nut1, get_coef)
spec_619_all_nut1 = bind_rows(spec_619_all_nut1, .id = "column_label")
write.csv(spec_619_all_nut1, "results/spec_619_coefs_nut1.csv")

spec_struc_619_all_nut1 <- lapply(spec_struc_mods_619_nut1, get_coef)
spec_struc_619_all_nut1 = bind_rows(spec_struc_619_all_nut1, .id = "column_label")
write.csv(spec_struc_619_all_nut1, "results/spec_struc_619_coefs_nut1.csv")



struc_703_all_nut1 <- lapply(struc_mods_703_nut1, get_coef)
struc_703_all_nut1 = bind_rows(struc_703_all_nut1, .id = "column_label")
write.csv(struc_703_all_nut1, "results/struc_703_coefs_nut1.csv")

spec_703_all_nut1 <- lapply(spec_mods_703_nut1, get_coef)
spec_703_all_nut1 = bind_rows(spec_703_all_nut1, .id = "column_label")
write.csv(spec_703_all_nut1, "results/spec_703_coefs_nut1.csv")

spec_struc_703_all_nut1 <- lapply(spec_struc_mods_703_nut1, get_coef)
spec_struc_703_all_nut1 = bind_rows(spec_struc_703_all_nut1, .id = "column_label")
write.csv(spec_struc_703_all_nut1, "results/spec_struc_703_coefs_nut1.csv")



struc_703_all_nut2 <- lapply(struc_mods_703_nut2, get_coef)
struc_703_all_nut2 = bind_rows(struc_703_all_nut2, .id = "column_label")
write.csv(struc_703_all_nut2, "results/struc_703_coefs_nut2.csv")

spec_703_all_nut2 <- lapply(spec_mods_703_nut2, get_coef)
spec_703_all_nut2 = bind_rows(spec_703_all_nut2, .id = "column_label")
write.csv(spec_703_all_nut2, "results/spec_703_coefs_nut2.csv")

spec_struc_703_all_nut2 <- lapply(spec_struc_mods_703_nut2, get_coef)
spec_struc_703_all_nut2 = bind_rows(spec_struc_703_all_nut2, .id = "column_label")
write.csv(spec_struc_703_all_nut2, "results/spec_struc_703_coefs_nut2.csv")



struc_717_all_nut1 <- lapply(struc_mods_717_nut1, get_coef)
struc_717_all_nut1 = bind_rows(struc_717_all_nut1, .id = "column_label")
write.csv(struc_717_all_nut1, "results/struc_717_coefs_nut1.csv")

spec_717_all_nut1 <- lapply(spec_mods_717_nut1, get_coef)
spec_717_all_nut1 = bind_rows(spec_717_all_nut1, .id = "column_label")
write.csv(spec_717_all_nut1, "results/spec_717_coefs_nut1.csv")

spec_struc_717_all_nut1 <- lapply(spec_struc_mods_717_nut1, get_coef)
spec_struc_717_all_nut1 = bind_rows(spec_struc_717_all_nut1, .id = "column_label")
write.csv(spec_struc_717_all_nut1, "results/spec_struc_717_coefs_nut1.csv")


struc_717_all_nut2 <- lapply(struc_mods_717_nut2, get_coef)
struc_717_all_nut2 = bind_rows(struc_717_all_nut2, .id = "column_label")
write.csv(struc_717_all_nut2, "results/struc_717_coefs_nut2.csv")

spec_717_all_nut2 <- lapply(spec_mods_717_nut2, get_coef)
spec_717_all_nut2 = bind_rows(spec_717_all_nut2, .id = "column_label")
write.csv(spec_717_all_nut2, "results/spec_717_coefs_nut2.csv")

spec_struc_717_all_nut2 <- lapply(spec_struc_mods_717_nut2, get_coef)
spec_struc_717_all_nut2 = bind_rows(spec_struc_717_all_nut2, .id = "column_label")
write.csv(spec_struc_717_all_nut2, "results/spec_struc_717_coefs_nut2.csv")



struc_717_all_nut3 <- lapply(struc_mods_717_nut3, get_coef)
struc_717_all_nut3 = bind_rows(struc_717_all_nut3, .id = "column_label")
write.csv(struc_717_all_nut3, "results/struc_717_coefs_nut3.csv")

spec_717_all_nut3 <- lapply(spec_mods_717_nut3, get_coef)
spec_717_all_nut3 = bind_rows(spec_717_all_nut3, .id = "column_label")
write.csv(spec_717_all_nut3, "results/spec_717_coefs_nut3.csv")

spec_struc_717_all_nut3 <- lapply(spec_struc_mods_717_nut3, get_coef)
spec_struc_717_all_nut3 = bind_rows(spec_struc_717_all_nut3, .id = "column_label")
write.csv(spec_struc_717_all_nut3, "results/spec_struc_717_coefs_nut3.csv")



#-------------------------------------------- by date plots --------------------------------------------------
struc_619_nut1$date <- "June 19"
struc_703_nut1$date <- "July 3"
struc_703_nut2$date <- "July 3"
struc_717_nut1$date <- "July 17"
struc_717_nut2$date <- "July 17"
struc_717_nut3$date <- "July 17"

struc_619_nut1$nutrient <- 1
struc_703_nut1$nutrient <- 1
struc_703_nut2$nutrient <- 2
struc_717_nut1$nutrient <- 1
struc_717_nut2$nutrient <- 2
struc_717_nut3$nutrient <- 3

spec_619_nut1$date <- "June 19"
spec_703_nut1$date <- "July 3"
spec_703_nut2$date <- "July 3"
spec_717_nut1$date <- "July 17"
spec_717_nut2$date <- "July 17"
spec_717_nut3$date <- "July 17"

spec_619_nut1$nutrient <- 1
spec_703_nut1$nutrient <- 1
spec_703_nut2$nutrient <- 2
spec_717_nut1$nutrient <- 1
spec_717_nut2$nutrient <- 2
spec_717_nut3$nutrient <- 3

spec_struc_619_nut1$date <- "June 19"
spec_struc_703_nut1$date <- "July 3"
spec_struc_703_nut2$date <- "July 3"
spec_struc_717_nut1$date <- "July 17"
spec_struc_717_nut2$date <- "July 17"
spec_struc_717_nut3$date <- "July 17"

spec_struc_619_nut1$nutrient <- 1
spec_struc_703_nut1$nutrient <- 1
spec_struc_703_nut2$nutrient <- 2
spec_struc_717_nut1$nutrient <- 1
spec_struc_717_nut2$nutrient <- 2
spec_struc_717_nut3$nutrient <- 3


#619 nut 1

struc_619_nut1$type <- "structural"
spec_619_nut1$type <- "spectral"
spec_struc_619_nut1$type <- "spectral+structural"
all_619_nut1 <- rbind(struc_619_nut1[,13:17], spec_619_nut1[,8:12], spec_struc_619_nut1[,15:19])


all_619_nut1$yvar <- factor(all_619_nut1$yvar,levels = c("N","P","K","B"))
all_619_nut1$type <- factor(all_619_nut1$type, levels = c("structural", "spectral", "spectral+structural"))

g1 <- ggplot(data=all_619_nut1, aes(x=yvar, y=adj_r2, fill=type)) +
  geom_bar(stat="identity", width=0.7, position=position_dodge(width=0.8)) + ylim(0,0.85) +theme_hc() + 
  theme(legend.position = "none") +theme(
    plot.title = element_blank(),
    axis.title.x = element_blank(),
    axis.title.y = element_blank(), text = element_text(size=20, family="CM Roman"))

#703 nut 1

struc_703_nut1$type <- "structural"
spec_703_nut1$type <- "spectral"
spec_struc_703_nut1$type <- "spectral+structural"
all_703_nut1 <- rbind(struc_703_nut1[,13:17], spec_703_nut1[,8:12], spec_struc_703_nut1[,15:19])


all_703_nut1$yvar <- factor(all_703_nut1$yvar,levels = c("N","P","K","B"))
all_703_nut1$type <- factor(all_703_nut1$type, levels = c("structural", "spectral", "spectral+structural"))

g2 <- ggplot(data=all_703_nut1, aes(x=yvar, y=adj_r2, fill=type)) + labs(y=expression("adjusted" ~ R^{2})) +
  geom_bar(stat="identity", width=0.7, position=position_dodge(width=0.8)) + ylim(0,0.85) +theme_hc() + 
  theme(legend.position = "none") +theme(
    plot.title = element_blank(),
    axis.title.x = element_blank(),
    axis.title.y = element_blank(), text = element_text(size=20, family="CM Roman"))


#703 nut 2

struc_703_nut2$type <- "structural"
spec_703_nut2$type <- "spectral"
spec_struc_703_nut2$type <- "spectral+structural"
all_703_nut2 <- rbind(struc_703_nut2[,13:17], spec_703_nut2[,8:12], spec_struc_703_nut2[,15:19])


all_703_nut2$yvar <- factor(all_703_nut2$yvar,levels = c("N","P","K","B"))
all_703_nut2$type <- factor(all_703_nut2$type, levels = c("structural", "spectral", "spectral+structural"))

g3 <- ggplot(data=all_703_nut2, aes(x=yvar, y=adj_r2, fill=type)) +
  geom_bar(stat="identity", width=0.7, position=position_dodge(width=0.8)) + ylim(0,0.85) +theme_hc() + 
  theme(legend.position = "none") +theme(
    plot.title = element_blank(),
    axis.title.x = element_blank(),
    axis.title.y = element_blank(), text = element_text(size=20, family="CM Roman"))
#717 nut 1

struc_717_nut1$type <- "structural"
spec_717_nut1$type <- "spectral"
spec_struc_717_nut1$type <- "spectral+structural"
all_717_nut1 <- rbind(struc_717_nut1[,13:17], spec_717_nut1[,8:12], spec_struc_717_nut1[,15:19])


all_717_nut1$yvar <- factor(all_717_nut1$yvar,levels = c("N","P","K","B"))
all_717_nut1$type <- factor(all_717_nut1$type, levels = c("structural", "spectral", "spectral+structural"))

g4 <- ggplot(data=all_717_nut1, aes(x=yvar, y=adj_r2, fill=type)) +
  geom_bar(stat="identity", width=0.7, position=position_dodge(width=0.8)) + ylim(0,0.85) +theme_hc() + 
  theme(legend.position = "none") +theme(
    plot.title = element_blank(),
    axis.title.x = element_blank(),
    axis.title.y = element_blank(), text = element_text(size=20, family="CM Roman"))

#717 nut 2

struc_717_nut2$type <- "structural"
spec_717_nut2$type <- "spectral"
spec_struc_717_nut2$type <- "spectral+structural"
all_717_nut2 <- rbind(struc_717_nut2[,13:17], spec_717_nut2[,8:12], spec_struc_717_nut2[,15:19])


all_717_nut2$yvar <- factor(all_717_nut2$yvar,levels = c("N","P","K","B"))
all_717_nut2$type <- factor(all_717_nut2$type, levels = c("structural", "spectral", "spectral+structural"))

g5 <- ggplot(data=all_717_nut2, aes(x=yvar, y=adj_r2, fill=type)) +
  geom_bar(stat="identity", width=0.7, position=position_dodge(width=0.8)) + ylim(0,0.85) +theme_hc() + 
  theme(legend.position = "none") +theme(
    plot.title = element_blank(),
    axis.title.x = element_blank(),
    axis.title.y = element_blank(), text = element_text(size=20, family="CM Roman"))

#717 nut 3

struc_717_nut3$type <- "structural"
spec_717_nut3$type <- "spectral"
spec_struc_717_nut3$type <- "spectral+structural"
all_717_nut3 <- rbind(struc_717_nut3[,13:17], spec_717_nut3[,8:12], spec_struc_717_nut3[,15:19])


all_717_nut3$yvar <- factor(all_717_nut3$yvar,levels = c("N","P","K","B"))
all_717_nut3$type <- factor(all_717_nut3$type, levels = c("structural", "spectral", "spectral+structural"))

g6 <- ggplot(data=all_717_nut3, aes(x=yvar, y=adj_r2, fill=type)) +
  geom_bar(stat="identity", width=0.7, position=position_dodge(width=0.8)) + ylim(0,0.85) + theme_hc() + 
  theme(legend.position = "none") +theme(
    plot.title = element_blank(),
    axis.title.x = element_blank(),
    axis.title.y = element_blank(), text = element_text(size=20, family="CM Roman"))

#xlab('nutrient') + labs(y=expression(paste("adjusted r"^"2"))) + ylim(0,0.85) + ggtitle("June 19")

# All dates
all_dates <- rbind(all_619_nut1, all_703_nut1, all_703_nut2, all_717_nut1, all_717_nut2, all_717_nut3)
all_dates$date <- factor(all_dates$date, levels = c("June 19", "July 3", "July 17"))

ggplot(data=all_dates, aes(x=yvar, y=adj_r2, fill=type)) +
  geom_bar(stat="identity", width=0.6, position=position_dodge(width=0.8))+
  xlab('response variable') + labs(y=expression(paste("adjusted R "^"2"))) + ylim(0,0.81) + facet_wrap(~date+nutrient) +
  theme(strip.text.x = element_text(size = 10, face = "bold"), text = element_text(size=10, family="CM Roman")) + theme_hc() + 
  scale_fill_discrete(name = "", labels = c("stuctural metrics", "spectral metrics", "spectral + structural metrics"))

ylab <- expression(paste("adjusted R "^"2"))
grid.arrange(g1, g2, g3, g4, g5, g6, layout_matrix = rbind(c(1, NA, NA), c(2,3, NA), c(4,5,6)))

#-------------------------------------------- by nutrient plots --------------------------------------------------
# 
# #N
# N_619 <- rbind(struc_619[1,13:15], spec_619[1,8:10], spec_struc_703[1,15:17])
# N_619$Date <- "June 13"
# N_703 <- rbind(struc_703[1,13:15], spec_703[1,8:10], spec_struc_703[1,15:17])
# N_703$Date <- "July 3"
# N_717 <- rbind(struc_717[1,13:15], spec_717[1,8:10], spec_struc_717[1,15:17])
# N_717$Date <- "July 17"
# N_all <- rbind(N_613, N_703, N_717)
# 
# N_all$Date <- factor(N_all$Date, levels = c("June 13", "July 3", "July 17"))
# N_all$type <- factor(N_all$type, levels = c("structural", "spectral", "spectral+structural"))
# 
# ggplot(data=N_all, aes(x=Date, y=adj_r2, fill=type)) +
#   geom_bar(stat="identity", width=0.7, position=position_dodge(width=0.8))+
#   xlab('date') + labs(y=expression(paste("adjusted r"^"2"))) + ylim(0,0.85) + ggtitle("Nitrogen")
# 
# 
# #P
# P_613 <- rbind(struc_619[2,13:15], spec_619[2,8:10], spec_struc_703[2,15:17])
# P_613$Date <- "June 13"
# P_703 <- rbind(struc_703[2,13:15], spec_703[2,8:10], spec_struc_703[2,15:17])
# P_703$Date <- "July 3"
# P_717 <- rbind(struc_717[2,13:15], spec_717[2,8:10], spec_struc_717[2,15:17])
# P_717$Date <- "July 17"
# P_all <- rbind(P_613, P_703, P_717)
# 
# P_all$Date <- factor(P_all$Date, levels = c("June 13", "July 3", "July 17"))
# P_all$type <- factor(P_all$type, levels = c("structural", "spectral", "spectral+structural"))
# 
# ggplot(data=P_all, aes(x=Date, y=adj_r2, fill=type)) +
#   geom_bar(stat="identity", width=0.7, position=position_dodge(width=0.8))+
#   xlab('date') + labs(y=expression(paste("adjusted r"^"2"))) + ylim(0,0.85) + ggtitle("Phosphorous")
# 
# 
# #K
# K_613 <- rbind(struc_619[3,13:15], spec_619[3,8:10], spec_struc_703[3,15:17])
# K_613$Date <- "June 13"
# K_703 <- rbind(struc_703[3,13:15], spec_703[3,8:10], spec_struc_703[3,15:17])
# K_703$Date <- "July 3"
# K_717 <- rbind(struc_717[3,13:15], spec_717[3,8:10], spec_struc_717[3,15:17])
# K_717$Date <- "July 17"
# K_all <- rbind(K_613, K_703, K_717)
# 
# K_all$Date <- factor(K_all$Date, levels = c("June 13", "July 3", "July 17"))
# K_all$type <- factor(K_all$type, levels = c("structural", "spectral", "spectral+structural"))
# 
# ggplot(data=K_all, aes(x=Date, y=adj_r2, fill=type)) +
#   geom_bar(stat="identity", width=0.7, position=position_dodge(width=0.8))+
#   xlab('date') + labs(y=expression(paste("adjusted r"^"2"))) + ylim(0,0.85) + ggtitle("Potassium")
# 
# 
# #Mg
# Mg_613 <- rbind(struc_619[4,13:15], spec_619[4,8:10], spec_struc_703[4,15:17])
# Mg_613$Date <- "June 13"
# Mg_703 <- rbind(struc_703[4,13:15], spec_703[4,8:10], spec_struc_703[4,15:17])
# Mg_703$Date <- "July 3"
# Mg_717 <- rbind(struc_717[4,13:15], spec_717[4,8:10], spec_struc_717[4,15:17])
# Mg_717$Date <- "July 17"
# Mg_all <- rbind(Mg_613, Mg_703, Mg_717)
# 
# Mg_all$Date <- factor(Mg_all$Date, levels = c("June 13", "July 3", "July 17"))
# Mg_all$type <- factor(Mg_all$type, levels = c("structural", "spectral", "spectral+structural"))
# 
# ggplot(data=Mg_all, aes(x=Date, y=adj_r2, fill=type)) +
#   geom_bar(stat="identity", width=0.7, position=position_dodge(width=0.8))+
#   xlab('date') + labs(y=expression(paste("adjusted r"^"2"))) + ylim(0,0.85) + ggtitle("Magnesium")
# 
# 
# #Ca
# Ca_613 <- rbind(struc_619[5,13:15], spec_619[5,8:10], spec_struc_703[5,15:17])
# Ca_613$Date <- "June 13"
# Ca_703 <- rbind(struc_703[5,13:15], spec_703[5,8:10], spec_struc_703[5,15:17])
# Ca_703$Date <- "July 3"
# Ca_717 <- rbind(struc_717[5,13:15], spec_717[5,8:10], spec_struc_717[5,15:17])
# Ca_717$Date <- "July 17"
# Ca_all <- rbind(Ca_613, Ca_703, Ca_717)
# 
# Ca_all$Date <- factor(Ca_all$Date, levels = c("June 13", "July 3", "July 17"))
# Ca_all$type <- factor(Ca_all$type, levels = c("structural", "spectral", "spectral+structural"))
# 
# ggplot(data=Ca_all, aes(x=Date, y=adj_r2, fill=type)) +
#   geom_bar(stat="identity", width=0.7, position=position_dodge(width=0.8))+
#   xlab('date') + labs(y=expression(paste("adjusted r"^"2"))) + ylim(0,0.85) + ggtitle("Calcium")
# 
# 
# #S
# S_613 <- rbind(struc_619[6,13:15], spec_619[6,8:10], spec_struc_703[6,15:17])
# S_613$Date <- "June 13"
# S_703 <- rbind(struc_703[6,13:15], spec_703[6,8:10], spec_struc_703[6,15:17])
# S_703$Date <- "July 3"
# S_717 <- rbind(struc_717[6,13:15], spec_717[6,8:10], spec_struc_717[6,15:17])
# S_717$Date <- "July 17"
# S_all <- rbind(S_613, S_703, S_717)
# 
# S_all$Date <- factor(S_all$Date, levels = c("June 13", "July 3", "July 17"))
# S_all$type <- factor(S_all$type, levels = c("structural", "spectral", "spectral+structural"))
# 
# ggplot(data=S_all, aes(x=Date, y=adj_r2, fill=type)) +
#   geom_bar(stat="identity", width=0.7, position=position_dodge(width=0.8))+
#   xlab('date') + labs(y=expression(paste("adjusted r"^"2"))) + ylim(0,0.85) + ggtitle("Sulfur")
# 
# 
# #B
# B_613 <- rbind(struc_619[10,13:15], spec_619[10,8:10], spec_struc_703[10,15:17])
# B_613$Date <- "June 13"
# B_703 <- rbind(struc_703[10,13:15], spec_703[10,8:10], spec_struc_703[10,15:17])
# B_703$Date <- "July 3"
# B_717 <- rbind(struc_717[10,13:15], spec_717[10,8:10], spec_struc_717[10,15:17])
# B_717$Date <- "July 17"
# B_all <- rbind(B_613, B_703, B_717)
# 
# B_all$Date <- factor(B_all$Date, levels = c("June 13", "July 3", "July 17"))
# B_all$type <- factor(B_all$type, levels = c("structural", "spectral", "spectral+structural"))
# 
# ggplot(data=B_all, aes(x=Date, y=adj_r2, fill=type)) +
#   geom_bar(stat="identity", width=0.7, position=position_dodge(width=0.8))+
#   xlab('date') + labs(y=expression(paste("adjusted r"^"2"))) + ylim(0,0.85) + ggtitle("Boron")

# -----------------------------tables-----------------------------

# Regression model tables for 619 nutrients ONLY

n_struc_mods <- list(struc_mods_619[[1]], struc_mods_703[[1]], struc_mods_717[[1]])
stargazer(n_struc_mods, omit.stat = c("f"), order = xvars_struc_fixed, single.row = F, no.space = T, align = T,
          covariate.labels = c("crop height mean (m)", "crop height histogram kurtosis", "Rumple index", "canopy relief ratio median", 
                               "canopy relief ratio IQR", "Moran's I", "water accumulation (m)"))

n_spec_mods <- list(spec_mods_619[[1]], spec_mods_703[[1]], spec_mods_717[[1]])
stargazer(n_spec_mods, omit.stat = c("f"), order = xvars_spec_fixed, single.row = F, no.space = T, align = T,
          covariate.labels = c("TGI mean", "TGI standard deviation"))

n_spec_struc_mods <- list(spec_struc_mods_619[[1]], spec_struc_mods_703[[1]], spec_struc_mods_717[[1]])
stargazer(n_spec_struc_mods, omit.stat = c("f"), order = xvars_all_fixed, single.row = F, no.space = T, align = T,
          covariate.labels = c("crop height mean (m)", "TGI mean", "crop height histogram kurtosis", "Rumple index", "canopy relief ratio median", 
                               "Moran's I", "water accumulation (m)", "TGI standard deviation"))


p_struc_mods <- list(struc_mods_619[[2]], struc_mods_703[[2]], struc_mods_717[[2]])
stargazer(p_struc_mods, omit.stat = c("f"), order = xvars_struc_fixed, single.row = F, no.space = T, align = T,
          covariate.labels = c("crop height mean (m)", "crop height histogram kurtosis", "canopy relief ratio median"))

p_spec_mods <- list(spec_mods_619[[2]], spec_mods_703[[2]], spec_mods_717[[2]])
stargazer(p_spec_mods, omit.stat = c("f"), order = xvars_spec_fixed, single.row = F, no.space = T, align = T,
          covariate.labels = c("TGI mean", "TGI standard deviation"))

p_spec_struc_mods <- list(spec_struc_mods_619[[2]], spec_struc_mods_703[[2]], spec_struc_mods_717[[2]])
stargazer(p_spec_struc_mods, omit.stat = c("f"), order = xvars_all_fixed, single.row = F, no.space = T, align = T,
          covariate.labels = c("crop height mean (m)", "TGI mean", "crop height histogram kurtosis", "Rumple index",
                               "canopy relief ratio IQR", "water accumulation (m)", "TGI standard deviation"))

k_struc_mods <- list(struc_mods_619[[3]], struc_mods_703[[3]], struc_mods_717[[3]])
stargazer(k_struc_mods, omit.stat = c("f"), order = xvars_struc_fixed, single.row = F, no.space = T, align = T,
          covariate.labels = c("Rumple index", "canopy relief ratio median"))

k_spec_mods <- list(spec_mods_619[[3]], spec_mods_703[[3]], spec_mods_717[[3]])
stargazer(k_spec_mods, omit.stat = c("f"), order = xvars_spec_fixed, single.row = F, no.space = T, align = T,
          covariate.labels = c("TGI mean", "TGI standard deviation"))

k_spec_struc_mods <- list(spec_struc_mods_619[[3]], spec_struc_mods_703[[3]], spec_struc_mods_717[[3]])
stargazer(k_spec_struc_mods, omit.stat = c("f"), order = xvars_all_fixed, single.row = F, no.space = T, align = T,
          covariate.labels = c("TGI mean", "canopy relief ratio median",
                               "canopy relief ratio IQR","Moran's I", "water accumulation (m)", "TGI standard deviation"))

b_struc_mods <- list(struc_mods_619[[10]], struc_mods_703[[10]], struc_mods_717[[10]])
stargazer(b_struc_mods, omit.stat = c("f"), order = xvars_struc_fixed, single.row = F, no.space = T, align = T,
          covariate.labels = c("crop height mean (m)", "crop height histogram kurtosis", "Rumple index", "canopy relief ratio median",
                               "canopy relief ratio IQR","Moran's I"))

b_spec_mods <- list(spec_mods_619[[10]], spec_mods_703[[10]], spec_mods_717[[10]])
stargazer(b_spec_mods, omit.stat = c("f"), order = xvars_spec_fixed, single.row = F, no.space = T, align = T,
          covariate.labels = c("TGI mean", "TGI standard deviation"))

b_spec_struc_mods <- list(spec_struc_mods_619[[10]], spec_struc_mods_703[[10]], spec_struc_mods_717[[10]])
stargazer(b_spec_struc_mods, omit.stat = c("f"), order = xvars_all_fixed, single.row = F, no.space = T, align = T,
          covariate.labels = c("crop height mean (m)", "TGI mean", "crop height histogram kurtosis", "Rumple index","canopy relief ratio median",
                               "canopy relief ratio IQR","Moran's I"))



# Regression model tables USING NUTRIENTS BY DATE

n_struc_mods <- list(struc_mods_619[[1]], struc_mods_703[[1]], struc_mods_717[[1]])
stargazer(n_struc_mods, omit.stat = c("f"), order = xvars_struc_fixed, single.row = F, no.space = T, align = T,
          covariate.labels = c("crop height mean (m)", "crop height histogram kurtosis", "Rumple index", "canopy relief ratio median", 
                               "canopy relief ratio IQR", "Moran's I", "water accumulation (m)"))

n_spec_mods <- list(spec_mods_619[[1]], spec_mods_703[[1]], spec_mods_717[[1]])
stargazer(n_spec_mods, omit.stat = c("f"), order = xvars_spec_fixed, single.row = F, no.space = T, align = T,
          covariate.labels = c("TGI mean", "TGI standard deviation"))

n_spec_struc_mods <- list(spec_struc_mods_619[[1]], spec_struc_mods_703[[1]], spec_struc_mods_717[[1]])
stargazer(n_spec_struc_mods, omit.stat = c("f"), order = xvars_all_fixed, single.row = F, no.space = T, align = T,
          covariate.labels = c("TGI mean", "Rumple index", "canopy relief ratio median", 
                               "Moran's I", "water accumulation (m)", "TGI standard deviation"))


p_struc_mods <- list(struc_mods_619[[2]], struc_mods_703[[2]], struc_mods_717[[2]])
stargazer(p_struc_mods, omit.stat = c("f"), order = xvars_struc_fixed, single.row = F, no.space = T, align = T,
          covariate.labels = c("crop height mean (m)", "crop height histogram kurtosis", "canopy relief ratio median"))

p_spec_mods <- list(spec_mods_619[[2]], spec_mods_703[[2]], spec_mods_717[[2]])
stargazer(p_spec_mods, omit.stat = c("f"), order = xvars_spec_fixed, single.row = F, no.space = T, align = T,
          covariate.labels = c("TGI mean", "TGI standard deviation"))

p_spec_struc_mods <- list(spec_struc_mods_619[[2]], spec_struc_mods_703[[2]], spec_struc_mods_717[[2]])
stargazer(p_spec_struc_mods, omit.stat = c("f"), order = xvars_all_fixed, single.row = F, no.space = T, align = T,
          covariate.labels = c("crop height mean (m)", "canopy relief ratio median",
                               "canopy relief ratio IQR", "water accumulation (m)", "TGI standard deviation"))

k_struc_mods <- list(struc_mods_619[[3]], struc_mods_703[[3]], struc_mods_717[[3]])
stargazer(k_struc_mods, omit.stat = c("f"), order = xvars_struc_fixed, single.row = F, no.space = T, align = T,
          covariate.labels = c("Rumple index", "canopy relief ratio median", "canopy relief ratio IQR", "Moran's I"))

k_spec_mods <- list(spec_mods_619[[3]], spec_mods_703[[3]], spec_mods_717[[3]])
stargazer(k_spec_mods, omit.stat = c("f"), order = xvars_spec_fixed, single.row = F, no.space = T, align = T,
          covariate.labels = c("TGI mean", "TGI standard deviation"))

k_spec_struc_mods <- list(spec_struc_mods_619[[3]], spec_struc_mods_703[[3]], spec_struc_mods_717[[3]])
stargazer(k_spec_struc_mods, omit.stat = c("f"), order = xvars_all_fixed, single.row = F, no.space = T, align = T,
          covariate.labels = c("TGI mean", "canopy relief ratio median",
                               "canopy relief ratio IQR","Moran's I", "water accumulation (m)", "TGI standard deviation"))

b_struc_mods <- list(struc_mods_619[[10]], struc_mods_703[[10]], struc_mods_717[[10]])
stargazer(b_struc_mods, omit.stat = c("f"), order = xvars_struc_fixed, single.row = F, no.space = T, align = T,
          covariate.labels = c("crop height mean (m)", "Rumple index", "canopy relief ratio median"))

b_spec_mods <- list(spec_mods_619[[10]], spec_mods_703[[10]], spec_mods_717[[10]])
stargazer(b_spec_mods, omit.stat = c("f"), order = xvars_spec_fixed, single.row = F, no.space = T, align = T)

b_spec_struc_mods <- list(spec_struc_mods_619[[10]], spec_struc_mods_703[[10]], spec_struc_mods_717[[10]])
stargazer(b_spec_struc_mods, omit.stat = c("f"), order = xvars_all_fixed, single.row = F, no.space = T, align = T,
          covariate.labels = c("TGI mean", "crop height mean (m)", "Rumple index","canopy relief ratio median",
                               "Moran's I",  "water accumulation (m)"))


# Nutrient summary tables

stargazer(agg_619[,c(2:4, 12)])
stargazer(agg_703[,c(2:4, 12)])
stargazer(agg_717[,c(2:4, 12)])

# Metric summary tables

stargazer(plots_619@data[,xvars_all_fixed])
stargazer(plots_703@data[,xvars_all_fixed])
stargazer(plots_717@data[,xvars_all_fixed])

stargazer(plots_619@data[,23:67])
stargazer(plots_703@data[,23:67])
stargazer(plots_717@data[,23:67])


# ------------------------------------ model summary chart -------------------------------------------------------------

# structural


struc_619_nut1 <- as.data.table(struc_619_nut1)
struc_703_nut1 <- as.data.table(struc_703_nut1)
struc_703_nut2 <- as.data.table(struc_703_nut2)
struc_717_nut1 <- as.data.table(struc_717_nut1)
struc_717_nut2 <- as.data.table(struc_717_nut2)
struc_717_nut3 <- as.data.table(struc_717_nut3)

struc_all <- rbind(struc_619_nut1, struc_703_nut1, struc_703_nut2, struc_717_nut1, struc_717_nut2, struc_717_nut3, fill = T)
struc_all <- struc_all[struc_all$yvar %in% c("N", "P", "K", "B"),]
struc_all$model <- paste(struc_all$date, "_", struc_all$nutrient)

struc_N <- struc_all[yvar == "N",][,c(1:7,13:14,18)]
struc_P <- struc_all[yvar == "P",][,c(1:7,13:14,18)]
struc_K <- struc_all[yvar == "K",][,c(1:7,13:14,18)]
struc_B <- struc_all[yvar == "B",][,c(1:7,13:14,18)]

write.csv(struc_N, "results/struc_N.csv")
write.csv(struc_P, "results/struc_P.csv")
write.csv(struc_K, "results/struc_K.csv")
write.csv(struc_B, "results/struc_B.csv")

# mycolor=gradientColor(low="red",mid="white",high="blue",n=20,plot=TRUE)
# ztable(struc_N[,c(1:7)]) %>% 
#   makeHeatmap(mycolor = mycolor) %>%
#   print(caption="Table 5. Heatmap table with 'Blue' palette")


# struc_N_heatmap <- data.frame(as.vector(t(struc_N[,c(1:7)])))
# struc_N_heatmap$model <- c(rep("June 19 _ 1", 7), rep("July 3 _ 1", 7), rep("July 3 _ 2", 7), rep("July 17 _ 1", 7), rep("July 17 _ 2", 7), rep("July 17 _ 3", 7))
# struc_N_heatmap$model <- factor(struc_N_heatmap$model,levels = rev(c("June 19 _ 1","July 3 _ 1","July 3 _ 2", "July 17 _ 1", "July 17 _ 2", "July 17 _ 3")))
# struc_N_heatmap$xvar <- names(struc_N)[1:7]
# struc_N_heatmap$xvar <- factor(struc_N_heatmap$xvar, levels = c("mean", "kurt", "rumple", "crrmedian", "crriqr", "moran", "watersum"))
# names(struc_N_heatmap)[1] <- "coef"
# 
# library(tikzDevice)
# 
# tikz(file = "plot_test.tex")
# library(ggplot2)
# 
# plot <- ggplot(struc_N_heatmap, aes(xvar, model)) + # x and y axes => Var1 and Var2
#   geom_tile(aes(fill = coef)) + # background colours are mapped according to the value column
#   geom_text(aes(fill = struc_N_heatmap$coef, label = round(struc_N_heatmap$coef, 4))) + # write the values
#   scale_fill_gradient2(low = muted("darkred"), 
#                        mid = "white", 
#                        high = muted("midnightblue"), 
#                        midpoint = 0)
# 
# print(plot)
# 
# #Necessary to close or the tikxDevice .tex file will not be written
# dev.off()
# 
# library(condformat)
# library(formattable)
# 
# condformat(struc_N[,c(1:7)]) %>%
# rule_fill_gradient2(crriqr, midpoint = 0) %>%
#   rule_fill_gradient2(crrmedian, midpoint = 0)%>%
#   rule_fill_gradient2(kurt, midpoint = 0)%>%
#   rule_fill_gradient2(mean, midpoint = 0)%>%
#   rule_fill_gradient2(moran, midpoint = 0)%>%
#   rule_fill_gradient2(rumple, midpoint = 0)%>%
#   rule_fill_gradient2(watersum, midpoint = 0)
# 
# condformat(struc_N[,c(1:7)]) %>%
#   rule_fill_gradient2(xvars_struc_fixed, expression = , midpoint = 0)
  


# spectral

spec_619_nut1 <- as.data.table(spec_619_nut1)
spec_703_nut1 <- as.data.table(spec_703_nut1)
spec_703_nut2 <- as.data.table(spec_703_nut2)
spec_717_nut1 <- as.data.table(spec_717_nut1)
spec_717_nut2 <- as.data.table(spec_717_nut2)
spec_717_nut3 <- as.data.table(spec_717_nut3)

spec_all <- rbind(spec_619_nut1, spec_703_nut1, spec_703_nut2, spec_717_nut1, spec_717_nut2, spec_717_nut3, fill = T)
spec_all <- spec_all[spec_all$yvar %in% c("N", "P", "K", "B"),]
spec_all$model <- paste(spec_all$date, "_", spec_all$nutrient)

spec_N <- spec_all[yvar == "N",][,c(1:2,8,9,13)]
spec_P <- spec_all[yvar == "P",][,c(1:2,8,9,13)]
spec_K <- spec_all[yvar == "K",][,c(1:2,8,9,13)]
spec_B <- spec_all[yvar == "B",][,c(1:2,8,9,13)]


write.csv(spec_N, "results/spec_N.csv")
write.csv(spec_P, "results/spec_P.csv")
write.csv(spec_K, "results/spec_K.csv")
write.csv(spec_B, "results/spec_B.csv")


# spectral + structural

spec_struc_619_nut1 <- as.data.table(spec_struc_619_nut1)
spec_struc_703_nut1 <- as.data.table(spec_struc_703_nut1)
spec_struc_703_nut2 <- as.data.table(spec_struc_703_nut2)
spec_struc_717_nut1 <- as.data.table(spec_struc_717_nut1)
spec_struc_717_nut2 <- as.data.table(spec_struc_717_nut2)
spec_struc_717_nut3 <- as.data.table(spec_struc_717_nut3)

spec_struc_all <- rbind(spec_struc_619_nut1, spec_struc_703_nut1, spec_struc_703_nut2, spec_struc_717_nut1, spec_struc_717_nut2, spec_struc_717_nut3, fill = T)
spec_struc_all <- spec_struc_all[spec_struc_all$yvar %in% c("N", "P", "K", "B"),]
spec_struc_all$model <- paste(spec_struc_all$date, "_", spec_struc_all$nutrient)

spec_struc_N <- spec_struc_all[yvar == "N",][,c(1:9,15,16,20)]
spec_struc_P <- spec_struc_all[yvar == "P",][,c(1:9,15,16,20)]
spec_struc_K <- spec_struc_all[yvar == "K",][,c(1:9,15,16,20)]
spec_struc_B <- spec_struc_all[yvar == "B",][,c(1:9,15,16,20)]


write.csv(spec_struc_N, "results/spec_struc_N.csv")
write.csv(spec_struc_P, "results/spec_struc_P.csv")
write.csv(spec_struc_K, "results/spec_struc_K.csv")
write.csv(spec_struc_B, "results/spec_struc_B.csv")

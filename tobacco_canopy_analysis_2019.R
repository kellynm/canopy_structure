library(raster)
library(rgdal)
library(data.table)
library(plyr)
library(lattice)
library(gridExtra)
library(velox)
library(rgeos)
library(gstat)
library(usdm)
library(viridis)
library(ggpubr)
library(foreign)
library(car)
library(PerformanceAnalytics)
library(rcompanion)
library(olsrr)
library(fastDummies)
library(tidyr)
library(graphics)
library(lme4)
library(MuMIn)
library(spdep)
#Load raster data for all dates
#setwd("/media/Kellyn/F20E17B40E177139/kpmontgo@ncsu.edu/LkWheeler_Sorghum/LkWheeler_Fusarium_Sorghum")
setwd("Q:/My Drive/Research/Canopy_Morphology/Tobacco_Project/2019")

tobacco_area <- readOGR("Q:/My Drive/Research/Canopy_Morphology/Tobacco_Project/2019/layers/boundary", "wilson19_boundary", stringsAsFactors = F)
tobacco_area <- spTransform(tobacco_area, CRS("+proj=lcc +lat_1=36.16666666666666 +lat_2=34.33333333333334 +lat_0=33.75 +lon_0=-79 +x_0=609601.22 +y_0=0 +datum=NAD83 +units=m +no_defs +ellps=GRS80 +towgs84=0,0,0"))

csm_619 <- crop(raster('Q:/My Drive/Research/Canopy_Morphology/Tobacco_Project/2019/csm/wilson19_619_csm.tif'), tobacco_area)
csm_703 <- crop(raster('Q:/My Drive/Research/Canopy_Morphology/Tobacco_Project/2019/csm/wilson19_703_csm.tif'), tobacco_area)
csm_717 <- crop(raster('Q:/My Drive/Research/Canopy_Morphology/Tobacco_Project/2019/csm/wilson19_717_csm.tif'), tobacco_area)
dem <- crop(raster('Q:/My Drive/Research/Canopy_Morphology/Tobacco_Project/2019/lidar/wilson19_dem.tif'), tobacco_area)

nosoil_619 <- crop(raster('Q:/My Drive/Research/Canopy_Morphology/Tobacco_Project/2019/nosoil/csm_nosoil_619.tif'), tobacco_area)
nosoil_703 <- crop(raster('Q:/My Drive/Research/Canopy_Morphology/Tobacco_Project/2019/nosoil/csm_nosoil_703.tif'), tobacco_area)
nosoil_717 <- crop(raster('Q:/My Drive/Research/Canopy_Morphology/Tobacco_Project/2019/nosoil/csm_nosoil_717.tif'), tobacco_area)

csm_masked_619 <- mask(csm_619, nosoil_619)
csm_masked_703 <- mask(csm_619, nosoil_703)
csm_masked_717 <- mask(csm_619, nosoil_717)

csm_619_velox <- velox(csm_masked_619)
csm_703_velox <- velox(csm_masked_703)
csm_717_velox <- velox(csm_masked_717)

plots <- readOGR("Q:/My Drive/Research/Canopy_Morphology/Tobacco_Project/2019/layers/plots", "wilson19_plots", stringsAsFactors = F)
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

nutrient <- fread("Q:/My Drive/Research/Canopy_Morphology/Tobacco_Project/2019/wilson19_nutrients.csv", stringsAsFactors = F)
#names(nutrient) <- c("Plot","treatment", "Sample #",  "N","P","K","Mg", "Ca","S","Zn", "Mn", "Cu", "Fe", "B", "Al", "Na") 

plots_619 <- plots
plots_619 <- merge(plots_619, nutrient[Date==617], by.x="join_id", by.y='Join_ID')
plots_703 <- plots
plots_703 <- merge(plots_703, nutrient[Date==617], by.x="join_id", by.y='Join_ID')
plots_717 <- plots
plots_717 <- merge(plots_717, nutrient[Date==617], by.x="join_id", by.y='Join_ID')

spad_617 <- fread("Q:/My Drive/Research/Canopy_Morphology/Tobacco_Project/2019/wilson19_617_spad.csv", stringsAsFactors = F)
spad_701 <- fread("Q:/My Drive/Research/Canopy_Morphology/Tobacco_Project/2019/wilson19_617_spad.csv", stringsAsFactors = F)
spad_718 <- fread("Q:/My Drive/Research/Canopy_Morphology/Tobacco_Project/2019/wilson19_617_spad.csv", stringsAsFactors = F)

plots_619 <- merge(plots_619, spad_617, by.x="spad_id", by.y='JOIN')
plots_703 <- merge(plots_703, spad_701, by.x="spad_id", by.y='JOIN')
plots_717 <- merge(plots_717, spad_718, by.x="spad_id", by.y='JOIN')


red_619 <- crop(raster('Q:/My Drive/Research/Canopy_Morphology/Tobacco_Project/2019/ortho/wilson/wilson19_619_red_georef.tif'), tobacco_area)
green_619 <- crop(raster('Q:/My Drive/Research/Canopy_Morphology/Tobacco_Project/2019/ortho/wilson/wilson19_619_green_georef.tif'), tobacco_area)
blue_619 <- crop(raster('Q:/My Drive/Research/Canopy_Morphology/Tobacco_Project/2019/ortho/wilson/wilson19_619_blue_georef.tif'), tobacco_area)
red_703 <- crop(raster('Q:/My Drive/Research/Canopy_Morphology/Tobacco_Project/2019/ortho/wilson/wilson19_703_red_georef.tif'), tobacco_area)
green_703 <- crop(raster('Q:/My Drive/Research/Canopy_Morphology/Tobacco_Project/2019/ortho/wilson/wilson19_703_green_georef.tif'), tobacco_area)
blue_703 <- crop(raster('Q:/My Drive/Research/Canopy_Morphology/Tobacco_Project/2019/ortho/wilson/wilson19_703_blue_georef.tif'), tobacco_area)
red_717 <- crop(raster('Q:/My Drive/Research/Canopy_Morphology/Tobacco_Project/2019/ortho/wilson/wilson19_717_red_georef.tif'), tobacco_area)
green_717 <- crop(raster('Q:/My Drive/Research/Canopy_Morphology/Tobacco_Project/2019/ortho/wilson/wilson19_717_green_georef.tif'), tobacco_area)
blue_717 <- crop(raster('Q:/My Drive/Research/Canopy_Morphology/Tobacco_Project/2019/ortho/wilson/wilson19_717_blue_georef.tif'), tobacco_area)

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
varInd_703 <- (green_masked_703-red_masked_703)/(green_masked_703+red_masked_703-blue_masked_703)
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


plots_619 <- merge(plots_619, varInd_df_list[[1]], by.x="plotid", by.y="plots")
plots_703 <- merge(plots_703, varInd_df_list[[2]], by.x="plotid", by.y="plots")
plots_717 <- merge(plots_717, varInd_df_list[[3]], by.x="plotid", by.y="plots")

plots_619 <- merge(plots_619, tgi_df_list[[1]], by.x="plotid", by.y="plots")
plots_703 <- merge(plots_703, tgi_df_list[[2]], by.x="plotid", by.y="plots")
plots_717 <- merge(plots_717, tgi_df_list[[3]], by.x="plotid", by.y="plots")


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

plots_619 <- merge(plots_619, waterDF, by.x="plotid", by.y="plots")
plots_703 <- merge(plots_703, waterDF, by.x="plotid", by.y="plots")
plots_717 <- merge(plots_717, waterDF, by.x="plotid", by.y="plots")

#---------------------------------------------------------- Canopy relief ratio----------------------------------------------------------------
crr <- function(x){
  (mean(x, na.rm=T)-min(x, na.rm=T))/(max(x, na.rm=T)-min(x, na.rm=T))
}

crrRast_619 <- focal(csm_619, w=matrix(1/81, nc=9, nr=9), crr)
crrRast_masked_619 <- mask(crrRast_619, nosoil_619)
crrRast_703 <- focal(csm_703, w=matrix(1/81, nc=9, nr=9), crr)
crrRast_masked_703 <- mask(crrRast_703, nosoil_703)
crrRast_717 <- focal(csm_717, w=matrix(1/81, nc=9, nr=9), crr)
crrRast_masked_717 <- mask(crrRast_717, nosoil_717)
crrRast_list <- list(crrRast_619,crrRast_703,crrRast_717)

writeRaster(crrRast_619, "Q:/My Drive/Research/Canopy_Morphology/Tobacco_Project/2019/results/crr_619.tif", format="GTiff", overwrite = T)
writeRaster(crrRast_703, "Q:/My Drive/Research/Canopy_Morphology/Tobacco_Project/2019/results/crr_703.tif", format="GTiff", overwrite = T)
writeRaster(crrRast_717, "Q:/My Drive/Research/Canopy_Morphology/Tobacco_Project/2019/results/crr_717.tif", format="GTiff", overwrite = T)

# Check if crr within plot is normally distributed
# crr_619_df <- as.data.frame(crrRast_619)
# shapiro.test(sample(crr_619_df$layer, 5000))
# histogram(crr_619_df$layer)
# velox_crr_619 <- velox(crrRast_619)
# crr_extract_619 <- velox_crr_619$extract(sp=plots)
# names(crr_extract_619) <- plots$Id

#histogram(crr_extract_619$`502-E`)
#shapiro.test(sample(crr_extract_619$`502-E`, 5000))
#ggqqplot(sample(crr_extract_619$`502-E`, 5000))

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
  # Normalize data
  # norm_crr <- as.data.frame(scale(crrDF[,2:4]))
  # norm_crr$plots <- crrDF$plots
  # norm_crr
}

crr_df_list <- lapply(crrRast_list, crr_metrics)

plots_619 <- merge(plots_619, crr_df_list[[1]], by.x="plotid", by.y="plots")
plots_703 <- merge(plots_703, crr_df_list[[2]], by.x="plotid", by.y="plots")
plots_717 <- merge(plots_717, crr_df_list[[3]], by.x="plotid", by.y="plots")


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

#merge all metrics into spatial polygons
plots_619 <- merge(plots_619, ch_df_list[[1]], by.x="plotid", by.y="plots")
plots_703 <- merge(plots_703, ch_df_list[[2]], by.x="plotid", by.y="plots")
plots_717 <- merge(plots_717, ch_df_list[[3]], by.x="plotid", by.y="plots")

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


csm_list <- list(csm_masked_619, csm_masked_703, csm_masked_717)
rumple_df_list <- lapply(csm_list, rumple)

#merge all metrics into spatial polygons
plots_619 <- merge(plots_619, rumple_df_list[[1]], by.x="plotid", by.y="plots")
plots_703 <- merge(plots_703, rumple_df_list[[2]], by.x="plotid", by.y="plots")
plots_717 <- merge(plots_717, rumple_df_list[[3]], by.x="plotid", by.y="plots")


# ------------------------------------------------ spatial autocorrelation -----------------------------------------------------

autocor_metrics <- function(x, w){
  plotNums <- plots$plotid
  morans_all <- numeric(length(plotNums))
  #moran_local_all <- numeric(length(plotNums))
  geary_all <- numeric(length(plotNums))
  count <- 1
  
  for (num in plotNums){
    tmp <- plots@data
    position <- match(num, tmp$plotid)
    extract <- crop(x, plots[position,])
    moran <- Moran(extract, w)
    morans_all[count] <- moran
    #moran_local <- MoranLocal(extract, w)
    #moran_local_all[count] <- moran_local
    geary <- Geary(extract, w)
    geary_all[count] <- geary
    count <- count+1
  }
  
  names(morans_all) <- plots$plotid
  #names(morans_local_all) <- plots$plotid
  names(geary_all) <- plots$plotid
  autocor_df <- data.frame(plots=plots$plotid, moran= morans_all, geary=geary_all)
  autocor_df
}

#Using approx 65 cm x 65 cm neighborhood matrix (approx in row plant spacing)
autocor_619 <- autocor_metrics(csm_619, matrix(c(rep.int(1,84), 0, rep.int(1,84)), nc=13, nr=13))
autocor_703 <- autocor_metrics(csm_703, matrix(c(rep.int(1,84), 0, rep.int(1,84)), nc=13, nr=13))
autocor_717 <- autocor_metrics(csm_717, matrix(c(rep.int(1,84), 0, rep.int(1,84)), nc=13, nr=13))

plots_619 <- merge(plots_619, autocor_619, by.x="plotid", by.y="plots")
plots_703 <- merge(plots_703, autocor_703, by.x="plotid", by.y="plots")
plots_717 <- merge(plots_717, autocor_717, by.x="plotid", by.y="plots")


# Create rasters representing local Moran's I for visualization
autocor_rasters <- function(x, w){
  plotNums <- plots$plotid
  moran_local_all <- as.list(x)
  count <- 1

  for (num in plotNums){
    tmp <- plots@data
    position <- match(num, tmp$plotid)
    cropped <- crop(x, plots[position,])
    #masked <- mask(cropped, plots[position,])
    moran_local <- GearyLocal(cropped, w)
    #NAvalue(moran_local) <- 0
    moran_local_all <- append(moran_local_all, moran_local)
    count <- count+1
  }

  moran_local_all <- moran_local_all[-1]
  names(moran_local_all) <- plots$plotid
  moran_local_all
}

local_moran_plots_619 <- autocor_rasters(csm_masked_619, matrix(c(rep.int(1,84), 0, rep.int(1,84)), nc=13, nr=13))
names(local_moran_plots_619) <- NULL
local_moran_plots_619$filename <- 'gearylocal_619.tif'
local_moran_plots_619$overwrite <- TRUE
local_moran_plots_619$fun <- sum

merged_local_moran_plots_619 <- do.call(raster::merge, local_moran_plots_619)
#writeRaster(merged_local_moran_plots_619, "Q:/My Drive/Research/Canopy_Morphology/Tobacco_Project/2019/results/local_moran_plots_619.tif", format="GTiff", overwrite = T)

local_moran_plots_703 <- autocor_rasters(csm_703, matrix(c(1,1,1,1,1,1,1,1,1,1,1,1,0,1,1,1,1,1,1,1,1,1,1,1,1), nc=5, nr=5))
names(local_moran_plots_703) <- NULL
local_moran_plots_703$filename <- 'moran_703.tif'
local_moran_plots_703$overwrite <- TRUE
local_moran_plots_703$fun <- sum
merged_local_moran_plots_703 <- do.call(raster::merge, local_moran_plots_703)
#writeRaster(merged_local_moran_plots_703, "Q:/My Drive/Research/Canopy_Morphology/Tobacco_Project/canopy_analysis/layers/moran/local_moran_plots_703.tif", format="GTiff", overwrite = T)

local_moran_plots_717 <- autocor_rasters(csm_717, matrix(c(1,1,1,1,1,1,1,1,1,1,1,1,0,1,1,1,1,1,1,1,1,1,1,1,1), nc=5, nr=5))
names(local_moran_plots_717) <- NULL
local_moran_plots_717$filename <- 'moran_717.tif'
local_moran_plots_717$overwrite <- TRUE
local_moran_plots_717$fun <- sum
merged_local_moran_plots_717 <- do.call(raster::merge, local_moran_plots_717)
#writeRaster(merged_local_moran_plots_717, "Q:/My Drive/Research/Canopy_Morphology/Tobacco_Project/canopy_analysis/layers/moran/local_moran_plots_717.tif", format="GTiff", overwrite = T)

# ------------------------------------- Aggregate plots to match nutrient observations --------------------------------------

agg_619 <- aggregate(plots_619@data[,23:67], by = list(plots_619$join_id), FUN = mean)
agg_619 <- agg_619[,c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,40,41,42,44,45,46)]
names(agg_619)[15] <- "spad"

agg_703 <- aggregate(plots_703@data[,23:67], by = list(plots_703$join_id), FUN = mean)
agg_703 <- agg_703[,c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,40,41,42,44,45,46)]
names(agg_703)[15] <- "spad"

agg_717 <- aggregate(plots_717@data[,23:67], by = list(plots_717$join_id), FUN = mean)
agg_717 <- agg_717[,c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,40,41,42,44,45,46)]
names(agg_717)[15] <- "spad"

# ------------------------------------Plot variables ------------------------------------------
library(ggpubr)

pairs(agg_619[,27:39], 
      main="Simple Scatterplot Matrix")

cor(agg_619[,2:40])

#plot(plots_619$mean, plots_619$BuAc)
plot(plots_619$N, plots_619$K, xlab = "Median Height (m)", ylab = "Yield (bu/ac)", cex = 1,cex.lab = 2, cex.axis = 2, pch=16)
lm_median <- lm(BuAc ~ median, data = plots_1cm)
text(par()$usr[1], par()$usr[4], bquote(R^2 ~ "=" ~ .(round(summary(lm_median)$r.squared, 2))), 
     cex = 2, col = 2, adj = c(-2.1, 1.6))
shapiro.test(plots_1cm$BuAc)

cor.test(plots_619$N, plots_619$K)
ggscatter(plots_1cm@data, x = "median", y = "BuAc", 
          add = "reg.line", conf.int = TRUE, 
          cor.coef = TRUE, cor.method = "pearson",
          xlab = "Median Height (m)", ylab = "Yield (bu/ac)")

#plot(plots_619$sum, plots_619$BuAc)
plot(plots_619$iqr, plots_619$BuAc, xlab = "Interquartile Range (m)", ylab = "Yield (bu/ac)", cex = 1,cex.lab = 2, cex.axis = 2, pch=16)
lm_iqr <- lm(BuAc ~ iqr, data = plots_1cm)
text(par()$usr[1], par()$usr[4], bquote(R^2 ~ "=" ~ .(round(summary(lm_iqr)$r.squared, 2))), 
     cex = 2, col = 2, adj = c(-2.1, 1.6))

plot(plots_619$rumple, plots_619$BuAc, xlab = "Rumple Index", ylab = "Yield (bu/ac)", cex = 1,cex.lab = 2, cex.axis = 2, pch=16)
lm_rumple <- lm(BuAc ~ rumple, data = plots_1cm)
text(par()$usr[1], par()$usr[4], bquote(R^2 ~ "=" ~ .(round(summary(lm_rumple)$r.squared, 2))), 
     cex = 2, col = 2, adj = c(-2.1, 1.6))

plot(plots_619$skew, plots_619$BuAc, xlab = "Skewness", ylab = "Yield (bu/ac)", cex = 1,cex.lab = 2, cex.axis = 2, pch=16)
lm_skew <- lm(BuAc ~ skew, data = plots_1cm)
text(par()$usr[1], par()$usr[4], bquote(R^2 ~ "=" ~ .(round(summary(lm_skew)$r.squared, 2))), 
     cex = 2, col = 2, adj = c(-2.1, 1.6))

plot(plots_619$kurt, plots_619$BuAc, xlab = "Kurtosis", ylab = "Yield (bu/ac)", cex = 1,cex.lab = 2, cex.axis = 2, pch=16)
lm_kurt <- lm(BuAc ~ kurt, data = plots_1cm)
text(par()$usr[1], par()$usr[4], bquote(R^2 ~ "=" ~ .(round(summary(lm_kurt)$r.squared, 2))), 
     cex = 2, col = 2, adj = c(-2.1, 1.6))

#plot(plots_619$sd, plots_619$BuAc)
plot(plots_619$crrmean, plots_619$BuAc, xlab = "CRR Mean", ylab = "Yield (bu/ac)", cex = 1,cex.lab = 2, cex.axis = 2, pch=16)
lm_crrmean <- lm(BuAc ~ crrmean, data = plots_1cm)
text(par()$usr[1], par()$usr[4], bquote(R^2 ~ "=" ~ .(round(summary(lm_crrmean)$r.squared, 2))), 
     cex = 2, col = 2, adj = c(-2.1, 1.6))

plot(plots_619$crrsd, plots_619$BuAc, xlab = "CRR Standard Deviation", ylab = "Yield (bu/ac)", cex = 1,cex.lab = 2, cex.axis = 2, pch=16)
lm_crrsd <- lm(BuAc ~ crrsd, data = plots_1cm)
text(par()$usr[1], par()$usr[4], bquote(R^2 ~ "=" ~ .(round(summary(lm_crrsd)$r.squared, 2))), 
     cex = 2, col = 2, adj = c(-3.4, 1.6))

#plot(plots_619$PlotCRR, plots_619$BuAc)
# plot(plots_619$moran, plots_619$BuAc)
# lm_moran <- lm(BuAc ~ moran, data = plots_1cm)
# text(par()$usr[1], par()$usr[4], bquote(R^2 ~ "=" ~ .(round(summary(lm_moran)$r.squared, 2))), 
#      cex = 2, col = 2, adj = c(-2.1, 1.6))

plot(plots_1cm$geary, plots_619$BuAc, xlab = "Geary's C", ylab = "Yield (bu/ac)", cex = 1,cex.lab = 2, cex.axis = 2, pch=16)
lm_geary <- lm(BuAc ~ geary, data = plots_1cm)
text(par()$usr[1], par()$usr[4], bquote(R^2 ~ "=" ~ .(round(summary(lm_geary)$r.squared, 2))), 
     cex = 2, col = 2, adj = c(-3.4, 1.6))

#tmp <- lm(BuAc~median, data=plots_619)
#summary(tmp)
#plot(plots_619$sill, plots_619$BuAc)
#plot(plots_619$range, plots_619$BuAc)
#plot(plots_619$Trt, plots_619$BuAc)
plot(plots_619$AUDPC, plots_619$BuAc, xlab = "AUDPC", ylab = "Yield (bu/ac)")
lm_audpc <- lm(BuAc ~ AUDPC, data = plots_1cm)
text(par()$usr[1], par()$usr[4], bquote(R^2 ~ "=" ~ .(round(summary(lm_audpc)$r.squared, 2))), 
     cex = 1.5, col = 2, adj = c(-6.7, 1.4))
#plot(plots_619$watersum, plots_619$BuAc)
#plot(plots_619$watersd, plots_619$BuAc)

par(mfrow=c(1,1))

# Plot histogram distribution of each variable
histogram(plots_619$N, nint=52)
histogram(plots_619$P, nint=52)
histogram(plots_619$K, nint=52)
histogram(plots_619$Mg, nint=52)
histogram(plots_619$Ca, nint=52)
histogram(plots_619$S, nint=52)
histogram(plots_619$Zn, nint=52)
histogram(plots_619$Mn, nint=52)
histogram(plots_619$Cu, nint=30)
histogram(plots_619$Fe, nint=30)
histogram(plots_619$B, nint=30)
histogram(plots_619$Al, nint=30)

# -------------------------------------------------- Normal transformation-------------------------------------------

#Check for normal distribution of independent variables and transform 

shapiro.test(agg_619$mean)
shapiro.test(agg_619$median)
shapiro.test(agg_619$sd)
shapiro.test(agg_619$iqr)
shapiro.test(agg_619$crrmedian)
shapiro.test(agg_619$crrmean)
shapiro.test(agg_619$crrsd) # not normal
agg_619$crrsd_ln <- transformTukey(agg_619$crrsd)
shapiro.test(agg_619$crriqr) # not normal
agg_619$crriqr_ln <- transformTukey(agg_619$crriqr)
shapiro.test(agg_619$skew)
shapiro.test(agg_619$kurt)
shapiro.test(agg_619$moran)
shapiro.test(agg_619$geary)
shapiro.test(agg_619$PlotCRR)
shapiro.test(agg_619$watersum) # not normal
agg_619$watersum_ln <- transformTukey(agg_619$watersum)
shapiro.test(agg_619$varInd_median)
shapiro.test(agg_619$varInd_mean)
shapiro.test(agg_619$varInd_sd) # not normal
agg_619$varInd_sd_ln <- transformTukey(agg_619$varInd_sd)
shapiro.test(agg_619$varInd_iqr)
shapiro.test(agg_619$tgi_mean) # not normal
agg_619$tgi_mean_ln <- transformTukey(agg_619$tgi_mean)
shapiro.test(agg_619$tgi_median) # not normal
agg_619$tgi_median_ln <- transformTukey(agg_619$tgi_median)
shapiro.test(agg_619$tgi_sd)
shapiro.test(agg_619$tgi_iqr)

shapiro.test(agg_703$mean)
shapiro.test(agg_703$median)
shapiro.test(agg_703$sd)
shapiro.test(agg_703$iqr) # not normal
agg_703$iqr_ln <- transformTukey(agg_703$iqr)
shapiro.test(agg_703$crrmedian)
shapiro.test(agg_703$crrmean)
shapiro.test(agg_703$crrsd)
shapiro.test(agg_703$crriqr)
shapiro.test(agg_703$skew)
shapiro.test(agg_703$kurt) # not normal
agg_703$kurt_ln <- transformTukey(agg_703$kurt)
shapiro.test(agg_703$moran)
shapiro.test(agg_703$geary)
shapiro.test(agg_703$PlotCRR)
shapiro.test(agg_703$watersum) # not normal
agg_703$watersum_ln <- transformTukey(agg_703$watersum)
shapiro.test(agg_703$varInd_median) # not normal
agg_703$varInd_median_ln <- transformTukey(agg_703$varInd_median)
shapiro.test(agg_703$varInd_mean) # not normal
agg_703$varInd_mean_ln <- transformTukey(agg_703$varInd_mean)
shapiro.test(agg_703$varInd_sd)
shapiro.test(agg_703$varInd_iqr)
shapiro.test(agg_703$tgi_mean) # not normal
agg_703$tgi_mean_ln <- transformTukey(agg_703$tgi_mean)
shapiro.test(agg_703$tgi_median) # not normal
agg_703$tgi_median_ln <- transformTukey(agg_703$tgi_median)
shapiro.test(agg_703$tgi_sd) # not normal
agg_703$tgi_sd_ln <- transformTukey(agg_703$tgi_sd)
shapiro.test(agg_703$tgi_iqr) # not normal
agg_703$tgi_iqr_ln <- transformTukey(agg_703$tgi_iqr)

shapiro.test(agg_717$mean)
shapiro.test(agg_717$median)
shapiro.test(agg_717$sd)
shapiro.test(agg_717$iqr)
shapiro.test(agg_717$crrmedian)
shapiro.test(agg_717$crrmean)
shapiro.test(agg_717$crrsd)
shapiro.test(agg_717$crriqr)
shapiro.test(agg_717$skew)
shapiro.test(agg_717$kurt) # not normal, cannot transform
shapiro.test(agg_717$moran)
shapiro.test(agg_717$geary) # not normal
agg_717$geary_ln <- transformTukey(agg_717$geary)
shapiro.test(agg_717$PlotCRR)
shapiro.test(agg_717$watersum) # not normal
agg_717$watersum_ln <- transformTukey(agg_717$watersum)
shapiro.test(agg_717$varInd_median) # not normal
agg_717$varInd_median_ln <- transformTukey(agg_717$varInd_median)
shapiro.test(agg_717$varInd_mean)
shapiro.test(agg_717$varInd_sd)
shapiro.test(agg_717$varInd_iqr)
shapiro.test(agg_717$tgi_mean)
shapiro.test(agg_717$tgi_median)
shapiro.test(agg_717$tgi_sd) # not normal
agg_717$tgi_sd_ln <- transformTukey(agg_717$tgi_sd)
shapiro.test(agg_717$tgi_iqr) # not normal
agg_717$tgi_iqr_ln <- transformTukey(agg_717$tgi_iqr)

# Check for normal distribution of dependent variables
shapiro.test(agg_619$N)



# --------------------------------------- Stats data prep --------------------------------------------------------------------------------
agg_619$treatment <- tstrsplit(agg_619$Group.1, "_")[[2]]
agg_619$fert <- tstrsplit(agg_619$Group.1, "_")[[1]]
agg_703$treatment <- tstrsplit(agg_703$Group.1, "_")[[2]]
agg_703$fert <- tstrsplit(agg_703$Group.1, "_")[[1]]
agg_717$treatment <- tstrsplit(agg_717$Group.1, "_")[[2]]
agg_717$fert <- tstrsplit(agg_717$Group.1, "_")[[1]]


normalize <- function(x)
{
  return((x- min(x)) /(max(x)-min(x)))
}


#minmax normalization of data
norm_619 <- as.data.frame(normalize(agg_619[,2:46]))
norm_619$Id <- agg_619$Group.1
norm_619$treatment <- agg_619$treatment
norm_619$fert <- agg_619$fert

norm_703 <- as.data.frame(normalize(agg_703[,2:49]))
norm_703$Id <- agg_703$Group.1
norm_703$treatment <- agg_703$treatment
norm_703$fert <- agg_703$fert

norm_717 <- as.data.frame(normalize(agg_717[,2:45]))
norm_717$Id <- agg_717$Group.1
norm_717$treatment <- agg_717$treatment
norm_717$fert <- agg_717$fert


# ---------------------------------------- MANOVA ------------------------------------------------------
med_619 <- norm_619$median
moran_619 <- norm_619$moran
crrmed_619 <- norm_619$crrmedian
crriqr_619 <- norm_619$crriqr
skew_619 <- norm_619$skew
kurt_619 <- norm_619$kurt
rumple_619 <- norm_619$rumple

norm_619$Treatment <- as.factor(norm_619$Treatment)
man_619 <- manova(cbind(med_619, moran_619, crrmed_619, crriqr_619, skew_619, kurt_619, rumple_619)~Treatment, data=norm_619)
summary(man_619)
summary.aov(man_619)


med_703 <- norm_703$median_ln
moran_703 <- norm_703$moran_ln
crrmed_703 <- norm_703$crrmedian_ln
skew_703 <- norm_703$skew_ln
kurt_703 <- norm_703$kurt_ln
rumple_703 <- norm_703$rumple

man_703 <- manova(cbind(med_703, moran_703, crrmed_703, skew_703, kurt_703, rumple_703)~Treatment, data=norm_703)
summary(man_703)
summary.aov(man_703)


med_717 <- norm_717$median_ln
moran_717 <- norm_717$moran_ln
crrmed_717 <- norm_717$crrmedian_ln
skew_717 <- norm_717$skew_ln
kurt_717 <- norm_717$kurt_ln
rumple_717 <- norm_717$rumple

man_717 <- manova(cbind(med_717, moran_717, crrmed_717, skew_717, kurt_717, rumple_717)~Trt_agg, data=norm_717)
summary(man_717)
summary.aov(man_717)

#---------------------------------------- LM structural only (partial mask) ------------------------------------------

lm_N_619 <- lm(N ~ crrmedian+crriqr_ln+kurt+rumple+moran+geary+watersum_ln, data = norm_619, na.action=na.fail)
vif(lm_N_619)
modsel_N_619 <- dredge(lm_N_619)
topmod_N_619 <- lm(N ~ crrmedian+moran+kurt, data = norm_619,na.action=na.fail)
vif(topmod_N_619)
summary(topmod_N_619)
N_619_pred <- predict(topmod_N_619)
N_619_res <- resid(topmod_N_619)
n_619_table <- cbind(norm_619$N, N_619_pred, N_619_res)
head(n_619_table)

lm_P_619 <- lm(P ~ crrmedian+crriqr_ln+kurt+rumple+moran+geary+watersum_ln, data = norm_619, na.action=na.fail)
summary(lm_P_619)
vif(lm_P_619)
modsel_P_619 <- dredge(lm_P_619)
topmod_P_619 <- lm(P ~ crrmedian+geary+watersum_ln, data = norm_619,na.action=na.fail)
summary(topmod_P_619)

lm_K_619 <- lm(K ~ crrmedian+crriqr_ln+kurt+rumple+moran+geary+watersum_ln, data = norm_619, na.action=na.fail)
summary(lm_K_619)
modsel_K_619 <- dredge(lm_K_619)
topmod_K_619 <- lm(K ~ crrmedian+watersum_ln, data = norm_619,na.action=na.fail)
summary(topmod_K_619)

lm_Ca_619 <- lm(Ca ~crrmedian+crriqr_ln+kurt+rumple+moran+geary+watersum_ln, data = norm_619, na.action=na.fail)
summary(lm_Ca_619)
modsel_Ca_619 <- dredge(lm_Ca_619)
topmod_Ca_619 <- lm(Ca ~ watersum_ln+moran, data = norm_619,na.action=na.fail)
summary(topmod_Ca_619)

lm_S_619 <- lm(S ~ crrmedian+crriqr_ln+kurt+rumple+moran+geary+watersum_ln, data = norm_619, na.action=na.fail)
summary(lm_S_619)
modsel_S_619 <- dredge(lm_S_619)
topmod_S_619 <- lm(S ~ crrmedian, data = norm_619,na.action=na.fail)
summary(topmod_S_619)

lm_B_619 <- lm(B ~ crrmedian+crriqr_ln+kurt+rumple+moran+geary+watersum_ln, data = norm_619, na.action=na.fail)
summary(lm_B_619)
modsel_B_619 <- dredge(lm_B_619)
topmod_B_619 <- lm(B ~ geary+moran, data = norm_619,na.action=na.fail)
vif(topmod_B_619)
summary(topmod_B_619)

lm_Zn_619 <- lm(Zn ~ crrmedian+crriqr_ln+kurt+rumple+moran+geary+watersum_ln, data = norm_619, na.action=na.fail)
summary(lm_Zn_619)
modsel_Zn_619 <- dredge(lm_Zn_619)
topmod_Zn_619 <- lm(Zn ~ crrmedian, data = norm_619,na.action=na.fail)
summary(topmod_Zn_619)

lm_Mg_619 <- lm(Mg ~  crrmedian+crriqr_ln+kurt+rumple+moran+geary+watersum_ln, data = norm_619, na.action=na.fail)
summary(lm_Mg_619)
modsel_Mg_619 <- dredge(lm_Mg_619)
topmod_Mg_619 <- lm(Mg ~ crrmedian+kurt+moran, data = norm_619,na.action=na.fail)
summary(topmod_Mg_619)

# function to get lm pvalue
lmp <- function (modelobject) {
  if (class(modelobject) != "lm") stop("Not an object of class 'lm' ")
  f <- summary(modelobject)$fstatistic
  p <- pf(f[1],f[2],f[3],lower.tail=F)
  attributes(p) <- NULL
  return(p)
}

res_N_619 <- as.data.frame(modsel_N_619[1])
res_N_619$adj_r2 <- summary(topmod_N_619)$adj.r.squared
res_N_619$pvalue <- lmp(topmod_N_619)

res_P_619 <- as.data.frame(modsel_P_619[1])
res_P_619$adj_r2 <- summary(topmod_P_619)$adj.r.squared
res_P_619$pvalue <- lmp(topmod_P_619)

res_K_619 <- as.data.frame(modsel_K_619[1])
res_K_619$adj_r2 <- summary(topmod_K_619)$adj.r.squared
res_K_619$pvalue <- lmp(topmod_K_619)

res_B_619 <- as.data.frame(modsel_B_619[1])
res_B_619$adj_r2 <- summary(topmod_B_619)$adj.r.squared
res_B_619$pvalue <- lmp(topmod_B_619)

res_Ca_619 <- as.data.frame(modsel_Ca_619[1])
res_Ca_619$adj_r2 <- summary(topmod_Ca_619)$adj.r.squared
res_Ca_619$pvalue <- lmp(topmod_Ca_619)

res_Mg_619 <- as.data.frame(modsel_Mg_619[1])
res_Mg_619$adj_r2 <- summary(topmod_Mg_619)$adj.r.squared
res_Mg_619$pvalue <- lmp(topmod_Mg_619)

res_S_619 <- as.data.frame(modsel_S_619[1])
res_S_619$adj_r2 <- summary(topmod_S_619)$adj.r.squared
res_S_619$pvalue <- lmp(topmod_S_619)

res_Zn_619 <- as.data.frame(modsel_Zn_619[1])
res_Zn_619$adj_r2 <- summary(topmod_Zn_619)$adj.r.squared
res_Zn_619$pvalue <- lmp(topmod_Zn_619)

res_619 <- rbind(res_N_619, res_P_619, res_K_619, res_B_619, res_Ca_619, res_Mg_619, res_S_619, res_Zn_619)
res_619$dependent <- c("N", "P", "K", "B", "Ca", "Mg", "S", "Zn")
write.csv(res_619, "results/res_struc_partialmask_619.csv")

# --------------------------------------------- All variables ---------------------------------------------------------

mean+median+sd+iqr+crrmean+crrmedian+crrsd+crriqr+skew+kurt+geary+moran+watersum+tgi_mean+tgi_median+tgi_sd+tgi_iqr+varInd_mean+varInd_median+varInd_sd+varInd_iqr

#------------------------------------------- LM with spectral+structural --------------------------------------------------------
lm_N_619 <- lm(N ~ crriqr_ln+crrmedian+kurt+rumple+moran+geary+watersum_ln+tgi_mean_ln+tgi_iqr, data = norm_619, na.action=na.fail)
vif(lm_N_619)
modsel_N_619 <- dredge(lm_N_619)
topmod_N_619 <- lm(N ~ geary+kurt+rumple+tgi_median_ln+watersum_ln, data = norm_619,na.action=na.fail)
summary(topmod_N_619)
N_619_pred <- predict(topmod_N_619)
N_619_res <- resid(topmod_N_619)
n_619_table <- cbind(norm_619$N, N_619_pred, N_619_res)
head(n_619_table)

lm_P_619 <- lm(P ~ crriqr_ln+crrmedian+kurt+rumple+moran+geary+watersum_ln+tgi_median_ln+tgi_iqr, data = norm_619, na.action=na.fail)
summary(lm_P_619)
modsel_P_619 <- dredge(lm_P_619)
topmod_P_619 <- lm(P ~ geary+tgi_iqr+watersum_ln, data = norm_619,na.action=na.fail)
summary(topmod_P_619)

lm_K_619 <- lm(K ~ crriqr_ln+crrmedian+kurt+rumple+moran+geary+watersum_ln+tgi_median_ln+tgi_iqr, data = norm_619, na.action=na.fail)
summary(lm_K_619)
modsel_K_619 <- dredge(lm_K_619)
topmod_K_619 <- lm(K ~ tgi_median_ln+watersum_ln, data = norm_619,na.action=na.fail)
summary(topmod_K_619)

lm_Ca_619 <- lm(Ca ~ crriqr_ln+crrmedian+kurt+rumple+moran+geary+watersum_ln+tgi_median_ln+tgi_iqr, data = norm_619, na.action=na.fail)
summary(lm_Ca_619)
modsel_Ca_619 <- dredge(lm_Ca_619)
topmod_Ca_619 <- lm(Ca ~ moran+watersum_ln, data = norm_619,na.action=na.fail)
summary(topmod_Ca_619)

lm_S_619 <- lm(S ~ crriqr_ln+skew+kurt+rumple+moran+geary+watersum_ln+tgi_median_ln+tgi_iqr, data = norm_619, na.action=na.fail)
summary(lm_S_619)
modsel_S_619 <- dredge(lm_S_619)
topmod_S_619 <- lm(S ~ geary+moran+tgi_iqr, data = norm_619,na.action=na.fail)
summary(topmod_S_619)

lm_B_619 <- lm(B ~ crriqr_ln+skew+kurt+rumple+moran+geary+watersum_ln+tgi_median_ln+tgi_iqr, data = norm_619, na.action=na.fail)
summary(lm_B_619)
modsel_B_619 <- dredge(lm_B_619)
topmod_B_619 <- lm(B ~ geary+moran, data = norm_619,na.action=na.fail)
summary(topmod_B_619)

lm_Zn_619 <- lm(Zn ~ crriqr_ln+skew+kurt+rumple+moran+geary+watersum_ln+tgi_median_ln+tgi_iqr, data = norm_619, na.action=na.fail)
summary(lm_Zn_619)
modsel_Zn_619 <- dredge(lm_Zn_619)
topmod_Zn_619 <- lm(Zn ~ rumple+tgi_iqr, data = norm_619,na.action=na.fail)
summary(topmod_Zn_619)

lm_Mg_619 <- lm(Mg ~ crriqr_ln+skew+kurt+rumple+moran+geary+watersum_ln+tgi_median_ln+tgi_iqr, data = norm_619, na.action=na.fail)
summary(lm_Mg_619)
modsel_Mg_619 <- dredge(lm_Mg_619)
topmod_Mg_619 <- lm(Mg ~ crriqr_ln+moran+tgi_median_ln+watersum_ln, data = norm_619,na.action=na.fail)
summary(topmod_Mg_619)

res_N_619 <- as.data.frame(modsel_N_619[1])
res_N_619$adj_r2 <- summary(topmod_N_619)$adj.r.squared
res_N_619$pvalue <- lmp(topmod_N_619)

res_P_619 <- as.data.frame(modsel_P_619[1])
res_P_619$adj_r2 <- summary(topmod_P_619)$adj.r.squared
res_P_619$pvalue <- lmp(topmod_P_619)

res_K_619 <- as.data.frame(modsel_K_619[1])
res_K_619$adj_r2 <- summary(topmod_K_619)$adj.r.squared
res_K_619$pvalue <- lmp(topmod_K_619)

res_B_619 <- as.data.frame(modsel_B_619[1])
res_B_619$adj_r2 <- summary(topmod_B_619)$adj.r.squared
res_B_619$pvalue <- lmp(topmod_B_619)

res_Ca_619 <- as.data.frame(modsel_Ca_619[1])
res_Ca_619$adj_r2 <- summary(topmod_Ca_619)$adj.r.squared
res_Ca_619$pvalue <- lmp(topmod_Ca_619)

res_Mg_619 <- as.data.frame(modsel_Mg_619[1])
res_Mg_619$adj_r2 <- summary(topmod_Mg_619)$adj.r.squared
res_Mg_619$pvalue <- lmp(topmod_Mg_619)

res_S_619 <- as.data.frame(modsel_S_619[1])
res_S_619$adj_r2 <- summary(topmod_S_619)$adj.r.squared
res_S_619$pvalue <- lmp(topmod_S_619)

res_Zn_619 <- as.data.frame(modsel_Zn_619[1])
res_Zn_619$adj_r2 <- summary(topmod_Zn_619)$adj.r.squared
res_Zn_619$pvalue <- lmp(topmod_Zn_619)

res_619_spectal <- rbind(res_N_619, res_P_619, res_K_619, res_B_619, res_Ca_619, res_Mg_619, res_S_619, res_Zn_619)
res_619_spectal$dependent <- c("N", "P", "K", "B", "Ca", "Mg", "S", "Zn")
write.csv(res_619_spectal, "results/res_619_spectral_struc.csv")

#------------------------------------------- LM with spectral only --------------------------------------------------------
lm_N_619 <- lm(N ~ tgi_mean_ln+tgi_sd+tgi_iqr+watersum_ln, data = norm_619, na.action=na.fail)
summary(lm_N_619)
vif(lm_N_619)
modsel_N_619 <- dredge(lm_N_619)
topmod_N_619 <- lm(N ~ tgi_mean_ln+tgi_iqr+watersum_ln, data = norm_619,na.action=na.fail)
summary(topmod_N_619)
N_619_pred <- predict(topmod_N_619)
N_619_res <- resid(topmod_N_619)
n_619_table <- cbind(norm_619$N, N_619_pred, N_619_res)
head(n_619_table)

lm_P_619 <- lm(P ~ tgi_mean_ln+tgi_sd+tgi_iqr+watersum_ln, data = norm_619, na.action=na.fail)
summary(lm_P_619)
modsel_P_619 <- dredge(lm_P_619)
topmod_P_619 <- lm(P ~ tgi_iqr+watersum_ln, data = norm_619,na.action=na.fail)
summary(topmod_P_619)

lm_K_619 <- lm(K ~ tgi_mean_ln+tgi_sd+tgi_iqr+watersum_ln, data = norm_619, na.action=na.fail)
summary(lm_K_619)
modsel_K_619 <- dredge(lm_K_619)
topmod_K_619 <- lm(K ~ tgi_mean_ln+tgi_sd+watersum_ln, data = norm_619,na.action=na.fail)
summary(topmod_K_619)

lm_Ca_619 <- lm(Ca ~ tgi_mean_ln+tgi_sd+tgi_iqr+watersum_ln, data = norm_619, na.action=na.fail)
summary(lm_Ca_619)
modsel_Ca_619 <- dredge(lm_Ca_619)
topmod_Ca_619 <- lm(Ca ~ watersum_ln, data = norm_619,na.action=na.fail)
summary(topmod_Ca_619)

lm_S_619 <- lm(S ~ tgi_mean_ln+tgi_sd+tgi_iqr+watersum_ln, data = norm_619, na.action=na.fail)
summary(lm_S_619)
modsel_S_619 <- dredge(lm_S_619)
topmod_S_619 <- lm(S ~ tgi_mean_ln, data = norm_619,na.action=na.fail)
summary(topmod_S_619)

lm_B_619 <- lm(B ~ tgi_mean_ln+tgi_sd+tgi_iqr+watersum_ln, data = norm_619, na.action=na.fail)
summary(lm_B_619)
modsel_B_619 <- dredge(lm_B_619)
topmod_B_619 <- lm(B ~ 1, data = norm_619,na.action=na.fail)
summary(topmod_B_619)

lm_Zn_619 <- lm(Zn ~ tgi_mean_ln+tgi_sd+tgi_iqr+watersum_ln, data = norm_619, na.action=na.fail)
summary(lm_Zn_619)
modsel_Zn_619 <- dredge(lm_Zn_619)
topmod_Zn_619 <- lm(Zn ~ tgi_iqr, data = norm_619,na.action=na.fail)
summary(topmod_Zn_619)

lm_Mg_619 <- lm(Mg ~ tgi_mean_ln+tgi_sd+tgi_iqr+watersum_ln, data = norm_619, na.action=na.fail)
summary(lm_Mg_619)
modsel_Mg_619 <- dredge(lm_Mg_619)
topmod_Mg_619 <- lm(Mg ~ tgi_mean_ln+watersum_ln, data = norm_619,na.action=na.fail)
summary(topmod_Mg_619)

res_N_619 <- as.data.frame(modsel_N_619[1])
res_N_619$adj_r2 <- summary(topmod_N_619)$adj.r.squared
res_N_619$pvalue <- lmp(topmod_N_619)

res_P_619 <- as.data.frame(modsel_P_619[1])
res_P_619$adj_r2 <- summary(topmod_P_619)$adj.r.squared
res_P_619$pvalue <- lmp(topmod_P_619)

res_K_619 <- as.data.frame(modsel_K_619[1])
res_K_619$adj_r2 <- summary(topmod_K_619)$adj.r.squared
res_K_619$pvalue <- lmp(topmod_K_619)

res_B_619 <- as.data.frame(modsel_B_619[1])
res_B_619$adj_r2 <- summary(topmod_B_619)$adj.r.squared
res_B_619$pvalue <- lmp(topmod_B_619)

res_Ca_619 <- as.data.frame(modsel_Ca_619[1])
res_Ca_619$adj_r2 <- summary(topmod_Ca_619)$adj.r.squared
res_Ca_619$pvalue <- lmp(topmod_Ca_619)

res_Mg_619 <- as.data.frame(modsel_Mg_619[1])
res_Mg_619$adj_r2 <- summary(topmod_Mg_619)$adj.r.squared
res_Mg_619$pvalue <- lmp(topmod_Mg_619)

res_S_619 <- as.data.frame(modsel_S_619[1])
res_S_619$adj_r2 <- summary(topmod_S_619)$adj.r.squared
res_S_619$pvalue <- lmp(topmod_S_619)

res_Zn_619 <- as.data.frame(modsel_Zn_619[1])
res_Zn_619$adj_r2 <- summary(topmod_Zn_619)$adj.r.squared
res_Zn_619$pvalue <- lmp(topmod_Zn_619)

res_619_spectal <- rbind(res_N_619, res_P_619, res_K_619, res_Ca_619, res_Mg_619, res_S_619, res_Zn_619)
res_619_spectal$dependent <- c("N", "P", "K", "Ca", "Mg", "S", "Zn")
write.csv(res_619_spectal, "results/res_619_spectralOnly.csv")

#---------------------------------------------------Repeat for 703 -----------------------------------------------
lm_N_703 <- lm(N ~ median+crrmedian+crrsd+kurt_ln+moran+rumple+watersum_ln, data = norm_703, na.action=na.fail)
vif(lm_N_703)
modsel_N_703 <- dredge(lm_N_703)
topmod_N_703 <- lm(N ~ moran+median, data = norm_703,na.action=na.fail)
summary(topmod_N_703)
N_703_pred <- predict(topmod_N_703)
N_703_res <- resid(topmod_N_703)
n_703_table <- cbind(norm_703$N, N_703_pred, N_703_res)
head(n_703_table)

lm_P_703 <- lm(P ~ median+crrmedian+crrsd+kurt_ln+moran+rumple+watersum_ln, data = norm_703, na.action=na.fail)
summary(lm_P_703)
modsel_P_703 <- dredge(lm_P_703)
topmod_P_703 <- lm(P ~ median+crrmedian, data = norm_703,na.action=na.fail)
summary(topmod_P_703)

lm_K_703 <- lm(K ~ median+crrmedian+crrsd+kurt_ln+moran+rumple+watersum_ln, data = norm_703, na.action=na.fail)
summary(lm_K_703)
modsel_K_703 <- dredge(lm_K_703)
topmod_K_703 <- lm(K ~ crrmedian, data = norm_703,na.action=na.fail)
summary(topmod_K_703)

lm_Ca_703 <- lm(Ca ~ median+crrmedian+crrsd+kurt_ln+moran+rumple+watersum_ln, data = norm_703, na.action=na.fail)
summary(lm_Ca_703)
modsel_Ca_703 <- dredge(lm_Ca_703)
topmod_Ca_703 <- lm(Ca ~ watersum_ln, data = norm_703,na.action=na.fail)
summary(topmod_Ca_703)

lm_S_703 <- lm(S ~ median+crrmedian+crrsd+kurt_ln+moran+rumple+watersum_ln, data = norm_703, na.action=na.fail)
summary(lm_S_703)
modsel_S_703 <- dredge(lm_S_703)
topmod_S_703 <- lm(S ~ crrmedian+rumple, data = norm_703,na.action=na.fail)
summary(topmod_S_703)

lm_B_703 <- lm(B ~ median+crrmedian+crrsd+kurt_ln+moran+rumple+watersum_ln, data = norm_703, na.action=na.fail)
summary(lm_B_703)
modsel_B_703 <- dredge(lm_B_703)
topmod_B_703 <- lm(B ~ crrmedian+moran+rumple+watersum_ln, data = norm_703,na.action=na.fail)
summary(topmod_B_703)

lm_Zn_703 <- lm(Zn ~ median+crrmedian+crrsd+kurt_ln+moran+rumple+watersum_ln, data = norm_703, na.action=na.fail)
summary(lm_Zn_703)
modsel_Zn_703 <- dredge(lm_Zn_703)
topmod_Zn_703 <- lm(Zn ~ 1, data = norm_703,na.action=na.fail)
summary(topmod_Zn_703)

lm_Mg_703 <- lm(Mg ~ median+crrmedian+crrsd+kurt_ln+moran+rumple+watersum_ln, data = norm_703, na.action=na.fail)
summary(lm_Mg_703)
modsel_Mg_703 <- dredge(lm_Mg_703)
topmod_Mg_703 <- lm(Mg ~ 1, data = norm_703,na.action=na.fail)
summary(topmod_Mg_703)

# function to get lm pvalue
lmp <- function (modelobject) {
  if (class(modelobject) != "lm") stop("Not an object of class 'lm' ")
  f <- summary(modelobject)$fstatistic
  p <- pf(f[1],f[2],f[3],lower.tail=F)
  attributes(p) <- NULL
  return(p)
}

res_N_703 <- as.data.frame(modsel_N_703[1])
res_N_703$adj_r2 <- summary(topmod_N_703)$adj.r.squared
res_N_703$pvalue <- lmp(topmod_N_703)

res_P_703 <- as.data.frame(modsel_P_703[1])
res_P_703$adj_r2 <- summary(topmod_P_703)$adj.r.squared
res_P_703$pvalue <- lmp(topmod_P_703)

res_K_703 <- as.data.frame(modsel_K_703[1])
res_K_703$adj_r2 <- summary(topmod_K_703)$adj.r.squared
res_K_703$pvalue <- lmp(topmod_K_703)

res_B_703 <- as.data.frame(modsel_B_703[1])
res_B_703$adj_r2 <- summary(topmod_B_703)$adj.r.squared
res_B_703$pvalue <- lmp(topmod_B_703)

res_Ca_703 <- as.data.frame(modsel_Ca_703[1])
res_Ca_703$adj_r2 <- summary(topmod_Ca_703)$adj.r.squared
res_Ca_703$pvalue <- lmp(topmod_Ca_703)

res_Mg_703 <- as.data.frame(modsel_Mg_703[1])
res_Mg_703$adj_r2 <- summary(topmod_Mg_703)$adj.r.squared
res_Mg_703$pvalue <- lmp(topmod_Mg_703)

res_S_703 <- as.data.frame(modsel_S_703[1])
res_S_703$adj_r2 <- summary(topmod_S_703)$adj.r.squared
res_S_703$pvalue <- lmp(topmod_S_703)

res_Zn_703 <- as.data.frame(modsel_Zn_703[1])
res_Zn_703$adj_r2 <- summary(topmod_Zn_703)$adj.r.squared
res_Zn_703$pvalue <- lmp(topmod_Zn_703)

res_703 <- rbind(res_N_703, res_P_703, res_K_703, res_B_703, res_Ca_703, res_S_703)
res_703$dependent <- c("N", "P", "K", "B", "Ca", "S")
write.csv(res_703, "results/res_struct_703.csv")

#------------------------------------------- LM with spectral and structural --------------------------------------------------------
lm_N_703 <- lm(N ~ median+sd+crrsd+kurt_ln+moran+watersum_ln+tgi_median_ln+tgi_sd_ln, data = norm_703, na.action=na.fail)
vif(lm_N_703)
modsel_N_703 <- dredge(lm_N_703)
topmod_N_703 <- lm(N ~ kurt+sd+tgi_median_ln+watersum_ln, data = norm_703,na.action=na.fail)
summary(topmod_N_703)
N_703_pred <- predict(topmod_N_703)
N_703_res <- resid(topmod_N_703)
n_703_table <- cbind(norm_703$N, N_703_pred, N_703_res)
head(n_703_table)

lm_P_703 <- lm(P ~ median+sd+crrsd+kurt_ln+moran+watersum_ln+tgi_median_ln+tgi_sd_ln, data = norm_703, na.action=na.fail)
summary(lm_P_703)
modsel_P_703 <- dredge(lm_P_703)
topmod_P_703 <- lm(P ~ sd, data = norm_703,na.action=na.fail)
summary(topmod_P_703)

lm_K_703 <- lm(K ~ median+sd+crrsd+kurt_ln+moran+watersum_ln+tgi_median_ln+tgi_sd_ln, data = norm_703, na.action=na.fail)
summary(lm_K_703)
modsel_K_703 <- dredge(lm_K_703)
topmod_K_703 <- lm(K ~ tgi_median_ln, data = norm_703,na.action=na.fail)
summary(topmod_K_703)

lm_Ca_703 <- lm(Ca ~ median+sd+crrsd+kurt_ln+moran+watersum_ln+tgi_median_ln+tgi_sd_ln, data = norm_703, na.action=na.fail)
summary(lm_Ca_703)
modsel_Ca_703 <- dredge(lm_Ca_703)
topmod_Ca_703 <- lm(Ca ~ watersum_ln, data = norm_703,na.action=na.fail)
summary(topmod_Ca_703)

lm_S_703 <- lm(S ~ median+sd+crrsd+kurt_ln+moran+watersum_ln+tgi_median_ln+tgi_sd_ln, data = norm_703, na.action=na.fail)
summary(lm_S_703)
modsel_S_703 <- dredge(lm_S_703)
topmod_S_703 <- lm(S ~ sd+tgi_median_ln, data = norm_703,na.action=na.fail)
summary(topmod_S_703)

lm_B_703 <- lm(B ~ median+sd+crrsd+kurt_ln+moran+watersum_ln+tgi_median_ln+tgi_sd_ln, data = norm_703, na.action=na.fail)
summary(lm_B_703)
modsel_B_703 <- dredge(lm_B_703)
topmod_B_703 <- lm(B ~ sd+tgi_median_ln, data = norm_703,na.action=na.fail)
summary(topmod_B_703)

lm_Zn_703 <- lm(Zn ~ median+sd+crrsd+kurt_ln+moran+watersum_ln+tgi_median_ln+tgi_sd_ln, data = norm_703, na.action=na.fail)
summary(lm_Zn_703)
modsel_Zn_703 <- dredge(lm_Zn_703)
topmod_Zn_703 <- lm(Zn ~ 1, data = norm_703,na.action=na.fail)
summary(topmod_Zn_703)

lm_Mg_703 <- lm(Mg ~ median+sd+crrsd+kurt_ln+moran+watersum_ln+tgi_median_ln+tgi_sd_ln, data = norm_703, na.action=na.fail)
summary(lm_Mg_703)
modsel_Mg_703 <- dredge(lm_Mg_703)
topmod_Mg_703 <- lm(Mg ~ tgi_sd_ln, data = norm_703,na.action=na.fail)
summary(topmod_Mg_703)

res_N_703 <- as.data.frame(modsel_N_703[1])
res_N_703$adj_r2 <- summary(topmod_N_703)$adj.r.squared
res_N_703$pvalue <- lmp(topmod_N_703)

res_P_703 <- as.data.frame(modsel_P_703[1])
res_P_703$adj_r2 <- summary(topmod_P_703)$adj.r.squared
res_P_703$pvalue <- lmp(topmod_P_703)

res_K_703 <- as.data.frame(modsel_K_703[1])
res_K_703$adj_r2 <- summary(topmod_K_703)$adj.r.squared
res_K_703$pvalue <- lmp(topmod_K_703)

res_B_703 <- as.data.frame(modsel_B_703[1])
res_B_703$adj_r2 <- summary(topmod_B_703)$adj.r.squared
res_B_703$pvalue <- lmp(topmod_B_703)

res_Ca_703 <- as.data.frame(modsel_Ca_703[1])
res_Ca_703$adj_r2 <- summary(topmod_Ca_703)$adj.r.squared
res_Ca_703$pvalue <- lmp(topmod_Ca_703)

res_Mg_703 <- as.data.frame(modsel_Mg_703[1])
res_Mg_703$adj_r2 <- summary(topmod_Mg_703)$adj.r.squared
res_Mg_703$pvalue <- lmp(topmod_Mg_703)

res_S_703 <- as.data.frame(modsel_S_703[1])
res_S_703$adj_r2 <- summary(topmod_S_703)$adj.r.squared
res_S_703$pvalue <- lmp(topmod_S_703)

res_Zn_703 <- as.data.frame(modsel_Zn_703[1])
res_Zn_703$adj_r2 <- summary(topmod_Zn_703)$adj.r.squared
res_Zn_703$pvalue <- lmp(topmod_Zn_703)

res_703_spec_struc <- rbind(res_N_703, res_P_703, res_K_703, res_B_703, res_Ca_703, res_Mg_703, res_S_703)
res_703_spec_struc$dependent <- c("N", "P", "K", "B", "Ca", "Mg", "S")
write.csv(res_703_spec_struc, "results/res_703_spec_struc.csv")


#------------------------------------------- LM with spectral only --------------------------------------------------------
lm_N_703 <- lm(N ~tgi_mean_ln+tgi_sd_ln+watersum_ln, data = norm_703, na.action=na.fail)
vif(lm_N_703)
modsel_N_703 <- dredge(lm_N_703)
topmod_N_703 <- lm(N ~ tgi_mean_ln, data = norm_703,na.action=na.fail)
summary(topmod_N_703)
N_703_pred <- predict(topmod_N_703)
N_703_res <- resid(topmod_N_703)
n_703_table <- cbind(norm_703$N, N_703_pred, N_703_res)
head(n_703_table)

lm_P_703 <- lm(P ~tgi_mean_ln+tgi_sd_ln+watersum_ln, data = norm_703, na.action=na.fail)
summary(lm_P_703)
modsel_P_703 <- dredge(lm_P_703)
topmod_P_703 <- lm(P ~ 1, data = norm_703,na.action=na.fail)
summary(topmod_P_703)

lm_K_703 <- lm(K ~tgi_mean_ln+tgi_sd_ln+watersum_ln, data = norm_703, na.action=na.fail)
summary(lm_K_703)
modsel_K_703 <- dredge(lm_K_703)
topmod_K_703 <- lm(K ~ tgi_mean_ln, data = norm_703,na.action=na.fail)
summary(topmod_K_703)

lm_Ca_703 <- lm(Ca ~tgi_mean_ln+tgi_sd_ln+watersum_ln, data = norm_703, na.action=na.fail)
summary(lm_Ca_703)
modsel_Ca_703 <- dredge(lm_Ca_703)
topmod_Ca_703 <- lm(Ca ~ watersum_ln, data = norm_703,na.action=na.fail)
summary(topmod_Ca_703)

lm_S_703 <- lm(S ~tgi_mean_ln+tgi_sd_ln+watersum_ln, data = norm_703, na.action=na.fail)
summary(lm_S_703)
modsel_S_703 <- dredge(lm_S_703)
topmod_S_703 <- lm(S ~ tgi_mean_ln, data = norm_703,na.action=na.fail)
summary(topmod_S_703)

lm_B_703 <- lm(B ~tgi_mean_ln+tgi_sd_ln+watersum_ln, data = norm_703, na.action=na.fail)
summary(lm_B_703)
modsel_B_703 <- dredge(lm_B_703)
topmod_B_703 <- lm(B ~ 1, data = norm_703,na.action=na.fail)
summary(topmod_B_703)

lm_Zn_703 <- lm(Zn ~tgi_mean_ln+tgi_sd_ln+watersum_ln, data = norm_703, na.action=na.fail)
summary(lm_Zn_703)
modsel_Zn_703 <- dredge(lm_Zn_703)
topmod_Zn_703 <- lm(Zn ~ 1, data = norm_703,na.action=na.fail)
summary(topmod_Zn_703)

lm_Mg_703 <- lm(Mg ~tgi_mean_ln+tgi_sd_ln+watersum_ln, data = norm_703, na.action=na.fail)
summary(lm_Mg_703)
modsel_Mg_703 <- dredge(lm_Mg_703)
topmod_Mg_703 <- lm(Mg ~ tgi_sd_ln, data = norm_703,na.action=na.fail)
summary(topmod_Mg_703)

res_N_703 <- as.data.frame(modsel_N_703[1])
res_N_703$adj_r2 <- summary(topmod_N_703)$adj.r.squared
res_N_703$pvalue <- lmp(topmod_N_703)

res_P_703 <- as.data.frame(modsel_P_703[1])
res_P_703$adj_r2 <- summary(topmod_P_703)$adj.r.squared
res_P_703$pvalue <- lmp(topmod_P_703)

res_K_703 <- as.data.frame(modsel_K_703[1])
res_K_703$adj_r2 <- summary(topmod_K_703)$adj.r.squared
res_K_703$pvalue <- lmp(topmod_K_703)

res_B_703 <- as.data.frame(modsel_B_703[1])
res_B_703$adj_r2 <- summary(topmod_B_703)$adj.r.squared
res_B_703$pvalue <- lmp(topmod_B_703)

res_Ca_703 <- as.data.frame(modsel_Ca_703[1])
res_Ca_703$adj_r2 <- summary(topmod_Ca_703)$adj.r.squared
res_Ca_703$pvalue <- lmp(topmod_Ca_703)

res_Mg_703 <- as.data.frame(modsel_Mg_703[1])
res_Mg_703$adj_r2 <- summary(topmod_Mg_703)$adj.r.squared
res_Mg_703$pvalue <- lmp(topmod_Mg_703)

res_S_703 <- as.data.frame(modsel_S_703[1])
res_S_703$adj_r2 <- summary(topmod_S_703)$adj.r.squared
res_S_703$pvalue <- lmp(topmod_S_703)

res_Zn_703 <- as.data.frame(modsel_Zn_703[1])
res_Zn_703$adj_r2 <- summary(topmod_Zn_703)$adj.r.squared
res_Zn_703$pvalue <- lmp(topmod_Zn_703)

res_703_spectal <- rbind(res_N_703, res_K_703, res_Ca_703, res_Mg_703, res_S_703)
res_703_spectal$dependent <- c("N", "K", "Ca", "Mg", "S")
write.csv(res_703_spectal, "results/res_703_spectralOnly.csv")


#---------------------------------------------------Repeat for 717 -----------------------------------------------

lm_N_717 <- lm(N ~ mean+sd+crrmedian+crrsd+kurt+moran+watersum_ln, data = norm_717, na.action=na.fail)
vif(lm_N_717)
modsel_N_717 <- dredge(lm_N_717)
topmod_N_717 <- lm(N ~ sd+crrsd+kurt+mean, data = norm_717,na.action=na.fail)
summary(topmod_N_717)
N_717_pred <- predict(topmod_N_717)
N_717_res <- resid(topmod_N_717)
n_717_table <- cbind(norm_717$N, N_717_pred, N_717_res)
head(n_717_table)

lm_P_717 <- lm(P ~ mean+sd+crrmedian+crrsd+kurt+moran+watersum_ln, data = norm_717, na.action=na.fail)
summary(lm_P_717)
modsel_P_717 <- dredge(lm_P_717)
topmod_P_717 <- lm(P ~ mean, data = norm_717,na.action=na.fail)
summary(topmod_P_717)

lm_K_717 <- lm(K ~ mean+sd+crrmedian+crrsd+kurt+moran+watersum_ln, data = norm_717, na.action=na.fail)
summary(lm_K_717)
modsel_K_717 <- dredge(lm_K_717)
topmod_K_717 <- lm(K ~ crrmedian+crrsd, data = norm_717,na.action=na.fail)
summary(topmod_K_717)

lm_Ca_717 <- lm(Ca ~ mean+sd+crrmedian+crrsd+kurt+moran+watersum_ln, data = norm_717, na.action=na.fail)
summary(lm_Ca_717)
modsel_Ca_717 <- dredge(lm_Ca_717)
topmod_Ca_717 <- lm(Ca ~ crrmedian+kurt, data = norm_717,na.action=na.fail)
summary(topmod_Ca_717)

lm_S_717 <- lm(S ~ mean+sd+crrmedian+crrsd+kurt+moran+watersum_ln, data = norm_717, na.action=na.fail)
summary(lm_S_717)
modsel_S_717 <- dredge(lm_S_717)
topmod_S_717 <- lm(S ~ mean, data = norm_717,na.action=na.fail)
summary(topmod_S_717)

lm_B_717 <- lm(B ~ mean+sd+crrmedian+crrsd+kurt+moran+watersum_ln, data = norm_717, na.action=na.fail)
summary(lm_B_717)
modsel_B_717 <- dredge(lm_B_717)
topmod_B_717 <- lm(B ~ crrmedian+sd, data = norm_717,na.action=na.fail)
summary(topmod_B_717)

lm_Zn_717 <- lm(Zn ~ mean+sd+crrmedian+crrsd+kurt+moran+watersum_ln, data = norm_717, na.action=na.fail)
summary(lm_Zn_717)
modsel_Zn_717 <- dredge(lm_Zn_717)
topmod_Zn_717 <- lm(Zn ~ kurt, data = norm_717,na.action=na.fail)
summary(topmod_Zn_717)

lm_Mg_717 <- lm(Mg ~ mean+sd+crrmedian+crrsd+kurt+moran+watersum_ln, data = norm_717, na.action=na.fail)
summary(lm_Mg_717)
modsel_Mg_717 <- dredge(lm_Mg_717)
topmod_Mg_717 <- lm(Mg ~ kurt, data = norm_717,na.action=na.fail)
summary(topmod_Mg_717)

# function to get lm pvalue
lmp <- function (modelobject) {
  if (class(modelobject) != "lm") stop("Not an object of class 'lm' ")
  f <- summary(modelobject)$fstatistic
  p <- pf(f[1],f[2],f[3],lower.tail=F)
  attributes(p) <- NULL
  return(p)
}

res_N_717 <- as.data.frame(modsel_N_717[1])
res_N_717$adj_r2 <- summary(topmod_N_717)$adj.r.squared
res_N_717$pvalue <- lmp(topmod_N_717)

res_P_717 <- as.data.frame(modsel_P_717[1])
res_P_717$adj_r2 <- summary(topmod_P_717)$adj.r.squared
res_P_717$pvalue <- lmp(topmod_P_717)

res_K_717 <- as.data.frame(modsel_K_717[1])
res_K_717$adj_r2 <- summary(topmod_K_717)$adj.r.squared
res_K_717$pvalue <- lmp(topmod_K_717)

res_B_717 <- as.data.frame(modsel_B_717[1])
res_B_717$adj_r2 <- summary(topmod_B_717)$adj.r.squared
res_B_717$pvalue <- lmp(topmod_B_717)

res_Ca_717 <- as.data.frame(modsel_Ca_717[1])
res_Ca_717$adj_r2 <- summary(topmod_Ca_717)$adj.r.squared
res_Ca_717$pvalue <- lmp(topmod_Ca_717)

res_Mg_717 <- as.data.frame(modsel_Mg_717[1])
res_Mg_717$adj_r2 <- summary(topmod_Mg_717)$adj.r.squared
res_Mg_717$pvalue <- lmp(topmod_Mg_717)

res_S_717 <- as.data.frame(modsel_S_717[1])
res_S_717$adj_r2 <- summary(topmod_S_717)$adj.r.squared
res_S_717$pvalue <- lmp(topmod_S_717)

res_Zn_717 <- as.data.frame(modsel_Zn_717[1])
res_Zn_717$adj_r2 <- summary(topmod_Zn_717)$adj.r.squared
res_Zn_717$pvalue <- lmp(topmod_Zn_717)

res_717 <- rbind(res_N_717, res_P_717, res_K_717, res_B_717, res_Ca_717, res_Mg_717, res_S_717, res_Zn_717)
res_717$dependent <- c("N", "P", "K", "B", "Ca", "Mg", "S", "Zn")
write.csv(res_717, "results/res_717_struc.csv")

#------------------------------------------- LM spectral and structural --------------------------------------------------------

lm_N_717 <- lm(N ~ sd+crrmedian+crrsd+kurt+moran+watersum_ln+tgi_median+tgi_sd_ln, data = norm_717, na.action=na.fail)
vif(lm_N_717)
modsel_N_717 <- dredge(lm_N_717)
topmod_N_717 <- lm(N ~ tgi_median+tgi_sd_ln, data = norm_717,na.action=na.fail)
summary(topmod_N_717)
N_717_pred <- predict(topmod_N_717)
N_717_res <- resid(topmod_N_717)
n_717_table <- cbind(norm_717$N, N_717_pred, N_717_res)
head(n_717_table)

lm_P_717 <- lm(P ~ sd+crrmedian+crrsd+kurt+moran+watersum_ln+tgi_median+tgi_sd_ln, data = norm_717, na.action=na.fail)
summary(lm_P_717)
modsel_P_717 <- dredge(lm_P_717)
topmod_P_717 <- lm(P ~ 1, data = norm_717,na.action=na.fail)
summary(topmod_P_717)

lm_K_717 <- lm(K ~ sd+crrmedian+crrsd+kurt+moran+watersum_ln+tgi_median+tgi_sd_ln, data = norm_717, na.action=na.fail)
summary(lm_K_717)
modsel_K_717 <- dredge(lm_K_717)
topmod_K_717 <- lm(K ~ crrmedian+tgi_sd_ln, data = norm_717,na.action=na.fail)
summary(topmod_K_717)

lm_Ca_717 <- lm(Ca ~ sd+crrmedian+crrsd+kurt+moran+watersum_ln+tgi_median+tgi_sd_ln, data = norm_717, na.action=na.fail)
summary(lm_Ca_717)
modsel_Ca_717 <- dredge(lm_Ca_717)
topmod_Ca_717 <- lm(Ca ~ crrmedian+kurt, data = norm_717,na.action=na.fail)
summary(topmod_Ca_717)

lm_S_717 <- lm(S ~ sd+crrmedian+crrsd+kurt+moran+watersum_ln+tgi_median+tgi_sd_ln, data = norm_717, na.action=na.fail)
summary(lm_S_717)
modsel_S_717 <- dredge(lm_S_717)
topmod_S_717 <- lm(S ~ tgi_sd_ln, data = norm_717,na.action=na.fail)
summary(topmod_S_717)

lm_B_717 <- lm(B ~ sd+crrmedian+crrsd+kurt+moran+watersum_ln+tgi_median+tgi_sd_ln, data = norm_717, na.action=na.fail)
summary(lm_B_717)
modsel_B_717 <- dredge(lm_B_717)
topmod_B_717 <- lm(B ~ crrmedian+sd, data = norm_717,na.action=na.fail)
summary(topmod_B_717)

lm_Zn_717 <- lm(Zn ~ sd+crrmedian+crrsd+kurt+moran+watersum_ln+tgi_median+tgi_sd_ln, data = norm_717, na.action=na.fail)
summary(lm_Zn_717)
modsel_Zn_717 <- dredge(lm_Zn_717)
topmod_Zn_717 <- lm(Zn ~ kurt, data = norm_717,na.action=na.fail)
summary(topmod_Zn_717)

lm_Mg_717 <- lm(Mg ~ sd+crrmedian+crrsd+kurt+moran+watersum_ln+tgi_median+tgi_sd_ln, data = norm_717, na.action=na.fail)
summary(lm_Mg_717)
modsel_Mg_717 <- dredge(lm_Mg_717)
topmod_Mg_717 <- lm(Mg ~ tgi_median, data = norm_717,na.action=na.fail)
summary(topmod_Mg_717)

res_N_717 <- as.data.frame(modsel_N_717[1])
res_N_717$adj_r2 <- summary(topmod_N_717)$adj.r.squared
res_N_717$pvalue <- lmp(topmod_N_717)

res_P_717 <- as.data.frame(modsel_P_717[1])
res_P_717$adj_r2 <- summary(topmod_P_717)$adj.r.squared
res_P_717$pvalue <- lmp(topmod_P_717)

res_K_717 <- as.data.frame(modsel_K_717[1])
res_K_717$adj_r2 <- summary(topmod_K_717)$adj.r.squared
res_K_717$pvalue <- lmp(topmod_K_717)

res_B_717 <- as.data.frame(modsel_B_717[1])
res_B_717$adj_r2 <- summary(topmod_B_717)$adj.r.squared
res_B_717$pvalue <- lmp(topmod_B_717)

res_Ca_717 <- as.data.frame(modsel_Ca_717[1])
res_Ca_717$adj_r2 <- summary(topmod_Ca_717)$adj.r.squared
res_Ca_717$pvalue <- lmp(topmod_Ca_717)

res_Mg_717 <- as.data.frame(modsel_Mg_717[1])
res_Mg_717$adj_r2 <- summary(topmod_Mg_717)$adj.r.squared
res_Mg_717$pvalue <- lmp(topmod_Mg_717)

res_S_717 <- as.data.frame(modsel_S_717[1])
res_S_717$adj_r2 <- summary(topmod_S_717)$adj.r.squared
res_S_717$pvalue <- lmp(topmod_S_717)

res_Zn_717 <- as.data.frame(modsel_Zn_717[1])
res_Zn_717$adj_r2 <- summary(topmod_Zn_717)$adj.r.squared
res_Zn_717$pvalue <- lmp(topmod_Zn_717)

res_717_spec_struc <- rbind(res_N_717, res_K_717, res_B_717, res_Ca_717, res_Mg_717, res_S_717, res_Zn_717)
res_717_spec_struc$dependent <- c("N", "K", "B", "Ca", "Mg", "S", "Zn")
write.csv(res_717_spec_struc, "results/res_717_spec_struc.csv")

#------------------------------------------- LM with spectral only --------------------------------------------------------
lm_N_717 <- lm(N ~ tgi_mean+tgi_iqr_ln+watersum_ln, data = norm_717, na.action=na.fail)
vif(lm_N_717)
modsel_N_717 <- dredge(lm_N_717)
topmod_N_717 <- lm(N ~ tgi_mean+tgi_iqr_ln, data = norm_717,na.action=na.fail)
summary(topmod_N_717)
N_717_pred <- predict(topmod_N_717)
N_717_res <- resid(topmod_N_717)
n_717_table <- cbind(norm_717$N, N_717_pred, N_717_res)
head(n_717_table)

lm_P_717 <- lm(P ~tgi_mean+tgi_iqr_ln+watersum_ln, data = norm_717, na.action=na.fail)
summary(lm_P_717)
modsel_P_717 <- dredge(lm_P_717)
topmod_P_717 <- lm(P ~ 1, data = norm_717,na.action=na.fail)
summary(topmod_P_717)

lm_K_717 <- lm(K ~tgi_mean+tgi_iqr_ln+watersum_ln, data = norm_717, na.action=na.fail)
summary(lm_K_717)
modsel_K_717 <- dredge(lm_K_717)
topmod_K_717 <- lm(K ~ tgi_iqr_ln, data = norm_717,na.action=na.fail)
summary(topmod_K_717)

lm_Ca_717 <- lm(Ca ~tgi_mean+tgi_iqr_ln+watersum_ln, data = norm_717, na.action=na.fail)
summary(lm_Ca_717)
modsel_Ca_717 <- dredge(lm_Ca_717)
topmod_Ca_717 <- lm(Ca ~ watersum_ln, data = norm_717,na.action=na.fail)
summary(topmod_Ca_717)

lm_S_717 <- lm(S ~tgi_mean+tgi_iqr_ln+watersum_ln, data = norm_717, na.action=na.fail)
summary(lm_S_717)
modsel_S_717 <- dredge(lm_S_717)
topmod_S_717 <- lm(S ~ tgi_iqr_ln, data = norm_717,na.action=na.fail)
summary(topmod_S_717)

lm_B_717 <- lm(B ~tgi_mean+tgi_iqr_ln+watersum_ln, data = norm_717, na.action=na.fail)
summary(lm_B_717)
modsel_B_717 <- dredge(lm_B_717)
topmod_B_717 <- lm(B ~ 1, data = norm_717,na.action=na.fail)
summary(topmod_B_717)

lm_Zn_717 <- lm(Zn ~tgi_mean+tgi_iqr_ln+watersum_ln, data = norm_717, na.action=na.fail)
summary(lm_Zn_717)
modsel_Zn_717 <- dredge(lm_Zn_717)
topmod_Zn_717 <- lm(Zn ~ 1, data = norm_717,na.action=na.fail)
summary(topmod_Zn_717)

lm_Mg_717 <- lm(Mg ~tgi_mean+tgi_iqr_ln+watersum_ln, data = norm_717, na.action=na.fail)
summary(lm_Mg_717)
modsel_Mg_717 <- dredge(lm_Mg_717)
topmod_Mg_717 <- lm(Mg ~ tgi_mean, data = norm_717,na.action=na.fail)
summary(topmod_Mg_717)

res_N_717 <- as.data.frame(modsel_N_717[1])
res_N_717$adj_r2 <- summary(topmod_N_717)$adj.r.squared
res_N_717$pvalue <- lmp(topmod_N_717)

res_P_717 <- as.data.frame(modsel_P_717[1])
res_P_717$adj_r2 <- summary(topmod_P_717)$adj.r.squared
res_P_717$pvalue <- lmp(topmod_P_717)

res_K_717 <- as.data.frame(modsel_K_717[1])
res_K_717$adj_r2 <- summary(topmod_K_717)$adj.r.squared
res_K_717$pvalue <- lmp(topmod_K_717)

res_B_717 <- as.data.frame(modsel_B_717[1])
res_B_717$adj_r2 <- summary(topmod_B_717)$adj.r.squared
res_B_717$pvalue <- lmp(topmod_B_717)

res_Ca_717 <- as.data.frame(modsel_Ca_717[1])
res_Ca_717$adj_r2 <- summary(topmod_Ca_717)$adj.r.squared
res_Ca_717$pvalue <- lmp(topmod_Ca_717)

res_Mg_717 <- as.data.frame(modsel_Mg_717[1])
res_Mg_717$adj_r2 <- summary(topmod_Mg_717)$adj.r.squared
res_Mg_717$pvalue <- lmp(topmod_Mg_717)

res_S_717 <- as.data.frame(modsel_S_717[1])
res_S_717$adj_r2 <- summary(topmod_S_717)$adj.r.squared
res_S_717$pvalue <- lmp(topmod_S_717)

res_Zn_717 <- as.data.frame(modsel_Zn_717[1])
res_Zn_717$adj_r2 <- summary(topmod_Zn_717)$adj.r.squared
res_Zn_717$pvalue <- lmp(topmod_Zn_717)

res_717_spectal <- rbind(res_N_717, res_K_717, res_Ca_717, res_Mg_717, res_S_717)
res_717_spectal$dependent <- c("N", "K", "Ca", "Mg", "S")
write.csv(res_717_spectal, "results/res_717_spectralOnly.csv")


#---------------------------------------------------Use first nutrient data for 703 and 717 -----------------------------------------------
# ----------------------------------------- 703 struct only -------------------------------------
lm_N_703 <- lm(N ~ median+crrmedian+crrsd+kurt_ln+moran+rumple+watersum_ln, data = norm_703, na.action=na.fail)
vif(lm_N_703)
modsel_N_703 <- dredge(lm_N_703)
topmod_N_703 <- lm(N ~ moran+median+watersum_ln, data = norm_703,na.action=na.fail)
summary(topmod_N_703)
N_703_pred <- predict(topmod_N_703)
N_703_res <- resid(topmod_N_703)
n_703_table <- cbind(norm_703$N, N_703_pred, N_703_res)
head(n_703_table)

lm_P_703 <- lm(P ~ median+crrmedian+crrsd+kurt_ln+moran+rumple+watersum_ln, data = norm_703, na.action=na.fail)
summary(lm_P_703)
modsel_P_703 <- dredge(lm_P_703)
topmod_P_703 <- lm(P ~ median+moran+watersum_ln, data = norm_703,na.action=na.fail)
summary(topmod_P_703)

lm_K_703 <- lm(K ~ median+crrmedian+crrsd+kurt_ln+moran+rumple+watersum_ln, data = norm_703, na.action=na.fail)
summary(lm_K_703)
modsel_K_703 <- dredge(lm_K_703)
topmod_K_703 <- lm(K ~ crrmedian, data = norm_703,na.action=na.fail)
summary(topmod_K_703)

lm_Ca_703 <- lm(Ca ~ median+crrmedian+crrsd+kurt_ln+moran+rumple+watersum_ln, data = norm_703, na.action=na.fail)
summary(lm_Ca_703)
modsel_Ca_703 <- dredge(lm_Ca_703)
topmod_Ca_703 <- lm(Ca ~ crrmedian+moran+watersum_ln, data = norm_703,na.action=na.fail)
summary(topmod_Ca_703)

lm_S_703 <- lm(S ~ median+crrmedian+crrsd+kurt_ln+moran+rumple+watersum_ln, data = norm_703, na.action=na.fail)
summary(lm_S_703)
modsel_S_703 <- dredge(lm_S_703)
topmod_S_703 <- lm(S ~ crrmedian+rumple, data = norm_703,na.action=na.fail)
summary(topmod_S_703)

lm_B_703 <- lm(B ~ median+crrmedian+crrsd+kurt_ln+moran+rumple+watersum_ln, data = norm_703, na.action=na.fail)
summary(lm_B_703)
modsel_B_703 <- dredge(lm_B_703)
topmod_B_703 <- lm(B ~ crrmedian+rumple, data = norm_703,na.action=na.fail)
summary(topmod_B_703)

lm_Zn_703 <- lm(Zn ~ median+crrmedian+crrsd+kurt_ln+moran+rumple+watersum_ln, data = norm_703, na.action=na.fail)
summary(lm_Zn_703)
modsel_Zn_703 <- dredge(lm_Zn_703)
topmod_Zn_703 <- lm(Zn ~ crrmedian+median, data = norm_703,na.action=na.fail)
summary(topmod_Zn_703)

lm_Mg_703 <- lm(Mg ~ median+crrmedian+crrsd+kurt_ln+moran+rumple+watersum_ln, data = norm_703, na.action=na.fail)
summary(lm_Mg_703)
modsel_Mg_703 <- dredge(lm_Mg_703)
topmod_Mg_703 <- lm(Mg ~ median+moran+watersum_ln, data = norm_703,na.action=na.fail)
summary(topmod_Mg_703)

# function to get lm pvalue
lmp <- function (modelobject) {
  if (class(modelobject) != "lm") stop("Not an object of class 'lm' ")
  f <- summary(modelobject)$fstatistic
  p <- pf(f[1],f[2],f[3],lower.tail=F)
  attributes(p) <- NULL
  return(p)
}

res_N_703 <- as.data.frame(modsel_N_703[1])
res_N_703$adj_r2 <- summary(topmod_N_703)$adj.r.squared
res_N_703$pvalue <- lmp(topmod_N_703)

res_P_703 <- as.data.frame(modsel_P_703[1])
res_P_703$adj_r2 <- summary(topmod_P_703)$adj.r.squared
res_P_703$pvalue <- lmp(topmod_P_703)

res_K_703 <- as.data.frame(modsel_K_703[1])
res_K_703$adj_r2 <- summary(topmod_K_703)$adj.r.squared
res_K_703$pvalue <- lmp(topmod_K_703)

res_B_703 <- as.data.frame(modsel_B_703[1])
res_B_703$adj_r2 <- summary(topmod_B_703)$adj.r.squared
res_B_703$pvalue <- lmp(topmod_B_703)

res_Ca_703 <- as.data.frame(modsel_Ca_703[1])
res_Ca_703$adj_r2 <- summary(topmod_Ca_703)$adj.r.squared
res_Ca_703$pvalue <- lmp(topmod_Ca_703)

res_Mg_703 <- as.data.frame(modsel_Mg_703[1])
res_Mg_703$adj_r2 <- summary(topmod_Mg_703)$adj.r.squared
res_Mg_703$pvalue <- lmp(topmod_Mg_703)

res_S_703 <- as.data.frame(modsel_S_703[1])
res_S_703$adj_r2 <- summary(topmod_S_703)$adj.r.squared
res_S_703$pvalue <- lmp(topmod_S_703)

res_Zn_703 <- as.data.frame(modsel_Zn_703[1])
res_Zn_703$adj_r2 <- summary(topmod_Zn_703)$adj.r.squared
res_Zn_703$pvalue <- lmp(topmod_Zn_703)

res_703 <- rbind(res_N_703, res_P_703, res_K_703, res_B_703, res_Ca_703, res_S_703, res_Mg_703)
res_703$dependent <- c("N", "P", "K", "B", "Ca", "S", "Mg")
write.csv(res_703, "results/res_703_struc_nut1.csv")

#------------------------------------------- LM with spectral and structural --------------------------------------------------------
lm_N_703 <- lm(N ~ median+sd+crrsd+kurt_ln+moran+watersum_ln+tgi_median_ln+tgi_sd_ln, data = norm_703, na.action=na.fail)
vif(lm_N_703)
modsel_N_703 <- dredge(lm_N_703)
topmod_N_703 <- lm(N ~ median+tgi_median_ln+watersum_ln, data = norm_703,na.action=na.fail)
summary(topmod_N_703)
N_703_pred <- predict(topmod_N_703)
N_703_res <- resid(topmod_N_703)
n_703_table <- cbind(norm_703$N, N_703_pred, N_703_res)
head(n_703_table)

lm_P_703 <- lm(P ~ median+sd+crrsd+kurt_ln+moran+watersum_ln+tgi_median_ln+tgi_sd_ln, data = norm_703, na.action=na.fail)
summary(lm_P_703)
modsel_P_703 <- dredge(lm_P_703)
topmod_P_703 <- lm(P ~ median+tgi_median_ln+watersum_ln, data = norm_703,na.action=na.fail)
summary(topmod_P_703)

lm_K_703 <- lm(K ~ median+sd+crrsd+kurt_ln+moran+watersum_ln+tgi_median_ln+tgi_sd_ln, data = norm_703, na.action=na.fail)
summary(lm_K_703)
modsel_K_703 <- dredge(lm_K_703)
topmod_K_703 <- lm(K ~ tgi_median_ln, data = norm_703,na.action=na.fail)
summary(topmod_K_703)

lm_Ca_703 <- lm(Ca ~ median+sd+crrsd+kurt_ln+moran+watersum_ln+tgi_median_ln+tgi_sd_ln, data = norm_703, na.action=na.fail)
summary(lm_Ca_703)
modsel_Ca_703 <- dredge(lm_Ca_703)
topmod_Ca_703 <- lm(Ca ~ moran+tgi_median_ln+watersum_ln, data = norm_703,na.action=na.fail)
summary(topmod_Ca_703)

lm_S_703 <- lm(S ~ median+sd+crrsd+kurt_ln+moran+watersum_ln+tgi_median_ln+tgi_sd_ln, data = norm_703, na.action=na.fail)
summary(lm_S_703)
modsel_S_703 <- dredge(lm_S_703)
topmod_S_703 <- lm(S ~ median+tgi_median_ln, data = norm_703,na.action=na.fail)
summary(topmod_S_703)

lm_B_703 <- lm(B ~ median+sd+crrsd+kurt_ln+moran+watersum_ln+tgi_median_ln+tgi_sd_ln, data = norm_703, na.action=na.fail)
summary(lm_B_703)
modsel_B_703 <- dredge(lm_B_703)
topmod_B_703 <- lm(B ~ median, data = norm_703,na.action=na.fail)
summary(topmod_B_703)

lm_Zn_703 <- lm(Zn ~ median+sd+crrsd+kurt_ln+moran+watersum_ln+tgi_median_ln+tgi_sd_ln, data = norm_703, na.action=na.fail)
summary(lm_Zn_703)
modsel_Zn_703 <- dredge(lm_Zn_703)
topmod_Zn_703 <- lm(Zn ~ median+tgi_median_ln+tgi_sd_ln, data = norm_703,na.action=na.fail)
summary(topmod_Zn_703)

lm_Mg_703 <- lm(Mg ~ median+sd+crrsd+kurt_ln+moran+watersum_ln+tgi_median_ln+tgi_sd_ln, data = norm_703, na.action=na.fail)
summary(lm_Mg_703)
modsel_Mg_703 <- dredge(lm_Mg_703)
topmod_Mg_703 <- lm(Mg ~ median+tgi_median_ln, data = norm_703,na.action=na.fail)
summary(topmod_Mg_703)

res_N_703 <- as.data.frame(modsel_N_703[1])
res_N_703$adj_r2 <- summary(topmod_N_703)$adj.r.squared
res_N_703$pvalue <- lmp(topmod_N_703)

res_P_703 <- as.data.frame(modsel_P_703[1])
res_P_703$adj_r2 <- summary(topmod_P_703)$adj.r.squared
res_P_703$pvalue <- lmp(topmod_P_703)

res_K_703 <- as.data.frame(modsel_K_703[1])
res_K_703$adj_r2 <- summary(topmod_K_703)$adj.r.squared
res_K_703$pvalue <- lmp(topmod_K_703)

res_B_703 <- as.data.frame(modsel_B_703[1])
res_B_703$adj_r2 <- summary(topmod_B_703)$adj.r.squared
res_B_703$pvalue <- lmp(topmod_B_703)

res_Ca_703 <- as.data.frame(modsel_Ca_703[1])
res_Ca_703$adj_r2 <- summary(topmod_Ca_703)$adj.r.squared
res_Ca_703$pvalue <- lmp(topmod_Ca_703)

res_Mg_703 <- as.data.frame(modsel_Mg_703[1])
res_Mg_703$adj_r2 <- summary(topmod_Mg_703)$adj.r.squared
res_Mg_703$pvalue <- lmp(topmod_Mg_703)

res_S_703 <- as.data.frame(modsel_S_703[1])
res_S_703$adj_r2 <- summary(topmod_S_703)$adj.r.squared
res_S_703$pvalue <- lmp(topmod_S_703)

res_Zn_703 <- as.data.frame(modsel_Zn_703[1])
res_Zn_703$adj_r2 <- summary(topmod_Zn_703)$adj.r.squared
res_Zn_703$pvalue <- lmp(topmod_Zn_703)

res_703_spec_struc <- rbind(res_N_703, res_P_703, res_K_703, res_B_703, res_Ca_703, res_Mg_703, res_S_703, res_Zn_703)
res_703_spec_struc$dependent <- c("N", "P", "K", "B", "Ca", "Mg", "S", "Zn")
write.csv(res_703_spec_struc, "results/res_703_spec_struc_nut1.csv")


#------------------------------------------- LM with spectral only --------------------------------------------------------
lm_N_703 <- lm(N ~tgi_mean_ln+tgi_sd_ln+watersum_ln, data = norm_703, na.action=na.fail)
vif(lm_N_703)
modsel_N_703 <- dredge(lm_N_703)
topmod_N_703 <- lm(N ~ tgi_mean_ln+watersum_ln, data = norm_703,na.action=na.fail)
summary(topmod_N_703)
N_703_pred <- predict(topmod_N_703)
N_703_res <- resid(topmod_N_703)
n_703_table <- cbind(norm_703$N, N_703_pred, N_703_res)
head(n_703_table)

lm_P_703 <- lm(P ~tgi_mean_ln+tgi_sd_ln+watersum_ln, data = norm_703, na.action=na.fail)
summary(lm_P_703)
modsel_P_703 <- dredge(lm_P_703)
topmod_P_703 <- lm(P ~ tgi_mean_ln+watersum_ln, data = norm_703,na.action=na.fail)
summary(topmod_P_703)

lm_K_703 <- lm(K ~tgi_mean_ln+tgi_sd_ln+watersum_ln, data = norm_703, na.action=na.fail)
summary(lm_K_703)
modsel_K_703 <- dredge(lm_K_703)
topmod_K_703 <- lm(K ~ tgi_mean_ln, data = norm_703,na.action=na.fail)
summary(topmod_K_703)

lm_Ca_703 <- lm(Ca ~tgi_mean_ln+tgi_sd_ln+watersum_ln, data = norm_703, na.action=na.fail)
summary(lm_Ca_703)
modsel_Ca_703 <- dredge(lm_Ca_703)
topmod_Ca_703 <- lm(Ca ~ watersum_ln, data = norm_703,na.action=na.fail)
summary(topmod_Ca_703)

lm_S_703 <- lm(S ~tgi_mean_ln+tgi_sd_ln+watersum_ln, data = norm_703, na.action=na.fail)
summary(lm_S_703)
modsel_S_703 <- dredge(lm_S_703)
topmod_S_703 <- lm(S ~ tgi_mean_ln+tgi_sd_ln, data = norm_703,na.action=na.fail)
summary(topmod_S_703)

lm_B_703 <- lm(B ~tgi_mean_ln+tgi_sd_ln+watersum_ln, data = norm_703, na.action=na.fail)
summary(lm_B_703)
modsel_B_703 <- dredge(lm_B_703)
topmod_B_703 <- lm(B ~ tgi_mean_ln, data = norm_703,na.action=na.fail)
summary(topmod_B_703)

lm_Zn_703 <- lm(Zn ~tgi_mean_ln+tgi_sd_ln+watersum_ln, data = norm_703, na.action=na.fail)
summary(lm_Zn_703)
modsel_Zn_703 <- dredge(lm_Zn_703)
topmod_Zn_703 <- lm(Zn ~ tgi_mean_ln+tgi_sd_ln, data = norm_703,na.action=na.fail)
summary(topmod_Zn_703)

lm_Mg_703 <- lm(Mg ~tgi_mean_ln+tgi_sd_ln+watersum_ln, data = norm_703, na.action=na.fail)
summary(lm_Mg_703)
modsel_Mg_703 <- dredge(lm_Mg_703)
topmod_Mg_703 <- lm(Mg ~ tgi_mean_ln, data = norm_703,na.action=na.fail)
summary(topmod_Mg_703)

res_N_703 <- as.data.frame(modsel_N_703[1])
res_N_703$adj_r2 <- summary(topmod_N_703)$adj.r.squared
res_N_703$pvalue <- lmp(topmod_N_703)

res_P_703 <- as.data.frame(modsel_P_703[1])
res_P_703$adj_r2 <- summary(topmod_P_703)$adj.r.squared
res_P_703$pvalue <- lmp(topmod_P_703)

res_K_703 <- as.data.frame(modsel_K_703[1])
res_K_703$adj_r2 <- summary(topmod_K_703)$adj.r.squared
res_K_703$pvalue <- lmp(topmod_K_703)

res_B_703 <- as.data.frame(modsel_B_703[1])
res_B_703$adj_r2 <- summary(topmod_B_703)$adj.r.squared
res_B_703$pvalue <- lmp(topmod_B_703)

res_Ca_703 <- as.data.frame(modsel_Ca_703[1])
res_Ca_703$adj_r2 <- summary(topmod_Ca_703)$adj.r.squared
res_Ca_703$pvalue <- lmp(topmod_Ca_703)

res_Mg_703 <- as.data.frame(modsel_Mg_703[1])
res_Mg_703$adj_r2 <- summary(topmod_Mg_703)$adj.r.squared
res_Mg_703$pvalue <- lmp(topmod_Mg_703)

res_S_703 <- as.data.frame(modsel_S_703[1])
res_S_703$adj_r2 <- summary(topmod_S_703)$adj.r.squared
res_S_703$pvalue <- lmp(topmod_S_703)

res_Zn_703 <- as.data.frame(modsel_Zn_703[1])
res_Zn_703$adj_r2 <- summary(topmod_Zn_703)$adj.r.squared
res_Zn_703$pvalue <- lmp(topmod_Zn_703)

res_703_spectal <- rbind(res_N_703, res_P_703, res_K_703, res_Ca_703, res_Mg_703, res_S_703, res_Zn_703)
res_703_spectal$dependent <- c("N", "P", "K", "Ca", "Mg", "S", "Zn")
write.csv(res_703_spectal, "results/res_703_spectralOnly_nut1.csv")

#---------------------------------------------------Repeat for 717 -----------------------------------------------

lm_N_717 <- lm(N ~ mean+sd+crrmedian+crrsd+kurt+moran+watersum_ln, data = norm_717, na.action=na.fail)
vif(lm_N_717)
modsel_N_717 <- dredge(lm_N_717)
topmod_N_717 <- lm(N ~ sd+crrsd+watersum_ln+mean, data = norm_717,na.action=na.fail)
summary(topmod_N_717)
N_717_pred <- predict(topmod_N_717)
N_717_res <- resid(topmod_N_717)
n_717_table <- cbind(norm_717$N, N_717_pred, N_717_res)
head(n_717_table)

lm_P_717 <- lm(P ~ mean+sd+crrmedian+crrsd+kurt+moran+watersum_ln, data = norm_717, na.action=na.fail)
summary(lm_P_717)
modsel_P_717 <- dredge(lm_P_717)
topmod_P_717 <- lm(P ~ sd+crrsd+watersum_ln+mean, data = norm_717,na.action=na.fail)
summary(topmod_P_717)

lm_K_717 <- lm(K ~ mean+sd+crrmedian+crrsd+kurt+moran+watersum_ln, data = norm_717, na.action=na.fail)
summary(lm_K_717)
modsel_K_717 <- dredge(lm_K_717)
topmod_K_717 <- lm(K ~ crrmedian+crrsd, data = norm_717,na.action=na.fail)
summary(topmod_K_717)

lm_Ca_717 <- lm(Ca ~ mean+sd+crrmedian+crrsd+kurt+moran+watersum_ln, data = norm_717, na.action=na.fail)
summary(lm_Ca_717)
modsel_Ca_717 <- dredge(lm_Ca_717)
topmod_Ca_717 <- lm(Ca ~ moran, data = norm_717,na.action=na.fail)
summary(topmod_Ca_717)

lm_S_717 <- lm(S ~ mean+sd+crrmedian+crrsd+kurt+moran+watersum_ln, data = norm_717, na.action=na.fail)
summary(lm_S_717)
modsel_S_717 <- dredge(lm_S_717)
topmod_S_717 <- lm(S ~ mean+crrmedian, data = norm_717,na.action=na.fail)
summary(topmod_S_717)

lm_B_717 <- lm(B ~ mean+sd+crrmedian+crrsd+kurt+moran+watersum_ln, data = norm_717, na.action=na.fail)
summary(lm_B_717)
modsel_B_717 <- dredge(lm_B_717)
topmod_B_717 <- lm(B ~ crrmedian+sd+kurt+moran, data = norm_717,na.action=na.fail)
summary(topmod_B_717)

lm_Zn_717 <- lm(Zn ~ mean+sd+crrmedian+crrsd+kurt+moran+watersum_ln, data = norm_717, na.action=na.fail)
summary(lm_Zn_717)
modsel_Zn_717 <- dredge(lm_Zn_717)
topmod_Zn_717 <- lm(Zn ~ mean+sd, data = norm_717,na.action=na.fail)
summary(topmod_Zn_717)

lm_Mg_717 <- lm(Mg ~ mean+sd+crrmedian+crrsd+kurt+moran+watersum_ln, data = norm_717, na.action=na.fail)
summary(lm_Mg_717)
modsel_Mg_717 <- dredge(lm_Mg_717)
topmod_Mg_717 <- lm(Mg ~ crrsd+kurt+mean+sd+watersum_ln, data = norm_717,na.action=na.fail)
summary(topmod_Mg_717)

# function to get lm pvalue
lmp <- function (modelobject) {
  if (class(modelobject) != "lm") stop("Not an object of class 'lm' ")
  f <- summary(modelobject)$fstatistic
  p <- pf(f[1],f[2],f[3],lower.tail=F)
  attributes(p) <- NULL
  return(p)
}

res_N_717 <- as.data.frame(modsel_N_717[1])
res_N_717$adj_r2 <- summary(topmod_N_717)$adj.r.squared
res_N_717$pvalue <- lmp(topmod_N_717)

res_P_717 <- as.data.frame(modsel_P_717[1])
res_P_717$adj_r2 <- summary(topmod_P_717)$adj.r.squared
res_P_717$pvalue <- lmp(topmod_P_717)

res_K_717 <- as.data.frame(modsel_K_717[1])
res_K_717$adj_r2 <- summary(topmod_K_717)$adj.r.squared
res_K_717$pvalue <- lmp(topmod_K_717)

res_B_717 <- as.data.frame(modsel_B_717[1])
res_B_717$adj_r2 <- summary(topmod_B_717)$adj.r.squared
res_B_717$pvalue <- lmp(topmod_B_717)

res_Ca_717 <- as.data.frame(modsel_Ca_717[1])
res_Ca_717$adj_r2 <- summary(topmod_Ca_717)$adj.r.squared
res_Ca_717$pvalue <- lmp(topmod_Ca_717)

res_Mg_717 <- as.data.frame(modsel_Mg_717[1])
res_Mg_717$adj_r2 <- summary(topmod_Mg_717)$adj.r.squared
res_Mg_717$pvalue <- lmp(topmod_Mg_717)

res_S_717 <- as.data.frame(modsel_S_717[1])
res_S_717$adj_r2 <- summary(topmod_S_717)$adj.r.squared
res_S_717$pvalue <- lmp(topmod_S_717)

res_Zn_717 <- as.data.frame(modsel_Zn_717[1])
res_Zn_717$adj_r2 <- summary(topmod_Zn_717)$adj.r.squared
res_Zn_717$pvalue <- lmp(topmod_Zn_717)

res_717 <- rbind(res_N_717, res_P_717, res_K_717, res_B_717, res_Ca_717, res_Mg_717, res_S_717, res_Zn_717)
res_717$dependent <- c("N", "P", "K", "B", "Ca", "Mg", "S", "Zn")
write.csv(res_717, "results/res_717_struc_nut1.csv")

#------------------------------------------- LM spectral and structural --------------------------------------------------------

lm_N_717 <- lm(N ~ sd+crrmedian+crrsd+kurt+moran+watersum_ln+tgi_median+tgi_sd_ln, data = norm_717, na.action=na.fail)
vif(lm_N_717)
modsel_N_717 <- dredge(lm_N_717)
topmod_N_717 <- lm(N ~ tgi_median+tgi_sd_ln+watersum_ln, data = norm_717,na.action=na.fail)
summary(topmod_N_717)
N_717_pred <- predict(topmod_N_717)
N_717_res <- resid(topmod_N_717)
n_717_table <- cbind(norm_717$N, N_717_pred, N_717_res)
head(n_717_table)

lm_P_717 <- lm(P ~ sd+crrmedian+crrsd+kurt+moran+watersum_ln+tgi_median+tgi_sd_ln, data = norm_717, na.action=na.fail)
summary(lm_P_717)
modsel_P_717 <- dredge(lm_P_717)
topmod_P_717 <- lm(P ~ tgi_sd_ln+watersum_ln, data = norm_717,na.action=na.fail)
summary(topmod_P_717)

lm_K_717 <- lm(K ~ sd+crrmedian+crrsd+kurt+moran+watersum_ln+tgi_median+tgi_sd_ln, data = norm_717, na.action=na.fail)
summary(lm_K_717)
modsel_K_717 <- dredge(lm_K_717)
topmod_K_717 <- lm(K ~ crrmedian+tgi_sd_ln, data = norm_717,na.action=na.fail)
summary(topmod_K_717)

lm_Ca_717 <- lm(Ca ~ sd+crrmedian+crrsd+kurt+moran+watersum_ln+tgi_median+tgi_sd_ln, data = norm_717, na.action=na.fail)
summary(lm_Ca_717)
modsel_Ca_717 <- dredge(lm_Ca_717)
topmod_Ca_717 <- lm(Ca ~ moran, data = norm_717,na.action=na.fail)
summary(topmod_Ca_717)

lm_S_717 <- lm(S ~ sd+crrmedian+crrsd+kurt+moran+watersum_ln+tgi_median+tgi_sd_ln, data = norm_717, na.action=na.fail)
summary(lm_S_717)
modsel_S_717 <- dredge(lm_S_717)
topmod_S_717 <- lm(S ~ crrmedian, data = norm_717,na.action=na.fail)
summary(topmod_S_717)

lm_B_717 <- lm(B ~ sd+crrmedian+crrsd+kurt+moran+watersum_ln+tgi_median+tgi_sd_ln, data = norm_717, na.action=na.fail)
summary(lm_B_717)
modsel_B_717 <- dredge(lm_B_717)
topmod_B_717 <- lm(B ~ crrmedian+kurt+moran+sd+tgi_sd_ln, data = norm_717,na.action=na.fail)
summary(topmod_B_717)

lm_Zn_717 <- lm(Zn ~ sd+crrmedian+crrsd+kurt+moran+watersum_ln+tgi_median+tgi_sd_ln, data = norm_717, na.action=na.fail)
summary(lm_Zn_717)
modsel_Zn_717 <- dredge(lm_Zn_717)
topmod_Zn_717 <- lm(Zn ~ crrmedian+kurt+tgi_median, data = norm_717,na.action=na.fail)
summary(topmod_Zn_717)

lm_Mg_717 <- lm(Mg ~ sd+crrmedian+crrsd+kurt+moran+watersum_ln+tgi_median+tgi_sd_ln, data = norm_717, na.action=na.fail)
summary(lm_Mg_717)
modsel_Mg_717 <- dredge(lm_Mg_717)
topmod_Mg_717 <- lm(Mg ~ tgi_median+tgi_sd_ln, data = norm_717,na.action=na.fail)
summary(topmod_Mg_717)

res_N_717 <- as.data.frame(modsel_N_717[1])
res_N_717$adj_r2 <- summary(topmod_N_717)$adj.r.squared
res_N_717$pvalue <- lmp(topmod_N_717)

res_P_717 <- as.data.frame(modsel_P_717[1])
res_P_717$adj_r2 <- summary(topmod_P_717)$adj.r.squared
res_P_717$pvalue <- lmp(topmod_P_717)

res_K_717 <- as.data.frame(modsel_K_717[1])
res_K_717$adj_r2 <- summary(topmod_K_717)$adj.r.squared
res_K_717$pvalue <- lmp(topmod_K_717)

res_B_717 <- as.data.frame(modsel_B_717[1])
res_B_717$adj_r2 <- summary(topmod_B_717)$adj.r.squared
res_B_717$pvalue <- lmp(topmod_B_717)

res_Ca_717 <- as.data.frame(modsel_Ca_717[1])
res_Ca_717$adj_r2 <- summary(topmod_Ca_717)$adj.r.squared
res_Ca_717$pvalue <- lmp(topmod_Ca_717)

res_Mg_717 <- as.data.frame(modsel_Mg_717[1])
res_Mg_717$adj_r2 <- summary(topmod_Mg_717)$adj.r.squared
res_Mg_717$pvalue <- lmp(topmod_Mg_717)

res_S_717 <- as.data.frame(modsel_S_717[1])
res_S_717$adj_r2 <- summary(topmod_S_717)$adj.r.squared
res_S_717$pvalue <- lmp(topmod_S_717)

res_Zn_717 <- as.data.frame(modsel_Zn_717[1])
res_Zn_717$adj_r2 <- summary(topmod_Zn_717)$adj.r.squared
res_Zn_717$pvalue <- lmp(topmod_Zn_717)

res_717_spec_struc <- rbind(res_N_717, res_P_717, res_K_717, res_B_717, res_Ca_717, res_Mg_717, res_S_717, res_Zn_717)
res_717_spec_struc$dependent <- c("N", "P", "K", "B", "Ca", "Mg", "S", "Zn")
write.csv(res_717_spec_struc, "results/res_717_spec_struc_nut1.csv")

#------------------------------------------- LM with spectral only --------------------------------------------------------
lm_N_717 <- lm(N ~ tgi_mean+tgi_iqr_ln+watersum_ln, data = norm_717, na.action=na.fail)
vif(lm_N_717)
modsel_N_717 <- dredge(lm_N_717)
topmod_N_717 <- lm(N ~ tgi_mean+tgi_iqr_ln+watersum_ln, data = norm_717,na.action=na.fail)
summary(topmod_N_717)
N_717_pred <- predict(topmod_N_717)
N_717_res <- resid(topmod_N_717)
n_717_table <- cbind(norm_717$N, N_717_pred, N_717_res)
head(n_717_table)

lm_P_717 <- lm(P ~tgi_mean+tgi_iqr_ln+watersum_ln, data = norm_717, na.action=na.fail)
summary(lm_P_717)
modsel_P_717 <- dredge(lm_P_717)
topmod_P_717 <- lm(P ~ tgi_mean+watersum_ln, data = norm_717,na.action=na.fail)
summary(topmod_P_717)

lm_K_717 <- lm(K ~tgi_mean+tgi_iqr_ln+watersum_ln, data = norm_717, na.action=na.fail)
summary(lm_K_717)
modsel_K_717 <- dredge(lm_K_717)
topmod_K_717 <- lm(K ~ tgi_mean+tgi_iqr_ln, data = norm_717,na.action=na.fail)
summary(topmod_K_717)

lm_Ca_717 <- lm(Ca ~tgi_mean+tgi_iqr_ln+watersum_ln, data = norm_717, na.action=na.fail)
summary(lm_Ca_717)
modsel_Ca_717 <- dredge(lm_Ca_717)
topmod_Ca_717 <- lm(Ca ~ tgi_mean, data = norm_717,na.action=na.fail)
summary(topmod_Ca_717)

lm_S_717 <- lm(S ~tgi_mean+tgi_iqr_ln+watersum_ln, data = norm_717, na.action=na.fail)
summary(lm_S_717)
modsel_S_717 <- dredge(lm_S_717)
topmod_S_717 <- lm(S ~ tgi_iqr_ln+watersum_ln, data = norm_717,na.action=na.fail)
summary(topmod_S_717)

lm_B_717 <- lm(B ~tgi_mean+tgi_iqr_ln+watersum_ln, data = norm_717, na.action=na.fail)
summary(lm_B_717)
modsel_B_717 <- dredge(lm_B_717)
topmod_B_717 <- lm(B ~ tgi_iqr_ln, data = norm_717,na.action=na.fail)
summary(topmod_B_717)

lm_Zn_717 <- lm(Zn ~tgi_mean+tgi_iqr_ln+watersum_ln, data = norm_717, na.action=na.fail)
summary(lm_Zn_717)
modsel_Zn_717 <- dredge(lm_Zn_717)
topmod_Zn_717 <- lm(Zn ~ tgi_iqr_ln, data = norm_717,na.action=na.fail)
summary(topmod_Zn_717)

lm_Mg_717 <- lm(Mg ~tgi_mean+tgi_iqr_ln+watersum_ln, data = norm_717, na.action=na.fail)
summary(lm_Mg_717)
modsel_Mg_717 <- dredge(lm_Mg_717)
topmod_Mg_717 <- lm(Mg ~ tgi_mean+tgi_iqr_ln, data = norm_717,na.action=na.fail)
summary(topmod_Mg_717)

res_N_717 <- as.data.frame(modsel_N_717[1])
res_N_717$adj_r2 <- summary(topmod_N_717)$adj.r.squared
res_N_717$pvalue <- lmp(topmod_N_717)

res_P_717 <- as.data.frame(modsel_P_717[1])
res_P_717$adj_r2 <- summary(topmod_P_717)$adj.r.squared
res_P_717$pvalue <- lmp(topmod_P_717)

res_K_717 <- as.data.frame(modsel_K_717[1])
res_K_717$adj_r2 <- summary(topmod_K_717)$adj.r.squared
res_K_717$pvalue <- lmp(topmod_K_717)

res_B_717 <- as.data.frame(modsel_B_717[1])
res_B_717$adj_r2 <- summary(topmod_B_717)$adj.r.squared
res_B_717$pvalue <- lmp(topmod_B_717)

res_Ca_717 <- as.data.frame(modsel_Ca_717[1])
res_Ca_717$adj_r2 <- summary(topmod_Ca_717)$adj.r.squared
res_Ca_717$pvalue <- lmp(topmod_Ca_717)

res_Mg_717 <- as.data.frame(modsel_Mg_717[1])
res_Mg_717$adj_r2 <- summary(topmod_Mg_717)$adj.r.squared
res_Mg_717$pvalue <- lmp(topmod_Mg_717)

res_S_717 <- as.data.frame(modsel_S_717[1])
res_S_717$adj_r2 <- summary(topmod_S_717)$adj.r.squared
res_S_717$pvalue <- lmp(topmod_S_717)

res_Zn_717 <- as.data.frame(modsel_Zn_717[1])
res_Zn_717$adj_r2 <- summary(topmod_Zn_717)$adj.r.squared
res_Zn_717$pvalue <- lmp(topmod_Zn_717)

res_717_spectal <- rbind(res_N_717, res_P_717, res_K_717, res_B_717, res_Ca_717, res_Mg_717, res_S_717, res_Zn_717)
res_717_spectal$dependent <- c("N", "P", "K", "B", "Ca", "Mg", "S", "Zn")
write.csv(res_717_spectal, "results/res_717_spectralOnly_nut1.csv")


#--------------------------------------------mixed effect models --------------------------------------------------

# add treatment/fert column
norm_619$group <- paste(norm_619$treatment, norm_619$fertilizer, sep="_")

mix_med_717 <- lmer(B ~ median+moran+(1|treatment), data = norm_717, na.action = na.fail)
summary(mix_med_717)
r.squaredGLMM(mix_med_717)
# Of course using treatment will better explain B levels... How is this helpful?


# ------------------------------ Try lm with nutrients as independent ---------------------------------------------

# All nutrients: 
#all nutrients: `N`+`P`+`K`+`Mg`+`Ca`+`S`+`Zn`+`Mn`+`Cu`+`Fe`+`B`+`Al`+`Na`
# nutrients and soil: "`N`+`P`+`K`+`Mg`+`Ca`+`S`+`Zn`+`Mn`+`Cu`+`Fe`+`B`+`Al`+`Na`+`soil_GoA`+`soil_Ra`+`soil_MaB`"

# # Function for calculating r^2 in SAR model, null model will be same for each date
# 
# w <- knn2nb(knearneigh(coordinates(plots_619), k=8))
# 
# null <- spautolm(median ~ 1, data = agg_619, listw=nb2listw(w), zero.policy = T, na.action = na.fail, family="SAR")
# #null <- lm(BuAc ~ 1, data=norm_1cm, na.action=na.fail)
# r2 <- function(model){
#   r.squaredLR(model, null = null)[1]
# }
# 
# adj_r2 <- function(model){
#   r2 <- r.squaredLR(model, null = null)[1]
#   n <- model$fit[[5]]
#   p <- model$parameters
#   1-(1-r2)*((n - 1)/(n - p - 1))
# }
# 
# # 619 OLS
# lm_619 <- lm(median ~ `N`+`Ca`+`S`+`Cu`+`B`+`Al`, data = norm_619, na.action = na.fail)
# summary(lm_619)
# qqnorm(lm_619$residuals)
# sd(lm_619$residuals)
# vif(lm_619) # Colinearity: removed Na, P, Fe, Mn, K
# dredge(lm_619, trace = T)
# 
# # 619 check for spatial autocorrelation
# std_619$ols_res <- lm_619$residuals
# # Weights matrix
# w <- knn2nb(knearneigh(coordinates(plots), k=8))
# moran.test(std_619$ols_res, nb2listw(w))
# 
# #Apply SAR model
# sar_619 <- spautolm(median ~ `N`+`Ca`+`S`+`Cu`+`B`+`Al`, data = std_619, listw=nb2listw(w), zero.policy = T, na.action = na.fail, family="SAR")
# summary(sar_619)
# sar_resi <- sar_619$fit[[9]]
# std_619$sar_resi <- sar_resi
# moran.test(std_619$sar_resi, nb2listw(w))
# 
# 
# lm_717 <- lm(median ~ `N`+`P`+`Ca`+`S`+`Cu`+`Fe`+`B`, data = std_717, na.action = na.fail)
# summary(lm_717)
# vif(lm_717) # Removed Na, Mg, Mn, K, Zn, Al
# dredge(lm_717, trace = T)
# 
# # 717 check for spatial autocorrelation
# std_717$ols_res <- lm_717$residuals
# moran.test(std_717$ols_res, nb2listw(w))
# 
# #Apply SAR model
# sar_717 <- spautolm(median ~ `N`+`P`+`Ca`+`S`+`Cu`+`Fe`+`B`, data = std_717, listw=nb2listw(w), zero.policy = T, na.action = na.fail, family="SAR")
# summary(sar_717)
# sar_resi <- sar_717$fit[[9]]
# std_717$sar_resi <- sar_resi
# moran.test(std_717$sar_resi, nb2listw(w))
# 
# lm_717 <- lm(median ~ `N`+`P`+`Mg`+`S`+`Mn`+`Fe`+`B`+`Al`, data = std_717, na.action = na.fail)
# summary(lm_717)
# vif(lm_717) # Eliminated Na, Cu, Ca, Zn, K
# dredge(lm_717, trace = T)
# 
# 
# # ---------------------------------SAR models with structural metric as response for each time step--------------------------------
# 
# r2 <- function(model){
#   r.squaredLR(model, null = null)[1]
# }
# 
# adj_r2 <- function(model){
#   r2 <- r.squaredLR(model, null = null)[1]
#   n <- model$fit[[5]]
#   p <- model$parameters
#   1-(1-r2)*((n - 1)/(n - p - 1))
# }
# 
# w <- knn2nb(knearneigh(coordinates(plots_619), k=8))
# 
# # ---------------------------------------SAR 619--------------------------------------------
# variables <- "`N`+`Ca`+`S`+`Cu`+`B`+`Al`+`soil_GoA`+`soil_Ra`+`soil_MaB`"
# dat <- std_619
# 
# 
# metric <- "median"
# null <- spautolm( as.formula(paste(metric, 1, sep="~")), data = dat, 
#                   listw=nb2listw(w), zero.policy = T, na.action = na.fail, family="SAR")
# sar <- spautolm(as.formula(paste(metric, variables, sep = " ~ ")), 
#                 data = dat, listw=nb2listw(w), zero.policy = T, na.action = na.fail, 
#                 family="SAR")
# summary(sar)
# adj_r2(sar)
# modsel <- dredge(sar, beta = "none", trace = T, extra = c('r2', 'adj_r2'))
# med_619 <- modsel[1]
# med_619$depvar <- "median"
# 
# metric <- "iqr"
# null <- spautolm( as.formula(paste(metric, 1, sep="~")), data = dat, 
#                   listw=nb2listw(w), zero.policy = T, na.action = na.fail, family="SAR")
# sar <- spautolm(as.formula(paste(metric, variables, sep = " ~ ")), 
#                 data = dat, listw=nb2listw(w), zero.policy = T, na.action = na.fail, 
#                 family="SAR")
# summary(sar)
# adj_r2(sar)
# modsel <- dredge(sar, beta = "none", trace = T, extra = c('r2', 'adj_r2'))
# iqr_619 <- modsel[1]
# iqr_619$depvar <- "iqr"
# 
# metric <- "crrmedian"
# null <- spautolm( as.formula(paste(metric, 1, sep="~")), data = dat, 
#                   listw=nb2listw(w), zero.policy = T, na.action = na.fail, family="SAR")
# sar <- spautolm(as.formula(paste(metric, variables, sep = " ~ ")), 
#                 data = dat, listw=nb2listw(w), zero.policy = T, na.action = na.fail, 
#                 family="SAR")
# summary(sar)
# adj_r2(sar)
# modsel <- dredge(sar, beta = "none", trace = T, extra = c('r2', 'adj_r2'))
# crrmed_619 <- modsel[1]
# crrmed_619$depvar <- "crrmed"
# 
# metric <- "crriqr"
# null <- spautolm( as.formula(paste(metric, 1, sep="~")), data = dat, 
#                   listw=nb2listw(w), zero.policy = T, na.action = na.fail, family="SAR")
# sar <- spautolm(as.formula(paste(metric, variables, sep = " ~ ")), 
#                 data = dat, listw=nb2listw(w), zero.policy = T, na.action = na.fail, 
#                 family="SAR")
# summary(sar)
# adj_r2(sar)
# modsel <- dredge(sar, beta = "none", trace = T, extra = c('r2', 'adj_r2'))
# crriqr_619 <- modsel[1]
# crriqr_619$depvar <- "crriqr"
# 
# metric <- "skew"
# null <- spautolm( as.formula(paste(metric, 1, sep="~")), data = dat, 
#                   listw=nb2listw(w), zero.policy = T, na.action = na.fail, family="SAR")
# sar <- spautolm(as.formula(paste(metric, variables, sep = " ~ ")), 
#                 data = dat, listw=nb2listw(w), zero.policy = T, na.action = na.fail, 
#                 family="SAR")
# summary(sar)
# adj_r2(sar)
# modsel <- dredge(sar, beta = "none", trace = T, extra = c('r2', 'adj_r2'))
# skew_619 <- modsel[1]
# skew_619$depvar <- "skew"
# 
# metric <- "kurt"
# null <- spautolm( as.formula(paste(metric, 1, sep="~")), data = dat, 
#                   listw=nb2listw(w), zero.policy = T, na.action = na.fail, family="SAR")
# sar <- spautolm(as.formula(paste(metric, variables, sep = " ~ ")), 
#                 data = dat, listw=nb2listw(w), zero.policy = T, na.action = na.fail, 
#                 family="SAR")
# summary(sar)
# adj_r2(sar)
# modsel <- dredge(sar, beta = "none", trace = T, extra = c('r2', 'adj_r2'))
# kurt_619 <- modsel[1]
# 
# metric <- "PlotCRR"
# null <- spautolm( as.formula(paste(metric, 1, sep="~")), data = dat, 
#                   listw=nb2listw(w), zero.policy = T, na.action = na.fail, family="SAR")
# sar <- spautolm(as.formula(paste(metric, variables, sep = " ~ ")), 
#                 data = dat, listw=nb2listw(w), zero.policy = T, na.action = na.fail, 
#                 family="SAR")
# summary(sar)
# adj_r2(sar)
# modsel <- dredge(sar, beta = "none", trace = T, extra = c('r2', 'adj_r2'))
# plotcrr_619 <- modsel[1]
# 
# metric <- "rumple"
# null <- spautolm( as.formula(paste(metric, 1, sep="~")), data = dat, 
#                   listw=nb2listw(w), zero.policy = T, na.action = na.fail, family="SAR")
# sar <- spautolm(as.formula(paste(metric, variables, sep = " ~ ")), 
#                 data = dat, listw=nb2listw(w), zero.policy = T, na.action = na.fail, 
#                 family="SAR")
# summary(sar)
# adj_r2(sar)
# modsel <- dredge(sar, beta = "none", trace = T, extra = c('r2', 'adj_r2'))
# rumple_619 <- modsel[1]
# 
# metric <- "moran"
# null <- spautolm( as.formula(paste(metric, 1, sep="~")), data = dat, 
#                   listw=nb2listw(w), zero.policy = T, na.action = na.fail, family="SAR")
# sar <- spautolm(as.formula(paste(metric, variables, sep = " ~ ")), 
#                 data = dat, listw=nb2listw(w), zero.policy = T, na.action = na.fail, 
#                 family="SAR")
# summary(sar)
# adj_r2(sar)
# modsel <- dredge(sar, beta = "none", trace = T, extra = c('r2', 'adj_r2'))
# moran_619 <- modsel[1]
# 
# 
# # ---------------------------------------SAR 717--------------------------------------------
# variables <- "`N`+`P`+`Ca`+`S`+`Cu`+`Fe`+`B`+`soil_GoA`+`soil_Ra`+`soil_MaB`"
# 
# dat <- std_717
# 
# metric <- "median"
# null <- spautolm( as.formula(paste(metric, 1, sep="~")), data = dat, 
#                   listw=nb2listw(w), zero.policy = T, na.action = na.fail, family="SAR")
# sar <- spautolm(as.formula(paste(metric, variables, sep = " ~ ")), 
#                 data = dat, listw=nb2listw(w), zero.policy = T, na.action = na.fail, 
#                 family="SAR")
# summary(sar)
# adj_r2(sar)
# modsel <- dredge(sar, beta = "none", trace = T, extra = c('r2', 'adj_r2'))
# med_717 <- modsel[1]
# 
# metric <- "iqr"
# null <- spautolm( as.formula(paste(metric, 1, sep="~")), data = dat, 
#                   listw=nb2listw(w), zero.policy = T, na.action = na.fail, family="SAR")
# sar <- spautolm(as.formula(paste(metric, variables, sep = " ~ ")), 
#                 data = dat, listw=nb2listw(w), zero.policy = T, na.action = na.fail, 
#                 family="SAR")
# summary(sar)
# adj_r2(sar)
# modsel <- dredge(sar, beta = "none", trace = T, extra = c('r2', 'adj_r2'))
# iqr_717 <- modsel[1]
# 
# metric <- "crrmedian"
# null <- spautolm( as.formula(paste(metric, 1, sep="~")), data = dat, 
#                   listw=nb2listw(w), zero.policy = T, na.action = na.fail, family="SAR")
# sar <- spautolm(as.formula(paste(metric, variables, sep = " ~ ")), 
#                 data = dat, listw=nb2listw(w), zero.policy = T, na.action = na.fail, 
#                 family="SAR")
# summary(sar)
# adj_r2(sar)
# modsel <- dredge(sar, beta = "none", trace = T, extra = c('r2', 'adj_r2'))
# crrmed_717 <- modsel[1]
# 
# metric <- "crriqr"
# null <- spautolm( as.formula(paste(metric, 1, sep="~")), data = dat, 
#                   listw=nb2listw(w), zero.policy = T, na.action = na.fail, family="SAR")
# sar <- spautolm(as.formula(paste(metric, variables, sep = " ~ ")), 
#                 data = dat, listw=nb2listw(w), zero.policy = T, na.action = na.fail, 
#                 family="SAR")
# summary(sar)
# adj_r2(sar)
# modsel <- dredge(sar, beta = "none", trace = T, extra = c('r2', 'adj_r2'))
# crriqr_717 <- modsel[1]
# 
# metric <- "skew"
# null <- spautolm( as.formula(paste(metric, 1, sep="~")), data = dat, 
#                   listw=nb2listw(w), zero.policy = T, na.action = na.fail, family="SAR")
# sar <- spautolm(as.formula(paste(metric, variables, sep = " ~ ")), 
#                 data = dat, listw=nb2listw(w), zero.policy = T, na.action = na.fail, 
#                 family="SAR")
# summary(sar)
# adj_r2(sar)
# modsel <- dredge(sar, beta = "none", trace = T, extra = c('r2', 'adj_r2'))
# skew_717 <- modsel[1]
# 
# metric <- "kurt"
# null <- spautolm( as.formula(paste(metric, 1, sep="~")), data = dat, 
#                   listw=nb2listw(w), zero.policy = T, na.action = na.fail, family="SAR")
# sar <- spautolm(as.formula(paste(metric, variables, sep = " ~ ")), 
#                 data = dat, listw=nb2listw(w), zero.policy = T, na.action = na.fail, 
#                 family="SAR")
# summary(sar)
# adj_r2(sar)
# modsel <- dredge(sar, beta = "none", trace = T, extra = c('r2', 'adj_r2'))
# kurt_717 <- modsel[1]
# 
# metric <- "PlotCRR"
# null <- spautolm( as.formula(paste(metric, 1, sep="~")), data = dat, 
#                   listw=nb2listw(w), zero.policy = T, na.action = na.fail, family="SAR")
# sar <- spautolm(as.formula(paste(metric, variables, sep = " ~ ")), 
#                 data = dat, listw=nb2listw(w), zero.policy = T, na.action = na.fail, 
#                 family="SAR")
# summary(sar)
# adj_r2(sar)
# modsel <- dredge(sar, beta = "none", trace = T, extra = c('r2', 'adj_r2'))
# plotcrr_717 <- modsel[1]
# 
# metric <- "rumple"
# null <- spautolm( as.formula(paste(metric, 1, sep="~")), data = dat, 
#                   listw=nb2listw(w), zero.policy = T, na.action = na.fail, family="SAR")
# sar <- spautolm(as.formula(paste(metric, variables, sep = " ~ ")), 
#                 data = dat, listw=nb2listw(w), zero.policy = T, na.action = na.fail, 
#                 family="SAR")
# summary(sar)
# adj_r2(sar)
# modsel <- dredge(sar, beta = "none", trace = T, extra = c('r2', 'adj_r2'))
# rumple_717 <- modsel[1]
# 
# metric <- "moran"
# null <- spautolm( as.formula(paste(metric, 1, sep="~")), data = dat, 
#                   listw=nb2listw(w), zero.policy = T, na.action = na.fail, family="SAR")
# sar <- spautolm(as.formula(paste(metric, variables, sep = " ~ ")), 
#                 data = dat, listw=nb2listw(w), zero.policy = T, na.action = na.fail, 
#                 family="SAR")
# summary(sar)
# adj_r2(sar)
# modsel <- dredge(sar, beta = "none", trace = T, extra = c('r2', 'adj_r2'))
# moran_717 <- modsel[1]
# 
# # ---------------------------------------SAR 717--------------------------------------------
# variables <- "`N`+`P`+`Mg`+`S`+`Mn`+`Fe`+`B`+`Al`"
# 
# dat <- std_717
# 
# metric <- "median"
# null <- spautolm( as.formula(paste(metric, 1, sep="~")), data = dat, 
#                   listw=nb2listw(w), zero.policy = T, na.action = na.fail, family="SAR")
# sar <- spautolm(as.formula(paste(metric, variables, sep = " ~ ")), 
#                 data = dat, listw=nb2listw(w), zero.policy = T, na.action = na.fail, 
#                 family="SAR")
# summary(sar)
# adj_r2(sar)
# modsel <- dredge(sar, beta = "none", trace = T, extra = c('r2', 'adj_r2'))
# med_717 <- modsel[1]
# 
# metric <- "iqr"
# null <- spautolm( as.formula(paste(metric, 1, sep="~")), data = dat, 
#                   listw=nb2listw(w), zero.policy = T, na.action = na.fail, family="SAR")
# sar <- spautolm(as.formula(paste(metric, variables, sep = " ~ ")), 
#                 data = dat, listw=nb2listw(w), zero.policy = T, na.action = na.fail, 
#                 family="SAR")
# summary(sar)
# adj_r2(sar)
# modsel <- dredge(sar, beta = "none", trace = T, extra = c('r2', 'adj_r2'))
# iqr_717 <- modsel[1]
# 
# metric <- "crrmedian"
# null <- spautolm( as.formula(paste(metric, 1, sep="~")), data = dat, 
#                   listw=nb2listw(w), zero.policy = T, na.action = na.fail, family="SAR")
# sar <- spautolm(as.formula(paste(metric, variables, sep = " ~ ")), 
#                 data = dat, listw=nb2listw(w), zero.policy = T, na.action = na.fail, 
#                 family="SAR")
# summary(sar)
# adj_r2(sar)
# modsel <- dredge(sar, beta = "none", trace = T, extra = c('r2', 'adj_r2'))
# crrmed_717 <- modsel[1]
# 
# metric <- "crriqr"
# null <- spautolm( as.formula(paste(metric, 1, sep="~")), data = dat, 
#                   listw=nb2listw(w), zero.policy = T, na.action = na.fail, family="SAR")
# sar <- spautolm(as.formula(paste(metric, variables, sep = " ~ ")), 
#                 data = dat, listw=nb2listw(w), zero.policy = T, na.action = na.fail, 
#                 family="SAR")
# summary(sar)
# adj_r2(sar)
# modsel <- dredge(sar, beta = "none", trace = T, extra = c('r2', 'adj_r2'))
# crriqr_717 <- modsel[1]
# 
# metric <- "skew"
# null <- spautolm( as.formula(paste(metric, 1, sep="~")), data = dat, 
#                   listw=nb2listw(w), zero.policy = T, na.action = na.fail, family="SAR")
# sar <- spautolm(as.formula(paste(metric, variables, sep = " ~ ")), 
#                 data = dat, listw=nb2listw(w), zero.policy = T, na.action = na.fail, 
#                 family="SAR")
# summary(sar)
# adj_r2(sar)
# modsel <- dredge(sar, beta = "none", trace = T, extra = c('r2', 'adj_r2'))
# skew_717 <- modsel[1]
# 
# metric <- "kurt"
# null <- spautolm( as.formula(paste(metric, 1, sep="~")), data = dat, 
#                   listw=nb2listw(w), zero.policy = T, na.action = na.fail, family="SAR")
# sar <- spautolm(as.formula(paste(metric, variables, sep = " ~ ")), 
#                 data = dat, listw=nb2listw(w), zero.policy = T, na.action = na.fail, 
#                 family="SAR")
# summary(sar)
# adj_r2(sar)
# modsel <- dredge(sar, beta = "none", trace = T, extra = c('r2', 'adj_r2'))
# kurt_717 <- modsel[1]
# 
# metric <- "PlotCRR"
# null <- spautolm( as.formula(paste(metric, 1, sep="~")), data = dat, 
#                   listw=nb2listw(w), zero.policy = T, na.action = na.fail, family="SAR")
# sar <- spautolm(as.formula(paste(metric, variables, sep = " ~ ")), 
#                 data = dat, listw=nb2listw(w), zero.policy = T, na.action = na.fail, 
#                 family="SAR")
# summary(sar)
# adj_r2(sar)
# modsel <- dredge(sar, beta = "none", trace = T, extra = c('r2', 'adj_r2'))
# plotcrr_717 <- modsel[1]
# 
# metric <- "rumple"
# null <- spautolm( as.formula(paste(metric, 1, sep="~")), data = dat, 
#                   listw=nb2listw(w), zero.policy = T, na.action = na.fail, family="SAR")
# sar <- spautolm(as.formula(paste(metric, variables, sep = " ~ ")), 
#                 data = dat, listw=nb2listw(w), zero.policy = T, na.action = na.fail, 
#                 family="SAR")
# summary(sar)
# adj_r2(sar)
# modsel <- dredge(sar, beta = "none", trace = T, extra = c('r2', 'adj_r2'))
# rumple_717 <- modsel[1]
# 
# metric <- "moran"
# null <- spautolm( as.formula(paste(metric, 1, sep="~")), data = dat, 
#                   listw=nb2listw(w), zero.policy = T, na.action = na.fail, family="SAR")
# sar <- spautolm(as.formula(paste(metric, variables, sep = " ~ ")), 
#                 data = dat, listw=nb2listw(w), zero.policy = T, na.action = na.fail, 
#                 family="SAR")
# summary(sar)
# adj_r2(sar)
# modsel <- dredge(sar, beta = "none", trace = T, extra = c('r2', 'adj_r2'))
# moran_717 <- modsel[1]

#---------------------------------------SAR nutrient results tables --------------------------------------
# 
# med_619$depvar <- "median"
# iqr_619$depvar <- "iqr"
# crrmed_619$depvar <- "crrmed"
# crriqr_619$depvar <- "crriqr"
# skew_619$depvar <- "skew"
# kurt_619$depvar <- "kurt"
# plotcrr_619$depvar <- "plotcrr"
# rumple_619$depvar <- "rumple"
# moran_619$depvar <- "moran"
# 
# med_717$depvar <- "median"
# iqr_717$depvar <- "iqr"
# crrmed_717$depvar <- "crrmed"
# crriqr_717$depvar <- "crriqr"
# skew_717$depvar <- "skew"
# kurt_717$depvar <- "kurt"
# plotcrr_717$depvar <- "plotcrr"
# rumple_717$depvar <- "rumple"
# moran_717$depvar <- "moran"
# 
# med_717$depvar <- "median"
# iqr_717$depvar <- "iqr"
# crrmed_717$depvar <- "crrmed"
# crriqr_717$depvar <- "crriqr"
# skew_717$depvar <- "skew"
# kurt_717$depvar <- "kurt"
# plotcrr_717$depvar <- "plotcrr"
# rumple_717$depvar <- "rumple"
# moran_717$depvar <- "moran"
# 
# sarresults_nutr_619 <- rbind(med_619, iqr_619, crrmed_619, crriqr_619, skew_619, kurt_619, plotcrr_619, rumple_619, moran_619)
# sarresults_nutr_619$depvar <- c("median", "iqr", "crrmed", "crriqr", "skew", "kurt", "plotcrr", "rumple", "moran")
# 
# sarresults_nutr_717 <- rbind(med_717, iqr_717, crrmed_717, crriqr_717, skew_717, kurt_717, plotcrr_717, rumple_717, moran_717)
# sarresults_nutr_717$depvar <- c("median", "iqr", "crrmed", "crriqr", "skew", "kurt", "plotcrr", "rumple", "moran")
# 
# sarresults_nutr_717 <- rbind(med_717, iqr_717, crrmed_717, crriqr_717, skew_717, kurt_717, plotcrr_717, rumple_717, moran_717)
# sarresults_nutr_717$depvar <- c("median", "iqr", "crrmed", "crriqr", "skew", "kurt", "plotcrr", "rumple", "moran")
# 
# write.csv(sarresults_nutr_619, "Q:/My Drive/Research/Canopy_Morphology/Tobacco_Project/2018/canopy_analysis/sarresults_nutr_619.csv")
# write.csv(sarresults_nutr_717, "Q:/My Drive/Research/Canopy_Morphology/Tobacco_Project/2018/canopy_analysis/sarresults_nutr_717.csv")
# write.csv(sarresults_nutr_717, "Q:/My Drive/Research/Canopy_Morphology/Tobacco_Project/2018/canopy_analysis/sarresults_nutr_717.csv")
# 

# # ---------------------------------SAR models with nutrients as response for each time step--------------------------------
# 
# variables <- "crrmedian+crriqr+median+skew+iqr+kurt+rumple+moran+soil_GoA+soil_Ra+soil_MaB+watersum"
# 
# 
# 
# r2 <- function(model){
#   r.squaredLR(model, null = null)[1]
# }
# 
# adj_r2 <- function(model){
#   r2 <- r.squaredLR(model, null = null)[1]
#   n <- model$fit[[5]]
#   p <- model$parameters
#   1-(1-r2)*((n - 1)/(n - p - 1))
# }
# 
# w <- knn2nb(knearneigh(coordinates(plots_619), k=8))
# 
# # ---------------------------------------SAR 619--------------------------------------------
# dat <- std_619
# 
# 
# metric <- "N"
# null <- spautolm( as.formula(paste(metric, 1, sep="~")), data = dat, 
#                   listw=nb2listw(w), zero.policy = T, na.action = na.fail, family="SAR")
# sar <- spautolm(as.formula(paste(metric, variables, sep = " ~ ")), 
#                 data = dat, listw=nb2listw(w), zero.policy = T, na.action = na.fail, 
#                 family="SAR")
# summary(sar)
# adj_r2(sar)
# modsel <- dredge(sar, beta = "none", trace = T, extra = c('r2', 'adj_r2'))
# n_619 <- modsel[1]
# 
# tmp <-spautolm(N~median+moran, 
#          data = dat, listw=nb2listw(w), zero.policy = T, na.action = na.fail, 
#          family="SAR")
# sar_resi <- tmp$fit[[9]]
# std_619$sar_resi <- sar_resi
# plot(std_619$sar_resi)
# 
# metric <- "P"
# null <- spautolm( as.formula(paste(metric, 1, sep="~")), data = dat, 
#                   listw=nb2listw(w), zero.policy = T, na.action = na.fail, family="SAR")
# sar <- spautolm(as.formula(paste(metric, variables, sep = " ~ ")), 
#                 data = dat, listw=nb2listw(w), zero.policy = T, na.action = na.fail, 
#                 family="SAR")
# summary(sar)
# adj_r2(sar)
# modsel <- dredge(sar, beta = "none", trace = T, extra = c('r2', 'adj_r2'))
# p_619 <- modsel[1]
# 
# metric <- "K"
# null <- spautolm( as.formula(paste(metric, 1, sep="~")), data = dat, 
#                   listw=nb2listw(w), zero.policy = T, na.action = na.fail, family="SAR")
# sar <- spautolm(as.formula(paste(metric, variables, sep = " ~ ")), 
#                 data = dat, listw=nb2listw(w), zero.policy = T, na.action = na.fail, 
#                 family="SAR")
# summary(sar)
# adj_r2(sar)
# modsel <- dredge(sar, beta = "none", trace = T, extra = c('r2', 'adj_r2'))
# k_619 <- modsel[1]
# 
# 
# metric <- "Ca"
# null <- spautolm( as.formula(paste(metric, 1, sep="~")), data = dat, 
#                   listw=nb2listw(w), zero.policy = T, na.action = na.fail, family="SAR")
# sar <- spautolm(as.formula(paste(metric, variables, sep = " ~ ")), 
#                 data = dat, listw=nb2listw(w), zero.policy = T, na.action = na.fail, 
#                 family="SAR")
# summary(sar)
# adj_r2(sar)
# modsel <- dredge(sar, beta = "none", trace = T, extra = c('r2', 'adj_r2'))
# ca_619 <- modsel[1]
# 
# metric <- "Mg"
# null <- spautolm( as.formula(paste(metric, 1, sep="~")), data = dat, 
#                   listw=nb2listw(w), zero.policy = T, na.action = na.fail, family="SAR")
# sar <- spautolm(as.formula(paste(metric, variables, sep = " ~ ")), 
#                 data = dat, listw=nb2listw(w), zero.policy = T, na.action = na.fail, 
#                 family="SAR")
# summary(sar)
# adj_r2(sar)
# modsel <- dredge(sar, beta = "none", trace = T, extra = c('r2', 'adj_r2'))
# mg_619 <- modsel[1]
# 
# metric <- "B"
# null <- spautolm( as.formula(paste(metric, 1, sep="~")), data = dat, 
#                   listw=nb2listw(w), zero.policy = T, na.action = na.fail, family="SAR")
# sar <- spautolm(as.formula(paste(metric, variables, sep = " ~ ")), 
#                 data = dat, listw=nb2listw(w), zero.policy = T, na.action = na.fail, 
#                 family="SAR")
# summary(sar)
# adj_r2(sar)
# modsel <- dredge(sar, beta = "none", trace = T, extra = c('r2', 'adj_r2'))
# b_619 <- modsel[1]
# 
# metric <- "Mn"
# null <- spautolm( as.formula(paste(metric, 1, sep="~")), data = dat, 
#                   listw=nb2listw(w), zero.policy = T, na.action = na.fail, family="SAR")
# sar <- spautolm(as.formula(paste(metric, variables, sep = " ~ ")), 
#                 data = dat, listw=nb2listw(w), zero.policy = T, na.action = na.fail, 
#                 family="SAR")
# summary(sar)
# adj_r2(sar)
# modsel <- dredge(sar, beta = "none", trace = T, extra = c('r2', 'adj_r2'))
# mn_619 <- modsel[1]
# 
# metric <- "S"
# null <- spautolm( as.formula(paste(metric, 1, sep="~")), data = dat, 
#                   listw=nb2listw(w), zero.policy = T, na.action = na.fail, family="SAR")
# sar <- spautolm(as.formula(paste(metric, variables, sep = " ~ ")), 
#                 data = dat, listw=nb2listw(w), zero.policy = T, na.action = na.fail, 
#                 family="SAR")
# summary(sar)
# adj_r2(sar)
# modsel <- dredge(sar, beta = "none", trace = T, extra = c('r2', 'adj_r2'))
# s_619 <- modsel[1]
# 
# metric <- "Zn"
# null <- spautolm( as.formula(paste(metric, 1, sep="~")), data = dat, 
#                   listw=nb2listw(w), zero.policy = T, na.action = na.fail, family="SAR")
# sar <- spautolm(as.formula(paste(metric, variables, sep = " ~ ")), 
#                 data = dat, listw=nb2listw(w), zero.policy = T, na.action = na.fail, 
#                 family="SAR")
# summary(sar)
# adj_r2(sar)
# modsel <- dredge(sar, beta = "none", trace = T, extra = c('r2', 'adj_r2'))
# zn_619 <- modsel[1]
# 
# metric <- "Cu"
# null <- spautolm( as.formula(paste(metric, 1, sep="~")), data = dat, 
#                   listw=nb2listw(w), zero.policy = T, na.action = na.fail, family="SAR")
# sar <- spautolm(as.formula(paste(metric, variables, sep = " ~ ")), 
#                 data = dat, listw=nb2listw(w), zero.policy = T, na.action = na.fail, 
#                 family="SAR")
# summary(sar)
# adj_r2(sar)
# modsel <- dredge(sar, beta = "none", trace = T, extra = c('r2', 'adj_r2'))
# cu_619 <- modsel[1]
# 
# metric <- "Fe"
# null <- spautolm( as.formula(paste(metric, 1, sep="~")), data = dat, 
#                   listw=nb2listw(w), zero.policy = T, na.action = na.fail, family="SAR")
# sar <- spautolm(as.formula(paste(metric, variables, sep = " ~ ")), 
#                 data = dat, listw=nb2listw(w), zero.policy = T, na.action = na.fail, 
#                 family="SAR")
# summary(sar)
# adj_r2(sar)
# modsel <- dredge(sar, beta = "none", trace = T, extra = c('r2', 'adj_r2'))
# fe_619 <- modsel[1]
# 
# metric <- "Al"
# null <- spautolm( as.formula(paste(metric, 1, sep="~")), data = dat, 
#                   listw=nb2listw(w), zero.policy = T, na.action = na.fail, family="SAR")
# sar <- spautolm(as.formula(paste(metric, variables, sep = " ~ ")), 
#                 data = dat, listw=nb2listw(w), zero.policy = T, na.action = na.fail, 
#                 family="SAR")
# summary(sar)
# adj_r2(sar)
# modsel <- dredge(sar, beta = "none", trace = T, extra = c('r2', 'adj_r2'))
# al_619 <- modsel[1]
# 
# # ---------------------------------------SAR 717--------------------------------------------
# dat <- std_717
# 
# 
# metric <- "N"
# null <- spautolm( as.formula(paste(metric, 1, sep="~")), data = dat, 
#                   listw=nb2listw(w), zero.policy = T, na.action = na.fail, family="SAR")
# sar <- spautolm(as.formula(paste(metric, variables, sep = " ~ ")), 
#                 data = dat, listw=nb2listw(w), zero.policy = T, na.action = na.fail, 
#                 family="SAR")
# summary(sar)
# adj_r2(sar)
# modsel <- dredge(sar, beta = "none", trace = T, extra = c('r2', 'adj_r2'))
# n_717 <- modsel[1]
# 
# metric <- "P"
# null <- spautolm( as.formula(paste(metric, 1, sep="~")), data = dat, 
#                   listw=nb2listw(w), zero.policy = T, na.action = na.fail, family="SAR")
# sar <- spautolm(as.formula(paste(metric, variables, sep = " ~ ")), 
#                 data = dat, listw=nb2listw(w), zero.policy = T, na.action = na.fail, 
#                 family="SAR")
# summary(sar)
# adj_r2(sar)
# modsel <- dredge(sar, beta = "none", trace = T, extra = c('r2', 'adj_r2'))
# p_717 <- modsel[1]
# 
# metric <- "K"
# null <- spautolm( as.formula(paste(metric, 1, sep="~")), data = dat, 
#                   listw=nb2listw(w), zero.policy = T, na.action = na.fail, family="SAR")
# sar <- spautolm(as.formula(paste(metric, variables, sep = " ~ ")), 
#                 data = dat, listw=nb2listw(w), zero.policy = T, na.action = na.fail, 
#                 family="SAR")
# summary(sar)
# adj_r2(sar)
# modsel <- dredge(sar, beta = "none", trace = T, extra = c('r2', 'adj_r2'))
# k_717 <- modsel[1]
# 
# metric <- "Ca"
# null <- spautolm( as.formula(paste(metric, 1, sep="~")), data = dat, 
#                   listw=nb2listw(w), zero.policy = T, na.action = na.fail, family="SAR")
# sar <- spautolm(as.formula(paste(metric, variables, sep = " ~ ")), 
#                 data = dat, listw=nb2listw(w), zero.policy = T, na.action = na.fail, 
#                 family="SAR")
# summary(sar)
# adj_r2(sar)
# modsel <- dredge(sar, beta = "none", trace = T, extra = c('r2', 'adj_r2'))
# ca_717 <- modsel[1]
# 
# metric <- "Mg"
# null <- spautolm( as.formula(paste(metric, 1, sep="~")), data = dat, 
#                   listw=nb2listw(w), zero.policy = T, na.action = na.fail, family="SAR")
# sar <- spautolm(as.formula(paste(metric, variables, sep = " ~ ")), 
#                 data = dat, listw=nb2listw(w), zero.policy = T, na.action = na.fail, 
#                 family="SAR")
# summary(sar)
# adj_r2(sar)
# modsel <- dredge(sar, beta = "none", trace = T, extra = c('r2', 'adj_r2'))
# mg_717 <- modsel[1]
# 
# metric <- "B"
# null <- spautolm( as.formula(paste(metric, 1, sep="~")), data = dat, 
#                   listw=nb2listw(w), zero.policy = T, na.action = na.fail, family="SAR")
# sar <- spautolm(as.formula(paste(metric, variables, sep = " ~ ")), 
#                 data = dat, listw=nb2listw(w), zero.policy = T, na.action = na.fail, 
#                 family="SAR")
# summary(sar)
# adj_r2(sar)
# modsel <- dredge(sar, beta = "none", trace = T, extra = c('r2', 'adj_r2'))
# b_717 <- modsel[1]
# 
# metric <- "Mn"
# null <- spautolm( as.formula(paste(metric, 1, sep="~")), data = dat, 
#                   listw=nb2listw(w), zero.policy = T, na.action = na.fail, family="SAR")
# sar <- spautolm(as.formula(paste(metric, variables, sep = " ~ ")), 
#                 data = dat, listw=nb2listw(w), zero.policy = T, na.action = na.fail, 
#                 family="SAR")
# summary(sar)
# adj_r2(sar)
# modsel <- dredge(sar, beta = "none", trace = T, extra = c('r2', 'adj_r2'))
# mn_717 <- modsel[1]
# 
# metric <- "S"
# null <- spautolm( as.formula(paste(metric, 1, sep="~")), data = dat, 
#                   listw=nb2listw(w), zero.policy = T, na.action = na.fail, family="SAR")
# sar <- spautolm(as.formula(paste(metric, variables, sep = " ~ ")), 
#                 data = dat, listw=nb2listw(w), zero.policy = T, na.action = na.fail, 
#                 family="SAR")
# summary(sar)
# adj_r2(sar)
# modsel <- dredge(sar, beta = "none", trace = T, extra = c('r2', 'adj_r2'))
# s_717 <- modsel[1]
# 
# metric <- "Zn"
# null <- spautolm( as.formula(paste(metric, 1, sep="~")), data = dat, 
#                   listw=nb2listw(w), zero.policy = T, na.action = na.fail, family="SAR")
# sar <- spautolm(as.formula(paste(metric, variables, sep = " ~ ")), 
#                 data = dat, listw=nb2listw(w), zero.policy = T, na.action = na.fail, 
#                 family="SAR")
# summary(sar)
# adj_r2(sar)
# modsel <- dredge(sar, beta = "none", trace = T, extra = c('r2', 'adj_r2'))
# zn_717 <- modsel[1]
# 
# metric <- "Cu"
# null <- spautolm( as.formula(paste(metric, 1, sep="~")), data = dat, 
#                   listw=nb2listw(w), zero.policy = T, na.action = na.fail, family="SAR")
# sar <- spautolm(as.formula(paste(metric, variables, sep = " ~ ")), 
#                 data = dat, listw=nb2listw(w), zero.policy = T, na.action = na.fail, 
#                 family="SAR")
# summary(sar)
# adj_r2(sar)
# modsel <- dredge(sar, beta = "none", trace = T, extra = c('r2', 'adj_r2'))
# cu_717 <- modsel[1]
# 
# metric <- "Fe"
# null <- spautolm( as.formula(paste(metric, 1, sep="~")), data = dat, 
#                   listw=nb2listw(w), zero.policy = T, na.action = na.fail, family="SAR")
# sar <- spautolm(as.formula(paste(metric, variables, sep = " ~ ")), 
#                 data = dat, listw=nb2listw(w), zero.policy = T, na.action = na.fail, 
#                 family="SAR")
# summary(sar)
# adj_r2(sar)
# modsel <- dredge(sar, beta = "none", trace = T, extra = c('r2', 'adj_r2'))
# fe_717 <- modsel[1]
# 
# metric <- "Al"
# null <- spautolm( as.formula(paste(metric, 1, sep="~")), data = dat, 
#                   listw=nb2listw(w), zero.policy = T, na.action = na.fail, family="SAR")
# sar <- spautolm(as.formula(paste(metric, variables, sep = " ~ ")), 
#                 data = dat, listw=nb2listw(w), zero.policy = T, na.action = na.fail, 
#                 family="SAR")
# summary(sar)
# adj_r2(sar)
# modsel <- dredge(sar, beta = "none", trace = T, extra = c('r2', 'adj_r2'))
# al_717 <- modsel[1]
# 
# 
# # ---------------------------------------SAR 717--------------------------------------------
# dat <- std_717
# 
# 
# metric <- "N"
# null <- spautolm( as.formula(paste(metric, 1, sep="~")), data = dat, 
#                   listw=nb2listw(w), zero.policy = T, na.action = na.fail, family="SAR")
# sar <- spautolm(as.formula(paste(metric, variables, sep = " ~ ")), 
#                 data = dat, listw=nb2listw(w), zero.policy = T, na.action = na.fail, 
#                 family="SAR")
# summary(sar)
# adj_r2(sar)
# modsel <- dredge(sar, beta = "none", trace = T, extra = c('r2', 'adj_r2'))
# n_717 <- modsel[1]
# 
# metric <- "P"
# null <- spautolm( as.formula(paste(metric, 1, sep="~")), data = dat, 
#                   listw=nb2listw(w), zero.policy = T, na.action = na.fail, family="SAR")
# sar <- spautolm(as.formula(paste(metric, variables, sep = " ~ ")), 
#                 data = dat, listw=nb2listw(w), zero.policy = T, na.action = na.fail, 
#                 family="SAR")
# summary(sar)
# adj_r2(sar)
# modsel <- dredge(sar, beta = "none", trace = T, extra = c('r2', 'adj_r2'))
# p_717 <- modsel[1]
# 
# metric <- "K"
# null <- spautolm( as.formula(paste(metric, 1, sep="~")), data = dat, 
#                   listw=nb2listw(w), zero.policy = T, na.action = na.fail, family="SAR")
# sar <- spautolm(as.formula(paste(metric, variables, sep = " ~ ")), 
#                 data = dat, listw=nb2listw(w), zero.policy = T, na.action = na.fail, 
#                 family="SAR")
# summary(sar)
# adj_r2(sar)
# modsel <- dredge(sar, beta = "none", trace = T, extra = c('r2', 'adj_r2'))
# k_717 <- modsel[1]
# 
# metric <- "Ca"
# null <- spautolm( as.formula(paste(metric, 1, sep="~")), data = dat, 
#                   listw=nb2listw(w), zero.policy = T, na.action = na.fail, family="SAR")
# sar <- spautolm(as.formula(paste(metric, variables, sep = " ~ ")), 
#                 data = dat, listw=nb2listw(w), zero.policy = T, na.action = na.fail, 
#                 family="SAR")
# summary(sar)
# adj_r2(sar)
# modsel <- dredge(sar, beta = "none", trace = T, extra = c('r2', 'adj_r2'))
# ca_717 <- modsel[1]
# 
# metric <- "Mg"
# null <- spautolm( as.formula(paste(metric, 1, sep="~")), data = dat, 
#                   listw=nb2listw(w), zero.policy = T, na.action = na.fail, family="SAR")
# sar <- spautolm(as.formula(paste(metric, variables, sep = " ~ ")), 
#                 data = dat, listw=nb2listw(w), zero.policy = T, na.action = na.fail, 
#                 family="SAR")
# summary(sar)
# adj_r2(sar)
# modsel <- dredge(sar, beta = "none", trace = T, extra = c('r2', 'adj_r2'))
# mg_717 <- modsel[1]
# 
# metric <- "B"
# null <- spautolm( as.formula(paste(metric, 1, sep="~")), data = dat, 
#                   listw=nb2listw(w), zero.policy = T, na.action = na.fail, family="SAR")
# sar <- spautolm(as.formula(paste(metric, variables, sep = " ~ ")), 
#                 data = dat, listw=nb2listw(w), zero.policy = T, na.action = na.fail, 
#                 family="SAR")
# summary(sar)
# adj_r2(sar)
# modsel <- dredge(sar, beta = "none", trace = T, extra = c('r2', 'adj_r2'))
# b_717 <- modsel[1]
# 
# metric <- "Mn"
# null <- spautolm( as.formula(paste(metric, 1, sep="~")), data = dat, 
#                   listw=nb2listw(w), zero.policy = T, na.action = na.fail, family="SAR")
# sar <- spautolm(as.formula(paste(metric, variables, sep = " ~ ")), 
#                 data = dat, listw=nb2listw(w), zero.policy = T, na.action = na.fail, 
#                 family="SAR")
# summary(sar)
# adj_r2(sar)
# modsel <- dredge(sar, beta = "none", trace = T, extra = c('r2', 'adj_r2'))
# mn_717 <- modsel[1]
# 
# metric <- "S"
# null <- spautolm( as.formula(paste(metric, 1, sep="~")), data = dat, 
#                   listw=nb2listw(w), zero.policy = T, na.action = na.fail, family="SAR")
# sar <- spautolm(as.formula(paste(metric, variables, sep = " ~ ")), 
#                 data = dat, listw=nb2listw(w), zero.policy = T, na.action = na.fail, 
#                 family="SAR")
# summary(sar)
# adj_r2(sar)
# modsel <- dredge(sar, beta = "none", trace = T, extra = c('r2', 'adj_r2'))
# s_717 <- modsel[1]
# 
# metric <- "Zn"
# null <- spautolm( as.formula(paste(metric, 1, sep="~")), data = dat, 
#                   listw=nb2listw(w), zero.policy = T, na.action = na.fail, family="SAR")
# sar <- spautolm(as.formula(paste(metric, variables, sep = " ~ ")), 
#                 data = dat, listw=nb2listw(w), zero.policy = T, na.action = na.fail, 
#                 family="SAR")
# summary(sar)
# adj_r2(sar)
# modsel <- dredge(sar, beta = "none", trace = T, extra = c('r2', 'adj_r2'))
# zn_717 <- modsel[1]
# 
# metric <- "Cu"
# null <- spautolm( as.formula(paste(metric, 1, sep="~")), data = dat, 
#                   listw=nb2listw(w), zero.policy = T, na.action = na.fail, family="SAR")
# sar <- spautolm(as.formula(paste(metric, variables, sep = " ~ ")), 
#                 data = dat, listw=nb2listw(w), zero.policy = T, na.action = na.fail, 
#                 family="SAR")
# summary(sar)
# adj_r2(sar)
# modsel <- dredge(sar, beta = "none", trace = T, extra = c('r2', 'adj_r2'))
# cu_717 <- modsel[1]
# 
# metric <- "Fe"
# null <- spautolm( as.formula(paste(metric, 1, sep="~")), data = dat, 
#                   listw=nb2listw(w), zero.policy = T, na.action = na.fail, family="SAR")
# sar <- spautolm(as.formula(paste(metric, variables, sep = " ~ ")), 
#                 data = dat, listw=nb2listw(w), zero.policy = T, na.action = na.fail, 
#                 family="SAR")
# summary(sar)
# adj_r2(sar)
# modsel <- dredge(sar, beta = "none", trace = T, extra = c('r2', 'adj_r2'))
# fe_717 <- modsel[1]
# 
# metric <- "Al"
# null <- spautolm( as.formula(paste(metric, 1, sep="~")), data = dat, 
#                   listw=nb2listw(w), zero.policy = T, na.action = na.fail, family="SAR")
# sar <- spautolm(as.formula(paste(metric, variables, sep = " ~ ")), 
#                 data = dat, listw=nb2listw(w), zero.policy = T, na.action = na.fail, 
#                 family="SAR")
# summary(sar)
# adj_r2(sar)
# modsel <- dredge(sar, beta = "none", trace = T, extra = c('r2', 'adj_r2'))
# al_717 <- modsel[1]
# 
# 
# # --------------------------------------------SAR structure results tables -------------------------------------
# sarresults_struct_619 <- rbind(n_619, p_619, k_619, ca_619, mg_619, b_619, mn_619, s_619, zn_619, cu_619, fe_619, al_619)
# sarresults_struct_619$depvar <- c("N", "P", "K", "Ca", "Mg", "B", "Mn", "S", "Zn", "Cu", "Fe", "Al")
# 
# sarresults_struct_717 <- rbind(n_717, p_717, k_717, ca_717, mg_717, b_717, mn_717, s_717, zn_717, cu_717, fe_717, al_717)
# sarresults_struct_717$depvar <- c("N", "P", "K", "Ca", "Mg", "B", "Mn", "S", "Zn", "Cu", "Fe", "Al")
# 
# sarresults_struct_717 <- rbind(n_717, p_717, k_717, ca_717, mg_717, b_717, mn_717, s_717, zn_717, cu_717, fe_717, al_717)
# sarresults_struct_717$depvar <- c("N", "P", "K", "Ca", "Mg", "B", "Mn", "S", "Zn", "Cu", "Fe", "Al")
# 
# write.csv(sarresults_struct_619, "Q:/My Drive/Research/Canopy_Morphology/Tobacco_Project/2019/sarresults_struct_619.csv")
# write.csv(sarresults_struct_717, "Q:/My Drive/Research/Canopy_Morphology/Tobacco_Project/2019/sarresults_struct_717.csv")
# write.csv(sarresults_struct_717, "Q:/My Drive/Research/Canopy_Morphology/Tobacco_Project/2019/sarresults_struct_717.csv")



# #--------------------------------------------Try NN regression---------------------------------------------------------
# library(neuralnet)
# set.seed(42)
# train_data_indices <- rep(FALSE, nrow(agg_619))
# train_data_indices[sample(1:nrow(agg_619), round(0.7 * nrow(agg_619)))] <- TRUE # randomly select 70% of the data for training
# nn <- neuralnet(N ~ median++iqr+crrmedian+crriqr+skew+kurt+moran+watersum,data=agg_619[train_data_indices, ],hidden=c(3,2),linear.output=T)
# pred_nn <- predict(nn, newdata=agg_619[-train_data_indices, ])
# 
# pred_nn_res <- pred_nn*(max(agg_619$median)-min(agg_619$median))+min(agg_619$median)
# test_res <- (agg_619[-train_data_indices, ]$median)*(max(agg_619$median)-min(agg_619$median))+min(agg_619$median)
# MSE_nn <- sum((test_res - pred_nn_res)^2)/nrow(agg_619[-train_data_indices, ])
# 
# # compare to SAR
# 
# w <- knn2nb(knearneigh(coordinates(plots_619[train_data_indices, ]), k=8))
# 
# variables <- "`N`+`Ca`+`S`+`Cu`+`B`+`Al`+`soil_GoA`+`soil_Ra`+`soil_MaB`"
# dat <- agg_619[train_data_indices, ]
# 
# metric <- "median"
# null <- errorsarlm( as.formula(paste(metric, 1, sep="~")), data = dat, 
#                   listw=nb2listw(w), zero.policy = T, na.action = na.fail)
# sar <- errorsarlm(as.formula(paste(metric, variables, sep = " ~ ")), 
#                 data = dat, listw=nb2listw(w), zero.policy = T, na.action = na.fail, 
#                )
# summary(sar)
# 
# predict_sar <- predict.sarlm(sar,newdata=agg_619[-train_data_indices, ] )
# MSE_sar <- sum((predict_sar - agg_619[-train_data_indices, ]$median)^2)/nrow(agg_619[-train_data_indices, ])
# 
# 
# # ---------------------------------------- Try RF regression --------------------------------------------------
# 
# #all nutrients: `N` + `P` + `K`+`Mg` +`Ca` + `S` + `Zn`+`Mn`+`Cu`+`Fe`+`B`+`Al`+`Na`
# 
# #--------------------------------RF 619 ---------------------------------
# 
# set.seed(42)
# train_data_indices <- rep(FALSE, nrow(agg_619))
# train_data_indices[sample(1:nrow(agg_619), round(0.7 * nrow(agg_619)))] <- TRUE # randomly select 70% of the data for training
# rf_regression<- randomForest(N ~ median+iqr+crrmedian+crriqr+skew+kurt+moran+watersum, data = agg_619[train_data_indices, ], importance=T)
# rf_regression
# varImpPlot(rf_regression)
# pred_struc <- predict(rf_regression, agg_619[!train_data_indices,]) # predict the rings
# #table(pred_struc, agg_619$median[!train_data_indices])
# plot(agg_619$skew[!train_data_indices], pred_struc, xlab="Observed", ylab="Predicted")
# 
# #--------------------------------RF 717 ---------------------------------
# 
# set.seed(42)
# train_data_indices <- rep(FALSE, nrow(agg_717))
# train_data_indices[sample(1:nrow(agg_717), round(0.7 * nrow(agg_717)))] <- TRUE # randomly select 70% of the data for training
# rf_regression<- randomForest(skew ~ `N` + `P` + `K`+`Mg` +`Ca` + `Zn`+`Cu`+`Fe`+`B`+`Al`+`Na`, data = agg_717[train_data_indices, ], importance=T)
# rf_regression
# varImpPlot(rf_regression)
# pred_struc <- predict(rf_regression, agg_717[!train_data_indices,]) # predict the rings
# #table(pred_struc, agg_717$median[!train_data_indices])
# plot(agg_717$skew[!train_data_indices], pred_struc, xlab="Observed", ylab="Predicted")
# 
# #--------------------------------RF 717 ---------------------------------
# 
# set.seed(42)
# train_data_indices <- rep(FALSE, nrow(agg_717))
# train_data_indices[sample(1:nrow(agg_717), round(0.7 * nrow(agg_717)))] <- TRUE # randomly select 70% of the data for training
# rf_regression<- randomForest(median ~ `watersum`+`N` + `P` + `K`+`Mg` +`Ca` + `Zn`+`Cu`+`Fe`+`B`+`Al`+`Na`, data = agg_717[train_data_indices, ], importance=T)
# rf_regression
# varImpPlot(rf_regression)
# pred_struc <- predict(rf_regression, agg_717[!train_data_indices,]) # predict the rings
# #table(pred_struc, agg_717$median[!train_data_indices])
# plot(agg_717$median[!train_data_indices], pred_struc, xlab="Observed", ylab="Predicted")
# 
# 
# # ---------------------------------------Try multinominal logisitic regression with nnet package --------------------------------
# 
# # 619
# 
# library(nnet)
# plots_619_df$Trt_agg <- relevel(plots_619_df$Trt_agg, ref = "100% Control")
# mltnom_619 <- multinom(Trt_agg ~ crrmedian+median+kurt+rumple+moran+dip+crriqr, data = plots_619_df, na.action = na.fail)
# summary(mltnom_619)
# 
# modsel_619 <- dredge(mltnom_619)
# modsel_619
# 
# topmod_619 <- multinom(Trt_agg ~ crrmedian+median+kurt+rumple+moran, data = plots_619_df, na.action = na.fail)
# summary(topmod_619)
# 
# library(DescTools)
# PseudoR2(topmod_619, which = "all")
# z <- summary(topmod_619)$coefficients/summary(topmod_619)$standard.errors
# z
# p <- (1 - pnorm(abs(z), 0, 1)) * 2
# p
# 
# coefs_619 <- coef(topmod_619)
# exp(coefs_619) #transformed to odds ratio
# (exp(coefs_619)-1)*100 # percent change in the odds for a one unit increase in the independent variable
# 
# # Try binary logistic regression with 100% control as baseline
# # def_100 <- as.data.frame(dt_619[Trt_agg == "100% Control" | Trt_agg == "deficiency",])
# # def_100 <- droplevels(def_100)
# # tox_100 <- as.data.frame(dt_619[Trt_agg == "100% Control" | Trt_agg == "toxicity",])
# # tox_100 <- droplevels(tox_100)
# # control_0_100 <- as.data.frame(dt_619[Trt_agg == "100% Control" | Trt_agg == "0% Control",])
# # control_0_100 <- droplevels(control_0_100)
# # 
# # binary_619 <- glm(Trt_agg ~ crrmedian+median+kurt+rumple+moran+crriqr+skew+iqr, data = control_0_100, family = binomial, na.action = na.fail)
# # modsel_619_binary <- dredge(binary_619)
# # topmod_619_binary <- glm(Trt_agg ~ crrmedian+iqr+moran+rumple+skew, data = tox_100, family = binomial)
# # summary(topmod_619_binary)
# # 
# # library(pscl)
# # 
# # pR2(topmod_619_binary)
# 
# # Try multinomial logisitc regression with mlogit package
# # library(foreign)
# # library(mlogit)
# # 
# # mlogit_data_619 <- mlogit.data(plots_619_df, choice = "Trt_agg", shape = "wide")
# # mlogit_619 <- mlogit(Trt_agg~crrmedian+kurt+median+skew_ln+rumple+moran, data = mlogit_data_619,
# #                      method = "nr", print.level = 0)
# # summary(mlogit_619)
# # 
# # # Weights matrix
# # w <- knn2nb(knearneigh(coordinates(plots), k=8))
# 
# 
# # 717
# plots_717_df$Trt_agg <- relevel(plots_717_df$Trt_agg, ref = "100% Control")
# mltnom_717 <- multinom(Trt_agg ~ crrmedian+median+kurt+rumple+moran+dip+crriqr, data = plots_717_df, na.action = na.fail)
# summary(mltnom_717)
# 
# modsel_717 <- dredge(mltnom_717)
# modsel_717
# 
# topmod_717 <- multinom(Trt_agg ~ crrmedian+kurt+median+moran, data = plots_717_df, na.action = na.fail)
# summary(topmod_717)
# PseudoR2(topmod_717, which = "all")
# z <- summary(topmod_717)$coefficients/summary(topmod_717)$standard.errors
# z
# p <- (1 - pnorm(abs(z), 0, 1)) * 2
# p
# 
# coefs_717 <- coef(topmod_717)
# exp(coefs_717) #transformed to odds ratio
# (exp(coefs_717)-1)*100 # percent change in the odds for a one unit increase in the independent variable
# 
# # 717
# plots_717_df$Trt_agg <- relevel(plots_717_df$Trt_agg, ref = "100% Control")
# mltnom_717 <- multinom(Trt_agg ~ crrmedian+median+kurt+rumple+moran+dip+crriqr, data = plots_717_df, na.action = na.fail)
# summary(mltnom_717)
# 
# modsel_717 <- dredge(mltnom_717)
# modsel_717
# 
# topmod_717 <- multinom(Trt_agg ~ median+rumple+moran, data = plots_717_df, na.action = na.fail)
# summary(topmod_717)
# PseudoR2(topmod_717, which = "all")
# z <- summary(topmod_717)$coefficients/summary(topmod_717)$standard.errors
# z
# p <- (1 - pnorm(abs(z), 0, 1)) * 2
# p
# 
# coefs_717 <- coef(topmod_717)
# exp(coefs_717) #transformed to odds ratio
# (exp(coefs_717)-1)*100 # percent change in the odds for a one unit increase in the independent variable
# 
# 
# # Try combining tox and def to see if model improves
# plots_717_df$Trt_agg_v2 <- as.character(plots_619_df$Trt_agg)
# plots_717_df[Trt_agg == "deficiency"]$Trt_agg_v2 <- "stressed"
# plots_717_df[Trt_agg == "toxicity"]$Trt_agg_v2 <- "stressed"
# plots_717_df$Trt_agg_v2 <- as.factor(plots_717_df$Trt_agg_v2) 
# 
# plots_717_df$Trt_agg <- relevel(plots_717_df$Trt_agg, ref = "100% Control")
# 
# mltnom_717_v2 <- multinom(Trt_agg_v2 ~ crrmedian+median+kurt+rumple+moran+dip+crriqr, data = plots_717_df, na.action = na.fail)
# summary(mltnom_717_v2)
# 
# modsel_717_v2 <- dredge(mltnom_717_v2)
# modsel_717_v2
# 
# topmod_717_v2 <- multinom(Trt_agg_v2 ~ median+rumple+moran+crrmedian, data = plots_717_df, na.action = na.fail)
# summary(topmod_717_v2)
# PseudoR2(topmod_717_v2, which = "all")
# z <- summary(topmod_717)$coefficients/summary(topmod_717)$standard.errors
# z
# p <- (1 - pnorm(abs(z), 0, 1)) * 2
# p
# 
# 
# # ----------------------------------------------Result tables------------------------------------------------------
# 
# vars <- topmod_619$coefnames[-1]
# sd_all <- numeric(length(vars))
# mean_all <- numeric(length(vars))
# count <- 1
# for (x in vars){
#   values <- as.data.frame(plots_619_df)[x]
#   sd_all[count] <- sd(values[,1])
#   mean_all[count] <- mean(values[,1])
#   count <- count+1
# }
# 
# vars_619 <- as.data.frame(vars)
# vars_619$mean <- mean_all
# vars_619$sd <- sd_all
# 
# #write.csv(vars_619, "Q:/My Drive/Research/Canopy_Morphology/Tobacco_Project/canopy_analysis/tobacco_results/vars_619.csv")
# 
# levels <- topmod_619$lev[-1]
# #df_all <- data.frame(var=character(), coefs=list(), SE=double(), Wald=double(), p=double(), exp=double())
# coefs <- numeric(length(levels))
# std_coefs <- numeric(length(levels))
# SE <- numeric(length(levels))
# Wald <- numeric(length(levels))
# p <- numeric(length(levels))
# exp <- numeric(length(levels))
# count <- 1
# 
# for (level in levels){
#   coefs[count] <- list(coef(topmod_619)[count,])
#   SE[count] <- list(summary(topmod_619)$standard.errors[count,])
#   z <- (summary(topmod_619)$coefficients/summary(topmod_619)$standard.errors)
#   Wald[count] <- list(z[count,])
#   p[count] <- list(((1 - pnorm(abs(z), 0, 1)) * 2)[count,])
#   std_coefs[count] <- list((sqrt(3)/pi) * (coefs[[count]][-1] * sd_all))
#   exp[count] <- list(exp(std_coefs[[count]]))
#   count <- count+1
# }
# 
# vars <- topmod_619$coefnames
# 
# results_0_619 <- as.data.frame(vars)
# results_0_619$level <- "0% Control"
# results_0_619$coefs <- coefs[[1]]
# results_0_619$std_coefs <- c(0,std_coefs[[1]])
# results_0_619$exp <- c(0,exp[[1]])
# results_0_619$SE <- SE[[1]]
# results_0_619$p <- p[[1]]
# results_0_619$wald <- Wald[[1]]
# 
# results_tox_619 <- as.data.frame(vars)
# results_tox_619$level <- "toxicity"
# results_tox_619$coefs <- coefs[[3]]
# results_tox_619$std_coefs <- c(0,std_coefs[[3]])
# results_tox_619$exp <- c(0,exp[[3]])
# results_tox_619$SE <- SE[[3]]
# results_tox_619$p <- p[[3]]
# results_tox_619$wald <- Wald[[3]]
# 
# results_def_619 <- as.data.frame(vars)
# results_def_619$level <- "deficiency"
# results_def_619$coefs <- coefs[[2]]
# results_def_619$std_coefs <- c(0,std_coefs[[2]])
# results_def_619$exp <- c(0,exp[[2]])
# results_def_619$SE <- SE[[2]]
# results_def_619$p <- p[[2]]
# results_def_619$wald <- Wald[[2]]
# 
# 
# 
# results_619 <- rbind(results_0_619, results_def_619, results_tox_619)
# 
# 
# # mlogit_results <- function(x){
# #   levels <- x$lev[-1]
# #   coefs <- numeric(length(levels))
# #   SE <- numeric(length(levels))
# #   Wald <- numeric(length(levels))
# #   p <- numeric(length(levels))
# #   exp <- numeric(length(levels))
# #   count <- 1
# #   for (level in levels){
# #     coefs[count] <- list(coef(x)[count,])
# #     SE[count] <- list(summary(x)$standard.errors[count,])
# #     Wald[count] <- list(summary(x)$coefficients/summary(topmod_619)$standard.errors[count,])
# #     p[count] <- list(((1 - pnorm(abs(z), 0, 1)) * 2)[count,])
# #     exp[count] <- list(exp(coefs[[count]]))
# #     count <- count+1
# #   }
# #   levels
# #   coefs
# #   SE
# #   Wald
# #   p
# #   exp
# # }
# # 
# # mlogit_results_619 <- mlogit_results(topmod_619)
# 
# 
# #write.csv(results_619, "Q:/My Drive/Research/Canopy_Morphology/Tobacco_Project/canopy_analysis/tobacco_results/results_619.csv")
# 
# 
# vars <- topmod_717$coefnames[-1]
# sd_all <- numeric(length(vars))
# mean_all <- numeric(length(vars))
# count <- 1
# for (x in vars){
#   values <- as.data.frame(plots_717_df)[x]
#   sd_all[count] <- sd(values[,1])
#   mean_all[count] <- mean(values[,1])
#   count <- count+1
# }
# 
# vars_717 <- as.data.frame(vars)
# vars_717$mean <- mean_all
# vars_717$sd <- sd_all
# 
# #write.csv(vars_717, "Q:/My Drive/Research/Canopy_Morphology/Tobacco_Project/canopy_analysis/tobacco_results/vars_717.csv")
# 
# levels <- topmod_717$lev[-1]
# #df_all <- data.frame(var=character(), coefs=list(), SE=double(), Wald=double(), p=double(), exp=double())
# coefs <- numeric(length(levels))
# std_coefs <- numeric(length(levels))
# SE <- numeric(length(levels))
# Wald <- numeric(length(levels))
# p <- numeric(length(levels))
# exp <- numeric(length(levels))
# count <- 1
# 
# for (level in levels){
#   coefs[count] <- list(coef(topmod_717)[count,])
#   SE[count] <- list(summary(topmod_717)$standard.errors[count,])
#   z <- (summary(topmod_717)$coefficients/summary(topmod_717)$standard.errors)
#   Wald[count] <- list(z[count,])
#   p[count] <- list(((1 - pnorm(abs(z), 0, 1)) * 2)[count,])
#   std_coefs[count] <- list((sqrt(3)/pi) * (coefs[[count]][-1] * sd_all))
#   exp[count] <- list(exp(std_coefs[[count]]))
#   count <- count+1
# }
# 
# vars <- topmod_717$coefnames
# 
# results_0_717 <- as.data.frame(vars)
# results_0_717$level <- "0% Control"
# results_0_717$coefs <- coefs[[1]]
# results_0_717$std_coefs <- c(0,std_coefs[[1]])
# results_0_717$exp <- c(0,exp[[1]])
# results_0_717$SE <- SE[[1]]
# results_0_717$p <- p[[1]]
# results_0_717$wald <- Wald[[1]]
# 
# results_tox_717 <- as.data.frame(vars)
# results_tox_717$level <- "toxicity"
# results_tox_717$coefs <- coefs[[3]]
# results_tox_717$std_coefs <- c(0,std_coefs[[3]])
# results_tox_717$exp <- c(0,exp[[3]])
# results_tox_717$SE <- SE[[3]]
# results_tox_717$p <- p[[3]]
# results_tox_717$wald <- Wald[[3]]
# 
# results_def_717 <- as.data.frame(vars)
# results_def_717$level <- "deficiency"
# results_def_717$coefs <- coefs[[2]]
# results_def_717$std_coefs <- c(0,std_coefs[[2]])
# results_def_717$exp <- c(0,exp[[2]])
# results_def_717$SE <- SE[[2]]
# results_def_717$p <- p[[2]]
# results_def_717$wald <- Wald[[2]]
# 
# 
# results_717 <- rbind(results_0_717, results_def_717, results_tox_717)
# 
# #write.csv(results_717, "Q:/My Drive/Research/Canopy_Morphology/Tobacco_Project/canopy_analysis/tobacco_results/results_717.csv")
# 
# 
# vars <- topmod_717$coefnames[-1]
# sd_all <- numeric(length(vars))
# mean_all <- numeric(length(vars))
# count <- 1
# for (x in vars){
#   values <- as.data.frame(plots_717_df)[x]
#   sd_all[count] <- sd(values[,1])
#   mean_all[count] <- mean(values[,1])
#   count <- count+1
# }
# 
# vars_717 <- as.data.frame(vars)
# vars_717$mean <- mean_all
# vars_717$sd <- sd_all
# 
# #write.csv(vars_717, "Q:/My Drive/Research/Canopy_Morphology/Tobacco_Project/canopy_analysis/tobacco_results/vars_717.csv")
# 
# levels <- topmod_717$lev[-1]
# #df_all <- data.frame(var=character(), coefs=list(), SE=double(), Wald=double(), p=double(), exp=double())
# coefs <- numeric(length(levels))
# std_coefs <- numeric(length(levels))
# SE <- numeric(length(levels))
# Wald <- numeric(length(levels))
# p <- numeric(length(levels))
# exp <- numeric(length(levels))
# count <- 1
# 
# for (level in levels){
#   coefs[count] <- list(coef(topmod_717)[count,])
#   SE[count] <- list(summary(topmod_717)$standard.errors[count,])
#   z <- (summary(topmod_717)$coefficients/summary(topmod_717)$standard.errors)
#   Wald[count] <- list(z[count,])
#   p[count] <- list(((1 - pnorm(abs(z), 0, 1)) * 2)[count,])
#   std_coefs[count] <- list((sqrt(3)/pi) * (coefs[[count]][-1] * sd_all))
#   exp[count] <- list(exp(std_coefs[[count]]))
#   count <- count+1
# }
# 
# vars <- topmod_717$coefnames
# 
# results_0_717 <- as.data.frame(vars)
# results_0_717$level <- "0% Control"
# results_0_717$coefs <- coefs[[1]]
# results_0_717$std_coefs <- c(0,std_coefs[[1]])
# results_0_717$exp <- c(0,exp[[1]])
# results_0_717$SE <- SE[[1]]
# results_0_717$p <- p[[1]]
# results_0_717$wald <- Wald[[1]]
# 
# results_tox_717 <- as.data.frame(vars)
# results_tox_717$level <- "toxicity"
# results_tox_717$coefs <- coefs[[3]]
# results_tox_717$std_coefs <- c(0,std_coefs[[3]])
# results_tox_717$exp <- c(0,exp[[3]])
# results_tox_717$SE <- SE[[3]]
# results_tox_717$p <- p[[3]]
# results_tox_717$wald <- Wald[[3]]
# 
# results_def_717 <- as.data.frame(vars)
# results_def_717$level <- "deficiency"
# results_def_717$coefs <- coefs[[2]]
# results_def_717$std_coefs <- c(0,std_coefs[[2]])
# results_def_717$exp <- c(0,exp[[2]])
# results_def_717$SE <- SE[[2]]
# results_def_717$p <- p[[2]]
# results_def_717$wald <- Wald[[2]]
# 
# results_717 <- rbind(results_0_717, results_def_717, results_tox_717)
# 
# #write.csv(results_717, "Q:/My Drive/Research/Canopy_Morphology/Tobacco_Project/canopy_analysis/tobacco_results/results_717.csv")
# 
# 
# # ------------------------------------- Random Forest Classification ------------------------------------------------------
# library(randomForest)
# set.seed(42)
# train_data_indices <- rep(FALSE, nrow(norm_1cm))
# train_data_indices[sample(1:nrow(norm_1cm), round(0.7 * nrow(norm_1cm)))] <- TRUE # randomly select 80% of the data for training
# rf_regression_1cm<- randomForest(BuAc_ln ~ median+skew+PlotCRR_ln+crrmedian_ln+crriqr+Trt+AUDPC_ln+geary+iqr_ln, data=norm_1cm[train_data_indices, ], importance=T)
# rf_regression_1cm
# varImpPlot(rf_regression_1cm)
# pred_trt <- predict(rf_regression_1cm, norm_1cm[!train_data_indices,]) # predict the rings
# table(pred_trt, norm_1cm$BuAc_ln[!train_data_indices])
# plot(norm_1cm$BuAc_ln[!train_data_indices], pred_trt, xlab="Observed", ylab="Predicted")
# #abline(a=0, b=1, lty=2, col=2)
# 
# set.seed(42)
# train_data_indices <- rep(FALSE, nrow(norm_5cm))
# train_data_indices[sample(1:nrow(norm_5cm), round(0.7 * nrow(norm_5cm)))] <- TRUE # randomly select 80% of the data for training
# rf_regression_5cm<- randomForest(BuAc_ln ~ median+skew+PlotCRR_ln+crrmedian_ln+crriqr+Trt+AUDPC_ln+geary+iqr_ln, data=norm_5cm[train_data_indices, ], importance=T)
# rf_regression_5cm
# varImpPlot(rf_regression_5cm)
# pred_trt <- predict(rf_regression_5cm, norm_5cm[!train_data_indices,]) # predict the rings
# table(pred_trt, norm_5cm$BuAc_ln[!train_data_indices])
# plot(norm_5cm$BuAc_ln[!train_data_indices], pred_trt, xlab="Observed", ylab="Predicted")
# #abline(a=0, b=1, lty=2, col=2)
# 
# set.seed(42)
# train_data_indices <- rep(FALSE, nrow(norm_10cm))
# train_data_indices[sample(1:nrow(norm_10cm), round(0.7 * nrow(norm_10cm)))] <- TRUE # randomly select 80% of the data for training
# rf_regression_10cm<- randomForest(BuAc_ln ~ median+skew+PlotCRR_ln+crrmedian_ln+crriqr+Trt+AUDPC_ln+geary+iqr_ln, data=norm_10cm[train_data_indices, ], importance=T)
# rf_regression_10cm
# varImpPlot(rf_regression_10cm)
# pred_trt <- predict(rf_regression_10cm, norm_10cm[!train_data_indices,]) # predict the rings
# table(pred_trt, norm_10cm$BuAc_ln[!train_data_indices])
# plot(norm_10cm$BuAc_ln[!train_data_indices], pred_trt, xlab="Observed", ylab="Predicted")
# #abline(a=0, b=1, lty=2, col=2)
# 
# set.seed(42)
# train_data_indices <- rep(FALSE, nrow(norm_20cm))
# train_data_indices[sample(1:nrow(norm_20cm), round(0.7 * nrow(norm_20cm)))] <- TRUE # randomly select 80% of the data for training
# rf_regression_20cm<- randomForest(BuAc_ln ~ median+skew+PlotCRR_ln+crrmedian_ln+crriqr+Trt+AUDPC_ln+geary+iqr_ln, data=norm_20cm[train_data_indices, ], importance=T)
# rf_regression_20cm
# varImpPlot(rf_regression_20cm)
# pred_trt <- predict(rf_regression_20cm, norm_20cm[!train_data_indices,]) # predict the rings
# table(pred_trt, norm_20cm$BuAc_ln[!train_data_indices])
# plot(norm_20cm$BuAc_ln[!train_data_indices], pred_trt, xlab="Observed", ylab="Predicted")
# #abline(a=0, b=1, lty=2, col=2)

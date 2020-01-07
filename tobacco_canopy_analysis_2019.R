library(raster)
library(rgdal)
library(data.table)
library(plyr)
library(lattice)
library(gridExtra)
library(spdep)
library(velox)
library(moments)
library(rgeos)
library(gstat)
library(usdm)
library(viridis)
library(randomForest)
library(diptest)
library(ggpubr)
library(foreign)
library(mlogit)
library(car)
library(nnet)
library(MuMIn)
library(data.table)
library(PerformanceAnalytics)
library(rcompanion)
library(olsrr)
library(spatialreg)
library(fastDummies)

#Load raster data for all dates
#setwd("/media/Kellyn/F20E17B40E177139/kpmontgo@ncsu.edu/LkWheeler_Sorghum/LkWheeler_Fusarium_Sorghum")
setwd("Q:/My Drive/Research/Canopy_Morphology/Tobacco_Project/2019")

tobacco_area <- readOGR("Q:/My Drive/Research/Canopy_Morphology/Tobacco_Project/2019/layers/boundary", "wilson19_boundary", stringsAsFactors = F)
tobacco_area <- spTransform(tobacco_area, CRS("+proj=lcc +lat_1=36.16666666666666 +lat_2=34.33333333333334 +lat_0=33.75 +lon_0=-79 +x_0=609601.22 +y_0=0 +datum=NAD83 +units=m +no_defs +ellps=GRS80 +towgs84=0,0,0"))

csm_619 <- crop(raster('Q:/My Drive/Research/Canopy_Morphology/Tobacco_Project/2019/csm/wilson19_619_csm.tif'), tobacco_area) #convert to cm
csm_703 <- crop(raster('Q:/My Drive/Research/Canopy_Morphology/Tobacco_Project/2019/csm/wilson19_703_csm.tif'), tobacco_area)
csm_717 <- crop(raster('Q:/My Drive/Research/Canopy_Morphology/Tobacco_Project/2019/csm/wilson19_717_csm.tif'), tobacco_area)
dem <- crop(raster('Q:/My Drive/Research/Canopy_Morphology/Tobacco_Project/2019/lidar/wilson19_dem.tif'), tobacco_area)

csm_619_velox <- velox(csm_619)
csm_703_velox <- velox(csm_703)
csm_717_velox <- velox(csm_717)

plots <- readOGR("Q:/My Drive/Research/Canopy_Morphology/Tobacco_Project/2019/layers/plots", "wilson19_plots", stringsAsFactors = F)
#plots$fertilizer <- as.factor(plots$fertilizer)
#plots$treatment <- as.factor(plots$treatment)
#plots$group <- as.factor(plots$group)
plots <- spTransform(plots, CRS("+proj=lcc +lat_1=36.16666666666666 +lat_2=34.33333333333334 +lat_0=33.75 +lon_0=-79 +x_0=609601.22 +y_0=0 +datum=NAD83 +units=m +no_defs +ellps=GRS80 +towgs84=0,0,0"))
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
plots_703 <- merge(plots_703, nutrient[Date==701], by.x="join_id", by.y='Join_ID')
plots_717 <- plots
plots_717 <- merge(plots_717, nutrient[Date==718], by.x="join_id", by.y='Join_ID')

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
crrRast_703 <- focal(csm_703, w=matrix(1/81, nc=9, nr=9), crr)
crrRast_717 <- focal(csm_717, w=matrix(1/81, nc=9, nr=9), crr)
crrRast_list <- list(crrRast_619,crrRast_703,crrRast_717)

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
  crriqr <- sapply(extract, IQR, na.rm=T)
  crrDF <- data.frame(plots=plots$plotid, crrmedian= crrmedian, crrmean=crrmean, crriqr = crriqr)
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

# # Check for multimodality
# dip.test(plot_extract_619$`502-E`)
# histogram(plot_extract_619$`502-E`)
# 
# multimodality_test <- function(x){
#   plotNums <- plots$plotid
#   dip_all <- numeric(length(plotNums))
#   names(dip_all) <- "dip"
#   count <- 1
#   for (plot in x){
#     test <- dip.test(plot)
#     p <- test$p.value
#     dip_all[count] <- p
#     count <- count+1
#   }
#   dip_df <- data.frame(plots=plots$plotid, dip= dip_all)
#   dip_df
# }
# 
# dip_df_list <- lapply(extract_list, multimodality_test)
# 
# #merge all metrics into spatial polygons
# plots_619 <- merge(plots_619, dip_df_list[[1]], by.x="plotid", by.y="plots")
# plots_703 <- merge(plots_703, dip_df_list[[2]], by.x="plotid", by.y="plots")
# plots_717 <- merge(plots_717, dip_df_list[[3]], by.x="plotid", by.y="plots")


# Check for normality
#shapiro.test(plot_extract_619$`502-E`)

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



min(plot_extract_619$`501-E`)

csm_list <- list(csm_619, csm_703, csm_717)
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

#Using approx 30 cm x 30 cm neighborhood matrix
autocor_619 <- autocor_metrics(csm_619, matrix(c(1,1,1,1,1,1,1,1,1,1,1,1,0,1,1,1,1,1,1,1,1,1,1,1,1), nc=5, nr=5))
autocor_703 <- autocor_metrics(csm_703, matrix(c(1,1,1,1,1,1,1,1,1,1,1,1,0,1,1,1,1,1,1,1,1,1,1,1,1), nc=5, nr=5))
autocor_717 <- autocor_metrics(csm_717, matrix(c(1,1,1,1,1,1,1,1,1,1,1,1,0,1,1,1,1,1,1,1,1,1,1,1,1), nc=5, nr=5))

plots_619 <- merge(plots_619, autocor_619, by.x="plotid", by.y="plots")
plots_703 <- merge(plots_703, autocor_703, by.x="plotid", by.y="plots")
plots_717 <- merge(plots_717, autocor_717, by.x="plotid", by.y="plots")


# Create rasters representing local Moran's I for visualization
# autocor_rasters <- function(x, w){
#   plotNums <- plots$plotid
#   moran_local_all <- as.list(x)
#   count <- 1
#     
#   for (num in plotNums){
#     tmp <- plots@data
#     position <- match(num, tmp$Id)
#     extract <- crop(x, plots[position,])
#     moran_local <- MoranLocal(extract, w)
#     moran_local_all <- append(moran_local_all, moran_local)
#     count <- count+1
#   }
#   
#   moran_local_all <- moran_local_all[-1]
#   names(moran_local_all) <- plots$plotid
#   moran_local_all
# }
# 
# local_moran_plots_619 <- autocor_rasters(csm_619, matrix(c(1,1,1,1,1,1,1,1,1,1,1,1,0,1,1,1,1,1,1,1,1,1,1,1,1), nc=5, nr=5))
# names(local_moran_plots_619) <- NULL
# local_moran_plots_619$filename <- 'moran_619.tif'
# local_moran_plots_619$overwrite <- TRUE
# local_moran_plots_619$fun <- sum
# merged_local_moran_plots_619 <- do.call(raster::merge, local_moran_plots_619)
# #writeRaster(merged_local_moran_plots_619, "Q:/My Drive/Research/Canopy_Morphology/Tobacco_Project/canopy_analysis/layers/moran/local_moran_plots_619.tif", format="GTiff", overwrite = T)
# 
# local_moran_plots_703 <- autocor_rasters(csm_703, matrix(c(1,1,1,1,1,1,1,1,1,1,1,1,0,1,1,1,1,1,1,1,1,1,1,1,1), nc=5, nr=5))
# names(local_moran_plots_703) <- NULL
# local_moran_plots_703$filename <- 'moran_703.tif'
# local_moran_plots_703$overwrite <- TRUE
# local_moran_plots_703$fun <- sum
# merged_local_moran_plots_703 <- do.call(raster::merge, local_moran_plots_703)
# #writeRaster(merged_local_moran_plots_703, "Q:/My Drive/Research/Canopy_Morphology/Tobacco_Project/canopy_analysis/layers/moran/local_moran_plots_703.tif", format="GTiff", overwrite = T)
# 
# local_moran_plots_717 <- autocor_rasters(csm_717, matrix(c(1,1,1,1,1,1,1,1,1,1,1,1,0,1,1,1,1,1,1,1,1,1,1,1,1), nc=5, nr=5))
# names(local_moran_plots_717) <- NULL
# local_moran_plots_717$filename <- 'moran_717.tif'
# local_moran_plots_717$overwrite <- TRUE
# local_moran_plots_717$fun <- sum
# merged_local_moran_plots_717 <- do.call(raster::merge, local_moran_plots_717)
# #writeRaster(merged_local_moran_plots_717, "Q:/My Drive/Research/Canopy_Morphology/Tobacco_Project/canopy_analysis/layers/moran/local_moran_plots_717.tif", format="GTiff", overwrite = T)

# ------------------------------------Plot variables ------------------------------------------
library(ggpubr)

imp_vars <- c("mean", "median", "skew", "iqr", "kurt", "sd", "crrsd", "crrmean", "moran", "Trt", "AUDPC", "watersum", "PlotCRR", "rumple") 

pairs(median~`N`+`P`+`K`+`Mg`+`Ca`+`S`+`Zn`+`Mn`+`Cu`+`Fe`+`B`+`Al`+`Na`,data=plots_619, lower.panel = NULL, 
      main="Simple Scatterplot Matrix")

par(mfrow=c(2,4), mar=c(5, 6, 2, 3))

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

#-------------------------------------Normalize variables---------------------------------------

library(rcompanion)

shapiro.test(plots_619$crrmedian)
# shapiro.test(plots_619$crriqr)
# shapiro.test(plots_619$dip)
# shapiro.test(plots_619$median)
shapiro.test(plots_619$iqr) #normal
# shapiro.test(plots_619$skew)
# shapiro.test(plots_619$kurt)
shapiro.test(plots_619$PlotCRR) #normal
# shapiro.test(plots_619$rumple)
# shapiro.test(plots_619$moran)
# shapiro.test(plots_619$geary)

shapiro.test(plots_619$N)
shapiro.test(plots_619$P)
shapiro.test(plots_619$K)
shapiro.test(plots_619$Mg)
shapiro.test(plots_619$Ca)
shapiro.test(plots_619$S)
shapiro.test(plots_619$Zn)
shapiro.test(plots_619$`Mn ppm`)
shapiro.test(plots_619$`Cu ppm`)
shapiro.test(plots_619$`Fe ppm`)
shapiro.test(plots_619$`B ppm`)
shapiro.test(plots_619$`Al ppm`)
shapiro.test(plots_619$`Na %`)

# 
# # NOTE the transformations below are not all natural log. Tukey transformation chooses the best transformation method to maximize Shapiro-Wilks W
plots_619$crrmedian_ln <- transformTukey(plots_619$crrmedian)
plots_619$crriqr_ln <- transformTukey(plots_619$crriqr) # still not normal
plots_619$dip_ln <- transformTukey(plots_619$dip) # still not normal
plots_619$median_ln <- transformTukey(plots_619$median)
plots_619$skew_ln <- transformTukey(plots_619$skew) # still not normal
plots_619$kurt_ln <- transformTukey(plots_619$kurt) # still not normal
plots_619$moran_ln <- transformTukey(plots_619$moran)
plots_619$geary_ln <- transformTukey(plots_619$geary)

library(MVTests)
mvnorm_test <- mvShapiro(plots_619@data[,18:30])
summary(mvnorm_test)

`N_ln` + `P_ln` + `K_ln`+`Mg_ln` +`Ca_ln` + `S_ln` + `Zn ppm_ln`+`Mn ppm_ln`+`Cu ppm_ln`+`Fe ppm_ln`+`B ppm_ln`+`Al ppm_ln`+`Na_ln`

plots_619$N_ln <- transformTukey(plots_619$N)
plots_619$`P_ln` <- transformTukey(plots_619$'P')
plots_619$`K_ln` <- transformTukey(plots_619$`K`)
plots_619$`Mg_ln` <- transformTukey(plots_619$`Mg`)
plots_619$`Ca_ln` <- transformTukey(plots_619$`Ca`)
plots_619$`S_ln` <- transformTukey(plots_619$`S`)
plots_619$`Zn_ln` <- transformTukey(plots_619$`Zn`)
plots_619$`Mn_ln` <- transformTukey(plots_619$`Mn`)
plots_619$`Cu_ln` <- transformTukey(plots_619$`Cu`)
plots_619$`Fe_ln` <- transformTukey(plots_619$`Fe`)
plots_619$`B_ln` <- transformTukey(plots_619$`B`)
plots_619$`Al_ln` <- transformTukey(plots_619$`Al`)
plots_619$`Na_ln` <- transformTukey(plots_619$`Na`)

plots_703$crrmedian_ln <- transformTukey(plots_703$crrmedian)
plots_703$crriqr_ln <- transformTukey(plots_703$crriqr)
plots_703$dip_ln <- transformTukey(plots_703$dip) # still not normal
plots_703$median_ln <- transformTukey(plots_703$median)
plots_703$skew_ln <- transformTukey(plots_703$skew) # still not normal
plots_703$kurt_ln <- transformTukey(plots_703$kurt) # still not normal
plots_703$moran_ln <- transformTukey(plots_703$moran)
plots_703$geary_ln <- transformTukey(plots_703$geary)
plots_703$`N_ln` <- transformTukey(plots_703$`N`)
plots_703$`P_ln` <- transformTukey(plots_703$'P')
plots_703$`K_ln` <- transformTukey(plots_703$`K`)
plots_703$`Mg_ln` <- transformTukey(plots_703$`Mg`)
plots_703$`Ca_ln` <- transformTukey(plots_703$`Ca`)
plots_703$`S_ln` <- transformTukey(plots_703$`S`)
plots_703$`Zn_ln` <- transformTukey(plots_703$`Zn`)
plots_703$`Mn_ln` <- transformTukey(plots_703$`Mn`)
plots_703$`Cu_ln` <- transformTukey(plots_703$`Cu`)
plots_703$`Fe_ln` <- transformTukey(plots_703$`Fe`)
plots_703$`B_ln` <- transformTukey(plots_703$`B`)
plots_703$`Al_ln` <- transformTukey(plots_703$`Al`)
plots_703$`Na_ln` <- transformTukey(plots_703$`Na`)

plots_717$crrmedian_ln <- transformTukey(plots_717$crrmedian)
plots_717$crriqr_ln <- transformTukey(plots_717$crriqr)
plots_717$dip_ln <- transformTukey(plots_717$dip) # still not normal
plots_717$median_ln <- transformTukey(plots_717$median)
plots_717$skew_ln <- transformTukey(plots_717$skew)
plots_717$kurt_ln <- transformTukey(plots_717$kurt)
plots_717$moran_ln <- transformTukey(plots_717$moran) # still not normal
plots_717$geary_ln <- transformTukey(plots_717$geary)
plots_717$`N_ln` <- transformTukey(plots_717$`N`)
plots_717$`P_ln` <- transformTukey(plots_717$'P')
plots_717$`K_ln` <- transformTukey(plots_717$`K`)
plots_717$`Mg_ln` <- transformTukey(plots_717$`Mg`)
plots_717$`Ca_ln` <- transformTukey(plots_717$`Ca`)
plots_717$`S_ln` <- transformTukey(plots_717$`S`)
plots_717$`Zn_ln` <- transformTukey(plots_717$`Zn`)
plots_717$`Mn_ln` <- transformTukey(plots_717$`Mn`)
plots_717$`Cu_ln` <- transformTukey(plots_717$`Cu`)
plots_717$`Fe_ln` <- transformTukey(plots_717$`Fe`)
plots_717$`B_ln` <- transformTukey(plots_717$`B`)
plots_717$`Al_ln` <- transformTukey(plots_717$`Al`)
plots_717$`Na_ln` <- transformTukey(plots_717$`Na`)


# --------------------------------------- Stats data prep --------------------------------------------------------------------------------

#standardize data
std_619 <- as.data.frame(scale(plots_619@data[,c(22:53)]))
std_619$Id <- plots_619@data$plotid
#std_619$Trt_agg <- plots_619@data$Trt_agg
std_619$Treatment <- plots_619@data$treatment
std_619$soil <- plots_619$soil
std_619$soil_GoA <- plots_619$soil_GoA
std_619$soil_MaB <- plots_619$soil_MaB
std_619$soil_Ra <- plots_619$soil_Ra

std_703 <- as.data.frame(scale(plots_703@data[,c(22:53)]))
std_703$Id <- plots_703$plotid
#std_703$Trt_agg <- plots_703$Trt_agg
std_703$Treatment <- plots_703@data$treatment
std_703$soil <- plots_703$soil
std_703$soil_GoA <- plots_703$soil_GoA
std_703$soil_MaB <- plots_703$soil_MaB
std_703$soil_Ra <- plots_703$soil_Ra

std_717 <- as.data.frame(scale(plots_717@data[,c(22:53)]))
std_717$Id <- plots_717$plotid
#std_717$Trt_agg <- plots_717$Trt_agg
std_717$Treatment <- plots_717@data$treatment
std_703$soil <- plots_703$soil
std_717$soil_GoA <- plots_717$soil_GoA
std_717$soil_MaB <- plots_717$soil_MaB
std_717$soil_Ra <- plots_717$soil_Ra



# ---------------------------------------- MANOVA ------------------------------------------------------
med_619 <- std_619$median
moran_619 <- std_619$moran
crrmed_619 <- std_619$crrmedian
crriqr_619 <- std_619$crriqr
skew_619 <- std_619$skew
kurt_619 <- std_619$kurt
rumple_619 <- std_619$rumple

std_619$Treatment <- as.factor(std_619$Treatment)
man_619 <- manova(cbind(med_619, moran_619, crrmed_619, crriqr_619, skew_619, kurt_619, rumple_619)~Treatment, data=std_619)
summary(man_619)
summary.aov(man_619)


med_703 <- std_703$median_ln
moran_703 <- std_703$moran_ln
crrmed_703 <- std_703$crrmedian_ln
skew_703 <- std_703$skew_ln
kurt_703 <- std_703$kurt_ln
rumple_703 <- std_703$rumple

man_703 <- manova(cbind(med_703, moran_703, crrmed_703, skew_703, kurt_703, rumple_703)~Treatment, data=std_703)
summary(man_703)
summary.aov(man_703)


med_717 <- std_717$median_ln
moran_717 <- std_717$moran_ln
crrmed_717 <- std_717$crrmedian_ln
skew_717 <- std_717$skew_ln
kurt_717 <- std_717$kurt_ln
rumple_717 <- std_717$rumple

man_717 <- manova(cbind(med_717, moran_717, crrmed_717, skew_717, kurt_717, rumple_717)~Trt_agg, data=std_717)
summary(man_717)
summary.aov(man_717)

# ------------------------------ Try lm with nutrients as independent ---------------------------------------------

# All nutrients: 
#all nutrients: `N`+`P`+`K`+`Mg`+`Ca`+`S`+`Zn`+`Mn`+`Cu`+`Fe`+`B`+`Al`+`Na`
# nutrients and soil: "`N`+`P`+`K`+`Mg`+`Ca`+`S`+`Zn`+`Mn`+`Cu`+`Fe`+`B`+`Al`+`Na`+`soil_GoA`+`soil_Ra`+`soil_MaB`"

# Function for calculating r^2 in SAR model, null model will be same for each date

null <- spautolm(median ~ 1, data = std_619, listw=nb2listw(w), zero.policy = T, na.action = na.fail, family="SAR")
#null <- lm(BuAc ~ 1, data=norm_1cm, na.action=na.fail)
r2 <- function(model){
  r.squaredLR(model, null = null)[1]
}

adj_r2 <- function(model){
  r2 <- r.squaredLR(model, null = null)[1]
  n <- model$fit[[5]]
  p <- model$parameters
  1-(1-r2)*((n - 1)/(n - p - 1))
}

# 619 OLS
lm_619 <- lm(median ~ `N`+`Ca`+`S`+`Cu`+`B`+`Al`, data = std_619, na.action = na.fail)
summary(lm_619)
qqnorm(lm_619$residuals)
sd(lm_619$residuals)
vif(lm_619) # Colinearity: removed Na, P, Fe, Mn, K
dredge(lm_619, trace = T)

# 619 check for spatial autocorrelation
std_619$ols_res <- lm_619$residuals
# Weights matrix
w <- knn2nb(knearneigh(coordinates(plots), k=8))
moran.test(std_619$ols_res, nb2listw(w))

#Apply SAR model
sar_619 <- spautolm(median ~ `N`+`Ca`+`S`+`Cu`+`B`+`Al`, data = std_619, listw=nb2listw(w), zero.policy = T, na.action = na.fail, family="SAR")
summary(sar_619)
sar_resi <- sar_619$fit[[9]]
std_619$sar_resi <- sar_resi
moran.test(std_619$sar_resi, nb2listw(w))


lm_703 <- lm(median ~ `N`+`P`+`Ca`+`S`+`Cu`+`Fe`+`B`, data = std_703, na.action = na.fail)
summary(lm_703)
vif(lm_703) # Removed Na, Mg, Mn, K, Zn, Al
dredge(lm_703, trace = T)

# 703 check for spatial autocorrelation
std_703$ols_res <- lm_703$residuals
moran.test(std_703$ols_res, nb2listw(w))

#Apply SAR model
sar_703 <- spautolm(median ~ `N`+`P`+`Ca`+`S`+`Cu`+`Fe`+`B`, data = std_703, listw=nb2listw(w), zero.policy = T, na.action = na.fail, family="SAR")
summary(sar_703)
sar_resi <- sar_703$fit[[9]]
std_703$sar_resi <- sar_resi
moran.test(std_703$sar_resi, nb2listw(w))

lm_717 <- lm(median ~ `N`+`P`+`Mg`+`S`+`Mn`+`Fe`+`B`+`Al`, data = std_717, na.action = na.fail)
summary(lm_717)
vif(lm_717) # Eliminated Na, Cu, Ca, Zn, K
dredge(lm_717, trace = T)


# ---------------------------------SAR models with structural metric as response for each time step--------------------------------

r2 <- function(model){
  r.squaredLR(model, null = null)[1]
}

adj_r2 <- function(model){
  r2 <- r.squaredLR(model, null = null)[1]
  n <- model$fit[[5]]
  p <- model$parameters
  1-(1-r2)*((n - 1)/(n - p - 1))
}

w <- knn2nb(knearneigh(coordinates(plots_619), k=8))

# ---------------------------------------SAR 619--------------------------------------------
variables <- "`N`+`Ca`+`S`+`Cu`+`B`+`Al`+`soil_GoA`+`soil_Ra`+`soil_MaB`"
dat <- std_619


metric <- "median"
null <- spautolm( as.formula(paste(metric, 1, sep="~")), data = dat, 
                  listw=nb2listw(w), zero.policy = T, na.action = na.fail, family="SAR")
sar <- spautolm(as.formula(paste(metric, variables, sep = " ~ ")), 
                data = dat, listw=nb2listw(w), zero.policy = T, na.action = na.fail, 
                family="SAR")
summary(sar)
adj_r2(sar)
modsel <- dredge(sar, beta = "none", trace = T, extra = c('r2', 'adj_r2'))
med_619 <- modsel[1]
med_619$depvar <- "median"

metric <- "iqr"
null <- spautolm( as.formula(paste(metric, 1, sep="~")), data = dat, 
                  listw=nb2listw(w), zero.policy = T, na.action = na.fail, family="SAR")
sar <- spautolm(as.formula(paste(metric, variables, sep = " ~ ")), 
                data = dat, listw=nb2listw(w), zero.policy = T, na.action = na.fail, 
                family="SAR")
summary(sar)
adj_r2(sar)
modsel <- dredge(sar, beta = "none", trace = T, extra = c('r2', 'adj_r2'))
iqr_619 <- modsel[1]
iqr_619$depvar <- "iqr"

metric <- "crrmedian"
null <- spautolm( as.formula(paste(metric, 1, sep="~")), data = dat, 
                  listw=nb2listw(w), zero.policy = T, na.action = na.fail, family="SAR")
sar <- spautolm(as.formula(paste(metric, variables, sep = " ~ ")), 
                data = dat, listw=nb2listw(w), zero.policy = T, na.action = na.fail, 
                family="SAR")
summary(sar)
adj_r2(sar)
modsel <- dredge(sar, beta = "none", trace = T, extra = c('r2', 'adj_r2'))
crrmed_619 <- modsel[1]
crrmed_619$depvar <- "crrmed"

metric <- "crriqr"
null <- spautolm( as.formula(paste(metric, 1, sep="~")), data = dat, 
                  listw=nb2listw(w), zero.policy = T, na.action = na.fail, family="SAR")
sar <- spautolm(as.formula(paste(metric, variables, sep = " ~ ")), 
                data = dat, listw=nb2listw(w), zero.policy = T, na.action = na.fail, 
                family="SAR")
summary(sar)
adj_r2(sar)
modsel <- dredge(sar, beta = "none", trace = T, extra = c('r2', 'adj_r2'))
crriqr_619 <- modsel[1]
crriqr_619$depvar <- "crriqr"

metric <- "skew"
null <- spautolm( as.formula(paste(metric, 1, sep="~")), data = dat, 
                  listw=nb2listw(w), zero.policy = T, na.action = na.fail, family="SAR")
sar <- spautolm(as.formula(paste(metric, variables, sep = " ~ ")), 
                data = dat, listw=nb2listw(w), zero.policy = T, na.action = na.fail, 
                family="SAR")
summary(sar)
adj_r2(sar)
modsel <- dredge(sar, beta = "none", trace = T, extra = c('r2', 'adj_r2'))
skew_619 <- modsel[1]
skew_619$depvar <- "skew"

metric <- "kurt"
null <- spautolm( as.formula(paste(metric, 1, sep="~")), data = dat, 
                  listw=nb2listw(w), zero.policy = T, na.action = na.fail, family="SAR")
sar <- spautolm(as.formula(paste(metric, variables, sep = " ~ ")), 
                data = dat, listw=nb2listw(w), zero.policy = T, na.action = na.fail, 
                family="SAR")
summary(sar)
adj_r2(sar)
modsel <- dredge(sar, beta = "none", trace = T, extra = c('r2', 'adj_r2'))
kurt_619 <- modsel[1]

metric <- "PlotCRR"
null <- spautolm( as.formula(paste(metric, 1, sep="~")), data = dat, 
                  listw=nb2listw(w), zero.policy = T, na.action = na.fail, family="SAR")
sar <- spautolm(as.formula(paste(metric, variables, sep = " ~ ")), 
                data = dat, listw=nb2listw(w), zero.policy = T, na.action = na.fail, 
                family="SAR")
summary(sar)
adj_r2(sar)
modsel <- dredge(sar, beta = "none", trace = T, extra = c('r2', 'adj_r2'))
plotcrr_619 <- modsel[1]

metric <- "rumple"
null <- spautolm( as.formula(paste(metric, 1, sep="~")), data = dat, 
                  listw=nb2listw(w), zero.policy = T, na.action = na.fail, family="SAR")
sar <- spautolm(as.formula(paste(metric, variables, sep = " ~ ")), 
                data = dat, listw=nb2listw(w), zero.policy = T, na.action = na.fail, 
                family="SAR")
summary(sar)
adj_r2(sar)
modsel <- dredge(sar, beta = "none", trace = T, extra = c('r2', 'adj_r2'))
rumple_619 <- modsel[1]

metric <- "moran"
null <- spautolm( as.formula(paste(metric, 1, sep="~")), data = dat, 
                  listw=nb2listw(w), zero.policy = T, na.action = na.fail, family="SAR")
sar <- spautolm(as.formula(paste(metric, variables, sep = " ~ ")), 
                data = dat, listw=nb2listw(w), zero.policy = T, na.action = na.fail, 
                family="SAR")
summary(sar)
adj_r2(sar)
modsel <- dredge(sar, beta = "none", trace = T, extra = c('r2', 'adj_r2'))
moran_619 <- modsel[1]


# ---------------------------------------SAR 703--------------------------------------------
variables <- "`N`+`P`+`Ca`+`S`+`Cu`+`Fe`+`B`+`soil_GoA`+`soil_Ra`+`soil_MaB`"

dat <- std_703

metric <- "median"
null <- spautolm( as.formula(paste(metric, 1, sep="~")), data = dat, 
                  listw=nb2listw(w), zero.policy = T, na.action = na.fail, family="SAR")
sar <- spautolm(as.formula(paste(metric, variables, sep = " ~ ")), 
                data = dat, listw=nb2listw(w), zero.policy = T, na.action = na.fail, 
                family="SAR")
summary(sar)
adj_r2(sar)
modsel <- dredge(sar, beta = "none", trace = T, extra = c('r2', 'adj_r2'))
med_703 <- modsel[1]

metric <- "iqr"
null <- spautolm( as.formula(paste(metric, 1, sep="~")), data = dat, 
                  listw=nb2listw(w), zero.policy = T, na.action = na.fail, family="SAR")
sar <- spautolm(as.formula(paste(metric, variables, sep = " ~ ")), 
                data = dat, listw=nb2listw(w), zero.policy = T, na.action = na.fail, 
                family="SAR")
summary(sar)
adj_r2(sar)
modsel <- dredge(sar, beta = "none", trace = T, extra = c('r2', 'adj_r2'))
iqr_703 <- modsel[1]

metric <- "crrmedian"
null <- spautolm( as.formula(paste(metric, 1, sep="~")), data = dat, 
                  listw=nb2listw(w), zero.policy = T, na.action = na.fail, family="SAR")
sar <- spautolm(as.formula(paste(metric, variables, sep = " ~ ")), 
                data = dat, listw=nb2listw(w), zero.policy = T, na.action = na.fail, 
                family="SAR")
summary(sar)
adj_r2(sar)
modsel <- dredge(sar, beta = "none", trace = T, extra = c('r2', 'adj_r2'))
crrmed_703 <- modsel[1]

metric <- "crriqr"
null <- spautolm( as.formula(paste(metric, 1, sep="~")), data = dat, 
                  listw=nb2listw(w), zero.policy = T, na.action = na.fail, family="SAR")
sar <- spautolm(as.formula(paste(metric, variables, sep = " ~ ")), 
                data = dat, listw=nb2listw(w), zero.policy = T, na.action = na.fail, 
                family="SAR")
summary(sar)
adj_r2(sar)
modsel <- dredge(sar, beta = "none", trace = T, extra = c('r2', 'adj_r2'))
crriqr_703 <- modsel[1]

metric <- "skew"
null <- spautolm( as.formula(paste(metric, 1, sep="~")), data = dat, 
                  listw=nb2listw(w), zero.policy = T, na.action = na.fail, family="SAR")
sar <- spautolm(as.formula(paste(metric, variables, sep = " ~ ")), 
                data = dat, listw=nb2listw(w), zero.policy = T, na.action = na.fail, 
                family="SAR")
summary(sar)
adj_r2(sar)
modsel <- dredge(sar, beta = "none", trace = T, extra = c('r2', 'adj_r2'))
skew_703 <- modsel[1]

metric <- "kurt"
null <- spautolm( as.formula(paste(metric, 1, sep="~")), data = dat, 
                  listw=nb2listw(w), zero.policy = T, na.action = na.fail, family="SAR")
sar <- spautolm(as.formula(paste(metric, variables, sep = " ~ ")), 
                data = dat, listw=nb2listw(w), zero.policy = T, na.action = na.fail, 
                family="SAR")
summary(sar)
adj_r2(sar)
modsel <- dredge(sar, beta = "none", trace = T, extra = c('r2', 'adj_r2'))
kurt_703 <- modsel[1]

metric <- "PlotCRR"
null <- spautolm( as.formula(paste(metric, 1, sep="~")), data = dat, 
                  listw=nb2listw(w), zero.policy = T, na.action = na.fail, family="SAR")
sar <- spautolm(as.formula(paste(metric, variables, sep = " ~ ")), 
                data = dat, listw=nb2listw(w), zero.policy = T, na.action = na.fail, 
                family="SAR")
summary(sar)
adj_r2(sar)
modsel <- dredge(sar, beta = "none", trace = T, extra = c('r2', 'adj_r2'))
plotcrr_703 <- modsel[1]

metric <- "rumple"
null <- spautolm( as.formula(paste(metric, 1, sep="~")), data = dat, 
                  listw=nb2listw(w), zero.policy = T, na.action = na.fail, family="SAR")
sar <- spautolm(as.formula(paste(metric, variables, sep = " ~ ")), 
                data = dat, listw=nb2listw(w), zero.policy = T, na.action = na.fail, 
                family="SAR")
summary(sar)
adj_r2(sar)
modsel <- dredge(sar, beta = "none", trace = T, extra = c('r2', 'adj_r2'))
rumple_703 <- modsel[1]

metric <- "moran"
null <- spautolm( as.formula(paste(metric, 1, sep="~")), data = dat, 
                  listw=nb2listw(w), zero.policy = T, na.action = na.fail, family="SAR")
sar <- spautolm(as.formula(paste(metric, variables, sep = " ~ ")), 
                data = dat, listw=nb2listw(w), zero.policy = T, na.action = na.fail, 
                family="SAR")
summary(sar)
adj_r2(sar)
modsel <- dredge(sar, beta = "none", trace = T, extra = c('r2', 'adj_r2'))
moran_703 <- modsel[1]

# ---------------------------------------SAR 717--------------------------------------------
variables <- "`N`+`P`+`Mg`+`S`+`Mn`+`Fe`+`B`+`Al`"

dat <- std_717

metric <- "median"
null <- spautolm( as.formula(paste(metric, 1, sep="~")), data = dat, 
                  listw=nb2listw(w), zero.policy = T, na.action = na.fail, family="SAR")
sar <- spautolm(as.formula(paste(metric, variables, sep = " ~ ")), 
                data = dat, listw=nb2listw(w), zero.policy = T, na.action = na.fail, 
                family="SAR")
summary(sar)
adj_r2(sar)
modsel <- dredge(sar, beta = "none", trace = T, extra = c('r2', 'adj_r2'))
med_717 <- modsel[1]

metric <- "iqr"
null <- spautolm( as.formula(paste(metric, 1, sep="~")), data = dat, 
                  listw=nb2listw(w), zero.policy = T, na.action = na.fail, family="SAR")
sar <- spautolm(as.formula(paste(metric, variables, sep = " ~ ")), 
                data = dat, listw=nb2listw(w), zero.policy = T, na.action = na.fail, 
                family="SAR")
summary(sar)
adj_r2(sar)
modsel <- dredge(sar, beta = "none", trace = T, extra = c('r2', 'adj_r2'))
iqr_717 <- modsel[1]

metric <- "crrmedian"
null <- spautolm( as.formula(paste(metric, 1, sep="~")), data = dat, 
                  listw=nb2listw(w), zero.policy = T, na.action = na.fail, family="SAR")
sar <- spautolm(as.formula(paste(metric, variables, sep = " ~ ")), 
                data = dat, listw=nb2listw(w), zero.policy = T, na.action = na.fail, 
                family="SAR")
summary(sar)
adj_r2(sar)
modsel <- dredge(sar, beta = "none", trace = T, extra = c('r2', 'adj_r2'))
crrmed_717 <- modsel[1]

metric <- "crriqr"
null <- spautolm( as.formula(paste(metric, 1, sep="~")), data = dat, 
                  listw=nb2listw(w), zero.policy = T, na.action = na.fail, family="SAR")
sar <- spautolm(as.formula(paste(metric, variables, sep = " ~ ")), 
                data = dat, listw=nb2listw(w), zero.policy = T, na.action = na.fail, 
                family="SAR")
summary(sar)
adj_r2(sar)
modsel <- dredge(sar, beta = "none", trace = T, extra = c('r2', 'adj_r2'))
crriqr_717 <- modsel[1]

metric <- "skew"
null <- spautolm( as.formula(paste(metric, 1, sep="~")), data = dat, 
                  listw=nb2listw(w), zero.policy = T, na.action = na.fail, family="SAR")
sar <- spautolm(as.formula(paste(metric, variables, sep = " ~ ")), 
                data = dat, listw=nb2listw(w), zero.policy = T, na.action = na.fail, 
                family="SAR")
summary(sar)
adj_r2(sar)
modsel <- dredge(sar, beta = "none", trace = T, extra = c('r2', 'adj_r2'))
skew_717 <- modsel[1]

metric <- "kurt"
null <- spautolm( as.formula(paste(metric, 1, sep="~")), data = dat, 
                  listw=nb2listw(w), zero.policy = T, na.action = na.fail, family="SAR")
sar <- spautolm(as.formula(paste(metric, variables, sep = " ~ ")), 
                data = dat, listw=nb2listw(w), zero.policy = T, na.action = na.fail, 
                family="SAR")
summary(sar)
adj_r2(sar)
modsel <- dredge(sar, beta = "none", trace = T, extra = c('r2', 'adj_r2'))
kurt_717 <- modsel[1]

metric <- "PlotCRR"
null <- spautolm( as.formula(paste(metric, 1, sep="~")), data = dat, 
                  listw=nb2listw(w), zero.policy = T, na.action = na.fail, family="SAR")
sar <- spautolm(as.formula(paste(metric, variables, sep = " ~ ")), 
                data = dat, listw=nb2listw(w), zero.policy = T, na.action = na.fail, 
                family="SAR")
summary(sar)
adj_r2(sar)
modsel <- dredge(sar, beta = "none", trace = T, extra = c('r2', 'adj_r2'))
plotcrr_717 <- modsel[1]

metric <- "rumple"
null <- spautolm( as.formula(paste(metric, 1, sep="~")), data = dat, 
                  listw=nb2listw(w), zero.policy = T, na.action = na.fail, family="SAR")
sar <- spautolm(as.formula(paste(metric, variables, sep = " ~ ")), 
                data = dat, listw=nb2listw(w), zero.policy = T, na.action = na.fail, 
                family="SAR")
summary(sar)
adj_r2(sar)
modsel <- dredge(sar, beta = "none", trace = T, extra = c('r2', 'adj_r2'))
rumple_717 <- modsel[1]

metric <- "moran"
null <- spautolm( as.formula(paste(metric, 1, sep="~")), data = dat, 
                  listw=nb2listw(w), zero.policy = T, na.action = na.fail, family="SAR")
sar <- spautolm(as.formula(paste(metric, variables, sep = " ~ ")), 
                data = dat, listw=nb2listw(w), zero.policy = T, na.action = na.fail, 
                family="SAR")
summary(sar)
adj_r2(sar)
modsel <- dredge(sar, beta = "none", trace = T, extra = c('r2', 'adj_r2'))
moran_717 <- modsel[1]

#---------------------------------------SAR nutrient results tables --------------------------------------

med_619$depvar <- "median"
iqr_619$depvar <- "iqr"
crrmed_619$depvar <- "crrmed"
crriqr_619$depvar <- "crriqr"
skew_619$depvar <- "skew"
kurt_619$depvar <- "kurt"
plotcrr_619$depvar <- "plotcrr"
rumple_619$depvar <- "rumple"
moran_619$depvar <- "moran"

med_703$depvar <- "median"
iqr_703$depvar <- "iqr"
crrmed_703$depvar <- "crrmed"
crriqr_703$depvar <- "crriqr"
skew_703$depvar <- "skew"
kurt_703$depvar <- "kurt"
plotcrr_703$depvar <- "plotcrr"
rumple_703$depvar <- "rumple"
moran_703$depvar <- "moran"

med_717$depvar <- "median"
iqr_717$depvar <- "iqr"
crrmed_717$depvar <- "crrmed"
crriqr_717$depvar <- "crriqr"
skew_717$depvar <- "skew"
kurt_717$depvar <- "kurt"
plotcrr_717$depvar <- "plotcrr"
rumple_717$depvar <- "rumple"
moran_717$depvar <- "moran"

sarresults_nutr_619 <- rbind(med_619, iqr_619, crrmed_619, crriqr_619, skew_619, kurt_619, plotcrr_619, rumple_619, moran_619)
sarresults_nutr_619$depvar <- c("median", "iqr", "crrmed", "crriqr", "skew", "kurt", "plotcrr", "rumple", "moran")

sarresults_nutr_703 <- rbind(med_703, iqr_703, crrmed_703, crriqr_703, skew_703, kurt_703, plotcrr_703, rumple_703, moran_703)
sarresults_nutr_703$depvar <- c("median", "iqr", "crrmed", "crriqr", "skew", "kurt", "plotcrr", "rumple", "moran")

sarresults_nutr_717 <- rbind(med_717, iqr_717, crrmed_717, crriqr_717, skew_717, kurt_717, plotcrr_717, rumple_717, moran_717)
sarresults_nutr_717$depvar <- c("median", "iqr", "crrmed", "crriqr", "skew", "kurt", "plotcrr", "rumple", "moran")

write.csv(sarresults_nutr_619, "Q:/My Drive/Research/Canopy_Morphology/Tobacco_Project/2018/canopy_analysis/sarresults_nutr_619.csv")
write.csv(sarresults_nutr_703, "Q:/My Drive/Research/Canopy_Morphology/Tobacco_Project/2018/canopy_analysis/sarresults_nutr_703.csv")
write.csv(sarresults_nutr_717, "Q:/My Drive/Research/Canopy_Morphology/Tobacco_Project/2018/canopy_analysis/sarresults_nutr_717.csv")


# ---------------------------------SAR models with nutrients as response for each time step--------------------------------

variables <- "crrmedian+crriqr+median+skew+iqr+kurt+rumple+moran+soil_GoA+soil_Ra+soil_MaB+watersum"



r2 <- function(model){
  r.squaredLR(model, null = null)[1]
}

adj_r2 <- function(model){
  r2 <- r.squaredLR(model, null = null)[1]
  n <- model$fit[[5]]
  p <- model$parameters
  1-(1-r2)*((n - 1)/(n - p - 1))
}

w <- knn2nb(knearneigh(coordinates(plots_619), k=8))

# ---------------------------------------SAR 619--------------------------------------------
dat <- std_619


metric <- "N"
null <- spautolm( as.formula(paste(metric, 1, sep="~")), data = dat, 
                  listw=nb2listw(w), zero.policy = T, na.action = na.fail, family="SAR")
sar <- spautolm(as.formula(paste(metric, variables, sep = " ~ ")), 
                data = dat, listw=nb2listw(w), zero.policy = T, na.action = na.fail, 
                family="SAR")
summary(sar)
adj_r2(sar)
modsel <- dredge(sar, beta = "none", trace = T, extra = c('r2', 'adj_r2'))
n_619 <- modsel[1]

tmp <-spautolm(N~median+moran, 
         data = dat, listw=nb2listw(w), zero.policy = T, na.action = na.fail, 
         family="SAR")
sar_resi <- tmp$fit[[9]]
std_619$sar_resi <- sar_resi
plot(std_619$sar_resi)

metric <- "P"
null <- spautolm( as.formula(paste(metric, 1, sep="~")), data = dat, 
                  listw=nb2listw(w), zero.policy = T, na.action = na.fail, family="SAR")
sar <- spautolm(as.formula(paste(metric, variables, sep = " ~ ")), 
                data = dat, listw=nb2listw(w), zero.policy = T, na.action = na.fail, 
                family="SAR")
summary(sar)
adj_r2(sar)
modsel <- dredge(sar, beta = "none", trace = T, extra = c('r2', 'adj_r2'))
p_619 <- modsel[1]

metric <- "K"
null <- spautolm( as.formula(paste(metric, 1, sep="~")), data = dat, 
                  listw=nb2listw(w), zero.policy = T, na.action = na.fail, family="SAR")
sar <- spautolm(as.formula(paste(metric, variables, sep = " ~ ")), 
                data = dat, listw=nb2listw(w), zero.policy = T, na.action = na.fail, 
                family="SAR")
summary(sar)
adj_r2(sar)
modsel <- dredge(sar, beta = "none", trace = T, extra = c('r2', 'adj_r2'))
k_619 <- modsel[1]


metric <- "Ca"
null <- spautolm( as.formula(paste(metric, 1, sep="~")), data = dat, 
                  listw=nb2listw(w), zero.policy = T, na.action = na.fail, family="SAR")
sar <- spautolm(as.formula(paste(metric, variables, sep = " ~ ")), 
                data = dat, listw=nb2listw(w), zero.policy = T, na.action = na.fail, 
                family="SAR")
summary(sar)
adj_r2(sar)
modsel <- dredge(sar, beta = "none", trace = T, extra = c('r2', 'adj_r2'))
ca_619 <- modsel[1]

metric <- "Mg"
null <- spautolm( as.formula(paste(metric, 1, sep="~")), data = dat, 
                  listw=nb2listw(w), zero.policy = T, na.action = na.fail, family="SAR")
sar <- spautolm(as.formula(paste(metric, variables, sep = " ~ ")), 
                data = dat, listw=nb2listw(w), zero.policy = T, na.action = na.fail, 
                family="SAR")
summary(sar)
adj_r2(sar)
modsel <- dredge(sar, beta = "none", trace = T, extra = c('r2', 'adj_r2'))
mg_619 <- modsel[1]

metric <- "B"
null <- spautolm( as.formula(paste(metric, 1, sep="~")), data = dat, 
                  listw=nb2listw(w), zero.policy = T, na.action = na.fail, family="SAR")
sar <- spautolm(as.formula(paste(metric, variables, sep = " ~ ")), 
                data = dat, listw=nb2listw(w), zero.policy = T, na.action = na.fail, 
                family="SAR")
summary(sar)
adj_r2(sar)
modsel <- dredge(sar, beta = "none", trace = T, extra = c('r2', 'adj_r2'))
b_619 <- modsel[1]

metric <- "Mn"
null <- spautolm( as.formula(paste(metric, 1, sep="~")), data = dat, 
                  listw=nb2listw(w), zero.policy = T, na.action = na.fail, family="SAR")
sar <- spautolm(as.formula(paste(metric, variables, sep = " ~ ")), 
                data = dat, listw=nb2listw(w), zero.policy = T, na.action = na.fail, 
                family="SAR")
summary(sar)
adj_r2(sar)
modsel <- dredge(sar, beta = "none", trace = T, extra = c('r2', 'adj_r2'))
mn_619 <- modsel[1]

metric <- "S"
null <- spautolm( as.formula(paste(metric, 1, sep="~")), data = dat, 
                  listw=nb2listw(w), zero.policy = T, na.action = na.fail, family="SAR")
sar <- spautolm(as.formula(paste(metric, variables, sep = " ~ ")), 
                data = dat, listw=nb2listw(w), zero.policy = T, na.action = na.fail, 
                family="SAR")
summary(sar)
adj_r2(sar)
modsel <- dredge(sar, beta = "none", trace = T, extra = c('r2', 'adj_r2'))
s_619 <- modsel[1]

metric <- "Zn"
null <- spautolm( as.formula(paste(metric, 1, sep="~")), data = dat, 
                  listw=nb2listw(w), zero.policy = T, na.action = na.fail, family="SAR")
sar <- spautolm(as.formula(paste(metric, variables, sep = " ~ ")), 
                data = dat, listw=nb2listw(w), zero.policy = T, na.action = na.fail, 
                family="SAR")
summary(sar)
adj_r2(sar)
modsel <- dredge(sar, beta = "none", trace = T, extra = c('r2', 'adj_r2'))
zn_619 <- modsel[1]

metric <- "Cu"
null <- spautolm( as.formula(paste(metric, 1, sep="~")), data = dat, 
                  listw=nb2listw(w), zero.policy = T, na.action = na.fail, family="SAR")
sar <- spautolm(as.formula(paste(metric, variables, sep = " ~ ")), 
                data = dat, listw=nb2listw(w), zero.policy = T, na.action = na.fail, 
                family="SAR")
summary(sar)
adj_r2(sar)
modsel <- dredge(sar, beta = "none", trace = T, extra = c('r2', 'adj_r2'))
cu_619 <- modsel[1]

metric <- "Fe"
null <- spautolm( as.formula(paste(metric, 1, sep="~")), data = dat, 
                  listw=nb2listw(w), zero.policy = T, na.action = na.fail, family="SAR")
sar <- spautolm(as.formula(paste(metric, variables, sep = " ~ ")), 
                data = dat, listw=nb2listw(w), zero.policy = T, na.action = na.fail, 
                family="SAR")
summary(sar)
adj_r2(sar)
modsel <- dredge(sar, beta = "none", trace = T, extra = c('r2', 'adj_r2'))
fe_619 <- modsel[1]

metric <- "Al"
null <- spautolm( as.formula(paste(metric, 1, sep="~")), data = dat, 
                  listw=nb2listw(w), zero.policy = T, na.action = na.fail, family="SAR")
sar <- spautolm(as.formula(paste(metric, variables, sep = " ~ ")), 
                data = dat, listw=nb2listw(w), zero.policy = T, na.action = na.fail, 
                family="SAR")
summary(sar)
adj_r2(sar)
modsel <- dredge(sar, beta = "none", trace = T, extra = c('r2', 'adj_r2'))
al_619 <- modsel[1]

# ---------------------------------------SAR 703--------------------------------------------
dat <- std_703


metric <- "N"
null <- spautolm( as.formula(paste(metric, 1, sep="~")), data = dat, 
                  listw=nb2listw(w), zero.policy = T, na.action = na.fail, family="SAR")
sar <- spautolm(as.formula(paste(metric, variables, sep = " ~ ")), 
                data = dat, listw=nb2listw(w), zero.policy = T, na.action = na.fail, 
                family="SAR")
summary(sar)
adj_r2(sar)
modsel <- dredge(sar, beta = "none", trace = T, extra = c('r2', 'adj_r2'))
n_703 <- modsel[1]

metric <- "P"
null <- spautolm( as.formula(paste(metric, 1, sep="~")), data = dat, 
                  listw=nb2listw(w), zero.policy = T, na.action = na.fail, family="SAR")
sar <- spautolm(as.formula(paste(metric, variables, sep = " ~ ")), 
                data = dat, listw=nb2listw(w), zero.policy = T, na.action = na.fail, 
                family="SAR")
summary(sar)
adj_r2(sar)
modsel <- dredge(sar, beta = "none", trace = T, extra = c('r2', 'adj_r2'))
p_703 <- modsel[1]

metric <- "K"
null <- spautolm( as.formula(paste(metric, 1, sep="~")), data = dat, 
                  listw=nb2listw(w), zero.policy = T, na.action = na.fail, family="SAR")
sar <- spautolm(as.formula(paste(metric, variables, sep = " ~ ")), 
                data = dat, listw=nb2listw(w), zero.policy = T, na.action = na.fail, 
                family="SAR")
summary(sar)
adj_r2(sar)
modsel <- dredge(sar, beta = "none", trace = T, extra = c('r2', 'adj_r2'))
k_703 <- modsel[1]

metric <- "Ca"
null <- spautolm( as.formula(paste(metric, 1, sep="~")), data = dat, 
                  listw=nb2listw(w), zero.policy = T, na.action = na.fail, family="SAR")
sar <- spautolm(as.formula(paste(metric, variables, sep = " ~ ")), 
                data = dat, listw=nb2listw(w), zero.policy = T, na.action = na.fail, 
                family="SAR")
summary(sar)
adj_r2(sar)
modsel <- dredge(sar, beta = "none", trace = T, extra = c('r2', 'adj_r2'))
ca_703 <- modsel[1]

metric <- "Mg"
null <- spautolm( as.formula(paste(metric, 1, sep="~")), data = dat, 
                  listw=nb2listw(w), zero.policy = T, na.action = na.fail, family="SAR")
sar <- spautolm(as.formula(paste(metric, variables, sep = " ~ ")), 
                data = dat, listw=nb2listw(w), zero.policy = T, na.action = na.fail, 
                family="SAR")
summary(sar)
adj_r2(sar)
modsel <- dredge(sar, beta = "none", trace = T, extra = c('r2', 'adj_r2'))
mg_703 <- modsel[1]

metric <- "B"
null <- spautolm( as.formula(paste(metric, 1, sep="~")), data = dat, 
                  listw=nb2listw(w), zero.policy = T, na.action = na.fail, family="SAR")
sar <- spautolm(as.formula(paste(metric, variables, sep = " ~ ")), 
                data = dat, listw=nb2listw(w), zero.policy = T, na.action = na.fail, 
                family="SAR")
summary(sar)
adj_r2(sar)
modsel <- dredge(sar, beta = "none", trace = T, extra = c('r2', 'adj_r2'))
b_703 <- modsel[1]

metric <- "Mn"
null <- spautolm( as.formula(paste(metric, 1, sep="~")), data = dat, 
                  listw=nb2listw(w), zero.policy = T, na.action = na.fail, family="SAR")
sar <- spautolm(as.formula(paste(metric, variables, sep = " ~ ")), 
                data = dat, listw=nb2listw(w), zero.policy = T, na.action = na.fail, 
                family="SAR")
summary(sar)
adj_r2(sar)
modsel <- dredge(sar, beta = "none", trace = T, extra = c('r2', 'adj_r2'))
mn_703 <- modsel[1]

metric <- "S"
null <- spautolm( as.formula(paste(metric, 1, sep="~")), data = dat, 
                  listw=nb2listw(w), zero.policy = T, na.action = na.fail, family="SAR")
sar <- spautolm(as.formula(paste(metric, variables, sep = " ~ ")), 
                data = dat, listw=nb2listw(w), zero.policy = T, na.action = na.fail, 
                family="SAR")
summary(sar)
adj_r2(sar)
modsel <- dredge(sar, beta = "none", trace = T, extra = c('r2', 'adj_r2'))
s_703 <- modsel[1]

metric <- "Zn"
null <- spautolm( as.formula(paste(metric, 1, sep="~")), data = dat, 
                  listw=nb2listw(w), zero.policy = T, na.action = na.fail, family="SAR")
sar <- spautolm(as.formula(paste(metric, variables, sep = " ~ ")), 
                data = dat, listw=nb2listw(w), zero.policy = T, na.action = na.fail, 
                family="SAR")
summary(sar)
adj_r2(sar)
modsel <- dredge(sar, beta = "none", trace = T, extra = c('r2', 'adj_r2'))
zn_703 <- modsel[1]

metric <- "Cu"
null <- spautolm( as.formula(paste(metric, 1, sep="~")), data = dat, 
                  listw=nb2listw(w), zero.policy = T, na.action = na.fail, family="SAR")
sar <- spautolm(as.formula(paste(metric, variables, sep = " ~ ")), 
                data = dat, listw=nb2listw(w), zero.policy = T, na.action = na.fail, 
                family="SAR")
summary(sar)
adj_r2(sar)
modsel <- dredge(sar, beta = "none", trace = T, extra = c('r2', 'adj_r2'))
cu_703 <- modsel[1]

metric <- "Fe"
null <- spautolm( as.formula(paste(metric, 1, sep="~")), data = dat, 
                  listw=nb2listw(w), zero.policy = T, na.action = na.fail, family="SAR")
sar <- spautolm(as.formula(paste(metric, variables, sep = " ~ ")), 
                data = dat, listw=nb2listw(w), zero.policy = T, na.action = na.fail, 
                family="SAR")
summary(sar)
adj_r2(sar)
modsel <- dredge(sar, beta = "none", trace = T, extra = c('r2', 'adj_r2'))
fe_703 <- modsel[1]

metric <- "Al"
null <- spautolm( as.formula(paste(metric, 1, sep="~")), data = dat, 
                  listw=nb2listw(w), zero.policy = T, na.action = na.fail, family="SAR")
sar <- spautolm(as.formula(paste(metric, variables, sep = " ~ ")), 
                data = dat, listw=nb2listw(w), zero.policy = T, na.action = na.fail, 
                family="SAR")
summary(sar)
adj_r2(sar)
modsel <- dredge(sar, beta = "none", trace = T, extra = c('r2', 'adj_r2'))
al_703 <- modsel[1]


# ---------------------------------------SAR 717--------------------------------------------
dat <- std_717


metric <- "N"
null <- spautolm( as.formula(paste(metric, 1, sep="~")), data = dat, 
                  listw=nb2listw(w), zero.policy = T, na.action = na.fail, family="SAR")
sar <- spautolm(as.formula(paste(metric, variables, sep = " ~ ")), 
                data = dat, listw=nb2listw(w), zero.policy = T, na.action = na.fail, 
                family="SAR")
summary(sar)
adj_r2(sar)
modsel <- dredge(sar, beta = "none", trace = T, extra = c('r2', 'adj_r2'))
n_717 <- modsel[1]

metric <- "P"
null <- spautolm( as.formula(paste(metric, 1, sep="~")), data = dat, 
                  listw=nb2listw(w), zero.policy = T, na.action = na.fail, family="SAR")
sar <- spautolm(as.formula(paste(metric, variables, sep = " ~ ")), 
                data = dat, listw=nb2listw(w), zero.policy = T, na.action = na.fail, 
                family="SAR")
summary(sar)
adj_r2(sar)
modsel <- dredge(sar, beta = "none", trace = T, extra = c('r2', 'adj_r2'))
p_717 <- modsel[1]

metric <- "K"
null <- spautolm( as.formula(paste(metric, 1, sep="~")), data = dat, 
                  listw=nb2listw(w), zero.policy = T, na.action = na.fail, family="SAR")
sar <- spautolm(as.formula(paste(metric, variables, sep = " ~ ")), 
                data = dat, listw=nb2listw(w), zero.policy = T, na.action = na.fail, 
                family="SAR")
summary(sar)
adj_r2(sar)
modsel <- dredge(sar, beta = "none", trace = T, extra = c('r2', 'adj_r2'))
k_717 <- modsel[1]

metric <- "Ca"
null <- spautolm( as.formula(paste(metric, 1, sep="~")), data = dat, 
                  listw=nb2listw(w), zero.policy = T, na.action = na.fail, family="SAR")
sar <- spautolm(as.formula(paste(metric, variables, sep = " ~ ")), 
                data = dat, listw=nb2listw(w), zero.policy = T, na.action = na.fail, 
                family="SAR")
summary(sar)
adj_r2(sar)
modsel <- dredge(sar, beta = "none", trace = T, extra = c('r2', 'adj_r2'))
ca_717 <- modsel[1]

metric <- "Mg"
null <- spautolm( as.formula(paste(metric, 1, sep="~")), data = dat, 
                  listw=nb2listw(w), zero.policy = T, na.action = na.fail, family="SAR")
sar <- spautolm(as.formula(paste(metric, variables, sep = " ~ ")), 
                data = dat, listw=nb2listw(w), zero.policy = T, na.action = na.fail, 
                family="SAR")
summary(sar)
adj_r2(sar)
modsel <- dredge(sar, beta = "none", trace = T, extra = c('r2', 'adj_r2'))
mg_717 <- modsel[1]

metric <- "B"
null <- spautolm( as.formula(paste(metric, 1, sep="~")), data = dat, 
                  listw=nb2listw(w), zero.policy = T, na.action = na.fail, family="SAR")
sar <- spautolm(as.formula(paste(metric, variables, sep = " ~ ")), 
                data = dat, listw=nb2listw(w), zero.policy = T, na.action = na.fail, 
                family="SAR")
summary(sar)
adj_r2(sar)
modsel <- dredge(sar, beta = "none", trace = T, extra = c('r2', 'adj_r2'))
b_717 <- modsel[1]

metric <- "Mn"
null <- spautolm( as.formula(paste(metric, 1, sep="~")), data = dat, 
                  listw=nb2listw(w), zero.policy = T, na.action = na.fail, family="SAR")
sar <- spautolm(as.formula(paste(metric, variables, sep = " ~ ")), 
                data = dat, listw=nb2listw(w), zero.policy = T, na.action = na.fail, 
                family="SAR")
summary(sar)
adj_r2(sar)
modsel <- dredge(sar, beta = "none", trace = T, extra = c('r2', 'adj_r2'))
mn_717 <- modsel[1]

metric <- "S"
null <- spautolm( as.formula(paste(metric, 1, sep="~")), data = dat, 
                  listw=nb2listw(w), zero.policy = T, na.action = na.fail, family="SAR")
sar <- spautolm(as.formula(paste(metric, variables, sep = " ~ ")), 
                data = dat, listw=nb2listw(w), zero.policy = T, na.action = na.fail, 
                family="SAR")
summary(sar)
adj_r2(sar)
modsel <- dredge(sar, beta = "none", trace = T, extra = c('r2', 'adj_r2'))
s_717 <- modsel[1]

metric <- "Zn"
null <- spautolm( as.formula(paste(metric, 1, sep="~")), data = dat, 
                  listw=nb2listw(w), zero.policy = T, na.action = na.fail, family="SAR")
sar <- spautolm(as.formula(paste(metric, variables, sep = " ~ ")), 
                data = dat, listw=nb2listw(w), zero.policy = T, na.action = na.fail, 
                family="SAR")
summary(sar)
adj_r2(sar)
modsel <- dredge(sar, beta = "none", trace = T, extra = c('r2', 'adj_r2'))
zn_717 <- modsel[1]

metric <- "Cu"
null <- spautolm( as.formula(paste(metric, 1, sep="~")), data = dat, 
                  listw=nb2listw(w), zero.policy = T, na.action = na.fail, family="SAR")
sar <- spautolm(as.formula(paste(metric, variables, sep = " ~ ")), 
                data = dat, listw=nb2listw(w), zero.policy = T, na.action = na.fail, 
                family="SAR")
summary(sar)
adj_r2(sar)
modsel <- dredge(sar, beta = "none", trace = T, extra = c('r2', 'adj_r2'))
cu_717 <- modsel[1]

metric <- "Fe"
null <- spautolm( as.formula(paste(metric, 1, sep="~")), data = dat, 
                  listw=nb2listw(w), zero.policy = T, na.action = na.fail, family="SAR")
sar <- spautolm(as.formula(paste(metric, variables, sep = " ~ ")), 
                data = dat, listw=nb2listw(w), zero.policy = T, na.action = na.fail, 
                family="SAR")
summary(sar)
adj_r2(sar)
modsel <- dredge(sar, beta = "none", trace = T, extra = c('r2', 'adj_r2'))
fe_717 <- modsel[1]

metric <- "Al"
null <- spautolm( as.formula(paste(metric, 1, sep="~")), data = dat, 
                  listw=nb2listw(w), zero.policy = T, na.action = na.fail, family="SAR")
sar <- spautolm(as.formula(paste(metric, variables, sep = " ~ ")), 
                data = dat, listw=nb2listw(w), zero.policy = T, na.action = na.fail, 
                family="SAR")
summary(sar)
adj_r2(sar)
modsel <- dredge(sar, beta = "none", trace = T, extra = c('r2', 'adj_r2'))
al_717 <- modsel[1]


# --------------------------------------------SAR structure results tables -------------------------------------
sarresults_struct_619 <- rbind(n_619, p_619, k_619, ca_619, mg_619, b_619, mn_619, s_619, zn_619, cu_619, fe_619, al_619)
sarresults_struct_619$depvar <- c("N", "P", "K", "Ca", "Mg", "B", "Mn", "S", "Zn", "Cu", "Fe", "Al")

sarresults_struct_703 <- rbind(n_703, p_703, k_703, ca_703, mg_703, b_703, mn_703, s_703, zn_703, cu_703, fe_703, al_703)
sarresults_struct_703$depvar <- c("N", "P", "K", "Ca", "Mg", "B", "Mn", "S", "Zn", "Cu", "Fe", "Al")

sarresults_struct_717 <- rbind(n_717, p_717, k_717, ca_717, mg_717, b_717, mn_717, s_717, zn_717, cu_717, fe_717, al_717)
sarresults_struct_717$depvar <- c("N", "P", "K", "Ca", "Mg", "B", "Mn", "S", "Zn", "Cu", "Fe", "Al")

write.csv(sarresults_struct_619, "Q:/My Drive/Research/Canopy_Morphology/Tobacco_Project/2019/sarresults_struct_619.csv")
write.csv(sarresults_struct_703, "Q:/My Drive/Research/Canopy_Morphology/Tobacco_Project/2019/sarresults_struct_703.csv")
write.csv(sarresults_struct_717, "Q:/My Drive/Research/Canopy_Morphology/Tobacco_Project/2019/sarresults_struct_717.csv")

#--------------------------------------------Linear mixed effect models--------------------------------------------------
library(lme4)

lme_med_619 <- lmer(median ~ N + K + P + Ca + S + B | soil)


#----------------------------------------------Try canonical correlation---------------------
nut_619 <- plots_619@data[, 22:34]
struc_619 <- plots_619@data[,35:50]
cancor_619 <- cancor(struc_619, nut_619)
cancor_619

#--------------------------------------------Try NN regression---------------------------------------------------------
library(neuralnet)
set.seed(42)
train_data_indices <- rep(FALSE, nrow(std_619))
train_data_indices[sample(1:nrow(std_619), round(0.7 * nrow(std_619)))] <- TRUE # randomly select 70% of the data for training
nn <- neuralnet(median ~ `N`+`Ca`+`S`+`Cu`+`B`+`Al`+`soil_GoA`+`soil_Ra`+`soil_MaB`+`watersum`,data=std_619[train_data_indices, ],hidden=c(3,2),linear.output=T)
pred_nn <- predict(nn, newdata=std_619[-train_data_indices, ])

pred_nn_res <- pred_nn*(max(std_619$median)-min(std_619$median))+min(std_619$median)
test_res <- (std_619[-train_data_indices, ]$median)*(max(std_619$median)-min(std_619$median))+min(std_619$median)
MSE_nn <- sum((test_res - pred_nn_res)^2)/nrow(std_619[-train_data_indices, ])

# compare to SAR

w <- knn2nb(knearneigh(coordinates(plots_619[train_data_indices, ]), k=8))

variables <- "`N`+`Ca`+`S`+`Cu`+`B`+`Al`+`soil_GoA`+`soil_Ra`+`soil_MaB`"
dat <- std_619[train_data_indices, ]

metric <- "median"
null <- errorsarlm( as.formula(paste(metric, 1, sep="~")), data = dat, 
                  listw=nb2listw(w), zero.policy = T, na.action = na.fail)
sar <- errorsarlm(as.formula(paste(metric, variables, sep = " ~ ")), 
                data = dat, listw=nb2listw(w), zero.policy = T, na.action = na.fail, 
               )
summary(sar)

predict_sar <- predict.sarlm(sar,newdata=std_619[-train_data_indices, ] )
MSE_sar <- sum((predict_sar - std_619[-train_data_indices, ]$median)^2)/nrow(std_619[-train_data_indices, ])


# ---------------------------------------- Try RF regression --------------------------------------------------

#all nutrients: `N` + `P` + `K`+`Mg` +`Ca` + `S` + `Zn`+`Mn`+`Cu`+`Fe`+`B`+`Al`+`Na`

#--------------------------------RF 619 ---------------------------------

set.seed(42)
train_data_indices <- rep(FALSE, nrow(std_619))
train_data_indices[sample(1:nrow(std_619), round(0.7 * nrow(std_619)))] <- TRUE # randomly select 70% of the data for training
rf_regression<- randomForest(median ~ `N`+`P`+`K`+`Mg`+`Ca`+`S`+`Zn`+`Mn`+`Cu`+`Fe`+`B`+`Al`+`Na`+`soil_GoA`+`soil_Ra`+`soil_MaB`+`watersum`, data = std_619[train_data_indices, ], importance=T)
rf_regression
varImpPlot(rf_regression)
pred_struc <- predict(rf_regression, std_619[!train_data_indices,]) # predict the rings
#table(pred_struc, std_619$median[!train_data_indices])
plot(std_619$skew[!train_data_indices], pred_struc, xlab="Observed", ylab="Predicted")

#--------------------------------RF 703 ---------------------------------

set.seed(42)
train_data_indices <- rep(FALSE, nrow(std_703))
train_data_indices[sample(1:nrow(std_703), round(0.7 * nrow(std_703)))] <- TRUE # randomly select 70% of the data for training
rf_regression<- randomForest(skew ~ `N` + `P` + `K`+`Mg` +`Ca` + `Zn`+`Cu`+`Fe`+`B`+`Al`+`Na`, data = std_703[train_data_indices, ], importance=T)
rf_regression
varImpPlot(rf_regression)
pred_struc <- predict(rf_regression, std_703[!train_data_indices,]) # predict the rings
#table(pred_struc, std_703$median[!train_data_indices])
plot(std_703$skew[!train_data_indices], pred_struc, xlab="Observed", ylab="Predicted")

#--------------------------------RF 717 ---------------------------------

set.seed(42)
train_data_indices <- rep(FALSE, nrow(std_703))
train_data_indices[sample(1:nrow(std_717), round(0.7 * nrow(std_717)))] <- TRUE # randomly select 70% of the data for training
rf_regression<- randomForest(median ~ `watersum`+`N` + `P` + `K`+`Mg` +`Ca` + `Zn`+`Cu`+`Fe`+`B`+`Al`+`Na`, data = std_717[train_data_indices, ], importance=T)
rf_regression
varImpPlot(rf_regression)
pred_struc <- predict(rf_regression, std_717[!train_data_indices,]) # predict the rings
#table(pred_struc, std_717$median[!train_data_indices])
plot(std_717$median[!train_data_indices], pred_struc, xlab="Observed", ylab="Predicted")


# ---------------------------------------Try multinominal logisitic regression with nnet package --------------------------------

# 619

library(nnet)
plots_619_df$Trt_agg <- relevel(plots_619_df$Trt_agg, ref = "100% Control")
mltnom_619 <- multinom(Trt_agg ~ crrmedian+median+kurt+rumple+moran+dip+crriqr, data = plots_619_df, na.action = na.fail)
summary(mltnom_619)

modsel_619 <- dredge(mltnom_619)
modsel_619

topmod_619 <- multinom(Trt_agg ~ crrmedian+median+kurt+rumple+moran, data = plots_619_df, na.action = na.fail)
summary(topmod_619)

library(DescTools)
PseudoR2(topmod_619, which = "all")
z <- summary(topmod_619)$coefficients/summary(topmod_619)$standard.errors
z
p <- (1 - pnorm(abs(z), 0, 1)) * 2
p

coefs_619 <- coef(topmod_619)
exp(coefs_619) #transformed to odds ratio
(exp(coefs_619)-1)*100 # percent change in the odds for a one unit increase in the independent variable

# Try binary logistic regression with 100% control as baseline
# def_100 <- as.data.frame(dt_619[Trt_agg == "100% Control" | Trt_agg == "deficiency",])
# def_100 <- droplevels(def_100)
# tox_100 <- as.data.frame(dt_619[Trt_agg == "100% Control" | Trt_agg == "toxicity",])
# tox_100 <- droplevels(tox_100)
# control_0_100 <- as.data.frame(dt_619[Trt_agg == "100% Control" | Trt_agg == "0% Control",])
# control_0_100 <- droplevels(control_0_100)
# 
# binary_619 <- glm(Trt_agg ~ crrmedian+median+kurt+rumple+moran+crriqr+skew+iqr, data = control_0_100, family = binomial, na.action = na.fail)
# modsel_619_binary <- dredge(binary_619)
# topmod_619_binary <- glm(Trt_agg ~ crrmedian+iqr+moran+rumple+skew, data = tox_100, family = binomial)
# summary(topmod_619_binary)
# 
# library(pscl)
# 
# pR2(topmod_619_binary)

# Try multinomial logisitc regression with mlogit package
# library(foreign)
# library(mlogit)
# 
# mlogit_data_619 <- mlogit.data(plots_619_df, choice = "Trt_agg", shape = "wide")
# mlogit_619 <- mlogit(Trt_agg~crrmedian+kurt+median+skew_ln+rumple+moran, data = mlogit_data_619,
#                      method = "nr", print.level = 0)
# summary(mlogit_619)
# 
# # Weights matrix
# w <- knn2nb(knearneigh(coordinates(plots), k=8))


# 703
plots_703_df$Trt_agg <- relevel(plots_703_df$Trt_agg, ref = "100% Control")
mltnom_703 <- multinom(Trt_agg ~ crrmedian+median+kurt+rumple+moran+dip+crriqr, data = plots_703_df, na.action = na.fail)
summary(mltnom_703)

modsel_703 <- dredge(mltnom_703)
modsel_703

topmod_703 <- multinom(Trt_agg ~ crrmedian+kurt+median+moran, data = plots_703_df, na.action = na.fail)
summary(topmod_703)
PseudoR2(topmod_703, which = "all")
z <- summary(topmod_703)$coefficients/summary(topmod_703)$standard.errors
z
p <- (1 - pnorm(abs(z), 0, 1)) * 2
p

coefs_703 <- coef(topmod_703)
exp(coefs_703) #transformed to odds ratio
(exp(coefs_703)-1)*100 # percent change in the odds for a one unit increase in the independent variable

# 717
plots_717_df$Trt_agg <- relevel(plots_717_df$Trt_agg, ref = "100% Control")
mltnom_717 <- multinom(Trt_agg ~ crrmedian+median+kurt+rumple+moran+dip+crriqr, data = plots_717_df, na.action = na.fail)
summary(mltnom_717)

modsel_717 <- dredge(mltnom_717)
modsel_717

topmod_717 <- multinom(Trt_agg ~ median+rumple+moran, data = plots_717_df, na.action = na.fail)
summary(topmod_717)
PseudoR2(topmod_717, which = "all")
z <- summary(topmod_717)$coefficients/summary(topmod_717)$standard.errors
z
p <- (1 - pnorm(abs(z), 0, 1)) * 2
p

coefs_717 <- coef(topmod_717)
exp(coefs_717) #transformed to odds ratio
(exp(coefs_717)-1)*100 # percent change in the odds for a one unit increase in the independent variable


# Try combining tox and def to see if model improves
plots_717_df$Trt_agg_v2 <- as.character(plots_619_df$Trt_agg)
plots_717_df[Trt_agg == "deficiency"]$Trt_agg_v2 <- "stressed"
plots_717_df[Trt_agg == "toxicity"]$Trt_agg_v2 <- "stressed"
plots_717_df$Trt_agg_v2 <- as.factor(plots_717_df$Trt_agg_v2) 

plots_717_df$Trt_agg <- relevel(plots_717_df$Trt_agg, ref = "100% Control")

mltnom_717_v2 <- multinom(Trt_agg_v2 ~ crrmedian+median+kurt+rumple+moran+dip+crriqr, data = plots_717_df, na.action = na.fail)
summary(mltnom_717_v2)

modsel_717_v2 <- dredge(mltnom_717_v2)
modsel_717_v2

topmod_717_v2 <- multinom(Trt_agg_v2 ~ median+rumple+moran+crrmedian, data = plots_717_df, na.action = na.fail)
summary(topmod_717_v2)
PseudoR2(topmod_717_v2, which = "all")
z <- summary(topmod_717)$coefficients/summary(topmod_717)$standard.errors
z
p <- (1 - pnorm(abs(z), 0, 1)) * 2
p


# ----------------------------------------------Result tables------------------------------------------------------

vars <- topmod_619$coefnames[-1]
sd_all <- numeric(length(vars))
mean_all <- numeric(length(vars))
count <- 1
for (x in vars){
  values <- as.data.frame(plots_619_df)[x]
  sd_all[count] <- sd(values[,1])
  mean_all[count] <- mean(values[,1])
  count <- count+1
}

vars_619 <- as.data.frame(vars)
vars_619$mean <- mean_all
vars_619$sd <- sd_all

#write.csv(vars_619, "Q:/My Drive/Research/Canopy_Morphology/Tobacco_Project/canopy_analysis/tobacco_results/vars_619.csv")

levels <- topmod_619$lev[-1]
#df_all <- data.frame(var=character(), coefs=list(), SE=double(), Wald=double(), p=double(), exp=double())
coefs <- numeric(length(levels))
std_coefs <- numeric(length(levels))
SE <- numeric(length(levels))
Wald <- numeric(length(levels))
p <- numeric(length(levels))
exp <- numeric(length(levels))
count <- 1

for (level in levels){
  coefs[count] <- list(coef(topmod_619)[count,])
  SE[count] <- list(summary(topmod_619)$standard.errors[count,])
  z <- (summary(topmod_619)$coefficients/summary(topmod_619)$standard.errors)
  Wald[count] <- list(z[count,])
  p[count] <- list(((1 - pnorm(abs(z), 0, 1)) * 2)[count,])
  std_coefs[count] <- list((sqrt(3)/pi) * (coefs[[count]][-1] * sd_all))
  exp[count] <- list(exp(std_coefs[[count]]))
  count <- count+1
}

vars <- topmod_619$coefnames

results_0_619 <- as.data.frame(vars)
results_0_619$level <- "0% Control"
results_0_619$coefs <- coefs[[1]]
results_0_619$std_coefs <- c(0,std_coefs[[1]])
results_0_619$exp <- c(0,exp[[1]])
results_0_619$SE <- SE[[1]]
results_0_619$p <- p[[1]]
results_0_619$wald <- Wald[[1]]

results_tox_619 <- as.data.frame(vars)
results_tox_619$level <- "toxicity"
results_tox_619$coefs <- coefs[[3]]
results_tox_619$std_coefs <- c(0,std_coefs[[3]])
results_tox_619$exp <- c(0,exp[[3]])
results_tox_619$SE <- SE[[3]]
results_tox_619$p <- p[[3]]
results_tox_619$wald <- Wald[[3]]

results_def_619 <- as.data.frame(vars)
results_def_619$level <- "deficiency"
results_def_619$coefs <- coefs[[2]]
results_def_619$std_coefs <- c(0,std_coefs[[2]])
results_def_619$exp <- c(0,exp[[2]])
results_def_619$SE <- SE[[2]]
results_def_619$p <- p[[2]]
results_def_619$wald <- Wald[[2]]



results_619 <- rbind(results_0_619, results_def_619, results_tox_619)


# mlogit_results <- function(x){
#   levels <- x$lev[-1]
#   coefs <- numeric(length(levels))
#   SE <- numeric(length(levels))
#   Wald <- numeric(length(levels))
#   p <- numeric(length(levels))
#   exp <- numeric(length(levels))
#   count <- 1
#   for (level in levels){
#     coefs[count] <- list(coef(x)[count,])
#     SE[count] <- list(summary(x)$standard.errors[count,])
#     Wald[count] <- list(summary(x)$coefficients/summary(topmod_619)$standard.errors[count,])
#     p[count] <- list(((1 - pnorm(abs(z), 0, 1)) * 2)[count,])
#     exp[count] <- list(exp(coefs[[count]]))
#     count <- count+1
#   }
#   levels
#   coefs
#   SE
#   Wald
#   p
#   exp
# }
# 
# mlogit_results_619 <- mlogit_results(topmod_619)


#write.csv(results_619, "Q:/My Drive/Research/Canopy_Morphology/Tobacco_Project/canopy_analysis/tobacco_results/results_619.csv")


vars <- topmod_703$coefnames[-1]
sd_all <- numeric(length(vars))
mean_all <- numeric(length(vars))
count <- 1
for (x in vars){
  values <- as.data.frame(plots_703_df)[x]
  sd_all[count] <- sd(values[,1])
  mean_all[count] <- mean(values[,1])
  count <- count+1
}

vars_703 <- as.data.frame(vars)
vars_703$mean <- mean_all
vars_703$sd <- sd_all

#write.csv(vars_703, "Q:/My Drive/Research/Canopy_Morphology/Tobacco_Project/canopy_analysis/tobacco_results/vars_703.csv")

levels <- topmod_703$lev[-1]
#df_all <- data.frame(var=character(), coefs=list(), SE=double(), Wald=double(), p=double(), exp=double())
coefs <- numeric(length(levels))
std_coefs <- numeric(length(levels))
SE <- numeric(length(levels))
Wald <- numeric(length(levels))
p <- numeric(length(levels))
exp <- numeric(length(levels))
count <- 1

for (level in levels){
  coefs[count] <- list(coef(topmod_703)[count,])
  SE[count] <- list(summary(topmod_703)$standard.errors[count,])
  z <- (summary(topmod_703)$coefficients/summary(topmod_703)$standard.errors)
  Wald[count] <- list(z[count,])
  p[count] <- list(((1 - pnorm(abs(z), 0, 1)) * 2)[count,])
  std_coefs[count] <- list((sqrt(3)/pi) * (coefs[[count]][-1] * sd_all))
  exp[count] <- list(exp(std_coefs[[count]]))
  count <- count+1
}

vars <- topmod_703$coefnames

results_0_703 <- as.data.frame(vars)
results_0_703$level <- "0% Control"
results_0_703$coefs <- coefs[[1]]
results_0_703$std_coefs <- c(0,std_coefs[[1]])
results_0_703$exp <- c(0,exp[[1]])
results_0_703$SE <- SE[[1]]
results_0_703$p <- p[[1]]
results_0_703$wald <- Wald[[1]]

results_tox_703 <- as.data.frame(vars)
results_tox_703$level <- "toxicity"
results_tox_703$coefs <- coefs[[3]]
results_tox_703$std_coefs <- c(0,std_coefs[[3]])
results_tox_703$exp <- c(0,exp[[3]])
results_tox_703$SE <- SE[[3]]
results_tox_703$p <- p[[3]]
results_tox_703$wald <- Wald[[3]]

results_def_703 <- as.data.frame(vars)
results_def_703$level <- "deficiency"
results_def_703$coefs <- coefs[[2]]
results_def_703$std_coefs <- c(0,std_coefs[[2]])
results_def_703$exp <- c(0,exp[[2]])
results_def_703$SE <- SE[[2]]
results_def_703$p <- p[[2]]
results_def_703$wald <- Wald[[2]]


results_703 <- rbind(results_0_703, results_def_703, results_tox_703)

#write.csv(results_703, "Q:/My Drive/Research/Canopy_Morphology/Tobacco_Project/canopy_analysis/tobacco_results/results_703.csv")


vars <- topmod_717$coefnames[-1]
sd_all <- numeric(length(vars))
mean_all <- numeric(length(vars))
count <- 1
for (x in vars){
  values <- as.data.frame(plots_717_df)[x]
  sd_all[count] <- sd(values[,1])
  mean_all[count] <- mean(values[,1])
  count <- count+1
}

vars_717 <- as.data.frame(vars)
vars_717$mean <- mean_all
vars_717$sd <- sd_all

#write.csv(vars_717, "Q:/My Drive/Research/Canopy_Morphology/Tobacco_Project/canopy_analysis/tobacco_results/vars_717.csv")

levels <- topmod_717$lev[-1]
#df_all <- data.frame(var=character(), coefs=list(), SE=double(), Wald=double(), p=double(), exp=double())
coefs <- numeric(length(levels))
std_coefs <- numeric(length(levels))
SE <- numeric(length(levels))
Wald <- numeric(length(levels))
p <- numeric(length(levels))
exp <- numeric(length(levels))
count <- 1

for (level in levels){
  coefs[count] <- list(coef(topmod_717)[count,])
  SE[count] <- list(summary(topmod_717)$standard.errors[count,])
  z <- (summary(topmod_717)$coefficients/summary(topmod_717)$standard.errors)
  Wald[count] <- list(z[count,])
  p[count] <- list(((1 - pnorm(abs(z), 0, 1)) * 2)[count,])
  std_coefs[count] <- list((sqrt(3)/pi) * (coefs[[count]][-1] * sd_all))
  exp[count] <- list(exp(std_coefs[[count]]))
  count <- count+1
}

vars <- topmod_717$coefnames

results_0_717 <- as.data.frame(vars)
results_0_717$level <- "0% Control"
results_0_717$coefs <- coefs[[1]]
results_0_717$std_coefs <- c(0,std_coefs[[1]])
results_0_717$exp <- c(0,exp[[1]])
results_0_717$SE <- SE[[1]]
results_0_717$p <- p[[1]]
results_0_717$wald <- Wald[[1]]

results_tox_717 <- as.data.frame(vars)
results_tox_717$level <- "toxicity"
results_tox_717$coefs <- coefs[[3]]
results_tox_717$std_coefs <- c(0,std_coefs[[3]])
results_tox_717$exp <- c(0,exp[[3]])
results_tox_717$SE <- SE[[3]]
results_tox_717$p <- p[[3]]
results_tox_717$wald <- Wald[[3]]

results_def_717 <- as.data.frame(vars)
results_def_717$level <- "deficiency"
results_def_717$coefs <- coefs[[2]]
results_def_717$std_coefs <- c(0,std_coefs[[2]])
results_def_717$exp <- c(0,exp[[2]])
results_def_717$SE <- SE[[2]]
results_def_717$p <- p[[2]]
results_def_717$wald <- Wald[[2]]

results_717 <- rbind(results_0_717, results_def_717, results_tox_717)

#write.csv(results_717, "Q:/My Drive/Research/Canopy_Morphology/Tobacco_Project/canopy_analysis/tobacco_results/results_717.csv")


# ------------------------------------- Random Forest Classification ------------------------------------------------------
library(randomForest)
set.seed(42)
train_data_indices <- rep(FALSE, nrow(norm_1cm))
train_data_indices[sample(1:nrow(norm_1cm), round(0.7 * nrow(norm_1cm)))] <- TRUE # randomly select 80% of the data for training
rf_regression_1cm<- randomForest(BuAc_ln ~ median+skew+PlotCRR_ln+crrmedian_ln+crriqr+Trt+AUDPC_ln+geary+iqr_ln, data=norm_1cm[train_data_indices, ], importance=T)
rf_regression_1cm
varImpPlot(rf_regression_1cm)
pred_trt <- predict(rf_regression_1cm, norm_1cm[!train_data_indices,]) # predict the rings
table(pred_trt, norm_1cm$BuAc_ln[!train_data_indices])
plot(norm_1cm$BuAc_ln[!train_data_indices], pred_trt, xlab="Observed", ylab="Predicted")
#abline(a=0, b=1, lty=2, col=2)

set.seed(42)
train_data_indices <- rep(FALSE, nrow(norm_5cm))
train_data_indices[sample(1:nrow(norm_5cm), round(0.7 * nrow(norm_5cm)))] <- TRUE # randomly select 80% of the data for training
rf_regression_5cm<- randomForest(BuAc_ln ~ median+skew+PlotCRR_ln+crrmedian_ln+crriqr+Trt+AUDPC_ln+geary+iqr_ln, data=norm_5cm[train_data_indices, ], importance=T)
rf_regression_5cm
varImpPlot(rf_regression_5cm)
pred_trt <- predict(rf_regression_5cm, norm_5cm[!train_data_indices,]) # predict the rings
table(pred_trt, norm_5cm$BuAc_ln[!train_data_indices])
plot(norm_5cm$BuAc_ln[!train_data_indices], pred_trt, xlab="Observed", ylab="Predicted")
#abline(a=0, b=1, lty=2, col=2)

set.seed(42)
train_data_indices <- rep(FALSE, nrow(norm_10cm))
train_data_indices[sample(1:nrow(norm_10cm), round(0.7 * nrow(norm_10cm)))] <- TRUE # randomly select 80% of the data for training
rf_regression_10cm<- randomForest(BuAc_ln ~ median+skew+PlotCRR_ln+crrmedian_ln+crriqr+Trt+AUDPC_ln+geary+iqr_ln, data=norm_10cm[train_data_indices, ], importance=T)
rf_regression_10cm
varImpPlot(rf_regression_10cm)
pred_trt <- predict(rf_regression_10cm, norm_10cm[!train_data_indices,]) # predict the rings
table(pred_trt, norm_10cm$BuAc_ln[!train_data_indices])
plot(norm_10cm$BuAc_ln[!train_data_indices], pred_trt, xlab="Observed", ylab="Predicted")
#abline(a=0, b=1, lty=2, col=2)

set.seed(42)
train_data_indices <- rep(FALSE, nrow(norm_20cm))
train_data_indices[sample(1:nrow(norm_20cm), round(0.7 * nrow(norm_20cm)))] <- TRUE # randomly select 80% of the data for training
rf_regression_20cm<- randomForest(BuAc_ln ~ median+skew+PlotCRR_ln+crrmedian_ln+crriqr+Trt+AUDPC_ln+geary+iqr_ln, data=norm_20cm[train_data_indices, ], importance=T)
rf_regression_20cm
varImpPlot(rf_regression_20cm)
pred_trt <- predict(rf_regression_20cm, norm_20cm[!train_data_indices,]) # predict the rings
table(pred_trt, norm_20cm$BuAc_ln[!train_data_indices])
plot(norm_20cm$BuAc_ln[!train_data_indices], pred_trt, xlab="Observed", ylab="Predicted")
#abline(a=0, b=1, lty=2, col=2)

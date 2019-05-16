library(raster)
library(rgdal)
library(data.table)
library(plyr)
library(lattice)
library(moments)
library(gridExtra)
library(spdep)
library(velox)
library(randomForest)
library(rgeos)
library(gstat)
library(usdm)
library(viridis)


#Load raster data for all dates
setwd("/media/Kellyn/F20E17B40E177139/kpmontgo@ncsu.edu/LkWheeler_Sorghum/LkWheeler_Fusarium_Sorghum")
#setwd("Q:/My Drive/LkWheeler_Sorghum/LkWheeler_Fusarium_Sorghum")

sorghum_area <- readOGR("sorghum_area", "sorghum_area", stringsAsFactors = F)

csm_1cm <- crop(raster('CSM/csm_adj_9_1_noDepth_1cm.tif'), sorghum_area)
csm_5cm <- crop(raster('CSM/csm_adj_9_1_noDepth_5cm.tif'), sorghum_area)
csm_10cm <- crop(raster('CSM/csm_adj_9_1_noDepth_10cm.tif'), sorghum_area)
csm_20cm <- crop(raster('CSM/csm_adj_9_1_noDepth_20cm.tif'), sorghum_area)

csm_1cm_velox <- velox(csm_1cm)
csm_5cm_velox <- velox(csm_5cm)
csm_10cm_velox <- velox(csm_10cm)
csm_20cm_velox <- velox(csm_20cm)

#ortho_9_1 <- crop(raster("orthos/ortho_9_1_18.tif"))

#------------------------------------------------------- Field data-------------------------------------------------------------------
field_data <- fread("Sorghum 2018 JB Yield AUDPC.csv")
field_data$Plot <- as.integer(field_data$Plot)

# If normalizing:
# field_data_norm <- field_data[,3]
# field_data_norm$BuAc <- field_data$BuAc
# field_data_norm$AUDPC <- field_data$AUDPC
# field_data_norm <- as.data.frame(scale(field_data_norm))
# field_data_norm$Trt <- field_data$Trt
# field_data_norm$Plot <- field_data$Plot

plots <- readOGR("sorghum_plots", "sorghum_plots_small", stringsAsFactors = F)
plots <- spTransform(plots, CRS("+proj=lcc +lat_1=36.16666666666666 +lat_2=34.33333333333334 +lat_0=33.75 +lon_0=-79 +x_0=609601.22 +y_0=0 +datum=NAD83 +units=m +no_defs +ellps=GRS80 +towgs84=0,0,0"))
plots$id <- as.integer(plots$id)

plots_1cm <- plots
plots_5cm <- plots
plots_10cm <- plots
plots_20cm <- plots

plots_1cm <- merge(plots_1cm, field_data, by.x="id", by.y="Plot")
plots_5cm <- merge(plots_5cm, field_data, by.x="id", by.y="Plot")
plots_10cm <- merge(plots_10cm, field_data, by.x="id", by.y="Plot")
plots_20cm <- merge(plots_20cm, field_data, by.x="id", by.y="Plot")

#--------------------------------------------------------- SIMWE water depth model---------------------------------------------------------------
water <- crop(raster("simwe/water_depth.tif"), sorghum_area)

water_velox <- velox(water)
water_extract <- water_velox$extract(sp=plots)
names(water_extract) <- plots$id

watermedian <- sapply(water_extract, median, na.rm=TRUE)
watermean <- sapply(water_extract, mean, na.rm=TRUE)
watersd <- sapply(water_extract, sd, na.rm=TRUE)
watersum <- sapply(water_extract, sum, na.rm=TRUE)

waterDF <- data.frame(plots=plots$id, watersum = watersum, watersd=watersd, watermean=watermean)
# norm_water <- as.data.frame(scale(waterDF[,2:4]))
# norm_water$plots <- waterDF$plots

plots_1cm <- merge(plots_1cm, waterDF, by.x="id", by.y="plots")
plots_5cm <- merge(plots_5cm, waterDF, by.x="id", by.y="plots")
plots_10cm <- merge(plots_10cm, waterDF, by.x="id", by.y="plots")
plots_20cm <- merge(plots_20cm, waterDF, by.x="id", by.y="plots")

#---------------------------------------------------------- Canopy relief ratio----------------------------------------------------------------
crr <- function(x){
  (mean(x, na.rm=T)-min(x, na.rm=T))/(max(x, na.rm=T)-min(x, na.rm=T))
}

crrRast_1cm <- focal(csm_1cm, w=matrix(1/2401, nc=49, nr=49), crr)
crrRast_5cm <- focal(csm_5cm, w=matrix(1/81, nc=9, nr=9), crr)
crrRast_10cm <- focal(csm_10cm, w=matrix(1/25, nc=5, nr=5), crr)
crrRast_20cm <- focal(csm_20cm, w=matrix(1/9, nc=3, nr=3), crr)
crrRast_list <- list(crrRast_1cm,crrRast_5cm,crrRast_10cm,crrRast_20cm)

crr_metrics <- function(x){
  veloxRast <- velox(x)
  extract <- veloxRast$extract(sp=plots)
  names(extract) <- plots$id
  crrmean <- sapply(extract, mean, na.rm=T)
  crrmedian <- sapply(extract, median, na.rm=T)
  crrsd <- sapply(extract, sd, na.rm=T)
  crrDF <- data.frame(plots=plots$id, crrmedian= crrmedian, crrmean=crrmean, crrsd = crrsd)
  crrDF
  # Normalize data
  # norm_crr <- as.data.frame(scale(crrDF[,2:4]))
  # norm_crr$plots <- crrDF$plots
  # norm_crr
}

crr_df_list <- lapply(crrRast_list, crr_metrics)

plots_1cm <- merge(plots_1cm, crr_df_list[[1]], by.x="id", by.y="plots")
plots_5cm <- merge(plots_5cm, crr_df_list[[2]], by.x="id", by.y="plots")
plots_10cm <- merge(plots_10cm, crr_df_list[[3]], by.x="id", by.y="plots")
plots_20cm <- merge(plots_20cm, crr_df_list[[4]], by.x="id", by.y="plots")


#---------------------------------------------------------- Crop height ----------------------------------------------------------
plot_extract_1cm <- csm_1cm_velox$extract(sp=plots)
plot_extract_5cm <- csm_5cm_velox$extract(sp=plots)
plot_extract_10cm <- csm_10cm_velox$extract(sp=plots)
plot_extract_20cm <- csm_20cm_velox$extract(sp=plots)
names(plot_extract_1cm) <- plots$id
names(plot_extract_5cm) <- plots$id
names(plot_extract_10cm) <- plots$id
names(plot_extract_20cm) <- plots$id

extract_list <- list(plot_extract_1cm, plot_extract_5cm, plot_extract_10cm, plot_extract_20cm)

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
  cropHeightDF <- data.frame(plots=plots$id, median= CHmedian, mean=CHmean, sd = CHsd, iqr = CHiqr, var = CHvar, sum = CHsum, skew = CHskew, kurt= CHkurt, PlotCRR = CHcrr, chRange = CHrange)
  cropHeightDF
  # Normalize data
  # norm_ch <- as.data.frame(scale(cropHeightDF[,2:10]))
  # norm_ch$plots <- cropHeightDF$plots
  # norm_ch
}

ch_df_list <- lapply(extract_list, ch_metrics)

#merge all metrics into spatial polygons
plots_1cm <- merge(plots_1cm, ch_df_list[[1]], by.x="id", by.y="plots")
plots_5cm <- merge(plots_5cm, ch_df_list[[2]], by.x="id", by.y="plots")
plots_10cm <- merge(plots_10cm, ch_df_list[[3]], by.x="id", by.y="plots")
plots_20cm <- merge(plots_20cm, ch_df_list[[4]], by.x="id", by.y="plots")

# 10th and 90th percentile yield
pct_10_90_df <- function(x){
  top_plots <- c('410', '401','407', '406', '408')
  top_yield <- x[top_plots]
  top_yield_all <- do.call(c, top_yield)
  top_yield_all <- as.data.frame(top_yield_all)
  names(top_yield_all) <- "height"
  top_yield_all$group <- "90th percentile yield"
  
  bottom_plots <- c('203', '202','205', '204', '212')
  bottom_yield <- x[bottom_plots]
  bottom_yield_all <- do.call(c, bottom_yield)
  bottom_yield_all <- as.data.frame(bottom_yield_all)
  names(bottom_yield_all) <- "height"
  bottom_yield_all$group <- "10th percentile yield"
  
  pct_10_90_yield <- rbind(top_yield_all, bottom_yield_all)
}

pct_10_90_df_list <- lapply(extract_list, pct_10_90_df)

# Histogram and denisty of top vs bottom yield
library(ggplot2)

hist_10_90 <- function(x){
  height_mean <- ddply(x, "group", summarise, height.mean=mean(height))
  height_median <- ddply(x, "group", summarise, height.median=median(height))
  ggplot(x, aes(x=height, fill=group, color=group)) +
    geom_histogram(binwidth=.03, position="identity", alpha=0.7) +
    geom_vline(data=height_median, aes(xintercept=height.median,  colour=group), linetype="solid", size=1) +
    scale_color_grey()+scale_fill_grey() +
    theme_classic() + labs(x="Crop Height (m)", y = "Count") + theme(legend.position= c(0.2,0.8), legend.spacing.x = unit(0.3, 'cm'), 
                                                                     legend.spacing.y = unit(0.2, 'cm'), legend.text = element_text(size=18),
                                                                     legend.title = element_blank(),
                                                                     axis.title.x = element_text(size=18, face="bold"),
                                                                     axis.text=element_text(size=14),
                                                                     axis.title.y = element_text(size=18, face="bold"))
}

lapply(pct_10_90_df_list, hist_10_90)


#--------------------------------------------------------- rumple index -----------------------------------------------------------------
library(lidR)

rumple <- function(x){
  plotNums <- plots$id
  rumple_all <- numeric(length(plotNums))
  names(rumple_all) <- "rumple"
  count <- 1
  for (num in plotNums){
    tmp <- plots@data
    position <- match(num, tmp$id)
    extract <- crop(x, plots[position,])
    rumple <- rumple_index(extract)
    rumple_all[count] <- rumple
    count <- count+1
  }
  rumple_df <- data.frame(plots=plots$id, rumple= rumple_all)
  rumple_df
  #Normalize data
  # norm_rumple <- as.data.frame(scale(rumple_df[,2]))
  # names(norm_rumple) <- "rumple"
  # norm_rumple$plots <- rumple_df$plots
  # norm_rumple
}

csm_list <- list(csm_1cm, csm_5cm, csm_10cm, csm_20cm)
rumple_df_list <- lapply(csm_list, rumple)

#merge all metrics into spatial polygons
plots_1cm <- merge(plots_1cm, rumple_df_list[[1]], by.x="id", by.y="plots")
plots_5cm <- merge(plots_5cm, rumple_df_list[[2]], by.x="id", by.y="plots")
plots_10cm <- merge(plots_10cm, rumple_df_list[[3]], by.x="id", by.y="plots")
plots_20cm <- merge(plots_20cm, rumple_df_list[[4]], by.x="id", by.y="plots")

# -------------------------------------------------------- variogram models ------------------------------------------------------------
# 
# variogram_metrics <- function(x){
#   plotNums <- plots$id
#   range_all <- numeric(length(plotNums))
#   sill_all <- numeric(length(plotNums))
#   count <- 1
#   
#   for (num in plotNums){
#     tmp <- plots@data
#     position <- match(num, tmp$id)
#     extract <- crop(x, plots[position,])
#     samp <- sampleRandom(extract, 200, sp=T, round(2)) #This sample number needs to be varied by resolution most likely
#     colnames(samp@data)[1] <- "height"
#     variogram <- variogram(samp$height~1, locations=samp)
#     fit_vgm <- fit.variogram(variogram, vgm(c("Exp", "Sph", "Mat", "Gau")))
#     range <- round(fit_vgm$range[2], 4)
#     sill <- round(fit_vgm$psill[2], 4)
#     range_all[count] <- range
#     sill_all[count] <- sill
#     count <- count+1
#   }
#   
#   names(sill_all) <- plots$id
#   names(range_all) <- plots$id
#   var_df <- data.frame(plots=plots$id, sill= sill_all, range=range_all)
#   # 
#   # Normalize data
#   norm_var <- as.data.frame(scale(var_df[,2:3]))
#   norm_var$plots <- var_df$plots
#   norm_var
# }
# 

# var_df_list <- lapply(csm_list, variogram_metrics)
# 
# plots_9_1 <- merge(plots_9_1, norm_var, by.x="id", by.y="plots")


# ------------------------------------------------ spatial autocorrelation -----------------------------------------------------

autocor_metrics <- function(x, w){
  plotNums <- plots$id
  morans_all <- numeric(length(plotNums))
  geary_all <- numeric(length(plotNums))
  count <- 1
  
  for (num in plotNums){
    tmp <- plots@data
    position <- match(num, tmp$id)
    extract <- crop(x, plots[position,])
    moran <- Moran(extract, w)
    morans_all[count] <- moran
    geary <- Geary(extract, w)
    geary_all[count] <- geary
    count <- count+1
  }
  
  names(morans_all) <- plots$id
  names(geary_all) <- plots$id
  autocor_df <- data.frame(plots=plots$id, moran= morans_all, geary=geary_all)
  autocor_df
}

#Using approx 30 cm x 30 cm weights matrix
autocor_1cm <- autocor_metrics(csm_1cm, matrix(c(1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,0,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1), nc=25, nr=25))
autocor_5cm <- autocor_metrics(csm_5cm, matrix(c(1,1,1,1,1,1,1,1,1,1,1,1,0,1,1,1,1,1,1,1,1,1,1,1,1), nc=5, nr=5))
autocor_10cm <- autocor_metrics(csm_10cm, matrix(c(1,1,1,1,0,1,1,1,1), nc=3, nr=3))
autocor_20cm <- autocor_metrics(csm_20cm, matrix(c(1,1,1,1,0,1,1,1,1), nc=3, nr=3))

plots_1cm <- merge(plots_1cm, autocor_1cm, by.x="id", by.y="plots")
plots_5cm <- merge(plots_5cm, autocor_5cm, by.x="id", by.y="plots")
plots_10cm <- merge(plots_10cm, autocor_10cm, by.x="id", by.y="plots")
plots_20cm <- merge(plots_20cm, autocor_20cm, by.x="id", by.y="plots")


# ------------------------------------Plot variables ------------------------------------------

imp_vars <- c("mean", "median", "skew", "iqr", "kurt", "sd", "crrsd", "crrmean", "moran", "Trt", "AUDPC", "watersum", "PlotCRR", "rumple") 

pairs(BuAc~median+skew+crrsd+crrmean+kurt+Trt+AUDPC+rumple+geary+iqr,data=plots_5cm, lower.panel = NULL, 
      main="Simple Scatterplot Matrix")

par(mfrow=c(2,4), mar=c(5, 6, 2, 3))

#plot(plots_5cm$mean, plots_5cm$BuAc)
plot(plots_5cm$median, plots_5cm$BuAc, xlab = "Median Height (m)", ylab = "Yield (bu/ac)", cex = 1,cex.lab = 2, cex.axis = 2, pch=16)
#plot(plots_5cm$sum, plots_5cm$BuAc)
plot(plots_5cm$iqr, plots_5cm$BuAc, xlab = "Interquartile Range (m)", ylab = "Yield (bu/ac)", cex = 1,cex.lab = 2, cex.axis = 2, pch=16)
plot(plots_5cm$rumple, plots_5cm$BuAc, xlab = "Rumple Index", ylab = "Yield (bu/ac)", cex = 1,cex.lab = 2, cex.axis = 2, pch=16)
plot(plots_5cm$skew, plots_5cm$BuAc, xlab = "Skewness", ylab = "Yield (bu/ac)", cex = 1,cex.lab = 2, cex.axis = 2, pch=16)
plot(plots_5cm$kurt, plots_5cm$BuAc, xlab = "Kurtosis", ylab = "Yield (bu/ac)", cex = 1,cex.lab = 2, cex.axis = 2, pch=16)
#plot(plots_5cm$sd, plots_5cm$BuAc)
plot(plots_5cm$crrmean, plots_5cm$BuAc, xlab = "Canopy Relief Ratio Mean", ylab = "Yield (bu/ac)", cex = 1,cex.lab = 2, cex.axis = 2, pch=16)
plot(plots_5cm$crrsd, plots_5cm$BuAc, xlab = "Canopy Relief Ratio Standard Deviation", ylab = "Yield (bu/ac)", cex = 1,cex.lab = 2, cex.axis = 2, pch=16)
#plot(plots_5cm$PlotCRR, plots_5cm$BuAc)
plot(plots_5cm$moran, plots_5cm$BuAc)
plot(plots_1cm$geary, plots_5cm$BuAc, xlab = "Geary's C", ylab = "Yield (bu/ac)", cex = 1,cex.lab = 2, cex.axis = 2, pch=16)
#tmp <- lm(BuAc~median, data=plots_5cm)
#summary(tmp)
#plot(plots_5cm$sill, plots_5cm$BuAc)
#plot(plots_5cm$range, plots_5cm$BuAc)
#plot(plots_5cm$Trt, plots_5cm$BuAc)
#plot(plots_5cm$AUDPC, plots_5cm$BuAc, xlab = "AUDPC", ylab = "Yield (bu/ac)")
#plot(plots_5cm$watersum, plots_5cm$BuAc)
#plot(plots_5cm$watersd, plots_5cm$BuAc)

par(mfrow=c(1,1))

# ggplot(data = plots_1cm, aes(x = carat, y = price)) + 
#   geom_point(aes(color = clarity), alpha = 0.3, position = "jitter")+
#   geom_smooth()+
#   facet_wrap(.~cut)+
#   theme_classic() + 
#   labs(
#     title = "Diamond (data) Mine",
#     y = "Price (US $)",
#     x = "Carat"
#   )

# Plot histogram distribution of each variable
histogram(plots_5cm$BuAc, nint=20)
histogram(plots_5cm$AUDPC, nint=20)
histogram(plots_5cm$watersum, nint=20)
histogram(plots_5cm$crrsd, nint=20)
histogram(plots_5cm$crrmean, nint=20)
histogram(plots_5cm$median, nint=20)
histogram(plots_5cm$sd, nint=20)
histogram(plots_5cm$skew, nint=20)
histogram(plots_5cm$kurt, nint=20)
histogram(plots_5cm$PlotCRR, nint=20)
histogram(plots_5cm$rumple, nint=20)
histogram(plots_5cm$geary, nint=20)

#----------------------------------------------- Statistical models ---------------------------------------------------------

# All : mean+median+sd+skew+PlotCRR+crrsd+crrmean+crrmedian+kurt+moran+Trt+AUDPC+watersum+watermean+watersd+rumple+geary+iqr

# Post colinearity 
# 1cm: median+PlotCRR+crrsd+crrmean+kurt+Trt+AUDPC+watersum+geary+iqr
# 5cm: median+PlotCRR+crrmean+Trt+AUDPC+watersum+geary+iqr
# 10cm: median+PlotCRR+crrsd+crrmean+moran+Trt+AUDPC+watersum

# normalize data
norm_1cm <- as.data.frame(scale(plots_1cm@data[,c(10, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35)]))
norm_1cm$Trt <- plots_1cm$Trt
norm_1cm$id <- plots_1cm$id

norm_5cm <- as.data.frame(scale(plots_5cm@data[,c(10, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35)]))
norm_5cm$Trt <- plots_5cm$Trt
norm_5cm$id <- plots_5cm$id

norm_10cm <- as.data.frame(scale(plots_10cm@data[,c(10, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35)]))
norm_10cm$Trt <- plots_10cm$Trt
norm_10cm$id <- plots_10cm$id

norm_20cm <- as.data.frame(scale(plots_20cm@data[,c(10, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35)]))
norm_20cm$Trt <- plots_20cm$Trt
norm_20cm$id <- plots_20cm$id

# 1cm
library(olsrr)
library(MASS)
library(leaps)
#median+PlotCRR+crrmean+watersum
ols_yield_1cm <- lm(BuAc ~  median+crrsd+crrmedian+skew+Trt+AUDPC+watersum+geary+iqr, data=norm_1cm)
ols_coll_diag(ols_yield_1cm)
step <- stepAIC(ols_yield_1cm, direction = "both")
step$anova

leaps<-regsubsets(BuAc ~  median+crrsd+crrmean+kurt+Trt+AUDPC+watersum+geary+iqr, data=norm_1cm,nbest=10)
# view results 
summary(leaps)
# plot a table of models showing variables in each model.
# models are ordered by the selection statistic.
plot(leaps,scale="r2")
# plot statistic by subset size 
library(car)
subsets(leaps, statistic="rsq")

norm_1cm$ols_yi_res <- ols_yield_1cm$residuals
summary(ols_yield_1cm)

w <- knn2nb(knearneigh(coordinates(plots), k=8))
moran.test(norm_1cm$ols_yi_res, nb2listw(w))

# backwards stepwise - median+AUDPC+geary
sp_err_yi_1cm <- errorsarlm(BuAc ~ mean+median+sd+skew+PlotCRR+crrsd+crrmean+crrmedian+kurt+moran+Trt+AUDPC+watersum+watermean+watersd+rumple+geary+iqr, data = norm_1cm, listw = nb2listw(w), zero.policy = T)
names(sp_err_yi_1cm)
summary(sp_err_yi_1cm)
norm_1cm$sp_err_resi <- sp_err_yi_1cm$residuals
moran.test(norm_1cm$sp_err_resi, nb2listw(w))

# Try RF regression to see if model can be improved

norm_1cm_df <- as.data.frame(norm_1cm)

set.seed(42)
train_data_indices <- rep(FALSE, nrow(norm_1cm_df))
train_data_indices[sample(1:nrow(norm_1cm_df), round(0.8 * nrow(norm_1cm_df)))] <- TRUE # randomly select 80% of the data for training
rf_regression_1cm<- randomForest(BuAc ~ median+skew+PlotCRR+crrsd+crrmean+moran+Trt+AUDPC+watersum+rumple, data=norm_1cm_df[train_data_indices, ], importance=T)
rf_regression_1cm
varImpPlot(rf_regression_1cm)
pred_yield <- predict(rf_regression_1cm, norm_1cm_df[!train_data_indices,]) # predict the rings
plot(norm_1cm_df$BuAc[!train_data_indices], pred_yield, xlab="Observed", ylab="Predicted")
abline(a=0, b=1, lty=2, col=2)

# 5cm

ols_yield_5cm <- lm(BuAc ~ median+skew+crrsd+crrmean+Trt+AUDPC+watersum+geary+rumple, data=norm_5cm)
ols_coll_diag(ols_yield_5cm)
norm_5cm$ols_yi_res <- ols_yield_5cm$residuals
summary(ols_yield_5cm)

w <- knn2nb(knearneigh(coordinates(plots), k=8))
moran.test(norm_5cm$ols_yi_res, nb2listw(w))

# backwards stepwise - median+Trt+geary
sp_err_yi_5cm <- errorsarlm(BuAc ~median+crrsd+crrmedian+skew+Trt+AUDPC+watersum+geary+iqr, data = norm_5cm, listw = nb2listw(w), zero.policy = T)
summary(sp_err_yi_5cm)
norm_5cm$sp_err_resi <- sp_err_yi_5cm$residuals
moran.test(norm_5cm$sp_err_resi, nb2listw(w))

# Try RF regression to see if model can be improved

norm_5cm_df <- as.data.frame(norm_5cm)

set.seed(42)
train_data_indices <- rep(FALSE, nrow(norm_5cm_df))
train_data_indices[sample(1:nrow(norm_5cm_df), round(0.8 * nrow(norm_5cm_df)))] <- TRUE # randomly select 80% of the data for training
rf_regression_5cm<- randomForest(BuAc ~ mean+sd+PlotCRR+crrmean+kurt+Trt+watersum, data=norm_5cm_df[train_data_indices, ], importance=T)
rf_regression_5cm
varImpPlot(rf_regression_5cm)
pred_yield <- predict(rf_regression_5cm, norm_5cm_df[!train_data_indices,]) # predict the rings
plot(norm_5cm_df$BuAc[!train_data_indices], pred_yield, xlab="Observed", ylab="Predicted")
abline(a=0, b=1, lty=2, col=2)

# 10cm

ols_yield_10cm <- lm(BuAc ~median+PlotCRR+crrsd+crrmean+moran+Trt+AUDPC+watersum, data=norm_10cm)
ols_coll_diag(ols_yield_10cm)
norm_10cm$ols_yi_res <- ols_yield_10cm$residuals
summary(ols_yield_10cm)

w <- knn2nb(knearneigh(coordinates(plots), k=8))
moran.test(norm_10cm$ols_yi_res, nb2listw(w))

# backwards stepwise - median+crrsd+Trt+iqr
sp_err_yi_10cm <- errorsarlm(BuAc ~ median+crrsd+crrmedian+skew+Trt+AUDPC+watersum+geary+iqr, data = norm_10cm, listw = nb2listw(w), zero.policy = T)
summary(sp_err_yi_10cm)
norm_10cm$sp_err_resi <- sp_err_yi_10cm$residuals
moran.test(norm_10cm$sp_err_resi, nb2listw(w))

# Try RF regression to see if model can be improved

norm_10cm_df <- as.data.frame(norm_10cm)

set.seed(42)
train_data_indices <- rep(FALSE, nrow(norm_10cm_df))
train_data_indices[sample(1:nrow(norm_10cm_df), round(0.8 * nrow(norm_10cm_df)))] <- TRUE # randomly select 80% of the data for training
rf_regression_10cm<- randomForest(BuAc ~ mean+sd+PlotCRR+crrmean+kurt+Trt+watersum, data=norm_10cm_df[train_data_indices, ], importance=T)
rf_regression_10cm
varImpPlot(rf_regression_10cm)
pred_yield <- predict(rf_regression_10cm, norm_10cm_df[!train_data_indices,]) # predict the rings
plot(norm_10cm_df$BuAc[!train_data_indices], pred_yield, xlab="Observed", ylab="Predicted")
abline(a=0, b=1, lty=2, col=2)

# 20cm

ols_yield_20cm <- lm(BuAc ~ mean+sd+skew+PlotCRR+crrsd+crrmean+kurt+moran+Trt+AUDPC+watersum, data=norm_20cm)
norm_20cm$ols_yi_res <- ols_yield_20cm$residuals
summary(ols_yield_20cm)

w <- knn2nb(knearneigh(coordinates(plots), k=8))
moran.test(norm_20cm$ols_yi_res, nb2listw(w))

# backwards stepwise - median+Trt
sp_err_yi_20cm <- errorsarlm(BuAc ~ median+crrsd+crrmedian+skew+Trt+AUDPC+watersum+geary+iqr, data = norm_20cm, listw = nb2listw(w), zero.policy = T)
summary(sp_err_yi_20cm)
norm_20cm$sp_err_resi <- sp_err_yi_20cm$residuals
moran.test(norm_20cm$sp_err_resi, nb2listw(w))

# Try RF regression to see if model can be improved

norm_20cm_df <- as.data.frame(norm_20cm)

set.seed(42)
train_data_indices <- rep(FALSE, nrow(norm_20cm_df))
train_data_indices[sample(1:nrow(norm_20cm_df), round(0.8 * nrow(norm_20cm_df)))] <- TRUE # randomly select 80% of the data for training
rf_regression_20cm<- randomForest(BuAc ~ mean+sd+PlotCRR+crrmean+kurt+Trt+watersum, data=norm_20cm_df[train_data_indices, ], importance=T)
rf_regression_20cm
varImpPlot(rf_regression_20cm)
pred_yield <- predict(rf_regression_20cm, norm_20cm_df[!train_data_indices,]) # predict the rings
plot(norm_20cm_df$BuAc[!train_data_indices], pred_yield, xlab="Observed", ylab="Predicted")
abline(a=0, b=1, lty=2, col=2)

# ------------------------------------- Table of regression results --------------------------------------------------------

extractResults <- function(x){
  sp_err_se <- x$rest.se
  sp_err_p <- summary(x)$Coef[,4]
  results <- cbind(sp_err_se, sp_err_p)
  results
}

coeff_5cm<-sp_err_yi_5cm$coefficients
coeff_1cm<-sp_err_yi_1cm$coefficients
coeff_10cm<-sp_err_yi_10cm$coefficients
coeff_20cm<-sp_err_yi_20cm$coefficients

sar_results_1cm <- extractResults(sp_err_yi_1cm)
sar_results_1cm <- cbind(coeff_1cm, sar_results_1cm)
sar_results_5cm <- extractResults(sp_err_yi_5cm)
sar_results_5cm <- cbind(coeff_5cm, sar_results_5cm)
sar_results_10cm <- extractResults(sp_err_yi_10cm)
sar_results_10cm <- cbind(coeff_10cm, sar_results_10cm)
sar_results_20cm <- extractResults(sp_err_yi_20cm)
sar_results_20cm <- cbind(coeff_20cm, sar_results_20cm)

write.csv(sar_results_1cm, "/media/Kellyn/F20E17B40E177139/kpmontgo@ncsu.edu/Research/Canopy_Morphology/results/sar_results_1cm.csv")
write.csv(sar_results_5cm, "/media/Kellyn/F20E17B40E177139/kpmontgo@ncsu.edu/Research/Canopy_Morphology/results/sar_results_5cm.csv")
write.csv(sar_results_10cm, "/media/Kellyn/F20E17B40E177139/kpmontgo@ncsu.edu/Research/Canopy_Morphology/results/sar_results_10cm.csv")
write.csv(sar_results_20cm, "/media/Kellyn/F20E17B40E177139/kpmontgo@ncsu.edu/Research/Canopy_Morphology/results/sar_results_20cm.csv")

# ------------------------------------------------------------ Figures -------------------------------------------------------------------

spplot(plots_10cm, zcol="BuAc")
spplot(plots_10cm, zcol="AUDPC")

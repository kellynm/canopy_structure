library(raster)
library(rgdal)
library(RColorBrewer)
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
library(classInt)


#Load raster data for all dates
setwd("/media/Kellyn/F20E17B40E177139/kpmontgo@ncsu.edu/LkWheeler_Sorghum/LkWheeler_Fusarium_Sorghum")
#setwd("Q:/My Drive/LkWheeler_Sorghum/LkWheeler_Fusarium_Sorghum")

sorghum_area <- readOGR("sorghum_area", "sorghum_area", stringsAsFactors = F)

csm_adj_9_1_1cm <- crop(raster('CSM/csm_adj_9_1_noDepth_1cm.tif'), sorghum_area)
csm_adj_9_1_5cm <- crop(raster('CSM/csm_adj_9_1_noDepth_5cm.tif'), sorghum_area)
csm_adj_9_1_10cm <- crop(raster('CSM/csm_adj_9_1_noDepth_10cm.tif'), sorghum_area)
csm_adj_9_1_20cm <- crop(raster('CSM/csm_adj_9_1_noDepth_20cm.tif'), sorghum_area)

csm_adj_9_1_1cm_velox <- velox(csm_adj_9_1_1cm)
csm_adj_9_1_5cm_velox <- velox(csm_adj_9_1_5cm)
csm_adj_9_1_10cm_velox <- velox(csm_adj_9_1_10cm)
csm_adj_9_1_20cm_velox <- velox(csm_adj_9_1_20cm)

#ortho_9_1 <- crop(raster("orthos/ortho_9_1_18.tif"))

#------------------------------------------------------- Field data-------------------------------------------------------------------
field_data <- fread("Sorghum 2018 JB Yield AUDPC.csv")
field_data$Plot <- as.integer(field_data$Plot)

# If normalizing:
field_data_tmp <- field_data[,3]
field_data_tmp$BuAc <- field_data$BuAc
field_data_tmp$AUDPC <- field_data$AUDPC
field_data_tmp <- as.data.frame(scale(field_data_tmp))
field_data_tmp$Trt <- field_data$Trt
field_data_tmp$Plot <- field_data$Plot
field_data <- field_data_tmp

plots <- readOGR("sorghum_plots", "sorghum_plots_small", stringsAsFactors = F)
plots <- spTransform(plots, CRS("+proj=lcc +lat_1=36.16666666666666 +lat_2=34.33333333333334 +lat_0=33.75 +lon_0=-79 +x_0=609601.22 +y_0=0 +datum=NAD83 +units=m +no_defs +ellps=GRS80 +towgs84=0,0,0"))
plots$id <- as.integer(plots$id)

plots_9_1_1cm <- plots
plots_9_1_5cm <- plots
plots_9_1_10cm <- plots
plots_9_1_20cm <- plots

plots_9_1_1cm <- merge(plots_9_1_1cm, field_data, by.x="id", by.y="Plot")
plots_9_1_5cm <- merge(plots_9_1_5cm, field_data, by.x="id", by.y="Plot")
plots_9_1_10cm <- merge(plots_9_1_10cm, field_data, by.x="id", by.y="Plot")
plots_9_1_20cm <- merge(plots_9_1_20cm, field_data, by.x="id", by.y="Plot")

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
norm_water <- as.data.frame(scale(waterDF[,2:4]))
norm_water$plots <- waterDF$plots

plots_9_1_1cm <- merge(plots_9_1_1cm, norm_water, by.x="id", by.y="plots")
plots_9_1_5cm <- merge(plots_9_1_5cm, norm_water, by.x="id", by.y="plots")
plots_9_1_10cm <- merge(plots_9_1_10cm, norm_water, by.x="id", by.y="plots")
plots_9_1_20cm <- merge(plots_9_1_20cm, norm_water, by.x="id", by.y="plots")

#---------------------------------------------------------- Canopy relief ratio----------------------------------------------------------------
crr <- function(x){
  (mean(x, na.rm=T)-min(x, na.rm=T))/(max(x, na.rm=T)-min(x, na.rm=T))
}

crrRast_9_1_1cm <- focal(csm_adj_9_1_1cm, w=matrix(1/2401, nc=49, nr=49), crr)
crrRast_9_1_5cm <- focal(csm_adj_9_1_5cm, w=matrix(1/81, nc=9, nr=9), crr)
crrRast_9_1_10cm <- focal(csm_adj_9_1_10cm, w=matrix(1/25, nc=5, nr=5), crr)
crrRast_9_1_20cm <- focal(csm_adj_9_1_20cm, w=matrix(1/9, nc=3, nr=3), crr)
crrRast_list <- list(crrRast_9_1_1cm,crrRast_9_1_5cm,crrRast_9_1_10cm,crrRast_9_1_20cm)

crr_metrics <- function(x){
  veloxRast <- velox(x)
  extract <- veloxRast$extract(sp=plots)
  names(extract) <- plots$id
  crrmean <- sapply(extract, mean, na.rm=T)
  crrmedian <- sapply(extract, median, na.rm=T)
  crrsd <- sapply(extract, sd, na.rm=T)
  crrDF <- data.frame(plots=plots$id, crrmedian= crrmedian, crrmean=crrmean, crrsd = crrsd)
  # Normalize data
  norm_crr <- as.data.frame(scale(crrDF[,2:4]))
  norm_crr$plots <- crrDF$plots
  norm_crr
}

crr_df_list <- lapply(crrRast_list, crr_metrics)

plots_9_1_1cm <- merge(plots_9_1_1cm, crr_df_list[[1]], by.x="id", by.y="plots")
plots_9_1_5cm <- merge(plots_9_1_5cm, crr_df_list[[2]], by.x="id", by.y="plots")
plots_9_1_10cm <- merge(plots_9_1_10cm, crr_df_list[[3]], by.x="id", by.y="plots")
plots_9_1_20cm <- merge(plots_9_1_20cm, crr_df_list[[4]], by.x="id", by.y="plots")


#---------------------------------------------------------- Crop height ----------------------------------------------------------
plot_extract_9_1_1cm <- csm_adj_9_1_1cm_velox$extract(sp=plots)
plot_extract_9_1_5cm <- csm_adj_9_1_5cm_velox$extract(sp=plots)
plot_extract_9_1_10cm <- csm_adj_9_1_10cm_velox$extract(sp=plots)
plot_extract_9_1_20cm <- csm_adj_9_1_20cm_velox$extract(sp=plots)
names(plot_extract_9_1_1cm) <- plots$id
names(plot_extract_9_1_5cm) <- plots$id
names(plot_extract_9_1_10cm) <- plots$id
names(plot_extract_9_1_20cm) <- plots$id

extract_list <- list(plot_extract_9_1_1cm, plot_extract_9_1_5cm, plot_extract_9_1_10cm, plot_extract_9_1_20cm)

ch_metrics <- function(x){
  CHmedian <- sapply(x, median, na.rm=TRUE)
  CHmean <- sapply(x, mean, na.rm=TRUE)
  CHsd <- sapply(x, sd, na.rm=TRUE)
  CHvar <- sapply(x, var, na.rm=TRUE)
  CHsum <- sapply(x, sum, na.rm=TRUE)
  CHskew <- sapply(x, skewness, na.rm=TRUE)
  CHkurt <- sapply(x, kurtosis, na.rm=TRUE)
  CHcrr <- sapply(x, crr)
  CHrange <- sapply(x, range, na.rm=TRUE)
  CHrange <- CHrange[2,]-CHrange[1,]
  cropHeightDF <- data.frame(plots=plots$id, median= CHmedian, mean=CHmean, sd = CHsd, var = CHvar, sum = CHsum, skew = CHskew, kurt= CHkurt, PlotCRR = CHcrr, chRange = CHrange)
  # Normalize data
  norm_ch <- as.data.frame(scale(cropHeightDF[,2:10]))
  norm_ch$plots <- cropHeightDF$plots
  norm_ch
}

ch_df_list <- lapply(extract_list, ch_metrics)

#merge all metrics into spatial polygons
plots_9_1_1cm <- merge(plots_9_1_1cm, ch_df_list[[1]], by.x="id", by.y="plots")
plots_9_1_5cm <- merge(plots_9_1_5cm, ch_df_list[[2]], by.x="id", by.y="plots")
plots_9_1_10cm <- merge(plots_9_1_10cm, ch_df_list[[3]], by.x="id", by.y="plots")
plots_9_1_20cm <- merge(plots_9_1_20cm, ch_df_list[[4]], by.x="id", by.y="plots")

# 10th and 90th percentile yield
pct_10_90_df <- function(x){
  top_plots <- c('410', '401','407', '406', '408')
  top_yield <- x[top_plots]
  top_yield_all <- do.call(c, top_yield)
  top_yield_all <- as.data.frame(top_yield_all)
  names(top_yield_all) <- "height"
  top_yield_all$group <- "top"
  
  bottom_plots <- c('203', '202','205', '204', '212')
  bottom_yield <- x[bottom_plots]
  bottom_yield_all <- do.call(c, bottom_yield)
  bottom_yield_all <- as.data.frame(bottom_yield_all)
  names(bottom_yield_all) <- "height"
  bottom_yield_all$group <- "bottom"
  
  pct_10_90_yield <- rbind(top_yield_all, bottom_yield_all)
}

pct_10_90_df_list <- lapply(extract_list, pct_10_90_df)

# Histogram and denisty of top vs bottom yield
library(ggplot2)

hist_10_90 <- function(x){
  height_mean <- ddply(x, "group", summarise, height.mean=mean(height))
  
  ggplot(x, aes(x=height, fill=group)) +
    geom_histogram(binwidth=.03, alpha=.5, position="identity") +
    geom_vline(data=height_mean, aes(xintercept=height.mean,  colour=group),
               linetype="dashed", size=1)
}

lapply(pct_10_90_df_list, hist_10_90)

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
csm_list <- list(csm_adj_9_1_1cm, csm_adj_9_1_5cm, csm_adj_9_1_10cm, csm_adj_9_1_20cm)
# var_df_list <- lapply(csm_list, variogram_metrics)
# 
# plots_9_1 <- merge(plots_9_1, norm_var, by.x="id", by.y="plots")


# ------------------------------------------------ spatial autocorrelation -----------------------------------------------------

autocor_metrics <- function(x){
  plotNums <- plots$id
  morans_all <- numeric(length(plotNums))
  geary_all <- numeric(length(plotNums))
  count <- 1
  
  for (num in plotNums){
    tmp <- plots@data
    position <- match(num, tmp$id)
    extract <- crop(x, plots[position,])
    moran <- Moran(extract)
    morans_all[count] <- moran
    geary <- Geary(extract)
    geary_all[count] <- geary
    count <- count+1
  }
  
  names(morans_all) <- plots$id
  names(geary_all) <- plots$id
  autocor_df <- data.frame(plots=plots$id, moran= morans_all, geary=geary_all)
  autocor_df
}

autocor_df_list <- lapply(csm_list, autocor_metrics)
plots_9_1_1cm <- merge(plots_9_1_1cm, autocor_df_list[[1]], by.x="id", by.y="plots")
plots_9_1_5cm <- merge(plots_9_1_5cm, autocor_df_list[[2]], by.x="id", by.y="plots")
plots_9_1_10cm <- merge(plots_9_1_10cm, autocor_df_list[[3]], by.x="id", by.y="plots")
plots_9_1_20cm <- merge(plots_9_1_20cm, autocor_df_list[[4]], by.x="id", by.y="plots")


# ------------------------------------Plot variables ------------------------------------------

imp_vars <- c("mean", "skew", "kurt", "sd", "crrsd", "crrmean", "moran", "Trt", "AUDPC", "watersum", "PlotCRR") 

plot(plots_9_1_1cm$mean, plots_9_1_1cm$BuAc)
plot(plots_9_1_1cm$sum, plots_9_1_1cm$BuAc)
plot(plots_9_1_1cm$skew, plots_9_1_1cm$BuAc)
plot(plots_9_1_1cm$kurt, plots_9_1_1cm$BuAc)
plot(plots_9_1_1cm$sd, plots_9_1_1cm$BuAc)
plot(plots_9_1_1cm$crrmean, plots_9_1_1cm$BuAc)
plot(plots_9_1_1cm$crrsd, plots_9_1_1cm$BuAc)
plot(plots_9_1_1cm$PlotCRR, plots_9_1_1cm$BuAc)
plot(plots_9_1_1cm$moran, plots_9_1_1cm$BuAc)
plot(plots_9_1_1cm$geary, plots_9_1_1cm$BuAc)
#plot(plots_9_1_1cm$sill, plots_9_1_1cm$BuAc)
#plot(plots_9_1_1cm$range, plots_9_1_1cm$BuAc)
plot(plots_9_1_1cm$Trt, plots_9_1_1cm$BuAc)
plot(plots_9_1_1cm$AUDPC, plots_9_1_1cm$BuAc)
plot(plots_9_1_1cm$watersum, plots_9_1_1cm$BuAc)
plot(plots_9_1_1cm$watersd, plots_9_1_1cm$BuAc)

# ggplot(data = plots_9_1_1cm, aes(x = carat, y = price)) + 
#   geom_point(aes(color = clarity), alpha = 0.3, position = "jitter")+
#   geom_smooth()+
#   facet_wrap(.~cut)+
#   theme_classic() + 
#   labs(
#     title = "Diamond (data) Mine",
#     y = "Price (US $)",
#     x = "Carat"
#   )


#----------------------------------------------- Statistical models ---------------------------------------------------------

#mean+sd+skew+PlotCRR+crrsd+crrmean+kurt+moran+Trt+AUDPC+watersum

# 1cm

ols_yield_9_1_1cm <- lm(BuAc ~ mean+sd+skew+PlotCRR+crrsd+crrmean+kurt+moran+Trt+AUDPC+watersum, data=plots_9_1_1cm)
plots_9_1_1cm$ols_yi_res <- ols_yield_9_1_1cm$residuals
summary(ols_yield_9_1_1cm)

w <- knn2nb(knearneigh(coordinates(plots), k=8))
moran.test(plots_9_1_1cm$ols_yi_res, nb2listw(w))

sp_err_yi_9_1_1cm <- errorsarlm(BuAc ~ mean+skew+crrmean+kurt+AUDPC, data = plots_9_1_1cm, listw = nb2listw(w), zero.policy = T)
summary(sp_err_yi_9_1_1cm)
plots_9_1_1cm$sp_err_resi <- sp_err_yi_9_1_1cm$residuals
moran.test(plots_9_1_1cm$sp_err_resi, nb2listw(w))

# Try RF regression to see if model can be improved

plots_9_1_1cm_df <- as.data.frame(plots_9_1_1cm)

set.seed(42)
train_data_indices <- rep(FALSE, nrow(plots_9_1_1cm_df))
train_data_indices[sample(1:nrow(plots_9_1_1cm_df), round(0.8 * nrow(plots_9_1_1cm_df)))] <- TRUE # randomly select 80% of the data for training
rf_regression_9_1_1cm<- randomForest(BuAc ~ mean+sd+PlotCRR+crrmean+kurt+Trt+watersum, data=plots_9_1_1cm_df[train_data_indices, ], importance=T)
rf_regression_9_1_1cm
varImpPlot(rf_regression_9_1_1cm)
pred_yield <- predict(rf_regression_9_1_1cm, plots_9_1_1cm_df[!train_data_indices,]) # predict the rings
plot(plots_9_1_1cm_df$BuAc[!train_data_indices], pred_yield, xlab="Observed", ylab="Predicted")
abline(a=0, b=1, lty=2, col=2)

# 5cm

ols_yield_9_1_5cm <- lm(BuAc ~ mean+sd+skew+PlotCRR+crrsd+crrmean+kurt+moran+Trt+AUDPC+watersum, data=plots_9_1_5cm)
plots_9_1_5cm$ols_yi_res <- ols_yield_9_1_5cm$residuals
summary(ols_yield_9_1_5cm)

w <- knn2nb(knearneigh(coordinates(plots), k=8))
moran.test(plots_9_1_5cm$ols_yi_res, nb2listw(w))

sp_err_yi_9_1_5cm <- errorsarlm(BuAc ~ mean+sd+PlotCRR+crrmean+kurt+AUDPC, data = plots_9_1_5cm, listw = nb2listw(w), zero.policy = T)
summary(sp_err_yi_9_1_5cm)
plots_9_1_5cm$sp_err_resi <- sp_err_yi_9_1_5cm$residuals
moran.test(plots_9_1_5cm$sp_err_resi, nb2listw(w))

# Try RF regression to see if model can be improved

plots_9_1_5cm_df <- as.data.frame(plots_9_1_5cm)

set.seed(42)
train_data_indices <- rep(FALSE, nrow(plots_9_1_5cm_df))
train_data_indices[sample(1:nrow(plots_9_1_5cm_df), round(0.8 * nrow(plots_9_1_5cm_df)))] <- TRUE # randomly select 80% of the data for training
rf_regression_9_1_5cm<- randomForest(BuAc ~ mean+sd+PlotCRR+crrmean+kurt+Trt+watersum, data=plots_9_1_5cm_df[train_data_indices, ], importance=T)
rf_regression_9_1_5cm
varImpPlot(rf_regression_9_1_5cm)
pred_yield <- predict(rf_regression_9_1_5cm, plots_9_1_5cm_df[!train_data_indices,]) # predict the rings
plot(plots_9_1_5cm_df$BuAc[!train_data_indices], pred_yield, xlab="Observed", ylab="Predicted")
abline(a=0, b=1, lty=2, col=2)

# 10cm

ols_yield_9_1_10cm <- lm(BuAc ~ mean+sd+skew+PlotCRR+crrsd+crrmean+kurt+moran+Trt+AUDPC+watersum, data=plots_9_1_10cm)
plots_9_1_10cm$ols_yi_res <- ols_yield_9_1_10cm$residuals
summary(ols_yield_9_1_10cm)

w <- knn2nb(knearneigh(coordinates(plots), k=8))
moran.test(plots_9_1_10cm$ols_yi_res, nb2listw(w))

sp_err_yi_9_1_10cm <- errorsarlm(BuAc ~ mean+sd+crrsd+crrmean+AUDPC, data = plots_9_1_10cm, listw = nb2listw(w), zero.policy = T)
summary(sp_err_yi_9_1_10cm)
plots_9_1_10cm$sp_err_resi <- sp_err_yi_9_1_10cm$residuals
moran.test(plots_9_1_10cm$sp_err_resi, nb2listw(w))

# Try RF regression to see if model can be improved

plots_9_1_10cm_df <- as.data.frame(plots_9_1_10cm)

set.seed(42)
train_data_indices <- rep(FALSE, nrow(plots_9_1_10cm_df))
train_data_indices[sample(1:nrow(plots_9_1_10cm_df), round(0.8 * nrow(plots_9_1_10cm_df)))] <- TRUE # randomly select 80% of the data for training
rf_regression_9_1_10cm<- randomForest(BuAc ~ mean+sd+PlotCRR+crrmean+kurt+Trt+watersum, data=plots_9_1_10cm_df[train_data_indices, ], importance=T)
rf_regression_9_1_10cm
varImpPlot(rf_regression_9_1_10cm)
pred_yield <- predict(rf_regression_9_1_10cm, plots_9_1_10cm_df[!train_data_indices,]) # predict the rings
plot(plots_9_1_10cm_df$BuAc[!train_data_indices], pred_yield, xlab="Observed", ylab="Predicted")
abline(a=0, b=1, lty=2, col=2)

# 20cm

ols_yield_9_1_20cm <- lm(BuAc ~ mean+sd+skew+PlotCRR+crrsd+crrmean+kurt+moran+Trt+AUDPC+watersum, data=plots_9_1_20cm)
plots_9_1_20cm$ols_yi_res <- ols_yield_9_1_20cm$residuals
summary(ols_yield_9_1_20cm)

w <- knn2nb(knearneigh(coordinates(plots), k=8))
moran.test(plots_9_1_20cm$ols_yi_res, nb2listw(w))

sp_err_yi_9_1_20cm <- errorsarlm(BuAc ~ mean+sd+PlotCRR+kurt+moran+AUDPC, data = plots_9_1_20cm, listw = nb2listw(w), zero.policy = T)
summary(sp_err_yi_9_1_20cm)
plots_9_1_20cm$sp_err_resi <- sp_err_yi_9_1_20cm$residuals
moran.test(plots_9_1_20cm$sp_err_resi, nb2listw(w))

# Try RF regression to see if model can be improved

plots_9_1_20cm_df <- as.data.frame(plots_9_1_20cm)

set.seed(42)
train_data_indices <- rep(FALSE, nrow(plots_9_1_20cm_df))
train_data_indices[sample(1:nrow(plots_9_1_20cm_df), round(0.8 * nrow(plots_9_1_20cm_df)))] <- TRUE # randomly select 80% of the data for training
rf_regression_9_1_20cm<- randomForest(BuAc ~ mean+sd+PlotCRR+crrmean+kurt+Trt+watersum, data=plots_9_1_20cm_df[train_data_indices, ], importance=T)
rf_regression_9_1_20cm
varImpPlot(rf_regression_9_1_20cm)
pred_yield <- predict(rf_regression_9_1_20cm, plots_9_1_20cm_df[!train_data_indices,]) # predict the rings
plot(plots_9_1_20cm_df$BuAc[!train_data_indices], pred_yield, xlab="Observed", ylab="Predicted")
abline(a=0, b=1, lty=2, col=2)

# ------------------------------------- Table of regression results --------------------------------------------------------

extractResults <- function(x){
  sp_err_se <- x$rest.se
  sp_err_p <- summary(x)$Coef[,4]
  results <- cbind(sp_err_se, sp_err_p)
  results
}

coeff_5cm<-sp_err_yi_9_1_5cm$coefficients
coeff_1cm<-sp_err_yi_9_1_1cm$coefficients
coeff_10cm<-sp_err_yi_9_1_10cm$coefficients
coeff_20cm<-sp_err_yi_9_1_20cm$coefficients

sar_results_1cm <- extractResults(sp_err_yi_9_1_1cm)
sar_results_1cm <- cbind(coeff_1cm, sar_results_1cm)
sar_results_5cm <- extractResults(sp_err_yi_9_1_5cm)
sar_results_5cm <- cbind(coeff_5cm, sar_results_5cm)
sar_results_10cm <- extractResults(sp_err_yi_9_1_10cm)
sar_results_10cm <- cbind(coeff_10cm, sar_results_10cm)
sar_results_20cm <- extractResults(sp_err_yi_9_1_20cm)
sar_results_20cm <- cbind(coeff_20cm, sar_results_20cm)

write.csv(sar_results_1cm, "/media/Kellyn/F20E17B40E177139/kpmontgo@ncsu.edu/Research/Canopy_Morphology/results/sar_results_1cm.csv")
write.csv(sar_results_5cm, "/media/Kellyn/F20E17B40E177139/kpmontgo@ncsu.edu/Research/Canopy_Morphology/results/sar_results_5cm.csv")
write.csv(sar_results_10cm, "/media/Kellyn/F20E17B40E177139/kpmontgo@ncsu.edu/Research/Canopy_Morphology/results/sar_results_10cm.csv")
write.csv(sar_results_20cm, "/media/Kellyn/F20E17B40E177139/kpmontgo@ncsu.edu/Research/Canopy_Morphology/results/sar_results_20cm.csv")

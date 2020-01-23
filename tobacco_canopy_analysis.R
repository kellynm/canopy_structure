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
library(caret)

#Load raster data for all dates
#setwd("/media/Kellyn/F20E17B40E177139/kpmontgo@ncsu.edu/LkWheeler_Sorghum/LkWheeler_Fusarium_Sorghum")
setwd("Q:/My Drive/Research/Canopy_Morphology/Tobacco_Project/2018/canopy_analysis")

tobacco_area <- readOGR("Q:/My Drive/Research/Canopy_Morphology/Tobacco_Project/2018/canopy_analysis/layers/Boundary", "Boundary", stringsAsFactors = F)
tobacco_area <- spTransform(tobacco_area, CRS("+proj=lcc +lat_1=36.16666666666666 +lat_2=34.33333333333334 +lat_0=33.75 +lon_0=-79 +x_0=609601.22 +y_0=0 +datum=NAD83 +units=m +no_defs +ellps=GRS80 +towgs84=0,0,0"))

csm_618 <- crop(raster('Q:/My Drive/Research/Canopy_Morphology/Tobacco_Project/2018/canopy_analysis/layers/CSM/CSM_618.tif'), tobacco_area) #convert to cm
csm_716 <- crop(raster('Q:/My Drive/Research/Canopy_Morphology/Tobacco_Project/2018/canopy_analysis/layers/CSM/CSM_716.tif'), tobacco_area)
csm_814 <- crop(raster('Q:/My Drive/Research/Canopy_Morphology/Tobacco_Project/2018/canopy_analysis/layers/CSM/CSM_814.tif'), tobacco_area)

csm_618_velox <- velox(csm_618)
csm_716_velox <- velox(csm_716)
csm_814_velox <- velox(csm_814)


plots <- readOGR("Q:/My Drive/Research/Canopy_Morphology/Tobacco_Project/2018/canopy_analysis/layers/plots", "plots_reshape", stringsAsFactors = F)
plots$Treatment <- as.factor(plots$Treatment)
plots <- spTransform(plots, CRS("+proj=lcc +lat_1=36.16666666666666 +lat_2=34.33333333333334 +lat_0=33.75 +lon_0=-79 +x_0=609601.22 +y_0=0 +datum=NAD83 +units=m +no_defs +ellps=GRS80 +towgs84=0,0,0"))

fert_plot_ids <- plots$Id[grep("-C", c(plots$Id, "-C"))]
fert_plots <- plots$Id %in% fert_plot_ids

plots$nutr_id <- substr(plots$Id, 1, 3)
plots$nutr_id[fert_plots] <- fert_plot_ids

nutrient <- fread("Q:/My Drive/Research/Canopy_Morphology/Tobacco_Project/2018/2018 Tobacco Field Tissue Data.csv")
names(nutrient) <- c("Plot","Treatment", "Sample #",  "N","P","K","Mg", "Ca","S","Zn", "Mn", "Cu", "Fe", "B", "Al", "Na") 

plots <- merge(plots, nutrient, by.x="nutr_id", by.y='Plot')

plots_618 <- plots
plots_716 <- plots
plots_814 <- plots

# ADD WATER MODEL?

#---------------------------------------------------------- Canopy relief ratio----------------------------------------------------------------
crr <- function(x){
  (mean(x, na.rm=T)-min(x, na.rm=T))/(max(x, na.rm=T)-min(x, na.rm=T))
}

crrRast_618 <- focal(csm_618, w=matrix(1/81, nc=9, nr=9), crr)
crrRast_716 <- focal(csm_716, w=matrix(1/81, nc=9, nr=9), crr)
crrRast_814 <- focal(csm_814, w=matrix(1/81, nc=9, nr=9), crr)
crrRast_list <- list(crrRast_618,crrRast_716,crrRast_814)

# Check if crr within plot is normally distributed
# crr_618_df <- as.data.frame(crrRast_618)
# shapiro.test(sample(crr_618_df$layer, 5000))
# histogram(crr_618_df$layer)
# velox_crr_618 <- velox(crrRast_618)
# crr_extract_618 <- velox_crr_618$extract(sp=plots)
# names(crr_extract_618) <- plots$Id

#histogram(crr_extract_618$`502-E`)
#shapiro.test(sample(crr_extract_618$`502-E`, 5000))
#ggqqplot(sample(crr_extract_618$`502-E`, 5000))

crr_metrics <- function(x){
  veloxRast <- velox(x)
  extract <- veloxRast$extract(sp=plots)
  names(extract) <- plots$Id
  crrmean <- sapply(extract, mean, na.rm=T)
  crrmedian <- sapply(extract, median, na.rm=T)
  crriqr <- sapply(extract, IQR, na.rm=T)
  crrDF <- data.frame(plots=plots$Id, crrmedian= crrmedian, crrmean=crrmean, crriqr = crriqr)
  crrDF
  # Normalize data
  # norm_crr <- as.data.frame(scale(crrDF[,2:4]))
  # norm_crr$plots <- crrDF$plots
  # norm_crr
}

crr_df_list <- lapply(crrRast_list, crr_metrics)

plots_618 <- merge(plots_618, crr_df_list[[1]], by.x="Id", by.y="plots")
plots_716 <- merge(plots_716, crr_df_list[[2]], by.x="Id", by.y="plots")
plots_814 <- merge(plots_814, crr_df_list[[3]], by.x="Id", by.y="plots")


#---------------------------------------------------------- Crop height ----------------------------------------------------------
plot_extract_618 <- csm_618_velox$extract(sp=plots)
plot_extract_716 <- csm_716_velox$extract(sp=plots)
plot_extract_814 <- csm_814_velox$extract(sp=plots)
names(plot_extract_618) <- plots$Id
names(plot_extract_716) <- plots$Id
names(plot_extract_814) <- plots$Id

extract_list <- list(plot_extract_618, plot_extract_716, plot_extract_814)

# Check for multimodality
# dip.test(plot_extract_618$`502-E`)
# histogram(plot_extract_618$`502-E`)
# 
# multimodality_test <- function(x){
#   plotNums <- plots$Id
#   dip_all <- numeric(length(plotNums))
#   names(dip_all) <- "dip"
#   count <- 1
#   for (plot in x){
#     test <- dip.test(plot)
#     p <- test$p.value
#     dip_all[count] <- p
#     count <- count+1
#   }
#   dip_df <- data.frame(plots=plots$Id, dip= dip_all)
#   dip_df
# }
# 
# dip_df_list <- lapply(extract_list, multimodality_test)
# 
# #merge all metrics into spatial polygons
# plots_618 <- merge(plots_618, dip_df_list[[1]], by.x="Id", by.y="plots")
# plots_716 <- merge(plots_716, dip_df_list[[2]], by.x="Id", by.y="plots")
# plots_814 <- merge(plots_814, dip_df_list[[3]], by.x="Id", by.y="plots")


# Check for normality
#shapiro.test(plot_extract_618$`502-E`)

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
  cropHeightDF <- data.frame(plots=plots$Id, median= CHmedian, mean=CHmean, sd = CHsd, iqr = CHiqr, var = CHvar, sum = CHsum, skew = CHskew, kurt= CHkurt, PlotCRR = CHcrr, chRange = CHrange)
  cropHeightDF
}

ch_df_list <- lapply(extract_list, ch_metrics)

#merge all metrics into spatial polygons
plots_618 <- merge(plots_618, ch_df_list[[1]], by.x="Id", by.y="plots")
plots_716 <- merge(plots_716, ch_df_list[[2]], by.x="Id", by.y="plots")
plots_814 <- merge(plots_814, ch_df_list[[3]], by.x="Id", by.y="plots")

#--------------------------------------------------------- rumple index -----------------------------------------------------------------
library(lidR)

rumple <- function(x){
  plotNums <- plots$Id
  rumple_all <- numeric(length(plotNums))
  names(rumple_all) <- "rumple"
  count <- 1
  for (num in plotNums){
    tmp <- plots@data
    position <- match(num, tmp$Id)
    extract <- mask(x, plots[position,])
    extract_df <- rasterToPoints(extract)
    rumple <- rumple_index(extract_df[,1], extract_df[,2], extract_df[,3])
    rumple_all[count] <- rumple
    count <- count+1
  }
  rumple_df <- data.frame(plots=plots$Id, rumple= rumple_all)
  rumple_df
}



min(plot_extract_618$`501-E`)

csm_list <- list(csm_618, csm_716, csm_814)
rumple_df_list <- lapply(csm_list, rumple)

#merge all metrics into spatial polygons
plots_618 <- merge(plots_618, rumple_df_list[[1]], by.x="Id", by.y="plots")
plots_716 <- merge(plots_716, rumple_df_list[[2]], by.x="Id", by.y="plots")
plots_814 <- merge(plots_814, rumple_df_list[[3]], by.x="Id", by.y="plots")


# ------------------------------------------------ spatial autocorrelation -----------------------------------------------------

autocor_metrics <- function(x, w){
  plotNums <- plots$Id
  morans_all <- numeric(length(plotNums))
  #moran_local_all <- numeric(length(plotNums))
  geary_all <- numeric(length(plotNums))
  count <- 1
  
  for (num in plotNums){
    tmp <- plots@data
    position <- match(num, tmp$Id)
    extract <- crop(x, plots[position,])
    moran <- Moran(extract, w)
    morans_all[count] <- moran
    #moran_local <- MoranLocal(extract, w)
    #moran_local_all[count] <- moran_local
    geary <- Geary(extract, w)
    geary_all[count] <- geary
    count <- count+1
  }
  
  names(morans_all) <- plots$Id
  #names(morans_local_all) <- plots$Id
  names(geary_all) <- plots$Id
  autocor_df <- data.frame(plots=plots$Id, moran= morans_all, geary=geary_all)
  autocor_df
}

#Using approx 30 cm x 30 cm neighborhood matrix
autocor_618 <- autocor_metrics(csm_618, matrix(c(1,1,1,1,1,1,1,1,1,1,1,1,0,1,1,1,1,1,1,1,1,1,1,1,1), nc=5, nr=5))
autocor_716 <- autocor_metrics(csm_716, matrix(c(1,1,1,1,1,1,1,1,1,1,1,1,0,1,1,1,1,1,1,1,1,1,1,1,1), nc=5, nr=5))
autocor_814 <- autocor_metrics(csm_814, matrix(c(1,1,1,1,1,1,1,1,1,1,1,1,0,1,1,1,1,1,1,1,1,1,1,1,1), nc=5, nr=5))

plots_618 <- merge(plots_618, autocor_618, by.x="Id", by.y="plots")
plots_716 <- merge(plots_716, autocor_716, by.x="Id", by.y="plots")
plots_814 <- merge(plots_814, autocor_814, by.x="Id", by.y="plots")


# Create rasters representing local Moran's I for visualization
# autocor_rasters <- function(x, w){
#   plotNums <- plots$Id
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
#   names(moran_local_all) <- plots$Id
#   moran_local_all
# }
# 
# local_moran_plots_618 <- autocor_rasters(csm_618, matrix(c(1,1,1,1,1,1,1,1,1,1,1,1,0,1,1,1,1,1,1,1,1,1,1,1,1), nc=5, nr=5))
# names(local_moran_plots_618) <- NULL
# local_moran_plots_618$filename <- 'moran_618.tif'
# local_moran_plots_618$overwrite <- TRUE
# local_moran_plots_618$fun <- sum
# merged_local_moran_plots_618 <- do.call(raster::merge, local_moran_plots_618)
# #writeRaster(merged_local_moran_plots_618, "Q:/My Drive/Research/Canopy_Morphology/Tobacco_Project/canopy_analysis/layers/moran/local_moran_plots_618.tif", format="GTiff", overwrite = T)
# 
# local_moran_plots_716 <- autocor_rasters(csm_716, matrix(c(1,1,1,1,1,1,1,1,1,1,1,1,0,1,1,1,1,1,1,1,1,1,1,1,1), nc=5, nr=5))
# names(local_moran_plots_716) <- NULL
# local_moran_plots_716$filename <- 'moran_716.tif'
# local_moran_plots_716$overwrite <- TRUE
# local_moran_plots_716$fun <- sum
# merged_local_moran_plots_716 <- do.call(raster::merge, local_moran_plots_716)
# #writeRaster(merged_local_moran_plots_716, "Q:/My Drive/Research/Canopy_Morphology/Tobacco_Project/canopy_analysis/layers/moran/local_moran_plots_716.tif", format="GTiff", overwrite = T)
# 
# local_moran_plots_814 <- autocor_rasters(csm_814, matrix(c(1,1,1,1,1,1,1,1,1,1,1,1,0,1,1,1,1,1,1,1,1,1,1,1,1), nc=5, nr=5))
# names(local_moran_plots_814) <- NULL
# local_moran_plots_814$filename <- 'moran_814.tif'
# local_moran_plots_814$overwrite <- TRUE
# local_moran_plots_814$fun <- sum
# merged_local_moran_plots_814 <- do.call(raster::merge, local_moran_plots_814)
# #writeRaster(merged_local_moran_plots_814, "Q:/My Drive/Research/Canopy_Morphology/Tobacco_Project/canopy_analysis/layers/moran/local_moran_plots_814.tif", format="GTiff", overwrite = T)

# ------------------------------------- Aggregate plots to match nutrient observations --------------------------------------

agg_618 <- aggregate(plots_618@data[,6:34], by = list(plots_618$nutr_id), FUN = mean)

agg_716 <- aggregate(plots_716@data[,6:34], by = list(plots_716$nutr_id), FUN = mean)

agg_814 <- aggregate(plots_814@data[,6:34], by = list(plots_814$nutr_id), FUN = mean)

 
# --------------------------------------- Stats data prep --------------------------------------------------------------------------------

normalize <- function(x)
{
  return((x- min(x)) /(max(x)-min(x)))
}

#minmax normalization of data
norm_618 <- as.data.frame(normalize(agg_618[,2:30]))
norm_618$Id <- agg_618$Group.1

norm_716 <- as.data.frame(normalize(agg_716[,2:30]))
norm_716$Id <- agg_716$Group.1

norm_814 <- as.data.frame(normalize(agg_814[,2:30]))
norm_814$Id <- agg_814$Group.1


#---------------------------------------- LM without spectral  ------------------------------------------

lm_N_618 <- lm(N ~ median+iqr+crrmedian+crriqr+skew+kurt+moran+rumple, data = norm_618, na.action=na.fail)
summary(lm_N_618)
vif(lm_N_618)
modsel_N_618 <- dredge(lm_N_618)
topmod_N_618 <- lm(N ~ crriqr+iqr+rumple, data = norm_618,na.action=na.fail)
summary(topmod_N_618)
N_618_pred <- predict(topmod_N_618)
N_618_res <- resid(topmod_N_618)
n_618_table <- cbind(norm_618$N, N_618_pred, N_618_res)
head(n_618_table)

lm_P_618 <- lm(P ~ median+iqr+crrmedian+crriqr+skew+kurt+moran+rumple, data = norm_618, na.action=na.fail)
summary(lm_P_618)
modsel_P_618 <- dredge(lm_P_618)
topmod_P_618 <- lm(P ~ rumple+skew, data = norm_618,na.action=na.fail)
summary(topmod_P_618)

lm_K_618 <- lm(K ~ median+iqr+crrmedian+crriqr+skew+kurt+moran+rumple, data = norm_618, na.action=na.fail)
summary(lm_K_618)
modsel_K_618 <- dredge(lm_K_618)
topmod_K_618 <- lm(K ~ iqr+median+rumple+skew, data = norm_618,na.action=na.fail)
summary(topmod_K_618)

lm_Ca_618 <- lm(Ca ~ median+iqr+crrmedian+crriqr+skew+kurt+moran+rumple, data = norm_618, na.action=na.fail)
summary(lm_Ca_618)
modsel_Ca_618 <- dredge(lm_Ca_618)
topmod_Ca_618 <- lm(Ca ~ median+skew, data = norm_618,na.action=na.fail)
summary(topmod_Ca_618)

lm_S_618 <- lm(S ~ median+iqr+crrmedian+crriqr+skew+kurt+moran+rumple, data = norm_618, na.action=na.fail)
summary(lm_S_618)
modsel_S_618 <- dredge(lm_S_618)
topmod_S_618 <- lm(S ~ crriqr+median+iqr+skew, data = norm_618,na.action=na.fail)
summary(topmod_S_618)

lm_B_618 <- lm(B ~ median+iqr+crrmedian+crriqr+skew+kurt+moran+rumple, data = norm_618, na.action=na.fail)
summary(lm_B_618)
modsel_B_618 <- dredge(lm_B_618)
topmod_B_618 <- lm(B ~ median+moran+skew, data = norm_618,na.action=na.fail)
summary(topmod_B_618)

lm_Zn_618 <- lm(Zn ~ median+iqr+crrmedian+crriqr+skew+kurt+moran+rumple, data = norm_618, na.action=na.fail)
summary(lm_Zn_618)
modsel_Zn_618 <- dredge(lm_Zn_618)
topmod_Zn_618 <- lm(Zn ~ iqr+median+skew, data = norm_618,na.action=na.fail)
summary(topmod_Zn_618)

lm_Mg_618 <- lm(Mg ~ median+iqr+crrmedian+crriqr+skew+kurt+moran+rumple, data = norm_618, na.action=na.fail)
summary(lm_Mg_618)
modsel_Mg_618 <- dredge(lm_Mg_618)
topmod_Mg_618 <- lm(Mg ~ crriqr+iqr+median+rumple, data = norm_618,na.action=na.fail)
summary(topmod_Mg_618)

# function to get lm pvalue
lmp <- function (modelobject) {
  if (class(modelobject) != "lm") stop("Not an object of class 'lm' ")
  f <- summary(modelobject)$fstatistic
  p <- pf(f[1],f[2],f[3],lower.tail=F)
  attributes(p) <- NULL
  return(p)
}

res_N_618 <- as.data.frame(modsel_N_618[1])
res_N_618$adj_r2 <- summary(topmod_N_618)$adj.r.squared
res_N_618$pvalue <- lmp(topmod_N_618)

res_P_618 <- as.data.frame(modsel_P_618[1])
res_P_618$adj_r2 <- summary(topmod_P_618)$adj.r.squared
res_P_618$pvalue <- lmp(topmod_P_618)

res_K_618 <- as.data.frame(modsel_K_618[1])
res_K_618$adj_r2 <- summary(topmod_K_618)$adj.r.squared
res_K_618$pvalue <- lmp(topmod_K_618)

res_B_618 <- as.data.frame(modsel_B_618[1])
res_B_618$adj_r2 <- summary(topmod_B_618)$adj.r.squared
res_B_618$pvalue <- lmp(topmod_B_618)

res_Ca_618 <- as.data.frame(modsel_Ca_618[1])
res_Ca_618$adj_r2 <- summary(topmod_Ca_618)$adj.r.squared
res_Ca_618$pvalue <- lmp(topmod_Ca_618)

res_Mg_618 <- as.data.frame(modsel_Mg_618[1])
res_Mg_618$adj_r2 <- summary(topmod_Mg_618)$adj.r.squared
res_Mg_618$pvalue <- lmp(topmod_Mg_618)

res_S_618 <- as.data.frame(modsel_S_618[1])
res_S_618$adj_r2 <- summary(topmod_S_618)$adj.r.squared
res_S_618$pvalue <- lmp(topmod_S_618)

res_Zn_618 <- as.data.frame(modsel_Zn_618[1])
res_Zn_618$adj_r2 <- summary(topmod_Zn_618)$adj.r.squared
res_Zn_618$pvalue <- lmp(topmod_Zn_618)

res_618 <- rbind(res_N_618, res_P_618, res_K_618, res_B_618, res_Ca_618, res_Mg_618, res_S_618, res_Zn_618)
res_618$dependent <- c("N", "P", "K", "B", "Ca", "Mg", "S", "Zn")
write.csv(res_618, "results/res_618.csv")

#------------------------------------------- LM with spectral --------------------------------------------------------
lm_N_618 <- lm(N ~ median+iqr+crrmedian+crriqr+skew+kurt+moran+rumple+spad+varInd_median+tgi_median, data = agg_618, na.action=na.fail)
summary(lm_N_618)
vif(lm_N_618)
modsel_N_618 <- dredge(lm_N_618)
topmod_N_618 <- lm(N ~ crrmedian+iqr+kurt+spad+varInd_median, data = agg_618,na.action=na.fail)
summary(topmod_N_618)
N_618_pred <- predict(topmod_N_618)
N_618_res <- resid(topmod_N_618)
n_618_table <- cbind(agg_618$N, N_618_pred, N_618_res)
head(n_618_table)

lm_P_618 <- lm(P ~ median+iqr+crrmedian+crriqr+skew+kurt+moran+rumple+spad+varInd_median+tgi_median, data = agg_618, na.action=na.fail)
summary(lm_P_618)
modsel_P_618 <- dredge(lm_P_618)
topmod_P_618 <- lm(P ~ crriqr+iqr+median+moran+skew+spad, data = agg_618,na.action=na.fail)
summary(topmod_P_618)

lm_K_618 <- lm(K ~ median+iqr+crrmedian+crriqr+skew+kurt+moran+rumple+spad+varInd_median+tgi_median, data = agg_618, na.action=na.fail)
summary(lm_K_618)
modsel_K_618 <- dredge(lm_K_618)
topmod_K_618 <- lm(K ~ varInd_median+tgi_median, data = agg_618,na.action=na.fail)
summary(topmod_K_618)

lm_Ca_618 <- lm(Ca ~ median+iqr+crrmedian+crriqr+skew+kurt+moran+rumple+spad+varInd_median+tgi_median, data = agg_618, na.action=na.fail)
summary(lm_Ca_618)
modsel_Ca_618 <- dredge(lm_Ca_618)
topmod_Ca_618 <- lm(Ca ~ iqr+moran, data = agg_618,na.action=na.fail)
summary(topmod_Ca_618)

lm_S_618 <- lm(S ~ median+iqr+crrmedian+crriqr+skew+kurt+moran+rumple+spad+varInd_median+tgi_median, data = agg_618, na.action=na.fail)
summary(lm_S_618)
modsel_S_618 <- dredge(lm_S_618)
topmod_S_618 <- lm(S ~ median, data = agg_618,na.action=na.fail)
summary(topmod_S_618)

lm_B_618 <- lm(B ~ median+iqr+crrmedian+crriqr+skew+kurt+moran+rumple+spad+varInd_median+tgi_median, data = agg_618, na.action=na.fail)
summary(lm_B_618)
modsel_B_618 <- dredge(lm_B_618)
topmod_B_618 <- lm(B ~ varInd_median+moran, data = agg_618,na.action=na.fail)
summary(topmod_B_618)

lm_Zn_618 <- lm(Zn ~ median+iqr+crrmedian+crriqr+skew+kurt+moran+rumple+spad+varInd_median+tgi_median, data = agg_618, na.action=na.fail)
summary(lm_Zn_618)
modsel_Zn_618 <- dredge(lm_Zn_618)
topmod_Zn_618 <- lm(Zn ~ crrmedian+spad, data = agg_618,na.action=na.fail)
summary(topmod_Zn_618)

lm_Mg_618 <- lm(Mg ~ median+iqr+crrmedian+crriqr+skew+kurt+moran+rumple+spad+varInd_median+tgi_median, data = agg_618, na.action=na.fail)
summary(lm_Mg_618)
modsel_Mg_618 <- dredge(lm_Mg_618)
topmod_Mg_618 <- lm(Mg ~ crrmedian+kurt+median+spad+tgi_median, data = agg_618,na.action=na.fail)
summary(topmod_Mg_618)

res_N_618 <- as.data.frame(modsel_N_618[1])
res_N_618$adj_r2 <- summary(topmod_N_618)$adj.r.squared
res_N_618$pvalue <- lmp(topmod_N_618)

res_P_618 <- as.data.frame(modsel_P_618[1])
res_P_618$adj_r2 <- summary(topmod_P_618)$adj.r.squared
res_P_618$pvalue <- lmp(topmod_P_618)

res_K_618 <- as.data.frame(modsel_K_618[1])
res_K_618$adj_r2 <- summary(topmod_K_618)$adj.r.squared
res_K_618$pvalue <- lmp(topmod_K_618)

res_B_618 <- as.data.frame(modsel_B_618[1])
res_B_618$adj_r2 <- summary(topmod_B_618)$adj.r.squared
res_B_618$pvalue <- lmp(topmod_B_618)

res_Ca_618 <- as.data.frame(modsel_Ca_618[1])
res_Ca_618$adj_r2 <- summary(topmod_Ca_618)$adj.r.squared
res_Ca_618$pvalue <- lmp(topmod_Ca_618)

res_Mg_618 <- as.data.frame(modsel_Mg_618[1])
res_Mg_618$adj_r2 <- summary(topmod_Mg_618)$adj.r.squared
res_Mg_618$pvalue <- lmp(topmod_Mg_618)

res_S_618 <- as.data.frame(modsel_S_618[1])
res_S_618$adj_r2 <- summary(topmod_S_618)$adj.r.squared
res_S_618$pvalue <- lmp(topmod_S_618)

res_Zn_618 <- as.data.frame(modsel_Zn_618[1])
res_Zn_618$adj_r2 <- summary(topmod_Zn_618)$adj.r.squared
res_Zn_618$pvalue <- lmp(topmod_Zn_618)

res_618_spectal <- rbind(res_N_618, res_P_618, res_K_618, res_B_618, res_Ca_618, res_Mg_618, res_S_618, res_Zn_618)
res_618_spectal$dependent <- c("N", "P", "K", "B", "Ca", "Mg", "S", "Zn")
write.csv(res_618_spectal, "results/res_618_spectral.csv")

#---------------------------------------------------Repeat for 716 -----------------------------------------------

lm_N_716 <- lm(N ~ median+iqr+crrmedian+crriqr+skew+kurt+moran+rumple, data = agg_716, na.action=na.fail)
summary(lm_N_716)
vif(lm_N_716)
modsel_N_716 <- dredge(lm_N_716)
topmod_N_716 <- lm(N ~ kurt+rumple+skew+moran, data = agg_716,na.action=na.fail)
summary(topmod_N_716)
N_716_pred <- predict(topmod_N_716)
N_716_res <- resid(topmod_N_716)
n_716_table <- cbind(agg_716$N, N_716_pred, N_716_res)
head(n_716_table)

lm_P_716 <- lm(P ~ median+iqr+crrmedian+crriqr+skew+kurt+moran+rumple, data = agg_716, na.action=na.fail)
summary(lm_P_716)
modsel_P_716 <- dredge(lm_P_716)
topmod_P_716 <- lm(P ~ median+crrmedian, data = agg_716,na.action=na.fail)
summary(topmod_P_716)

lm_K_716 <- lm(K ~ median+iqr+crrmedian+crriqr+skew+kurt+moran+rumple, data = agg_716, na.action=na.fail)
summary(lm_K_716)
modsel_K_716 <- dredge(lm_K_716)
topmod_K_716 <- lm(K ~ crrmedian+median, data = agg_716,na.action=na.fail)
summary(topmod_K_716)

lm_Ca_716 <- lm(Ca ~ median+iqr+crrmedian+crriqr+skew+kurt+moran+rumple, data = agg_716, na.action=na.fail)
summary(lm_Ca_716)
modsel_Ca_716 <- dredge(lm_Ca_716)
topmod_Ca_716 <- lm(Ca ~ watersum, data = agg_716,na.action=na.fail)
summary(topmod_Ca_716)

lm_S_716 <- lm(S ~ median+iqr+crrmedian+crriqr+skew+kurt+moran+rumple, data = agg_716, na.action=na.fail)
summary(lm_S_716)
modsel_S_716 <- dredge(lm_S_716)
topmod_S_716 <- lm(S ~ 1, data = agg_716,na.action=na.fail)
summary(topmod_S_716)

lm_B_716 <- lm(B ~ median+iqr+crrmedian+crriqr+skew+kurt+moran+rumple, data = agg_716, na.action=na.fail)
summary(lm_B_716)
modsel_B_716 <- dredge(lm_B_716)
topmod_B_716 <- lm(B ~ crriqr+kurt, data = agg_716,na.action=na.fail)
summary(topmod_B_716)

lm_Zn_716 <- lm(Zn ~ median+iqr+crrmedian+crriqr+skew+kurt+moran+rumple, data = agg_716, na.action=na.fail)
summary(lm_Zn_716)
modsel_Zn_716 <- dredge(lm_Zn_716)
topmod_Zn_716 <- lm(Zn ~ crrmedian+kurt+median, data = agg_716,na.action=na.fail)
summary(topmod_Zn_716)

lm_Mg_716 <- lm(Mg ~ median+iqr+crrmedian+crriqr+skew+kurt+moran+rumple, data = agg_716, na.action=na.fail)
summary(lm_Mg_716)
modsel_Mg_716 <- dredge(lm_Mg_716)
topmod_Mg_716 <- lm(Mg ~ crrmedian+median+rumple, data = agg_716,na.action=na.fail)
summary(topmod_Mg_716)

# function to get lm pvalue
lmp <- function (modelobject) {
  if (class(modelobject) != "lm") stop("Not an object of class 'lm' ")
  f <- summary(modelobject)$fstatistic
  p <- pf(f[1],f[2],f[3],lower.tail=F)
  attributes(p) <- NULL
  return(p)
}

res_N_716 <- as.data.frame(modsel_N_716[1])
res_N_716$adj_r2 <- summary(topmod_N_716)$adj.r.squared
res_N_716$pvalue <- lmp(topmod_N_716)

res_P_716 <- as.data.frame(modsel_P_716[1])
res_P_716$adj_r2 <- summary(topmod_P_716)$adj.r.squared
res_P_716$pvalue <- lmp(topmod_P_716)

res_K_716 <- as.data.frame(modsel_K_716[1])
res_K_716$adj_r2 <- summary(topmod_K_716)$adj.r.squared
res_K_716$pvalue <- lmp(topmod_K_716)

res_B_716 <- as.data.frame(modsel_B_716[1])
res_B_716$adj_r2 <- summary(topmod_B_716)$adj.r.squared
res_B_716$pvalue <- lmp(topmod_B_716)

res_Ca_716 <- as.data.frame(modsel_Ca_716[1])
res_Ca_716$adj_r2 <- summary(topmod_Ca_716)$adj.r.squared
res_Ca_716$pvalue <- lmp(topmod_Ca_716)

res_Mg_716 <- as.data.frame(modsel_Mg_716[1])
res_Mg_716$adj_r2 <- summary(topmod_Mg_716)$adj.r.squared
res_Mg_716$pvalue <- lmp(topmod_Mg_716)

res_S_716 <- as.data.frame(modsel_S_716[1])
res_S_716$adj_r2 <- summary(topmod_S_716)$adj.r.squared
res_S_716$pvalue <- lmp(topmod_S_716)

res_Zn_716 <- as.data.frame(modsel_Zn_716[1])
res_Zn_716$adj_r2 <- summary(topmod_Zn_716)$adj.r.squared
res_Zn_716$pvalue <- lmp(topmod_Zn_716)

res_716 <- rbind(res_N_716, res_P_716, res_K_716, res_B_716, res_Ca_716, res_Mg_716, res_Zn_716)
res_716$dependent <- c("N", "P", "K", "B", "Ca", "Mg", "Zn")
write.csv(res_716, "results/res_716.csv")

#------------------------------------------- LM with spectral --------------------------------------------------------
lm_N_716 <- lm(N ~ median+iqr+crrmedian+crriqr+skew+kurt+moran+rumple+spad+varInd_median+tgi_median, data = agg_716, na.action=na.fail)
summary(lm_N_716)
vif(lm_N_716)
modsel_N_716 <- dredge(lm_N_716)
topmod_N_716 <- lm(N ~ crriqr+crrmedian+iqr+skew+spad+tgi_median, data = agg_716,na.action=na.fail)
summary(topmod_N_716)
N_716_pred <- predict(topmod_N_716)
N_716_res <- resid(topmod_N_716)
n_716_table <- cbind(agg_716$N, N_716_pred, N_716_res)
head(n_716_table)

lm_P_716 <- lm(P ~ median+iqr+crrmedian+crriqr+skew+kurt+moran+rumple+spad+varInd_median+tgi_median, data = agg_716, na.action=na.fail)
summary(lm_P_716)
modsel_P_716 <- dredge(lm_P_716)
topmod_P_716 <- lm(P ~ crrmedian+median+spad, data = agg_716,na.action=na.fail)
summary(topmod_P_716)

lm_K_716 <- lm(K ~ median+iqr+crrmedian+crriqr+skew+kurt+moran+rumple+spad+varInd_median+tgi_median, data = agg_716, na.action=na.fail)
summary(lm_K_716)
modsel_K_716 <- dredge(lm_K_716)
topmod_K_716 <- lm(K ~ crrmedian+median+spad, data = agg_716,na.action=na.fail)
summary(topmod_K_716)

lm_Ca_716 <- lm(Ca ~ median+iqr+crrmedian+crriqr+skew+kurt+moran+rumple+spad+varInd_median+tgi_median, data = agg_716, na.action=na.fail)
summary(lm_Ca_716)
modsel_Ca_716 <- dredge(lm_Ca_716)
topmod_Ca_716 <- lm(Ca ~ crrmedian+varInd_median, data = agg_716,na.action=na.fail)
summary(topmod_Ca_716)

lm_S_716 <- lm(S ~ median+iqr+crrmedian+crriqr+skew+kurt+moran+rumple+spad+varInd_median+tgi_median, data = agg_716, na.action=na.fail)
summary(lm_S_716)
modsel_S_716 <- dredge(lm_S_716)
topmod_S_716 <- lm(S ~ spad, data = agg_716,na.action=na.fail)
summary(topmod_S_716)

lm_B_716 <- lm(B ~ median+iqr+crrmedian+crriqr+skew+kurt+moran+rumple+spad+varInd_median+tgi_median, data = agg_716, na.action=na.fail)
summary(lm_B_716)
modsel_B_716 <- dredge(lm_B_716)
topmod_B_716 <- lm(B ~ crriqr+kurt, data = agg_716,na.action=na.fail)
summary(topmod_B_716)

lm_Zn_716 <- lm(Zn ~ median+iqr+crrmedian+crriqr+skew+kurt+moran+rumple+spad+varInd_median+tgi_median, data = agg_716, na.action=na.fail)
summary(lm_Zn_716)
modsel_Zn_716 <- dredge(lm_Zn_716)
topmod_Zn_716 <- lm(Zn ~ crrmedian+median+spad, data = agg_716,na.action=na.fail)
summary(topmod_Zn_716)

lm_Mg_716 <- lm(Mg ~ median+iqr+crrmedian+crriqr+skew+kurt+moran+rumple+spad+varInd_median+tgi_median, data = agg_716, na.action=na.fail)
summary(lm_Mg_716)
modsel_Mg_716 <- dredge(lm_Mg_716)
topmod_Mg_716 <- lm(Mg ~ crrmedian+varInd_median+spad+tgi_median, data = agg_716,na.action=na.fail)
summary(topmod_Mg_716)

res_N_716 <- as.data.frame(modsel_N_716[1])
res_N_716$adj_r2 <- summary(topmod_N_716)$adj.r.squared
res_N_716$pvalue <- lmp(topmod_N_716)

res_P_716 <- as.data.frame(modsel_P_716[1])
res_P_716$adj_r2 <- summary(topmod_P_716)$adj.r.squared
res_P_716$pvalue <- lmp(topmod_P_716)

res_K_716 <- as.data.frame(modsel_K_716[1])
res_K_716$adj_r2 <- summary(topmod_K_716)$adj.r.squared
res_K_716$pvalue <- lmp(topmod_K_716)

res_B_716 <- as.data.frame(modsel_B_716[1])
res_B_716$adj_r2 <- summary(topmod_B_716)$adj.r.squared
res_B_716$pvalue <- lmp(topmod_B_716)

res_Ca_716 <- as.data.frame(modsel_Ca_716[1])
res_Ca_716$adj_r2 <- summary(topmod_Ca_716)$adj.r.squared
res_Ca_716$pvalue <- lmp(topmod_Ca_716)

res_Mg_716 <- as.data.frame(modsel_Mg_716[1])
res_Mg_716$adj_r2 <- summary(topmod_Mg_716)$adj.r.squared
res_Mg_716$pvalue <- lmp(topmod_Mg_716)

res_S_716 <- as.data.frame(modsel_S_716[1])
res_S_716$adj_r2 <- summary(topmod_S_716)$adj.r.squared
res_S_716$pvalue <- lmp(topmod_S_716)

res_Zn_716 <- as.data.frame(modsel_Zn_716[1])
res_Zn_716$adj_r2 <- summary(topmod_Zn_716)$adj.r.squared
res_Zn_716$pvalue <- lmp(topmod_Zn_716)

res_716_spectal <- rbind(res_N_716, res_P_716, res_K_716, res_B_716, res_Ca_716, res_Mg_716, res_S_716, res_Zn_716)
res_716_spectal$dependent <- c("N", "P", "K", "B", "Ca", "Mg", "S", "Zn")
write.csv(res_716_spectal, "results/res_716_spectral.csv")

#---------------------------------------------------Repeat for 814 -----------------------------------------------

lm_N_814 <- lm(N ~ median+iqr+crrmedian+crriqr+skew+kurt+moran+rumple, data = agg_814, na.action=na.fail)
summary(lm_N_814)
vif(lm_N_814)
modsel_N_814 <- dredge(lm_N_814)
topmod_N_814 <- lm(N ~ median+rumple, data = agg_814,na.action=na.fail)
summary(topmod_N_814)
N_814_pred <- predict(topmod_N_814)
N_814_res <- resid(topmod_N_814)
n_814_table <- cbind(agg_814$N, N_814_pred, N_814_res)
head(n_814_table)

lm_P_814 <- lm(P ~ median+iqr+crrmedian+crriqr+skew+kurt+moran+rumple, data = agg_814, na.action=na.fail)
summary(lm_P_814)
modsel_P_814 <- dredge(lm_P_814)
topmod_P_814 <- lm(P ~ iqr+kurt, data = agg_814,na.action=na.fail)
summary(topmod_P_814)

lm_K_814 <- lm(K ~ median+iqr+crrmedian+crriqr+skew+kurt+moran+rumple, data = agg_814, na.action=na.fail)
summary(lm_K_814)
modsel_K_814 <- dredge(lm_K_814)
topmod_K_814 <- lm(K ~ crrmedian+crriqr, data = agg_814,na.action=na.fail)
summary(topmod_K_814)

lm_Ca_814 <- lm(Ca ~ median+iqr+crrmedian+crriqr+skew+kurt+moran+rumple, data = agg_814, na.action=na.fail)
summary(lm_Ca_814)
modsel_Ca_814 <- dredge(lm_Ca_814)
topmod_Ca_814 <- lm(Ca ~ crrmedian+kurt+skew, data = agg_814,na.action=na.fail)
summary(topmod_Ca_814)

lm_S_814 <- lm(S ~ median+iqr+crrmedian+crriqr+skew+kurt+moran+rumple, data = agg_814, na.action=na.fail)
summary(lm_S_814)
modsel_S_814 <- dredge(lm_S_814)
topmod_S_814 <- lm(S ~ crriqr+median, data = agg_814,na.action=na.fail)
summary(topmod_S_814)

lm_B_814 <- lm(B ~ median+iqr+crrmedian+crriqr+skew+kurt+moran+rumple, data = agg_814, na.action=na.fail)
summary(lm_B_814)
modsel_B_814 <- dredge(lm_B_814)
topmod_B_814 <- lm(B ~ crriqr+moran, data = agg_814,na.action=na.fail)
summary(topmod_B_814)

lm_Zn_814 <- lm(Zn ~ median+iqr+crrmedian+crriqr+skew+kurt+moran+rumple, data = agg_814, na.action=na.fail)
summary(lm_Zn_814)
modsel_Zn_814 <- dredge(lm_Zn_814)
topmod_Zn_814 <- lm(Zn ~ 1, data = agg_814,na.action=na.fail)
summary(topmod_Zn_814)

lm_Mg_814 <- lm(Mg ~ median+iqr+crrmedian+crriqr+skew+kurt+moran+rumple, data = agg_814, na.action=na.fail)
summary(lm_Mg_814)
modsel_Mg_814 <- dredge(lm_Mg_814)
topmod_Mg_814 <- lm(Mg ~ crrmedian+kurt+skew, data = agg_814,na.action=na.fail)
summary(topmod_Mg_814)

# function to get lm pvalue
lmp <- function (modelobject) {
  if (class(modelobject) != "lm") stop("Not an object of class 'lm' ")
  f <- summary(modelobject)$fstatistic
  p <- pf(f[1],f[2],f[3],lower.tail=F)
  attributes(p) <- NULL
  return(p)
}

res_N_814 <- as.data.frame(modsel_N_814[1])
res_N_814$adj_r2 <- summary(topmod_N_814)$adj.r.squared
res_N_814$pvalue <- lmp(topmod_N_814)

res_P_814 <- as.data.frame(modsel_P_814[1])
res_P_814$adj_r2 <- summary(topmod_P_814)$adj.r.squared
res_P_814$pvalue <- lmp(topmod_P_814)

res_K_814 <- as.data.frame(modsel_K_814[1])
res_K_814$adj_r2 <- summary(topmod_K_814)$adj.r.squared
res_K_814$pvalue <- lmp(topmod_K_814)

res_B_814 <- as.data.frame(modsel_B_814[1])
res_B_814$adj_r2 <- summary(topmod_B_814)$adj.r.squared
res_B_814$pvalue <- lmp(topmod_B_814)

res_Ca_814 <- as.data.frame(modsel_Ca_814[1])
res_Ca_814$adj_r2 <- summary(topmod_Ca_814)$adj.r.squared
res_Ca_814$pvalue <- lmp(topmod_Ca_814)

res_Mg_814 <- as.data.frame(modsel_Mg_814[1])
res_Mg_814$adj_r2 <- summary(topmod_Mg_814)$adj.r.squared
res_Mg_814$pvalue <- lmp(topmod_Mg_814)

res_S_814 <- as.data.frame(modsel_S_814[1])
res_S_814$adj_r2 <- summary(topmod_S_814)$adj.r.squared
res_S_814$pvalue <- lmp(topmod_S_814)

res_Zn_814 <- as.data.frame(modsel_Zn_814[1])
res_Zn_814$adj_r2 <- summary(topmod_Zn_814)$adj.r.squared
res_Zn_814$pvalue <- lmp(topmod_Zn_814)

res_814 <- rbind(res_N_814, res_P_814, res_K_814, res_B_814, res_Ca_814, res_Mg_814, res_S_814)
res_814$dependent <- c("N", "P", "K", "B", "Ca", "Mg", "S")
write.csv(res_814, "results/res_814.csv")

#------------------------------------------- LM with spectral --------------------------------------------------------

lm_N_814 <- lm(N ~ median+iqr+crrmedian+crriqr+skew+kurt+moran+rumple+spad+varInd_median+tgi_median, data = agg_814, na.action=na.fail)
summary(lm_N_814)
vif(lm_N_814)
modsel_N_814 <- dredge(lm_N_814)
topmod_N_814 <- lm(N ~ spad+varInd_median, data = agg_814,na.action=na.fail)
summary(topmod_N_814)
N_814_pred <- predict(topmod_N_814)
N_814_res <- resid(topmod_N_814)
n_814_table <- cbind(agg_814$N, N_814_pred, N_814_res)
head(n_814_table)

lm_P_814 <- lm(P ~ median+iqr+crrmedian+crriqr+skew+kurt+moran+rumple+spad+varInd_median+tgi_median, data = agg_814, na.action=na.fail)
summary(lm_P_814)
modsel_P_814 <- dredge(lm_P_814)
topmod_P_814 <- lm(P ~ iqr+kurt+spad, data = agg_814,na.action=na.fail)
summary(topmod_P_814)

lm_K_814 <- lm(K ~ median+iqr+crrmedian+crriqr+skew+kurt+moran+rumple+spad+varInd_median+tgi_median, data = agg_814, na.action=na.fail)
summary(lm_K_814)
modsel_K_814 <- dredge(lm_K_814)
topmod_K_814 <- lm(K ~ crrmedian+spad, data = agg_814,na.action=na.fail)
summary(topmod_K_814)

lm_Ca_814 <- lm(Ca ~ median+iqr+crrmedian+crriqr+skew+kurt+moran+rumple+spad+varInd_median+tgi_median, data = agg_814, na.action=na.fail)
summary(lm_Ca_814)
modsel_Ca_814 <- dredge(lm_Ca_814)
topmod_Ca_814 <- lm(Ca ~ crrmedian+kurt+skew, data = agg_814,na.action=na.fail)
summary(topmod_Ca_814)

lm_S_814 <- lm(S ~ median+iqr+crrmedian+crriqr+skew+kurt+moran+rumple+spad+varInd_median+tgi_median, data = agg_814, na.action=na.fail)
summary(lm_S_814)
modsel_S_814 <- dredge(lm_S_814)
topmod_S_814 <- lm(S ~ crriqr+median+spad+varInd_median, data = agg_814,na.action=na.fail)
summary(topmod_S_814)

lm_B_814 <- lm(B ~ median+iqr+crrmedian+crriqr+skew+kurt+moran+rumple+spad+varInd_median+tgi_median, data = agg_814, na.action=na.fail)
summary(lm_B_814)
modsel_B_814 <- dredge(lm_B_814)
topmod_B_814 <- lm(B ~ crriqr+kurt, data = agg_814,na.action=na.fail)
summary(topmod_B_814)

lm_Zn_814 <- lm(Zn ~ median+iqr+crrmedian+crriqr+skew+kurt+moran+rumple+spad+varInd_median+tgi_median, data = agg_814, na.action=na.fail)
summary(lm_Zn_814)
modsel_Zn_814 <- dredge(lm_Zn_814)
topmod_Zn_814 <- lm(Zn ~ spad, data = agg_814,na.action=na.fail)
summary(topmod_Zn_814)

lm_Mg_814 <- lm(Mg ~ median+iqr+crrmedian+crriqr+skew+kurt+moran+rumple+spad+varInd_median+tgi_median, data = agg_814, na.action=na.fail)
summary(lm_Mg_814)
modsel_Mg_814 <- dredge(lm_Mg_814)
topmod_Mg_814 <- lm(Mg ~ crrmedian+varInd_median, data = agg_814,na.action=na.fail)
summary(topmod_Mg_814)

res_N_814 <- as.data.frame(modsel_N_814[1])
res_N_814$adj_r2 <- summary(topmod_N_814)$adj.r.squared
res_N_814$pvalue <- lmp(topmod_N_814)

res_P_814 <- as.data.frame(modsel_P_814[1])
res_P_814$adj_r2 <- summary(topmod_P_814)$adj.r.squared
res_P_814$pvalue <- lmp(topmod_P_814)

res_K_814 <- as.data.frame(modsel_K_814[1])
res_K_814$adj_r2 <- summary(topmod_K_814)$adj.r.squared
res_K_814$pvalue <- lmp(topmod_K_814)

res_B_814 <- as.data.frame(modsel_B_814[1])
res_B_814$adj_r2 <- summary(topmod_B_814)$adj.r.squared
res_B_814$pvalue <- lmp(topmod_B_814)

res_Ca_814 <- as.data.frame(modsel_Ca_814[1])
res_Ca_814$adj_r2 <- summary(topmod_Ca_814)$adj.r.squared
res_Ca_814$pvalue <- lmp(topmod_Ca_814)

res_Mg_814 <- as.data.frame(modsel_Mg_814[1])
res_Mg_814$adj_r2 <- summary(topmod_Mg_814)$adj.r.squared
res_Mg_814$pvalue <- lmp(topmod_Mg_814)

res_S_814 <- as.data.frame(modsel_S_814[1])
res_S_814$adj_r2 <- summary(topmod_S_814)$adj.r.squared
res_S_814$pvalue <- lmp(topmod_S_814)

res_Zn_814 <- as.data.frame(modsel_Zn_814[1])
res_Zn_814$adj_r2 <- summary(topmod_Zn_814)$adj.r.squared
res_Zn_814$pvalue <- lmp(topmod_Zn_814)

res_814_spectal <- rbind(res_N_814, res_P_814, res_K_814, res_B_814, res_Ca_814, res_Mg_814, res_S_814, res_Zn_814)
res_814_spectal$dependent <- c("N", "P", "K", "B", "Ca", "Mg", "S", "Zn")
write.csv(res_814_spectal, "results/res_814_spectal.csv")




# ------------------------------ Try lm with nutrients as independent ---------------------------------------------

# # All nutrients: 
# #all nutrients: `N`+`P`+`K`+`Mg`+`Ca`+`S`+`Zn`+`Mn`+`Cu`+`Fe`+`B`+`Al`+`Na`
# 
# # Function for calculating r^2 in SAR model, null model will be same for each date
# 
# null <- spautolm(median_ln ~ 1, data = std_618, listw=nb2listw(w), zero.policy = T, na.action = na.fail, family="SAR")
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
# # 618 OLS
# lm_618 <- lm(median ~ `N`+`K`+`Mg`+`Ca`+`Mn`+`Cu`+`Fe`+`B`+`Al`+`Na`, data = plots_618, na.action = na.fail)
# summary(lm_618)
# qqnorm(lm_618$residuals)
# sd(lm_618$residuals)
# vif(lm_618) # Colinearity: removed Zn, P, and S for VIF > 4
# dredge(lm_618, trace = T)
# 
# # 618 check for spatial autocorrelation
# std_618$ols_res <- lm_618$residuals
# # Weights matrix
# w <- knn2nb(knearneigh(coordinates(plots), k=8))
# moran.test(std_618$ols_res, nb2listw(w))
# 
# #Apply SAR model
# sar_618 <- spautolm(median ~ `N`+`K`+`Mg`+`Ca`+`Mn`+`Cu`+`Fe`+`B`+`Al`+`Na`, data = std_618, listw=nb2listw(w), zero.policy = T, na.action = na.fail, family="SAR")
# summary(sar_618)
# sar_resi <- sar_618$fit[[9]]
# std_618$sar_resi <- sar_resi
# moran.test(std_618$sar_resi, nb2listw(w))
# 
# 
# lm_716 <- lm(median ~ `N`+`K`+`Mg`+`Ca`+`Mn`+`Cu`+`Fe`+`B`+`Al`+`Na`, data = std_716, na.action = na.fail)
# summary(lm_716)
# vif(lm_716)
# dredge(lm_716, trace = T)
# 
# lm_814 <- lm(rumple ~ `N` + `P` + `K`+`Mg` +`Ca` + `S` + `Cu` +`Fe`+`B`+`Al`+`Na`, data = std_814, na.action = na.fail)
# summary(lm_814)
# vif(lm_814)
# dredge(lm_814, trace = T)


# ---------------------------------SAR models with structural metric as response for each time step--------------------------------

# variables <- "`N`+`K`+`Mg`+`Ca`+`Mn`+`Cu`+`Fe`+`B`+`Al`+`Na`"
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
# w <- knn2nb(knearneigh(coordinates(plots_618), k=8))
# 
# # ---------------------------------------SAR 618--------------------------------------------
# dat <- std_618
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
# med_618 <- modsel[1]
# med_618$depvar <- "median"
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
# iqr_618 <- modsel[1]
# iqr_618$depvar <- "iqr"
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
# crrmed_618 <- modsel[1]
# crrmed_618$depvar <- "crrmed"
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
# crriqr_618 <- modsel[1]
# crriqr_618$depvar <- "crriqr"
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
# skew_618 <- modsel[1]
# skew_618$depvar <- "skew"
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
# kurt_618 <- modsel[1]
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
# plotcrr_618 <- modsel[1]
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
# rumple_618 <- modsel[1]
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
# moran_618 <- modsel[1]
# 
# 
# # ---------------------------------------SAR 716--------------------------------------------
# dat <- std_716
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
# med_716 <- modsel[1]
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
# iqr_716 <- modsel[1]
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
# crrmed_716 <- modsel[1]
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
# crriqr_716 <- modsel[1]
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
# skew_716 <- modsel[1]
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
# kurt_716 <- modsel[1]
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
# plotcrr_716 <- modsel[1]
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
# rumple_716 <- modsel[1]
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
# moran_716 <- modsel[1]
# 
# # ---------------------------------------SAR 814--------------------------------------------
# dat <- std_814
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
# med_814 <- modsel[1]
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
# iqr_814 <- modsel[1]
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
# crrmed_814 <- modsel[1]
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
# crriqr_814 <- modsel[1]
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
# skew_814 <- modsel[1]
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
# kurt_814 <- modsel[1]
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
# plotcrr_814 <- modsel[1]
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
# rumple_814 <- modsel[1]
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
# moran_814 <- modsel[1]
# 
# #---------------------------------------SAR nutrient results tables --------------------------------------
# 
# med_618$depvar <- "median"
# iqr_618$depvar <- "iqr"
# crrmed_618$depvar <- "crrmed"
# crriqr_618$depvar <- "crriqr"
# skew_618$depvar <- "skew"
# kurt_618$depvar <- "kurt"
# plotcrr_618$depvar <- "plotcrr"
# rumple_618$depvar <- "rumple"
# moran_618$depvar <- "moran"
# 
# med_716$depvar <- "median"
# iqr_716$depvar <- "iqr"
# crrmed_716$depvar <- "crrmed"
# crriqr_716$depvar <- "crriqr"
# skew_716$depvar <- "skew"
# kurt_716$depvar <- "kurt"
# plotcrr_716$depvar <- "plotcrr"
# rumple_716$depvar <- "rumple"
# moran_716$depvar <- "moran"
# 
# med_814$depvar <- "median"
# iqr_814$depvar <- "iqr"
# crrmed_814$depvar <- "crrmed"
# crriqr_814$depvar <- "crriqr"
# skew_814$depvar <- "skew"
# kurt_814$depvar <- "kurt"
# plotcrr_814$depvar <- "plotcrr"
# rumple_814$depvar <- "rumple"
# moran_814$depvar <- "moran"
# 
# sarresults_nutr_618 <- rbind(med_618, iqr_618, crrmed_618, crriqr_618, skew_618, kurt_618, plotcrr_618, rumple_618, moran_618)
# sarresults_nutr_618$depvar <- c("median", "iqr", "crrmed", "crriqr", "skew", "kurt", "plotcrr", "rumple", "moran")
# 
# sarresults_nutr_716 <- rbind(med_716, iqr_716, crrmed_716, crriqr_716, skew_716, kurt_716, plotcrr_716, rumple_716, moran_716)
# sarresults_nutr_716$depvar <- c("median", "iqr", "crrmed", "crriqr", "skew", "kurt", "plotcrr", "rumple", "moran")
# 
# sarresults_nutr_814 <- rbind(med_814, iqr_814, crrmed_814, crriqr_814, skew_814, kurt_814, plotcrr_814, rumple_814, moran_814)
# sarresults_nutr_814$depvar <- c("median", "iqr", "crrmed", "crriqr", "skew", "kurt", "plotcrr", "rumple", "moran")
# 
# write.csv(sarresults_nutr_618, "Q:/My Drive/Research/Canopy_Morphology/Tobacco_Project/2018/canopy_analysis/sarresults_nutr_618.csv")
# write.csv(sarresults_nutr_716, "Q:/My Drive/Research/Canopy_Morphology/Tobacco_Project/2018/canopy_analysis/sarresults_nutr_716.csv")
# write.csv(sarresults_nutr_814, "Q:/My Drive/Research/Canopy_Morphology/Tobacco_Project/2018/canopy_analysis/sarresults_nutr_814.csv")
# 
# 
# # ---------------------------------SAR models with nutrients as response for each time step--------------------------------
# 
# variables <- "crrmedian+crriqr+median+kurt+rumple+moran"
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
# w <- knn2nb(knearneigh(coordinates(plots_618), k=8))
# 
# # ---------------------------------------SAR 618--------------------------------------------
# dat <- std_618
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
# n_618 <- modsel[1]
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
# p_618 <- modsel[1]
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
# k_618 <- modsel[1]
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
# ca_618 <- modsel[1]
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
# mg_618 <- modsel[1]
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
# b_618 <- modsel[1]
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
# mn_618 <- modsel[1]
# 
# # ---------------------------------------SAR 716--------------------------------------------
# dat <- std_716
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
# n_716 <- modsel[1]
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
# p_716 <- modsel[1]
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
# k_716 <- modsel[1]
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
# ca_716 <- modsel[1]
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
# mg_716 <- modsel[1]
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
# b_716 <- modsel[1]
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
# mn_716 <- modsel[1]
# 
# # ---------------------------------------SAR 814--------------------------------------------
# dat <- std_814
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
# n_814 <- modsel[1]
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
# p_814 <- modsel[1]
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
# k_814 <- modsel[1]
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
# ca_814 <- modsel[1]
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
# mg_814 <- modsel[1]
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
# b_814 <- modsel[1]
# 
# metric <- "Mn"
# null <- spautolm( as.formula(paste(metric, 1, sep="~")), data = dat, 
#                   listw=nb2listw(w), zero.policy = T, na.action = na.fail, family="SAR")
# sar <- spautolm(as.formula(paste(metric, variables, sep = " ~ ")), 
#                 data = dat, listw=nb2listw(w), zero.policy = T, na.action = na.fail, 
#                 family="SAR")
# tmp <- spautolm(B ~ moran, data = std_618, listw=nb2listw(w), zero.policy = T, na.action = na.fail, 
#          family="SAR")
# 
# summary(sar)
# adj_r2(sar)
# modsel <- dredge(sar, beta = "none", trace = T, extra = c('r2', 'adj_r2'))
# mn_814 <- modsel[1]
# 
# # --------------------------------------------SAR structure results tables -------------------------------------
# sarresults_struct_618 <- rbind(n_618, p_618, k_618, ca_618, mg_618, b_618, mn_618)
# sarresults_struct_618$depvar <- c("N", "P", "K", "Ca", "Mg", "B", "Mn")
# 
# sarresults_struct_716 <- rbind(n_716, p_716, k_716, ca_716, mg_716, b_716, mn_716)
# sarresults_struct_716$depvar <- c("N", "P", "K", "Ca", "Mg", "B", "Mn")
# 
# sarresults_struct_814 <- rbind(n_814, p_814, k_814, ca_814, mg_814, b_814, mn_814)
# sarresults_struct_814$depvar <- c("N", "P", "K", "Ca", "Mg", "B", "Mn")
# 
# write.csv(sarresults_struct_618, "Q:/My Drive/Research/Canopy_Morphology/Tobacco_Project/2018/canopy_analysis/sarresults_struct_618.csv")
# write.csv(sarresults_struct_716, "Q:/My Drive/Research/Canopy_Morphology/Tobacco_Project/2018/canopy_analysis/sarresults_struct_716.csv")
# write.csv(sarresults_struct_814, "Q:/My Drive/Research/Canopy_Morphology/Tobacco_Project/2018/canopy_analysis/sarresults_struct_814.csv")
# 
# # ---------------------------------------- Try RF regression --------------------------------------------------
# 
# #all nutrients: `N` + `P` + `K`+`Mg` +`Ca` + `S` + `Zn`+`Mn`+`Cu`+`Fe`+`B`+`Al`+`Na`
# 
# #--------------------------------RF 618 ---------------------------------
# 
# set.seed(42)
# train_data_indices <- rep(FALSE, nrow(std_618))
# train_data_indices[sample(1:nrow(std_618), round(0.7 * nrow(std_618)))] <- TRUE # randomly select 70% of the data for training
# rf_regression<- randomForest(skew ~ `N` + `P` + `K`+`Mg` +`Ca` + `Zn`+`Cu`+`Fe`+`B`+`Al`+`Na`, data = std_618[train_data_indices, ], importance=T)
# rf_regression
# varImpPlot(rf_regression)
# pred_struc <- predict(rf_regression, std_618[!train_data_indices,]) # predict the rings
# #table(pred_struc, std_618$median[!train_data_indices])
# plot(std_618$skew[!train_data_indices], pred_struc, xlab="Observed", ylab="Predicted")
# 
# #--------------------------------RF 716 ---------------------------------
# 
# set.seed(42)
# train_data_indices <- rep(FALSE, nrow(std_716))
# train_data_indices[sample(1:nrow(std_716), round(0.7 * nrow(std_716)))] <- TRUE # randomly select 70% of the data for training
# rf_regression<- randomForest(skew ~ `N` + `P` + `K`+`Mg` +`Ca` + `Zn`+`Cu`+`Fe`+`B`+`Al`+`Na`, data = std_716[train_data_indices, ], importance=T)
# rf_regression
# varImpPlot(rf_regression)
# pred_struc <- predict(rf_regression, std_716[!train_data_indices,]) # predict the rings
# #table(pred_struc, std_716$median[!train_data_indices])
# plot(std_716$skew[!train_data_indices], pred_struc, xlab="Observed", ylab="Predicted")
# 
# #--------------------------------RF 814 ---------------------------------
# 
# set.seed(42)
# train_data_indices <- rep(FALSE, nrow(std_716))
# train_data_indices[sample(1:nrow(std_814), round(0.7 * nrow(std_814)))] <- TRUE # randomly select 70% of the data for training
# rf_regression<- randomForest(median ~ `N` + `P` + `K`+`Mg` +`Ca` + `Zn`+`Cu`+`Fe`+`B`+`Al`+`Na`, data = std_814[train_data_indices, ], importance=T)
# rf_regression
# varImpPlot(rf_regression)
# pred_struc <- predict(rf_regression, std_814[!train_data_indices,]) # predict the rings
# #table(pred_struc, std_814$median[!train_data_indices])
# plot(std_814$median[!train_data_indices], pred_struc, xlab="Observed", ylab="Predicted")
# 
# 
# # ---------------------------------------Try multinominal logisitic regression with nnet package --------------------------------
# 
# # 618
# 
# library(nnet)
# plots_618_df$Trt_agg <- relevel(plots_618_df$Trt_agg, ref = "100% Control")
# mltnom_618 <- multinom(Trt_agg ~ crrmedian+median+kurt+rumple+moran+dip+crriqr, data = plots_618_df, na.action = na.fail)
# summary(mltnom_618)
# 
# modsel_618 <- dredge(mltnom_618)
# modsel_618
# 
# topmod_618 <- multinom(Trt_agg ~ crrmedian+median+kurt+rumple+moran, data = plots_618_df, na.action = na.fail)
# summary(topmod_618)
# 
# library(DescTools)
# PseudoR2(topmod_618, which = "all")
# z <- summary(topmod_618)$coefficients/summary(topmod_618)$standard.errors
# z
# p <- (1 - pnorm(abs(z), 0, 1)) * 2
# p
# 
# coefs_618 <- coef(topmod_618)
# exp(coefs_618) #transformed to odds ratio
# (exp(coefs_618)-1)*100 # percent change in the odds for a one unit increase in the independent variable

# Try binary logistic regression with 100% control as baseline
# def_100 <- as.data.frame(dt_618[Trt_agg == "100% Control" | Trt_agg == "deficiency",])
# def_100 <- droplevels(def_100)
# tox_100 <- as.data.frame(dt_618[Trt_agg == "100% Control" | Trt_agg == "toxicity",])
# tox_100 <- droplevels(tox_100)
# control_0_100 <- as.data.frame(dt_618[Trt_agg == "100% Control" | Trt_agg == "0% Control",])
# control_0_100 <- droplevels(control_0_100)
# 
# binary_618 <- glm(Trt_agg ~ crrmedian+median+kurt+rumple+moran+crriqr+skew+iqr, data = control_0_100, family = binomial, na.action = na.fail)
# modsel_618_binary <- dredge(binary_618)
# topmod_618_binary <- glm(Trt_agg ~ crrmedian+iqr+moran+rumple+skew, data = tox_100, family = binomial)
# summary(topmod_618_binary)
# 
# library(pscl)
# 
# pR2(topmod_618_binary)

# Try multinomial logisitc regression with mlogit package
# library(foreign)
# library(mlogit)
# 
# mlogit_data_618 <- mlogit.data(plots_618_df, choice = "Trt_agg", shape = "wide")
# mlogit_618 <- mlogit(Trt_agg~crrmedian+kurt+median+skew_ln+rumple+moran, data = mlogit_data_618,
#                      method = "nr", print.level = 0)
# summary(mlogit_618)
# 
# # Weights matrix
# w <- knn2nb(knearneigh(coordinates(plots), k=8))

# 
# # 716
# plots_716_df$Trt_agg <- relevel(plots_716_df$Trt_agg, ref = "100% Control")
# mltnom_716 <- multinom(Trt_agg ~ crrmedian+median+kurt+rumple+moran+dip+crriqr, data = plots_716_df, na.action = na.fail)
# summary(mltnom_716)
# 
# modsel_716 <- dredge(mltnom_716)
# modsel_716
# 
# topmod_716 <- multinom(Trt_agg ~ crrmedian+kurt+median+moran, data = plots_716_df, na.action = na.fail)
# summary(topmod_716)
# PseudoR2(topmod_716, which = "all")
# z <- summary(topmod_716)$coefficients/summary(topmod_716)$standard.errors
# z
# p <- (1 - pnorm(abs(z), 0, 1)) * 2
# p
# 
# coefs_716 <- coef(topmod_716)
# exp(coefs_716) #transformed to odds ratio
# (exp(coefs_716)-1)*100 # percent change in the odds for a one unit increase in the independent variable
# 
# # 814
# plots_814_df$Trt_agg <- relevel(plots_814_df$Trt_agg, ref = "100% Control")
# mltnom_814 <- multinom(Trt_agg ~ crrmedian+median+kurt+rumple+moran+dip+crriqr, data = plots_814_df, na.action = na.fail)
# summary(mltnom_814)
# 
# modsel_814 <- dredge(mltnom_814)
# modsel_814
# 
# topmod_814 <- multinom(Trt_agg ~ median+rumple+moran, data = plots_814_df, na.action = na.fail)
# summary(topmod_814)
# PseudoR2(topmod_814, which = "all")
# z <- summary(topmod_814)$coefficients/summary(topmod_814)$standard.errors
# z
# p <- (1 - pnorm(abs(z), 0, 1)) * 2
# p
# 
# coefs_814 <- coef(topmod_814)
# exp(coefs_814) #transformed to odds ratio
# (exp(coefs_814)-1)*100 # percent change in the odds for a one unit increase in the independent variable
# 
# 
# # Try combining tox and def to see if model improves
# plots_814_df$Trt_agg_v2 <- as.character(plots_618_df$Trt_agg)
# plots_814_df[Trt_agg == "deficiency"]$Trt_agg_v2 <- "stressed"
# plots_814_df[Trt_agg == "toxicity"]$Trt_agg_v2 <- "stressed"
# plots_814_df$Trt_agg_v2 <- as.factor(plots_814_df$Trt_agg_v2) 
# 
# plots_814_df$Trt_agg <- relevel(plots_814_df$Trt_agg, ref = "100% Control")
# 
# mltnom_814_v2 <- multinom(Trt_agg_v2 ~ crrmedian+median+kurt+rumple+moran+dip+crriqr, data = plots_814_df, na.action = na.fail)
# summary(mltnom_814_v2)
# 
# modsel_814_v2 <- dredge(mltnom_814_v2)
# modsel_814_v2
# 
# topmod_814_v2 <- multinom(Trt_agg_v2 ~ median+rumple+moran+crrmedian, data = plots_814_df, na.action = na.fail)
# summary(topmod_814_v2)
# PseudoR2(topmod_814_v2, which = "all")
# z <- summary(topmod_814)$coefficients/summary(topmod_814)$standard.errors
# z
# p <- (1 - pnorm(abs(z), 0, 1)) * 2
# p
# 
# 
# # ----------------------------------------------Result tables------------------------------------------------------
# 
# vars <- topmod_618$coefnames[-1]
# sd_all <- numeric(length(vars))
# mean_all <- numeric(length(vars))
# count <- 1
# for (x in vars){
#   values <- as.data.frame(plots_618_df)[x]
#   sd_all[count] <- sd(values[,1])
#   mean_all[count] <- mean(values[,1])
#   count <- count+1
# }
# 
# vars_618 <- as.data.frame(vars)
# vars_618$mean <- mean_all
# vars_618$sd <- sd_all
# 
# #write.csv(vars_618, "Q:/My Drive/Research/Canopy_Morphology/Tobacco_Project/canopy_analysis/tobacco_results/vars_618.csv")
# 
# levels <- topmod_618$lev[-1]
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
#   coefs[count] <- list(coef(topmod_618)[count,])
#   SE[count] <- list(summary(topmod_618)$standard.errors[count,])
#   z <- (summary(topmod_618)$coefficients/summary(topmod_618)$standard.errors)
#   Wald[count] <- list(z[count,])
#   p[count] <- list(((1 - pnorm(abs(z), 0, 1)) * 2)[count,])
#   std_coefs[count] <- list((sqrt(3)/pi) * (coefs[[count]][-1] * sd_all))
#   exp[count] <- list(exp(std_coefs[[count]]))
#   count <- count+1
# }
# 
# vars <- topmod_618$coefnames
# 
# results_0_618 <- as.data.frame(vars)
# results_0_618$level <- "0% Control"
# results_0_618$coefs <- coefs[[1]]
# results_0_618$std_coefs <- c(0,std_coefs[[1]])
# results_0_618$exp <- c(0,exp[[1]])
# results_0_618$SE <- SE[[1]]
# results_0_618$p <- p[[1]]
# results_0_618$wald <- Wald[[1]]
# 
# results_tox_618 <- as.data.frame(vars)
# results_tox_618$level <- "toxicity"
# results_tox_618$coefs <- coefs[[3]]
# results_tox_618$std_coefs <- c(0,std_coefs[[3]])
# results_tox_618$exp <- c(0,exp[[3]])
# results_tox_618$SE <- SE[[3]]
# results_tox_618$p <- p[[3]]
# results_tox_618$wald <- Wald[[3]]
# 
# results_def_618 <- as.data.frame(vars)
# results_def_618$level <- "deficiency"
# results_def_618$coefs <- coefs[[2]]
# results_def_618$std_coefs <- c(0,std_coefs[[2]])
# results_def_618$exp <- c(0,exp[[2]])
# results_def_618$SE <- SE[[2]]
# results_def_618$p <- p[[2]]
# results_def_618$wald <- Wald[[2]]
# 
# 
# 
# results_618 <- rbind(results_0_618, results_def_618, results_tox_618)
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
# #     Wald[count] <- list(summary(x)$coefficients/summary(topmod_618)$standard.errors[count,])
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
# # mlogit_results_618 <- mlogit_results(topmod_618)
# 
# 
# #write.csv(results_618, "Q:/My Drive/Research/Canopy_Morphology/Tobacco_Project/canopy_analysis/tobacco_results/results_618.csv")
# 
# 
# vars <- topmod_716$coefnames[-1]
# sd_all <- numeric(length(vars))
# mean_all <- numeric(length(vars))
# count <- 1
# for (x in vars){
#   values <- as.data.frame(plots_716_df)[x]
#   sd_all[count] <- sd(values[,1])
#   mean_all[count] <- mean(values[,1])
#   count <- count+1
# }
# 
# vars_716 <- as.data.frame(vars)
# vars_716$mean <- mean_all
# vars_716$sd <- sd_all
# 
# #write.csv(vars_716, "Q:/My Drive/Research/Canopy_Morphology/Tobacco_Project/canopy_analysis/tobacco_results/vars_716.csv")
# 
# levels <- topmod_716$lev[-1]
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
#   coefs[count] <- list(coef(topmod_716)[count,])
#   SE[count] <- list(summary(topmod_716)$standard.errors[count,])
#   z <- (summary(topmod_716)$coefficients/summary(topmod_716)$standard.errors)
#   Wald[count] <- list(z[count,])
#   p[count] <- list(((1 - pnorm(abs(z), 0, 1)) * 2)[count,])
#   std_coefs[count] <- list((sqrt(3)/pi) * (coefs[[count]][-1] * sd_all))
#   exp[count] <- list(exp(std_coefs[[count]]))
#   count <- count+1
# }
# 
# vars <- topmod_716$coefnames
# 
# results_0_716 <- as.data.frame(vars)
# results_0_716$level <- "0% Control"
# results_0_716$coefs <- coefs[[1]]
# results_0_716$std_coefs <- c(0,std_coefs[[1]])
# results_0_716$exp <- c(0,exp[[1]])
# results_0_716$SE <- SE[[1]]
# results_0_716$p <- p[[1]]
# results_0_716$wald <- Wald[[1]]
# 
# results_tox_716 <- as.data.frame(vars)
# results_tox_716$level <- "toxicity"
# results_tox_716$coefs <- coefs[[3]]
# results_tox_716$std_coefs <- c(0,std_coefs[[3]])
# results_tox_716$exp <- c(0,exp[[3]])
# results_tox_716$SE <- SE[[3]]
# results_tox_716$p <- p[[3]]
# results_tox_716$wald <- Wald[[3]]
# 
# results_def_716 <- as.data.frame(vars)
# results_def_716$level <- "deficiency"
# results_def_716$coefs <- coefs[[2]]
# results_def_716$std_coefs <- c(0,std_coefs[[2]])
# results_def_716$exp <- c(0,exp[[2]])
# results_def_716$SE <- SE[[2]]
# results_def_716$p <- p[[2]]
# results_def_716$wald <- Wald[[2]]
# 
# 
# results_716 <- rbind(results_0_716, results_def_716, results_tox_716)
# 
# #write.csv(results_716, "Q:/My Drive/Research/Canopy_Morphology/Tobacco_Project/canopy_analysis/tobacco_results/results_716.csv")
# 
# 
# vars <- topmod_814$coefnames[-1]
# sd_all <- numeric(length(vars))
# mean_all <- numeric(length(vars))
# count <- 1
# for (x in vars){
#   values <- as.data.frame(plots_814_df)[x]
#   sd_all[count] <- sd(values[,1])
#   mean_all[count] <- mean(values[,1])
#   count <- count+1
# }
# 
# vars_814 <- as.data.frame(vars)
# vars_814$mean <- mean_all
# vars_814$sd <- sd_all
# 
# #write.csv(vars_814, "Q:/My Drive/Research/Canopy_Morphology/Tobacco_Project/canopy_analysis/tobacco_results/vars_814.csv")
# 
# levels <- topmod_814$lev[-1]
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
#   coefs[count] <- list(coef(topmod_814)[count,])
#   SE[count] <- list(summary(topmod_814)$standard.errors[count,])
#   z <- (summary(topmod_814)$coefficients/summary(topmod_814)$standard.errors)
#   Wald[count] <- list(z[count,])
#   p[count] <- list(((1 - pnorm(abs(z), 0, 1)) * 2)[count,])
#   std_coefs[count] <- list((sqrt(3)/pi) * (coefs[[count]][-1] * sd_all))
#   exp[count] <- list(exp(std_coefs[[count]]))
#   count <- count+1
# }
# 
# vars <- topmod_814$coefnames
# 
# results_0_814 <- as.data.frame(vars)
# results_0_814$level <- "0% Control"
# results_0_814$coefs <- coefs[[1]]
# results_0_814$std_coefs <- c(0,std_coefs[[1]])
# results_0_814$exp <- c(0,exp[[1]])
# results_0_814$SE <- SE[[1]]
# results_0_814$p <- p[[1]]
# results_0_814$wald <- Wald[[1]]
# 
# results_tox_814 <- as.data.frame(vars)
# results_tox_814$level <- "toxicity"
# results_tox_814$coefs <- coefs[[3]]
# results_tox_814$std_coefs <- c(0,std_coefs[[3]])
# results_tox_814$exp <- c(0,exp[[3]])
# results_tox_814$SE <- SE[[3]]
# results_tox_814$p <- p[[3]]
# results_tox_814$wald <- Wald[[3]]
# 
# results_def_814 <- as.data.frame(vars)
# results_def_814$level <- "deficiency"
# results_def_814$coefs <- coefs[[2]]
# results_def_814$std_coefs <- c(0,std_coefs[[2]])
# results_def_814$exp <- c(0,exp[[2]])
# results_def_814$SE <- SE[[2]]
# results_def_814$p <- p[[2]]
# results_def_814$wald <- Wald[[2]]
# 
# results_814 <- rbind(results_0_814, results_def_814, results_tox_814)
# 
# #write.csv(results_814, "Q:/My Drive/Research/Canopy_Morphology/Tobacco_Project/canopy_analysis/tobacco_results/results_814.csv")
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

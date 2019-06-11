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


#Load raster data for all dates
#setwd("/media/Kellyn/F20E17B40E177139/kpmontgo@ncsu.edu/LkWheeler_Sorghum/LkWheeler_Fusarium_Sorghum")
setwd("Q:/My Drive/Research/Canopy_Morphology/Tobacco_Project/canopy_analysis/layers")

tobacco_area <- readOGR("Boundary", "Boundary", stringsAsFactors = F)
tobacco_area <- spTransform(tobacco_area, CRS("+proj=lcc +lat_1=36.16666666666666 +lat_2=34.33333333333334 +lat_0=33.75 +lon_0=-79 +x_0=609601.22 +y_0=0 +datum=NAD83 +units=m +no_defs +ellps=GRS80 +towgs84=0,0,0"))

csm_618 <- crop(raster('CSM/CSM_618.tif'), tobacco_area)
csm_716 <- crop(raster('CSM/CSM_716.tif'), tobacco_area)
csm_814 <- crop(raster('CSM/CSM_814.tif'), tobacco_area)

csm_618_velox <- velox(csm_618)
csm_716_velox <- velox(csm_716)
csm_814_velox <- velox(csm_814)

#ortho_9_1 <- crop(raster("orthos/ortho_9_1_18.tif"))

#------------------------------------------------------- Plots -------------------------------------------------------------------

plots <- readOGR("plots", "plots_reshape", stringsAsFactors = F)
plots$Treatment <- as.factor(plots$Treatment)
plots <- spTransform(plots, CRS("+proj=lcc +lat_1=36.16666666666666 +lat_2=34.33333333333334 +lat_0=33.75 +lon_0=-79 +x_0=609601.22 +y_0=0 +datum=NAD83 +units=m +no_defs +ellps=GRS80 +towgs84=0,0,0"))

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

crr_metrics <- function(x){
  veloxRast <- velox(x)
  extract <- veloxRast$extract(sp=plots)
  names(extract) <- plots$Id
  crrmean <- sapply(extract, mean, na.rm=T)
  crrmedian <- sapply(extract, median, na.rm=T)
  crrsd <- sapply(extract, sd, na.rm=T)
  crrDF <- data.frame(plots=plots$Id, crrmedian= crrmedian, crrmean=crrmean, crrsd = crrsd)
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

extract <- mask(csm_618, plots[position,])
extract_df <- rasterToPoints(extract)
rumple <- rumple_index(extract_df[,1], extract_df[,2], extract_df[,3])

tmp <- plots@data
position <- match("606-E", tmp$Id)

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
  geary_all <- numeric(length(plotNums))
  count <- 1
  
  for (num in plotNums){
    tmp <- plots@data
    position <- match(num, tmp$Id)
    extract <- crop(x, plots[position,])
    moran <- Moran(extract, w)
    morans_all[count] <- moran
    geary <- Geary(extract, w)
    geary_all[count] <- geary
    count <- count+1
  }
  
  names(morans_all) <- plots$Id
  names(geary_all) <- plots$Id
  autocor_df <- data.frame(plots=plots$Id, moran= morans_all, geary=geary_all)
  autocor_df
}

#Using approx 30 cm x 30 cm weights matrix
autocor_618 <- autocor_metrics(csm_618, matrix(c(1,1,1,1,1,1,1,1,1,1,1,1,0,1,1,1,1,1,1,1,1,1,1,1,1), nc=5, nr=5))
autocor_716 <- autocor_metrics(csm_716, matrix(c(1,1,1,1,1,1,1,1,1,1,1,1,0,1,1,1,1,1,1,1,1,1,1,1,1), nc=5, nr=5))
autocor_814 <- autocor_metrics(csm_814, matrix(c(1,1,1,1,1,1,1,1,1,1,1,1,0,1,1,1,1,1,1,1,1,1,1,1,1), nc=5, nr=5))

plots_618 <- merge(plots_618, autocor_618, by.x="Id", by.y="plots")
plots_716 <- merge(plots_716, autocor_716, by.x="Id", by.y="plots")
plots_814 <- merge(plots_814, autocor_814, by.x="Id", by.y="plots")

# Aggregate treatments

library(data.table)

deficiency <- c("-K", "-N", "-P")
toxicity <- c("+B 2 lb", "+B 4 lb", "+B 8 lb")

dt_618 <- as.data.table(plots_618)
dt_618$Trt_agg <- as.character(plots_618$Treatment)
dt_618[Treatment %in% deficiency]$Trt_agg <- "deficiency"
dt_618[Treatment %in% toxicity]$Trt_agg <- "toxicity"
dt_618$Trt_agg <- as.factor(dt_618$Trt_agg)
plots_618@data <- dt_618

dt_716 <- as.data.table(plots_716)
dt_716$Trt_agg <- as.character(plots_716$Treatment)
dt_716[Treatment %in% deficiency]$Trt_agg <- "deficiency"
dt_716[Treatment %in% toxicity]$Trt_agg <- "toxicity"
dt_716$Trt_agg <- as.factor(dt_716$Trt_agg)
plots_716@data <- dt_716

dt_814 <- as.data.table(plots_814)
dt_814$Trt_agg <- as.character(plots_814$Treatment)
dt_814[Treatment %in% deficiency]$Trt_agg <- "deficiency"
dt_814[Treatment %in% toxicity]$Trt_agg <- "toxicity"
dt_814$Trt_agg <- as.factor(dt_814$Trt_agg)
plots_814@data <- dt_814

# --------------------------------------- Regression --------------------------------------------------------------------------------

library(PerformanceAnalytics)

plots_618_df <- as.data.frame(plots_618)
plots_716_df <- as.data.frame(plots_716)
plots_814_df <- as.data.frame(plots_814)

chart.Correlation(plots_618_df[3:18], 
                  method="spearman",
                  histogram=TRUE,
                  pch=16)

# normalize data
norm_618 <- as.data.frame(scale(plots_618_df[,3:18]))
norm_618$Treatment <- plots_618_df$Treatment
norm_618$Id <- plots_618_df$Id
norm_618$Trt_agg <- plots_618_df$Trt_agg

norm_716 <- as.data.frame(scale(plots_716_df[,3:18]))
norm_716$Treatment <- plots_716_df$Treatment
norm_716$Id <- plots_716_df$Id
norm_716$Trt_agg <- plots_716_df$Trt_agg

norm_814 <- as.data.frame(scale(plots_814_df[,3:18]))
norm_814$Treatment <- plots_814_df$Treatment
norm_814$Id <- plots_814_df$Id
norm_814$Trt_agg <- plots_814_df$Trt_agg

# all variables - crrmedian+crrmean+crrsd+median+mean+sd+iqr+var+sum+skew+kurt+PlotCRR+chRange+rumple+moran+geary
# 618 variables (post collinearity analysis) - 
#     all trt classes: crrmedian+crrsd+mean+sd+kurt+geary
#     agg trt classes: crrmedian+crrsd+skew+kurt+rumple+moran
# 716 variables (post collinearity analysis) - 
#     all trt classes: crrmean+crrsd+median+kurt+rumple+geary
#     agg trt classes: 
# 814 variables (post collinearity analysis) - 
#     all trt classes: crrmean+crrsd+median+iqr+kurt+geary
#     agg trt classes:

library(car)
glm_618 <- glm(Treatment ~ crrmedian+crrsd+median+kurt+rumple+geary, data = norm_618, 
              family = binomial(link="logit"), na.action = na.fail)
vif(glm_618)

glm_716 <- glm(Treatment ~ crrmean+crrsd+median+kurt+rumple+moran, data = norm_716, 
               family = binomial(link="logit"), na.action = na.fail)
vif(glm_716)

glm_814 <- glm(Treatment ~ crrmedian+median+iqr+kurt+PlotCRR+geary, data = norm_814, 
               family = binomial(link="logit"), na.action = na.fail)
vif(glm_814)

# Weights matrix
w <- knn2nb(knearneigh(coordinates(plots), k=8))

library(MuMIn)

train_data_indices <- rep(FALSE, nrow(norm_618))
train_data_indices[sample(1:nrow(norm_618), round(0.8 * nrow(norm_618)))] <- TRUE # randomly select 80% of the data for training

train_618 <- norm_618[train_data_indices, ]

glm_618 = glm(Treatment ~ crrmedian+crrsd+mean+sd+kurt+geary,
              data=train_618,
              family = binomial(link="logit"), na.action = na.fail
)

modsel_618 <- dredge(glm_618, beta = "none", trace = T)

glm_best_618 = glm(Treatment ~ crrsd,
                   data=norm_618,
                   family = binomial(link="logit"), na.action = na.fail
)

glm_best_618_agg = glm(Trt_agg ~ crrmedian+sd+sum+skew+kurt+moran,
                   data=norm_618,
                   family = binomial(link="logit"), na.action = na.fail
)

summary.glm(glm_best_618)
summary.glm(glm_best_618_agg)
# norm_618$resi <- glm_best_618$residuals
# moran.test(norm_618$resi, nb2listw(w))
# 
# durbinWatsonTest(glm_best_618)
# 
# library(ROCR)
# pred_data_618 <- norm_618[!train_data_indices, ]
# pred_618 = predict(glm_best_618,newdata=pred_data_618,type="response")
# fitpred = prediction(fitpreds,data$val$income)
# fitperf = performance(fitpred,"tpr","fpr")
# plot(fitperf,col="green",lwd=2,main="ROC Curve for Logistic:  Adult")
# abline(a=0,b=1,lwd=2,lty=2,col="gray")

modsel_716 <- dredge(glm_716, beta = "none", trace = T)

glm_best_716 = glm(Treatment ~ crrsd+median,
                   data=norm_716,
                   family = binomial(link="logit"), na.action = na.fail
)

summary.glm(glm_best_716)

glm_814 = glm(Treatment ~ crrmedian+median+iqr+kurt+PlotCRR+geary,
              data=norm_814,
              family = binomial(link="logit"), na.action = na.fail
)

modsel_814 <- dredge(glm_814, beta = "none", trace = T)
modsel_814

glm_best_814 = glm(Treatment ~ geary+iqr+kurt,
                   data=norm_814,
                   family = binomial(link="logit"), na.action = na.fail
)

summary.glm(glm_best_814)

# ------------------------------------- Random Forest Classification ------------------------------------------------------------
library(randomForest)
set.seed(42)
train_data_indices <- rep(FALSE, nrow(plots_618_df))
train_data_indices[sample(1:nrow(plots_618_df), round(0.8 * nrow(plots_618_df)))] <- TRUE # randomly select 80% of the data for training
rf_regression_618<- randomForest(Trt_agg ~ crrmedian+crrsd+skew+kurt+rumple+moran, data=plots_618_df[train_data_indices, ], importance=T)
rf_regression_618
varImpPlot(rf_regression_618)
pred_trt <- predict(rf_regression_618, plots_618_df[!train_data_indices,]) # predict the rings
plot(plots_618_df$Treatment[!train_data_indices], pred_trt, xlab="Observed", ylab="Predicted")
abline(a=0, b=1, lty=2, col=2)

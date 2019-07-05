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

#Load raster data for all dates
#setwd("/media/Kellyn/F20E17B40E177139/kpmontgo@ncsu.edu/LkWheeler_Sorghum/LkWheeler_Fusarium_Sorghum")
setwd("Q:/My Drive/Research/Canopy_Morphology/Tobacco_Project/canopy_analysis/layers")

tobacco_area <- readOGR("Boundary", "Boundary", stringsAsFactors = F)
tobacco_area <- spTransform(tobacco_area, CRS("+proj=lcc +lat_1=36.16666666666666 +lat_2=34.33333333333334 +lat_0=33.75 +lon_0=-79 +x_0=609601.22 +y_0=0 +datum=NAD83 +units=m +no_defs +ellps=GRS80 +towgs84=0,0,0"))

csm_618 <- crop(raster('CSM/CSM_618.tif'), tobacco_area)*100 #convert to cm
csm_716 <- crop(raster('CSM/CSM_716.tif'), tobacco_area)*100
csm_814 <- crop(raster('CSM/CSM_814.tif'), tobacco_area)*100

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

# Check if crr within plot is normally distributed
crr_618_df <- as.data.frame(crrRast_618)
shapiro.test(sample(crr_618_df$layer, 5000))
histogram(crr_618_df$layer)
velox_crr_618 <- velox(crrRast_618)
crr_extract_618 <- velox_crr_618$extract(sp=plots)
names(crr_extract_618) <- plots$Id

histogram(crr_extract_618$`502-E`)
shapiro.test(sample(crr_extract_618$`502-E`, 5000))
ggqqplot(sample(crr_extract_618$`502-E`, 5000))

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
dip.test(plot_extract_618$`502-E`)
histogram(plot_extract_618$`502-E`)

multimodality_test <- function(x){
  plotNums <- plots$Id
  dip_all <- numeric(length(plotNums))
  names(dip_all) <- "dip"
  count <- 1
  for (plot in x){
    test <- dip.test(plot)
    p <- test$p.value
    dip_all[count] <- p
    count <- count+1
  }
  dip_df <- data.frame(plots=plots$Id, dip= dip_all)
  dip_df
}

dip_df_list <- lapply(extract_list, multimodality_test)

#merge all metrics into spatial polygons
plots_618 <- merge(plots_618, dip_df_list[[1]], by.x="Id", by.y="plots")
plots_716 <- merge(plots_716, dip_df_list[[2]], by.x="Id", by.y="plots")
plots_814 <- merge(plots_814, dip_df_list[[3]], by.x="Id", by.y="plots")


# Check for normality
shapiro.test(plot_extract_618$`502-E`)

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

#Using approx 30 cm x 30 cm neighborhood matrix
autocor_618 <- autocor_metrics(csm_618, matrix(c(1,1,1,1,1,1,1,1,1,1,1,1,0,1,1,1,1,1,1,1,1,1,1,1,1), nc=5, nr=5))
autocor_716 <- autocor_metrics(csm_716, matrix(c(1,1,1,1,1,1,1,1,1,1,1,1,0,1,1,1,1,1,1,1,1,1,1,1,1), nc=5, nr=5))
autocor_814 <- autocor_metrics(csm_814, matrix(c(1,1,1,1,1,1,1,1,1,1,1,1,0,1,1,1,1,1,1,1,1,1,1,1,1), nc=5, nr=5))

plots_618 <- merge(plots_618, autocor_618, by.x="Id", by.y="plots")
plots_716 <- merge(plots_716, autocor_716, by.x="Id", by.y="plots")
plots_814 <- merge(plots_814, autocor_814, by.x="Id", by.y="plots")

#-------------------------------------Normalize variables---------------------------------------

library(rcompanion)

shapiro.test(plots_618$crrmedian)
shapiro.test(plots_618$crriqr)
shapiro.test(plots_618$dip)
shapiro.test(plots_618$median)
shapiro.test(plots_618$iqr) #normal
shapiro.test(plots_618$skew)
shapiro.test(plots_618$kurt)
shapiro.test(plots_618$PlotCRR) #normal
shapiro.test(plots_618$rumple)
shapiro.test(plots_618$moran)
shapiro.test(plots_618$geary)

# NOTE the transformations below are not all natural log. Tukey transformation choses the best transformation method.
plots_618$crrmedian_ln <- transformTukey(plots_618$crrmedian)
plots_618$crriqr_ln <- transformTukey(plots_618$crriqr) # still not normal
plots_618$dip_ln <- transformTukey(plots_618$dip) # still not normal
plots_618$median_ln <- transformTukey(plots_618$median)
plots_618$skew_ln <- transformTukey(plots_618$skew) # still not normal
plots_618$kurt_ln <- transformTukey(plots_618$kurt) # still not normal
plots_618$moran_ln <- transformTukey(plots_618$moran)
plots_618$geary_ln <- transformTukey(plots_618$geary)

plots_716$crrmedian_ln <- transformTukey(plots_716$crrmedian)
plots_716$crriqr_ln <- transformTukey(plots_716$crriqr)
plots_716$dip_ln <- transformTukey(plots_716$dip) # still not normal
plots_716$median_ln <- transformTukey(plots_716$median)
plots_716$skew_ln <- transformTukey(plots_716$skew) # still not normal
plots_716$kurt_ln <- transformTukey(plots_716$kurt) # still not normal
plots_716$moran_ln <- transformTukey(plots_716$moran)
plots_716$geary_ln <- transformTukey(plots_716$geary)

plots_814$crrmedian_ln <- transformTukey(plots_814$crrmedian)
plots_814$crriqr_ln <- transformTukey(plots_814$crriqr)
plots_814$dip_ln <- transformTukey(plots_814$dip) # still not normal
plots_814$median_ln <- transformTukey(plots_814$median)
plots_814$skew_ln <- transformTukey(plots_814$skew) 
plots_814$kurt_ln <- transformTukey(plots_814$kurt)
plots_814$moran_ln <- transformTukey(plots_814$moran) # still not normal
plots_814$geary_ln <- transformTukey(plots_814$geary)

# --------------------------------------- Regression --------------------------------------------------------------------------------
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


library(PerformanceAnalytics)

plots_618_df <- as.data.frame(plots_618)
plots_716_df <- as.data.frame(plots_716)
plots_814_df <- as.data.frame(plots_814)

chart.Correlation(plots_618_df[,c(3, 7, 14, 17, 18)], 
                  method="pearson",
                  histogram=TRUE,
                  pch=16)

# normalize data
norm_618 <- as.data.frame(scale(plots_618@data[,c(3, 5,6, 7, 10, 13, 14, 15, 17, 18)]))
norm_618$Treatment <- plots_618@data$Treatment
norm_618$Id <- plots_618@data$Id
norm_618$Trt_agg <- plots_618@data$Trt_agg

norm_716 <- as.data.frame(scale(plots_716@data[,3:27]))
norm_716$Treatment <- plots_716$Treatment
norm_716$Id <- plots_716$Id
norm_716$Trt_agg <- plots_716$Trt_agg

norm_814 <- as.data.frame(scale(plots_814@data[,3:27]))
norm_814$Treatment <- plots_814$Treatment
norm_814$Id <- plots_814$Id
norm_814$Trt_agg <- plots_814_df$Trt_agg

# all variables - crrmedian_ln+crriqr_ln+median_ln+iqr+skew_ln+kurt_ln+PlotCRR+rumple+moran_ln+geary_ln+dip_ln
# 618 variables (post collinearity analysis) - 
#     all trt classes: crrmedian_ln+median_ln+iqr+PlotCRR+rumple+moran_ln+dip_ln
#     agg trt classes: crrmedian_ln+crriqr_ln+median_ln+skew_ln+rumple+moran_ln+dip_ln
# 716 variables (post collinearity analysis) - 
#     all trt classes: crrmedian_ln+crriqr_ln+median_ln+skew_ln+rumple+geary_ln+dip_ln
#     agg trt classes: crrmedian_ln+crriqr_ln+median_ln+skew_ln+rumple+moran_ln+dip_ln
# 814 variables (post collinearity analysis) - 
#     all trt classes: crrmedian_ln+crriqr_ln+median_ln+iqr+kurt_ln+geary_ln+dip_ln
#     agg trt classes: crrmedian_ln+crriqr_ln+median_ln+iqr+skew_ln+moran_ln+geary_ln+dip_ln

library(car)
glm_618 <- glm(Trt_agg ~ crrmedian+median+kurt+rumple+moran+dip, data = plots_618, 
              family = binomial(link="logit"), na.action = na.fail)
vif(glm_618)

glm_716 <- glm(Trt_agg ~ crrmedian_ln+crriqr_ln+median_ln+skew_ln+rumple+moran_ln+dip_ln, data = norm_716, 
               family = binomial(link="logit"), na.action = na.fail)
vif(glm_716)

glm_814 <- glm(Trt_agg ~ crrmedian_ln+crriqr_ln+median_ln+iqr+skew_ln+moran_ln+geary_ln+dip_ln, data = norm_814, 
               family = binomial(link="logit"), na.action = na.fail)
vif(glm_814)

# Try multinominal logisitic regression with nnet package

library(nnet)
plots_618_df$Trt_agg <- relevel(plots_618_df$Trt_agg, ref = "100% Control")
mltnom_618 <- multinom(Trt_agg ~ crrmedian+median+kurt+rumple+moran+crriqr+skew+iqr, data = plots_618_df, na.action = na.fail)
summary(mltnom_618)

modsel_618 <- dredge(mltnom_618)
modsel_618

topmod_618 <- multinom(Trt_agg ~ median+kurt+rumple+moran+skew+iqr, data = plots_618_df, na.action = na.fail)
summary(topmod_618)
z <- summary(topmod_618)$coefficients/summary(topmod_618)$standard.errors
z
p <- (1 - pnorm(abs(z), 0, 1)) * 2
p

coefs_618 <- coef(topmod_618)
exp(coefs_618) #transformed to odds ratio
(exp(coefs_618)-1)*100 # percent change in the odds for a one unit increase in the independent variable

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



modsel_716 <- dredge(glm_716, beta = "none", trace = T)

glm_best_716 = glm(Trt_agg ~ crrmedian_ln+median_ln+dip_ln+moran_ln+skew_ln,
                   data=norm_716,
                   family = binomial(link="logit"), na.action = na.fail
)

summary.glm(glm_best_716)

plots_716$glm_trtag_resi <-glm_best_618$residuals
spplot(plots_716, zcol = "glm_trtag_resi")

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
train_data_indices[sample(1:nrow(plots_618_df), round(0.7 * nrow(plots_618_df)))] <- TRUE # randomly select 80% of the data for training
rf_regression_618<- randomForest(Trt_agg ~ crrmedian_ln+median_ln+kurt_ln+rumple+geary_ln+dip_ln, data=plots_618_df[train_data_indices, ], importance=T)
rf_regression_618
varImpPlot(rf_regression_618)
pred_trt <- predict(rf_regression_618, plots_618_df[!train_data_indices,]) # predict the rings
table(pred_trt, plots_618_df$Trt_agg[!train_data_indices])
plot(plots_618_df$Trt_agg[!train_data_indices], pred_trt, xlab="Observed", ylab="Predicted")
abline(a=0, b=1, lty=2, col=2)

set.seed(42)
train_data_indices <- rep(FALSE, nrow(plots_716_df))
train_data_indices[sample(1:nrow(plots_716_df), round(0.7 * nrow(plots_716_df)))] <- TRUE # randomly select 80% of the data for training
rf_regression_716<- randomForest(Trt_agg ~ crrmedian_ln+median_ln+kurt_ln+rumple+moran_ln+dip_ln, data=plots_716_df[train_data_indices, ], importance=T)
rf_regression_716
varImpPlot(rf_regression_716)
pred_trt <- predict(rf_regression_716, plots_716_df[!train_data_indices,]) # predict the rings
table(pred_trt, plots_716_df$Trt_agg[!train_data_indices])
plot(plots_716_df$Trt_agg[!train_data_indices], pred_trt, xlab="Observed", ylab="Predicted")
abline(a=0, b=1, lty=2, col=2)

set.seed(42)
train_data_indices <- rep(FALSE, nrow(plots_814_df))
train_data_indices[sample(1:nrow(plots_814_df), round(0.7 * nrow(plots_814_df)))] <- TRUE # randomly select 80% of the data for training
rf_regression_814<- randomForest(Trt_agg ~ crrmedian_ln+median_ln+kurt_ln+rumple+moran_ln+dip_ln, data=plots_814_df[train_data_indices, ], importance=T)
rf_regression_814
varImpPlot(rf_regression_814)
pred_trt <- predict(rf_regression_814, plots_814_df[!train_data_indices,]) # predict the rings
table(pred_trt, plots_814_df$Trt_agg[!train_data_indices])
plot(plots_814_df$Trt_agg[!train_data_indices], pred_trt, xlab="Observed", ylab="Predicted")
abline(a=0, b=1, lty=2, col=2)

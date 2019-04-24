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
# roughly equal to one plant per polygon
#sorghum_polygons <- readOGR("sorghum_polygons", "sorghum_polygons", stringsAsFactors = F)

# areas below are null, min and max values across all polygons
# csm_vegOnly_9_1 <- crop(raster("CSM/vegOnly_9_1.tif"), sorghum_area)
# glob_min_9_1 <- 0.0000000264168
# glob_max_9_1 <- 1.40754
# glob_range_9_1 <- glob_max_9_1-glob_min_9_1

csm_adj_9_1 <- crop(raster('CSM/csm_adj_9_1_noDepth_20cm.tif'), sorghum_area)
#csm_adj_9_7 <- crop(raster("CSM/csm_adj_9_7_18_5sm.tif"), sorghum_area)
# csm_adj_9_21 <- crop(raster("CSM/csm_adj_9_1_18_0sm_noDepth.tif"), sorghum_area)
# #csm_adj_10_2 <- crop(raster("CSM/csm_adj_10_2_18_5sm.tif"), sorghum_area)
# csm_adj_10_12 <- crop(raster("CSM/csm_adj_9_1_18_0sm_noDepth.tif"), sorghum_area)
# 
# weight_slope <- crop(raster("series_weighted_slope.tif"), sorghum_area)

#ortho_9_1 <- crop(raster("orthos/ortho_9_1_18.tif"))

# SIMWE water depth model
water <- crop(raster("simwe/water_depth.tif"), sorghum_area)

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

crr <- function(x){
  (mean(x, na.rm=T)-min(x, na.rm=T))/(max(x, na.rm=T)-min(x, na.rm=T))
}

# crr_globRange <- function(x){
#   (mean(x, na.rm=T)-min(x, na.rm=T))/(glob_max_9_1-glob_min_9_1)
# }
# 
# pred_r_squared <- function(linear.model) {
#   lm.anova <- anova(linear.model)
#   tss <- sum(lm.anova$"Sum Sq")
#   # predictive R^2
#   pred.r.squared <- 1 - PRESS(linear.model)/(tss)
#   return(pred.r.squared)
# }
# 
# PRESS <- function(linear.model) {
#   pr <- residuals(linear.model)/(1 - lm.influence(linear.model)$hat)
#   PRESS <- sum(pr^2)
#   return(PRESS)
# }

csm_adj_9_1_velox <- velox(csm_adj_9_1)
# csm_adj_9_21_velox <- velox(csm_adj_9_21)
# csm_adj_10_12_velox <- velox(csm_adj_10_12)
# 
# weight_slope_velox <- velox(weight_slope)

#veg_9_1_velox<- velox(csm_vegOnly_9_1)



# Canopy relief ratio

#crrRast_9_1 <- focal(csm_adj_9_1, w=matrix(1,3,3), crr)
crrRast_9_1 <- crop(raster('CRR/crr_15_neighbor.tif'), sorghum_area)
crrRast_9_1_velox <- velox(crrRast_9_1)
crr_extract_9_1 <- crrRast_9_1_velox$extract(sp=plots)
names(crr_extract_9_1) <- plots$id

# histogram(crr_extract_9_1$`410`, nint=28)
# histogram(crr_extract_9_1$`401`, nint=28)
# histogram(crr_extract_9_1$`407`, nint=28)
# histogram(crr_extract_9_1$`203`, nint=28)
# histogram(crr_extract_9_1$`202`, nint=28)
# histogram(crr_extract_9_1$`205`, nint=28)

crrmedian <- sapply(crr_extract_9_1, median, na.rm=TRUE)
crrmean <- sapply(crr_extract_9_1, mean, na.rm=TRUE)
crrsd <- sapply(crr_extract_9_1, sd, na.rm=TRUE)

crrDF <- data.frame(plots=plots$id, crrmedian= crrmedian, crrmean=crrmean, crrsd = crrsd)
# Normalize data
norm_crr <- as.data.frame(scale(crrDF[,2:4]))
norm_crr$plots <- crrDF$plots

#merge all metrics into spatial polygons
plots_9_1 <- plots
plots_9_1 <- merge(plots_9_1, norm_crr, by.x="id", by.y="plots")
plots_9_1 <- merge(plots_9_1, field_data, by.x="id", by.y="Plot")

# spplot(plots_9_1, zcol="BuAc")
# 
# # weighted slope
# 
# slope_extract <-weight_slope_velox$extract(sp=plots)
# names(slope_extract) <- plots$id
# 
# slopemedian <- sapply(slope_extract, median, na.rm=TRUE)
# slopemean <- sapply(slope_extract, mean, na.rm=TRUE)
# slopesd <- sapply(slope_extract, sd, na.rm=TRUE)
# slopesum <- sapply(slope_extract, sum, na.rm=TRUE)
# 
# slopeDF <- data.frame(plots=plots$id, slopesum = slopesum, slopesd=slopesd, slopemean=slopemean)
# plots_9_1 <- merge(plots_9_1, slopeDF, by.x="id", by.y="plots")

# SIMWE water depth

water_velox <- velox(water)
water_extract <- water_velox$extract(sp=plots)
names(water_extract) <- plots$id

# histogram(water_extract$`410`, nint=28)


watermedian <- sapply(water_extract, median, na.rm=TRUE)
watermean <- sapply(water_extract, mean, na.rm=TRUE)
watersd <- sapply(water_extract, sd, na.rm=TRUE)
watersum <- sapply(water_extract, sum, na.rm=TRUE)

waterDF <- data.frame(plots=plots$id, watersum = watersum, watersd=watersd, watermean=watermean)
# Normalize data
norm_water <- as.data.frame(scale(waterDF[,2:4]))
norm_water$plots <- waterDF$plots

plots_9_1 <- merge(plots_9_1, norm_water, by.x="id", by.y="plots")

# Crop height data exploration
plot_extract_9_1 <- csm_adj_9_1_velox$extract(sp=plots)
names(plot_extract_9_1) <- plots$id


#library("sm")

hist(plot_extract_9_1$`410`, breaks=45, xlab = "Height (m)", col = "green") #highest yield
histogram(plot_extract_9_1$`401`, nint=45) #2nd highest yield, descending order
histogram(plot_extract_9_1$`407`, nint=45)
histogram(plot_extract_9_1$`406`, nint=45)
histogram(plot_extract_9_1$`408`, nint=45)

histogram(plot_extract_9_1$`203`, nint=45, xlab = "Lowest yield plot") #lowest yield
histogram(plot_extract_9_1$`202`, nint=45)
histogram(plot_extract_9_1$`205`, nint=45)
histogram(plot_extract_9_1$`204`, nint=45)
histogram(plot_extract_9_1$`212`, nint=45)

# Cumulative distribution functions
# extract_unlist_9_1 <- lapply(plot_extract_9_1, unlist)
# tmp <- lapply(plot_extract_9_1, ecdf)
# ch_ecdf_9_1 <- lapply(plot_extract_9_1, ecdf)
# par(mar=c(5, 6, 4, 2) + 0.1)
# plot(ch_ecdf_9_1[["410"]], add=T)#main="Hypsographs - 90th Percentile Yield Plots", xlab="Crop Height (m)", ylab="Fraction of Data", cex.axis=2, lwd=5, cex.lab=2, col="blue") #highest yield
# plot(ch_ecdf_9_1[["401"]], add=T)
# plot(ch_ecdf_9_1[["407"]], add=T)
# plot(ch_ecdf_9_1[["406"]], add=T)
# plot(ch_ecdf_9_1[["408"]], add=T)
# plot(ch_ecdf_9_1[["405"]])
# plot(ch_ecdf_9_1[["402"]])
# plot(ch_ecdf_9_1[["309"]])
# plot(ch_ecdf_9_1[["404"]])
# plot(ch_ecdf_9_1[["308"]])
# 
# plot(ch_ecdf_9_1[["203"]], add=T) #main="Hypsographs - 10th Percentile Yield Plots", xlab="crop height (m)", ylab="Fraction of Data", cex.axis=2, lwd=5, cex.lab=2, col="blue")#lowest yield
# plot(ch_ecdf_9_1[["202"]], add=T)
# plot(ch_ecdf_9_1[["205"]], add=T)
# plot(ch_ecdf_9_1[["204"]], add=T)
# plot(ch_ecdf_9_1[["212"]], add=T)
# plot(ch_ecdf_9_1[["102"]])
# plot(ch_ecdf_9_1[["305"]])
# plot(ch_ecdf_9_1[["206"]])
# plot(ch_ecdf_9_1[["101"]])
# plot(ch_ecdf_9_1[["306"]])
#
# Top 10% yield avg hypsograph
top_plots <- c('410', '401','407', '406', '408')
top_yield_9_1 <- plot_extract_9_1[top_plots]
top_yield_9_1_all <- do.call(c, top_yield_9_1)
# top_ecdf_9_1 <- ecdf(top_yield_9_1_all)
# par(mar=c(5, 6, 4, 2) + 0.1)
# plot(top_ecdf_9_1, main="Hypsographs - 90th Percentile Yield", xlab="Crop Height (m)", ylab="Fraction of Data", cex.axis=2, lwd=5, cex.lab=2, col="blue")

# # Bottom 10% yield avg hypsograph
bottom_plots <- c('203', '202','205', '204', '212')
bottom_yield_9_1 <- plot_extract_9_1[bottom_plots]
bottom_yield_9_1_all <- do.call(c, bottom_yield_9_1)
# bottom_ecdf_9_1 <- ecdf(bottom_yield_9_1_all)
# par(mar=c(5, 6, 4, 2) + 0.1, mgp=c(3, 1, 0))
# plot(bottom_ecdf_9_1, main="Hypsographs - 10th Percentile Yield", xlab="Crop Height (m)", ylab="Fraction of Data", cex.axis=2, lwd=5, cex.lab=2, col="blue")

# # #Hypsographs plotted together
# par(mar=c(6, 6.5, 4, 2) + 0.1, mgp=c(4, 1, 0))
# plot(bottom_ecdf_9_1, main="Hypsograph", xlab="Crop Height (m)", ylab="Fraction of Data", cex.main=4, cex.axis=2.5, lwd=5, cex.lab=3, col="red")
# plot(top_ecdf_9_1,lwd=5, cex.lab=2, col="blue", add=T)
# legend(0,0.9, legend=c("10th percentile yield", "90th percentile yield"), col=(c("red", "blue")), lty=1, lwd=5, cex=2.5, bty="n", y.intersp = 0.6)

# Histogram and denisty of top vs bottom yield
library(ggplot2)

bottom_yield_9_1_all <- as.data.frame(bottom_yield_9_1_all)
names(bottom_yield_9_1_all) <- "height"
bottom_yield_9_1_all$group <- "bottom"
top_yield_9_1_all <- as.data.frame(top_yield_9_1_all)
names(top_yield_9_1_all) <- "height"
top_yield_9_1_all$group <- "top"
pct_10_90_yield <- rbind(top_yield_9_1_all, bottom_yield_9_1_all)

height_mean <- ddply(pct_10_90_yield, "group", summarise, height.mean=mean(height))

ggplot(pct_10_90_yield, aes(x=height, fill=group)) +
  geom_histogram(binwidth=.03, alpha=.5, position="identity") +
  geom_vline(data=height_mean, aes(xintercept=height.mean,  colour=group),
             linetype="dashed", size=1)

# Get sill and range values from fitted variogram models

plotNums <- plots$id
range_all <- numeric(length(plotNums))
sill_all <- numeric(length(plotNums))
count <- 1

for (num in plotNums){
  tmp <- plots@data
  position <- match(num, tmp$id)
  extract <- crop(csm_adj_9_1, plots[position,])
  samp <- sampleRandom(extract, 200, sp=T, round(2))
  colnames(samp@data)[1] <- "height"
  variogram <- variogram(samp$height~1, locations=samp)
  fit_vgm <- fit.variogram(variogram, vgm(c("Exp", "Sph", "Mat", "Gau")))
  range <- round(fit_vgm$range[2], 4)
  sill <- round(fit_vgm$psill[2], 4)
  range_all[count] <- range
  sill_all[count] <- sill
  count <- count+1
}

names(sill_all) <- plots$id
names(range_all) <- plots$id
var_df <- data.frame(plots=plots$id, sill= sill_all, range=range_all)
# 
# Normalize data
norm_var <- as.data.frame(scale(var_df[,2:3]))
norm_var$plots <- var_df$plots
plots_9_1 <- merge(plots_9_1, norm_var, by.x="id", by.y="plots")

# breaks.qt <- classIntervals(plots_9_1$range, n = 16, style = "quantile", intervalClosure = "right")
# spplot(plots_9_1, zcol="range", main=list(label="Variogram Sill", cex=1.5), par.settings=list(fontsize=list(text=25)), col.regions = viridis(16), at=breaks.qt$brks)


# Variogram with sample points
#samp_410 <- sampleRandom(plot_410, 4000, sp=T, round(2))
#colnames(samp_410@data)[1] <- "height"

# variogram_410 <- variogram(samp_410$height~1, locations=samp_410)
# fit_vgm_410 <- fit.variogram(variogram_410, vgm(c("Exp", "Sph", "Mat", "Gau")))
# range <- round(fit_vgm_410$range[2], 4)
# print(range)
# sill <- round(fit_vgm_410$psill[2], 4)
# print(sill)
# modelAnnotation <- paste("Range =", range, "Sill =", sill, sep=" ")
# 
# mypanel = function(x,y,...) { 
#   vgm.panel.xyplot(x,y,...)
#   panel.abline(h=fit_vgm_410$psill[2], color = 'red')
#   panel.abline(v=fit_vgm_410$range[2], color = 'red')
#   panel.text(4,0.01,label=modelAnnotation, font=1, cex=2)
# }
# plot(variogram_410, fit_vgm_410, panel=mypanel, main=list(label="Plot 410 Variogram and Fitted Model"), xlab=list(label="Distance (m)",cex=2), ylab=list(label="Semivariance",cex=2), scale=list(cex=2), lwd=5, col="red", pch=21, cex=2)


# Crop height metrics
set.seed(42)
CHmedian <- sapply(plot_extract_9_1, median, na.rm=TRUE)
CHmean <- sapply(plot_extract_9_1, mean, na.rm=TRUE)
CHsd <- sapply(plot_extract_9_1, sd, na.rm=TRUE)
CHvar <- sapply(plot_extract_9_1, var, na.rm=TRUE)
CHsum <- sapply(plot_extract_9_1, sum, na.rm=TRUE)
CHskew <- sapply(plot_extract_9_1, skewness, na.rm=TRUE)
CHkurt <- sapply(plot_extract_9_1, kurtosis, na.rm=TRUE)
CHcrr <- sapply(plot_extract_9_1, crr)
CHrange <- sapply(plot_extract_9_1, range, na.rm=TRUE)
CHrange <- CHrange[2,]-CHrange[1,]

cropHeightDF <- data.frame(plots=plots$id, median= CHmedian, mean=CHmean, sd = CHsd, var = CHvar, sum = CHsum, skew = CHskew, kurt= CHkurt, PlotCRR = CHcrr, chRange = CHrange)

# Normalize data
norm_ch <- as.data.frame(scale(cropHeightDF[,2:10]))
norm_ch$plots <- cropHeightDF$plots
# norm_classified <- as.data.frame(scale(classifiedCounts_9_1[,1:4]))
# norm_classified$plots <- cropHeightDF$plots
# norm_geomorph <- as.data.frame(scale(geomorphCounts_9_1[,1:10]))
# norm_geomorph$plots <- cropHeightDF$plots
# norm_varInd <- as.data.frame(scale(varIndDF[,2:7]))
# norm_varInd$plots <- cropHeightDF$plots
# norm_tgi <- as.data.frame(scale(tgiDF[,2:7]))
# norm_tgi$plots <- cropHeightDF$plots

#merge all metrics into spatial polygons
plots_9_1 <- merge(plots_9_1, norm_ch, by.x="id", by.y="plots")

# plots_9_1 <- merge(plots_9_1, norm_classified, by.x="id", by.y="plots")
# plots_9_1 <- merge(plots_9_1, norm_geomorph, by.x="id", by.y="plots")
# plots_9_1 <- merge(plots_9_1, norm_varInd, by.x="id", by.y="plots")
# plots_9_1 <- merge(plots_9_1, norm_tgi, by.x="id", by.y="plots")

#my.palette <- brewer.pal(n = 8, name = "YlGnBu")

spplot(plots_9_1, zcol="crrsd", main=list(label="Canopy Relief Ratio", cex=2))
spplot(plots_9_1, zcol="mean", main=list(label="Mean Crop Height (m)", cex=2))
spplot(plots_9_1, zcol="sd", main=list(label="Canopy Rugosity", cex=2))
breaks.qt <- classIntervals(plots_9_1$range, n = 16, style = "quantile", intervalClosure = "right")
spplot(plots_9_1, zcol="range", main=list(label="Yield (bu/ac)", cex=1.5), par.settings=list(fontsize=list(text=25)), col.regions = viridis(16), at=breaks.qt$brks)
#grid.arrange(p1, p2, p3, nrow = 1)

plot(plots_9_1$BuAc, plots_9_1$PlotCRR)


# Use Moran's I to compute spatial autocorrelation for each plot
tmp <- MoranLocal(csm_adj_9_1)
plot(tmp)

plotNums <- plots$id
morans_all <- numeric(length(plotNums))
geary_all <- numeric(length(plotNums))
count <- 1

for (num in plotNums){
  tmp <- plots@data
  position <- match(num, tmp$id)
  extract <- crop(csm_adj_9_1, plots[position,])
  #samp <- sampleRandom(extract, 4000, sp=T, round(2))
  #colnames(samp@data)[1] <- "height"
  moran <- Moran(extract)
  morans_all[count] <- moran
  geary <- Geary(extract)
  geary_all[count] <- geary
  count <- count+1
}

names(morans_all) <- plots$id
names(geary_all) <- plots$id
autocor_df <- data.frame(plots=plots$id, moran= morans_all, geary=geary_all)

plots_9_1 <- merge(plots_9_1, autocor_df, by.x="id", by.y="plots")

spplot(plots_9_1, zcol="geary")


plot(plots_9_1$mean, plots_9_1$BuAc)
plot(plots_9_1$sum, plots_9_1$BuAc)
plot(plots_9_1$skew, plots_9_1$BuAc)
plot(plots_9_1$kurt, plots_9_1$BuAc)
plot(plots_9_1$sd, plots_9_1$BuAc)
plot(plots_9_1$crrmean, plots_9_1$BuAc)
plot(plots_9_1$moran, plots_9_1$BuAc)
plot(plots_9_1$geary, plots_9_1$BuAc)
plot(plots_9_1$sill, plots_9_1$BuAc)
plot(plots_9_1$range, plots_9_1$BuAc)
plot(plots_9_1$Trt, plots_9_1$BuAc)
plot(plots_9_1$AUDPC, plots_9_1$BuAc)
plot(plots_9_1$watersum, plots_9_1$BuAc)
plot(plots_9_1$watersd, plots_9_1$BuAc)

# library(diptest)
# dip.test(plot_extract_9_1$'203')

hist(plot_extract_9_1$`410`, breaks=45) #highest yield
#hist(ozone, breaks = 15, freq = F, xlab = 'Ozone (ppb)', yl)im = c(0, 0.025), ylab = 'Probability', main = 'Histogram of Ozone Pollution Data with Kernel Density Plot')
lines(density(plot_extract_9_1$'410', na.rm = T))
histogram(plot_extract_9_1$`401`, nint=45) #2nd highest yield, descending order
histogram(plot_extract_9_1$`407`, nint=45)
histogram(plot_extract_9_1$`406`, nint=45)
histogram(plot_extract_9_1$`408`, nint=45)

histogram(plot_extract_9_1$`203`, nint=45) #lowest yield
plot(density(plot_extract_9_1$'203'))
histogram(plot_extract_9_1$`202`, nint=45)
histogram(plot_extract_9_1$`205`, nint=45)
histogram(plot_extract_9_1$`204`, nint=45)
histogram(plot_extract_9_1$`212`, nint=45)


write.csv(plots_9_1, "plots_9_1_metrics_20cm.csv")

# repeat metrics by polygon
# polygon_extract_9_1 <- veg_9_1_velox$extract(sp=sorghum_polygons)
# # CRR using global range
# crr_plant <- sapply(polygon_extract_9_1, crr_globRange)
# 
# CHmedian <- sapply(polygon_extract_9_1, median, na.rm=TRUE)
# CHmean <- sapply(polygon_extract_9_1, mean, na.rm=TRUE)
# CHsd <- sapply(polygon_extract_9_1, sd, na.rm=TRUE)
# CHvar <- sapply(polygon_extract_9_1, var, na.rm=TRUE)
# CHsum <- sapply(polygon_extract_9_1, sum, na.rm=TRUE)
# CHskew <- sapply(polygon_extract_9_1, skewness, na.rm=TRUE)
# CHkurt <- sapply(polygon_extract_9_1, kurtosis, na.rm=TRUE)
# 
# polygons_9_1 <- sorghum_polygons
# polygons_9_1$crr_plant <- crr_plant
# #spplot(polygons_9_1, zcol = "crr_plant")
# polygons_9_1$med <- CHmedian
# polygons_9_1$mean <- CHmean
# polygons_9_1$sd <- CHsd
# polygons_9_1$var <- CHvar
# polygons_9_1$sum <- CHsum
# polygons_9_1$skew <- CHskew
# polygons_9_1$kurt <- CHkurt
# 
# breaks.qt <- classIntervals(polygons_9_1$mean, n = 14, style = "quantile", intervalClosure = "right")
# spplot(polygons_9_1, zcol="mean", main=list(label="Mean Crop Height (m)", cex=1.5),par.settings=list(fontsize=list(text=25)), col.regions=viridis(16),at=breaks.qt$brks)
# breaks.qt <- classIntervals(polygons_9_1$sd, n = 14, style = "quantile", intervalClosure = "right")
# spplot(polygons_9_1, zcol="sd", main=list(label="Canopy Rugosity",  cex=1.5),  col.regions=viridis(16) , par.settings=list(fontsize=list(text=25)), at=breaks.qt$brks)
# breaks.qt <- classIntervals(polygons_9_1$crr_plant, n = 14, style = "quantile", intervalClosure = "right")
# spplot(polygons_9_1, zcol="crr_plant", main=list(label="Canopy Relief Ratio",  cex=1.5), par.settings=list(fontsize=list(text=25)), col.regions=viridis(16) ,at=breaks.qt$brks)
# 
# 
# spplot(polygons_9_1[which(polygons_9_1$id_1=='410'),], zcol="mean", main=list(label="Plot 410 Mean Crop Height"))
# spplot(polygons_9_1[which(polygons_9_1$id_1=='203'),], zcol="mean", main=list(label="Plot 203 Mean Crop Height"))
# #grid.arrange(p1, p2, nrow=1)
# 
# spplot(polygons_9_1[which(polygons_9_1$id_1=='410'),], zcol="crr_plant", main=list(label="Plot 410 Canopy Relief Ratio"))
# spplot(polygons_9_1[which(polygons_9_1$id_1=='203'),], zcol="crr_plant", main=list(label="Plot 203 Canopy Relief Ratio"))
# 
# p1 <- spplot(polygons_9_1, zcol="crr_plant")
# p2 <- spplot(polygons_9_1, zcol="mean")
# p3 <- spplot(polygons_9_1, zcol="sd")
# grid.arrange(p1, p2, p3, nrow = 1)
# 
# # Aggregate by plant metrics to plot level
# crs.new <- CRS("+proj=lcc +lat_1=36.16666666666666 +lat_2=34.33333333333334
#                +lat_0=33.75 +lon_0=-79 +x_0=609601.22 +y_0=0 +datum=NAD83 +units=m
#                +no_defs +ellps=GRS80 +towgs84=0,0,0")
# proj4string(plots_9_1) <- crs.new
# proj4string(polygons_9_1) <- crs.new
# polycrr_avg <- over(plots_9_1, polygons_9_1[,"crr_plant"], fn = mean)
# plots_9_1$crr_poly <- polycrr_avg$crr_plant
# polycrr_sd <- over(plots_9_1, polygons_9_1[,"crr_plant"], fn = sd)
# plots_9_1$crr_sd_poly <- polycrr_sd$crr_plant
# 
# polyMean_avg <- over(plots_9_1, polygons_9_1[,"mean"], fn = mean)
# plots_9_1$mean_poly <- polyAvg_mean$mean
# polyMean_sd <- over(plots_9_1, polygons_9_1[,"mean"], fn = sd)
# plots_9_1$mean_sd_poly <- polyMean_sd$mean
# 
# polySd_avg <- over(plots_9_1, polygons_9_1[,"sd"], fn = mean)
# plots_9_1$sd_poly <- polyAvg_sd$sd
# polySd_sd <- over(plots_9_1, polygons_9_1[,"sd"], fn = sd)
# plots_9_1$sd_sd_poly <- polySd_sd$sd
# 
# top_bottom <- c(top_plots, bottom_plots)
# plots_9_1[plots_9_1$id %in% top_bottom,c(1,25:30)]@data
# 
# 

# Statistical models

ols_yield_9_1 <- lm(BuAc ~ mean+sd+skew+PlotCRR, data=plots_9_1)
plots_9_1$ols_yi_res <- ols_yield_9_1$residuals
summary(ols_yield_9_1)
plot(ols_yield_9_1)
spplot(plots_9_1, zcol="ols_yi_res")
print(pred_r2_9_1 <- pred_r_squared(ols_yield_9_1))

#w <- 1/as.matrix(dist(coordinates(plots_9_1)))
#diag(w) <- 0
w <- knn2nb(knearneigh(coordinates(plots), k=8))
moran.test(plots_9_1$ols_yi_res, nb2listw(w))
#Residuals are not correlated

sp_err_yi_9_1 <- errorsarlm(BuAc ~ mean+skew+PlotCRR+sd, data = plots_9_1, listw = nb2listw(w), zero.policy = T)
summary(sp_err_yi_9_1)
plots_9_1$sp_err_resi <- sp_err_yi_9_1$residuals
spplot(plots_9_1, zcol="crrsd")
print(pred_r2_9_1_sperr <- pred_r_squared(sp_err_yi_9_1))
moran.test(plots_9_1$sp_err_resi, nb2listw(w))

# sp_lag_yi_9_1 <- lagsarlm(BuAc ~ mean+sd+skew+PlotCRR, data = plots_9_1, listw = nb2listw(w), zero.policy = T)
# summary(sp_lag_yi_9_1)
# plots_9_1$sp_lag_resi <- sp_lag_yi_9_1$residuals
# spplot(plots_9_1, zcol="sp_lag_resi")
# print(pred_r2_9_1_splag <- pred_r_squared(sp_lag_yi_9_1))
# moran.test(plots_9_1$sp_lag_resi, nb2listw(w))

# Try RF regression to see if model can be improved

plots_9_1_df <- as.data.frame(plots_9_1)

set.seed(42)
train_data_indices <- rep(FALSE, nrow(plots_9_1_df))
train_data_indices[sample(1:nrow(plots_9_1_df), round(0.8 * nrow(plots_9_1_df)))] <- TRUE # randomly select 80% of the data for training
rf_regression_9_1<- randomForest(BuAc ~ mean+watersum+watermean+sd, data=plots_9_1_df[train_data_indices, ], importance=T)
rf_regression_9_1
varImpPlot(rf_regression_9_1)
pred_yield <- predict(rf_regression_9_1, plots_9_1_df[!train_data_indices,]) # predict the rings
plot(plots_9_1_df$BuAc[!train_data_indices], pred_yield, xlab="Observed", ylab="Predicted")
abline(a=0, b=1, lty=2, col=2)

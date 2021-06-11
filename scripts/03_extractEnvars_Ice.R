## Extract sea-ice variables

library(raadtools)
library(lubridate)
library(ggplot2)
library(spatialEco)
library(pals)
library(rworldmap)
library(rgeos)

setwd("~/humpbacks/sasiHumpbacks/")

# Get fitted data

dat <- readRDS("./output/fitted_tracks_w_env_01.RDS")


# Lon, lat, date for convenient extractions
xy <- dat[, c("lon", "lat")]
xyt <- dat[, c("lon", "lat", "date")]

#------------------------
# ICE
#------------------------

# Create a mask to flag values outside the NSIDC raster
# to update extracted values

# Get an ice file
msk <- readice_daily()
msk[] <- 1

# Study extent on nsidc grid for calculating distance to ice edge
study_area <- raster(ext = extent(-50, +30, -90, -35), crs = "+proj=longlat +datum=WGS84 +no_defs")
study_area <- projectRaster(study_area, res = c(25000, 25000), crs = crs(msk))

# Convert lonlat to polar stereographic
xystereo <- xy
coordinates(xystereo) <- c("lon", "lat")
crs(xystereo) <- "+proj=longlat +datum=WGS84 +no_defs"
xystereo <- spTransform(xystereo, crs(msk))

dat$stereo_lon <- coordinates(xystereo)[,1]
dat$stereo_lat <- coordinates(xystereo)[,2]

# Extract using the stereo coords
flag <- raster::extract(msk, xystereo)
dat$NSIDCflag <- as.vector(flag)
rm(flag)

# Set NAs to zero
dat[which(is.na(dat$NSIDCflag)), "NSIDCflag"] <- 0

# Get world map
data(countriesLow, package = "rworldmap")
wrld <- spTransform(countriesLow, crs(msk))

rm(msk)
rm(xystereo)

#------------------------
# Test appropriate ice concentration lag
dat$ICE_0w <- raadtools::extract(x = readice_daily, y = xyt)

# Create lags
iceLag <- function(lonlatdate, lag = 0, fun = readice_daily, ...) {
  # lonlatdate = 3 column matrix with longitude, latitude and datetime
  # lag = Period in weeks to calculate lag
  # fun = function that should be used for extraction
  xyt_lag <- lonlatdate
  xyt_lag$date <- ymd_hms(xyt_lag$date)
  xyt_lag$date <- xyt_lag$date %m-% weeks(lag) # Using %m-% instead of minus rolls back to last day of month
  icelag <- raadtools::extract(x = fun, y = xyt_lag, ...)
  return(icelag)
}

# Extract
# 1 week
dat$ICE_1w <- iceLag(lonlatdate = xyt, lag = 1, fun = readice_daily)

# 2 weeks
dat$ICE_2w <- iceLag(lonlatdate = xyt, lag = 2, fun = readice_daily)

# 1 month
dat$ICE_4w <- iceLag(lonlatdate = xyt, lag = 4, fun = readice_daily)

# 2 months
dat$ICE_8w <- iceLag(lonlatdate = xyt, lag = 8, fun = readice_daily)

# 4 months
dat$ICE_16w <- iceLag(lonlatdate = xyt, lag = 16, fun = readice_daily)

# Update NA values that were off the NSIDC grid
# to have zero sea ice concentration
dat[dat$NSIDCflag == 0 & is.na(dat$ICE_0w), ]$ICE_0w <- 0
dat[dat$NSIDCflag == 0 & is.na(dat$ICE_1w), ]$ICE_1w <- 0
dat[dat$NSIDCflag == 0 & is.na(dat$ICE_2w), ]$ICE_2w <- 0
dat[dat$NSIDCflag == 0 & is.na(dat$ICE_4w), ]$ICE_4w <- 0
dat[dat$NSIDCflag == 0 & is.na(dat$ICE_8w), ]$ICE_8w <- 0
dat[dat$NSIDCflag == 0 & is.na(dat$ICE_16w), ]$ICE_16w <- 0

# Correlations
cor(dat$g, dat$ICE_0w, use = "complete.obs")
cor(dat$g, dat$ICE_1w, use = "complete.obs")
cor(dat$g, dat$ICE_2w, use = "complete.obs")
cor(dat$g, dat$ICE_4w, use = "complete.obs")
cor(dat$g, dat$ICE_8w, use = "complete.obs")
cor(dat$g, dat$ICE_16w, use = "complete.obs")

# GAM
library(mgcv)
m0 <- gam(g ~ s(ICE_0w), data = dat, family = "betar", REML = F)
m1 <- gam(g ~ s(ICE_1w), data = dat, family = "betar", REML = F)
m2 <- gam(g ~ s(ICE_2w), data = dat, family = "betar", REML = F)
m4 <- gam(g ~ s(ICE_4w), data = dat, family = "betar", REML = F)
m8 <- gam(g ~ s(ICE_8w), data = dat, family = "betar", REML = F)
m16 <- gam(g ~ s(ICE_16w), data = dat, family = "betar", REML = F)

aic <- AIC(m0, m1, m2, m4, m8)
aic

aic[order(aic$AIC, decreasing = F), ]

# Summary table, including
# Deviance explained
# R^2
# RMSE

hold <- aic
hold$model <- row.names(aic)

for (i in 1:length(hold$model)) {
  this.mod <- eval(parse(text = hold$model[i]))
  hold$r_squared[i] <- round(summary(this.mod)$r.sq, 2)
  hold$dev_explained[i] <- round(summary(this.mod)$dev.expl*100, 2)
  preds <- predict(this.mod, newdata = dat, type = "response")
  hold$rmse[i] <- sqrt(mean((dat$g - preds)^2, na.rm = T))
  rm(preds, this.mod)
}

hold <- hold[order(hold$AIC, decreasing = F), ]
hold$AIC <- round(hold$AIC, 2)
hold$df <- round(hold$df, 2)
hold$rmse <- round(hold$rmse, 3)
hold$dAIC <- hold$AIC - min(hold$AIC)
hold

write.csv(hold, "./output/gams/ice_gams_summary.csv", row.names = F)

# Plot
tiff("./output/figs/g_v_iceconc_4w.tiff", width = 5, height = 5, units = "in", res = 300)
ggplot(data = dat, aes(x = ICE_4w, y = g)) +
  geom_point() +
  geom_smooth() +
  geom_vline(xintercept = 15, colour = "red") +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 1)) +
  labs(x = "Ice concentration 1 month prior (%)", y = "Move persistence") +
  theme_bw()
dev.off()

tiff("./output/figs/g_v_iceconc_0w.tiff", width = 5, height = 5, units = "in", res = 300)
ggplot(data = dat, aes(x = ICE_0w, y = g)) +
  geom_point() +
  geom_smooth() +
  geom_vline(xintercept = 15, colour = "red") +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 1)) +
  labs(x = "Ice concentration (%)", y = "Move persistence") +
  theme_bw()
dev.off()


#------------------------
# Distance to ice-edge
# dat$ICEDIST_0w <- raadtools::extract(x = distance_to_ice_edge,
#                                   y = xyt)
# 
# dat$ICEDIST_4w <- iceLag(lonlatdate = xyt,
#                          lag = 4,
#                          fun = distance_to_ice_edge)
# 
# # Convert to km
# dat$ICEDIST_0w <- dat$ICEDIST_0w/1000
# dat$ICEDIST_4w <- dat$ICEDIST_4w/1000

# Plot
# tiff("./out/figs/g_v_icedist.tiff", width = 5, height = 5, units = "in", res = 300)
# ggplot(data = dat, aes(x = ICEDIST_0w, y = g)) +
#   geom_point() +
#   geom_smooth() +
#   geom_vline(xintercept = 100, colour = "red") +
#   scale_x_continuous(expand = c(0, 0)) +
#   scale_y_continuous(expand = c(0, 0), limits = c(0, 1)) +
#   labs(x = "Distance to ice edge (km)", y = "Move persistence") +
#   theme_bw()
# dev.off()

#------------------------
# Lagged ice params

# Assign index for later reassembly
dat$dx <- 1:nrow(dat)

# Date - day only
dat$date_day <- as_date(ymd_hms(dat$date))

# Unique days
dts <- sort(unique(dat$date_day))

# Create a reference grid for calculating distance to ice edge

#-----------
hold <- data.frame()

for (i in 1:length(dts)) {
  print(paste(i, "of", length(dts)), sep = " ")
  # Dates for the preceding month
  this.date <- dts[i]
  start.date <- this.date %m-% months(1)
  these.dates <- seq(start.date, this.date, by = "1 day")
  
  # Data for this date
  this.data <- dat[dat$date_day == this.date, ]
  
  if (sum(this.data$NSIDCflag) > 0) { # Run only if there are locs on the NSIDC grid 
    # Calculate & extract
    # Get daily ice for the preceding two months
    this.ice <- stack(readice_daily(date = these.dates))
    
    # Calculate cv
    this.cv <- calc(this.ice, fun = cv)
    this.data$ICECV <- raster::extract(this.cv, this.data[ , c("stereo_lon", "stereo_lat")])
    # Set NA values to zero
    this.data[which(is.na(this.data$ICECV)), "ICECV"] <- 0
    
    # Calculate slope of sea ice concentration using Kendall's test
    this.melt <- raster.kendall(this.ice)
    this.data$ICETREND <- raster::extract(this.melt, this.data[ , c("stereo_lon", "stereo_lat")])
    
    # Map for an animation
    ice.now <- readice_daily(this.date)
    edge <- rasterToContour(ice.now, levels = 15)
    
    #...& extract distance to ice contour beyond nsidc grid
    this.icedist <- study_area
    dd <- gDistance(edge, as(this.icedist, "SpatialPoints"), byid = TRUE)
    this.icedist[] = apply(dd,1,min)
    this.data$ICEDISTfull <- raster::extract(this.icedist, this.data[ , c("stereo_lon", "stereo_lat")])/1000
    
    # Plot
    png(paste0("./output/figs/icemaps/", this.date, ".png"), width = 1000, height = 1000) 
    plot(ice.now,
         col = rev(brewer.blues(100)),
         main = this.date,
         legend.shrink = 1.0,
         axes = F,
         box = F)
    plot(edge, add = T, col = "#33BBEE")
    points(this.data$stereo_lon, this.data$stereo_lat, pch = 19, col = "#EE7733", cex = 1.0)
    plot(wrld, col = "grey", add = T, border = FALSE)
    dev.off()
    
    rm(this.ice, this.cv, this.melt, ice.now, edge, this.icedist, dd)
    
  } else {
    this.data$ICECV <- 0
    this.data$ICETREND <- 0
    
    # Still need to extract icedist regardless
    ice.now <- readice_daily(this.date)
    edge <- rasterToContour(ice.now, levels = 15)
    this.icedist <- study_area
    dd <- gDistance(edge, as(this.icedist, "SpatialPoints"), byid = TRUE)
    this.icedist[] = apply(dd,1,min)
    this.data$ICEDISTfull <- raster::extract(this.icedist, this.data[ , c("stereo_lon", "stereo_lat")])/1000
    rm(ice.now, edge, this.icedist, dd)
  }
  
  hold <- rbind(hold, this.data)
  rm(this.data)
  
}

# To animate using ffmpeg in console:
# ffmpeg -y -r 1 -pattern_type glob -i "*.png" -vb 8192k -s 1100x1080 -vcodec mpeg4 ice.mp4

# Reorder
hold <- hold[order(hold$dx), ]

# Update values off the NSIDC grid
hold[hold$NSIDCflag == 0 & is.na(hold$ICECV), ]$ICECV <- 0
hold[hold$NSIDCflag == 0 & is.na(hold$ICETREND), ]$ICETREND <- 0

# Update ice distances to be negative inside the sea ice
hold[!is.na(hold$ICE_0w) & hold$ICE_0w >= 15, ]$ICEDISTfull <- hold[!is.na(hold$ICE_0w) & hold$ICE_0w >= 15, ]$ICEDISTfull*-1

# Save
saveRDS(hold, "./output/fitted_tracks_w_env_02.RDS")

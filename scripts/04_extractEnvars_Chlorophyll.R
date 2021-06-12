## Extract chlorophyll

library(raster)
library(lubridate)
library(ggplot2)
library(mgcv)

setwd("D:/UCSC/Analysis/sasiHumpbacks/")

# Source the function to get chlorophyll file
source("./scripts/99_build_chlorophyll_dataset.R")

# Get fitted data
dat <- readRDS("./output/fitted_tracks_w_env_02.RDS")

#------------------------
# Get chlorophyll

# Assign index for later reassembly
dat$dx <- 1:nrow(dat)

# Date - unique year & month
dat$ymd <- paste(year(dat$date), month(dat$date), 15, sep = "-")
dts <- sort(unique(dat$ymd))

#-----------
hold <- data.frame()

for (i in 1:length(dts)) {
  print(paste(i, "of", length(dts)), sep = " ")

  this.date <- dts[i]
  
  # Get this month's data
  this.data <- dat[dat$ymd == this.date, ]
  
  # Geth the chla
  this.chla <- get_chla(date = this.date, lag = 0, datadir = "../mega/data_chlorophyll_download/")
  this.chla.lag1 <- get_chla(date = this.date, lag = 1, datadir = "../mega/data_chlorophyll_download/")
  this.chla.lag2 <- get_chla(date = this.date, lag = 2, datadir = "../mega/data_chlorophyll_download/")
  this.chla.lag3 <- get_chla(date = this.date, lag = 3, datadir = "../mega/data_chlorophyll_download/")
  
  this.data$CHLA_0m <- raster::extract(this.chla, this.data[ , c("lon", "lat")])
  this.data$CHLA_1m <- raster::extract(this.chla.lag1, this.data[ , c("lon", "lat")])
  this.data$CHLA_2m <- raster::extract(this.chla.lag2, this.data[ , c("lon", "lat")])
  this.data$CHLA_3m <- raster::extract(this.chla.lag3, this.data[ , c("lon", "lat")])
  
  hold <- rbind(hold, this.data)
  
  rm(this.date, this.data, this.chla, this.chla.lag1, this.chla.lag2, this.chla.lag3)
  
}

#------------------------
# Reorganise and save output

# Reorder
hold <- hold[order(hold$dx), ]

# Drop unneccessary columns
hold$dx <- NULL
hold$ymd <- NULL

# Save
saveRDS(hold, "./output/fitted_tracks_w_env_03.RDS")

#------------------------
# Check the two lags
m0 <- gam(g ~ s(CHLA_0m), data = hold, family = "betar", REML = F)
m1 <- gam(g ~ s(CHLA_1m), data = hold, family = "betar", REML = F)
m2 <- gam(g ~ s(CHLA_2m), data = hold, family = "betar", REML = F)
m3 <- gam(g ~ s(CHLA_3m), data = hold, family = "betar", REML = F)

AIC(m0, m1, m2, m3)

#------------------------
# Plot chlrophyll by breeding stock
png("./output/figs/CHLA.png", height = 1200, width = 1200)
ggplot(data = hold, aes(x = CHLA_0m, y = g, group = stock)) +
  geom_point(alpha = 0.1) +
  scale_y_continuous(limits = c(0, 1)) +
  scale_x_log10() +
  facet_wrap(~ stock, ncol = 2) +
  geom_rug(alpha = 0.1) +
  geom_smooth() +
  labs(title = "CHLA", x = "Chlorophyll a concentration (mg/m^3)", y = "Move persistence")
dev.off()

#------------------------
# Plot example chlorophyll for supplement
library(SOmap)
library(pals)

# Get example chlorophyll
example_chla <- get_chla(date = "2016-03-15", lag = 0)

# Crop
example_chla <- crop(example_chla, extent(-180, +180, -90, -40))

# And log
example_chla <- log(example_chla)

# Basemap
mp <- SOmap(
  trim = -40,
      bathy_legend = F,
      border_col = c("white", "white"),
      border_width = 0.01,
      straight = TRUE,
      graticules = TRUE)

# Get rid of bathymetry and its legend
raster::values(mp$bathy[[1]]$plotargs$x) <- NA_real_
mp$bathy_legend <- NULL

# Plot chlorophyll to file
png("./out/figs/envar_examples/CHLA.png",
     height = 4,
     width = 6,
     units = "in",
    res = 300)

mp

SOplot(SOproj(example_chla), col = rev(ocean.algae(125)), add = T,
       legend.args = list(text = "log10(CHLA)"))

dev.off()

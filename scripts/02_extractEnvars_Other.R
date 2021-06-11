## Extract other environmental covariates

library(raster)
library(raadtools)
library(lubridate)
# library(boaR)
library(grec)

setwd("~/humpbacks/sasiHumpbacks/")

# -----------------------
# Get fitted data
dat <- readRDS("./output/fitted_tracks.RDS")

# Lon, lat, date for convenient extractions
xy <- dat[, c("lon", "lat")]
xyt <- dat[, c("lon", "lat", "date")]

#------------------------
# Bathymetry
#------------------------
dat$DEPTH <- raadtools::extract(x = readtopo,
                              y = xy,
                              topo = "gebco_19")

dat$DEPTH[dat$DEPTH > 0] <- 0

#------------------------
# Slope & Terrain Ruggedness
#------------------------
dep <- readtopo(topo = "gebco_19", xylim = c(-50, +30, -70, -35))
slp <- terrain(dep, opt = "slope", unit = "degrees")
tri <- terrain(dep, opt = "TRI")

dat$SLOPE <- raster::extract(slp, xy)
dat$TRI <- raster::extract(tri, xy)

rm(dep, slp, tri)

#------------------------
# Distance to shelf
#------------------------
shlf <- readderivaadc("distance_shelf")
dat$SHELFDIST <- raster::extract(shlf, xy)

rm(shlf)

#------------------------
# Distance to upper slope (Antarctica)
#------------------------
slp <- readderivaadc("distance_upper_slope")
dat$SLOPEDIST <- raster::extract(slp, xy)

rm(slp)

#------------------------
## SST
#------------------------
dat$SST <- raadtools::extract(x = readghrsst,
                              y = xyt,
                              latest = FALSE)
dat$SST <- dat$SST - 273.15 # k to c

#------------------------
## SSH
#------------------------

# Done below in loop, environmental data is in [0:360]
# dat$SSH <- raadtools::extract(x = readssh,
#                               y = xyt,
#                               latest = FALSE,
#                               ssha = TRUE)

#------------------------
## EKE
#------------------------
dat$CURU <- raadtools::extract(x = readcurr,
                              y = xyt,
                              latest = FALSE,
                              uonly = TRUE)

dat$CURV <- raadtools::extract(x = readcurr,
                               y = xyt,
                               latest = FALSE,
                               vonly = TRUE)

dat$EKE <- 0.5*(dat$CURU^2 + dat$CURV^2)

#------------------------
## WIND
#------------------------
dat$WINDU <- raadtools::extract(x = readwind,
                               y = xyt,
                               latest = FALSE,
                               uonly = TRUE)

dat$WINDV <- raadtools::extract(x = readwind,
                                y = xyt,
                                latest = FALSE,
                                vonly = TRUE)

#------------------------
## SSTGRAD and SSHGRAD
#------------------------

# Loop through days, calculate and extract spatial gradients of sst and ssh

#------------------------
# Assign index for later reassembly
dat$dx <- 1:nrow(dat)

# Date - day only
dat$date_day <- as_date(ymd_hms(dat$date))

# Unique days
dts <- sort(unique(dat$date_day))

#-----------

hold <- data.frame()

for (i in 1:length(dts)) {
  print(paste(i, "of", length(dts)), sep = " ")

  this.date <- dts[i]
  
  # Data for this date
  this.data <- dat[dat$date_day == this.date, ]
  
  # Extent
  this.ext <- extent(min(this.data$lon)-1, max(this.data$lon)+1, min(this.data$lat)-1, max(this.data$lat)+1)
  
  # Get sst and calculate gradient
  this.sst <- readghrsst(date = this.date, xylim = this.ext)
  this.sst <- this.sst - 273.15 # k to c
  this.sst.grad <- raster::terrain(this.sst, opt = "slope", unit = "radians")
  this.sst.grad <- calc(this.sst.grad, fun = tan)*1000
  this.data$SSTGRAD <- raster::extract(this.sst.grad, this.data[ , c("lon", "lat")])
  
  # Get ssh and calculate gradient
  this.ssha <- rotate(readssh(date = this.date, ssha = T))
  this.data$SSH <- raster::extract(this.ssha, this.data[ , c("lon", "lat")])
  
  this.ssh <- rotate(readssh(date = this.date, ssha = F))
  this.ssh <- raster::terrain(this.ssh, opt = "slope", unit = "radians")
  this.ssh <- calc(this.ssh, fun = tan)*1000
  this.data$SSHGRAD <- raster::extract(this.ssh, this.data[ , c("lon", "lat")])
  
  # sst fronts
  
  # Unsure why the strange flipping and transposition is needed to get boa output correct
  # sst.fronts <- boa(lon = xFromCol(this.sst),
  #                   lat = yFromRow(this.sst),
  #                   ingrid = t(apply(raster::as.matrix(this.sst), 2, rev)))
  
  # Use library 'grec' instead
  sst.fronts <- detectFronts(this.sst, method = "BelkinOReilly2009")

  this.data$SSTFRONT<- raster::extract(sst.fronts, this.data[ , c("lon", "lat")])
  
  # Bind result
  if (nrow(this.data) > 0) {
  hold <- rbind(hold, this.data)
  }
  
  # Clean up
  rm(this.date, this.data, this.ext, this.sst, this.sst.grad, this.ssh, this.ssha, sst.fronts)
  
}

dat <- hold

# Reorder
dat <- dat[order(dat$dx), ]

# Reorganize
dat$dx <- NULL
dat$date_day <- NULL

# Save
dat <- saveRDS(dat, "./output/fitted_tracks_w_env_01.RDS")

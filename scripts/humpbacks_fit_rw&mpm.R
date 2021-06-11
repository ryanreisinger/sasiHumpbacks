#FoieGras
# setwd("C:/R/Humpbackwhalemixing")

setwd("C:\\Users\\Ryan Reisinger\\Documents\\Academic\\UCSC\\Work\\Analysis\\sasiHumpbacks\\")

#library(remotes)
library(tidyverse)
#remotes::install_github("ianjonsen/foieGras")
library(foieGras)
library(sf)
#install.packages("rnaturalearth")
library(rnaturalearth)
#install.packages("rnaturalearthdata)
library(rnaturalearthdata)
library(wesanderson)
library(patchwork)

#---------------------------------------
# Load all three datasets and combine them

# Brazil BSA
this_data <- read.csv("./data/Brazil_HW_S_of_40_Mdu_5.csv", sep = ";", stringsAsFactors = F)

A <- data.frame("stock" = "A",
                 "id" = this_data$ï..ID,
                 "date" = strptime(paste(this_data$Date, this_data$Time, sep = " "), format = "%Y/%m/%d %H:%M:%S"),
                 "lc" = this_data$Quality,
                 "lon" = this_data$Longitude,
                 "lat" = this_data$Latitude
)

# South Africa BSB2
this_data <- read.csv("./data/BSB2.csv", sep = ";", stringsAsFactors = F)

B2 <- data.frame("stock" = "B2",
                 "id" = this_data$ID,
                 "date" = strptime(paste(this_data$Date, this_data$Time, sep = " "), format = "%Y/%m/%d %H:%M:%S"),
                 "lc" = this_data$Quality,
                 "lon" = this_data$Longitude,
                 "lat" = this_data$Latitude
)

# South Africa BSC
this_data <- read.csv("./data/BSC.csv", sep = ";", stringsAsFactors = F)

C <- data.frame("stock" = "C",
                 "id" = this_data$ID,
                 "date" = strptime(paste(this_data$Date, this_data$Time, sep = " "), format = "%Y/%m/%d %H:%M:%S"),
                 "lc" = this_data$Quality,
                 "lon" = this_data$Longitude,
                 "lat" = this_data$Latitude
)

tracks <- rbind(A, B2, C)

# Remove any invalid or missing dates
tracks <- filter(tracks, !is.na(tracks$date))

#--------------------------------------------------------------------------
## Split tracks with large gaps
## into segments

int.thresh <- 3 # Gap threshold in days

tracks$id_original <- tracks$id

ids <- unique(tracks$id)
all.d <- data.frame()

for (i in 1:length(ids)) {
  this.id <- ids[i]
  sub <- tracks[tracks$id == this.id, ]
  intervals <- diff(sub$date)
  units(intervals) <- "days"
  sub$int <- c(0, intervals)
  sub$dx <- 0
  sub[sub$int > int.thresh, "dx"] <- 1
  sub$dx2 <- cumsum(sub$dx)
  sub$id <- paste0(sub$id, "_segment", sub$dx2)
  sub$int <- NULL
  sub$dx <- NULL
  sub$dx2 <- NULL
  all.d <- rbind(all.d, sub)
}

tracks <- all.d
rm(all.d)

# ----------------------
## Filter again to remove fragments
nlocs <- tracks %>%
  group_by(id) %>%
  tally %>%
  filter(., n > 2)

tracks <- dplyr::filter(tracks, tracks$id %in% nlocs$id)

#---------------------------------------
# Fit the random walk model

# First try one by one to identify problems
ids <- unique(tracks$id)

these_fits <- data.frame()

for( i in ids) {
  print(i) # So we can see where it fails
  this_track <- dplyr::filter(tracks, id == i)
  this_fit <- fit_ssm(dplyr::select(this_track, id, date, lc, lon, lat),
                 model = "rw", # Which model
                 vmax = 10, # Max speed in m/s
                 ang = c(15, 25), # Angle
                 distlim = c(2500, 5000), # Spike distance in m
                 time.step = 12 # Time step in hours
                 
  )
  
  if(this_fit$converged) {
  this_grab <- grab(this_fit, what = "predicted", as_sf = FALSE)
  these_fits <- rbind(these_fits, this_grab)
  } else {
    print("FAILED")
  }
}

# Plot to check
world <- ne_countries(scale = "medium", returnclass = "sf") #define world data

ggplot(data = world) +
  geom_sf(color = "grey30", fill = "grey30") +
  geom_point(data = these_fits, aes(x = lon, y = lat)) +
  theme_bw() +
  coord_sf(crs= "+proj=longlat +datum=WGS84",
           xlim = c(min(these_fits$lon), max(these_fits$lon)),
           ylim = c(min(these_fits$lat), max(these_fits$lat)))
  
#---------------------------------------
# Fit the move persistence model
# mpm <- fit_mpm(x = these_fits,
#                model = "mpm",
#                verbose = 1)

# First try one by one to identify problems
ids <- unique(these_fits$id)

these_mpm <- data.frame()

for(i in ids) {
  print(i) # So we can see where it fails
  this <- dplyr::filter(these_fits, id == i)
  # Fit mpm, ignoring errrors
  tryCatch({
    print(i)
    this_mpm <- fit_mpm(x = this,
                        model = "mpm",
                        verbose = 1)
    
    this_mpm_out_dat <- grab(this_mpm, what = "data", as_sf = F)
    this_mpm_out_fit <- grab(this_mpm, what = "fitted", as_sf = F)
    this_grab <- cbind(this_mpm_out_dat[ , c("id", "date", "lon", "lat")], this_mpm_out_fit[ , c("g", "g.se")])
  }, error=function(e){})
  
  if (nrow(this_grab) > 0) {
    these_mpm <- rbind(these_mpm, this_grab)
  } else {
    print("FAILED")
  }
}

# Put the segments back together
these_mpm$id <- substr(x = these_mpm$id, start = 1, stop = nchar(these_mpm$id)-9)

# Plot
p1 <- ggplot(data = world) +
  geom_sf(color = "grey30", fill = "grey30") +
  geom_path(data = these_mpm, aes(x = lon, y = lat, group = id), colour = "grey") +
  geom_point(data = these_mpm, aes(x = lon, y = lat, group = id, colour = g)) +
  theme_bw() +
  coord_sf(crs= "+proj=longlat +datum=WGS84",
           xlim = c(min(these_mpm$lon), max(these_mpm$lon)),
           ylim = c(min(these_mpm$lat), max(these_mpm$lat))) +
  scale_colour_gradientn(colours = rev(wes_palette(name = "Zissou1", 
                                                   type = "continuous")), name = expression(gamma[t]), limits = c(0, 1))

p1

#---------------------------------------
# Add the breeding stock information
out <- dplyr::left_join(x = these_mpm, y = select(tracks, stock, id_original), by = c("id" = "id_original"))

# Plot by stock
p2 <- ggplot(data = world) +
  geom_sf(color = "grey30", fill = "grey30") +
  geom_path(data = out, aes(x = lon, y = lat, group = id), colour = "grey") +
  geom_point(data = out, aes(x = lon, y = lat, group = id, colour = stock)) +
  theme_bw() +
  coord_sf(crs= "+proj=longlat +datum=WGS84",
           xlim = c(min(out$lon), max(out$lon)),
           ylim = c(min(out$lat), max(out$lat)))

p2

# Save the two plots together
p <- p1/p2
png("./output/tracks.png", width = 8, height = 8, unit = "in", res = 300)
plot(p)
dev.off()

#---------------------------------------
# Save
saveRDS(out, "./output/fitted_tracks.RDS")
















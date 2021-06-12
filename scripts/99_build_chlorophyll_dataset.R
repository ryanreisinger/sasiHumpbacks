# Build chlorophyll dataset

## Ryan Reisinger
## April 2020


library(raster)

#-------------------------------
## Download chlorophyll data with bowerbird and blueant
# Set up bowerbird
if (FALSE) {
  
  library(bowerbird)
  library(blueant)
  
  my_directory <- "D:/UCSC/Analysis/mega/data_chlorophyll/"
  cf <- bb_config(local_file_root = my_directory)
  
  # Define the source
  # src <- sources(name="Oceandata SeaWiFS Level-3 mapped monthly 9km chl-a")
  src <- sources(name="Oceandata VIIRS Level-3 mapped 32-day 9km chl-a")
  src$user <- "ryan.r.reisinger@gmail.com"
  src$password <- "xrVW4b8S748f56n"
  cf <- bb_add(cf, src)
  
  # Download
  status <- bb_sync(cf, verbose = TRUE) # Sync
}

#-------------------------------
## List the files

# Directory
datadir <- "D:/UCSC/Analysis/mega/data_chlorophyll_download/"

# List the .nc files in the data directory
chl_files <- list.files(datadir)[which(grepl("*.nc", list.files(datadir)))]

#-------------------------------
# Construct a reference table
chl_ref <- data.frame("filename" = chl_files)

# Start and end date
chl_ref$start_date <- strptime(substr(chl_ref$filename, 2, 8), format = "%Y%j", tz = "UTC")
chl_ref$end_date <- strptime(substr(chl_ref$filename, 9, 15), format = "%Y%j", tz = "UTC")


#-------------------------------
# Function to read in the correct file, given a date

get_chla <- function(date = NULL, datadir = "D:/UCSC/Analysis/mega/data_chlorophyll_download/", lag = 0) {
  
  require(raster)
  require(lubridate)
  
  # date = date as text in ymd format
  # datadir = directory containing chlorophyll files
  # lag = lag in months, 0 if matching months are required
  
  dt <- strptime(date, format = "%Y-%m-%d", tz = "UTC")
  
  # If a lag is specified, subtract
  if (lag == 0) {
    dt <- dt
  } else {
    dt <- dt %m-% months(lag)
  }
  
  if(!is.na(dt)) {
    which.file <- as.character(chl_ref[dt >= chl_ref$start_date & dt <= chl_ref$end_date, ]$filename)
    if (length(which.file) == 0) {
      stop("No matching data found")
    }
    which.file <- paste0(datadir, which.file)
    this.chla <- raster(which.file)
    return(this.chla)
  } else {
    stop("Not a valid date")
  }
}

# tst <- get_chla("2015-06-01", lag = 1)



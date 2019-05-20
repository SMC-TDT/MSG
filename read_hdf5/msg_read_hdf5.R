#!/usr/bin/env Rscript

## msg_read_hdf5.R ####################################################################################

## Created : Carles Cayuela, June 2014; Anna del Moral, July 2015
## Debugged : Patricia Altube, May 2019
## Last modified : Patricia Altube, May 2019

## Description : Reads and displays MSG data in hdf5 format (returns: png image and geotiff file)

###################################################################################################

name_scr <- "msg_read_hdf5.R"
name_fconf <- "msg_read_hdf5_config.R"
name_ffuns <- "msg_read_hdf5_funs.R"

## INITIALIZATION #################################################################################

# Get the path of the current script
args <- commandArgs(trailingOnly = FALSE)
script_dir <- NULL

if (any(grep("--slave", args))){
  file <- gsub("--file=", "", args[grep("--file=", args)])
  script_dir <- gsub(name_scr, "", file)
}else if (any(grep("--interactive", args))){
  try(script_dir <- dirname(sys.frame(1)$ofile), silent=TRUE)
}

# Load functions and user configuration
file_funs <- paste(script_dir, name_ffuns, sep="/")
file_conf <- paste(script_dir, name_fconf, sep="/")

if (!file.exists(file_funs)){
  stop("Functions file not found")
}
if (!file.exists(file_conf)){
  stop("Configuration file not found")
}

source(file_funs)
source(file_conf)

## INPUT ARGUMENTS ################################################################################

# FOR TESTING from RStudio
# args <- list("file_in"=paste0("/home/pav/Desktop/MSG_hdf5/data/",
#                               "MSG3-SEVI-MSG15-0100-NA-20171018165740.354000000Z-20171018165757-1303506.h5"),
#              "path_out"="/home/pav/Desktop/MSG_hdf5/img")

args <- parse_cl()
file_in <- args$file_in
path_out <- args$path_out

if (is.null(file_in)){
  stop("Input file not specified")
}
if (!file.exists(file_in)){
  stop("Input file not found")
}

if (is.null(path_out)){
  path_out <- paste0(script_dir, "/img")
}

###################################################################################################

# Available color palettes
col_palettes <- list("greyscale"=gray.colors(255, start=1, end=0, gamma=1.2, alpha=NULL),
                     "rainbow"=rainbow(255), "terrain"=terrain.colors(255), 
                     "topo"=topo.colors(255))

## DATA ###########################################################################################

# Map from "rworldmap" library
map_pol <- getMap(resolution = "low")
e_map <- extent(ext_img)
suppressWarnings(map_pol <- crop(map_pol, e_map))

# Scan type from file name (TO BE MODIFIED):
name_f <- gsub("(\\/[[:graph:]]*\\/)*", "",file_in)
scan <- "NOM"
if (grepl("RSS", name_f)){
  scan <- "RSS"
}

# Read hdf5 file
data_all <- h5read(file=file_in, name=data_name_h5, compoundAsDataFrame=FALSE)

data_data <- data_all[["DATA"]]
data_meta <- data_all[["METADATA"]]
data_hdr <- data_meta[["HEADER"]]
data_sub <- data_meta[["SUBSET"]]

# Channel number
ch <- as.numeric(get_channel_msg(data_sub))
name_ch <- paste("Channel", sprintf("%02d", ch))

# Channel parameters
ch_pars <- channels[[as.character(ch)]]
prj <- projs[[scan]]
wnum <- 10^4/ch_pars[["wl"]] # wavenumber

# Necessary description lists
img_acq_descr <- list_entry_value(data_hdr[["ImageAcquisition"]][["ImageAcquisition_DESCR"]])
img_descr <- list_entry_value(data_hdr[["ImageDescription"]][["ImageDescription_DESCR"]])
geom_descr <- list_entry_value(data_hdr[["GeometricProcessing"]][["GeometricProcessing_DESCR"]])
subset_descr <- list_entry_value(data_sub)

# Coordinate transformation parameters
rect_off <- transf_cnt[[ch_pars[["id"]]]]
rect_pars <- subset_descr[grep(names(subset_descr), pattern=ch_pars[["id"]])]
grid_pars <- img_descr[grep(names(img_descr), pattern=ch_pars[["id"]])]

rect_lin_S <- as.numeric(rect_pars[grep("SouthLine", names(rect_pars))])
rect_col_E <- as.numeric(rect_pars[grep("EastColumn", names(rect_pars))])
lon_ssp <- img_descr[["ProjectionDescription-LongitudeOfSSP"]]*(pi/180) # Longitude of sub-satellite point

# Data for calculation of radiance
cal_arr <- data_hdr[["RadiometricProcessing"]][["Level15ImageCalibration_ARRAY"]]
data_img <- data_data[[name_ch]][["IMAGE_DATA"]]

# Calibration slope and offset
cal_slope <- cal_arr[["Cal_Slope"]][ch]
cal_offset <- cal_arr[["Cal_Offset"]][ch]

# Date and time
days <- img_acq_descr[["PlannedAcquisitionTime-TrueRepeatCycleStart-days"]]
ms <- img_acq_descr[["PlannedAcquisitionTime-TrueRepeatCycleStart-milliseconds"]]
dt <- get_datetime_msg(days, ms, date_ref=tref_CCSDS)

## COMPUTATIONS ###################################################################################

# Output title and text
tit_out <- paste(paste0(scan, sprintf("%02d", ch)), format(dt, format="%Y-%m-%d %H:%M:%S")) 
name_out <- paste("MSG", format(dt, format="%Y%m%d%H%M%S"), 
                  paste0(scan, sprintf("%02d", ch)), sep="_")

#Output files
path_out_tiff <- paste(path_out, "tiff", sep="/")
create_dir(path_out_tiff)
png_out <- paste(path_out, paste0(name_out, ".png"), sep="/")
tif_out <- paste(path_out_tiff, paste0(name_out, ".tiff"), sep="/")

# Radiance (transform digital number to radiance)
rad <- cal_offset + (cal_slope*data_img) 

# Pixel coordinates to lat/lon coordinates 
rad_rast <- t(raster(rad))
rad_arr <- as.matrix(rad_rast)

l <- row(rad_arr) + rect_lin_S - 1
c <- col(rad_arr) + rect_col_E - 1

coord_geo <- pixel2geo(col=c, row=l, coff=rect_off[["coff"]], roff=rect_off[["loff"]], 
                       cfac=rect_off[["cfac"]], rfac=rect_off[["lfac"]], 
                       lon_ssp=lon_ssp, h_sat=h_msg)

# Solar position
sun_pars <- solar_el_and_dist(dt=dt, coord_geo=coord_geo)
sun_dist <- sun_pars$d_avg
el_sun <- sun_pars[[paste0("el_", el_sun_test)]]
  
# Display variable (albedo for VIS and TB for IR):
if (ch_pars[["vis"]]){
  
  if(el_sun>=el_sun_min){
    irad_sun <- ch_pars[["irad"]]/(pi*(sun_dist^2))
    img_var <- 100*rad/irad_sun # reflectance (in percent)
    img_var[Which(img_var<0.1, cells=T)] <- 0.1
    
  }else{
    print(paste("Channel", ch, ": The Sun is too low to calculate reflectance"))
    stopifnot(sun_pos$altitude>=30)
  }
  
}else{
  img_var <- (c2*wnum)/log(1 + (c1*wnum^3/rad))
}

img_rast <- t(raster(img_var))
img_arr <- as.matrix(img_rast)

# Normalized Geostationary Projection
y <- (array(coord_geo$lat))
x <- (array(coord_geo$lon))
z <- rev(array(img_arr))

xyz <- data.frame(x, y, z)
xyz <- na.omit(xyz)

img_sp <- SpatialPoints(coords=cbind(xyz$x, xyz$y), proj4string=CRS("+proj=longlat"))
img_pr <- spTransform(img_sp, CRS(prj))
e <- extent(img_pr)
extent(img_rast) <- e
projection(img_rast) <- prj

map_pol <- spTransform(map_pol, CRSobj=CRS(prj))

# Output figure
png(png_out, width=ncol(img_rast), height=nrow(img_rast))
par(mar = c(0.1, 0.1, 0.1, 0.1))
image(img_rast, col=col_palettes[[col]], axes=FALSE, xlab=NA, ylab=NA)
plot(map_pol, border=col_shp, add=T)
text(x=xmax(img_rast)-450000, y=ymin(img_rast)+60000, labels=tit_out, col=col_txt, cex=1.6)
dev.off()

# Output tiff
writeRaster(img_rast, tif_out, format="GTiff", overwrite=TRUE)

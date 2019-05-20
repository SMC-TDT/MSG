
## MSG_read_hdf5_funs.R ###########################################################################

## Created : Patricia Altube, May 2019
## Debugged : Patricia Altube, May 2019
## Last modified : Patricia Altube, May 2019

## Description : Functions used by the MSG_read_hdf5.R process
# create_dir
# get_channel_msg
# get_datetime_msg
# list_entry_value
# parse_cl
# pixel2geo
# solar_pars

## LIBRARIES ######################################################################################

library("rhdf5",  quietly=TRUE, warn.conflicts=FALSE) # Read hdf5 data
library(raster, quietly=TRUE, warn.conflicts=FALSE) # Raster operations
library(rgdal, quietly=TRUE, warn.conflicts=FALSE) # Map/data projection
library(rworldmap, quietly=TRUE, warn.conflicts=FALSE) # Overlay map
library(oce, quietly=TRUE, warn.conflicts=FALSE) # Solar position
library(optparse, quietly=TRUE, warn.conflicts=FALSE) # Command line arguments

## FUNCTIONS ######################################################################################

create_dir <- function(path){
  # -----------------------------------------------------------------------------------------------
  # Creates directory (recursively) if it does not exist
  #
  # PARAMETERS
  # **********
  # path : str
  #   directory path
  #
  # -----------------------------------------------------------------------------------------------
  
  if (!file.exists(path)){
    dir.create(path, showWarnings =F, recursive=T)
  }
}

get_datetime_msg <- function(days, ms, date_ref=as.Date("1958-01-01"), time_ref="00:00:00"){
  # -----------------------------------------------------------------------------------------------
  # Calculates date and time from days and milliseconds 
  #
  # PARAMETERS
  # **********
  # days : int
  #   Number of days since reference date
  # ms : int/float
  #   Milliseconds since time_ref
  #
  # OTHER PARAMETERS
  # ****************
  # date_ref : Date
  #   Reference date object
  # time_ref : str
  #   Reference time
  #
  # RETURNS
  # *******
  # dt : POSIXlt
  #   Date and time object
  # -----------------------------------------------------------------------------------------------
  
  dt0 <- as.POSIXlt(paste(date_ref+days, time_ref), tz="UTC")
  dt <- as.POSIXlt(ms/1000, format="%S", origin=dt0, tz="UTC")
  
  return(dt)
}

get_channel_msg <- function(meta_subset){
  # -----------------------------------------------------------------------------------------------
  # (ad-hoc function) Retrieves channel number from hdf5 MSG data (METADATA>SUBSET) 
  #
  # PARAMETERS
  # **********
  # meta_subset : list
  #   Level 5 data subset description list with 2 elements: 'Entry name' and 'Value'
  #
  # RETURNS
  # *******
  # chan : int
  #   Channel number
  # -----------------------------------------------------------------------------------------------
  
  indx <- which(meta_subset[["Value"]]=="X")
  str <- meta_subset[["EntryName"]][indx]
  chan <- regmatches(str, regexpr("\\d+", str))
  
  return(chan)
}

list_entry_value <- function(ls, name_entry="EntryName", name_val="Value"){
  # -----------------------------------------------------------------------------------------------
  # Generates a list from a list of two elements in which one element contains 
  # the names and the other the corresponding values
  #
  # PARAMETERS
  # **********
  # ls : list
  #   Two element list: [[name_entry]] with names and [[name_val]] with values
  #
  # OTHER PARAMETERS
  # ****************
  # name_entry : str
  #   Name of the list entry containing the names
  # name_val : str
  #   Name of the list entry containing the values
  #
  # RETURNS
  # *******
  # new_ls : list
  #   List with name_entry=name_val entry pairs
  # -----------------------------------------------------------------------------------------------
  
  new_ls <- as.list(ls[[name_val]])
  names(new_ls) <- ls[[name_entry]]
  
  indx_num <- grep("^\\d+$", new_ls)
  new_ls[indx_num] <- as.numeric(new_ls[indx_num])
  
  return(new_ls)
  
}

parse_cl <- function(){
  # -----------------------------------------------------------------------------------------------
  # Command-line argument parser: 
  #   Input file: -i, --file_in
  #   Output folder: -p, --path_out
  #
  # PARAMETERS
  # **********
  #
  # RETURNS
  # *******
  # opt : list(file_in, path_out)
  #   Command line arguments
  # -----------------------------------------------------------------------------------------------
  
  option_list = list(make_option(c("-i", "--file_in"), type="character", default=NULL, 
                                 help="Input file (complete path)", metavar="character"),
                     make_option(c("-p", "--path_out"), type="character", default=NULL,
                                 help="Output folder [default= ./img/]", metavar="character"))
  
  opt_parser <- OptionParser(option_list=option_list);
  opt <- parse_args(opt_parser);
  
  return(opt)
}

pixel2geo <- function(col, row, coff, roff, cfac, rfac, lon_ssp, h_sat=42164){
  # -----------------------------------------------------------------------------------------------
  # (ad-hoc function) Transforms pixel coordinates to lat/lon coordinates 
  #
  # PARAMETERS
  # **********
  # col : int (num, arr or mat)
  #   Column number
  # row : int (num, arr or mat)
  #   Row number
  # coff : int
  #   Column offset
  # roff : int
  #   Row offset
  # cfac : int
  #   Column scaling factor
  # rfac : int
  #   Row scaling factor
  # lon_ssp : float
  #   Longitude of sub-satellite point
  #
  # OTHER PARAMETERS
  # ****************
  # h_sat : float
  #   Height of satellite [km] from center of Earth
  #
  # RETURNS
  # *******
  # coord_geo : list(lat=, lon=)
  #   List with geographical coordinates
  # -----------------------------------------------------------------------------------------------
  
  cnt0 <- 2^(-16)
  cnt1 <- 1.006803
  cnt2 <- 1737121856
  
  # Intermediate coordinates (rad):
  x <- (col - coff)/(cnt0*cfac)
  y <- (row - roff)/(cnt0*rfac)
  
  # Coordinate transformation function:
  a1 <- cos(x)*cos(y)
  a2 <- cos(y)*cos(y) + cnt1*sin(y)*sin(y)
  
  sa <- (h_sat*a1)^(2) - cnt2*a2
  sd <- sqrt(sa)
  sn <- (h_sat*a1 - sd)/a2
  s1 <- h_sat - sn*a1
  s2 <- sn*sin(x)*cos(y)
  s3 <- -sn*sin(y)
  sxy <- sqrt(s1^(2) + s2^(2))
  
  # Geographic Coordinates:
  lon <- atan(s2/s1) + lon_ssp
  lat <- atan((cnt1*s3)/sxy)
  
  coord_geo <- list("lat"=lat*180/pi, "lon"=lon*180/pi)
  return(coord_geo)
}

solar_el_and_dist <- function(dt, coord_geo){
  # -----------------------------------------------------------------------------------------------
  # (ad-hoc function) Calculates solar position for given time and coordinate set and returns 
  # ad-hoc statistics of solar elevation and distance
  #
  # PARAMETERS
  # **********
  # dt : POSIXlt
  #   Date and time object
  # coord_geo : list(lat=, lon=)
  #   List with geographical coordinates (see pixel2geo function)
  #
  # RETURNS
  # *******
  # sun_pos_stat : list(el_min=, el_max=, el_avg=, el_med=, d_avg=)
  #   List with statistics of solar elevation and distance
  # -----------------------------------------------------------------------------------------------
  
  df <- data.frame(lat=array(coord_geo$lat), lon=array(coord_geo$lon))
  df$dt <- dt
  
  sun_pos <- sunAngle(df$dt, df$lon, df$lat)
  
  elev_min <- min(sun_pos$altitude, na.rm=TRUE)
  elev_max <- max(sun_pos$altitude, na.rm=TRUE)
  elev_med <- median(sun_pos$altitude, na.rm=TRUE)
  elev_avg <- mean(sun_pos$altitude, na.rm=TRUE)
  dist_avg <- mean(sun_pos$distance, na.rm=TRUE)
  
  sun_pos_stat <- list(el_min=elev_min, el_max=elev_max, 
                       el_avg=elev_avg, el_med=elev_med, d_avg=dist_avg)
  return(sun_pos_stat)
  
}

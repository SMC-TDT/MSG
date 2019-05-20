
## USER SETTINGS ########################################################################
## Set the values of the variable (right hand side of the <- sign) in this section for 
## customization.  Do not change the variable names (the left hand side of the <- sign). 

# color palette, from R Palettes{grDevices}: "greyscale", "rainbow", "terrain", "topo" 
col <- "greyscale"

# Color of overlay map and text
col_shp <- "yellow"
col_txt <- "yellow"

# Minimum solar elevation required for VIS channel processing
el_sun_min <- 30

# Which elevation of the set is tested against el_sun_min:
# min: minimum
# med: median
# avg: mean
# max: maximum
el_sun_test <- "avg"

## DEFAULT METADATA #####################################################################
## Do not modify this section unless you know what you are doing

# Extent of the image: c(lon_min, lon_max, lat_min, lat_max)
ext_img <- c(-20, 65, 10, 65) 

# Planck constants:
c1 <- 1.19104e-05 # [mW cm^4 m^-2 sr^-1]
c2 <- 1.43877 # [K cm]

# Origin date
tref_CCSDS <- as.Date("1958-01-01")

# name of the dataset in the HDF5 file
data_name_h5 <- "/U-MARF/MSG/Level1.5"

# Distance from Earth centre to Satellite
h_msg <- 42164

# Scan types:
# "NOM": nominal
# "RSS": rapid scan

# Image projections depending on scan type
projs <- list("NOM"="+proj=geos +h=35785831 +a=6378169 +b=6356583 +lon_0=0.0", 
              "RSS"="+proj=geos +h=35785831 +a=6378169 +b=6356583 +lon_0=9.5")

# Spectral channel types:
# "VIS": visible 
# "IR" : infrared 
# "HRV": High Resolution Visible 

# Channel parameters
# id: identifier without wavelength
# wl: wavelength in micrometers
# irad: solar irradiance in cm
# vis: whether the channel is visible
# transf: coordinate transformation type

channels <- list("1"=list(id="VIS", wl=0.635, irad=65.2296, vis=TRUE, transf="GEN"),
                 "2"=list(id="VIS", wl=0.810, irad=73.0127, vis=TRUE, transf="GEN"),
                 "3"=list(id="IR", wl=1.640, irad=62.3715, vis=TRUE, transf="GEN"), 
                 "4"=list(id="IR", wl=3.920, irad=NA, vis=FALSE, transf="GEN"),
                 "5"=list(id="IR", wl=6.250, irad=NA, vis=FALSE, transf="GEN"),
                 "6"=list(id="IR", wl=7.350, irad=NA, vis=FALSE, transf="GEN"),
                 "7"=list(id="IR", wl=8.700, irad=NA, vis=FALSE, transf="GEN"), 
                 "8"=list(id="IR", wl=9.660, irad=NA, vis=FALSE, transf="GEN"),
                 "9"=list(id="IR", wl=10.800, irad=NA, vis=FALSE, transf="GEN"), 
                 "10"=list(id="IR", wl=12.000, irad=NA, vis=FALSE, transf="GEN"),
                 "11"=list(id="IR", wl=13.400, irad=NA, vis=FALSE, transf="GEN"), 
                 "12"=list(id="HRV", wl=0.750, irad=NA, vis=TRUE, transf="HRV"))

# Coordinate transformation constants depending on channel id
transf_cnt <- list("VIS"=list("coff"=1856, "loff"=1856, "cfac"=-781648343, "lfac"=-781648343),
                   "IR"=list("coff"=1856, "loff"=1856, "cfac"=-781648343, "lfac"=-781648343),
                   "HRV"=list("coff"=5567, "loff"=5567, "cfac"=-2344944937, "lfac"=-2344944937))


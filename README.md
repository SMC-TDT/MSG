# MSG - Meteosat Second Generation data processing

This repository contains code related to MSG data management which has been developed by the Remote Sensing team at the Meteorological Service of Catalonia.
Currently, the repository includes:
### **read_hdf5**: 
- Description: R code for reading and displaying MSG data in hdf5 format. 
- Input: MSG hdf5 files. 
- Output: Image in png format and a geotiff file.


## Dependencies

R libraries:
- oce
- optparse
- raster
- rgdal
- rgeos
- rhdf5
- rworldmap

## Usage
### read_hdf5:

First, modify, if desired, the configuration file: msg_read_hdf5_config.R

#### Single file processing: #### 

Execute msg_read_hdf5.R:
~~~~
Rscript msg_read_hdf5.R -i /input/path/inputfile.h5 -p /output/path/
~~~~
Command line arguments:

    -i, --file_in: input hdf5 file with its complete path

    -p, --path_out: output path (will be created if it does not exist already)

  

# MSG - Meteosat Second Generation data processing

This repository contains code related to MSG data management which has been developed by the Remote Sensing team at the Meteorological Service of Catalonia.
Currently, the repository includes:
- read_hdf5: R code for reading and displaying MSG data in hdf5 format

## Dependencies

R libraries:
- rhdf5
- raster
- rgdal
- rworldmap
- oce
- optparse

## Usage
### read_hdf5:
Modify, if desired, the configuration file: msg_read_hdf5_config.R

**Single file processing:** execute msg_read_hdf5.R:
~~~~
Rscript msg_read_hdf5.R -i /input/path/inputfile.h5 -p /output/path/
~~~~
Command line arguments:

--file_in, -i: input hdf5 file with its complete path

--path_out, -p: output path (will be created if it does not exist already)
  

![alt tag](https://github.com/UniToDSTGruppoClima/CondMerg/blob/master/logo/CondMerg.png) |
# CondMerg
CondMerg is an open source software developed in R language, cross platform and easily adaptable to different needing; it implements the conditional merging method, in addition to some other experimental variants of that approach (called bias field conditional merging and mean conditional merging). It is optimized for batch processing of multiple files.

#### PRE-REQUISITE
A Java Virtual Machine is required (https://www.java.com/it/download/).
In addition, following R packages have to be installed: *"rgdal", "xlsx", "automap", "raster", "plyr", "mblm", "doParallel"*.
The easiest way to do this is inside R environment, using the command `install.packages(c("rgdal","xlsx","automap","raster","plyr","mblm","doParallel"))`.
On some Linux systems (e.g.: Ubuntu 14.04 / Mint 17.3), this doesn't work and requires some dependencies that could be installed with shell command `apt-get install default-jdk r-base r-base-dev libgdal1-dev libproj-dev r-cran-plyr` followed by the previous install.packages command in R (without "plyr").

This code has been tested under R version 3.0.2, 3.1.3 and 3.2.4; for older or newer versions problems might occur.

#### INPUT
The program takes in input a folder and processes all CSV files contained in that folder that starts with "RainGauges_" and are followed by a string X defining a date and/or an hour (or everything else). Every single csv has 4 columns defining the weather station identification code (called cod_station in input file), its UTM x and y coordinates (in columns utm_x and utm_y) and the value of the measured precipitation (column precip). In addition, the program needs a geo-referred TIF file for each csv which contains radar data for the same time period used in the csv files. The tif file name has to  start with the suffix "Radar_" and has to continue with the same string X defining a date and/or an hour used in the referred csv file. It has to use the coordinates reference system based on WGS 84 / UTM Zone 32N (but it is  very easy to modify the code to adopt other CRS).
Summarizing, input files are:
  - RainGauges_X.csv -> Tabular information given by rain gauges
  - Radar_X.tif -> Raster information given by radar

#### OUTPUT
The program creates 3 folders where it saves computation results related to the merging between radar data and gauges data. Main results are geo-referred TIF files about the merging processes (ordinary kriging of rain gauges, conditional merging, bias field conditional merging and mean conditional merging) and about differences between raw radar data and single merging methods. Other important output are related to the cross-validation and are mainly tables about statistical information (indexes) and scatter plots about data compared with the rain gauges. There are also scatter plots about the cross-validated data compared with the not validated one. Remember that the X in the filenames is the string defining a date and/or an hour (or everything else) that is used in input files.
Summarizing, we have:
- Folder 1 - Main output images
  -   Krig_gauges_X.geo.tif -> Raster containing kriging of rain gauges
  -   Cond_merg_X.geo.tif -> Raster containing conditional merging
  -   Bias_field_cond_merg_X.geo.tif -> Raster containing bias field cond. merg.
  -   Mean_cond_merg_X.geo.tif -> Raster containing mean conditional merging
  -   Data_X_k-fold.csv -> Temp information used for k-fold cross validation
- Folder 2 - Differences between raw data and previous images
  -   Err_radar_menus_condmerg_X.geo.tif -> Raster of radar - cond. merg.
  -   Err_radar_menus_bfcondmerg_X.geo.tif -> Radar - b. f. cond. merg.
  -   Err_radar_menus_meancondmerg_X.geo.tif -> Radar - mean cond. merg.
- Folder 3 - Statistical cross-validation using k-fold method
  -   k-fold_condmerg_errors_X.png -> Validated data vs not validated
  -   k-fold_bfcondmerg_errors_X.png -> Validated data vs not validated
  -   k-fold_meancondmerg_errors_X.png -> Validated data vs not validated
  -  Indexes_k-fold_X.csv -> Detailed indexes on previous information
  -   scatter_raw_radar_X.png -> Scatter plot (SP) of raw radar vs gauges
  -   scatter_condmerg_k-fold_X.png -> SP of validated condmerg vs gauges
  -   scatter_condmerg_X.png -> SP of conditional merging vs gauges
  -  scatter_bfcondmerg_k-fold_X.png -> SP of validated bfcondmerg vs gauges
  -  scatter_bfcondmerg_X.png -> SP of b. f. conditional merging vs gauges
  -  scatter_meancondmerg_k-fold_X.png -> SP of validated meancond vs gauges
  -  scatter_meancondmerg_X.png -> SP of mean conditional merging vs gauges
  -  Indexes_vs_pluvio_X.xlsx -> Detailed indexes on previous information

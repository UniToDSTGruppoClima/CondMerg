###############################################################################
############################################################################### 
########################        CondMerg         ##############################
############  Conditional merging and k-fold cross validation  ################
###############################################################################
###############################################################################
###                                                                         ###
### Copyright (C) 2016 - Diego Guenzi                                       ###
###                                                                         ###
###                                                                         ###
### This program is free software: you can redistribute it and/or modify    ###
### it under the terms of the GNU General Public License as published by    ###
### the Free Software Foundation, either version 3 of the License, or       ###
### (at your option) any later version.                                     ###
###                                                                         ###
### This program is distributed in the hope that it will be useful,         ###
### but WITHOUT ANY WARRANTY; without even the implied warranty of          ###
### MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the           ###
### GNU General Public License for more details.                            ###
###                                                                         ###
### You should have received a copy of the GNU General Public License       ###
### along with this program.  If not, see <http://www.gnu.org/licenses/>    ###
###                                                                         ###
###############################################################################
###############################################################################


############################### CondMerg ######################################
# DESCRIPTION OF THE SOFTWARE
# CondMerg is an open source software developed in R language, cross platform 
# and easily adaptable to different needing; it implements the conditional 
# merging method, in addition to some other experimental variants of that 
# approach (called bias field conditional merging and mean conditional
# merging). It is optimized for batch processing of multiple files.
#
# PRE-REQUISITE
# A Java Virtual Machine is required (https://www.java.com/it/download/).
# In addition, following R-packages have to be installed:
#   "rgdal", "xlsx", "automap", "raster", "plyr", "mblm", "doParallel"
# The easiest way to do this is inside R environment, using the command
# install.packages(c("rgdal","xlsx","automap","raster","plyr","mblm","doParallel"))
# On some Linux systems (e.g.: Ubuntu 14.04 / Mint 17.3), this doesn't work
# and requires some dependencies that could be installed with shell command
# apt-get install default-jdk r-base r-base-dev libgdal1-dev libproj-dev r-cran-plyr
# followed by the previous install.packages command in R (without "plyr").
#
# INPUT
# The program takes in input a folder and processes all CSV files contained in
# that folder that starts with "RainGauges_" and are followed by a string X
# defining a date and/or an hour (or everything else). Every single csv has 4 
# columns defining the weather station identification code (called cod_station
# in input file), its UTM x and y coordinates (in columns utm_x and utm_y) and
# the value of the measured precipitation (column precip). In addition, the
# program needs a geo-referred TIF file for each csv which contains radar data
# for the same time period used in the csv files. The tif file name has to 
# start with the suffix "Radar_" and has to continue with the same string X
# defining a date and/or an hour used in the referred csv file. It has to use
# the coordinates reference system based on WGS 84 / UTM Zone 32N (but it is 
# very easy to modify the code to adopt other CRS).
# Summarizing, input files are:
#   RainGauges_X.csv -> Tabular information given by rain gauges
#   Radar_X.tif -> Raster information given by radar
#
# OUTPUT
# The program creates 3 folders where it saves computation results related to
# the merging between radar data and gauges data. Main results are 
# geo-referred TIF files about the merging processes (ordinary kriging of
# rain gauges, conditional merging, bias field conditional merging and mean
# conditional merging) and about differences between raw radar data and 
# single merging methods. Other important output are related to the 
# cross-validation and are mainly tables about statistical information
# (indexes) and scatter plots about data compared with the rain gauges.
# There are also scatter plots about the cross-validated data compared
# with the not validated one.
# Remember that the X in the filenames is the string defining a date 
# and/or an hour (or everything else) that is used in input files.
# Summarizing, we have:
# Folder 1 - Main output images
#   Krig_gauges_X.geo.tif -> Raster containing kriging of rain gauges
#   Cond_merg_X.geo.tif -> Raster containing conditional merging
#   Bias_field_cond_merg_X.geo.tif -> Raster containing bias field cond. merg.
#   Mean_cond_merg_X.geo.tif -> Raster containing mean conditional merging
#   Data_X_k-fold.csv -> Temp information used for k-fold cross validation
# Folder 2 - Differences between raw data and previous images
#   Err_radar_menus_condmerg_X.geo.tif -> Raster of radar - cond. merg.
#   Err_radar_menus_bfcondmerg_X.geo.tif -> Radar - b. f. cond. merg.
#   Err_radar_menus_meancondmerg_X.geo.tif -> Radar - mean cond. merg.
# Folder 3 - Statistical cross-validation using k-fold method
#   k-fold_condmerg_errors_X.png -> Validated data vs not validated
#   k-fold_bfcondmerg_errors_X.png -> Validated data vs not validated
#   k-fold_meancondmerg_errors_X.png -> Validated data vs not validated
#   Indexes_k-fold_X.csv -> Detailed indexes on previous information
#   scatter_raw_radar_X.png -> Scatter plot (SP) of raw radar vs gauges
#   scatter_condmerg_k-fold_X.png -> SP of validated condmerg vs gauges
#   scatter_condmerg_X.png -> SP of conditional merging vs gauges
#   scatter_bfcondmerg_k-fold_X.png -> SP of validated bfcondmerg vs gauges
#   scatter_bfcondmerg_X.png -> SP of b. f. conditional merging vs gauges
#   scatter_meancondmerg_k-fold_X.png -> SP of validated meancond vs gauges
#   scatter_meancondmerg_X.png -> SP of mean conditional merging vs gauges
#   Indexes_vs_pluvio_X.xlsx -> Detailed indexes on previous information
#
# This code has been tested under R-Version 3.0.2, 3.1.3 and 3.2.4; for older
# or newer versions problems might occur.
###############################################################################
# Versions:
# 
# v1.0 - 20160331: First public release after code cleaning and review
# 
###############################################################################


###############################################################################
####                         INITIALIZATION                                ####
###############################################################################

#### Hardcoded and editable parameters ####

if (interactive()) { # User has to choose the folder for input and output
  setwd(choose.dir(getwd(), "Choose the folder with your input files"))
} else { # Folder for input and output are hardcoded
  setwd("C:/temp/")
}
ignore_small_precipitation_in_bias = TRUE # If TRUE ignores denominator rains < 1mm/h for bias analysis; if FALSE, bias could grow a lot (i.e.: 1/0.001)
k = 2 # Iterations of k-fold cross validation on training sets and test sets (k=0 means do not execute k-fold)
percentile = "99%" # Definition of the percentile of data used for second regression line in plots (just for them, not used into elaborations)
results_folder = "conditional_merging_results" # Name of the first output folder
error_folder = "differences_results" # Name of the second output folder
kfold_folder = "k-fold_cross_validation_results" # Name of the third output folder

#### End of hardcoded and editable parameters ####

start.time = Sys.time()
library(xlsx)
library(automap)
library(raster)
library(plyr)
library(mblm)
library(doParallel)
if (k==1) { # Doesn't make sense
  k = 0 # Use it like "do not execute k-fold"
}
dir.create(paste("./",results_folder,"/",sep=""), showWarnings=FALSE)
dir.create(paste("./",error_folder,"/",sep=""), showWarnings=FALSE)
if (k!=0) {
  dir.create(paste("./",kfold_folder,"/",sep=""), showWarnings=FALSE)
}
runOnCores = detectCores()
if (runOnCores > 2) {
  runOnCores = runOnCores -1
}
mycluster = makeCluster(runOnCores)
registerDoParallel(mycluster, cores=runOnCores)

###############################################################################
####                       FUNCTION DEFINITION                             ####
###############################################################################

OrdinaryKriging <- function(f_formula, f_data, f_grid, f_filename) {
  ##########################################################
  # This function does the ordinary kriging
  #
  # INPUT
  # f_formula: has to be in the format "<name of the field in f_data you want to krige>~1"
  # f_data: points you want to spatialize using f_formula
  # f_grid: a grid used to define size and resolution of kriging
  # f_filename: if has a value, it's the geo tif filename (that could be used in QGIS) where kriging will be saved
  #
  # OUTPUT
  # Raster image containing the kriging
  ##########################################################
  
  kriging_result = autoKrige(f_formula, f_data, f_grid, fix.values=c(0,NA,NA))
  kriging_result$krige_output$var1.pred[kriging_result$krige_output$var1.pred < 0] = 0 # Removes negative values
  kriging_result$krige_output$var1.pred = round(kriging_result$krige_output$var1.pred, digits = 4)
  rain = raster(kriging_result$krige_output)
  crs(rain) = "+proj=utm +zone=32 +datum=WGS84 +units=m +no_defs" # WGS 84 / UTM Zone 32N
  if (!missing(f_filename)) {
    writeRaster(rain, filename=f_filename, format="GTiff", overwrite=TRUE, datatype="FLT4S", options=c("COMPRESS=NONE"))
  }
  
  return(rain)
} # end of function OrdinaryKriging

ScatterPlot <- function(title, label_x, my_x, label_y, my_y, my_filename) {
  ##########################################################
  # This function makes a scatter plot with two regression lines (whole dataset and < of a specific percentile)
  # Warning: it uses a globally set variable called "percentile"
  #
  # INPUT
  # title: name of the plot
  # label_x: x axis name
  # my_x: data for the x axis
  # label_y: y axis name
  # my_y: data for the y axis
  # my_filename: png plot filename
  #
  # OUTPUT
  # Dataframe containing the R squared, a and b coefficients for the regression line
  ##########################################################
  
  my_max = max(my_x, my_y, na.rm=TRUE) # Find maximum value in dataset
  my_max = ((my_max %/% 5) + 1) * 5 # Find the multiple of 5 just after my_max to plot both axis
  
  png(my_filename, width=7, height=7, units="in", res=120)
  plot(my_x, my_y, xlab=label_x, ylab=label_y, xlim=c(0,my_max), ylim=c(0,my_max), main=title, col="red") # Plots data
  
  lin_regression = lm(formula=my_y~my_x) # Computes linear regression on the whole dataset
  abline(lin_regression, col="red") # Plots linear regression
  plot_info = summary(lin_regression)
  r2_line = round(plot_info$r.squared, 5) # R squared
  par_a_line = round(plot_info$coefficients[2], 5) # Parameter a in equation y=ax+b
  par_b_line = round(plot_info$coefficients[1], 5) # Parameter b in equation y=ax+b
  text(par("usr")[2], par("usr")[3]+(2*my_max)/30, paste("R^2 = ",r2_line," ",sep=""), adj=c(1,0), col="red", cex=.75) # Plot R squared
  text(par("usr")[2], par("usr")[3]+(1*my_max)/30, paste("y = ",par_a_line," x + ",par_b_line," ",sep=""), adj=c(1,0), col="red", cex=.75) # Plot equation
  
  theilsen = mblm(formula=my_y~my_x) # Fits data with Theil-Sen single median method
  abline(theilsen, col="orange") # Plots Theil-Sen estimator
  
  quant_val = quantile(abs(my_x - my_y), probs=seq(0, 1, 0.01)) # Quantiles definition on differences
  data_up_to_quant = cbind(my_x, my_y, abs(my_x - my_y))
  data_up_to_quant[data_up_to_quant[,3] > quant_val[percentile]] = NA # Removes data greater than previously defined percentile
  points(data_up_to_quant[,1], data_up_to_quant[,2], col="black") # Overwrites all points < defined percentile with new color
  
  lin_regr_except_percentile = lm(formula=data_up_to_quant[,2]~data_up_to_quant[,1])
  abline(lin_regr_except_percentile, col="black")
  plot_info = summary(lin_regr_except_percentile)
  r2_percentile = round(plot_info$r.squared, 5)
  par_a_percentile = round(plot_info$coefficients[2], 5)
  par_b_percentile = round(plot_info$coefficients[1], 5)
  text(par("usr")[2], par("usr")[3]+(4*my_max)/30, paste("R^2 = ",r2_percentile," ",sep=""), adj=c(1,0), col="black", cex=.75)
  text(par("usr")[2], par("usr")[3]+(3*my_max)/30, paste("y = ",par_a_percentile," x + ",par_b_percentile," ",sep=""), adj=c(1,0), col="black", cex=.75)
  
  abline(0, 1, col="green", lty=2)
  legend('topleft', lty=1, col=c("black","red","orange","green"), legend=c(paste("Regression if >",percentile,"percentile excluded"),
                                                                           "Regression on the whole dataset","Theil-Sen estimator","Bisector"), bty='n', cex=.75)
  dev.off()
  
  return(data.frame(r2_line, par_a_line, par_b_line))
} # end of function ScatterPlot

###############################################################################
####                              MAIN PROGRAM                             ####
###############################################################################

#### Loop for each data file in current directory ####

my_list = list.files(pattern="^RainGauges_")
# SINGLE THREAD VERSION IS for (i in seq_len(length(my_list))) {
blackhole = foreach(i=icount(length(my_list)), .packages=c("automap","raster","plyr")) %dopar% {
  
  #### Obtain date/time from file ####
  
  filename = strsplit(my_list[i], ".csv$")[[1]]
  date_time = strsplit(filename, "RainGauges_")[[1]][2]
  
  #### Read information from csv ####
  
  csv = read.csv(my_list[i], header=TRUE, sep=",")
  my_csv = na.omit(csv)
  coordinates(my_csv) = ~ utm_x+utm_y
  
  #### Associates a group for k-fold cross validation purposes to csv information ####
  
  if (k==0) { # Execute only on the whole dataset, do not do k-fold
    my_csv$group = sample(1:1, nrow(my_csv), replace=TRUE) # Fake group
  } else { # Do k-fold
    my_csv$group = sample(1:k, nrow(my_csv), replace=TRUE) # Assign group information to data
  }
  group_list = 0:k # List of available groups: 0 for execution on complete dataset, 1..k for k-fold
  
  #### Obtain radar data from geo tif ####
  
  finalcsv = cbind(my_csv$cod_station, my_csv$utm_x, my_csv$utm_y, my_csv$precip, my_csv$group)
  colnames(finalcsv) = c("cod_station", "utm_x", "utm_y", "precip", "group")
  finalcsv = data.frame(finalcsv)
  backupcsv = finalcsv
  coordinates(finalcsv) = ~ utm_x+utm_y
  raw_rad_data = raster(paste("Radar_", date_time, ".tif", sep=""))
  crs(raw_rad_data) = "+proj=utm +zone=32 +datum=WGS84 +units=m +no_defs"
  rad_values = raw_rad_data[finalcsv] # Obtain radar values from tif file
  finalcsv = cbind(backupcsv, rad_values)
  coordinates(finalcsv) = ~ utm_x+utm_y
  
  #### Grid computation based on radar geo tif ####
  
  step = res(raw_rad_data)[1] # Cell edge size
  x_min = extent(raw_rad_data)[1] + step/2 # Starting x
  x_max = extent(raw_rad_data)[2] # Ending x
  y_min = extent(raw_rad_data)[3] + step/2 # Starting y
  y_max = extent(raw_rad_data)[4] # Ending y
  gridx = seq(x_min, x_max, by=step)
  gridy = seq(y_min, y_max, by=step) 
  grid = expand.grid(x=gridx, y=gridy)
  gridded(grid) = ~ x+y # Grid definition
  
  #### Loop for each group in the k-fold ####
  
  for (j in group_list) {
    
    if (j==0) { # Work on complete dataset
      training_set = finalcsv # Complete dataset
    } else { # Work on a subset (using a training set and a test set)
      groups_4_training = subset(group_list, !(group_list %in% c(0,j))) # Groups to be taken in training set (exclude 0 and current j, used for test set)
      training_set = subset(finalcsv, group %in% groups_4_training) # Training set is a subset of the complete dataset
    }
    
    #### Radar Kriging ####
    
    k_radar = OrdinaryKriging(rad_values~1, training_set, grid)
    
    #### Rain gauge Kriging ####
    
    if (j==0) {
      k_gauges = OrdinaryKriging(precip~1, training_set, grid, paste("./",results_folder,"/","krig_gauges_",date_time,".geo.tif",sep=""))
    } else {
      k_gauges = OrdinaryKriging(precip~1, training_set, grid)
    }
    
    #### Conditional merging ####
    
    cond_merg = raw_rad_data - k_radar + k_gauges
    cond_merg[cond_merg < 0] = 0 # Removes negative values
    if (j==0) {
      writeRaster(cond_merg, filename=paste("./",results_folder,"/","cond_merg_",date_time,".geo.tif",sep=""), 
                  format="GTiff", overwrite=TRUE, datatype="FLT4S", options=c("COMPRESS=NONE")) # Saves the geo tif that could be used in QGIS
    }
    
    #### EXPERIMENTAL (could be ignored) Bias field conditional merging ####
    
    bf_cond_merg = (raw_rad_data / k_radar) * k_gauges # Instead of using subtraction, bf_cond_merg uses the ratio (bias)
    bf_cond_merg[is.infinite(bf_cond_merg)] = 0 # Removes infinite values
    bf_cond_merg[is.na(bf_cond_merg)] = 0 # Removes NA values
    if (j==0) {
      writeRaster(bf_cond_merg, filename=paste("./",results_folder,"/","bias_field_cond_merg_",date_time,".geo.tif",sep=""), 
                  format="GTiff", overwrite=TRUE, datatype="FLT4S", options=c("COMPRESS=NONE")) # Saves the geo tif that could be used in QGIS
    }
    
    #### EXPERIMENTAL (could be ignored) Mean conditional merging ####
    
    mean_cond_merg = (bf_cond_merg + cond_merg) /2 # Mean between bias field conditional merging and conditional merging
    mean_cond_merg[mean_cond_merg>(2*cond_merg)] = cond_merg[mean_cond_merg>(2*cond_merg)] # If more than double, ignore it and use cond. merging
    if (j==0) {
      writeRaster(mean_cond_merg, filename=paste("./",results_folder,"/","mean_cond_merg_",date_time,".geo.tif",sep=""), 
                  format="GTiff", overwrite=TRUE, datatype="FLT4S", options=c("COMPRESS=NONE")) # Saves the geo tif that could be used in QGIS
    }
    
    #### Radar over/underestimation (Radar error) ####
    
    if (j==0) {
      ou_estim_condmerg = raw_rad_data - cond_merg
      ou_estim_bfcondmerg = raw_rad_data - bf_cond_merg
      ou_estim_meancondmerg = raw_rad_data - mean_cond_merg
      writeRaster(ou_estim_condmerg, filename=paste("./",error_folder,"/","err_radar_menus_condmerg_",date_time,".geo.tif",sep=""), 
                  format="GTiff", overwrite=TRUE, datatype="FLT4S", options=c("COMPRESS=NONE")) # Saves the geo tif that could be used in QGIS
      writeRaster(ou_estim_bfcondmerg, filename=paste("./",error_folder,"/","err_radar_menus_bfcondmerg_",date_time,".geo.tif",sep=""), 
                  format="GTiff", overwrite=TRUE, datatype="FLT4S", options=c("COMPRESS=NONE")) # Saves the geo tif that could be used in QGIS
      writeRaster(ou_estim_meancondmerg, filename=paste("./",error_folder,"/","err_radar_menus_meancondmerg_",date_time,".geo.tif",sep=""), 
                  format="GTiff", overwrite=TRUE, datatype="FLT4S", options=c("COMPRESS=NONE")) # Saves the geo tif that could be used in QGIS
    }
    
    #### Complete list of values from geo tif ####
    
    radar_coord = xyFromCell(raw_rad_data, 1:ncell(raw_rad_data)) # Extract coordinates from tif
    radar_data = extract(raw_rad_data, 1:ncell(raw_rad_data))
    condmerg_data = extract(cond_merg, 1:ncell(cond_merg))
    bfcondmerg_data = extract(bf_cond_merg, 1:ncell(bf_cond_merg))
    meancondmerg_data = extract(mean_cond_merg, 1:ncell(mean_cond_merg))
    combine = cbind(radar_coord, radar_data, condmerg_data, bfcondmerg_data, meancondmerg_data) # Combine all values
    combine=data.frame(combine)
    
    #### Values for the station points ####
    
    csv_align = finalcsv
    csv_align$X = (finalcsv$utm_x %/% step) * step + (x_min %% step) # Transforms x of stations on x for the grid
    csv_align$Y = (finalcsv$utm_y %/% step) * step + (y_min %% step) # Transforms y of stations on y for the grid
    tab = merge(rename(combine,c("x"="X","y"="Y")), csv_align, by=c("X","Y")) #Join between stations and grid
    tab = rename(tab, c("X"="grid_x","Y"="grid_y","radar_data"="grid_radar_data","rad_values"="station_radar_data",
                        "utm_x"="station_x","utm_y"="station_y","precip"="gauges_data"))
    if (j==0) { # Complete dataset
      ref_tab = tab # Backup of results table
      final_tab = data.frame() # Create empty container for results
    } else { # I'm on a partial training set -> let's start with real k-fold
      tab = subset(tab, (tab$group %in% c(j))) # Use the test set only
      temp = tab[,c(3:7)]
      tab = cbind(temp, station_radar_data=tab[,c(12)])
      tab = rename(tab, c("grid_radar_data"="grid_radar_data_TEST","condmerg_data"="condmerg_data_TEST","bfcondmerg_data"="bfcondmerg_data_TEST",
                          "meancondmerg_data"="meancondmerg_data_TEST","station_radar_data"="station_radar_data_TEST"))
      temp_tab = merge(ref_tab, tab, by=c("cod_station")) # Join between previous values and test set
      drop_cols = c("grid_radar_data_TEST", "station_radar_data_TEST")
      temp_tab = temp_tab[,!(names(temp_tab) %in% drop_cols)] # Drop useless columns
      final_tab = rbind(final_tab, temp_tab) # Prepares final result
      if (j==k) { # Last test set completed
        final_tab = final_tab[,!(names(final_tab) %in% drop_cols)] # Drop useless columns
        final_tab = cbind(final_tab[,1:3], final_tab[,8:9], final_tab[11], final_tab[10], final_tab[4], final_tab[12],
                          final_tab[5], final_tab[13], final_tab[6], final_tab[14], final_tab[7], final_tab[15])
        write.csv(final_tab, file=paste("./",results_folder,"/","Data_",date_time,"_k-fold.csv",sep=""),row.names=FALSE)
      }
    }
  }
}

#### Loop for each file with values from the stations ####

my_list = list.files(path=paste("./",results_folder,"/",sep=""),pattern="^Data_")
# SINGLE THREAD VERSION IS for (i in seq_len(length(my_list))) {
blackhole = foreach(i=icount(length(my_list)), .packages=c("xlsx","mblm")) %dopar% {
  csv = read.csv(paste("./",results_folder,"/",my_list[i],sep=""), header=TRUE, sep=",")
  is.na(csv) = sapply(csv, is.infinite) # Puts NA instead of infinite (obtained from x/0); NAs (from 0/0) will be ignored
  filename = strsplit(my_list[i],"_k-fold.csv$")[[1]]
  date_time = strsplit(filename,"Data_")[[1]][2]
  
  #### Prepares data for error and bias analysis (k-fold groups vs complete dataset) ####
  
  temp_res = cbind(csv[1], csv[7], csv[9]) # Basic info
  temp_res = cbind(temp_res, round(csv[8]-csv[9],5)) # Grid / radar displacement error
  temp_res = cbind(temp_res, round(csv[11]-csv[10],5)) # Conditional merging error
  temp_res = cbind(temp_res, round(csv[10]/csv[11],5)) # Conditional merging bias
  temp_res = cbind(temp_res, round(csv[13]-csv[12],5)) # Bias field conditional merging error
  temp_res = cbind(temp_res, round(csv[12]/csv[13],5)) # Bias field conditional merging bias
  temp_res = cbind(temp_res, round(csv[15]-csv[14],5)) # Mean conditional merging error
  temp_res = cbind(temp_res, round(csv[14]/csv[15],5)) # Mean conditional merging bias
  temp_res = rename(temp_res, c("grid_radar_data"="err_displac_grid","condmerg_data_TEST"="err_condmerg",
                                "condmerg_data"="bias_condmerg","bfcondmerg_data_TEST"="err_bfcondmerg","bfcondmerg_data"="bias_bfcondmerg",
                                "meancondmerg_data_TEST"="err_meancondmerg","meancondmerg_data"="bias_meancondmerg"))
  results = data.frame(Algorithm=c("Conditional merging","Bias field conditional merging","Mean conditional merging")) # Prepares results dataframe
  
  #### ME - Mean Error (k-fold groups vs complete dataset) ####
  
  condmerg = mean(temp_res$err_condmerg, na.rm=TRUE)
  bfcondmerg = mean(temp_res$err_bfcondmerg, na.rm=TRUE)
  meancondmerg = mean(temp_res$err_meancondmerg, na.rm=TRUE)
  results = cbind(results, Mean_error=c(condmerg,bfcondmerg,meancondmerg))
  
  #### MAE - Mean Absolute Error (k-fold groups vs complete dataset) ####
  
  condmerg = mean(abs(temp_res$err_condmerg), na.rm=TRUE)
  bfcondmerg = mean(abs(temp_res$err_bfcondmerg), na.rm=TRUE)
  meancondmerg = mean(abs(temp_res$err_meancondmerg), na.rm=TRUE)
  results = cbind(results, Mean_absolute_error=c(condmerg,bfcondmerg,meancondmerg))
  
  #### RMAE - Relative Mean Absolute Error (k-fold groups vs complete dataset) ####
  
  condmerg = results$Mean_absolute_error[1]/mean(csv[,10], na.rm=TRUE)
  bfcondmerg = results$Mean_absolute_error[2]/mean(csv[,12], na.rm=TRUE)
  meancondmerg = results$Mean_absolute_error[3]/mean(csv[,14], na.rm=TRUE)
  results = cbind(results, Relative_mean_absolute_error=c(condmerg,bfcondmerg,meancondmerg))
  
  #### MB - Mean Bias (k-fold groups vs complete dataset) ####
  
  is.na(temp_res$bias_condmerg) = sapply(temp_res$bias_condmerg, is.infinite) # Puts NA instead of infinite (from x/0); NAs (from 0/0) will be ignored
  is.na(temp_res$bias_bfcondmerg) = sapply(temp_res$bias_bfcondmerg, is.infinite)
  is.na(temp_res$bias_meancondmerg) = sapply(temp_res$bias_meancondmerg, is.infinite)
  if (ignore_small_precipitation_in_bias) {
    temp_res$bias_condmerg[csv[11] < 1] = NA # Puts NA where denominator < 1 mm (very small rainfall amounts could cause bias to explode)
    temp_res$bias_bfcondmerg[csv[13] < 1] = NA
    temp_res$bias_meancondmerg[csv[15] < 1] = NA
  }
  condmerg = mean(temp_res$bias_condmerg, na.rm=TRUE)
  bfcondmerg = mean(temp_res$bias_bfcondmerg, na.rm=TRUE)
  meancondmerg = mean(temp_res$bias_meancondmerg, na.rm=TRUE)
  results = cbind(results, Mean_bias=c(condmerg,bfcondmerg,meancondmerg))
  
  #### RMSE - Root Mean Square Error (k-fold groups vs complete dataset) ####
  
  condmerg = sqrt(mean((temp_res$err_condmerg)^2, na.rm=TRUE))
  bfcondmerg = sqrt(mean((temp_res$err_bfcondmerg)^2, na.rm=TRUE))
  meancondmerg = sqrt(mean((temp_res$err_meancondmerg)^2, na.rm=TRUE))
  results = cbind(results, Root_mean_square_error=c(condmerg,bfcondmerg,meancondmerg))
  
  #### SD - Standard Deviation (k-fold groups vs complete dataset) ####
  
  condmerg = sqrt(mean((results$Mean_error[1]-csv[,10])^2, na.rm=TRUE))
  bfcondmerg = sqrt(mean((results$Mean_error[2]-csv[,12])^2, na.rm=TRUE))
  meancondmerg = sqrt(mean((results$Mean_error[3]-csv[,14])^2, na.rm=TRUE))
  results = cbind(results, Standard_deviation=c(condmerg,bfcondmerg,meancondmerg))
  
  #### Conditional merging scatter plot (k-fold groups vs complete dataset) ####
  
  stat_condmerg = ScatterPlot("Conditional merging - k-fold groups vs complete dataset",
                              "Complete method", csv[,10], "Estimated by k-fold", csv[,11],
                              paste("./",kfold_folder,"/","k-fold_condmerg_errors_",date_time,".png",sep=""))
  
  #### Bias field conditional merging scatter plot (k-fold groups vs complete dataset) ####
  
  stat_bfcondmerg = ScatterPlot("Bias field conditional merging - k-fold groups vs complete dataset",
                                "Complete method", csv[,12], "Estimated by k-fold", csv[,13],
                                paste("./",kfold_folder,"/","k-fold_bfcondmerg_errors_",date_time,".png",sep=""))
  
  #### Mean conditional merging scatter plot (k-fold groups vs complete dataset) ####
  
  stat_meancondmerg = ScatterPlot("Mean conditional merging - k-fold groups vs complete dataset",
                                  "Complete method", csv[,14], "Estimated by k-fold", csv[,15],
                                  paste("./",kfold_folder,"/","k-fold_meancondmerg_errors_",date_time,".png",sep=""))
  
  #### Write results (k-fold groups vs complete dataset) ####
  
  results = cbind(results, R2=c(stat_condmerg[,1],stat_bfcondmerg[,1],stat_meancondmerg[,1]))
  results = cbind(results, A=c(stat_condmerg[,2],stat_bfcondmerg[,2],stat_meancondmerg[,2]))
  results = cbind(results, B=c(stat_condmerg[,3],stat_bfcondmerg[,3],stat_meancondmerg[,3]))
  write.csv(results,file=paste("./",kfold_folder,"/","Indexes_k-fold_",date_time,".csv",sep=""),row.names=FALSE)
  
  #### Prepares data for error and bias analysis (raw/computed data vs rain gauges) ####
  
  temp_res = cbind(csv[1], csv[7], csv[9]) # Basic info
  temp_res = rename(temp_res,c("gauges_data"="rgauge_data","station_radar_data"="rawdata_radar"))
  temp_res = cbind(temp_res,round(csv[9]-csv[7],5)) # Radar raw data error
  temp_res = rename(temp_res,c("station_radar_data"="err_rawradar"))
  temp_res = cbind(temp_res,round(csv[7]/csv[9],5)) # Radar raw data bias
  temp_res = rename(temp_res,c("gauges_data"="bias_rawradar"))
  temp_res = cbind(temp_res,round(csv[11]-csv[7],5)) # Conditional merging error (k-fold related)
  temp_res = rename(temp_res,c("condmerg_data_TEST"="err_condmerg_kfold"))
  temp_res = cbind(temp_res,round(csv[10]-csv[7],5)) # Conditional merging error
  temp_res = rename(temp_res,c("condmerg_data"="err_condmerg"))
  temp_res = cbind(temp_res,round(csv[7]/csv[11],5)) # Conditional merging bias (k-fold related)
  temp_res = rename(temp_res,c("gauges_data"="bias_condmerg_kfold"))
  temp_res = cbind(temp_res,round(csv[7]/csv[10],5)) # Conditional merging bias
  temp_res = rename(temp_res,c("gauges_data"="bias_condmerg"))
  temp_res = cbind(temp_res,round(csv[13]-csv[7],5)) # Bias field conditional merging error (k-fold related)
  temp_res = rename(temp_res,c("bfcondmerg_data_TEST"="err_bfcondmerg_kfold"))
  temp_res = cbind(temp_res,round(csv[12]-csv[7],5)) # Bias field conditional merging error
  temp_res = rename(temp_res,c("bfcondmerg_data"="err_bfcondmerg"))
  temp_res = cbind(temp_res,round(csv[7]/csv[13],5)) # Bias field conditional merging bias (k-fold related)
  temp_res = rename(temp_res,c("gauges_data"="bias_bfcondmerg_kfold")) 
  temp_res = cbind(temp_res,round(csv[7]/csv[12],5)) # Bias field conditional merging bias
  temp_res = rename(temp_res,c("gauges_data"="bias_bfcondmerg"))
  temp_res = cbind(temp_res,round(csv[15]-csv[7],5)) # Mean conditional merging error (k-fold related)
  temp_res = rename(temp_res,c("meancondmerg_data_TEST"="err_meancondmerg_kfold"))
  temp_res = cbind(temp_res,round(csv[14]-csv[7],5)) # Mean conditional merging error
  temp_res = rename(temp_res,c("meancondmerg_data"="err_meancondmerg"))
  temp_res = cbind(temp_res,round(csv[7]/csv[15],5)) # Mean conditional merging bias (k-fold related)
  temp_res = rename(temp_res,c("gauges_data"="bias_meancondmerg_kfold"))
  temp_res = cbind(temp_res,round(csv[7]/csv[14],5)) # Mean conditional merging bias
  temp_res = rename(temp_res,c("gauges_data"="bias_meancondmerg"))
  temp_res = rename(temp_res,c("rgauge_data"="gauges_data"))
  results = data.frame(Algorithm=c("Raw radar data","Conditional merging (k-fold)","Conditional merging","Bias field conditional merging (k-fold)",
                                   "Bias field conditional merging","Mean conditional merging (k-fold)","Mean conditional merging"))
  
  #### ME - Mean Error (raw/computed data vs rain gauges) ####
  
  rawradar = mean(temp_res$err_rawradar, na.rm=TRUE)
  condmerg_kfold = mean(temp_res$err_condmerg_kfold, na.rm=TRUE)
  condmerg = mean(temp_res$err_condmerg, na.rm=TRUE)
  bfcondmerg_kfold = mean(temp_res$err_bfcondmerg_kfold, na.rm=TRUE)
  bfcondmerg = mean(temp_res$err_bfcondmerg, na.rm=TRUE)
  meancondmerg_kfold = mean(temp_res$err_meancondmerg_kfold, na.rm=TRUE)
  meancondmerg = mean(temp_res$err_meancondmerg, na.rm=TRUE)
  results = cbind(results, Mean_error=c(rawradar,condmerg_kfold,condmerg,bfcondmerg_kfold,bfcondmerg,meancondmerg_kfold,meancondmerg))
  
  #### MAE - Mean Absolute Error (raw/computed data vs rain gauges) ####
  
  rawradar = mean(abs(temp_res$err_rawradar), na.rm=TRUE)
  condmerg_kfold = mean(abs(temp_res$err_condmerg_kfold), na.rm=TRUE)
  condmerg = mean(abs(temp_res$err_condmerg), na.rm=TRUE)
  bfcondmerg_kfold = mean(abs(temp_res$err_bfcondmerg_kfold), na.rm=TRUE)
  bfcondmerg = mean(abs(temp_res$err_bfcondmerg), na.rm=TRUE)
  meancondmerg_kfold = mean(abs(temp_res$err_meancondmerg_kfold), na.rm=TRUE)
  meancondmerg = mean(abs(temp_res$err_meancondmerg), na.rm=TRUE)
  results = cbind(results, Mean_absolute_error=c(rawradar,condmerg_kfold,condmerg,bfcondmerg_kfold,bfcondmerg,meancondmerg_kfold,meancondmerg))
  
  #### RMAE - Relative Mean Absolute Error (raw/computed data vs rain gauges) ####
  
  gauges_mean = mean(csv[,7], na.rm=TRUE)
  rawradar = results$Mean_absolute_error[1]/gauges_mean
  condmerg_kfold = results$Mean_absolute_error[2]/gauges_mean
  condmerg = results$Mean_absolute_error[3]/gauges_mean
  bfcondmerg_kfold = results$Mean_absolute_error[4]/gauges_mean
  bfcondmerg = results$Mean_absolute_error[5]/gauges_mean
  meancondmerg_kfold = results$Mean_absolute_error[6]/gauges_mean
  meancondmerg = results$Mean_absolute_error[7]/gauges_mean
  results = cbind(results, Relative_mean_absolute_error=c(rawradar,condmerg_kfold,condmerg,bfcondmerg_kfold,bfcondmerg,meancondmerg_kfold,meancondmerg))
  
  #### MB - Mean Bias (raw/computed data vs rain gauges) ####
  
  is.na(temp_res$bias_rawradar) = sapply(temp_res$bias_rawradar, is.infinite) # Puts NA instead of infinite (from x/0); NAs (from 0/0) will be ignored
  is.na(temp_res$bias_condmerg_kfold) = sapply(temp_res$bias_condmerg_kfold, is.infinite)
  is.na(temp_res$bias_condmerg) = sapply(temp_res$bias_condmerg, is.infinite)
  is.na(temp_res$bias_bfcondmerg_kfold) = sapply(temp_res$bias_bfcondmerg_kfold, is.infinite)
  is.na(temp_res$bias_bfcondmerg) = sapply(temp_res$bias_bfcondmerg, is.infinite)
  is.na(temp_res$bias_meancondmerg_kfold) = sapply(temp_res$bias_meancondmerg_kfold, is.infinite)
  is.na(temp_res$bias_meancondmerg) = sapply(temp_res$bias_meancondmerg, is.infinite)
  if (ignore_small_precipitation_in_bias) {
    temp_res$bias_rawradar[csv[9] < 1] = NA # Puts NA where denominator < 1 mm (very small rainfall amounts could cause bias to explode)
    temp_res$bias_condmerg_kfold[csv[11] < 1] = NA
    temp_res$bias_condmerg[csv[10] < 1] = NA
    temp_res$bias_bfcondmerg_kfold[csv[13] < 1] = NA
    temp_res$bias_bfcondmerg[csv[12] < 1] = NA
    temp_res$bias_meancondmerg_kfold[csv[15] < 1] = NA
    temp_res$bias_meancondmerg[csv[14] < 1] = NA
  }
  rawradar = mean(temp_res$bias_rawradar, na.rm=TRUE)
  condmerg_kfold = mean(temp_res$bias_condmerg_kfold, na.rm=TRUE)
  condmerg = mean(temp_res$bias_condmerg, na.rm=TRUE)
  bfcondmerg_kfold = mean(temp_res$bias_bfcondmerg_kfold, na.rm=TRUE)
  bfcondmerg = mean(temp_res$bias_bfcondmerg, na.rm=TRUE)
  meancondmerg_kfold = mean(temp_res$bias_meancondmerg_kfold, na.rm=TRUE)
  meancondmerg = mean(temp_res$bias_meancondmerg, na.rm=TRUE)
  results = cbind(results, Mean_bias=c(rawradar,condmerg_kfold,condmerg,bfcondmerg_kfold,bfcondmerg,meancondmerg_kfold,meancondmerg))
  
  #### RMSE - Root Mean Square Error (raw/computed data vs rain gauges) ####
  
  rawradar = sqrt(mean((temp_res$err_rawradar)^2, na.rm=TRUE))
  condmerg_kfold = sqrt(mean((temp_res$err_condmerg_kfold)^2, na.rm=TRUE))
  condmerg = sqrt(mean((temp_res$err_condmerg)^2, na.rm=TRUE))
  bfcondmerg_kfold = sqrt(mean((temp_res$err_bfcondmerg_kfold)^2, na.rm=TRUE))
  bfcondmerg = sqrt(mean((temp_res$err_bfcondmerg)^2, na.rm=TRUE))
  meancondmerg_kfold = sqrt(mean((temp_res$err_meancondmerg_kfold)^2, na.rm=TRUE))
  meancondmerg = sqrt(mean((temp_res$err_meancondmerg)^2, na.rm=TRUE))
  results = cbind(results, Root_mean_square_error=c(rawradar,condmerg_kfold,condmerg,bfcondmerg_kfold,bfcondmerg,meancondmerg_kfold,meancondmerg))
  
  #### SD - Standard Deviation (raw/computed data vs rain gauges) ####
  
  rawradar = sqrt(mean((results$Mean_error[1]-csv[,7])^2, na.rm=TRUE))
  condmerg_kfold = sqrt(mean((results$Mean_error[2]-csv[,7])^2, na.rm=TRUE))
  condmerg = sqrt(mean((results$Mean_error[3]-csv[,7])^2, na.rm=TRUE))
  bfcondmerg_kfold = sqrt(mean((results$Mean_error[4]-csv[,7])^2, na.rm=TRUE))
  bfcondmerg = sqrt(mean((results$Mean_error[5]-csv[,7])^2, na.rm=TRUE))
  meancondmerg_kfold = sqrt(mean((results$Mean_error[6]-csv[,7])^2, na.rm=TRUE))
  meancondmerg = sqrt(mean((results$Mean_error[7]-csv[,7])^2, na.rm=TRUE))
  results = cbind(results, Standard_deviation=c(rawradar,condmerg_kfold,condmerg,bfcondmerg_kfold,bfcondmerg,meancondmerg_kfold,meancondmerg))
  
  #### Raw radar data vs rain gauges scatter plot ####
  
  stat_rawradar = ScatterPlot("Raw radar data vs rain gauges", "Rain gauges", csv[,7], "Raw radar data", csv[,9],
                              paste("./",kfold_folder,"/","scatter_raw_radar_",date_time,".png",sep=""))
  
  #### Conditional merging (k-fold) vs rain gauges scatter plot ####
  
  stat_condmergk = ScatterPlot("Conditional merging (k-fold) vs rain gauges", "Rain gauges", csv[,7], "Conditional merging (k-fold)", csv[,11],
                               paste("./",kfold_folder,"/","scatter_condmerg_k-fold_",date_time,".png",sep=""))
  
  #### Conditional merging vs rain gauges scatter plot ####
  
  stat_condmerg = ScatterPlot("Conditional merging vs rain gauges", "Rain gauges", csv[,7], "Conditional merging", csv[,10],
                              paste("./",kfold_folder,"/","scatter_condmerg_",date_time,".png",sep=""))
  
  #### Bias field conditional merging (k-fold) vs rain gauges scatter plot ####
  
  stat_bfcondmergk = ScatterPlot("Bias field conditional merging (k-fold) vs rain gauges", "Rain gauges", csv[,7], "Bias field conditional merging (k-fold)", csv[,13],
                                 paste("./",kfold_folder,"/","scatter_bfcondmerg_k-fold_",date_time,".png",sep=""))
  
  #### Bias field conditional merging vs rain gauges scatter plot ####
  
  stat_bfcondmerg = ScatterPlot("Bias field conditional merging vs rain gauges", "Rain gauges", csv[,7], "Bias field conditional merging", csv[,12],
                                paste("./",kfold_folder,"/","scatter_bfcondmerg_",date_time,".png",sep=""))
  
  #### Mean conditional merging (k-fold) vs rain gauges scatter plot ####
  
  stat_meancondmergk = ScatterPlot("Mean conditional merging (k-fold) vs rain gauges", "Rain gauges", csv[,7], "Mean conditional merging (k-fold)", csv[,15],
                                   paste("./",kfold_folder,"/","scatter_meancondmerg_k-fold_",date_time,".png",sep=""))
  
  #### Mean conditional merging vs rain gauges scatter plot ####
  
  stat_meancondmerg = ScatterPlot("Mean conditional merging vs rain gauges", "Rain gauges", csv[,7], "Mean conditional merging", csv[,14],
                                  paste("./",kfold_folder,"/","scatter_meancondmerg_",date_time,".png",sep=""))
  
  #### Prepares results for writing (raw/computed data vs rain gauges) ####
  
  results = cbind(results, R2=c(stat_rawradar[,1],stat_condmergk[,1],stat_condmerg[,1],stat_bfcondmergk[,1],
                                stat_bfcondmerg[,1],stat_meancondmergk[,1],stat_meancondmerg[,1]))
  results = cbind(results, A=c(stat_rawradar[,2],stat_condmergk[,2],stat_condmerg[,2],stat_bfcondmergk[,2],
                               stat_bfcondmerg[,2],stat_meancondmergk[,2],stat_meancondmerg[,2]))
  results = cbind(results, B=c(stat_rawradar[,3],stat_condmergk[,3],stat_condmerg[,3],stat_bfcondmergk[,3],
                               stat_bfcondmerg[,3],stat_meancondmergk[,3],stat_meancondmerg[,3]))
  
  #### Write results in xlsx format (raw/computed data vs rain gauges) ####
  
  #write.csv(results,file=paste("./",kfold_folder,"/","Indexes_vs_pluvio_",date_time,".csv",sep=""),row.names=FALSE) # Simple alternative to xlsx (disabled)
  new_results = results
  myexcel = paste("./",kfold_folder,"/","Indexes_vs_pluvio_",date_time,".xlsx",sep="")
  write.xlsx(new_results, myexcel, sheetName="Summary", row.names=FALSE)
  wb = loadWorkbook(myexcel)
  fo_g = Fill(foregroundColor="green")
  cs_g = CellStyle(wb, fill=fo_g)
  fo_y = Fill(foregroundColor="yellow") # Highlight color
  cs_y = CellStyle(wb, fill=fo_y) # Cell style with highlight
  sheets = getSheets(wb)
  sheet = sheets[["Summary"]]
  
  highlight = NULL # Prepares vector for highlight
  rows = getRows(sheet, rowIndex=2:(nrow(new_results)+1)) # Number of rows
  for (w in 2:length(new_results[1,])) { # For each column
    cells = getCells(rows, colIndex = w) # Get all cells in column w
    values = lapply(cells, getCellValue) # Get all values from cells
    if (w==5 || w==8 || w==9) { # On these rows, index is good if near 1
      minimum = values[min(abs(1-as.numeric(values))) == abs(1-as.numeric(values))] # Take the nearest to 1 value
    } else { # On other rows, index is good if near 0
      minimum = values[min(abs(as.numeric(values))) == abs(as.numeric(values))] # Take the nearest to 0 value
    }
    for (j in names(values)) { # For every value in the column
      if (as.numeric(values[j]) == minimum) { # Value is minimal?
        highlight = c(highlight, j) # Write its position in highlight vector
      }
    }
  }
  cells = getCells(rows, colIndex = 2:length(new_results[1,])) # Take all cells
  lapply(names(cells[highlight]),function(ii)setCellStyle(cells[[ii]],cs_y)) # Set style at highlight cells
  
  highlight = NULL # Resets vector for highlight
  for (w in 2:length(new_results[1,])) { # For each column
    cells = getCells(rows, colIndex = w) # Get all cells in column w
    values = lapply(cells, getCellValue) # Get all values from cells
    values[3] = Inf
    values[5] = Inf
    values[7] = Inf # Exclude values from non k-fold rows
    if (w==5 || w==8 || w==9) { # On these rows, index is good if near 1
      minimum = values[min(abs(1-as.numeric(values))) == abs(1-as.numeric(values))] # Take the nearest to 1 value
    } else { # On other rows, index is good if near 0
      minimum = values[min(abs(as.numeric(values))) == abs(as.numeric(values))] # Take the nearest to 0 value
    }
    for (j in names(values)) { # For every value in the column
      if (as.numeric(values[j]) == minimum) { # Value is minimal?
        highlight = c(highlight, j) # Write its position in highlight vector
      }
    }
  }
  cells = getCells(rows, colIndex = 2:length(new_results[1,])) # Take all cells
  lapply(names(cells[highlight]),function(ii)setCellStyle(cells[[ii]],cs_g)) # Set style at highlight cells
  
  saveWorkbook(wb, myexcel) # Update previously created xlsx
}

stopCluster(mycluster)
end.time = Sys.time()
cat(capture.output(end.time - start.time))

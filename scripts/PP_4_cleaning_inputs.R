####################################################################################################
#### Format RADS input data
#### Author: Allie Rhea (Allison.Rhea@colostate.edu)
#### Date Created: 03/18/2021
####################################################################################################
#
#
# Summary: this script uniformly prepares all RADS input data. First, it projects all layers to UTM
# zone 13N and then it clips the extent to the county boundary with a 5 km buffer.
#
#
# Note: County boundary was buffered in ArcGIS pro.The Buffer tool was used with the following 
# specifications - distance = 5 km, Side type = full, method = planar, dissolve type = no dissolve.
#
####################################################################################################
#-> Get working and packages directory paths from command call
setwd('C:/Users/aerhe/Desktop/LC_RADS/scripts')
wd <- getwd()
pd <- paste(c(unlist(strsplit(wd,'/'))[1:(length(unlist(strsplit(wd,'/')))-1)],
              'Portable_R/Packages'),collapse='/')
####################################################################################################

###########################################START MESSAGE############################################
cat('Process fuel treatment data\n',sep='')
cat('Started at: ',as.character(Sys.time()),'\n\n',sep='')
cat('Messages, Errors, and Warnings (if they exist):\n')
####################################################################################################

############################################START SET UP############################################

#-> Load packages
.libPaths(pd)
packages <- c('raster','rgdal','rgeos','st')
for(package in packages){
  if(suppressMessages(!require(package,lib.loc=pd,character.only=T))){
    install.packages(package,lib=pd,repos='https://repo.miserver.it.umich.edu/cran/')
    suppressMessages(library(package,lib.loc=pd,character.only=T))
  }
}

#-> Maximize raster processing speed
rasterOptions(maxmemory=10^9)

#-> Load in buffered AOI extent
setwd(paste0(wd,'/INPUT/SPATIAL/Extent')) #===> Change directory
ma <- shapefile('ma_extent.shp')

#-> Load input data
setwd(paste0(wd,'/INPUT/SPATIAL/Test')) #===> Change directory
rds <- shapefile('CO_USFS_roads.shp')

#-> Relabel projection information so it plays well with Arc
proj <- '+proj=utm +zone=13 +ellps=GRS80 +datum=NAD83'
crs(ma) <- CRS(proj)
crs(rds) <- CRS(proj)
plot(rds)

#->L crop input layers to buffered AOI extent
roads <- crop(roads,cc)

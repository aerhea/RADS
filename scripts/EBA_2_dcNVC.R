####################################################################################################
#### Calculate delta conditional net value change (dcNVC) by treatment
#### Author: Ben Gannon (bengannon@gmail.com)
#### Date Created: 05/07/2019
#### Last Modified: 02/04/2021
####################################################################################################
# This is a generic workflow designed to calculate dcNVC by treatment from an ecological benefits
# table.
####################################################################################################
#-> Get working directory path from command call
#initial.options <- commandArgs(trailingOnly=F)
#setwd(dirname(sub('--file=','',initial.options[grep('--file=',initial.options)])))
setwd('C:/Users/semue/Documents/CFRI/Projects/JCOS/BOX_RADS/Modeling/RAD_RUNS/JCOS_RADS/JCOS_RADS_03052021/RADS_application/scripts')
wd <- getwd()
pd <- paste(c(unlist(strsplit(wd,'/'))[1:(length(unlist(strsplit(wd,'/')))-1)],
            'Portable_R/Packages'),collapse='/')
####################################################################################################

###########################################START MESSAGE############################################
cat('Calculate delta conditional net value change (dcNVC) by treatment\n',sep='')
log <- file(paste0(wd,'/OUTPUT/EBA_dcNVC/EBA2.log'))
sink(file=log,append=T,type=c('output','message'),split=T)
cat('Started at: ',as.character(Sys.time()),'\n\n',sep='')
cat('Messages, Errors, and Warnings (if they exist):\n')
####################################################################################################

############################################START SET UP############################################

#-> Load packages
.libPaths(pd)
packages <- c('raster','rgdal','rgeos','plyr')
for(package in packages){
	if(suppressMessages(!require(package,lib.loc=pd,character.only=T))){
		install.packages(package,lib=pd,repos='https://repo.miserver.it.umich.edu/cran/')
		suppressMessages(library(package,lib.loc=pd,character.only=T))
	}
}

#-> Maximize raster processing speed
rasterOptions(maxmemory = 10^9)

#-> Load in input tables
setwd(paste0(wd,'/INPUT')) #===> Change directory
rs <- read.csv('EBA_settings.csv',header=T)
rs <- rs[rs$Include==1,]
ms <- read.csv('map_settings.csv',header=T)
ri <- read.csv('EBA_relative_importance.csv',header=T)
tspecs <- read.csv('so_treatment_specs.csv',header=T)

#-> Load in relative extent raster
setwd(paste0(wd,'/OUTPUT/EBA_extents')) #===> Change directory
re <- read.csv('relative_area.csv',header=T)

#-> Load in feature classes
setwd(paste0(wd,'/INPUT/SPATIAL')) #===> Change directory
extent_fc <- readOGR(dsn='VECTOR_INPUT.gdb',layer='mapping_extent',verbose=F)
ra_extent_fc <- readOGR(dsn='VECTOR_INPUT.gdb',layer='risk_analysis_extent',verbose=F)
roads_fc <- readOGR(dsn='VECTOR_INPUT.gdb',layer='Major_roads',verbose=F)
rivers_fc <- suppressWarnings(readOGR(dsn='VECTOR_INPUT.gdb',layer='Major_streams',verbose=F))

#-> Set path to HVRA binary rasters
inrasters <- paste0(wd,'/OUTPUT/EBA_extents/GIS')

#-> Load in extent/template raster
setwd(paste0(wd,'/INPUT/SPATIAL/Extent')) #===> Change directory
ma_extent <- raster('ma_extent.tif')
ra_extent <- raster('ra_extent.tif')

#-> Load in hillshade
setwd(paste0(wd,'/INPUT/SPATIAL')) #===> Change directory
hillshade <- raster('hillcl.tif')
hillshade <- crop(hillshade,extent_fc)
SoG <- colorRampPalette(c('grey0','grey100'))(50)

#############################################END SET UP#############################################

#########################################START CLEAR OUTPUT#########################################

###---> Empty map output folder
e <- try(setwd(paste0(wd,'/OUTPUT/EBA_dcNVC'))) #===> Change directory

#-> Empty output folder...but be careful!
if(class(e)=='try-error'){
	stop('Remove files operation halted due to directory specification error')
}else{
	#-> Empty output folder
	d_list <- list.files()
	if(length(d_list) > 0){
		unlink(d_list[grep('EBA2.log',d_list,invert=T)],recursive=T)
	}
}

#-> Add GIS subdirectory
dir.create('GIS')

##########################################END CLEAR OUTPUT##########################################

#####################################START RELATIVE IMPORTANCE######################################
# Note that no assumptions are made that the weights add up to 100%. Instead, the category
# RIs are coverted to relative measures that sum to one. Then, the HVRA RIs are relativized
# by category and they are combined with multiplication.

#-> Check that all relative importance categories are represented
ri_ri <- unique(ri$Category)
ri_rs <- unique(rs$Category)
extra_ri <- ri_ri[!(ri_ri %in% ri_rs)]
extra_rs <- ri_rs[!(ri_rs %in% ri_ri)]
if(length(extra_ri) > 0){
	if(length(extra_ri)==1){
		cat('Warning:',paste0(extra_ri),
		    'is assigned relative importance but not included.\n')
	}
	if(length(extra_ri)>1){
		cat('Warning:',paste(extra_ri,collapse=', '),
		    'are assigned relative importance but not included.\n')
	}
}
if(length(extra_rs) > 0){
	if(length(extra_rs)==1){
		cat('Warning:',paste0(extra_rs),
		    'is included but not assigned relative importance.\n')
	}
	if(length(extra_rs)>1){
		cat('Warning:',paste(extra_rs,collapse=', '),
		    'are included but not assigned relative importance.\n')
	}
}

#-> Calculate Relative Importance/Relative Extent weights
ri$RI_cat <- ri$RI/sum(ri$RI) # Calculate total relative importance for category
w <- rs[,c('Layer','Category','RI_HVRA')] # Subset important information from rs table
ct <- ddply(w,.(Category),summarize,
            CatTot = sum(RI_HVRA))
w <- merge(w,ct,by='Category',all.x=T)	 
w$RI_HVRA <- w$RI_HVRA/w$CatTot
w$CatTot <- NULL
w <- merge(w,ri[,c('Category','RI_cat')],by='Category',all.x=T)
w$RI <- w$RI_HVRA*w$RI_cat
w <- merge(w,re[,c('Layer','RE')],by='Layer',all.x=T)
w$RIoRE <- w$RI/w$RE
rs <- merge(rs,w[,c('Layer','RIoRE')],by='Layer',all.x=T)

kfields <- c('Layer','HVRA','Category','RI_HVRA','RIoRE')
write.csv(rs[,kfields],'RelImp_over_RelExt_check.csv',row.names=F)

######################################END RELATIVE IMPORTANCE#######################################

##########################################START PROCESSING##########################################

#-> Get list of treatments
trts <- tspecs$FB_Code

#-> Get vector of HVRA categories
cats <- unique(rs$Category)

#-> Create empty list to store benefit totals
bTots <- list()

###---> Calculate delta cNVC by treatment
setwd(paste0(wd,'/OUTPUT/EBA_dcNVC/GIS')) #===> Change directory
for(i in 1:length(trts)){
	cat('#-> Calculating dcNVC for ',paste(trts[i]),'\n',sep='')
	
	#-> Create empty list to store category-level cNVC
	cNVCs <- list()

	#-> Process eNVC by category
	for(j in 1:length(cats)){
		
		cat('Processing ',paste(cats[j]),'\n',sep='')
		
		#-> Subset table
		catrs <- rs[rs$Category==paste(cats[j]),]
		
		#-> Create empty list to store HVRA-level cNVC
		hvra.l <- list()
		
		#-> Create list to store coverage rasters
		cov.l <- list()
		
		#-> Calculate weighted cNVC for each HVRA
		for(k in 1:nrow(catrs)){
			
			#-> Read in raster
			cNVC <- raster(paste0(inrasters,'/HVRA_',catrs$Layer[k],'.tif'))
			
			#-> Save coverage
			coverage <- cNVC; coverage[!is.na(coverage)] <- 1; coverage[is.na(coverage)] <- 0
			cov.l[[k]] <- coverage
			
			#-> Fill null with zero
			cNVC[is.na(cNVC)] <- 0
			
			#-> Apply treatment ecological benefits
			cNVC <- cNVC*catrs[k,paste(trts[i])]
			
			#-> Calculate cNVC with relative importance adjustment
			# RIoRE = relative importance/relative extent
			hvra.l[[k]] <- cNVC*catrs$RIoRE[k]
			
			#-> Get total risk associated with HVRA
			bTot <- catrs[k,c('Layer','HVRA','Category')]
			bTot$Treatment <- as.character(tspecs$Treatment[i])
			bTot$Tot_cNVC <- cellStats(cNVC*catrs$RIoRE[k]*ra_extent,'sum')
			bTots[[length(bTots)+1]] <- bTot
		
		}
		
		if(length(hvra.l) > 1){
		
			#-> Calculate cNVC for category
			cNVC <- sum(stack(hvra.l))
			
			#-> Get coverage
			coverage <- sum(stack(cov.l))
			
		}else{
			
			cNVC <- hvra.l[[1]]
			coverage <- cov.l[[1]]
			
		}
		
		#-> Set cNVC raster to null where there is no coverage
		cNVC[coverage==0] <- NA
		
		#-> Save category raster to file
		writeRaster(suppressWarnings(cNVC*ma_extent),paste0('dcNVC_trt_',i,'_category_',j,'.tif'),
		            format='GTiff',overwrite=T)
		
		#-> Save to category list
		cNVCs[[j]] <- cNVC

	}

	#-> Calculate integrated cNVC across categories
	if(length(cNVCs) > 1){
		ts <- stack(cNVCs) # Stack rasters
		ts[is.na(ts)] <- 0	# Fill nulls with zeros
		cNVC <- sum(ts)
	}else{
		cNVC <- cNVCs[[1]]
		cNVC[is.na(cNVC)] <- 0	# Fill nulls with zeros
	}
		
	#-> Save composite cNVC raster
	writeRaster(cNVC*ma_extent,paste0(trts[i],'_dcNVC.tif'),format='GTiff',overwrite=T)
	
	cat('\n')
}


#-> Save category key
Number <- seq(1,length(cats),1)
ckey <- data.frame(Number,Category=cats)
write.csv(ckey,'category_key.csv',row.names=F)

###---> Compile and save benefit report
setwd(paste0(wd,'/OUTPUT/EBA_dcNVC')) #===> Change directory

#-> Compile list to data frame
bdf <- do.call('rbind',bTots)

#-> Save
write.csv(bdf,'Total_benefit_report.csv',row.names=F)

###########################################END PROCESSING###########################################

####################################################################################################
cat('\nFinished at: ',as.character(Sys.time()),'\n\n',sep='')
sink()
cat('Contents saved to log file.\n',sep='')
cat('Close command window to proceed!\n',sep='')
############################################END MESSAGE#############################################

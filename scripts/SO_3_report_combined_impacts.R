####################################################################################################
#### Report on impacts of combined priorities
#### Date Created: 02/23/2021
#### Last Modified: 02/24/2021
####################################################################################################
# 
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
cat('Report on impacts of combined priorities\n',sep='')
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

#-> Load in run specifications files
setwd(paste0(wd,'/INPUT')) #===> Change directory
budgets <- read.csv('so_budgets.csv',header=T)
tspecs <- read.csv('so_treatment_specs.csv',header=T)
ptab <- read.csv('combined_priorities.csv',header=T)
MinArea <- read.csv('so_project_size.csv',header=T)[1,'MinArea']
MaxArea <- read.csv('so_project_size.csv',header=T)[1,'MaxArea']

#-> Load in treatment units
setwd(paste0(wd,'/INPUT/SPATIAL/Treatment_units')) #===> Change directory
TUs <- shapefile('treatment_units.shp')

#-> Load in individual priorities
setwd(paste0(wd,'/OUTPUT/Treatment_plan_riskreduction')) #===> Change directory
RR_priority <- data.frame(shapefile('treatment_priority.shp'))
RR_bigplan <- data.frame(shapefile(paste0('tplan_',max(budgets$Budget),'.shp')))
setwd(paste0(wd,'/OUTPUT/Treatment_plan_ecobenefits')) #===> Change directory
EB_priority <- data.frame(shapefile('treatment_priority.shp'))
EB_bigplan <- data.frame(shapefile(paste0('tplan_',max(budgets$Budget),'.shp')))

#-> Load in treatment unit metrics
setwd(paste0(wd,'/OUTPUT/Treatment_plan_riskreduction/Supplementary')) #===> Change directory
RR_metrics <- read.csv('TrtUnit_metics.csv',header=T)
setwd(paste0(wd,'/OUTPUT/Treatment_plan_ecobenefits/Supplementary')) #===> Change directory
EB_metrics <- read.csv('TrtUnit_metics.csv',header=T)

#-> Define projection
proj <- CRS('+proj=utm +zone=13 +datum=NAD83')

#-> Define Con function to match ArcGIS syntax
Con <- function(condition, trueValue, falseValue){
	return(condition*trueValue + (!condition)*falseValue)
}

#############################################END SET UP#############################################

###########################################START ANALYSIS###########################################
setwd(paste0(wd,'/OUTPUT/Treatment_plan_combined')) #===> Change directory

###---> Step 1: Merge priorities

#-> Transfer priority fields
RR_priority$RR_Priority <- RR_priority$Priority
EB_priority$EB_Priority <- EB_priority$Priority

#-> Combine
cpt <- data.frame(cbind(UID=TUs$UID,PARK_NAME=TUs$PARK_NAME))
cpt <- merge(cpt,RR_priority[,c('UID','RR_Priority')],by='UID',all.x=T)
cpt <- merge(cpt,EB_priority[,c('UID','EB_Priority')],by='UID',all.x=T)
cpt[is.na(cpt)] <- 'None'

#-> Merge with combined priorities
cpt <- merge(cpt,ptab,by=c('RR_Priority','EB_Priority'),all.x=T)

#-> Merge with feasibility and cost information
cpt <- merge(cpt,RR_metrics[,c('UID','T1_feas','T2_feas','T3_feas','Tot_feas','T1_Cost','T2_Cost',
             'T3_Cost','T1_RR','T2_RR','T3_RR')],by='UID',all.x=T)
cpt <- merge(cpt,EB_metrics[,c('UID','T1_EB','T2_EB','T3_EB')])

#-> Reorganize
cpt <- cpt[,c('UID','PARK_NAME', 'RR_Priority','EB_Priority','ComboRank','T1_feas','T2_feas','T3_feas',
           'Tot_feas','T1_Cost','T2_Cost','T3_Cost','T1_RR','T2_RR','T3_RR','T1_EB','T2_EB',
		   'T3_EB')]
cpt <- cpt[order(cpt$ComboRank),]


###---> Step 2: Simplify treatment plan

#-> Create treatment type fields
cpt$T1_ac <- 0; cpt$T2_ac <- 0; cpt$T3_ac <- 0

#-> Simplify by treatment unit
# Pulling columns for types by name in case order changes in future efforts
thin <- paste0('T',which(tspecs$Treatment=='Thin'),'_ac')
Rx <- paste0('T',which(tspecs$Treatment=='Rx Fire'),'_ac')
thinRx <- paste0('T',which(tspecs$Treatment=='Thin and Rx Fire'),'_ac')
UIDs <- cpt$UID[!is.na(cpt$ComboRank)]
for(i in 1:length(UIDs)){
	#-> Pull assigned acreage from relevant plan
	sRR <- RR_bigplan[RR_bigplan$UID==UIDs[i],]
	sEB <- EB_bigplan[EB_bigplan$UID==UIDs[i],]
	#-> Transfer max acres by type to table
	comp <- rbind(sRR,sEB)
	cpt$T1_ac[i] <- max(comp$T1_ac,na.rm=T)
	cpt$T2_ac[i] <- max(comp$T2_ac,na.rm=T)
	cpt$T3_ac[i] <- max(comp$T3_ac,na.rm=T)
	#-> Combine thin and Rx and transfer if appropriate
	if(sum(cpt[i,thin],cpt[i,Rx],cpt[i,thinRx]) > cpt[i,'Tot_feas']){
		if(cpt[i,thinRx] > 0){
			cpt[i,thinRx] <- cpt[i,paste0('T',which(tspecs$Treatment=='Thin and Rx Fire'),'_feas')]
			cpt[i,thin] <- 0
			rfeas <- cpt[i,'Tot_feas'] - cpt[i,thinRx]
			cpt[i,Rx] <- min(rfeas,cpt[i,paste0('T',which(tspecs$Treatment=='Rx Fire'),'_feas')])
			if(cpt[i,Rx] < MinArea){
				cpt[i,Rx] <- 0
			}
		}else{ # i.e., the combination of thin and Rx are the issue
			cpt[i,thinRx] <- cpt[i,paste0('T',which(tspecs$Treatment=='Thin and Rx Fire'),'_feas')]
			cpt[i,thin] <- 0
			rfeas <- cpt[i,'Tot_feas'] - cpt[i,thinRx]
			cpt[i,Rx] <- min(rfeas,cpt[i,paste0('T',which(tspecs$Treatment=='Rx Fire'),'_feas')])
			if(cpt[i,Rx] < MinArea){
				cpt[i,Rx] <- 0
			}
		}
	}

}

#-> Save treamtment unit summary table
write.csv(cpt,'treatment_unit_summary.csv',row.names=F)

###---> Step 3: Calculate performance metrics by priority level

#-> Calculate metrics
pmd <- ddply(cpt,.(ComboRank),summarize,
		     T1_USD = sum(T1_ac*T1_Cost),
			 T2_USD = sum(T2_ac*T2_Cost),
			 T3_USD = sum(T3_ac*T3_Cost),
			 RR = sum(T1_ac*T1_RR,T2_ac*T2_RR,T3_ac*T3_RR),
			 EB = sum(T1_ac*T1_EB,T2_ac*T2_EB,T3_ac*T3_EB),
			 T1_ac = sum(T1_ac),
			 T2_ac = sum(T2_ac),
			 T3_ac = sum(T3_ac),
			 Tot_ac = sum(T1_ac,T2_ac,T3_ac),			 
			 Tot_USD = sum(T1_USD,T2_USD,T3_USD))

#-> Merge to priority table
ptab <- merge(ptab,pmd,by='ComboRank',all.x=T)

#-> Calculate cumulative metrics
ptab[is.na(ptab)] <- 0
ptab <- ptab[order(ptab$ComboRank),]
ptab$cT1_ac <- cumsum(ptab$T1_ac)
ptab$cT2_ac <- cumsum(ptab$T2_ac)
ptab$cT3_ac <- cumsum(ptab$T3_ac)
ptab$cTot_ac <- cumsum(ptab$Tot_ac)
ptab$cTot_USD <- cumsum(ptab$Tot_USD)
ptab$cRR <- cumsum(ptab$RR)
ptab$cEB <- cumsum(ptab$EB)

#-> Save summary table
write.csv(ptab[ptab$ComboRank > 0,],'priority_level_summary.csv',row.names=F)

#-> Summary figure
tiff('Avoided_impact_curve.tif',width=1000,height=1250,compression='lzw',pointsize=22,
     type='windows')
par(mfrow=c(3,1))

cols <- c('forestgreen','red','orange','sienna','purple','cyan')
xtics <- seq(0,floor(max(ptab$cTot_USD)/10000000),1)*10000000

#-> Avoided impact curve
ylim <- c(0,max(ptab$cRR)*1.05)
plot(cRR~cTot_USD,data=ptab,type='l',xaxs='i',yaxs='i',ylim=ylim,lwd=2,axes=F,
     ylab='Risk Reduction',xlab='Fuel treatment budget ($)',main='Risk Reduction by Budget')
axis(1,at=xtics,labels=paste0('$',xtics/1000000,'M'))
axis(2,at=pretty(ptab$cRR),labels=format(pretty(ptab$cRR),big.mark=',',big.interval=3L,
     trim=T,scientific=F))
box()
 
#-> Ecological benefits curve
ylim <- c(0,max(ptab$cEB)*1.05)
plot(cEB~cTot_USD,data=ptab,type='l',xaxs='i',yaxs='i',ylim=ylim,lwd=2,axes=F,
     ylab='Ecological Benefit',xlab='Fuel treatment budget ($)',main='Ecological Benefit by Budget')
axis(1,at=xtics,labels=paste0('$',xtics/1000000,'M'))
axis(2,at=pretty(ptab$cEB),labels=format(pretty(ptab$cEB),big.mark=',',big.interval=3L,
     trim=T,scientific=F))
box()

#-> Treatment allocation curves
ylim <- c(0,max(ptab$cTot_ac,na.rm=T)*1.05)
plot(cTot_ac~cTot_USD,data=ptab,xaxs='i',yaxs='i',ylim=ylim,lwd=2,col=NA,axes=F,
     ylab='Treated Acres',xlab='Fuel treatment budget ($)',main='Treatment Allocation')
for(i in 1:nrow(tspecs)){
	lines(ptab$cTot_USD,ptab[,paste0('cT',i,'_ac')],lwd=2,col=cols[i])
}
lines(cTot_ac~cTot_USD,data=ptab,lwd=2,col='black')	 
axis(1,at=xtics,labels=paste0('$',xtics/1000000,'M'))
axis(2,at=pretty(c(0,max(ylim))),labels=format(pretty(c(0,max(ylim))),big.mark=',',
     big.interval=3L,trim=T,scientific=F))
box()

legend('topleft',c(as.character(tspecs$Treatment),'Total'),
       col=c(cols[1:nrow(tspecs)],'black'),lwd=2,bty='n')
g <- dev.off()
	
###---> Save plan as shapefile

#-> Merge priority information to treatment units
TUs <- merge(TUs,cpt[,c('UID','T1_ac','T2_ac','T3_ac','RR_Priority','EB_Priority','ComboRank')],
             by='UID',all.x=T)
TUs <- TUs[!is.na(TUs$ComboRank) & TUs$ComboRank > 0,]

#-> Save as shapefile
writeOGR(obj=TUs,dsn='.',layer='treatment_priority',driver='ESRI Shapefile',overwrite=T)

#-> Quick visual assessment in R
ComboRank <- seq(min(ptab$ComboRank),max(ptab$ComboRank),1)
Cols <- colorRampPalette(c('red','yellow','forestgreen'))(length(ComboRank))
ctab <- data.frame(ComboRank,Cols)
test <- merge(TUs,ctab,by='ComboRank')
plot(test,col=as.character(test$Cols))

############################################END ANALYSIS############################################

####################################################################################################
cat('\nFinished at: ',as.character(Sys.time()),'\n\n',sep='')
cat('Close command window to proceed!\n',sep='')
############################################END MESSAGE#############################################

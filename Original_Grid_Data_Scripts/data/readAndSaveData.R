rm(list=ls(all=TRUE))
Sys.setlocale('LC_ALL','C')
options(scipen = 999)

# ===== Settings =======

readIn=T #T for reading in and calculation
ncells=67420
ny=105
cluster=F
plotting=F

if (cluster){
  iFol="/Volumes/RachelExternal/Thesis/DataFabian/data/"
  oFol="/Volumes/RachelExternal/Thesis/DataFabian/plots/"
}else{
  iFol="/Volumes/RachelExternal/Thesis/DataFabian/data/"
  oFol="/Volumes/RachelExternal/Thesis/DataFabian/plots/"
}

parallel=F # no massive parallelization, can still run on 4 cpus of local machine
if(cluster) {
  .libPaths(paste0("/p/projects/lpjml/R.3.4.4/library"))
  if(require(Rmpi)) { # R implementation of MPI interface
    if(require(doMPI)) { # interface for foreach construct to run in MPI parallel mode
      cl <- startMPIcluster() # start cluster (link R instances together)
      num.cluster <- clusterSize(cl)
      parallel <- TRUE
      if(num.cluster>1) {
        # we are using more than 1 CPU, so really run in parallel mode
        registerDoMPI(cl) # tells foreach to use MPI parallel mode
        print(paste("Running in parallel mode on",num.cluster,"worker nodes."))
      } else {
        registerDoSEQ() # tells foreach to use sequential mode
        print("Running in sequential mode.")
      }
    } else {
      library(foreach)
      registerDoSEQ() # tells foreach to use sequential mode
      print("Running in sequential mode.")
    }
  } else {
    library(foreach)
    registerDoSEQ() # tells foreach to use sequential mode
    print("Running in sequential mode.")
  }
} else {
  library(foreach)
  registerDoSEQ() # tells foreach to use sequential mode
  print("Running in sequential mode.")
}
# ======== functions ==========
## reads a header from an LPJ input file                              ##
## tries to determine header version unless force_version is provided ##
## return value: list with 3 components:                              ##
## - header name, e.g. LPJGRID                                        ##
## - header values (11 in total), if header version is <3, partially  ##
##   filled with default values                                       ##
## - endian of file (little or big)                                   ##
readheader <- function(filename, force_version=NULL) {
  if(!file.exists(filename)) {
    stop(paste("Error in readheader:", filename, "does not exist"))
  }
  zz <- file(filename, "rb")
  headername <- rawToChar(readBin(zz, raw(), n=30), multiple=TRUE)
  headername <- headername[1:(min(which(!grepl("[[:alpha:]_]", headername)))-1)]
  headername <- paste(headername, collapse="")
  if(substr(headername, 1,3) != "LPJ") {
    close(zz)
    stop(paste("Error in readheader: invalid header name", headername))
  }
  seek(zz, nchar(headername))
  endian <- .Platform$endian
  version <- readBin(zz, integer(), size=4, n=1, endian=endian)
  if(bitwAnd(version, 0xff)==0) {
    endian <- ifelse(endian=="little", "big", "little")
    seek(zz, nchar(headername))
    version <- readBin(zz, integer(), size=4, n=1, endian=endian)
  }
  if(!is.null(force_version)) {
    print(paste("Forcing header version to", force_version))
    version <- force_version
  }
  headerdata <- readBin(zz, integer(), size=4, n=6, endian=endian)
  names(headerdata) <- c("order", "firstyear", "nyear", "firstcell", "ncell", "nbands")
  if(version == 2) {
    headerdata <- c(headerdata, readBin(zz, double(), size=4, n=2, endian=endian))
    names(headerdata) <- c(names(headerdata[1:6]), "cellsize_lon", "scalar")
  }
  if(version == 3) {
    headerdata <- c(headerdata, readBin(zz, double(), size=4, n=3, endian=endian))
    headerdata <- c(headerdata, readBin(zz, integer(), size=4, n=1, endian=endian))
    names(headerdata) <- c(names(headerdata[1:(length(headerdata)-4)]), "cellsize_lon", "scalar", "cellsize_lat", "datatype")
  } else {
    if(length(headerdata)==6) {
      headerdata <- c(headerdata, cellsize_lon=0.5, scalar=1, cellsize_lat=0.5, datatype=1)
      warning("Type 1 header. Adding default values for cellsize, scalar and datatype which may not be correct in all cases")
    }
    if(length(headerdata)==8) {
      headerdata <- c(headerdata, cellsize_lat=as.double(headerdata["cellsize_lon"]), datatype=1)
      warning("Type 2 header. Adding default value for datatype which may not be correct in all cases")
    }
  }
  close(zz)
  return(list(name=headername, header=c(version=version, headerdata), endian=endian))
}

#read LPJmL input with header
autoReadInput <- function(inFile,getyearstart=-1,getyearstop=-1,manu=F,msize=4){
  hdr=readheader(filename=inFile)$header
  print(hdr)
  startyear=hdr[3]
  stopyear=hdr[3]+hdr[4]-1
  ncells=hdr[6]
  nbands=hdr[7]
  if (getyearstart==-1){
    getyearstart=startyear
  }
  if (getyearstop==-1){
    getyearstop=stopyear
  }
  if (hdr[1]==1){#header version 1 
    headersize=36
  }else if (hdr[1]==2){#header version 2
    headersize=43
  }else{ #header version 3
    headersize=51
  }
  if (length(hdr)>10){
    if (hdr[11]==0){
      size=1
      inputType="char"
    }else if(hdr[11]==1){
      size=2
      inputType="integer"
    }else if(hdr[11]==2){
      size=4
      inputType="integer"
    }else if(hdr[11]==3){
      size=4
      inputType="double"
    }else if(hdr[11]==4){
      size=8
      inputType="double"
    }
  }
  if (manu){size=msize}
  if (getyearstop>stopyear){
    stop(paste("unexpected usage: getyearstop (",getyearstop,") larger than stopyear (",stopyear,") -- stopping"))
  }
  if (getyearstart<startyear){
    stop(paste("unexpected usage: getyearstart (",getyearstart,") smaller than startyear (",startyear,") -- stopping"))
  }
  nyears=getyearstop-getyearstart+1
  input <- file(inFile,"rb")
  seek(input, where = headersize+(getyearstart-startyear)*ncells*size*nbands, origin="start")
  if (inputType == "integer"){
    dataIn <- readBin(input,integer(),n = nyears*ncells*nbands, size=size)
  }else if (inputType == "double"){
    dataIn <- readBin(input,double(),n = nyears*ncells*nbands, size=size)
  }else{
    dataIn <- readBin(input,character(),n = nyears*ncells*nbands, size=size)
  }
  close(input)      #remove to save space
  #print(paste("nyears,nbands,ncells,size,headersize=",nyears,nbands,ncells,size,headersize))
  if (nyears==1){
    dim(dataIn) <- c(nbands,ncells)
  }else{
    dim(dataIn) <- c(nbands,ncells,nyears)
  }
  return(dataIn*hdr[["scalar"]])
}

#converting degrees to radians 
deg2rad <- function(deg){ 
  return(deg*pi/180)
}

#calculate distance of two points on a sphere
surfDistance <- function(lat1, long1, lat2, long2){ 
  #Convert the latitudes and longitudes from degree to radians.
  lat1 = deg2rad(lat1)
  long1 = deg2rad(long1)
  lat2 = deg2rad(lat2)
  long2 = deg2rad(long2)
  
  #Haversine Formula 
  dlong = long2-long1
  dlat = lat2-lat1
  ans = sin(dlat/2)^2 + cos(lat1) * cos(lat2) * sin(dlong/2)^2
  ans = 2 * asin(sqrt(ans))
  
  #return result in km
  return(ans * 6371)
} 

#returns a list of indizes of the cells in radius maxDist around the given cell
getVicinityCells <- function(cell,grid,getIndexFrame,maxDist,ires,prec=2){
  clat=grid[2,cell]
  clon=grid[1,cell]
  ny=floor(maxDist/surfDistance(lat1=clat,long1=clon,lat2=clat+1/ires,long2=clon))#lats
  nx=floor(maxDist/surfDistance(lat1=clat,long1=clon,lat2=clat,long2=clon+1/ires))#lons
  y=c(rev(seq(from = -1/ires,by = -1/ires,length.out = ny)),0,seq(from =1/ires,by=1/ires,length.out = ny))#lats
  x=c(rev(seq(from = -1/ires,by = -1/ires,length.out = nx)),0,seq(from =1/ires,by=1/ires,length.out = nx))#lons
  offset=expand.grid(x,y)
  coords=offset+matrix(rep(c(clon,clat),length(offset[,1])),nrow = length(offset[,1]),byrow = T)
  lara=round(range(grid[2,]),prec)#grid is in lon,lat
  lora=round(range(grid[1,]),prec)#grid is in lon,lat
  inds=array(0,dim=length(coords[,1]))
  for (i in 1:length(coords[,1])){
    loc=round(getNearestCellCoords(coords[i,],ires = ires),prec)
    if (loc[1]<lora[1] | loc[1]>lora[2]){inds[i]=NA}
    else if (loc[2]<lara[1] | loc[2]>lara[2]){inds[i]=NA}
    else{inds[i]=getIndexFrame[[as.character(loc[1]),as.character(loc[2])]]}
  }
  return(inds)
}

# returns an array - for each cell the distance to the closest 
#    grid cell with irrigation, 0 if that cell is already irrigated
# relies on global variables:
#   ncells,irrFrac,lon,lat
# and functions:
#   surfDistance,deg2rad
distanceToNearestIrrigPlot <- function(y){
  distNextIrr=array(99999,dim=c(ncells))#just for year y
  for (c in 1:ncells){
    if (c%%5000==0){write(paste0("y:",y,", c: ",c),stdout())}
    if (irrFrac[c,y]>0){
      distNextIrr[c]=0
      next
    }
    irrCells=which(irrFrac[,y]>0)
    for(n in irrCells){
      ncdist=surfDistance(lat1 = lat[c],long1 = lon[c],lat2 = lat[n],long2 = lon[n])
      if (ncdist<distNextIrr[c]){distNextIrr[c]=ncdist}
    }
  }
  return(distNextIrr)
}

#plot a global LPJmL array (67420 cells) with manual breaks and colorramp 
plotGlobalMan <- function(data,file,title,brks,palette="YlGnBu",legendtitle,legYes,eps){
  if (!length(palette)==(length(brks)-1)){colorRampPalette(RColorBrewer::brewer.pal(9,palette))(length(brks)-1)}
  ires=2
  legendticks=seq(from=0,to = 100,length.out = length(brks))
  data[data<brks[1]] <- brks[1]
  data[data>brks[length(brks)]] <- brks[length(brks)]
  if (eps){
    file=strsplit(file,".",fixed=TRUE)[[1]]
    file=paste(c(file[1:(length(file)-1)],"eps"),collapse=".")
    ps.options(family = c("Helvetica"), pointsize = 18)
    postscript(file,horizontal = FALSE, onefile = FALSE, width=22, height=8.5,paper="special")
  }else{
    png(file, width=800*ires, height=400*ires, units="px",res=400,pointsize = 4)
  }
  ra <- raster::raster(ncols=360*ires, nrows=180*ires)
  range <- range(data)
  ra[raster::cellFromXY(ra,cbind(lon,lat))] <-  data
  extent <- raster::extent(c(-180, 180, -60, 90))
  if (legYes){
    par(bty="n",oma=c(0,0,0,0),mar=c(0,0,0,3),xpd=T)
  }else{
    par(bty="n",oma=c(0,0,0,0),mar=c(0,0,0,0))
  }
  raster::plot(ra,ext=extent,breaks=brks,col=palette,main="",legend=FALSE,axes=FALSE,maxpixels=360*180*ires*ires)
  title(title,line=-1)
  if (legYes){
    fields::image.plot(legend.only=TRUE,zlim=range(brks),col = palette,useRaster=FALSE,breaks=legendticks,
                       lab.breaks=round(brks,2),legend.shrink = 0.8,legend.args=list(legendtitle,side=3, font=2, line=1))
  }
  maps::map('world',add=TRUE,res=0, lwd=0.1,ylim=c(-60,90))
  dev.off()
}

#plot global lpj array
plotGlobalWflex <- function(data,file,title,man=F,mbrk=-1,mpalette=-1,pow2max=10,pow2min=0,colPos="Blues",colNeg="Reds",legendtitle,legYes,onlyPos=F,lon,lat,printBorders=T,logscale=T){
  ext=c(-180, 180, -60, 90)
  xmagn=1
  ires=2
  if (!man){
    if (logscale){
      if (onlyPos){
        legendticks <- c(0,2^seq(pow2min,pow2max,1))
        brks <- c(seq(pow2min,pow2max,length.out = length(legendticks)))
        palette <- c("white",colorRampPalette(RColorBrewer::brewer.pal(9,colPos))(length(legendticks)-2))  
      }else{
        legendticks <- c(-(2^seq(pow2max,pow2min,-1)),2^seq(pow2min,pow2max,1))
        brks <- seq(-pow2max,pow2max,length.out = length(legendticks))
        palette <- c(rev(colorRampPalette(RColorBrewer::brewer.pal(9,colNeg))(length(legendticks)/2-1)),"white",colorRampPalette(RColorBrewer::brewer.pal(9,colPos))(length(legendticks)/2-1))
      }
    }else{ #linscale
        lout=1+2*7
        if (onlyPos){
          legendticks <- round(seq(0,(2^pow2max),length.out = lout),0)
          brks <- legendticks
          palette <- c("white",colorRampPalette(RColorBrewer::brewer.pal(9,colPos))(length(legendticks)-2))
        }else{
          legendticks <- round(seq(-(2^pow2max),(2^pow2max),length.out = lout),0)
          brks <- legendticks
          palette <- c(rev(colorRampPalette(RColorBrewer::brewer.pal(9,colNeg))(length(legendticks)/2-1)),"white",colorRampPalette(RColorBrewer::brewer.pal(9,colPos))(length(legendticks)/2-1))
        }
    }
  }else{ #man
        brks=mbrk
        legendticks=mbrk
        palette=mpalette
      }

  
  data[data<legendticks[1]] <- legendticks[1]
  data[data>legendticks[length(legendticks)]] <- legendticks[length(legendticks)]

  fileP=strsplit(file,".",fixed=TRUE)[[1]]
  expFormat=fileP[-1]
  if (expFormat=="eps"){
    postscript(file,horizontal = FALSE, onefile = FALSE, width=14, height=8,paper="special",family = c("Helvetica"), pointsize = 12)
  }else if (expFormat=="pdf"){  
    pdf(file,width=14,height=8,paper="special",family = c("Helvetica"),pointsize = 12)
  }else if (expFormat=="png"){ 
    png(file, width=7, height=3.5, units="in", res=400, pointsize=6,type="cairo")
  }else{
    print("Please specify export image format for image as one of '.pdf', '.png', or '.eps'")
  }
  ra <- raster::raster(ncols=360*ires, nrows=180*ires)
  range <- range(data)
  ra[raster::cellFromXY(ra,cbind(lon,lat))] <-  data
  extent <- raster::extent(ext)
  if (legYes){
    par(bty="n",oma=c(0,0,0,0),mar=c(0,0,0,3),xpd=T)
  }else{
    par(bty="n",oma=c(0,0,0,0),mar=c(0,0,0,0))
  }
  raster::plot(ra,ext=extent,breaks=legendticks,col=palette,main="",legend=FALSE,axes=FALSE,maxpixels=360*180*ires*ires)
  title(title,line=-1)
  if (legYes){
    fields::image.plot(legend.only=TRUE,zlim=range(data),col = palette,useRaster=FALSE,breaks=brks,lab.breaks=round(legendticks,2),legend.shrink = 0.8,#legend.width = 2*xmagn,
                       legend.args=list(legendtitle,side=3, font=2, line=1))
  }
  if(printBorders){maps::map('world',add=TRUE,res=0, lwd=0.1,ylim=c(-60,90))}
  dev.off()
}
#example usage
#plotGlobalWflex(data = wd,file = paste("/home/stenzel/withdrawals_lin.png"),title = "withdrawals in 2090",pow2max = 6,pow2min = -2,legendtitle = "km3/yr",legYes = T,onlyPos = T,printBorders = T,lon=lon,lat=lat,logscale = F)

# ========= main ============

#library(lpjmliotools)
source("/home/stenzel/Nextcloud_PIK/7_irrigation_History/data/lpjmliotools_20201104.R")
load(file=paste0(oFol,"latlon.RData"))
cellarea <- (111194.9/2)*(111194.9/2)*cos(lat/180*pi) # cellarea in m2 


if (readIn){
cft_landuse=readCFT(inFile = paste0(iFol,"landuse/cftfrac.bin"),startyear = 1901,stopyear = 2005,bands = 32,size = 4,getyearstart = 1901,headersize = 0,getyearstop = 2005)
cft_allcrops=readCFT(inFile = paste0(iFol,"allcrops/cftfrac.bin"),startyear = 1901,stopyear = 2005,bands = 32,size = 4,getyearstart = 1901,headersize = 0,getyearstop = 2005)
harvest_allcrops=readCFT(inFile = paste0(iFol,"allcrops/pft_harvest.pft.bin"),startyear = 1901,stopyear = 2005,bands = 32,size = 4,getyearstart = 1901,headersize = 0,getyearstop = 2005)
discharge_landuse=readMonthly(inFile = paste0(iFol,"landuse/mdischarge.bin"),startyear = 1901,stopyear = 2005,size = 4,getyearstart = 1901,headersize = 0,getyearstop = 2005)

#calculate total irr fraction
irrFrac=apply(cft_landuse[,17:32,],c(1,3),sum)
globalIrrigArea=colSums(irrFrac*cellarea,na.rm = T)*10^-10 #from m2 to Mha

#calculate cft specific harvest increase through (potential) irrigation
increase=(harvest_allcrops[,17:32,]-harvest_allcrops[,1:16,])/harvest_allcrops[,1:16,]
increase[which(harvest_allcrops[,17:32,]<=0)]=0 #remove cells which do not have harvests -> result in INF/NaN values
increase[which(harvest_allcrops[,1:16,]<=0)]=0
mean_increase=apply(increase,c(1,3),mean)
median_increase=apply(increase,c(1,3),median)

#add runoff_landuse, prec_landuse, transp_green_landuse, transp_blue_landuse, evap_landuse, interc_landuse
runoff_landuse=readMonthly(inFile = paste0(iFol,"landuse/mrunoff.bin"),startyear = 1901,stopyear = 2005,size = 4,getyearstart = 1901,headersize = 0,getyearstop = 2005)
prec_landuse=readMonthly(inFile = paste0(iFol,"landuse/mprec.bin"),startyear = 1901,stopyear = 2005,size = 4,getyearstart = 1901,headersize = 0,getyearstop = 2005)
transp_green_landuse=readMonthly(inFile = paste0(iFol,"landuse/mtransp.bin"),startyear = 1901,stopyear = 2005,size = 4,getyearstart = 1901,headersize = 0,getyearstop = 2005)
transp_blue_landuse=readMonthly(inFile = paste0(iFol,"landuse/mtransp_b.bin"),startyear = 1901,stopyear = 2005,size = 4,getyearstart = 1901,headersize = 0,getyearstop = 2005)
evap_landuse=readMonthly(inFile = paste0(iFol,"landuse/mevap.bin"),startyear = 1901,stopyear = 2005,size = 4,getyearstart = 1901,headersize = 0,getyearstop = 2005)
interc_landuse=readMonthly(inFile = paste0(iFol,"landuse/minterc.bin"),startyear = 1901,stopyear = 2005,size = 4,getyearstart = 1901,headersize = 0,getyearstop = 2005)

system.time({
  if (cluster){
    #require(parallel)
    distNextIrr <- foreach(y = 1:ny, .inorder=FALSE, .combine="c", .verbose=T) %dopar% {
      return(distanceToNearestIrrigPlot(y))
    }
  }else{
    require(parallel)
    numCores <- detectCores()
    ret3 <- mclapply(1:ny,distanceToNearestIrrigPlot,mc.cores = numCores)
    distNextIrr=array(unlist(ret3))
    dim(distNextIrr)=c(ncells,ny)
  }
})#end system.time

if (cluster){
  if(parallel) {
    closeCluster(cl)
    mpi.quit()
  }
}

for (y in 1:105){
  
}

#reading country-list
require(plyr)
country_codes=data.frame(t(autoReadInput(inFile = paste0(iFol,"landuse/cow_full_2018.bin"))))
colnames(country_codes)=c("COW","REG")
country_names<-read.csv(file="/media/Storage/irrigation/countrylist.csv",header=T)
lpj_COW_ISO<-read.csv(file="/media/Storage/irrigation/country_list_LPJmL_ISO_codes.csv",header=T)
colnames(lpj_COW_ISO)=c("Countryname","ISO")
lpj_countries=merge(lpj_COW_ISO,country_names,by="Countryname")
fullCountryData=join(country_codes,lpj_countries)
#save(fullCountryData,file="/media/Storage/LOCOMOTION/coutry_data.RData",version = 2)

#GDP from Maddison project
countryGDP <- read.csv(file = paste('/media/Storage/irrigation/GDP/mpd2018.csv',sep=""),stringsAsFactors = FALSE)
names(countryGDP)[names(countryGDP) == "countrycode"] <- "ISO"
lpjGDPpc=array(0,dim=c(ncells,ny))
for (y in 1901:2005){
  yearGDP=subset(countryGDP,year==y)
  test <- join(fullCountryData,yearGDP)
  lpjGDPpc[,(y-1900)]=test$cgdppc
}

#Population
lpjPOP=array(0,dim=c(ncells,ny))
popdens=autoReadInput(inFile="/media/Storage/irrigation/popdens_HYDE3_1901_2011_bi.clm")
dim(popdens)=c(ncells,111)
plotGlobalW(data = popdens[,99],file = paste0(oFol,"popDens_",2000,".png"),title = "",pow2max = 11,pow2min = 0,legendtitle = "kcap",legYes = T,eps = F,onlyPos = T)



save(runoff_landuse,prec_landuse,transp_green_landuse,transp_blue_landuse,evap_landuse,interc_landuse,popdens,distNextIrr,mean_increase,
     median_increase,globalIrrigArea,irrFrac,discharge_landuse,lpjGDPpc,lat,lon,cellarea,file=paste0(oFol,"data.RData"),version=2)
}
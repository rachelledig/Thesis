### ================ global variable definitions ============
ndays=c(31,28,31,30,31,30,31,31,30,31,30,31)

### ================ read/write routines ===================

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


#' Transform even grid to lpjml grid
#'
#' Transform even grid (e.g. from a netcdf) to lpjml grid
#'
#' @param gridIn even grid to be transformed (e.g. c(720,360) for 0.5 deg)
#' @param lonIn longitude vector of input grid
#' @param latIn latitude vector of input grid
#' @param res resolution (default 0.5)
#'
#' @return outlist transformed data in lpjml format
#'
#' @examples
#' evenGrid2lpjGrid(gridIn=climateCategories,lonIn=seq(-179.75,179.75,0.5),latIn=seq(-89.75,89.75,0.5),res=0.5)
#'
#' @export
evenGrid2lpjGrid <- function(gridIn,lonIn,latIn,ncells,res=0.5){
  l=length(dim(gridIn))
  if (l==2){
    outlist <- array(0,dim=c(ncells))
  }else if(l==3){
    outlist <- array(0,dim=c(ncells,length(gridIn[1,1,])))
  }
  if (l==2){
    for (i in 1:ncells){
      outlist[i] <- gridIn[which(round(lonIn,2)==lon[i]),which(round(latIn,2)==lat[i])] # minus lat[i] only if the grid is flipped in lat
    }
  }else if (l==3){
    for (y in c(1:length(gridIn[1,1,])) ){
      for (i in 1:ncells){
        outlist[i,y] <- gridIn[which(round(lonIn,2)==lon[i]),which(round(latIn,2)==lat[i]),y] # minus lat[i] only if the grid is flipped in lat
      }
    }  
  }
  return(outlist)
}


#' Get header of LPJmL input file
#'
#' Reads the header of a binary CLM(2) LPJmL input file and returns the content as data_frame
#'
#' @param f.in character string containing the file to read the header from
#'
#' @return data_frame of header data
#'     h$v : header version,
#'     h$o : order,
#'     h$fy : firstyear,
#'     h$ny : nyear,
#'     h$fc : firstcell,
#'     h$nc : ncells,
#'     h$nb : nbands,
#'     h$r : resolution,
#'     h$b : scaling
#'
#' @examples
#' freadheader(system.file("extdata", "grid.bin", package = "lpjmliotools", mustWork = TRUE))
#'
#' @export
freadheader <- function(f.in){
  h <- list()
  h$name <- readChar(f.in, 7)                  # header name
  #h$name <- readBin(f.in,character(), size=8)   # header name
  h$v <- readBin(f.in, integer(), n=1, size=4) # header version
  h$o <- readBin(f.in, integer(), n=1, size=4) # order
  h$fy <- readBin(f.in,integer(), n=1, size=4) # firstyear
  h$ny <- readBin(f.in,integer(), n=1, size=4) # nyear
  h$fc <- readBin(f.in, integer(),n=1,size=4)  # firstcell
  h$nc <- readBin(f.in,integer(),n=1,size=4)   # ncells
  h$nb <- readBin(f.in,integer(),n=1,size=4)   # nbands
  h$r <- readBin(f.in,double(),n=1, size=4)    # resolution
  h$b <- readBin(f.in,double(),n=1, size=4)    # scaling
  return(h)
}

#' Get country code list
#'
#' Reads and returns the list of country codes contained in
#'
#' @param fileName character string containing the file to read the cow from ()
#'
#' @return array of 2 IDs (countrycode,regioncode) for each of the 67420 gridcells
#'
#' @examples
#' lpjCOW(system.file("extdata", "cow_mg_2006.bin", package = "lpjmliotools", mustWork = TRUE))
#'
#' @export
lpjCOW <- function(fileName,ncells){
  zz    <- file(fileName,"rb")
  seek(zz, where = 43, origin="start")
  cow <- readBin(zz, integer(), 2*ncells, size=2)
  dim(cow)=c(2,ncells)
  close(zz)
  return(cow)
}

#' Read monthly LPJmL output file
#'
#' Returns a range of years from a monthly LPJmL output of dimension c(ncells,12,nyears)
#'
#' @param inFile character string containing the file to read the data from
#' @param startyear absolute startyear of output
#' @param stopyear absolute stopyear of output
#' @param size size of each cell's data (2 or 4 bytes)
#' @param headersize size of header data (defaults to 0 bytes)
#' @param getyearstart start of range to return
#' @param getyearstop end of range to return
#'
#' @return array of monthly data c(ncells,12,nyears) for each of the 67420 gridcells over requested range of years
#'
#' @examples
#' \dontrun{
#' readMonthly(inFile="mwateramount.bin",startyear=1861,stopyear=2005,size=4,headersize=0,
#'             getyearstart=1984,getyearstop=2005)
#' }
#'
#' @export
readMonthly <- function(inFile,startyear,stopyear,size,headersize=0,getyearstart,getyearstop,ncells=67420){
  if (getyearstop>stopyear){
    stop(paste("unexpected usage: getyearstop (",getyearstop,") larger than stopyear (",stopyear,") -- stopping"))
  }
  if (getyearstart<startyear){
    stop(paste("unexpected usage: getyearstart (",getyearstart,") smaller than startyear (",startyear,") -- stopping"))
  }
  nyears=getyearstop-getyearstart+1
  input <- file(inFile,"rb")
  seek(input, where = headersize+(getyearstart-startyear)*12*ncells*size, origin="start")
  monthly <- readBin(input,double(),n = nyears*ncells*12, size=size)
  close(input)      #remove to save space
  if (nyears==1){
    dim(monthly) <- c(ncells,12)
  }else{
    dim(monthly) <- c(ncells,12,nyears)
  }
  return(monthly)
}

#' Read yearly LPJmL input file
#'
#' Returns a range of years from a LPJmL input of dimension c(ncells,nyears)
#'
#' @param inFile character string containing the file to read the data from
#' @param startyear absolute startyear of output
#' @param stopyear absolute stopyear of output
#' @param size size of each cell's data (2 or 4 bytes)
#' @param inputType type of variable to read, integer/double
#' @param headersize size of header data (defaults to 43 bytes)
#' @param getyearstart start of range to return
#' @param getyearstop end of range to return
#'
#' @return array of data c(ncells,nyears) for each of the 67420 gridcells over requested range of years
#'
#' @examples
#' \dontrun{
#' readYearlyInput(inFile="wateruse_1900_2005.bin",startyear=1900,stopyear=2005,size=4,headersize=43,
#'            getyearstart=1984,getyearstop=2005)
#' }
#'
#' @export
readYearlyInput <- function(inFile,startyear,stopyear,size,inputType,headersize=43,getyearstart,getyearstop,ncells=67420){
  if (getyearstop>stopyear){
    stop(paste("unexpected usage: getyearstop (",getyearstop,") larger than stopyear (",stopyear,") -- stopping"))
  }
  if (getyearstart<startyear){
    stop(paste("unexpected usage: getyearstart (",getyearstart,") smaller than startyear (",startyear,") -- stopping"))
  }
  nyears=getyearstop-getyearstart+1
  input <- file(inFile,"rb")
  seek(input, where = headersize+(getyearstart-startyear)*ncells*size, origin="start")
  if (inputType == "integer"){
    dataIn <- readBin(input,integer(),n = nyears*ncells, size=size)
  }else if (inputType == "double"){
    dataIn <- readBin(input,double(),n = nyears*ncells, size=size)
  }else{
    print("unknown input type")
  }
  close(input)      #remove to save space
  if (nyears==1){
    dim(dataIn) <- c(ncells)
  }else{
    dim(dataIn) <- c(ncells,nyears)
  }
  return(dataIn)
}


#' Read yearly LPJmL output file
#'
#' Returns a range of years from a LPJmL output of dimension c(ncells,nyears)
#'
#' @param inFile character string containing the file to read the data from
#' @param startyear absolute startyear of output
#' @param stopyear absolute stopyear of output
#' @param size size of each cell's data (2 or 4 bytes)
#' @param headersize size of header data (defaults to 0 bytes)
#' @param getyearstart start of range to return
#' @param getyearstop end of range to return
#' @param ncells number of lpj cells (67420 for 30min res, 2298847 for 5min res)
#'
#' @return array of data c(ncells,nyears) for each of the 67420 gridcells over requested range of years
#'
#' @examples
#' \dontrun{
#' readYearly(inFile="hdates.bin",startyear=1901,stopyear=2005,size=2,headersize=0,
#'            getyearstart=1984,getyearstop=2005)
#' }
#'
#' @export
readYearly <- function(inFile,startyear,stopyear,size,headersize=0,getyearstart,getyearstop,ncells=67420){
  if (getyearstop>stopyear){
    stop(paste("unexpected usage: getyearstop (",getyearstop,") larger than stopyear (",stopyear,") -- stopping"))
  }
  if (getyearstart<startyear){
    stop(paste("unexpected usage: getyearstart (",getyearstart,") smaller than startyear (",startyear,") -- stopping"))
  }
  nyears=getyearstop-getyearstart+1
  input <- file(inFile,"rb")
  seek(input, where = headersize+(getyearstart-startyear)*ncells*size, origin="start")
  if (size==2){
    dataIn <- readBin(input,integer(),n = nyears*ncells, size=size)
  }else if (size==4){
    dataIn <- readBin(input,double(),n = nyears*ncells, size=size)
  }else{
    print("unknown data size")
  }
  close(input)      #remove to save space
  if (nyears==1){
    dim(dataIn) <- c(ncells)
  }else{
    dim(dataIn) <- c(ncells,nyears)
  }
  return(dataIn)
}

#' Read cft LPJmL output
#'
#' Returns a range of years from a cft LPJmL output of dimension c(ncells,bands,nyears).
#'
#' @param inFile character string containing the file to read the data from
#' @param startyear absolute startyear of output
#' @param stopyear absolute stopyear of output
#' @param bands number of bands (32 for standard LPJmL output, 64 for standard input)
#' @param size size of each cell's data (2 or 4 bytes)
#' @param headersize size of header data (defaults to 0 bytes)
#' @param getyearstart start of range to return
#' @param getyearstop end of range to return
#' @param ncells number of lpj cells (67420 for 30min res, 2298847 for 5min res)
#'
#' @return array of data c(ncells,bands,nyears) for each of the 67420 gridcells over requested range of years
#
#' @examples
#' \dontrun{
#' readCFT(inFile="cftfrac.bin",startyear=1861,stopyear=2005,bands=32,size=4,headersize=0,
#'         getyearstart=2005,getyearstop=2005,ncells=67420)
#' }
#' 
#' @export
readCFT <- function(inFile,startyear,stopyear,bands,size,headersize,getyearstart,getyearstop,ncells=67420){
  if (getyearstop>stopyear){
    stop(paste("unexpected usage: getyearstop (",getyearstop,") larger than stopyear (",stopyear,") -- stopping"))
  }
  if (getyearstart<startyear){
    stop(paste("unexpected usage: getyearstart (",getyearstart,") smaller than startyear (",startyear,") -- stopping"))
  }
  nyears=getyearstop-getyearstart+1
  input <- file(inFile,"rb")
  seek(input, where = headersize+(getyearstart-startyear)*bands*ncells*size, origin="start")
  if (size==2){
    cftfracs <- readBin(input,integer(),n = nyears*ncells*bands, size=size)
  }else if (size==4){
    cftfracs <- readBin(input,double(),n = nyears*ncells*bands, size=size)
  }else{
    print("unknown data size")
  }
  close(input)
  if (nyears==1){
    dim(cftfracs) <- c(ncells,bands)
  }else{
    dim(cftfracs) <- c(ncells,bands,nyears)
  }
  return(cftfracs)
}

#' Read cft LPJmL input
#'
#' Returns a range of years from a cft LPJmL input of dimension c(ncells,bands,nyears).
#'
#' @param inFile character string containing the file to read the data from
#' @param startyear absolute startyear of output
#' @param stopyear absolute stopyear of output
#' @param bands number of cft-bands (32 for standard LPJmL input, 64 for standard input)
#' @param size size of each cell's data (2 or 4 bytes)
#' @param headersize size of header data (defaults to 43 bytes)
#' @param getyearstart start of range to return
#' @param getyearstop end of range to return
#'
#' @return array of cft-fractions c(bands,ncells,nyears) for each of the 67420 gridcells over requested range of years
#
#' @examples
#' \dontrun{
#' readCFTinput(inFile="cftinput.bin",startyear=1901,stopyear=2015,bands=64,size=4,headersize=43,
#'         getyearstart=2005,getyearstop=2005)
#' }
#' 
#' @export
readCFTinput <- function(inFile,startyear,stopyear,bands,size,dtype,headersize=43,getyearstart,getyearstop,ncells=67420){
  if (getyearstop>stopyear){
    stop(paste("unexpected usage: getyearstop (",getyearstop,") larger than stopyear (",stopyear,") -- stopping"))
  }
  if (getyearstart<startyear){
    stop(paste("unexpected usage: getyearstart (",getyearstart,") smaller than startyear (",startyear,") -- stopping"))
  }
  nyears=getyearstop-getyearstart+1
  input <- file(inFile,"rb")
  seek(input, where = (getyearstart-startyear)*bands*ncells*size+headersize, origin="start")
  if (dtype=="integer"){
      cftfracs <- readBin(input,integer(),n = nyears*ncells*bands, size=size)
  }else if (dtype=="double"){
      cftfracs <- readBin(input,double(),n = nyears*ncells*bands, size=size)
  }else{
    print("unknown data type")
  }
  close(input)
  if (nyears==1){
    dim(cftfracs) <- c(bands,ncells)
  }else{
    dim(cftfracs) <- c(bands,ncells,nyears)
  }
  return(cftfracs)
}

#' Write header of LPJmL input file
#' 
#' Writes the header of a binary CLM(2) LPJmL input file.
#' 
#' @param file.out character string containing the name of the file to write the header in  
#' @param headername character string containing the name of the header 
#' @param version CLM Version (usually 2)
#' @param firstyear first year of data set 
#' @param lastyear last year covered in the data set
#' @param bands number of cft-bands (32 or 64 for standard input)
#' @param scalar for cft fracs set to 0.001
#' 
#' @return None
#' 
#' @examples
#' \dontrun{
#' fwriteheader(fwriteheader("cftfracs.clm",headername="LPJLUSE",version=2, firstyear = 2006,lastyear = 2100, bands=64, scalar=0.001)
#' }
#' 
#' @export
fwriteheader <- function(file.out,headername,version=2,firstyear,lastyear,bands,scalar){
  nyears=lastyear-firstyear+1
  writeChar(headername,file.out,eos=NULL)
  writeBin(as.integer(version),file.out,size=4,endian=.Platform$endian)   # CLIMATE VERSION
  writeBin(as.integer(1),file.out,size=4,endian=.Platform$endian)       # ORDER
  writeBin(as.integer(firstyear),file.out,size=4,endian=.Platform$endian) # FIRSTYEAR
  writeBin(as.integer(nyears),file.out,size=4,endian=.Platform$endian)       # NYEAR
  writeBin(as.integer(0),file.out,size=4,endian=.Platform$endian)       # FIRSTCELL
  writeBin(as.integer(ncells),file.out,size=4,endian=.Platform$endian)       # NCELL
  writeBin(as.integer(bands),file.out,size=4,endian=.Platform$endian)       # NBAND
  writeBin(0.5,file.out,size=4,endian=.Platform$endian)               # CELLSIZE
  writeBin(scalar,file.out,size=4,endian=.Platform$endian)               # SCALAR
}

#' Write cft LPJmL input file
#' 
#' Writes a binary CLM(2) LPJmL input file with header from an array of dimension c(ncells,bands,nyears) 
#' and checks whether the sums of all cft fracs per grid cell are always <= 1000. 
#' 
#' @param filename character string containing the name of the file to write the cft input in
#' @param input_array array of dimension c(ncells,bands,nyears)
#' @param firstyear first year in array
#' @param lastyear last year in array 
#' @param bands number of cft-bands (32 or 64 for standard input) 
#' 
#' @return boolean containing TRUE (the sums of all cft fracs per grid cell are always equal or below 1000) or FALSE (the sums of all cft fracs per grid cell are not always equal or below 1000)
#' 
#' @examples
#' \dontrun{
#' writeCFTinput(filename="cftfrac_rcp26_2006-2100_64bands.clm", input_array = CFT_2006_2100, firstyear=2006,lastyear=2100,bands=64)
#' }
#' 
writeCFTinput <- function(filename, input_array, firstyear, lastyear, bands){
  nyears=lastyear-firstyear+1
  #check if sum of lu shares in array always <= 1000
  lutotal<-apply(input_array,c(2,3),sum)
  check<-FALSE
  range<-range(lutotal)
  if (range[2]<=1000){check[1]<-TRUE}
  print("LU share sums always equal to or below 1000?");print(check)
  # open file for binary writing 
  f.out<- file(filename,"wb")   
  # write header
  fwriteheader(f.out,"LPJLUSE",2, firstyear,nyears, bands, 0.001)
  # write data
  for (i in firstyear:lastyear){
    writeBin(as.integer(as.vector(input_array[,,i-firstyear+1])),f.out,size=2)
  } 
  close(f.out)
}

#' Write wateruse LPJmL input file
#' 
#' Writes a binary CLM(2) LPJmL input file with header from an array of dimension c(ncells,bands,nyears) 
#' 
#' @param filename character string containing the name of the file to write the cft input in
#' @param input_array array of dimension c(ncells,bands,nyears) or c(ncells,nyears)
#' @param firstyear first year in array
#' @param lastyear last year in array 
#' @param bands number of cft-bands (32 or 64 for standard input) 
#' @param out_name LPJmL variable to write out e.g. "LPJLUSE" oder "LPJWUSE"
#' @param out_size size of integer to be writte out (2/4)
#' @param out_scaling optional scaling factor to be applied, defaults to 1
#' 
#' @examples
#' \dontrun{
#' writeWUinput(filename="wateruse_wd_2band_2006_2050.clm", input_array = waterinput, firstyear=2006, lastyear=2050, bands=2, out_name="LPJWUSE", out_size=4, out_scaling=1)
#' }
#' 
writeWUinput <- function(filename, input_array, fromyear, toyear, out_bands, out_name, out_size, out_scaling=1){
  # open file for binary writing 
  f.out<- file(filename,"wb")
  # write header
  print(paste("Writing Header to: ",filename, out_name, 2, fromyear, toyear, out_bands, out_scaling))
  
  fwriteheader(file.out = f.out, headername = out_name, version = 2, firstyear = fromyear, lastyear = toyear, bands = out_bands, scalar = out_scaling)
  # write data
  for (i in fromyear:toyear){
    if (out_bands==1){
      writeBin(as.integer(as.vector(input_array[,i-fromyear+1]*1/out_scaling)),f.out,size=out_size)
    }else{
      for (c in 1:ncells){
        writeBin(as.integer(as.vector(input_array[c,,i-fromyear+1]*1/out_scaling)),f.out,size=out_size)
      }
    }
  } 
  close(f.out)
}

### ================ plotting routines ===================
#plot a global LPJmL array (67420 cells) with manual breaks and a custom colorramp (if manual colorramp is not provided, YlGnBu is used) 
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

#plot an array of country data, e.g. obtained from aggregateLPJmLdataCountry()
plotCountryData <- function(data,cowList,file,sty="lin",title="",legendtitle=""){
  ra=range(data,na.rm = T)
  if (sty=="lin"){
    brks=seq(from=ra[1],to = ra[2],length.out = 12)
    palette=RColorBrewer::brewer.pal(11,"Spectral")
  }else if (sty=="log"){
    if (min(data)<0){
      brks=c(-(2^seq(from=log(ra[2],base=2),to = 1,length.out = 6)),2^seq(from=1,to = log(ra[2],base=2),length.out = 6))
      palette=RColorBrewer::brewer.pal(11,"RdBu")
    }else{
      palette=RColorBrewer::brewer.pal(11,"Spectral")
      brks=c(0,2^seq(from=2,to = log(ra[2],base=2),length.out = 11))
    }
  }else{
    print("Style not known, use 'log' or 'lin'")
  }
  ires=2
  legendticks=seq(from=0,to = 100,length.out = length(brks))
  lpjdata=array(NA,dim=67420)
  co=sort(unique(cowList))
  for (c in 1:length(co)){
    lpjdata[which(cowList==co[c])]=data[c]
  }
  png(file, width=800*ires, height=400*ires, units="px",res=400,pointsize = 4)
  ra <- raster::raster(ncols=360*ires, nrows=180*ires)
  ra[raster::cellFromXY(ra,cbind(lon,lat))] <-  lpjdata
  extent <- raster::extent(c(-180, 180, -60, 90))
  par(bty="n",oma=c(0,0,0,0),mar=c(0,0,0,3),xpd=T)
  raster::plot(ra,ext=extent,breaks=brks,col=palette,main="",legend=FALSE,axes=FALSE,maxpixels=360*180*ires*ires)
  title(title,line=-1)
  fields::image.plot(legend.only=TRUE,zlim=range(brks),col = palette,useRaster=FALSE,breaks=legendticks,
                     lab.breaks=round(brks,2),legend.shrink = 0.8,legend.args=list(legendtitle,side=3, font=2, line=1))
  maps::map('world',add=TRUE,res=0, lwd=0.1,ylim=c(-60,90))
  dev.off()
}

#aggregates all data in lpjml format to country level, using the aggregation Method provided (e.g. sum, mean)
aggregateLPJmLdata2Country <- function(input,cowList,aggMethod="sum"){
  clist=sort(unique(cowList))
  dataList=array(0,dim=length(clist))
  for (c in 1:length(dataList)){
    if (aggMethod=="sum"){
      dataList[c]=sum(input[which(cowList==clist[c])])
    }else if (aggMethod=="mean"){
      dataList[c]=mean(input[which(cowList==clist[c])])
    }else{print("Unknown aggregation Method aggMethod")}
  }
  rownames(dataList)=clist
  return(dataList)
}

#plot regional lpj array
plotRegionalW <- function(data,file,title,pow2max,pow2min,ext,colPos="Blues",colNeg="Reds",legendtitle,legYes,onlyPos=F,eps,ires=12,lon,lat,map=T){
  if (onlyPos){
    legendticks <- c(0,2^seq(pow2min,pow2max,1))
    brks <- c(seq(pow2min,pow2max,length.out = length(legendticks)))
    palette <- c("white",colorRampPalette(RColorBrewer::brewer.pal(9,colPos))(length(legendticks)-2))  
  }else{
    legendticks <- c(-(2^seq(pow2max,pow2min,-1)),2^seq(pow2min,pow2max,1))
    brks <- seq(-pow2max,pow2max,length.out = length(legendticks))
    palette <- c(rev(colorRampPalette(RColorBrewer::brewer.pal(9,colNeg))(length(legendticks)/2-1)),"white",colorRampPalette(RColorBrewer::brewer.pal(9,colPos))(length(legendticks)/2-1))
  }
  data[data<legendticks[1]] <- legendticks[1]
  data[data>legendticks[length(legendticks)]] <- legendticks[length(legendticks)]
  if (eps){
    file=strsplit(file,".",fixed=TRUE)[[1]]
    file=paste(c(file[1:(length(file)-1)],"eps"),collapse=".")
    ps.options(family = c("Helvetica"), pointsize = 18)
    postscript(file,horizontal = FALSE, onefile = FALSE, width=22, height=8.5,paper="special")
  }else{
    png(file, width=800*ires, height=400*ires, units="px",res=400,pointsize = 4)
  }
  #print(paste(ires,length(lon),length(lat)))
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
    fields::image.plot(legend.only=TRUE,zlim=c(-pow2max,pow2max),col = palette,useRaster=FALSE,breaks=brks,lab.breaks=round(legendticks,2),legend.shrink = 0.8,
                       legend.args=list(legendtitle,side=3, font=2, line=1))
  }
  if(map){maps::map('world',add=TRUE,res=0, lwd=0.1,ylim=c(-60,90))}
  dev.off()
}

plotNetcdf <- function(data,file,title,brks,palette,legendtitle,lon,lat,map=T){
  png(file, width=1600, height=900, units="px",res=400,pointsize = 4)
  image.plot(x=lon,y=lat,z=data,breaks = brks,col = palette,bty="n",xlab="",ylab="",axes=F)#,legend.args = c(legendtitle))
  title(title)
  if(map){maps::map('world',add=TRUE,res=0, lwd=0.1,ylim=c(-60,90))}
  dev.off()
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
plotGlobalWflex <- function(data,file,title,man=F,mbrk=-1,mpalette=-1,pow2max,pow2min,colPos="Blues",colNeg="Reds",legendtitle,legYes,onlyPos=F,eps,ires=12,lon,lat,map=T,ext=c(-180, 180, -60, 90),logscale=T){
  
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
  }else{
    if (!man){
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
    }else{
      brks=mbrk
      legendticks=mbrk
      palette=mpalette
    }
  }

  data[data<legendticks[1]] <- legendticks[1]
  data[data>legendticks[length(legendticks)]] <- legendticks[length(legendticks)]
  if (max(ext[1:2])>180){
    lon[lon<0]=lon[lon<0]+360
  }
  #print(paste0("range(lon)",range(lon)))
  xmagn=(ext[2]-ext[1])/360
  ymagn=(ext[4]-ext[3])/180
  if (eps){
    file=strsplit(file,".",fixed=TRUE)[[1]]
    file=paste(c(file[1:(length(file)-1)],"eps"),collapse=".")
    ps.options(family = c("Helvetica"), pointsize = 18)
    postscript(file,horizontal = FALSE, onefile = FALSE, width=22, height=8.5,paper="special")
  }else{
    png(file, width=800*ires*xmagn, height=360*ires*ymagn, units="px",res=400,pointsize = 4)
  }
  ra <- raster::raster(ncols=360*ires, nrows=180*ires)
  range <- range(data)
  ra[raster::cellFromXY(ra,cbind(lon,lat))] <-  data
  extent <- raster::extent(ext)
  if (legYes){
    par(bty="n",oma=c(0,0,0,6),mar=c(0,0,0,6),xpd=T)
  }else{
    par(bty="n",oma=c(0,0,0,0),mar=c(0,0,0,0))
  }
  raster::plot(ra,ext=extent,breaks=legendticks,col=palette,main="",legend=FALSE,axes=FALSE,maxpixels=360*180*ires*ires)
  title(title,line=-1,cex.main=4*xmagn*ires/12)
  if (legYes){
    fields::image.plot(legend.only=TRUE,zlim=range(data),col = palette,useRaster=FALSE,breaks=brks,lab.breaks=round(legendticks,2),legend.width = 6*xmagn,
                       legend.args=list(legendtitle,side=3, font=2, line=0),axis.args=list(cex.axis=4*xmagn*ires/12))
  }
  if(map){maps::map('world',add=TRUE,res=0, lwd=0.1,ylim=c(-60,90))}
  dev.off()
}

#' Plot global water stress array with histogram
#'
#' Creates a PNG/eps with a plot of a global water stress array
#'
#' @param data array with data to plot in LPJmL specific array c(67420)
#' @param belocs array with be plantation locations in LPJmL specific array c(67420)
#' @param file character string for location/file to save plot to
#' @param title character string title for plot
#' @param legendtitle character string legend title
#' @param legYes show legend (boolean)
#' @param beplot show small beplot (boolean)
#' @param eps write eps file instead of PNG (boolean)
#'
#' @return None
#'
#' @examples
#' \dontrun{
#' plotStressHist(data=irrigation2006,belocs=bearray[,2100],file=paste("~/","mwateramount_2005_06.png",sep=""),
#'             title = paste("irrigation amount 2006 in mm/yr",sep=""),
#'             legendtitle="legendtitle",legYes=TRUE, beplot=TRUE,eps=FALSE)
#' }
#' @export
plotStressHist <- function(data,pop,belocs,file,title,legendtitle,beplot=F,legYes=F,hist=F,areaLeg=T,eps){
  #legendticks <- c("extreme high","high","medium","low","no","low","medium","high","extreme high")  
  if (min(data,na.rm=T)>=0){
    brks <- c(0,0.1,20,40,70,100)  
    palette <- c("white","gold","coral","red","red4")
  }else{
    brks <- c(-100,-70,-40,-20,-0.1,0.1,20,40,70,100)
    palette <- c("darkblue","dodgerblue3","deepskyblue","cyan","white","gold","coral","red","red4")
  }
  data[!is.finite(data)] <- 0
  data[data<0] <- 0
  data[data>100] <- 100
  if (eps){
    file=strsplit(file,".",fixed=TRUE)[[1]]
    file=paste(c(file[1:(length(file)-1)],"eps"),collapse=".")
    ps.options(family = c("Helvetica"), pointsize = 18)
    postscript(file,horizontal = FALSE, onefile = FALSE, width=22, height=8.5,paper="special")
  }else{
    png(file, width=7.25, height=3.5, units="in", res=400, pointsize=6,type="cairo")
  }

  ra <- raster::raster(ncols=720, nrows=360)
  range <- range(data)
  ra[raster::cellFromXY(ra,cbind(lon,lat))] <-  data
  raBE <- raster::raster(ncols=720, nrows=360)
  raBE[raster::cellFromXY(raBE,cbind(lon,lat))] <-  belocs
  
  extent <- raster::extent(c(-180, 180, -60, 90))
  par(bty="n",oma=c(0,0,0,0),mar=c(0,0,0,0))
  raster::plot(ra,ext=extent,breaks=brks,col=palette,main="",legend=FALSE,axes=FALSE)
  title(title,line=-1)
  maps::map('world',add=TRUE,res=0.4, lwd=0.25,ylim=c(-60,90))
  area=cellarea/10^10#from m2 to Mha
  area0_01=sum(area[data>=0&data<=0.1],na.rm=T)
  area01_20=sum(area[data>0.1&data<=20],na.rm=T)
  area20_40=sum(area[data>20&data<=40],na.rm=T)
  area40_70=sum(area[data>40&data<=70],na.rm=T)
  area70_100=sum(area[data>70&data<=100],na.rm=T)
  ar=round(c(area0_01,area01_20,area20_40,area40_70,area70_100),0)
  
  pop0_01=sum(pop[data>=0&data<=0.1],na.rm=T)
  pop01_20=sum(pop[data>0.1&data<=20],na.rm=T)
  pop20_40=sum(pop[data>20&data<=40],na.rm=T)
  pop40_70=sum(pop[data>40&data<=70],na.rm=T)
  pop70_100=sum(pop[data>70&data<=100],na.rm=T)
  po=round(c(pop0_01,pop01_20,pop20_40,pop40_70,pop70_100),0)
  dev.off()
  
  #legend(s)
  file=strsplit(file,".",fixed=TRUE)[[1]]
  if (eps){
    file=paste(paste(c(file[1:(length(file)-1)]),collapse="."),"_legend.eps",sep="")
    ps.options(family = c("Helvetica"), pointsize = 18)
    postscript(file,horizontal = FALSE, onefile = FALSE, width=5, height=5,paper="special")
  }else{
    file=paste(paste(c(file[1:(length(file)-1)]),collapse="."),"_legend.png",sep="")
    png(file, width=2, height=2, units="in",bg="transparent", res=300, pointsize=6,type="cairo")
  }
  par(oma=c(0,0,0,0),mar=c(0,0,0,0))
  if (legYes){
    #plot 180,360?
    fields::image.plot(legend.only=TRUE,zlim=c(-100,100),col = palette,useRaster=FALSE,breaks=brks,lab.breaks=brks,legend.shrink = 0.8,legend.args=list(legendtitle,side=3, font=2, line=1))
  }else if (areaLeg){
    plot(NA,xlim=c(0,100),ylim=c(0,100),axes=F,type="n",main="",xlab="",ylab="")
    legend(x=0,y=100,paste(ar,"Mha"),pt.bg=palette[1:5],pch=22,title="Area",fill=NULL,bg="white",border=NA,col="black",pt.cex=4,lty=0,lwd=1,cex = 2.5)
  }else if (beplot){
    plot(NA,xlim=c(-100,100),ylim=c(-100,100),axes=F,main="")
    rect(-53,-77,40,-35,col="white",border="black")
    text(x=-10,y=-38,"bioenergy locations",cex=1)
    par(bty="n",fig=c(0.32,0.57,0.00,0.28),new=T,mar=c(0,0,2.5,0))
    raster::plot(raBE,ext=extent,breaks=c(-0.5,0.5,1.5),col=c("gray80","black"),main="",legend=F,axes=F)#is.args=list(tick=F,lwd=0.25,lwd.tick=0, labels=FALSE))
  }else if (hist){
    par(fig=c(0.02,0.28,0.04,0.5),new=T,mar=c(2,2,2,2))
    plot.new()
    if (min(data,na.rm=T)>=0){
      #data[data<0.1]=NA
      hist(data,breaks=seq(0,100,10),main="")
    }else{
      #data[data<0.1&data>-0.1]=NA
      hist(data,breaks=seq(-100,100,20),main="")
    }
  }
  dev.off()
  return(cbind(ar,po))
}


#' Plot global water stress array
#'
#' Creates a PNG/eps with a plot of a global water stress array
#'
#' @param data array with data to plot in LPJmL specific array c(67420)
#' @param file character string for location/file to save plot to
#' @param title character string title for plot
#' @param legendtitle character string legend title
#' @param legYes show legend (boolean)
#' @param eps write eps file instead of PNG (boolean)
#'
#' @return None
#'
#' @examples
#' \dontrun{
#' plotStress(data=irrigation2006,file=paste("~/","mwateramount_2005_06.png",sep=""),
#'             title = paste("irrigation amount 2006 in mm/yr",sep=""),
#'             legendtitle="legendtitle",legYes=TRUE,eps=FALSE)
#' }
#' @export
plotStress <- function(data,file,title,legendtitle,legYes,eps){
  #legendticks <- c("extreme high","high","medium","low","no","low","medium","high","extreme high")  
  brks <- c(-100,-70,-40,-20,-0.1,0.1,20,40,70,100)  
  
  data[data<brks[1]] <- brks[1]
  data[data>brks[length(brks)]] <- brks[length(brks)]
  
  palette <- c("darkblue","dodgerblue3","deepskyblue","cyan","white","gold","coral","red","red4")
  if (eps){
    file=strsplit(file,".",fixed=TRUE)[[1]]
    file=paste(c(file[1:(length(file)-1)],"eps"),collapse=".")
    ps.options(family = c("Helvetica"), pointsize = 18)
    postscript(file,horizontal = FALSE, onefile = FALSE, width=22, height=8.5,paper="special")
  }else{
    png(file, width=7.25, height=3.5, units="in", res=300, pointsize=6,type="cairo")
  }
  ra <- raster::raster(ncols=720, nrows=360)
  range <- range(data)
  ra[raster::cellFromXY(ra,cbind(lon,lat))] <-  data
  extent <- raster::extent(c(-180, 180, -60, 90))
  if (legYes){
    par(bty="n",oma=c(0,0,0,0),mar=c(0,0,0,3),xpd=T)
  }else{
    par(bty="n",oma=c(0,0,0,0),mar=c(0,0,0,0))
  }
  raster::plot(ra,ext=extent,breaks=brks,col=palette,main="",legend=FALSE,axes=FALSE)
  title(title,line=-1)
  if (legYes){
    fields::image.plot(legend.only=TRUE,zlim=c(-100,100),col = palette,useRaster=FALSE,breaks=brks,lab.breaks=brks,legend.shrink = 0.8,legend.args=list(legendtitle,side=3, font=2, line=1))
  }
  maps::map('world',add=TRUE,res=0.4, lwd=0.25,ylim=c(-60,90))
  dev.off()
}

#' Plot global month indices on LPJmL grid, for example for months with highest precipitation
#'
#' Creates a PNG/eps with a plot of a month index from 1:12
#'
#' @param data array with data to plot in LPJmL specific array c(67420)
#' @param file character string for location/file to save plot to
#' @param title character string title for plot
#' @param legendtitle character string legend title
#' @param legYes show legend (boolean)
#' @param eps write eps file instead of PNG (boolean)
#'
#' @return None
#'
#' @examples
#' plotMonth(data=irrigation2006,file=paste("~/","mwateramount_2005_06.png",sep=""),
#'             title = paste("irrigation amount 2006 in mm/yr",sep=""),
#'             legendtitle="legendtitle",legYes=TRUE,eps=FALSE)
#'
#' @export
plotMonth <- function(data,file,title,legendtitle,legYes,eps){
  #legendticks <- c("extreme high","high","medium","low","no","low","medium","high","extreme high")
  brks <- seq(-0.5,12.5,1)  
  
  data[data<0] <- 0
  data[data>12] <- 0 
  
  palette <- c("gray",RColorBrewer::brewer.pal(12,"Paired")[c(1,2,11,3,4,7,8,5,6,12,9,10)])
  if (eps){
    file=strsplit(file,".",fixed=TRUE)[[1]]
    file=paste(c(file[1:(length(file)-1)],"eps"),collapse=".")
    ps.options(family = c("Helvetica"), pointsize = 18)
    postscript(file,horizontal = FALSE, onefile = FALSE, width=22, height=8.5,paper="special")
  }else{
    png(file, width=7.25, height=3.5, units="in", res=300, pointsize=6,type="cairo")
  }
  ra <- raster::raster(ncols=720, nrows=360)
  range <- range(data)
  ra[raster::cellFromXY(ra,cbind(lon,lat))] <-  data
  extent <- raster::extent(c(-180, 180, -60, 90))
  if (legYes){
    par(bty="n",oma=c(0,0,0,0),mar=c(0,0,0,3),xpd=T)
  }else{
    par(bty="n",oma=c(0,0,0,0),mar=c(0,0,0,0))
  }
  raster::plot(ra,ext=extent,breaks=brks,col=palette,main="",legend=FALSE,axes=FALSE)
  title(title,line=-1)
  maps::map('world',add=TRUE,res=0.4, lwd=0.25,ylim=c(-60,90))
  if (legYes){
    legend("right",legend=c("NA","Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec"),fill = palette,title=legendtitle)
  }
  dev.off()
}

#' Plot global month indices on LPJmL grid with transparency for not highly stressed cells
#'
#' Creates a PNG/eps with a plot of a month index from 1:12
#'
#' @param data array with data to plot in LPJmL specific array c(67420)
#' @param stressf array with stresslevel [0,100] of cells, determining the transparency
#' @param file character string for location/file to save plot to
#' @param title character string title for plot
#' @param legendtitle character string legend title
#' @param legYes show legend (boolean)
#' @param eps write eps file instead of PNG (boolean)
#'
#' @return None
#'
#' @examples
#' plotMonth(data=irrigation2006,file=paste("~/","mwateramount_2005_06.png",sep=""),
#'             title = paste("irrigation amount 2006 in mm/yr",sep=""),
#'             legendtitle="legendtitle",legYes=TRUE,eps=FALSE)
#'
#' @export
plotMonthTransp <- function(data,stressf,file,title,legendtitle,legYes,eps){
  #legendticks <- c("extreme high","high","medium","low","no","low","medium","high","extreme high")
  brks <- seq(-0.5,12.5,1)  
  #data[!is.finite(data)]=0
  data[data<0] <- 0
  data[data>12] <- 0
  data20=data
  data20[stressf>=0.2]=NA
  data40=data
  data40[stressf>=0.4|stressf<0.2]=NA
  data70=data
  data70[stressf>=0.7|stressf<0.4]=NA
  data100=data
  data100[stressf<0.7]=NA
  palette <- c("gray",RColorBrewer::brewer.pal(12,"Paired")[c(1,2,11,3,4,7,8,5,6,12,9,10)])
  if (eps){
    file=strsplit(file,".",fixed=TRUE)[[1]]
    file=paste(c(file[1:(length(file)-1)],"eps"),collapse=".")
    ps.options(family = c("Helvetica"), pointsize = 18)
    postscript(file,horizontal = FALSE, onefile = FALSE, width=22, height=8.5,paper="special")
  }else{
    png(file, width=7.25, height=3.5, units="in", res=300, pointsize=6,type="cairo")
  }
  ra <- raster::raster(ncols=720, nrows=360)
  ra20 <- raster::raster(ncols=720, nrows=360)
  ra40 <- raster::raster(ncols=720, nrows=360)
  ra70 <- raster::raster(ncols=720, nrows=360)
  ra100 <- raster::raster(ncols=720, nrows=360)
  range <- range(data)
  ra[raster::cellFromXY(ra,cbind(lon,lat))] <-  data
  ra20[raster::cellFromXY(ra20,cbind(lon,lat))] <-  data20
  ra40[raster::cellFromXY(ra40,cbind(lon,lat))] <-  data40
  ra70[raster::cellFromXY(ra70,cbind(lon,lat))] <-  data70
  ra100[raster::cellFromXY(ra20,cbind(lon,lat))] <-  data100
  
  extent <- raster::extent(c(-180, 180, -60, 90))
  if (legYes){
    par(bty="n",oma=c(0,0,0,0),mar=c(0,0,0,3),xpd=T)
  }else{
    par(bty="n",oma=c(0,0,0,0),mar=c(0,0,0,0))
  }
  #raster::plot(ra,ext=extent,breaks=brks,col=palette,main="",legend=FALSE,axes=FALSE)#,alpha=0.25)
  raster::plot(ra20,ext=extent,breaks=brks,col=palette,main="",legend=FALSE,axes=FALSE,alpha=0.1)
  par(new=T)
  raster::plot(ra40,ext=extent,breaks=brks,col=palette,main="",legend=FALSE,axes=FALSE,alpha=0.3)
  par(new=T)
  raster::plot(ra70,ext=extent,breaks=brks,col=palette,main="",legend=FALSE,axes=FALSE,alpha=0.6)
  par(new=T)
  raster::plot(ra100,ext=extent,breaks=brks,col=palette,main="",legend=FALSE,axes=FALSE,alpha=1)
  
  title(title,line=-1)
  maps::map('world',add=TRUE,res=0.4, lwd=0.25,ylim=c(-60,90))
  if (legYes){
    legend("right",legend=c("NA","Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec"),fill = palette,title=legendtitle)
  }
  dev.off()
}


#' Plot global LPJmL array
#'
#' Creates a PNG/eps with a plot of a global LPJmL array
#'    Data is plotted in range: c(-2^pow2max,-2^-pow2min,0,2^-pow2min,2^pow2max)
#'    colors for pos and neg values can be given, default is Blues for the positive
#'    and Reds for the negative numbers
#'    0-range (from 2^-pow2min to 2^pow2min) is white.
#'    The negatives can be omitted by setting onlyPos=T, in case there are only pos values.  
#'
#' @param data array with data to plot in LPJmL specific array c(67420)
#' @param file character string for location/file to save plot to
#' @param title character string title for plot
#' @param pow2max upper (positive) end of data range to plot (2^pow2max)
#' @param pow2min smallest positive number to be distinguished from 0 (2^-pow2min)
#' @param colPos color palette for the positives
#' @param colNeg color palette for the negatives
#' @param legendtitle character string legend title
#' @param legYes show legend (boolean)
#' @param onlyPos show only positive (boolean)
#' @param eps write eps file instead of PNG (boolean)
#'
#' @return None
#'
#' @examples
#' plotGlobalW(data=irrigation2006,file=paste("~/","mwateramount_2005_06.png",sep=""),
#'             title = paste("irrigation amount 2006 in mm/yr",sep=""),pow2max=15,pow2min=0,
#'             legendtitle="legendtitle",legYes=TRUE,eps=FALSE)
#'
#' @export
plotGlobalW <- function(data,file,title,pow2max,pow2min,colPos="Blues",colNeg="Reds",legendtitle,legYes,onlyPos=F,eps){
  if (onlyPos){
    legendticks <- c(0,2^seq(pow2min,pow2max,1))
    brks <- c(seq(pow2min,pow2max,length.out = length(legendticks)))
    palette <- c("white",colorRampPalette(RColorBrewer::brewer.pal(9,colPos))(length(legendticks)-2))  
  }else{
    legendticks <- c(-(2^seq(pow2max,pow2min,-1)),2^seq(pow2min,pow2max,1))
    brks <- seq(-pow2max,pow2max,length.out = length(legendticks))
    palette <- c(rev(colorRampPalette(RColorBrewer::brewer.pal(9,colNeg))(length(legendticks)/2-1)),"white",colorRampPalette(RColorBrewer::brewer.pal(9,colPos))(length(legendticks)/2-1))
  }
  data[data<legendticks[1]] <- legendticks[1]
  data[data>legendticks[length(legendticks)]] <- legendticks[length(legendticks)]
  if (eps){
    file=strsplit(file,".",fixed=TRUE)[[1]]
    file=paste(c(file[1:(length(file)-1)],"eps"),collapse=".")
    ps.options(family = c("Helvetica"), pointsize = 18)
    postscript(file,horizontal = FALSE, onefile = FALSE, width=22, height=8.5,paper="special")
  }else{
    png(file, width=7.25, height=3.5, units="in", res=300, pointsize=6,type="cairo")
  }
  ra <- raster::raster(ncols=720, nrows=360)
  range <- range(data)
  ra[raster::cellFromXY(ra,cbind(lon,lat))] <-  data
  extent <- raster::extent(c(-180, 180, -60, 90))
  if (legYes){
    par(bty="n",oma=c(0,0,0,0),mar=c(0,0,0,3),xpd=T)
  }else{
    par(bty="n",oma=c(0,0,0,0),mar=c(0,0,0,0))
  }
  raster::plot(ra,ext=extent,breaks=legendticks,col=palette,main="",legend=FALSE,axes=FALSE)
  title(title,line=-1)
  if (legYes){
    fields::image.plot(legend.only=TRUE,zlim=c(-pow2max,pow2max),col = palette,useRaster=FALSE,breaks=brks,lab.breaks=round(legendticks,2),legend.shrink = 0.8,
                       legend.args=list(legendtitle,side=3, font=2, line=1))
  }
  maps::map('world',add=TRUE,res=0.4, lwd=0.25,ylim=c(-60,90))
  dev.off()
}

#' Plot global LPJmL array
#'
#' Creates a PNG/eps with a plot of a global LPJmL array
#'    Data is plotted linearly in range: c(-max,min,0,min,max)
#'    colors for pos and neg values can be given, default is Blues for the positive
#'    and Reds for the negative numbers
#'    0-range (from -min to min) is white.
#'    The negatives can be omitted by setting onlyPos=T, in case there are only pos values.  
#'
#' @param data array with data to plot in LPJmL specific array c(67420)
#' @param file character string for location/file to save plot to
#' @param title character string title for plot
#' @param max upper (positive) end of data range to plot
#' @param min smallest positive number to be distinguished from 0
#' @param colPos color palette for the positives
#' @param colNeg color palette for the negatives
#' @param legendtitle character string legend title
#' @param legYes show legend (boolean)
#' @param onlyPos show only positive (boolean)
#' @param eps write eps file instead of PNG (boolean)
#'
#' @return None
#'
#' @examples
#' plotGlobalWlin(data=irrigation2006,file=paste("~/","mwateramount_2005_06.png",sep=""),
#'             title = paste("irrigation amount 2006 in mm/yr",sep="")max=15,min=0,
#'             legendtitle="legendtitle",legYes=TRUE,eps=FALSE)
#'
#' @export
plotGlobalWlin <- function(data,file,title,max,min,colPos="Blues",colNeg="Reds",legendtitle,legYes,onlyPos=F,eps){
  if (onlyPos){
    lengthbrks <- 16
    if (max-min > 10){
      brks <- round(seq(min,max,length.out = lengthbrks),0)
    }else{
      brks <- round(seq(min,max,length.out = lengthbrks),1)
    }
    palette <- c("white",colorRampPalette(RColorBrewer::brewer.pal(9,colPos))(lengthbrks-2))  
  }else{
    lengthbrks <- 24
    if (max-min > 10){
      brks <- round(seq(-max,max,length.out = lengthbrks),0)
    }else{
      brks <- round(seq(-max,max,length.out = lengthbrks),1)
    }
    palette <- c(rev(colorRampPalette(RColorBrewer::brewer.pal(9,colNeg))(lengthbrks/2-1)),"white",colorRampPalette(RColorBrewer::brewer.pal(9,colPos))(lengthbrks/2-1))
  }
  data[data<brks[1]] <- brks[1]
  data[data>brks[length(brks)]] <- brks[length(brks)]
  if (eps){
    file=strsplit(file,".",fixed=TRUE)[[1]]
    file=paste(c(file[1:(length(file)-1)],"eps"),collapse=".")
    ps.options(family = c("Helvetica"), pointsize = 18)
    postscript(file,horizontal = FALSE, onefile = FALSE, width=22, height=8.5,paper="special")
  }else{
    png(file, width=7.25, height=3.5, units="in", res=300, pointsize=6,type="cairo")
  }
  ra <- raster::raster(ncols=720, nrows=360)
  range <- range(data)
  ra[raster::cellFromXY(ra,cbind(lon,lat))] <-  data
  extent <- raster::extent(c(-180, 180, -60, 90))
  if (legYes){
    par(bty="n",oma=c(0,0,0,0),mar=c(0,0,0,3),xpd=T)
  }else{
    par(bty="n",oma=c(0,0,0,0),mar=c(0,0,0,0))
  }
  raster::plot(ra,ext=extent,breaks=brks,col=palette,main="",legend=FALSE,axes=FALSE)
  title(title,line=-1)
  if (legYes){
    fields::image.plot(legend.only=TRUE,zlim=c(-max,max),col = palette,useRaster=FALSE,breaks=brks,lab.breaks=brks,legend.shrink = 0.8,
                       legend.args=list(legendtitle,side=3, font=2, line=1))
  }
  maps::map('world',add=TRUE,res=0.4, lwd=0.25,ylim=c(-60,90))
  dev.off()
}

#' Plot global LPJmL array with weighted sums of the pos and neg values
#'
#' Creates a PNG/eps with a plot of a global LPJmL array
#'    Data is plotted in range: c(-2^pow2max,-2^-pow2min,0,2^-pow2min,2^pow2max)
#'    colors for pos and neg values can be given, default is Blues for the positive
#'    and Reds for the negative numbers
#'    0-range (from 2^-pow2min to 2^pow2min) is white.
#'    Also globally weighted sums for all pos and neg values are plotted. 
#'    
#' @param data array with data to plot in LPJmL specific array c(67420)
#' @param file character string for location/file to save plot to
#' @param title character string title for plot
#' @param pow2max upper (positive) end of data range to plot (2^pow2max)
#' @param pow2min smallest positive number to be distinguished from 0 (2^-pow2min)
#' @param colPos color palette for the positives
#' @param colNeg color palette for the negatives
#' @param legendtitle character string legend title
#' @param legYes show legend (boolean)
#' @param eps write eps file instead of PNG (boolean)
#'
#' @return None
#'
#' @examples
#' plotGlobalWsum(data=irrigation2006,file=paste("~/","mwateramount_2005_06.png",sep=""),
#'             title = paste("irrigation amount 2006 in mm/yr",sep=""),pow2max=15,pow2min=0,
#'             legendtitle="legendtitle",legYes=TRUE,eps=FALSE)
#'
#' @export
plotGlobalWsum <- function(data,popData=0,file,title,pow2max,pow2min,colPos="Blues",colNeg="Reds",posLab="Pos",negLab="Neg",legendtitle,legYes,eps,pie=T){
  legendticks <- c(-(2^seq(pow2max,pow2min,-1)),2^seq(pow2min,pow2max,1))
  brks <- seq(-pow2max,pow2max,length.out = length(legendticks))
  data[which(!is.finite(data))] <- 0
  if (!length(popData)==1){
    pdata=data*popData
    ptot=sum(pdata[pdata>0])-sum(pdata[pdata<0])
    psump=sum(pdata[pdata>0])/ptot*100
    psumn=-sum(pdata[pdata<0])/ptot*100
    pp=round(cbind(psump,psumn),0)
  }
  gdata=data*cellarea
  gtot=sum(gdata[gdata>0])-sum(gdata[gdata<0])
  gsump=sum(gdata[gdata>0])/gtot*100
  gsumn=-sum(gdata[gdata<0])/gtot*100
  gp=round(cbind(gsump,gsumn),0)
  data[data<legendticks[1]] <- legendticks[1]
  data[data>legendticks[length(legendticks)]] <- legendticks[length(legendticks)]
  palette <- c(rev(colorRampPalette(RColorBrewer::brewer.pal(9,colNeg))(length(legendticks)/2-1)),"white",colorRampPalette(RColorBrewer::brewer.pal(9,colPos))(length(legendticks)/2-1))
  if (eps){
    file=strsplit(file,".",fixed=TRUE)[[1]]
    file=paste(c(file[1:(length(file)-1)],"eps"),collapse=".")
    ps.options(family = c("Helvetica"), pointsize = 18)
    postscript(file,horizontal = FALSE, onefile = FALSE, width=22, height=8.5,paper="special")
  }else{
    png(file, width=7.25, height=3.5, units="in", res=300, pointsize=6,type="cairo")
  }
  ra <- raster::raster(ncols=720, nrows=360)
  range <- range(data)
  ra[raster::cellFromXY(ra,cbind(lon,lat))] <-  data
  extent <- raster::extent(c(-180, 180, -60, 90))
  modLegendTicks=seq(0,length(legendticks)-1,1)
  if (legYes){
    par(bty="n",oma=c(0,0,0,0),mar=c(0,0,0,3),xpd=T)
  }else{
    par(bty="n",oma=c(0,0,0,0),mar=c(0,0,0,0))
  }
  raster::plot(ra,ext=extent,breaks=legendticks,col=palette,main="",legend=FALSE,axes=FALSE)
  title(title,line=-1)
  if (legYes){
    fields::image.plot(legend.only=TRUE,zlim=c(-pow2max,pow2max),col = palette,useRaster=FALSE,breaks=brks,lab.breaks=round(legendticks,2),legend.shrink = 0.8,
                       legend.args=list(legendtitle,side=3, font=2, line=1))
  }
  maps::map('world',add=TRUE,res=0.4, lwd=0.25,ylim=c(-60,90))
  dev.off()
  file=strsplit(file,".",fixed=TRUE)[[1]]
  if (eps){
    file=paste(paste(c(file[1:(length(file)-1)]),collapse="."),"_legend.eps",sep="")
    ps.options(family = c("Helvetica"), pointsize = 18)
    postscript(file,horizontal = FALSE, onefile = FALSE, width=5, height=5,paper="special")
  }else{
    file=paste(paste(c(file[1:(length(file)-1)]),collapse="."),"_legend.png",sep="")
    png(file, width=2.5, height=2.5, units="in",bg = "transparent",res=300, pointsize=6,type="cairo")
  }
  par(oma=c(0,0,0,0),mar=c(4,4,4,4),xpd=T)
  if (pie){
    #par(fig=c(0,0.2,0,0.2),oma=c(0,0,0,0),mar=c(0,0,0,0))
    pie(x = gp,labels = paste(gp,"%",sep=""),col = palette[c(length(palette)-3,4)],cex=2)
  }else{
    plot(0,0,xlim=c(-205,-120),ylim=c(-60,25),xaxt="n",type="n",yaxt="n")
    legend(-200,20,c("Area weighted:",paste(c(paste(posLab,": ",sep=""),paste(negLab,": ",sep="")),gp,"%",sep="")),text.col=c("black",palette[c(length(palette),1)]),pt.cex=3,lty=0,lwd=1,cex = 1.5,bty = "n")
    if (!length(popData)==1){
      legend(-200,-10,c("Population weighted:",paste(c(paste(posLab,": ",sep=""),paste(negLab,": ",sep="")),pp,"%",sep="")),text.col=c("black",palette[c(length(palette),1)]),pt.cex=3,lty=0,lwd=1,cex = 1.5,bty = "n")
    }
  }
  dev.off()
  
}

#' Plot global agreement of GCM based LPJmL array
#'
#' Creates a PNG/eps with a plot of a global agreement in the LPJmL array
#'    Data is plotted in range: c(-2^pow2max,-2^-pow2min,0,2^-pow2min,2^pow2max)
#'    colors for pos and neg values can be given, default is Blues for the positive
#'    and Reds for the negative numbers
#'    0-range (from 2^-pow2min to 2^pow2min) is white.
#'    Also globally weighted sums for all pos and neg values are plotted. 
#'    
#' @param data array with data from the 4 GCMs to plot in LPJmL specific array c(67420,4)
#' @param file character string for location/file to save plot to
#' @param title character string title for plot
#' @param legendtitle character string legend title
#' @param legYes show legend (boolean)
#' @param eps write eps file instead of PNG (boolean)
#'
#' @return None
#'
#' @examples
#' plotGlobalAgree()
#'
#' @export
plotGlobalAgree <- function(data,file,title,legendtitle,legYes,eps){
  occ <- function(vec){
    lv=length(vec[,1])
    most=array(0,dim=c(lv))
    for (i in 1:lv){
      most[i]=max(tabulate(vec[i,]))
    }
    return(most)
  }
  if (eps){
    file=strsplit(file,".",fixed=TRUE)[[1]]
    file=paste(c(file[1:(length(file)-1)],"eps"),collapse=".")
    ps.options(family = c("Helvetica"), pointsize = 18)
    postscript(file,horizontal = FALSE, onefile = FALSE, width=22, height=8.5,paper="special")
  }else{
    png(file, width=7.25, height=3.5, units="in", res=300, pointsize=6,type="cairo")
  }
  data=occ(data)
  brks <- seq(0.5,4.5,1)
  palette <- RColorBrewer::brewer.pal(9,"Set3")[c(4,3,5,7)]
  ra <- raster::raster(ncols=720, nrows=360)
  range <- range(data)
  ra[raster::cellFromXY(ra,cbind(lon,lat))] <-  data
  extent <- raster::extent(c(-180, 180, -60, 90))
  if (legYes){
    par(bty="n",oma=c(0,0,0,0),mar=c(0,0,0,3),xpd=T)
  }else{
    par(bty="n",oma=c(0,0,0,0),mar=c(0,0,0,0))
  }
  raster::plot(ra,ext=extent,breaks=brks,col=palette,main="",legend=FALSE,axes=FALSE)
  title(title,line=-1)
  if (legYes){
    legend(-179,-20,title="Agreeing GCMs",cex=1.3,legend = c("1","2","3","4"),fill=palette,bg="white",border=NA)
  }
  maps::map('world',add=TRUE,res=0.4, lwd=0.25,ylim=c(-60,90))
  dev.off()
}

#' Plot attribution of CC or BE values
#'    
#' @param data scenario array with RCP26 vs RCP 60 water stress diff
#' @param ccArray array with CC attribution 
#' @param beArray array with BE attribution 
#' @param file character string for location/file to save plot to
#' @param title character string title for plot
#' @param pow2max upper (positive) end of data range to plot (2^pow2max)
#' @param pow2min smallest positive number to be distinguished from 0 (2^-pow2min)
#' @param colPos color palette for the positives
#' @param colNeg color palette for the negatives
#' @param legendtitle character string legend title
#' @param legYes show legend (boolean)
#' @param eps write eps file instead of PNG (boolean)
#'
#' @return None
#'
#' @examples
#' plotAttrib(data=waterstressRCP26-waterstressRCP60,ccArray=ccArray,beArray=beArray,file=paste("~/","attribution.png",sep=""),
#'             title = "",pow2max=15,pow2min=0,
#'             legendtitle="legendtitle",legYes=TRUE,eps=FALSE)
#'
#' @export
plotAttrib <- function(data,ccArray,beArray,file,title,pow2max,pow2min,legendtitle,legYes,eps){
  legendticks <- c(0,2^seq(pow2min,pow2max,1))
  brks <- seq(0,pow2max,length.out = length(legendticks))
  
  data[which(!is.finite(data))] <- 0
  ccArray[which(!is.finite(ccArray))] <- 0
  beArray[which(!is.finite(beArray))] <- 0
  
  data_rcp60_cc = abs(data)
  data_rcp60_be = abs(data)
  data_rcp60_na = abs(data)
  data_rcp26_cc = abs(data)
  data_rcp26_be = abs(data)
  data_rcp26_na = abs(data)
  data_rcp60_cc[data>0]=0
  data_rcp60_na[data>0]=0
  data_rcp60_be[data>0]=0
  data_rcp26_cc[data<0]=0
  data_rcp26_na[data<0]=0
  data_rcp26_be[data<0]=0
  ccArray = abs(ccArray)
  beArray = abs(beArray)
  
  data_rcp60_cc[data_rcp60_cc>0 & ccArray<1.2*beArray]=0
  data_rcp60_be[data_rcp60_be>0 & beArray<1.2*ccArray]=0
  data_rcp60_na[data_rcp60_na>0 & (beArray>=1.2*ccArray & ccArray>=1.2*beArray)]=0
  
  data_rcp26_cc[data_rcp26_cc>0 & ccArray<1.2*beArray]=0
  data_rcp26_be[data_rcp26_be>0 & beArray<1.2*ccArray]=0
  data_rcp26_na[data_rcp26_na>0 & (beArray>=1.2*ccArray & ccArray>=1.2*beArray)]=0
  data=cbind(data_rcp60_cc,data_rcp60_be,data_rcp60_na,data_rcp26_cc,data_rcp26_be,data_rcp26_na)
  
  data[data<2^pow2min]=NA
  #View(data)
  #data[data<legendticks[1]] <- legendticks[1]
  #data[data>legendticks[length(legendticks)]] <- legendticks[length(legendticks)]
  
  palette6c <- c("white",colorRampPalette(RColorBrewer::brewer.pal(9,"Blues"))(length(legendticks)-2))
  palette6b <- c("white",colorRampPalette(RColorBrewer::brewer.pal(9,"Greens"))(length(legendticks)-2))
  palette6n <- c("white",colorRampPalette(RColorBrewer::brewer.pal(9,"Purples"))(length(legendticks)-2))
  palette2c <- c("white",colorRampPalette(RColorBrewer::brewer.pal(9,"Reds"))(length(legendticks)-2))
  palette2b <- c("white",colorRampPalette(RColorBrewer::brewer.pal(9,"Oranges"))(length(legendticks)-2))
  palette2n <- c("white",colorRampPalette(RColorBrewer::brewer.pal(9,"Greys"))(length(legendticks)-2))
  palette=cbind(palette6c,palette6b,palette6n,palette2c,palette2b,palette2n)
  if (eps){
    file=strsplit(file,".",fixed=TRUE)[[1]]
    file=paste(c(file[1:(length(file)-1)],"eps"),collapse=".")
    ps.options(family = c("Helvetica"), pointsize = 18)
    postscript(file,horizontal = FALSE, onefile = FALSE, width=22, height=8.5,paper="special")
  }else{
    png(file, width=10, height=3.5, units="in", res=300, pointsize=6,type="cairo")
  }
  write(paste("writing to ",file),stdout())
  range <- range(legendticks)
  extent <- raster::extent(c(-180, 180, -60, 90))
  par(fig=c(0,0.7,0,1), bty="n",oma=c(0,0,0,0),mar=c(0,0,0,0),xpd=T)
  for (i in 1:6){
    ra <- raster::raster(ncols=720, nrows=360)
    ra[raster::cellFromXY(ra,cbind(lon,lat))] <-  data[,i]
    if (i>1){
      par(fig=c(0,0.7,0,1), oma=c(0,0,0,0),new=TRUE,mar = c(0,0,0,0),bty="n",oma=c(0,0,0,0),mar=c(0,0,0,0),xpd=T)
    }
    raster::plot(ra,ext=extent,breaks=legendticks,col=palette[,i],main="",legend=FALSE,axes=FALSE)
  }
  par(fig=c(0,0.7,0,1), oma=c(0,0,0,0),mar = c(0,0,0,0),bty="n",xpd=T,new=TRUE)
  maps::map('world',res=0.4, lwd=0.25,ylim=c(-60,90),add=T)
  title(title,line=-1)
  if (legYes){
    for (i in 1:6){
      #write(paste(palette[,i]),stdout())
      par(fig=c(0.7+(0.3*(i-1)/6),0.7+(0.3*i/6),0,1), new=TRUE,xpd=T)
      #par(fig=c(0.7,0.9,0,1), new=TRUE,xpd=T)
      
      fields::image.plot(legend.only=TRUE,zlim=c(0,pow2max),col = palette[,i],useRaster=FALSE,breaks=brks,lab.breaks=round(legendticks,2),
                         legend.args=list(legendtitle,side=3, font=2, line=1))
    }
  }
  dev.off()
}

#' Plot attribution v2 of CC or BE values
#'    
#' @param data scenario array with RCP26 vs RCP 60 water stress diff
#' @param ccArray array with CC attribution 
#' @param beArray array with BE attribution 
#' @param file character string for location/file to save plot to
#' @param title character string title for plot
#' @param colPos color palette for the positives
#' @param colNeg color palette for the negatives
#' @param legendtitle character string legend title
#' @param legYes show legend (boolean)
#' @param eps write eps file instead of PNG (boolean)
#'
#' @return None
#'
#' @examples
#' plotAttrib2(data=waterstressRCP26-waterstressRCP60,ccArray=ccArray,beArray=beArray,file=paste("~/","attribution.png",sep=""),
#'             title = "",pow2max=15,pow2min=0,
#'             legendtitle="legendtitle",legYes=TRUE,eps=FALSE)
#'
#' @export
plotAttrib2 <- function(data,ccluArray,ccArray,luArray,beArray,popData,file,title,legendtitle,nscen="High climate change",pscen="Bioenergy",legYes=F,areaYes=F,pie=T,eps){
  if (legYes && areaYes){areaYes=F}
  if (pie && areaYes){areaYes=F}
  if (legYes && pie){pie=F}
  
  brks <- seq(0.5,8.5,1)
  data[which(!is.finite(data))] <- 0
  ccArray[which(!is.finite(ccArray))] <- 0
  luArray[which(!is.finite(luArray))] <- 0
  ccluArray[which(!is.finite(ccluArray))] <- 0
  beArray[which(!is.finite(beArray))] <- 0
  
  data_plot = array(NA,dim=c(ncells))
  #sign does not matter, but the absolute difference says sth about the similarity
  ccArray = abs(ccArray) 
  luArray = abs(luArray) 
  ccluArray = abs(ccluArray)
  beArray = abs(beArray)
  
  #stress is higher in rcp6.0 (data<0)
  data_plot[data<0 & ccluArray>1.2*beArray]=1 # higher stressdiff due to climlu-variation 
  data_plot[data_plot==1 & luArray>ccArray]=2
  data_plot[data<0 & beArray>1.2*ccluArray]=3
  data_plot[data<0 & (beArray<=1.2*ccluArray & ccluArray<=1.2*beArray)]=4
  
  #stress is higher in rcp2.6 (data>0)
  data_plot[data>0 & ccluArray>1.2*beArray]=5
  data_plot[data_plot==5 & luArray>ccArray]=6
  data_plot[data>0 & beArray>1.2*ccluArray]=7
  data_plot[data>0 & (beArray<=1.2*ccluArray & ccluArray<=1.2*beArray)]=8
  data_plot[data>-(2^-4) & data<(2^-4)]=NA
  data_area=data*cellarea/10000000
  data_pop=data*popData/100000
  total_area=c(-sum(data_area[data_plot==1],na.rm = T),-sum(data_area[data_plot==2],na.rm = T),
               -sum(data_area[data_plot==3],na.rm = T),-sum(data_area[data_plot==4],na.rm = T),
               sum(data_area[data_plot==5],na.rm = T),sum(data_area[data_plot==6],na.rm = T),
               sum(data_area[data_plot==7],na.rm = T),sum(data_area[data_plot==8],na.rm = T))
  total_pop=c(-sum(data_pop[data_plot==1],na.rm = T),-sum(data_pop[data_plot==2],na.rm = T),
              -sum(data_pop[data_plot==3],na.rm = T),-sum(data_pop[data_plot==4],na.rm = T),
              sum(data_pop[data_plot==5],na.rm = T),sum(data_pop[data_plot==6],na.rm = T),
              sum(data_pop[data_plot==7],na.rm = T),sum(data_pop[data_plot==8],na.rm = T))
  asum=sum(total_area)
  psum=sum(total_pop)
  palette <- c("orangered","orange","coral","hotpink","slateblue","lightseagreen","steelblue1","turquoise1")
  if (eps){
    file=strsplit(file,".",fixed=TRUE)[[1]]
    file=paste(c(file[1:(length(file)-1)],"eps"),collapse=".")
    ps.options(family = c("Helvetica"), pointsize = 18)
    postscript(file,horizontal = FALSE, onefile = FALSE, width=22, height=8.5,paper="special")
  }else{
    png(file, width=7.25, height=3.5, units="in", res=300, pointsize=6,type="cairo")
  }
  range <- range(brks)
  extent <- raster::extent(c(-180, 180, -60, 90))
  par(bty="n",oma=c(0,0,0,0),mar=c(0,0,0,0),xpd=T)
  ra <- raster::raster(ncols=720, nrows=360)
  ra[raster::cellFromXY(ra,cbind(lon,lat))] <- data_plot
  raster::plot(ra,ext=extent,breaks=brks,col=palette,main="",legend=FALSE,axes=FALSE)
  maps::map('world',res=0.4, lwd=0.25,ylim=c(-60,90),add=T)
  title(title,line=-1)
  dev.off()
  
  file=strsplit(file,".",fixed=TRUE)[[1]]
  if (eps){
    file=paste(paste(c(file[1:(length(file)-1)]),collapse="."),"_legend.eps",sep="")
    ps.options(family = c("Helvetica"), pointsize = 18)
    postscript(file,horizontal = FALSE, onefile = FALSE, width=5, height=5,paper="special")
  }else{
    file=paste(paste(c(file[1:(length(file)-1)]),collapse="."),"_legend.png",sep="")
    png(file, width=2.5, height=2.5, units="in",bg = "transparent", res=300, pointsize=6,type="cairo")
  }
  par(oma=c(0,0,0,0),mar=c(4,4,4,4),xpd=T)
  if (pie){
    shrs=round(total_area/asum*100,0)
    pie2(x = shrs,labels = paste(shrs,"%",sep=""),col = palette,cex=2)
  }else{
    plot(NA,xlim=c(-200,-120),ylim=c(-90,40),type="n",axes=F,xlab="",ylab="")
    if (legYes){
      rect(xleft = -184,ybottom = -113,xright = -115,ytop = 52,col = "white")
      legend(-180,30,title=paste("Higher WSI in scenario\n\"",nscen,"\"\nattributable to:",sep=""),cex=1.3,legend = c("Climate change","Landuse","irrigated bioenergy","undetermined"),
             fill=palette[1:4], horiz=F,border=NULL,bty="n",box.col="white",bg="white",ncol=1)
      legend(-180,-52,title=paste("Higher WSI in scenario\n\"",pscen,"\"\nattributable to:",sep=""),cex=1.3,legend = c("Climate change","Landuse","irrigated bioenergy","undetermined"),
             fill=palette[5:8], horiz=F,border=NULL,bty="n",box.col="white",bg="white",ncol=1)
    }else if (areaYes){
      #rect(-179,-43,-90,12,col="white",border="white")
      legend(-176,10,title="Shares area /    ",cex=1.3,legend = paste(round(total_area/asum*100,0)),
             fill=palette, horiz=F,border=NULL,bty="n",box.col="white",ncol=1)
      legend(-172,10,title="                    pop weighted [%]",cex=1.3,legend = paste(round(total_pop/psum*100,0)), bty="n")
    }
  }
  dev.off()
  return(data_plot)
}

#' Plot global LPJmL array inside RStudio
#'
#' Plot of a global LPJmL array inside RStudio
#'    Data is plotted in range: c(-2^pow2max,-2^-pow2min,0,2^-pow2min,2^pow2max)
#'    where the positive values are colored green to blue,
#'    0-range is white,
#'    and the negative ones red to yellow
#'
#' @param data array with data to plot in LPJmL specific array c(67420)
#' @param title character string title for plot
#' @param pow2max upper (positive) end of data range to plot (2^pow2max)
#' @param pow2min smallest positive number to be distinguished from 0 (2^-pow2min)
#' @param legendtitle character string legend title
#' @param legYes show legend (boolean)
#'
#' @return None
#
#' @examples
#' plotGlobalWtoscreen(data=irrigation2006,title = paste("irrigation amount 2006 in mm/yr",sep=""),
#'                     pow2max=15,pow2min=0,"legendtitle",legYes=TRUE)
#'
#' @export
plotGlobalWtoscreen <- function(data,title,pow2max,pow2min,legendtitle,legYes){
  legendticks <- c(-(2^seq(pow2max,-pow2min,-1)),2^seq(-pow2min,pow2max,1))
  brks <- seq(-pow2max,pow2max,length.out = length(legendticks))
  data[data<legendticks[1]] <- legendticks[1]
  data[data>legendticks[length(legendticks)]] <- legendticks[length(legendticks)]
  palette <- c(rev(colorRampPalette(RColorBrewer::brewer.pal(9,"YlOrRd"))(pow2max+pow2min)),"white",colorRampPalette(RColorBrewer::brewer.pal(9,"GnBu"))(pow2max+pow2min))
  ra <- raster::raster(ncols=720, nrows=360)
  range <- range(data)
  ra[raster::cellFromXY(ra,cbind(lon,lat))] <-  data
  extent <- raster::extent(c(-180, 180, -60, 90))
  modLegendTicks=seq(0,length(legendticks)-1,1)
  par(bty="n")
  raster::plot(ra,ext=extent,breaks=legendticks,col=palette,main=title,legend=FALSE,axes=FALSE)
  if (legYes){
    legendtitle=""
    fields::image.plot(legend.only=TRUE,zlim=c(-pow2max,pow2max),col = palette, breaks=brks,lab.breaks=legendticks,
                       legend.args=list(legendtitle,side=4, font=2, line=2.5))
  }
  maps::map('world',add=TRUE,res=0.4, lwd=0.25,ylim=c(-60,90))
}

#' Plot global LPJmL array with only positive values
#'
#' Creates a PNG/eps with a plot of a global LPJmL array
#'    Data is plotted in range: c(0,2^pow2max)
#'    positive values are colored accoding to chosen palette,
#'    0-range is white,
#'
#' @param data array with data to plot in LPJmL specific array c(67420)
#' @param file character string for location/file to save plot to
#' @param title character string title for plot
#' @param pow2max upper (positive) end of data range to plot (2^pow2max)
#' @param col color palette name from brewer.pal (character string)
#' @param legendtitle character string legend title
#' @param legYes show legend (boolean)
#' @param eps write eps file instead of PNG (boolean)
#'
#' @return None
#'
#' @examples
#' \dontrun{
#' plotGlobal_pos(data=irrigation2006,file=paste("~/","mwateramount_2005_06.png",sep=""),
#'             title = paste("irrigation amount 2006 in mm/yr",sep=""),pow2max=15,col="GnBu",
#'             legendtitle="legendtitle",legYes=TRUE,eps=FALSE)
#' }
#' @export
plotGlobal_pos <- function(data,file,title,pow2max,col,legendtitle,legYes,eps){
  legendticks <- c(0,2^seq(0,(pow2max+1),1))
  brks <- seq(0,(pow2max+1),length.out = length(legendticks))
  data[data>legendticks[length(legendticks)]] <- legendticks[length(legendticks)]
  palette <- c("white",colorRampPalette(brewer.pal(9,col))(length(legendticks)-2))
  if (eps){
    file=strsplit(file,".",fixed=T)[[1]]
    file=paste(c(file[1:(length(file)-1)],"eps"),collapse=".")
    ps.options(family = c("Helvetica"), pointsize = 18)
    postscript(file,horizontal = FALSE, onefile = FALSE, width=22, height=8.5,paper="special")
  }else{
    png(file, width=7.25, height=3.5, units="in", res=300, pointsize=6,type="cairo")
  }
  ra <- raster(ncols=720, nrows=360)
  range <- range(data)
  ra[cellFromXY(ra,cbind(lon,lat))] <-  data
  extent <- extent(c(-180, 180, -60, 90))
  if (legYes){
    par(bty="n",oma=c(0,0,0,0),mar=c(0,0,0,3),xpd=T)
  }else{
    par(bty="n",oma=c(0,0,0,0),mar=c(0,0,0,0))
  }
  plot(ra,ext=extent,breaks=legendticks,col=palette,main=title,legend=F,axes=F)
  if (legYes){
    image.plot(legend.only=T,zlim=c(0,pow2max+1),col = palette,useRaster=F,breaks=brks,legend.shrink = 0.8,
               #axis.args=list(at=brks,labels=legendticks),
               legend.args=list(legendtitle,side=3, font=2, line=1),axis.args=list(at=brks[1:(length(brks)-1)],labels=legendticks[1:(length(legendticks)-1)]))
  }
  map('world',add=T,res=0.4, lwd=0.25)
  dev.off()
}

#' Plot global LPJmL array with positive and negative values
#'
#' Creates a PNG/eps with a plot of a global LPJmL array
#'    Data is plotted in range: c(-2^negpow2max,2^pow2max)
#'    positive and negative values are colored accoding to chosen palettes,
#'    0-range is white,
#'
#' @param data array with data to plot in LPJmL specific array c(67420)
#' @param file character string for location/file to save plot to
#' @param title character string title for plot
#' @param negpow2max lower (negative) end of data range to plot (-2^negpow2max)
#' @param pow2max upper (positive) end of data range to plot (2^pow2max)
#' @param poscol color palette name from brewer.pal (character string) for positive values
#' @param negcol color palette name from brewer.pal (character string) for negative values
#' @param legendtitle character string legend title
#' @param legYes show legend (boolean)
#' @param eps write eps file instead of PNG (boolean)
#'
#' @return None
#'
#' @examples
#' \dontrun{
#' plotGlobal_pos_neg(data=EFR_deficits2006,file=paste("~/","EFRdeficits_06.png",sep=""),
#'             title = "",negpow2max=10, pow2max=15,negcol="GnBu",poscol="Reds"
#'             legendtitle="legendtitle",legYes=TRUE,eps=FALSE)
#' }
#' @export
plotGlobal_pos_neg <- function(data,file,title,pow2max,negpow2max,negcol, poscol,legendtitle,legYes,eps){
  legendticks <- c(-(2^seq(negpow2max+1,0,-1)),2^seq(0,pow2max+1,1))
  brks <- seq(-(negpow2max+1),(pow2max+1),length.out = length(legendticks))
  data[data<legendticks[1]] <- legendticks[1]
  data[data>legendticks[length(legendticks)]] <- legendticks[length(legendticks)]
  if (pow2max<negpow2max){
    palette  <- c(colorRampPalette(rev(RColorBrewer::brewer.pal(9,negcol)[1:9]))(length(seq(negpow2max+1,0,-1))-1),"white",colorRampPalette(RColorBrewer::brewer.pal(9,poscol)[1:9])(length(seq((negpow2max+1),0,-1))-1)[1:(pow2max+1)])
  } else if (pow2max>negpow2max){
    palette  <- c(colorRampPalette(rev(RColorBrewer::brewer.pal(9,negcol)[1:9]))(length(seq(0,pow2max+1,1))-1)[(length(seq(0,pow2max+1,1))-(negpow2max+1)):(length(seq(0,pow2max+1,1))-1)],"white",colorRampPalette(RColorBrewer::brewer.pal(9,poscol)[1:9])(length(seq(0,pow2max+1,1))-1))
  } else {
    palette  <- c(colorRampPalette(rev(RColorBrewer::brewer.pal(9,negcol)[1:9]))(length(seq(negpow2max+1,0,-1))-1),"white",colorRampPalette(RColorBrewer::brewer.pal(9,poscol)[1:9])(length(seq(0,pow2max+1,1))-1))
  }
  if (eps){
    file=strsplit(file,".",fixed=T)[[1]]
    file=paste(c(file[1:(length(file)-1)],"eps"),collapse=".")
    ps.options(family = c("Helvetica"), pointsize = 18)
    postscript(file,horizontal = FALSE, onefile = FALSE, width=22, height=8.5,paper="special")
  }else{
    png(file, width=7.25, height=3.5, units="in", res=300, pointsize=6,type="cairo")
  }
  ra <- raster(ncols=720, nrows=360)
  range <- range(data)
  ra[cellFromXY(ra,cbind(lon,lat))] <-  data
  extent <- extent(c(-180, 180, -60, 90))
  if (legYes){
    par(bty="n",oma=c(0,0,0,0),mar=c(0,0,0,3),xpd=T)
  }else{
    par(bty="n",oma=c(0,0,0,0),mar=c(0,0,0,0))
  }
  plot(ra,ext=extent,breaks=legendticks,col=palette,main="",legend=F,axes=F)
  title(title, line = -1)
  if (legYes){
    image.plot(legend.only=T,zlim=c(-(negpow2max+1),pow2max),col = palette,useRaster=F,breaks=brks,legend.shrink = 0.8,
               axis.args=list(at=brks[2:(length(brks)-1)],labels=legendticks[2:(length(legendticks)-1)]),
               legend.args=list(legendtitle,side=3, font=2, line=1))
  }
  map('world',add=T,res=0.4, lwd=0.25)
  dev.off()
}

#' Plot global LPJmL CFT array inside RStudio
#'
#' Plot of one band of a global LPJmL cft array inside RStudio
#'            with classification into 10% classes.
#'
#' @param lushares array with cft band data to plot in LPJmL specific array c(67420)
#' @param title character string title for plot
#' @param legendtitle character string legend title
#'
#' @return None
#
#' @examples
#' \dontrun{
#' plotLUsharesToScreen(lushares=cftfracs2005[,3],title = "cft fracs for maize",legendtitle="")
#' }
#'
#' @export

plotLUsharesToScreen <- function(lushares,title,legendtitle){
  par(mar=c(4, 0, 4, 15) + 0.1,oma=c(0,0,0,2))
  par(cex=1.2,cex.main=1.2)
  cols  <- c(RColorBrewer::brewer.pal(10,"RdYlBu"),"white")
  range <- c(0,1)
  brk   <- c(0,0.001,seq(0,1,0.1)[2:length(seq(0,1,0.1))])
  par(cex=1.2,cex.main=1.2)
  lu.ras <- raster::raster(ncols=720, nrows=360)
  lu.ras[raster::cellFromXY(lu.ras,cbind(lon,lat))] <-  1-lushares
  raster::plot(lu.ras,ylim=c(-60,90),xlim=c(-180,180),zlim=range,breaks=brk,col=(cols),legend=F,xaxt="n",yaxt="n",main=title,axes=F,box=F)
  maps::map('world',add=TRUE,res=0.4, lwd=0.25,ylim=c(-60,90))  #Legend
  par(xpd=T)
  legend(190,70,title=legendtitle,cex=1.6,
         rev(c("0","0-10%","11-20%","21-30%","31-40%","41-50%","51-60%","61-70%","71-80%","81-90%","91-100%")),
         fill=(cols), horiz=F,border=NULL,bty="o",box.col="white",bg="white",ncol=1)
  text(190,70,pos=4,"landuse intensity",cex=1.6)
  par(xpd=F)
}

#' Plot global LPJmL CFT array
#'
#' Plot of one band of a global LPJmL cft array
#'             classification into 10% classes.
#'
#' @param file character string for location/file to save plot to
#' @param lushares array with cft band data to plot in LPJmL specific array c(67420)
#' @param title character string title for plot
#' @param legendtitle character string legend title
#' @param eps write eps file instead of PNG (boolean)
#'
#' @return None
#
#' @examples
#' \dontrun{
#' plotLUshares(file="~/cftfracs2005_maize.png",lushares=cftfracs2005[,3],
#'              title = "cft fracs for maize",legendtitle="",eps=FALSE)
#' }
#'
#' @export
plotLUshares <- function(file,lushares,title,legendtitle,eps){
  if (eps){
    file=strsplit(file,".",fixed=T)[[1]]
    file=paste(c(file[1:(length(file)-1)],"eps"),collapse=".")
    ps.options(family = c("Helvetica"), pointsize = 18) 
    postscript(file,horizontal = FALSE, onefile = FALSE, width=20, height=10,paper="special")
  }else{
    png(file, width=7.2, height=3.6, units="in", res=300, pointsize=6,type="cairo")
  }
  par(mar=c(4, 0, 4, 15) + 0.1,oma=c(0,0,0,2))
  par(cex=1.2,cex.main=1.2)
  cols  <- c(RColorBrewer::brewer.pal(10,"RdYlBu"),"white")
  range <- c(0,1)
  brk   <- c(0,0.001,seq(0,1,0.1)[2:length(seq(0,1,0.1))])
  par(cex=1.2,cex.main=1.2)
  lu.ras <- raster::raster(ncols=720, nrows=360)
  lu.ras[raster::cellFromXY(lu.ras,cbind(lon,lat))] <-  1-lushares
  raster::plot(lu.ras,ylim=c(-60,90),xlim=c(-180,180),zlim=range,breaks=brk,col=(cols),legend=F,xaxt="n",yaxt="n",main=title,axes=F,box=F)
  maps::map('world',add=TRUE,res=0.4, lwd=0.25,ylim=c(-60,90))  #Legend
  par(xpd=T)
  legend(190,70,title=legendtitle,cex=1.6,
         rev(c("0","0-10%","11-20%","21-30%","31-40%","41-50%","51-60%","61-70%","71-80%","81-90%","91-100%")),
         fill=(cols), horiz=F,border=NULL,bty="o",box.col="white",bg="white",ncol=1)
  text(190,70,pos=4,"landuse intensity",cex=1.6)
  par(xpd=F)
  dev.off()
}

#' Plot global LPJmL CFT array
#'
#' Plot of one band of a global LPJmL cft array
#'             classification into 10% classes, except for first class: 0-1%, 1-10%, 
#'             colors according to chosen palette
#'
#' @param file character string for location/file to save plot to
#' @param lushares array with cft band data to plot in LPJmL specific array c(67420)
#' @param title character string title for plot
#' @param col color palette name from brewer.pal (character string)
#' @param legendtitle character string legend title
#' @param eps write eps file instead of PNG (boolean)
#'
#' @return None
#
#' @examples
#' \dontrun{
#' plotLUshares2(file="~/cftfracs2005_maize.png",lushares=cftfracs2005[,3],
#'              title = "cft fracs for maize",col="Reds",legendtitle="",eps=FALSE)
#' }
#'
#' @export

plotLUshares2<- function(file,lushares,title, col, legendtitle, eps){
  if (eps){
    file=strsplit(file,".",fixed=T)[[1]]
    file=paste(c(file[1:(length(file)-1)],"eps"),collapse=".")
    ps.options(family = c("Helvetica"), pointsize = 18) 
    postscript(file,horizontal = FALSE, onefile = FALSE, width=20, height=10,paper="special")
  }else{
    png(file, width=7.2, height=3.6, units="in", res=300, pointsize=6,type="cairo")
  }
  par(cex=1.2,cex.main=1.2)
  cols  <- c("white",RColorBrewer::brewer.pal(9,col)[2],colorRampPalette(RColorBrewer::brewer.pal(9,col)[3:9])(10))
  range <- c(0,1)
  brk   <- c(0,0.0001,0.01,seq(0,1,0.1)[2:length(seq(0,1,0.1))])
  lu.ras <- raster::raster(ncols=720, nrows=360)
  lu.ras[raster::cellFromXY(lu.ras,cbind(lon,lat))] <-  lushares
  extent <- raster::extent(c(-180, 180, -60, 90))
  legend.breaks<- c(0,0.05,seq(0.05,1,0.095)[2:length(seq(0.05,1,0.095))]) # damit auch 0 bis 1 in der Legende sichtbar ist
  legend.cols <-c(RColorBrewer::brewer.pal(9,col)[2],colorRampPalette(RColorBrewer::brewer.pal(9,col)[3:9])(10))
  par(bty="n",oma=c(0,0,0,0),mar=c(0,0,0,3),xpd=T)
  raster::plot(lu.ras,ext=extent,breaks=brk,col=(cols),legend=F,main=title,axes=F,box=F)
  maps::map('world',add=TRUE,res=0.4, lwd=0.25,ylim=c(-60,90))  #Legend
  fields::image.plot(legend.only=TRUE,zlim=range,col = legend.cols,useRaster=FALSE,breaks=legend.breaks,lab.breaks=c("0","1 %","10 %","20 %","30 %","40 %","50 %","60 %","70 %","80 %","90 %","100 %"),legend.shrink = 0.8,
                     legend.args=list(legendtitle,side=3, font=2, line=1))
  dev.off()
}

#' Plot bioenergy harvest over time
#'
#' Plot bioenergy harvest over time
#'
#' @param file character string for location/file to save plot to
#' @param beHarvest array yearly beharvest data -- array c(years,type,clim,rf:irr)
#' @param title character string title for plot
#' @param eps write eps file instead of PNG (boolean)
#'
#' @return None
#
#' @examples
#' \dontrun{
#' plotBEharvest(file="beHarvest.png",beHarvest=beHarvestYearly,
#'              title = "",eps=FALSE)
#' }
#'
#' @export
plotBEharvest <- function(file,beHarv,beHarv_cal,title,expFormat="png"){
  if (expFormat=="eps"){
    file=strsplit(file,".",fixed=TRUE)[[1]]
    file=paste(c(file[1:(length(file)-1)],"eps"),collapse=".")
    ps.options(family = c("Helvetica"), pointsize = 18)
    postscript(file,horizontal = FALSE, onefile = FALSE, width=18, height=18,paper="special")
  }else if (expFormat=="pdf"){  
    file=strsplit(file,".",fixed=TRUE)[[1]]
    file=paste(c(file[1:(length(file)-1)],"pdf"),collapse=".")
    pdf(file,width=14.4,height=7.2,family = c("Helvetica"),pointsize = 12,paper='special',version = "1.5")
  }else{
    png(file, width=7.2, height=3.6, units="in", res=500, pointsize=6,type="cairo")
  }
  par(fig=c(0,0.9,0,1),oma=par("oma")+c(0,0.1,0,0))
  plot(0,type="n", xlim=c(2006,2099), ylim=c(0,11.5),xlab="Year", ylab="harvest [GtC]",main=title,cex.axis=1.5,cex.lab=1.5,xaxs = "i", yaxs = "i")
  grid(ny=NULL,nx=NULL)
  cols=rep("white",12)
  ind=c(5,11,1,12,6,8,18,10,16,13,17,14,15)
  cols[ind]=alpha(c("gray30","orange","orangered","red3","red4","lawngreen","yellowgreen","limegreen","green4","turquoise1","steelblue1","royalblue","royalblue4"),0.7)#RColorBrewer::brewer.pal(10,"Paired")[c(7,8,5,6,3,4,9,1,2,10)])
  ltys=c("solid", "dashed", "dotted", "longdash")#, "dotdash")
  clims=c("HadGEM","MIROC5","GFDL","IPSL")
  scens=c("baseline 0","noefr 15","noefr 30","noefr 45","noefr 60","efr 30","efr 45","efr 60","efr 90","efrwm 30","efrwm 45","efrwm 60","efrwm 90")
  #plot sums
  ytop=11.2
  xpos=seq(2011,2046,length.out = 6)
  rect(xleft = 2006,ybottom = 3.9,xright = xpos[length(xpos)]+2,ytop = ytop+0.3,col = "white")
  text(x=xpos,y=ytop,labels=c("Sum [GtC]",c(clims,"mean")),cex=0.9)
  text(x=xpos,y=ytop-1*0.5,labels=c("ISIMIP2b demand",round(c(colSums(apply(beHarv_cal[,ind[1],1:4,],c(1,2),sum)),mean(colSums(apply(beHarv_cal[,ind[1],1:4,],c(1,2),sum)))),0)),cex=0.9)
  for (i in 1:(length(beHarv[1,,1,1])-5)){
    text(x=xpos,y=ytop-(i+1)*0.5,labels=c(scens[i],round(c(colSums(apply(beHarv[,ind[i],1:4,],c(1,2),sum)),mean(colSums(apply(beHarv[,ind[i],1:4,],c(1,2),sum)))),0)),cex=0.9)
  }
  #plot lines
  for (t in ind[1:(length(beHarv[1,,1,1])-5)]){
    for (c in 1:4){
      lines(x=2006:2099,y=rowSums(beHarv[,t,c,]),col=cols[t],lty=ltys[c])
    }
  }
  colb=alpha("black",0.7)
  for (c in 1:4){
    lines(x=2006:2099,y=rowSums(beHarv_cal[,5,c,]),col=colb,lty=ltys[c])
  }
  legend(2085,5.4,title="scenario",legend = c("ISIMIP2b demand",scens[1:(length(beHarv[1,,1,1])-5)]),col=c(colb,cols[ind]), lty="solid",cex=0.8)
  legend(2067,4,title="climate",legend = clims,col="black", lty=ltys,cex=0.8)
  par(fig=c(0.9,1,0.0,1), mar=c(5.1, 0, 4.1, 0),new=TRUE,xpd=T)
  plot(0,type="n", xlim=c(0,1), ylim=c(0,100),xlab="", ylab="",main="",cex.axis=1.5,cex.lab=1.5,axes=F)
  # plot new prod. increase axis and labels
  ys=c(mean(apply(beHarv_cal[,5,1:4,],c(2),sum)),rowMeans(apply(beHarv[,ind,1:4,],c(2,3),sum)))
  yp=round((ys-ys[2])/(ys[1]-ys[2])*100)
  #lines(x=rep(0.2,2),y=yp[1:2],col="black",lwd=0.7)
  text(x=0.65,y=yp,labels=paste(yp,"%"),adj=1,cex=1)
  text(x=0.85,y=(yp[1]+yp[2])/2,labels="total productivity increase",srt=90,cex=1)
  #plot means
  b=yp[1]
  yp=yp[2:length(yp)]
  for (i in 1:(length(beHarv[1,,1,1])-5)){
    xs=c(0.05,0.25)
    #if (i==11 || i==9){xs=c(0,0.15)}
    #else if(i==3 || i==4){xs=c(0.15,0.3)}
    lines(x=xs,y=rep(yp[i],2),col=cols[ind[i]],lwd=1.5,xaxt="n")
  }
  lines(x=c(0.05,0.25),y=rep(b,2),col=colb,lwd=1.5,xaxt="n")
  #lines(x=c(0,0.2),y=rep(mean(rowSums(beHarv_cal[94,5,,])),2),col=colb,lwd=1.5)
  dev.off()
}

#' Write out latex formatted data table 
#'
#' Write out latex formatted data table for paper
#' layout: scen|BEharvest|withdrawals|maxStressedArea|meanStressedArea
#'
#' @param outfile character string for location/file to save latex table to
#' @param BEdata array with bioenergy data
#' @param WDdata array with bioenergy irrigation water withdrawal data
#' @param maxAREAdata array with stressed area data
#' @param meanAREAdata array with stressed area data
#'
#' @return None
#
#' @examples
#' \dontrun{
#' write_latex_stress_table(file="~/table.tex",BEdata=beDATA,WDdata=waterwdDATA,maxAREAdata=MAXstressedAREAS,meanAREAdata=MEANstressedAREAS)
#' }
#'
#' @export
write_latex_stress_table <- function(outfile,BEdata,WDdata,maxAREAdata,meanAREAdata){
  ofile=file(outfile,"w")
  names=array("",dim=c(length(BEdata)))
  scenarios=c(7,3,5,11,1,12,6,8,18,10,16,13,17,14,15)
  names[scenarios]=c("today", "RCP6.0", "RCP2.6 0\\% (baseline 0)","RCP2.6 15\\%", "RCP2.6 30\\% (noefr 30)", "RCP2.6 45\\%", 
                                                     "RCP2.6 60\\%", "RCP2.6 30\\% EFR", "RCP2.6 45\\% EFR", "RCP2.6 60\\% EFR (efr 60)", "RCP2.6 90\\% EFR", 
                                                     "RCP2.6 30\\% EFR WM", "RCP2.6 45\\% EFR WM (efrwm 45)", "RCP2.6 60\\% EFR WM", "RCP2.6 90\\% EFR WM")
  postline=array("",dim=c(length(BEdata)))
  postline[scenarios]=c("\\hline\n\\\\\n\\hline", "\\hline\n\\\\", "","\\hline", "\\hline", "", "\\\\", "", 
                                                        "\\hline", "\\hline", "\\\\", "\\hline", "\\hline", "","")
  write(paste("Writing LaTeX Table of stress area data to: ",outfile,sep=""),stdout())
  write(paste("\\begin{tabular}{lcccc}\n\\toprule\n",
              "& \\textbf{bioenergy}&\\textbf{total}&\\textbf{Area experiencing}&\\textbf{Area experiencing} \\\\\n", 
"& \\textbf{harvest}&\\textbf{withdrawals}&\\textbf{$>$40\\% water stress in at}&\\textbf{$>$40\\% water stress in} \\\\\n",
"& \\textbf{[GtC]}&\\textbf{[km$^3$/yr]}&\\textbf{least one month [Mha]}&\\textbf{the yearly mean [Mha]}\\\\\n\\\\\n\\hline\n",sep=""),file=ofile)
  for (sc in scenarios){
    write(paste(names[sc]," & ",round(BEdata[sc],0)," & ",round(WDdata[sc],0)," & ",round(maxAREAdata[sc,1],0)," & ",round(maxAREAdata[sc,2]/10^6,0), " & ",round(meanAREAdata[sc,1],0) , " & ",round(meanAREAdata[sc,2]/10^6,0) ,"\\\\\n",postline[sc],sep=""), file = ofile)
  }
  write(paste("\\\\ \n \\bottomrule\n\\end{tabular}\n",sep=""),file=ofile)
  close(ofile)
}


#' Plot global relative differences
#'
#' Plot of relative differences derived from two LPJmL arrays
#'             classification into 10% classes
#'
#' @param file character string for location/file to save plot to
#' @param diff_data array with relative differences in LPJmL specific array c(67420)
#' @param min lower end of data range (between -1.5 and 0, with one digit allowed to the left of the decimal point); one additional break is added to display values lower than min
#' @param max upper end of data range (between 0 and 1.5, with one digit allowed to the left of the decimal point); one additional break is added to display values higher than max
#' @param title character string title for plot
#' @param poscol color palette name from brewer.pal (character string) for poitive values
#' @param negcol color palette name from brewer.pal (character string) for negative values
#' @param legendtitle character string legend title
#'
#' @return None
#
#' @examples
#' \dontrun{
#' plot_relDiff(file="~/difference in soil carbon_2005_2006.png",diff_data=diff_soilc_2005_2006, min=-0.5, max=0.7, negcol="Reds", poscol="Blues",
#'              title = "relativ difference in soil carbon between 2005 and 2006",legendtitle="diff. in %",eps=FALSE)
#' }
#'
#' @export
plot_relDiff <- function(file,diff_data, min, max, negcol, poscol, title,legendtitle){
  png(file, width=7.2, height=3.6, units="in", res=300, pointsize=6,type="cairo")
  par(cex=1.2,cex.main=1.2)
  if (max<abs(min)&&max!=0){
    cols  <- c(colorRampPalette(rev(RColorBrewer::brewer.pal(9,negcol)[3:9]))(length(seq(min-0.1,0,0.1))-1),"white",colorRampPalette(RColorBrewer::brewer.pal(9,poscol)[3:9])(length(seq(min-0.1,0,0.1))-1)[1:((max+0.1)*10)])
    legend.cols<-c(colorRampPalette(rev(RColorBrewer::brewer.pal(9,negcol)[3:9]))(length(seq(min-0.1,0,0.1))-1),colorRampPalette(RColorBrewer::brewer.pal(9,poscol)[3:9])(length(seq(min-0.1,0,0.1))-1)[1:((max+0.1)*10)])
  } else if (max>abs(min)&&min!=0){
    cols  <- c(colorRampPalette(rev(RColorBrewer::brewer.pal(9,negcol)[3:9]))(length(seq(0,max+0.1,0.1))-1)[((length(seq(0,max+0.1,0.1))-1)-(abs(min-0.1)*10)+1):(length(seq(0,max+0.1,0.1))-1)],"white",colorRampPalette(RColorBrewer::brewer.pal(9,poscol)[3:9])(length(seq(0,max+0.1,0.1))-1))
    legend.cols<-c(colorRampPalette(rev(RColorBrewer::brewer.pal(9,negcol)[3:9]))(length(seq(0,max+0.1,0.1))-1)[((length(seq(0,max+0.1,0.1))-1)-(abs(min-0.1)*10)+1):(length(seq(0,max+0.1,0.1))-1)],colorRampPalette(RColorBrewer::brewer.pal(9,poscol)[3:9])(length(seq(0,max+0.1,0.1))-1))
  } else if (max==0){
    cols  <- c(colorRampPalette(rev(RColorBrewer::brewer.pal(9,negcol)[3:9]))(length(seq(min-0.1,0,0.1))-1),"white",RColorBrewer::brewer.pal(3,poscol)[3])
    legend.cols<-c(colorRampPalette(rev(RColorBrewer::brewer.pal(9,negcol)[3:9]))(length(seq(min-0.1,0,0.1))-1),RColorBrewer::brewer.pal(3,poscol)[3])
  } else if (min==0){
    cols  <- c(RColorBrewer::brewer.pal(3,negcol)[3],"white",colorRampPalette(RColorBrewer::brewer.pal(9,poscol)[3:9])(length(seq(0,max+0.1,0.1))-1))
    legend.cols<-c(RColorBrewer::brewer.pal(3,negcol)[3],colorRampPalette(RColorBrewer::brewer.pal(9,poscol)[3:9])(length(seq(0,max+0.1,0.1))-1))
  } else {
    cols  <- c(colorRampPalette(rev(RColorBrewer::brewer.pal(9,negcol)[3:9]))(length(seq(min-0.1,0,0.1))-1),"white",colorRampPalette(RColorBrewer::brewer.pal(9,poscol)[3:9])(length(seq(0,max+0.1,0.1))-1))
    legend.cols<-c(colorRampPalette(rev(RColorBrewer::brewer.pal(9,negcol)[3:9]))(length(seq(min-0.1,0,0.1))-1),colorRampPalette(RColorBrewer::brewer.pal(9,poscol)[3:9])(length(seq(min-0.1,0,0.1))-1))
  }
  range <- c(min-0.1,max+0.1)
  brk   <- c(seq(min-0.1,0,0.1)[1:(length(seq(min-0.1,0,0.1))-1)],-0.0001,0.0001,seq(0,max+0.1,0.1)[2:length(seq(0,max+0.1,0.1))])
  diff_data[diff_data<brk[1]] <- brk[1] #smaller values than min also displayed 
  diff_data[diff_data>brk[length(brk)]] <- brk[length(brk)] #higher values than max also displayed 
  lu.ras <- raster::raster(ncols=720, nrows=360)
  lu.ras[raster::cellFromXY(lu.ras,cbind(lon,lat))] <-  diff_data
  extent <- raster::extent(c(-180, 180, -60, 90))
  if(max==0){
    legendlabel<-c(as.character(c(paste(round(brk*100,0),"%",sep=" ")))[2:(length(seq(min-0.1,0,0.1))-1)],0)
  }
  else if (min==0){
    legendlabel<-c(0, as.character(c(paste(round(brk*100),"%",sep=" ")))[(length(brk)-length(seq(0,max+0.1,0.1))+2):(length(brk)-1)])
  } else {
    legendlabel<-c(as.character(c(paste(round(brk*100,0),"%",sep=" ")))[2:(length(seq(min-0.1,0,0.1))-1)],0, as.character(c(paste(round(brk*100),"%",sep=" ")))[(length(brk)-length(seq(0,max+0.1,0.1))+2):(length(brk)-1)])
  }
  legend.breaks<- c(seq(min-0.1,-0.1,0.1),0,seq(0.1,max+0.1,0.1))
  par(bty="n",oma=c(0,0,0,0),mar=c(0,0,0,3),xpd=T)
  raster::plot(lu.ras,ext=extent,zlim=range,breaks=brk,col=(cols),legend=F,xaxt="n",yaxt="n",main=title,axes=F,box=F)
  maps::map('world',add=TRUE,res=0.4, lwd=0.25,ylim=c(-60,90))  #Legend
  fields::image.plot(legend.only=TRUE,zlim=range,col = legend.cols,useRaster=FALSE,breaks=legend.breaks,legend.shrink = 0.8,
                     legend.args=list(legendtitle,side=3, font=2, line=1),axis.args=list(at=legend.breaks[2:(length(legend.breaks)-1)],labels=legendlabel))
  dev.off()
}

#' Copy of the pie function with labels further away from pie
#'
#' Copy of the pie function with labels further away from pie
#'
#' @export
pie2 <- function (x, labels = names(x), edges = 200, radius = 0.8, clockwise = FALSE, 
          init.angle = if (clockwise) 90 else 0, density = NULL, angle = 45, 
          col = NULL, border = NULL, lty = NULL, main = NULL, ...) 
{
  if (!is.numeric(x) || any(is.na(x) | x < 0)) 
    stop("'x' values must be positive.")
  if (is.null(labels)) 
    labels <- as.character(seq_along(x))
  else labels <- as.graphicsAnnot(labels)
  x <- c(0, cumsum(x)/sum(x))
  dx <- diff(x)
  nx <- length(dx)
  plot.new()
  pin <- par("pin")
  xlim <- ylim <- c(-1, 1)
  if (pin[1L] > pin[2L]) 
    xlim <- (pin[1L]/pin[2L]) * xlim
  else ylim <- (pin[2L]/pin[1L]) * ylim
  dev.hold()
  on.exit(dev.flush())
  plot.window(xlim, ylim, "", asp = 1)
  if (is.null(col)) 
    col <- if (is.null(density)) 
      c("white", "lightblue", "mistyrose", "lightcyan", 
        "lavender", "cornsilk")
  else par("fg")
  if (!is.null(col)) 
    col <- rep_len(col, nx)
  if (!is.null(border)) 
    border <- rep_len(border, nx)
  if (!is.null(lty)) 
    lty <- rep_len(lty, nx)
  angle <- rep(angle, nx)
  if (!is.null(density)) 
    density <- rep_len(density, nx)
  twopi <- if (clockwise) 
    -2 * pi
  else 2 * pi
  t2xy <- function(t) {
    t2p <- twopi * t + init.angle * pi/180
    list(x = radius * cos(t2p), y = radius * sin(t2p))
  }
  for (i in 1L:nx) {
    n <- max(2, floor(edges * dx[i]))
    P <- t2xy(seq.int(x[i], x[i + 1], length.out = n))
    polygon(c(P$x, 0), c(P$y, 0), density = density[i], angle = angle[i], 
            border = border[i], col = col[i], lty = lty[i])
    P <- t2xy(mean(x[i + 0:1]))
    lab <- as.character(labels[i])
    if (!is.na(lab) && nzchar(lab)) {
      lines(c(1, 1.05) * P$x, c(1, 1.05) * P$y)
      text(1.1 * P$x, 1.1 * P$y, labels[i], xpd = TRUE, 
           adj = ifelse(P$x < 0, 1, 0), ...)
    }
  }
  title(main = main, ...)
  invisible(NULL)
}

### ===== end plotting ======

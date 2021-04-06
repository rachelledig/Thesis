# script taken from https://pjbartlein.github.io/REarthSysSci/netCDF.html#map-the-data
#data taken from D’Odorico, P., Chiarelli, D.D., Rosa L., Bini A., Zilberman D., Rulli, M.C. The value of water in agriculture in a global high-resolution analysis. PNAS. 2020

install.packages("ncdf4")
library(ncdf4)
install.packages("chron")
library(chron)
install.packages("lattice")
library(lattice)
install.packages("RColorBrewer")
library(RColorBrewer)
library(raster)
library(sp)


# set path and filename
ncpath <- "/Volumes/RachelExternal/Thesis/Data/"
ncname <- "Current value of water (Figure 5A)"  
ncfname <- paste(ncpath, ncname, ".nc", sep="")
dname <- "current2"  # note: tmp means temperature (not temporary)

# open a netCDF file
ncin <- nc_open(ncfname)
print(ncin)

v <- ncin$var[[1]]
size <- v$varsize
dims <- v$ndims
nt <- size[dims]              # length of time dimension
lat <- ncin$dim$y  # latitude position
lon <- ncin$dim$x 


r<-list()
for (i in 1:nt) {
  start <- rep(1,dims)     # begin with start=(1,1,...,1)
  start[dims] <- i             # change to start=(1,1,...,i) to read    timestep i
  count <- size                # begin with count=(nx,ny,...,nt), reads entire var
  count[dims] <- 1             # change to count=(nx,ny,...,1) to read 1 tstep
  
  dt<-ncvar_get(ncin, varid = 'current2', start = start, count = count)
  
  # convert to raster
  r[i]<-raster(dt)
}

# create layer stack with time dimension
r<-stack(r)

# transpose the raster to have correct orientation
rt<-t(r)
extent(rt)<-extent(c(range(lon), range(lat)))

# plot the result
spplot(rt)

# get longitude and latitude
lon = ncvar_get(ncin,"x")
nlon <- dim(lon)
head(lon)

lat = ncvar_get(ncin,"y")
nlat <- dim(lat)
head(lat)

print(c(nlon,nlat))


# get price
price_array = ncvar_get(ncin,dname)
dlname <- ncatt_get(ncin,dname,"long_name")
dim(price_array)

# get global attributes
Conventions <- ncatt_get(ncin,0,"Conventions")

# quick map
image(lon,lat,price_array, col=rev(brewer.pal(10,"RdBu")))
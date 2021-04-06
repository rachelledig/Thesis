library("sf")
library("ggplot2")
library("rnaturalearth")
#library("rnaturalearthdata")
library("rnaturalearthhires")
library("rgeos")

gplotspath="/Volumes/RachelExternal/Thesis/DataFabian/plots/"

countryGDP <- read.csv(file = paste(file.choose()),stringsAsFactors = FALSE)
names(countryGDP)[names(countryGDP) == "countrycode"] <- "adm0_a3"

world <- ne_countries(scale = "large", returnclass = "sf")
#plot GDP maps
for (y in c(1900:2000)){
  yearGDP=subset(countryGDP,year==y)
  test <- merge(world,yearGDP,by="adm0_a3")
  file=paste(gplotspath,"worldGDP_",y,".png",sep="")
  write(paste("Plotting map to:",file),stdout())
  png(file, width=7.25, height=3.5, units="in", res=400, pointsize=6,type="cairo")
  print(ggplot(data = test) +
    geom_sf(data=test,aes(fill=cgdppc),color="white",lwd=0.05) +
    scale_fill_viridis_c(limits=c(0,150000),trans = "sqrt", alpha = .4,na.value = "white") +
    theme(axis.ticks = element_blank(),
          axis.text = element_blank(),
          axis.line  = element_blank(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_blank(),
          panel.background = element_blank()) +
    coord_sf(crs = st_crs('+proj=moll')) )
  dev.off()
}
#plot area equipped for irrigation 
setClass("num.with.commas")
setAs("character", "num.with.commas", function(from) as.numeric(gsub(",", "", from) ) )
countryAEI <- read.csv(file = paste('/media/Storage/irrigation/HID_v10/Supplement_S3_national_1900-2005.csv',sep=""),colClasses=c(rep('character',3),rep('num.with.commas',14)))
names(countryAEI)[names(countryAEI) == "CODE"] <- "adm0_a3"
decs <- c(seq(1900,1970,10),seq(1980,2005,5))
for (y in decs[14]){
  aei=countryAEI[,c(2,which(decs==y)+3)]
  colnames(aei) <- c("adm0_a3","AEI")
  test <- merge(world,aei,by="adm0_a3")
  file=paste(gplotspath,"worldAEI_",y,".png",sep="")
  write(paste("Plotting map to:",file),stdout())
  png(file, width=7.25, height=3.5, units="in", res=400, pointsize=6,type="cairo")
  print(ggplot(data = test) +
          geom_sf(data=test,aes(fill=AEI),color="white",lwd=0.05) +
          #scale_colour_brewer(limits=c(0,300000000),palette = "GnBu",type = seq,na.value = "white") +
          scale_fill_viridis_c(limits=c(0,100000000),trans="sqrt",alpha = .4,na.value = "white") +
          theme(axis.ticks = element_blank(),
                axis.text = element_blank(),
                axis.line  = element_blank(),
                panel.grid.major = element_blank(),
                panel.grid.minor = element_blank(),
                panel.border = element_blank(),
                panel.background = element_blank()) +
          coord_sf(crs = st_crs('+proj=moll')) 
        )
  dev.off()
}
rm(list=ls(all=TRUE))
Sys.setlocale('LC_ALL','C')
options(scipen = 999)

# ===== Settings =======

ncells=67420
plotting=T

iFol="/Volumes/RachelExternal/Thesis/DataFabian/Data/"
oFol="/Volumes/RachelExternal/Thesis/DataFabian/plots/"
load(file=paste0(iFol,"data.RData"))

#contains:
## distNextIrr [1:67420, 1:105]: distance (in km) for each gridcell to next irrigated cell at that time
## irrFrac: fraction (0..1) of cell that is irrigated
## lat/lon: latitude (-90..90)/longitude (-180..180) of each lpjml gridcell
## lpjGDPpc: GDP per capita (2011 US dollars/cap) for each lpj-grid
## mean_increase: mean increase in productivity (delta gC/m2) in each gridcell across all crops
## median_increase: median increase in productivity (delta gC/m2) in each gridcell across all crops
## popdens [1:67420, 1:111]: population density (cap/km2) starting 1901 to 2011 
# discharge_landuse [1:67420, 1:12, 1:105]: (hm)3/d == (1.000.000 m3)/d
# globalIrrigArea: Mha
# prec_landuse [1:67420, 1:12, 1:105]: monthly precipitation (mm/month)
# runoff_landuse [1:67420, 1:12, 1:105]: amount of water (fed by precipitation) that is not consumed and is handed over to the next gridcell in the riverrouting (mm/month)
# transp_green_landuse [1:67420, 1:12, 1:105]: amount of green water (fed by precipitation) that is transpired by plants (mm/month)
# transp_blue_landuse [1:67420, 1:12, 1:105]: amount of blue water (fed by irrigation) that is transpired by plants (mm/month)
# evap_landuse [1:67420, 1:12, 1:105]: amount of water (fed by precipitation) that evaporates on bare soil (mm/month)
# interc_landuse [1:67420, 1:12, 1:105]: amount of water (fed by precipitation) that evaporates from plants leaves (is intercepted -- never reaches the ground) (mm/month)

load(file=paste0(iFol,"datacftfrac.RData"))

# cftfrac [1:67420, 1:32, 1:105]: fraction of gridcell cultivated with the given cft (crop functional type); cft 1-16 are the rainfed fractions, 17-32 the irrigated ones - the sum of which you have already in irrFrac
# 1 - Temperate cereals (wheat, rye, barley; wheat)
# 2 - Rice (paddy rice; rice)
# 3 - Maize (maize for food; maize)
# 4 - Tropical cereals (millet, sorghum; millet)
# 5 - Pulses (pulses; field peas)
# 6 - Temperate roots (sugar beet; sugar beet)
# 7 - Tropical roots (cassava; cassava)
# 8 - Sunflower (sunflower; sunflower)
# 9 - Soybean (soybean; soybean)
# 10 - Groundnuts (groundnuts; groundnuts)
# 11 - Rapeseed (rapeseed; rapeseed)
# 12 - Sugarcane (sugarcane: sugarcane)
# 13 - others (potatoes, oil palm, citrus, date palm, grapes/vine, cotton, cocoa, coffee, other perennial crops, other annual crops; managed grassland)
# 14 - managed grasslands (pastures; managed grasslands)
# 15 - bio-energy grass
# 16 - bio-energy tree

# https://www.dataquest.io/blog/statistical-learning-for-predictive-modeling-r/
# https://medium.com/@ODSC/hierarchical-bayesian-models-in-r-9a18e6acdf2b

load(file=paste0(iFol,"countryData_67420.RData"))
#contains:
#fullCountryData - data.frame, which contains for every grid cell the LPJmL country id (COW), 
#                  the LPJmL regional code (REG) - both defined in "include/managepar.h",
#                  the more of less official english Countryname, and the international 3 character ISO code
# example usage: countryIrrfrac2005=aggregateLPJmLdata2Country(input=irrFrac[,105]*cellarea,cowList = fullCountryData$ISO,aggMethod = "sum")
#                plotCountryData(data = countryIrrfrac2005/10^6,sty="log",cowList = fullCountryData$ISO,file = paste0(oFol,"countryIrrFrac2005.png"),title = "Country specific total irrigated area",legendtitle = "km2")



# ======== functions ==========
source(paste0(iFol,"lpjmliotools_20201104.R"))


# ========= main ============
if (plotting){
  years=c(2:105)#c(2,10,20,30,40,50,60,70,80,90,100,105)
  for (y in years){
    plotGlobalW(data = mean_increase[,y]*100,file = paste0(oFol,"irrigHarvestIncrease_mean_",y+1900,".png"),title = y+1900,pow2max = 12,pow2min = 0,legendtitle = "",legYes = T,eps = F)
    plotGlobalW(data = median_increase[,y]*100,file = paste0(oFol,"irrigHarvestIncrease_median_",y+1900,".png"),title = y+1900,pow2max = 12,pow2min = 0,legendtitle = "",legYes = T,eps = F)
    plotGlobalMan(data = irrFrac[,y]*100,file = paste0(oFol,"irrigFrac_",y+1900,".png"),title = y+1900,brks=c(-1,0,0.1,1,5,10,20,35,50,75,100),palette=c("#ffffff","#ffebbe","#f5ca7b","#c8d79e","#aacd65","#bdd1ff","#7b8ef5","#444e89","#d89e9d","#8a4444"),legendtitle = "",legYes = T,eps = F)
    plotGlobalW(data = distNextIrr[,y],file = paste0(oFol,"distNextIrr_",y+1900,".png"),title = y+1900,pow2max = 11,pow2min = 4,legendtitle = "",legYes = T,eps = F,onlyPos = T)
    plotGlobalWlin(data = lpjGDPpc[,y]/1000,file = paste0(oFol,"GDPpc_",y+1900,".png"),title = y+1900,max = 50,min = 0,legendtitle = "",legYes = T,eps = F,onlyPos = T)
  }
}

    #make some gifs 

#install.packages("magick")
library(magick)
dir_out <- file.path(oFol)

#lets make some gifs of these
imgs <- list.files(dir_out, pattern = "irrigHarvestIncrease_mean_", full.names = TRUE)
img_list <- lapply(imgs, image_read)

## join the images together
img_joined <- image_join(img_list)

## animate at 2 frames per second
img_animated <- image_animate(img_joined, fps = 2)

## view animated image
img_animated

## save to disk
image_write(image = img_animated,
            path = "meanincreaseyield.gif")

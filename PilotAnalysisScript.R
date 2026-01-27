### Pilot Analysis of Dark Bias presence in LPI and BioTIME. 
# Code written by Hannah Wauchope July-Dec 2025. Last updated 12th Jan 2026.
# See GitHUB repo README for full details. 

#Load Packages
library(plyr)
library(tidyverse)
library(terra)
library(tidyterra)
library(ncdf4)
library(pbapply)
library(ggplot2)
library(pbmcapply)
library(viridis)

DataFP <- "/Users/hannahwauchope/Dropbox/Work/Data/"
options(scipen = 999)

#### SECTION 1 Get and Simplify Biome Data, Get Continent Data ####
## Get Template Raster (for consistent gridcell size etc)
TempRas <- rast("/Users/hannahwauchope/Dropbox/Work/Data/HYDE/AnthromesPost1500/anthromes2000AD.asc") #This is just a template raster to define a spatial grid size (30 arc seconds, ~10km at equator)
values(TempRas) <- 1

### 1A BIOME DATA ###......................................................

## Get WWF Biomes
WWFBiome <- vect(paste0(DataFP, "WWF_Biomes/wwf_terr_ecos.shp"))
BiomeRas <- rasterize(WWFBiome, TempRas, field="BIOME") #Rasterize Biome shapefiles

## Bring Biomes into Broader Groupings:

# Tropical Forests = 1
# 1 = Tropical & Subtropical Moist Broadleaf Forests
# 2 = Tropical & Subtropical Dry Broadleaf Forests
# 3 = Tropical & Subtropical Coniferous Forests
values(BiomeRas)[values(BiomeRas)==1] <- 1
values(BiomeRas)[values(BiomeRas)==2] <- 1
values(BiomeRas)[values(BiomeRas)==3] <- 1


#Temperate Forests = 4
# 4 = Temperate Broadleaf & Mixed Forests
# 5 = Temperate Conifer Forests
values(BiomeRas)[values(BiomeRas)==4] <- 4
values(BiomeRas)[values(BiomeRas)==5] <- 4


#Boreal Forests/Taiga = 6
#Stays same

# Grasslands & Shrublands = 7
# 7 = Tropical & Subtropical Grasslands, Savannas & Shrublands
# 8 = Temperate Grasslands, Savannas & Shrublands
# 9 = Flooded Grasslands & Savannas
# 10 = Montane Grasslands & Shrublands
values(BiomeRas)[values(BiomeRas)==7] <- 7
values(BiomeRas)[values(BiomeRas)==8] <- 7
values(BiomeRas)[values(BiomeRas)==9] <- 7
values(BiomeRas)[values(BiomeRas)==10] <- 7


# Tundra = 11
#Stays same

#Other Forests = 12
# 12 = Mediterranean Forests, Woodlands & Scrub
# 14 = Mangroves
values(BiomeRas)[values(BiomeRas)==12] <- 12
values(BiomeRas)[values(BiomeRas)==14] <- 12


# Deserts & Xeric Shrublands = 13
#Stays same

# Ice/NA = 0 
values(BiomeRas)[values(BiomeRas)==99] <- 0
values(BiomeRas)[values(BiomeRas)==98] <- 0

## 1B CONTINENT DATA ##......................................................
#Load continent shapefile
Continents <- vect("/Users/hannahwauchope/Dropbox/Work/Data/GISLayers/CountryContinentWorldSHPS/World/TM_WORLD_BORDERS-0.3-SSudan.shp")

#Rasterize
ContinentRas <- rasterize(Continents, TempRas, field="REGION", touches=TRUE)

#Add South America ("REGION" in above line does not distinguish between North and South Americas)
SubContinentRas <- rasterize(Continents, TempRas, field="SUBREGION", touches=TRUE)
values(ContinentRas) <- ifelse(values(SubContinentRas)==5, 12, values(ContinentRas))

BiomeRas <- project(BiomeRas, ContinentRas)

#### SECTION 2 Get LPI and BioTIME Data ####
## 2A LPI DATA ##......................................................
LPI <- read.csv("/Users/hannahwauchope/Dropbox/Work/Data/LPI/lpi_pops_20191018.csv")

#Clean LPI Data, convert to long format, group Biomes as per above (Biome data already extracted from WWF shapefile by LPI)
LPI_Cleaned <- LPI %>% 
  select(-contains(c("__"))) %>% 
  rename_with(~ gsub("X", "", .x, fixed = TRUE)) %>% 
  gather(Year, Count, contains(c("19", "20"))) %>% 
  filter(Specific_location==1, Count!="NULL") %>% 
  mutate(Group=case_when(Class %in% c("Actinopterygii","Cephalaspidomorphi", "Elasmobranchii",
                                      "Amphibia", "Sarcopterygii", "Myxini", "Holocephali") ~ "Fish",
                         Class == "Mammalia" ~ "Mammals",
                         Class == "Aves" ~ "Birds",
                         Class == "Reptilia" ~ "Reptiles"))

## 2B CONTINENT DATA ##......................................................
BT <- data.table::fread(paste0(DataFP, "BioTIME/BioTIMEQuery_24_06_2021.csv"))

#Subset to just samples, and post 1950 to align with LPI/habitat data
BT_Samps <- unique(BT[,c("STUDY_ID", "PLOT", "YEAR", "LATITUDE", "LONGITUDE")]) %>% 
  rename(Year = YEAR, Study=STUDY_ID, Plot=PLOT, Latitude=LATITUDE, Longitude=LONGITUDE)
BT_Samps <- subset(BT_Samps, Year>1949)

#### SECTION 3 Get Habitat Loss Data ####
## 3A NETCDF EXTRACTION FUNCTIONS ##......................................................

#Function for extracting data from nc files
GetNCDat <- function(ncobject, VarNames, years){
  pbmclapply(VarNames, function(x){
    VarDat <- ncdf4::ncvar_get(ncobject, x)
    VarDat <- VarDat[,,years]
    return(VarDat)
  }, mc.cores=6)
}

#Function for converting this data into rasters by year, summing across all sub rasters (i.e. all categories of primary habitat loss)
GetNCRasts <- function(VarDat, VarNames, yearsConv){
  AllYearsRas <- pblapply(1:length(yearsConv), function(yr){
    AllVars <- lapply(1:length(VarDat), function(Var){
      VarYear <- VarDat[[Var]][,,yr]
      VarYear <- as.vector(VarYear) 
      VarYear_df <- data.frame(cbind(lonlat,VarYear)) #get biodata
      names(VarYear_df) <- c("lon","lat", VarNames[[Var]])
      Var_ras <- rast(VarYear_df)  #Convert first two columns as lon-lat and third as value
      return(Var_ras)
    })
    AllVars <- sum(rast(AllVars))
    names(AllVars) <- yearsConv[[yr]]
    return(AllVars)
  })
  AllYearsRas <- rast(AllYearsRas)
  return(AllYearsRas)
}

#Get the year data for the rasters
env_nc      <- ncdf4::nc_open(paste0(DataFP, "LUH/transitions.nc")) #open the netcdf
years <- ncdf4::ncvar_get(env_nc, "time") #get the avaliable years
years <- years+1
years <- years[(years+849)>1949]
yearsConv <- years+849

## 3B EXTRACT TRANSITION DATA ##......................................................
#Get transition data. I'm using the first year (i.e. 2014 = 2014-2015 transition)
if(!file.exists(paste0(paste0(DataFP, "LUH/SummedPrimaryHabitatLoss1950to2014.tif")))){
  env_nc      <- ncdf4::nc_open(paste0(DataFP, "LUH/transitions.nc")) #open the netcdf
  names(env_nc$var)
  longitude   <- ncdf4::ncvar_get(env_nc, "lon") #get coordinate data
  latitude    <- ncdf4::ncvar_get(env_nc, "lat") #get coordinate data
  lonlat <- as.matrix(expand.grid(longitude,latitude)) #Make a matrix of coordinates
  
  LossPrimNames <-  names(env_nc$var)[1:18] #1-18
  
  LossPrim <- GetNCDat(env_nc, LossPrimNames, years) #extract primary habitat loss data for all years
  
  AllPrimLoss <- GetNCRasts(LossPrim, LossPrimNames, yearsConv)
  
  writeRaster(AllPrimLoss, file=paste0(DataFP, "LUH/SummedPrimaryHabitatLoss1950to2014.tif"))
} else {
  AllPrimLoss <- rast(paste0(DataFP, "LUH/SummedPrimaryHabitatLoss1950to2014.tif"))
}

#### SECTION 4 Classify LPI and BioTIME by Biome and Continent ####
## 4A LPI CLASSIFY ##......................................................

#Biome Data already extracted for LPI - reclassify
LPI_Cleaned <- LPI_Cleaned %>%
  mutate(Biome = case_when(T_biome=="Tropical and subtropical moist broadleaf forests" ~ "Tropical Forest",
                  T_biome=="Tropical and subtropical dry broadleaf forests" ~ "Tropical Forest",
                  T_biome=="Tropical and subtropical coniferous forests" ~ "Tropical Forest",
                  T_biome=="Temperate coniferous forests" ~ "Temperate Forest",
                  T_biome=="Temperate broadleaf and mixed forests" ~ "Temperate Forest",
                  T_biome=="Boreal forests/taiga" ~ "Boreal Forest",
                  T_biome=="Deserts and xeric shrublands"~"Desert",
                  T_biome=="Tundra" ~ "Tundra",
                  T_biome=="Tropical and subtropical grasslands, savannas and shrublands" ~ "Grass & Shrubland",
                  T_biome=="Temperate grasslands, savannas and shrublands" ~ "Grass & Shrubland",
                  T_biome=="Flooded grasslands and savannas" ~ "Grass & Shrubland",
                  T_biome=="Montane grasslands and shrublands" ~ "Grass & Shrubland",
                  T_biome=="Mediterranean forests, woodlands and scrub" ~ "Other Forest",
                  T_biome=="Mangroves" ~ "Other Forest"))

#Convert into a vector
LPIVect <- vect(LPI_Cleaned, geom=c("Longitude", "Latitude"))

#Extract Continents
LPIVect <- cbind(LPIVect, extract(ContinentRas, LPIVect)[,2]) %>% 
  rename(Continent = y) %>% 
  mutate(Continent = case_when(Continent == 0 ~ "Antarctica",
                               Continent == 2 ~ "Africa",
                               Continent == 9 ~ "Oceania",
                               Continent == 12 ~ "South America",
                               Continent == 19 ~ "North America",
                               Continent == 142 ~ "Asia",
                               Continent == 150 ~ "Europe"))

## 4B BIOTIME CLASSIFY ##......................................................

#Convert into a vector
BTVect <- vect(BT_Samps, geom=c("Longitude", "Latitude"))

#Extract WWF Biomes
BTVect <- cbind(BTVect, extract(BiomeRas, BTVect)[,2])
BTVect <- BTVect %>% 
  rename(Biome = y) %>% 
  mutate(Biome = case_when(Biome==0 ~ "Ice",
                           Biome==1 ~ "Tropical Forest",
                           Biome==4 ~ "Temperate Forest",
                           Biome==6 ~ "Boreal Forest",
                           Biome==7 ~ "Grass & Shrubland",
                           Biome==11 ~ "Tundra",
                           Biome==12 ~ "Other Forest",
                           Biome==13 ~ "Desert"))

#Extract Continents
BTVect <- cbind(BTVect, extract(ContinentRas, BTVect)[,2]) %>% 
  rename(Continent = y) %>% 
  mutate(Continent = case_when(Continent == 0 ~ "Antarctica",
                               Continent == 2 ~ "Africa",
                               Continent == 9 ~ "Oceania",
                               Continent == 12 ~ "South America",
                               Continent == 19 ~ "North America",
                               Continent == 142 ~ "Asia",
                               Continent == 150 ~ "Europe"))



#### SECTION 5 Get Distributions of Habitat loss by Dataset, Continent and Biome ####

## 5A CREATE BIOME/CONTINENT MASKS ##......................................................
# Create raster masks of each Biome and each Continent, so that we can extract habitat loss data by each

#Biome
BiomeMaskNames <- c("TropicalForest", "TemperateForest", "BorealForest", "Grasslands")
BiomeMasks <- lapply(c(1,4,6,7), function(x){
  Mask <- BiomeRas 
  values(Mask) <- ifelse(values(Mask)==x,1,NA)
  Mask <- resample(Mask, AllPrimLoss[[1]])
  crs(Mask) <- crs(AllPrimLoss[[1]])
  return(Mask)
})
names(BiomeMasks) <- BiomeMaskNames

#Continent
ContinentMaskNames <- c("Antarctica", "Africa", "Oceania", "South America", "North America", "Asia", "Europe")
ContinentMasks <- lapply(c(0,2,9,12,19,142,150), function(x){
  Mask <- ContinentRas 
  values(Mask) <- ifelse(values(Mask)==x,1,NA)
  Mask <- resample(Mask, AllPrimLoss[[1]])
  crs(Mask) <- crs(AllPrimLoss[[1]])
  return(Mask)
})
names(ContinentMasks) <- ContinentMaskNames

## 5B HABITAT LOSS DISTRIBUTIONS ACROSS ALL GRID CELLS ##......................................................

#This loops through each year, and extracts the distribution of habitat loss in each biome x continent combination. Loss distributions are binned into 100,000 bins (i.e. # grid cells per bin) to reduce data size. 
GetPropsOfLostPrimary <- pbmclapply(as.character(yearsConv[1:length(yearsConv)-1]), function(yr){
  SetOfLoss <- AllPrimLoss[[yr]]
  LossMasks <- lapply(BiomeMaskNames, function(Mask){
    CLossMasks <- lapply(ContinentMaskNames, function(CMask){
      LossMask <- SetOfLoss*BiomeMasks[[Mask]]
      LossMask <- LossMask*ContinentMasks[[CMask]]
      LossDat <- data.frame("PropLoss" = values(LossMask))
      names(LossDat) <- "PropLoss"
      LossDat <- LossDat %>% 
        drop_na() %>% 
        group_by(group = cut(PropLoss, breaks = seq(0, 1, length.out=100000))) %>% #, Prot
        summarise(NCells = n()) %>% 
        mutate(PropLoss = gsub("[]]", "", str_split_fixed(group, ",", 2)[,2]),
               PropLoss = ifelse(is.na(PropLoss), 0, PropLoss),
               Year = yr,
               Biome = Mask,
               Continent = CMask) %>% 
        select(-group)
      if(nrow(LossDat)>1) {return(LossDat)} else {return(NULL)}
    })
    return(bind_rows(CLossMasks))
  })
  return(bind_rows(LossMasks))
}, mc.cores=6)

GetPropsOfLostPrimary <- bind_rows(GetPropsOfLostPrimary)
GetPropsOfLostPrimary$Data <- "World"

## 5C HABITAT LOSS DISTRIBUTIONS IN LPI DATA ##......................................................

#This loop operates as the one above, but first extracts habitat loss for every LPI sample, and then groups these into bins for comparability with the world data
GetLPIPropsofLostPrimary <- pbmclapply(as.character(yearsConv[1:length(yearsConv)-1]), function(yr){
  SetOfLoss <- AllPrimLoss[[yr]]
  LPIDat <- LPIVect[LPIVect$Year==as.numeric(yr)+1,]
  LossMasks <- lapply(BiomeMaskNames, function(Mask){
    LossMask <- SetOfLoss*BiomeMasks[[Mask]]
    LPIExtract <- cbind(LPIDat, extract(LossMask, LPIDat)[,2])
    names(LPIExtract)[names(LPIExtract)=="y"] <- "PropLoss"
    LPIExtract <- as.data.frame(LPIExtract) %>% 
      drop_na(PropLoss) %>% 
      group_by(group = cut(PropLoss, breaks = seq(0, 1, length.out=100000)), Continent) %>% 
      summarise(NCells = n()) %>% 
      mutate(PropLoss = as.numeric(gsub("[]]", "", str_split_fixed(group, ",", 2)[,2])),
             PropLoss = ifelse(is.na(PropLoss), 0, PropLoss),
             Year = yr,
             Biome = Mask) %>% 
      select(-group)
    return(LPIExtract)
  })
  return(bind_rows(LossMasks))
}, mc.cores=6)

#Clean
GetLPIPropsofLostPrimary <- bind_rows(GetLPIPropsofLostPrimary) %>%
  ungroup() %>% 
  select(Continent, NCells, PropLoss, Year, Biome) %>% 
  mutate(Data="LPI")

## 5D HABITAT LOSS DISTRIBUTIONS IN BIOTIME DATA ##......................................................

#As before, but BioTIME
GetBTPropsofLostPrimary <- pbmclapply(as.character(yearsConv[1:length(yearsConv)-1]), function(yr){
  SetOfLoss <- AllPrimLoss[[yr]]
  BTDat <- BTVect[BTVect$Year==as.numeric(yr)+1,]
  LossMasks <- lapply(BiomeMaskNames, function(Mask){
    LossMask <- SetOfLoss*BiomeMasks[[Mask]]
    BTExtract <- cbind(BTDat, extract(LossMask, BTDat)[,2])
    names(BTExtract)[names(BTExtract)=="y"] <- "PropLoss"
    BTExtract <- as.data.frame(BTExtract) %>% 
      drop_na(PropLoss) %>% 
      group_by(group = cut(PropLoss, breaks = seq(0, 1, length.out=100000)), Continent) %>% 
      summarise(NCells = n()) %>% 
      mutate(PropLoss = as.numeric(gsub("[]]", "", str_split_fixed(group, ",", 2)[,2])),
             PropLoss = ifelse(is.na(PropLoss), 0, PropLoss),
             Year = yr,
             Biome = Mask) %>% 
      select(-group)
    return(BTExtract)
  })
  return(bind_rows(LossMasks))
}, mc.cores=6)
GetBTPropsofLostPrimary <- bind_rows(GetBTPropsofLostPrimary)
GetBTPropsofLostPrimary$Data <- "BioTIME"

## 5E BRING TOGETHER DISTRIBUTION ESTIMATES ##......................................................
GetPropsOfLostPrimary$PropLoss <- as.numeric(GetPropsOfLostPrimary$PropLoss)
AllProp <- data.table::rbindlist(list(GetPropsOfLostPrimary, GetLPIPropsofLostPrimary,
                                      GetBTPropsofLostPrimary), use.names=TRUE, fill=TRUE) %>%
  select(-group) %>%
  drop_na(Continent) %>% 
  filter(!Continent %in% c("Antarctica")) %>% 
  group_by(Data, Biome, Continent) %>% 
  mutate(NCellsWeight = NCells/sum(NCells),
         PropLoss=PropLoss+0.00001)

#### SECTION 6 Plot results ####
#Create a dataframe of all Biome x Continent Combinations
AllPlots <- expand.grid(unique(AllProp$Biome), unique(AllProp$Continent))
names(AllPlots) <- c("Biome", "Continent")
AllPlots$SaveName <- paste0(AllPlots$Biome, "_", AllPlots$Continent)

#Loop through each Biome x Continent combination and create density plot comparing all grid cells with those sampled by LPI/BioTIME
lapply(1:nrow(AllPlots), function(x){
  DensDat <- AllProp %>% filter(Biome == AllPlots[x,1],
                                Continent == AllPlots[x,2])
  if(nrow(DensDat)==0){
    return(NULL)
  } #Skip if no data for this combo (e.g. African Boreal Forest)
  (DensityPlot <- ggplot(DensDat, aes(x=as.numeric(PropLoss), weight=NCellsWeight, group=Data, colour=Data, fill=Data))+ #weight=NCellsWeight, 
      geom_density(alpha=0.6, bw=0.2, adjust=2)+
      theme_classic()+
      ylab("Frequency")+
      scale_fill_manual(values=c("BioTIME"="#D7A305","LPI"="#589000","World"= "grey70"))+
      scale_colour_manual(values=c("BioTIME"="#D7A305","LPI"="#589000","World"= "grey70"))+
      scale_y_continuous(expand=c(0,0))+
      scale_x_continuous(trans='log10', limits=c(0.00001,1), breaks=c(0.00001, 0.85), labels=c(0, 100), name="Log % Lost Primary Habitat", expand=c(0,0))+
      theme(axis.text.y = element_blank(), axis.ticks.y= element_blank(), text = element_text(family="Gill Sans", size=20),
            legend.position="none", aspect.ratio = 0.55))
  ggsave(paste0("Figures/DensityPlot_", AllPlots[x,3], ".png"), DensityPlot, units="cm", height=8, width=14)
})



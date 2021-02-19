#Lukas DeFilippo 12/12/19, last updated 4/24/20
#Coho forecast script - read in data from csv, subset and re-structure as necessary, fit forecast model

#Load required packages####
library(reshape)
library(rstan)
library(shinystan)
library(sp)
library(rgdal)
library(car)
library(gstat)
library(tidyr)
library(dvmisc)
library(geoR)
library(boot)
library(invgamma)
#Mapping libraries
library("ggplot2")
theme_set(theme_bw())
library("sf")
library("rnaturalearth")
library("rnaturalearthdata")
library(mapdata)
library(lwgeom)
library(maptools)
library(ggrepel)
library(msm)
library(RColorBrewer)
library(emdbook)
library(mvtnorm)
library(bayesplot)
library(corrplot)
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())
set.seed(666)
#####


#Load functions####
#Function to calculate 95% quantiles
cred <- function(x){
  return(quantile(x, c(0.025, 0.975)))
  
  #Function to calculate 50% quantiles
}
cred.50 <- function(x){
  return(quantile(x, c(0.25, 0.75)))
}

cred.80 <- function(x){
  return(quantile(x, c(0.1, 0.9)))
}

#calculate 95% quantiles for columns
cred.2 <- function(x){
  return(apply( x, 2, cred))
}

#calculate 50% quantiles for columns
cred.3 <- function(x){
  return(apply( x, 2, cred.50))
}

#calculate 50% quantiles for columns
cred.4 <- function(x){
  return(apply( x, 2, cred.80))
}
#Function to calculate the medians of a column
colMedian <- function(x){
  return(apply(x, 2, median))
}

#Posterior plotting function for comparing simulated data to model estimates for single parameters
post_plot_func <- function(par, par_real, name){
  plot(median(par),pch=3, cex=1, col=rgb(238, 64, 0, max = 255, alpha = 200), xaxt='n', ylab=NA, xlab=NA, main=name,
       ylim=c(quantile(par, 0.01), quantile(par, 0.99)))
  arrows((1), quantile(par, 0.25), 1, quantile(par, 0.75), angle=90, code=3, length=0.0, lwd=5, lty=1, col=rgb(238, 64, 0, max = 255, alpha = 150))
  arrows((1), quantile(par, 0.025), 1, quantile(par, 0.975), angle=90, code=3, length=0.025, lwd=1, lty=1, col=rgb(238, 64, 0, max = 255, alpha = 150))
  points(par_real, col='black', cex=1, pch=16)
  #print(par_real > quantile(par, c(0.25, 0.75))[1] & par_real < quantile(par, c(0.25, 0.75))[2])
  #print(par_real > quantile(par, c(0.025, 0.975))[1] & par_real < quantile(par, c(0.025, 0.975))[2])
}

#Function for plotting exponential kernel
exp_func <- function(d, rho_dist, alpha_dist){ 
  dc <- exp(-(d^2)/(2*rho_dist^2))
  return(dc)
}
######

#Clean and prepare the data
dat_func <- function(coho_file='Coho data_3-30-20.csv', cwt_file='CWT_FRAM_Matches_complete20200303.csv', stream_file='Coho_KM_3.31.2020_2.csv', name_file='coho_names_3.csv'){
  coho_dat <- read.csv(file=coho_file)
  
  str(coho_dat)
  coho_dat$Managment.Unit..FRAM. <- as.factor(coho_dat$Managment.Unit..FRAM.)
  coho_dat$SubPopulation <- as.factor(coho_dat$SubPopulation)
  coho_dat$SaSI.Population  <- as.factor(coho_dat$SaSI.Population)
  coho_dat$Smolt.Abundance.Population <- as.factor(coho_dat$Smolt.Abundance.Population )
  str(coho_dat)
  
  
  #Remove blank entries for the SASI and non-SASI subpopulations
  coho_dat$SaSI.Population[coho_dat$SaSI.Population==""] <- NA
  coho_dat$SubPopulation[coho_dat$SubPopulation==""] <- NA
  
  #Create new columns that 1) identifies units treated as individual populations within the model and 2) defines the number of spawning adults (different for SASI versus FRAM pops)
  Population <- NA
  Spawners <- NA
  for(i in 1:nrow(coho_dat)){
    if(is.na(coho_dat$SaSI.Population[i]) & is.na(coho_dat$SubPopulation[i])){
      Population[i] <- as.character(coho_dat$Managment.Unit..FRAM.[i])
      Spawners[i] <- coho_dat$Age.3.Escapement..FRAM.[i]
    }
    if(is.na(coho_dat$SaSI.Population[i])==FALSE){
      Population[i] <- as.character(coho_dat$SaSI.Population[i])
      Spawners[i] <- coho_dat$SASI.Natural.Origin.Abundance[i]
    }
    if(is.na(coho_dat$SubPopulation[i])==FALSE & is.na(coho_dat$SaSI.Population[i])){
      Population[i] <- as.character(coho_dat$SubPopulation[i])
      Spawners[i] <- coho_dat$SupPopulation.Escapement[i]
    }
    if(is.na(coho_dat$SubPopulation[i])==FALSE & is.na(coho_dat$SaSI.Population[i])==FALSE){
      Population[i] <- as.character(coho_dat$SaSI.Population[i])
      Spawners[i] <- coho_dat$SASI.Natural.Origin.Abundance[i]
    }
  }
  coho_dat$Population <- as.factor(Population)
  coho_dat$Spawners <- Spawners
  
  #assign composite SASI spawners to discovery creek
  coho_dat$Spawners[coho_dat$Population=="Discovery Bay"] <- coho_dat$SASI.CompositeOrigin.Abundance[coho_dat$Population=="Discovery Bay"]
  
  #Calculate the harvest for both SASI and FRAM units using the FRAM harvest rate
  coho_dat$Harvest <- coho_dat$Spawners/(1/coho_dat$Harvest....FRAM.-1)
  
  #Check the population structure of the FRAM data 
  levels(coho_dat$Population)
  levels(coho_dat$Managment.Unit..FRAM.)
  
  #Remove the Queets and Clearwater individual subpopulations (FRAM unit Queets River Fall natural composed entirely of Queets and Clearwater)
  coho_dat <- subset(coho_dat, coho_dat$Population != "Queets" & coho_dat$Population != "Clearwater")
  coho_dat <- droplevels(coho_dat)
  levels(coho_dat$Population)
  
  #Negative count data?
  coho_dat[which(coho_dat$Spawners < 0),]
  coho_dat[which(coho_dat$Harvest < 0),]
  coho_dat[which(coho_dat$Smolt.Abundance < 0),]
  
  #Assume the negative spawner abundance count should be positive?
  coho_dat$Spawners[which(coho_dat$Spawners < 0)] <- 432.7124
  coho_dat[which(coho_dat$Spawners < 0),]
  
  #Filter out the really small streams
  aggregate(coho_dat$Spawners, by=list('Population'=coho_dat$Population), mean, na.rm=TRUE)
  unique(coho_dat$Population[coho_dat$Spawners < 1])
  #Bell and Johnson creek have chronically low spawner abundance and counts less than 0
  
  aggregate(coho_dat$Harvest, by=list('Population'=coho_dat$Population), mean, na.rm=TRUE)
  unique(coho_dat$Population[coho_dat$Harvest < 1])
  
  aggregate(coho_dat$Spawners, by=list('Population'=coho_dat$Population), mean, na.rm=TRUE)
  unique(coho_dat$Population[coho_dat$Smolt.Abundance < 1])
  
  #Bell, Johnson and Jimmy come latelt creek are pretty small, remove for now
  coho_dat <- subset(coho_dat, coho_dat$Population != 'Bell Creek' & coho_dat$Population !='Johnson Creek' & coho_dat$Population !='Jimmy Come Lately Creek')
  
  coho_dat <- subset(coho_dat, coho_dat$Population != 'Deep Creek' & coho_dat$Population !='McDonald Creek' & coho_dat$Population !='Siebert Creek' &
  coho_dat$Population != 'Salt Creek' & coho_dat$Population !='Discovery Bay' & coho_dat$Population !='East Twin Creek' &
  coho_dat$Population != 'West Twin Creek' & coho_dat$Population !='Northeast Hood Canal')
  
  coho_dat <- droplevels(coho_dat)
  
  #Create more succint coho data frame
  coho_dat_2 <- coho_dat[,c(1,6,14,15,22,23,24)]
  
  #Read in the stream distances and append to the main dataframe
  stream_dist <- read.csv(file=stream_file)
  stream_dist$SaSI.Population <- as.factor(stream_dist$SaSI.Population)
  stream_dist$Smolt.Population <- as.factor(stream_dist$Smolt.Population)
  stream_dist$Population <- as.factor(stream_dist$Population)
  
  stream_dist <- stream_dist[is.na(stream_dist$Population)==FALSE & stream_dist$Population !="",]
  stream_dist <- droplevels(stream_dist)
  
  levels(stream_dist$Population)[-which(levels(stream_dist$Population) %in% levels(coho_dat_2$Population))]
  levels(coho_dat_2$Population)[-which(levels(coho_dat_2$Population) %in% levels(stream_dist$Population))]
  
  coho_dat_2 <- merge(coho_dat_2, stream_dist[,c(4,8)], by='Population', all=FALSE)
  coho_dat_2 <- droplevels(coho_dat_2)
  
  #Read in CWT data file 
  CWT_dat <- read.csv(file=cwt_file)
  CWT_dat$Smolt.Ocean.Survival.Population <- as.factor(CWT_dat$Smolt.Ocean.Survival.Population)
  CWT_dat$Managment.Unit..FRAM. <- as.factor(CWT_dat$Managment.Unit..FRAM.)
  
  levels(CWT_dat$Managment.Unit..FRAM.)
  levels(CWT_dat$Smolt.Ocean.Survival.Population)
  
  #Add column to CWT data identifying Hatchery populations
  #Identify whether a population's marine survival estimates are based off of a hatchery
  x <- levels(CWT_dat$Smolt.Ocean.Survival.Population)[grep(' H',levels(CWT_dat$Smolt.Ocean.Survival.Population))]
  CWT_dat$Hatchery <- NA
  for(i in 1:length(levels(CWT_dat$Smolt.Ocean.Survival.Population))){
    if(unique(CWT_dat$Smolt.Ocean.Survival.Population)[i] %in% x){
      CWT_dat$Hatchery[CWT_dat$Smolt.Ocean.Survival.Population==unique(CWT_dat$Smolt.Ocean.Survival.Population)[i]] <- 1
    }else{
      CWT_dat$Hatchery[CWT_dat$Smolt.Ocean.Survival.Population==unique(CWT_dat$Smolt.Ocean.Survival.Population)[i]] <- 0
    }
  }
  
  #Identify which FRAM units have multiple marine survival time series in the CWT database
  x <- c()
  for(i in 1:length(unique(CWT_dat$Smolt.Ocean.Survival.Population))){
    x[i] <- unique(CWT_dat$Managment.Unit..FRAM.[CWT_dat$Smolt.Ocean.Survival.Population==unique(CWT_dat$Smolt.Ocean.Survival.Population)[i]])
  }
  
  x[duplicated(x)]
  which(x%in%x[duplicated(x)])    
  #CWT populations that are associated with the same FRAM unit
  unique(CWT_dat$Smolt.Ocean.Survival.Population)[which(x%in%x[duplicated(x)])]
  #Fram units with multiple CWT time-series
  CWT_dat$Managment.Unit..FRAM.[CWT_dat$Smolt.Ocean.Survival.Population %in% unique(CWT_dat$Smolt.Ocean.Survival.Population)[which(x%in%x[duplicated(x)])]]
  
  #For now, pick a single marine survival CWT time series for each FRAM unit for which there are duplicates
  CWT_dat <- subset(CWT_dat, CWT_dat$Smolt.Ocean.Survival.Population != 'Minter Crk H' & CWT_dat$Smolt.Ocean.Survival.Population != 'Baker H' & CWT_dat$Smolt.Ocean.Survival.Population != 'Satsop H')
  CWT_dat <- droplevels(CWT_dat)
  
  #Check all duplicates have been removed
  x <- c()
  for(i in 1:length(unique(CWT_dat$Smolt.Ocean.Survival.Population))){
    x[i] <- unique(CWT_dat$Managment.Unit..FRAM.[CWT_dat$Smolt.Ocean.Survival.Population==unique(CWT_dat$Smolt.Ocean.Survival.Population)[i]])
  }
  
  x[duplicated(x)]
  which(x%in%x[duplicated(x)])
  
  #Merge the marine survival and the core data
  coho_dat_full <- merge(coho_dat_2, CWT_dat[,c(1,2,5,6,7,9)], by.x=c('Population','Calendar.Year'), by.y=c('Managment.Unit..FRAM.', 'Calendar.Year'), all=TRUE)
  coho_dat_full <- droplevels(coho_dat_full)
  #Verify the structure data 
  str(coho_dat_full)
  
  #Correct spurrious coordinates for Green River marine entry point
  coho_dat_full$Longitude[coho_dat_full$Population=='Green River Wild'] <- -122.2145
  coho_dat_full$Latitude[coho_dat_full$Population=='Green River Wild'] <- 47.3519
  
  #Correct spurrious coordinates for Area 10E
  coho_dat_full$Longitude[coho_dat_full$Population=='Area 10E Miscellaneous Wild'] <- -122.8242
  coho_dat_full$Latitude[coho_dat_full$Population=='Area 10E Miscellaneous Wild'] <- 47.5896
  
  #Subset years to remove early years with little data
  coho_dat_full <- subset(coho_dat_full, coho_dat_full$Calendar.Year > 1985)
  coho_dat_full <- droplevels(coho_dat_full)
  
  #Filter out the Puyallup smolt data
  coho_dat_full$Smolt.Abundance[coho_dat_full$Population=='Puyallup River Wild'] <- NA
  
  #Filter out first couple years of the Quillayute time series
  coho_dat_full$Harvest[coho_dat_full$Population== 'Quillayute River Summer Natural' & coho_dat_full$Calendar.Year== 1986] <- NA
  coho_dat_full$Harvest[coho_dat_full$Population== 'Quillayute River Fall Natural' & coho_dat_full$Calendar.Year== 1986] <- NA
  coho_dat_full$Harvest[coho_dat_full$Population== 'Quillayute River Summer Natural' & coho_dat_full$Calendar.Year== 1987] <- NA
  coho_dat_full$Harvest[coho_dat_full$Population== 'Quillayute River Fall Natural' & coho_dat_full$Calendar.Year== 1987] <- NA
  
  #Include only puget sound pops
  #sound <- unique(na.omit(coho_dat_full$Population[coho_dat_full$Longitude > -123.5]))
  #sound <- unique(na.omit(coho_dat_full$Population[coho_dat_full$Longitude < -123.5]))
  
  #coho_dat_full <- subset(coho_dat_full, coho_dat_full$Population %in% sound)
  coho_dat_full <- droplevels(coho_dat_full)
  
  #add numeric code to populations
  coho_dat_full$pop <- as.numeric(coho_dat_full$Population)
  coho_dat_full$yr <- as.numeric(coho_dat_full$Calendar.Year) - min(as.numeric(coho_dat_full$Calendar.Year)) + 1
  
  #Population 10 is missing an entry for 2013 <- add in a dummy row
  new_row <- coho_dat_full[coho_dat_full$Calendar.Year==2012 & coho_dat_full$Population=='Area 7-7A Independent Wild',]
  new_row$Calendar.Year <- 2013
  new_row$yr <- 28
  new_row$Smolt.Abundance <- NA
  new_row$Spawners <- NA
  new_row$Harvest <- NA
  new_row$Fishery_Plus_Escapement <- NA
  new_row$Release_No <- NA
  
  which(coho_dat_full$Calendar.Year==2012 & coho_dat_full$Population=='Area 7-7A Independent Wild')
  coho_dat_full <- rbind(coho_dat_full[1:which(coho_dat_full$Calendar.Year==2012 & coho_dat_full$Population=='Area 7-7A Independent Wild'),], new_row, coho_dat_full[which(coho_dat_full$Calendar.Year==2014 & coho_dat_full$Population=='Area 7-7A Independent Wild'):nrow(coho_dat_full),])
  
  #Population 29 is missing an entry for 2004
  new_row <- coho_dat_full[coho_dat_full$Calendar.Year==2003 & coho_dat_full$Population=='Port Gamble Bay Wild',]
  new_row$Calendar.Year <- 2004
  new_row$yr <- 19
  new_row$Smolt.Abundance <- NA
  new_row$Spawners <- NA
  new_row$Harvest <- NA
  new_row$Fishery_Plus_Escapement <- NA
  new_row$Release_No <- NA
  
  which(coho_dat_full$Calendar.Year==2003 & coho_dat_full$Population=='Port Gamble Bay Wild')
  coho_dat_full <- rbind(coho_dat_full[1:which(coho_dat_full$Calendar.Year==2003 & coho_dat_full$Population=='Port Gamble Bay Wild'),], new_row, coho_dat_full[which(coho_dat_full$Calendar.Year==2005 & coho_dat_full$Population=='Port Gamble Bay Wild'):nrow(coho_dat_full),])
  
  #Population 20 is missing an entry for 2000 and 2001
  new_row_1 <- coho_dat_full[coho_dat_full$Calendar.Year==1999 & coho_dat_full$Population=='Grays Harbor Miscellaneous Wild',]
  new_row_2 <- coho_dat_full[coho_dat_full$Calendar.Year==1999 & coho_dat_full$Population=='Grays Harbor Miscellaneous Wild',]
  
  new_row_1$Calendar.Year <-2000
  new_row_1$yr <- 15
  new_row_1$Smolt.Abundance <- NA
  new_row_1$Spawners <- NA
  new_row_1$Harvest <- NA
  new_row_1$Fishery_Plus_Escapement <- NA
  new_row_1$Release_No <- NA
  
  new_row_2$Calendar.Year <-2001
  new_row_2$yr <- 16
  new_row_2$Smolt.Abundance <- NA
  new_row_2$Spawners <- NA
  new_row_2$Harvest <- NA
  new_row_2$Fishery_Plus_Escapement <- NA
  new_row_2$Release_No <- NA
  
  which(coho_dat_full$Calendar.Year==1999 & coho_dat_full$Population=='Grays Harbor Miscellaneous Wild')
  coho_dat_full <- rbind(coho_dat_full[1:which(coho_dat_full$Calendar.Year==1999 & coho_dat_full$Population=='Grays Harbor Miscellaneous Wild'),], new_row_1,new_row_2, coho_dat_full[which(coho_dat_full$Calendar.Year==2002 & coho_dat_full$Population=='Grays Harbor Miscellaneous Wild'):nrow(coho_dat_full),])
  
  #re-arrange data frame by population
  coho_dat_full <- coho_dat_full[order(coho_dat_full$pop),]
  
  #Convert latitude and longitude for each population to eastings and northings
  coho_space_Lat <- aggregate(c(coho_dat_full$Latitude), by=list('Population'= coho_dat_full$Population), FUN=median, na.rm=TRUE)
  coho_space_Long <- aggregate(c(coho_dat_full$Longitude), by=list('Population'= coho_dat_full$Population), FUN=median, na.rm=TRUE)
  
  coho_space <- merge(coho_space_Long, coho_space_Lat, by='Population')
  write.csv(coho_space, file='coho_coords.csv')
  colnames(coho_space) <- c('Population', 'long', 'lat')
  coordinates(coho_space) <-c("long", 'lat')
  proj4string(coho_space) <- CRS("+proj=longlat +datum=WGS84")
  utm_dat <- spTransform(coho_space, CRS("+proj=utm +zone=10T ellps=WGS84"))
  coord <-utm_dat@coords
  
  #Assign populations to either puget sound
  coho_dat_full$Basin <- NA
  coho_dat_full$Basin[coho_dat_full$Longitude < -123.80] <- 0
  coho_dat_full$Basin[coho_dat_full$Longitude >= -123.80] <- 1
  
  #Add column with more plot-friendly names to the data file
  coho_names <- read.csv(file=name_file)
  coho_names$Name <- as.factor(coho_names$Name)
  coho_dat_full <- merge(coho_dat_full, coho_names, by='pop', all=TRUE, no.dups=TRUE)
  coho_dat_full <- droplevels(coho_dat_full)
  
  
  #Create single list object to store coho data frame and coordinates
  dat_list <- list()
  dat_list$coho_dat_full <- coho_dat_full
  dat_list$coord <- coord
  
  return(dat_list)
}
dat_list <- dat_func(coho_file='Coho data_3-30-20.csv', cwt_file='CWT_FRAM_Matches_complete20200303.csv', stream_file='Coho_KM_3.31.2020_2.csv', name_file='coho_names_3_2.csv')
#dat_list <- dat_func(coho_file='Coho data_3-30-20.csv', cwt_file='CWT_FRAM_Matches_complete20200303.csv', stream_file='Coho_KM_3.31.2020_2.csv', name_file='coho_names_3.csv')

#Read in official state forecast data, change names to match other data and identify overlapping populations#####
state_cast <- read.csv(file='ForecastCompilation_SourceDataPre_MB080520-1.csv')
state_cast <- subset(state_cast, state_cast$StockLongName %in% c(
                       "Area 10 Miscellaneous Wild UnMarked",
                       "Area 10E Miscellaneous Wild UnMarked",
                       "Area 11 Miscellaneous Wild UnMarked",
                       "Area 12/12B Wild UnMarked",
                       "Area 12A Wild UnMarked",
                       "Area 12C/12D Wild UnMarked",
                       "Area 13 Miscellaneous Wild UnMarked",
                       "Area 13A Miscellaneous Wild UnMarked",
                       "Area 13B Miscellaneous Wild UnMarked",
                       "Area 7/7A Independent Wild UnMarked",
                       "Baker (Skagit) Wild UnMarked",
                       "Chehalis River Wild UnMarked",
                       "Deschutes River (WA) Wild UnMarked",
                       "Dungeness River Wild UnMarked",
                       "East JDF Miscellaneous Wild UnMarked",
                       "Elwha River Wild UnMarked",
                       "Grays Harbor Miscellaneous Wild UnMarked",
                       "Green River Wild UnMarked",
                       "Hoh River Wild UnMarked",
                       "Humptulips River Wild UnMarked",
                       "Lake Washington Wild UnMarked",
                       "Nisqually River Wild UnMarked",
                       "Nooksack River Wild UnMarked",
                       "Port Gamble Bay Wild UnMarked",
                       "Puyallup River Wild UnMarked",
                       "Queets River Fall Natural UnMarked",
                       "Quillayute River Fall Natural UnMarked",
                       "Quillayute River Summer Natural UnMarked",
                       "Quinault River Fall Natural UnMarked",
                       "Samish River Wild UnMarked",
                       "Skagit River Wild UnMarked",
                       "Skokomish River Wild UnMarked",
                       "Snohomish River Hatchery UnMarked",
                       "Stillaguamish River Wild UnMarked",
                       "West JDF Miscellaneous Wild UnMarked",
                       "Willapa Bay Natural UnMarked"
                       #"Wash Early Wild UnMarked",
                       #"Wash Late Wild UnMarked"
                       ))
state_cast <- droplevels(state_cast)
state_cast <- state_cast[-which(duplicated(state_cast[,c(1,3)])),]
state_cast$RunSize <- state_cast$RunSize/1.2317
state_cast$StockLongName <- as.factor(state_cast$StockLongName)
levels(state_cast$StockLongName)[levels(state_cast$StockLongName)=="Area 10 Miscellaneous Wild UnMarked"] <- "Area 10 Miscellaneous Wild"
levels(state_cast$StockLongName)[levels(state_cast$StockLongName)=="Area 10E Miscellaneous Wild UnMarked"] <- "Area 10E Miscellaneous Wild"
levels(state_cast$StockLongName)[levels(state_cast$StockLongName)=="Area 11 Miscellaneous Wild UnMarked"] <- "Area 11 Miscellaneous Wild"
levels(state_cast$StockLongName)[levels(state_cast$StockLongName)=="Area 12/12B Wild UnMarked"] <- "Area 12-12B Wild"
levels(state_cast$StockLongName)[levels(state_cast$StockLongName)=="Area 12A Wild UnMarked"] <- "Area 12A Wild"
levels(state_cast$StockLongName)[levels(state_cast$StockLongName)=="Area 12C/12D Wild UnMarked"] <- "Area 12C-12D Wild"
levels(state_cast$StockLongName)[levels(state_cast$StockLongName)=="Area 13 Miscellaneous Wild UnMarked"] <- "Area 13 Miscellaneous Wild"
levels(state_cast$StockLongName)[levels(state_cast$StockLongName)=="Area 13A Miscellaneous Wild UnMarked"] <- "Area 13A Miscellaneous Wild"
levels(state_cast$StockLongName)[levels(state_cast$StockLongName)=="Area 13B Miscellaneous Wild UnMarked"] <- "Area 13B Miscellaneous Wild"
levels(state_cast$StockLongName)[levels(state_cast$StockLongName)=="Area 7/7A Independent Wild UnMarked"] <- "Area 7-7A Independent Wild"
levels(state_cast$StockLongName)[levels(state_cast$StockLongName)=="Baker (Skagit) Wild UnMarked"] <- "Baker (Skagit) Wild"
levels(state_cast$StockLongName)[levels(state_cast$StockLongName)=="Chehalis River Wild UnMarked"] <- "Chehalis River Wild"
levels(state_cast$StockLongName)[levels(state_cast$StockLongName)=="Deschutes River (WA) Wild UnMarked"] <- "Deschutes River (WA) Wild"
levels(state_cast$StockLongName)[levels(state_cast$StockLongName)=="Dungeness River Wild UnMarked"] <- "Dungeness River Wild"
levels(state_cast$StockLongName)[levels(state_cast$StockLongName)=="East JDF Miscellaneous Wild UnMarked"] <- "East JDF Miscellaneous Wild"
levels(state_cast$StockLongName)[levels(state_cast$StockLongName)=="Elwha River Wild UnMarked"] <- "Elwha River Wild"
levels(state_cast$StockLongName)[levels(state_cast$StockLongName)=="Grays Harbor Miscellaneous Wild UnMarked"] <- "Grays Harbor Miscellaneous Wild"
levels(state_cast$StockLongName)[levels(state_cast$StockLongName)=="Green River Wild UnMarked"] <- "Green River Wild"
levels(state_cast$StockLongName)[levels(state_cast$StockLongName)=="Hoh River Wild UnMarked"] <- "Hoh River Wild" 
levels(state_cast$StockLongName)[levels(state_cast$StockLongName)=="Humptulips River Wild UnMarked"] <- "Humptulips River Wild"
levels(state_cast$StockLongName)[levels(state_cast$StockLongName)=="Lake Washington Wild UnMarked"] <- "Lake Washington Wild"
levels(state_cast$StockLongName)[levels(state_cast$StockLongName)=="Nisqually River Wild UnMarked"] <- "Nisqually River Wild"
levels(state_cast$StockLongName)[levels(state_cast$StockLongName)=="Nooksack River Wild UnMarked"] <- "Nooksack River Wild"
levels(state_cast$StockLongName)[levels(state_cast$StockLongName)=="Port Gamble Bay Wild UnMarked"] <- "Port Gamble Bay Wild"
levels(state_cast$StockLongName)[levels(state_cast$StockLongName)=="Puyallup River Wild UnMarked"] <- "Puyallup River Wild"
levels(state_cast$StockLongName)[levels(state_cast$StockLongName)=="Queets River Fall Natural UnMarked"] <- "Queets River Fall Natural"
levels(state_cast$StockLongName)[levels(state_cast$StockLongName)=="Quillayute River Fall Natural UnMarked"] <- "Quillayute River Fall Natural"
levels(state_cast$StockLongName)[levels(state_cast$StockLongName)=="Quillayute River Summer Natural UnMarked"] <- "Quillayute River Summer Natural"
levels(state_cast$StockLongName)[levels(state_cast$StockLongName)=="Quinault River Fall Natural UnMarked"] <- "Quinault River Fall Natural"
levels(state_cast$StockLongName)[levels(state_cast$StockLongName)=="Samish River Wild UnMarked"] <- "Samish River Wild" 
levels(state_cast$StockLongName)[levels(state_cast$StockLongName)=="Skagit River Wild UnMarked"] <- "Skagit River Wild"
levels(state_cast$StockLongName)[levels(state_cast$StockLongName)=="Skokomish River Wild UnMarked"] <- "Skokomish River Wild"
levels(state_cast$StockLongName)[levels(state_cast$StockLongName)=="Snohomish River Hatchery UnMarked"] <- "Snohomish River Wild" 
levels(state_cast$StockLongName)[levels(state_cast$StockLongName)=="Stillaguamish River Wild UnMarked"] <- "Stillaguamish River Wild"
levels(state_cast$StockLongName)[levels(state_cast$StockLongName)=="West JDF Miscellaneous Wild UnMarked"] <- "West JDF Miscellaneous Wild"
levels(state_cast$StockLongName)[levels(state_cast$StockLongName)=="Willapa Bay Natural UnMarked"] <- "Willapa Bay Natural" 
#levels(state_cast$StockLongName)[levels(state_cast$StockLongName)== "Wash Early Wild UnMarked"] <- 'Wash Early Wild'
#levels(state_cast$StockLongName)[levels(state_cast$StockLongName)== "Wash Late Wild UnMarked"] <- 'Wash Late Wild'

cast_pop <- c(which(levels(dat_list$coho_dat_full$Population) %in% levels(state_cast$StockLongName)))

state_cast <- state_cast[,c(1,3,5)]
state_cast <- spread(data=state_cast, key=YEAR, value=RunSize)
#####

length(cast_pop)
#Execute the stan model and plot output
exec_fun <- function(coho_dat_full=dat_list$coho_dat_full, coord=dat_list$coord, exec=FALSE, n_iter=1000, n_thin=1, n_adapt=0.8, n_tree=10, n_chain=3, plots=TRUE,
                     sigma_esc = 0.2, mod_filename='LD_coho_forecast_6_2_4_fit', cast_pop=c(1:36)){
  
  #Create vector of basin identifiers for coastal and puget sound populations
  basin <- aggregate(coho_dat_full$Basin, by=list('Population'=coho_dat_full$pop), mean, na.rm=TRUE)
  coast <- basin$Population[basin$x==0]
  sound <- basin$Population[basin$x==1]  
  
  #Compute number of populations
  n_pop <- length(unique(coho_dat_full$Population))
  
  #Compute populations with smolt data
  pop_smolt <- unique(coho_dat_full$pop[which(is.na(coho_dat_full$Smolt.Abundance)==FALSE)])
  #Compute populations with marine survival data
  pop_MS <- unique(coho_dat_full$pop[which(is.na(coho_dat_full$Release_No)==FALSE)])
  #Compute populations with harvest data
  pop_catch <- unique(coho_dat_full$pop[which(is.na(coho_dat_full$Harvest)==FALSE)])
  #Compute populations with escapement data
  pop_esc <- unique(coho_dat_full$pop[which(is.na(coho_dat_full$Spawners)==FALSE)])
  
  #Compute number of populations with each data type available
  n_pop_smolt <- length(pop_smolt)
  n_pop_MS <- length(pop_MS)
  n_pop_catch <- length(pop_catch)
  n_pop_esc <- length(pop_esc)
  
  #create vectors to store the number of years of smolt and escapement data available for each population being considered
  n_year_true_smolt <- vector(length=n_pop_smolt)
  n_year_true_esc <- vector(length=n_pop_esc)
  n_year_true_harvest <- vector(length=n_pop_catch)
  n_year_true_MS <- vector(length=n_pop_MS)
  
  #Calculate the number of years for which each data type is available for each population
  #Smolt data
  for (i in 1:n_pop_smolt){
    n_year_true_smolt[i] <- length(na.omit(coho_dat_full$Smolt.Abundance[coho_dat_full$pop==pop_smolt[i]]))
  }
  #escapement data
  for (i in 1:n_pop_esc){
    n_year_true_esc[i] <- length(na.omit(coho_dat_full$Spawners[coho_dat_full$pop==pop_esc[i]]))
  }
  #harvest data
  for (i in 1:n_pop_catch){
    n_year_true_harvest[i] <- length(na.omit(coho_dat_full$Harvest[coho_dat_full$pop==pop_catch[i]]))
  }
  #marine survival data
  for (i in 1:n_pop_MS){
    n_year_true_MS[i] <- length(na.omit(coho_dat_full$Fishery_Plus_Escapement[coho_dat_full$pop==pop_MS[i]]))
  }
  
  #Compile all the years across populations for which there is escapement and smolt data (not NAs)
  esc_true <- coho_dat_full$yr[which(is.na(coho_dat_full$Spawners)==FALSE)]
  smolt_true <- coho_dat_full$yr[which(is.na(coho_dat_full$Smolt.Abundance)==FALSE)]
  harvest_true <- coho_dat_full$yr[which(is.na(coho_dat_full$Harvest )==FALSE)]
  MS_true <- coho_dat_full$yr[which(is.na(coho_dat_full$Fishery_Plus_Escapement)==FALSE)]
  
  #Compile all the actual observations of smolt and escapement data that are not NA
  smolt_dat <- coho_dat_full$Smolt.Abundance[which(is.na(coho_dat_full$Smolt.Abundance)==FALSE)]                 
  esc_dat <- coho_dat_full$Spawners[which(is.na(coho_dat_full$Spawners)==FALSE)]  
  harvest_dat <- coho_dat_full$Harvest[which(is.na(coho_dat_full$Harvest )==FALSE)]      
  MS_dat_x <- coho_dat_full$Fishery_Plus_Escapement[which(is.na(coho_dat_full$Fishery_Plus_Escapement)==FALSE)]      
  MS_dat_N <- coho_dat_full$Release_No [which(is.na(coho_dat_full$Release_No)==FALSE)]      
  
  #Create slice points for cutting up the intact (non-NA) smolt and escapement data by population (stan does not accept ragged data structures)
  slice_smolt_start <- vector(length=n_pop_smolt)
  slice_smolt_start[1] <- 1
  slice_smolt_end <- vector(length=n_pop_smolt)
  slice_smolt_end[n_pop_smolt] <- length(smolt_true)
  
  slice_esc_start <- vector(length=n_pop_esc)
  slice_esc_start[1] <- 1
  slice_esc_end <- vector(length=n_pop_esc)
  slice_esc_end[n_pop_esc] <- length(esc_true)
  
  slice_harvest_start <- vector(length=n_pop_catch)
  slice_harvest_start[1] <- 1
  slice_harvest_end <- vector(length=n_pop_catch)
  slice_harvest_end[n_pop_catch] <- length(harvest_true)
  
  slice_MS_start <- vector(length=n_pop_MS)
  slice_MS_start[1] <- 1
  slice_MS_end <- vector(length=n_pop_MS)
  slice_MS_end[n_pop_MS] <- length(MS_true)
  
  for(p in 2:n_pop_smolt){
    slice_smolt_start[p] <- slice_smolt_start[p-1] + length(na.omit(coho_dat_full$Smolt.Abundance[coho_dat_full$pop==pop_smolt[p-1]]))
  }
  slice_smolt_end[1:(n_pop_smolt-1)] <- slice_smolt_start[2:n_pop_smolt]-1
  
  
  for(p in 2:n_pop_esc){
    slice_esc_start[p] <- slice_esc_start[p-1] + length(na.omit(coho_dat_full$Spawners[coho_dat_full$pop==pop_esc[p-1]]))
  }
  slice_esc_end[1:(n_pop_esc-1)] <- slice_esc_start[2:n_pop_esc]-1
  
  
  for(p in 2:n_pop_catch){
    slice_harvest_start[p] <- slice_harvest_start[p-1] + length(na.omit(coho_dat_full$Harvest[coho_dat_full$pop==pop_catch[p-1]]))
  }
  slice_harvest_end[1:(n_pop_catch-1)] <- slice_harvest_start[2:n_pop_catch]-1
  
  
  for(p in 2:n_pop_MS){
    slice_MS_start[p] <- slice_MS_start[p-1] + length(na.omit(coho_dat_full$Fishery_Plus_Escapement[coho_dat_full$pop==pop_MS[p-1]]))
  }
  slice_MS_end[1:(n_pop_MS-1)] <- slice_MS_start[2:n_pop_MS]-1
  
  #test if the cut points are correctly pulling the right populations and plot
  pdf(file='data.pdf')
  par(mfrow=c(3,3))
  for(i in 1:n_pop_smolt){
    print(all(smolt_dat[slice_smolt_start[i]:slice_smolt_end[i]] == na.omit(coho_dat_full$Smolt.Abundance[coho_dat_full$pop==pop_smolt[i]])))
    plot(coho_dat_full$Smolt.Abundance[coho_dat_full$pop==pop_smolt[i]], ylab='Smolt Abundance', xlab='Year', main=paste(unique(coho_dat_full$Name)[pop_smolt[i]]), type='o', pch=16, col='black')
  }
  for(i in 1:n_pop_esc){
    print(all(esc_dat[slice_esc_start[i]:slice_esc_end[i]] == na.omit(coho_dat_full$Spawners[coho_dat_full$pop==pop_esc[i]])))
    plot(coho_dat_full$Spawners[coho_dat_full$pop==pop_esc[i]], ylab='Escapement Abundance', xlab='Year', main=paste(unique(coho_dat_full$Name)[pop_esc[i]]), type='o', pch=16, col='black')
  }
  for(i in 1:n_pop_catch){
    print(all(harvest_dat[slice_harvest_start[i]:slice_harvest_end[i]] == na.omit(coho_dat_full$Harvest[coho_dat_full$pop==pop_catch[i]])))
    plot(coho_dat_full$Harvest[coho_dat_full$pop==pop_catch[i]], ylab='Harvest Abundance', xlab='Year', main=paste(unique(coho_dat_full$Name)[pop_catch[i]]), type='o', pch=16, col='black')
  }
  for(i in 1:n_pop_MS){
    print(all(MS_dat_x[slice_MS_start[i]:slice_MS_end[i]] == na.omit(coho_dat_full$Fishery_Plus_Escapement[coho_dat_full$pop==pop_MS[i]])))
    plot(coho_dat_full$Fishery_Plus_Escapement[coho_dat_full$pop==pop_MS[i]]/coho_dat_full$Release_No[coho_dat_full$pop==pop_MS[i]], 
         ylab='Marine survival', xlab='Year', main=paste(unique(coho_dat_full$Name)[pop_MS[i]]), type='o', pch=16, col='black')
  }
  dev.off()
  
  
  #####
  #Fit stan model####
  if(exec==TRUE){
    mod_fit <- stan(file = 'LD_coho_forecast_6_2_4.stan', data = list(
      n_year <- length(unique(coho_dat_full$Calendar.Year)),
      n_pop <- length(unique(coho_dat_full$pop)),
      pop_smolt <- pop_smolt,
      pop_esc <- pop_esc,
      pop_catch <- pop_catch,
      pop_MS <- pop_MS,
      n_pop_smolt <-n_pop_smolt,
      n_pop_esc <- n_pop_esc,
      n_pop_catch <- n_pop_catch,
      n_pop_MS <- n_pop_MS,
      stream_dist = aggregate(coho_dat_full$KM, by=list('Population'=coho_dat_full$Population), mean)$x,
      smolt_true <- smolt_true,
      esc_true <- esc_true,
      MS_true <- MS_true,
      harvest_true <- harvest_true,
      smolt_dat <- smolt_dat,
      esc_dat <- esc_dat,
      harvest_dat <- harvest_dat,
      MS_dat_x <- round(MS_dat_x),
      MS_dat_N <- round(MS_dat_N),
      n_smolt <- length(smolt_dat),
      n_esc <- length(esc_dat),
      n_harvest <- length(harvest_dat),
      n_MS <- length(MS_dat_x),
      #sigma_catch <- sigma_catch,
      sigma_esc <- sigma_esc,
      #sigma_smolt <- sigma_smolt,
      n_hatchery <- length(unique(na.omit(coho_dat_full$pop[coho_dat_full$Hatchery==1]))),
      hatchery <- which(pop_MS %in% unique(na.omit(coho_dat_full$pop[coho_dat_full$Hatchery==1]))),
      wild <- which(pop_MS %in% unique(na.omit(coho_dat_full$pop[coho_dat_full$Hatchery==0]))),
      slice_smolt_start <- slice_smolt_start,
      slice_smolt_end <- slice_smolt_end,
      slice_esc_start <- slice_esc_start,
      slice_esc_end <-  slice_esc_end,
      slice_harvest_start <- slice_harvest_start,
      slice_harvest_end <-  slice_harvest_end,
      slice_MS_start <- slice_MS_start,
      slice_MS_end <-  slice_MS_end,
      u=matrix(1, nrow=1, ncol=n_pop),
      dist = as.matrix(dist(coord, 'euclidean', diag=TRUE, upper=TRUE)/10000)
      #dist = as.matrix(dist(coord[sound,]/10000, 'euclidean', diag=TRUE, upper=TRUE))
    ),iter = n_iter, chains = n_chain, control=list(adapt_delta=n_adapt, max_treedepth=n_tree), thin=n_thin, seed=666)
    
    saveRDS(mod_fit, paste(mod_filename,'.rds',sep=''))
  }
  
  if(exec==FALSE){
    mod_fit <- readRDS(paste(mod_filename,'.rds',sep=''))
  }
  
  df_1 <- as.data.frame(mod_fit)
  
  if(plots==TRUE){
    #Diagnostic Plots####
    #Year adjustment for plotting the xaxes on graphs
    yr <- min(coho_dat_full$Calendar.Year)-1
    n_year <- length(unique(coho_dat_full$Calendar.Year))
    
    #Plot the spatial kernal
    pdf(file='spatial_kernel_II_AR.pdf')
    rho_dist <- df_1$rho_dist_MS
    alpha_dist <- df_1$alpha_dist_MS
    sigma_dist <- df_1$sigma_dist_MS
    
    dist <- seq(0,30,0.1)
    cor_vec_med <- vector(length=length(dist))
    cor_vec_lower_50 <- vector(length=length(dist))
    cor_vec_upper_50 <- vector(length=length(dist))
    cor_vec_lower_95 <- vector(length=length(dist))
    cor_vec_upper_95 <- vector(length=length(dist))
    
    #Puget sound populations
    for(i in 1:length(dist)){
      cor_vec_med[i] <- exp_func(d=dist[i], rho_dist=median(rho_dist), alpha_dist=median(alpha_dist))
      cor_vec_lower_50[i] <- exp_func(d=dist[i], rho_dist=cred.50(rho_dist)[1], alpha_dist=cred.50(alpha_dist)[1])
      cor_vec_upper_50[i] <- exp_func(d=dist[i], rho_dist=cred.50(rho_dist)[2], alpha_dist=cred.50(alpha_dist)[2])
      cor_vec_lower_95[i] <- exp_func(d=dist[i], rho_dist=cred(rho_dist)[1], alpha_dist=cred(alpha_dist)[1])
      cor_vec_upper_95[i] <- exp_func(d=dist[i], rho_dist=cred(rho_dist)[2], alpha_dist=cred(alpha_dist)[2])
    }
    par(mfrow=c(2,1), mar=c(0,0,0,0), oma=c(10,12,10,12))
    mat <- matrix(c(1,1,1,1,0,0,0,0,0,0,
                    2,2,2,2,0,3,0,4,0,5,
                    6,6,6,6,0,7,0,8,0,9),3,10, byrow=TRUE)
    layout(mat=mat, widths=rep(c(rep(1, 4), 0.15, 0.4, 0.5, 0.4, 0.5, 0.4),3), heights=rep(1,30))
    hist(as.vector(dist(coord, 'euclidean', diag=FALSE, upper=FALSE)/10000), xlim=c(0,30), xaxs='i', yaxs='i', main=NA, xaxt='n', col=rgb(13, 82, 144, max = 255, alpha = 195), breaks=30, border=rgb(13, 82, 144, max = 255, alpha = 195), yaxt='n')
    axis(side=2, at=c(0, 10,20,30,40), labels=c('',10,20,30,40))
    mtext(side=2, 'Frequency', cex=1, line=2.5)
    text(x=1, y=38, 'A', cex=1.1)
    plot(1:length(dist), cor_vec_med, col=rgb(13, 82, 144, max = 255, alpha = 255), type='l', ylim=c(0, 1), xaxt='n', xaxs='i', lwd=1.25)
    xx <- c(1:(length(dist)),(length(dist)):1)
    yy <- c(cor_vec_lower_50, rev(cor_vec_upper_50))
    polygon(x=xx, y=yy, density = -1, border = rgb(13, 82, 144, max = 255, alpha = 255), lwd=0.75, col=rgb(13, 82, 144, max = 255, alpha = 150))
    yy <- c(cor_vec_lower_95, rev(cor_vec_upper_95))
    polygon(x=xx, y=yy, density = -1, border = rgb(13, 82, 144, max = 255, alpha = 255), lwd=0.75, col=rgb(13, 82, 144, max = 255, alpha = 100))
    mtext(side=1, 'Euclidean distance (km)', cex=1, line=2.5)
    
    mtext(side=2, 'Marine survival',cex=1, line=3.5)
    mtext(side=2, 'correlation',cex=1, line=2.25)
    axis(side=1, at=seq(0, length(dist)-1, 10), labels=seq(0, length(dist)-1, 10))
    text(x=15, y=0.045, 'B', cex=1.1)
    
    plot(1, median(alpha_dist),yaxt='n', xaxt='n', ylim=c(0,quantile(df_1$alpha_dist, 0.98)), xlim=c(0.5,1.5), pch=16, col=rgb(13, 82, 144, max = 255, alpha = 255), cex=1.5)
    arrows((1), quantile(df_1$alpha_dist, 0.25), 1, quantile(df_1$alpha_dist, 0.75), angle=90, code=3, length=0.0, lwd=5, lty=1, col=rgb(13, 82, 144, max = 255, alpha = 175))
    arrows((1), quantile(df_1$alpha_dist, 0.025), 1, quantile(df_1$alpha_dist, 0.975), angle=90, code=3, length=0.01, lwd=1.5, lty=1, col=rgb(13, 82, 144, max = 255, alpha = 150))
    mtext(side=1, expression(gamma), line=1.0, at=1, cex=1.2)
    text(x=1, y=0.025, 'C', cex=1.1)
    axis(side=4)
    
    plot(1, median(sigma_dist),yaxt='n', xaxt='n', pch=16, col=rgb(13, 82, 144, max = 255, alpha = 255), cex=1.5, xlim=c(0.5,1.5), ylim=c(0,quantile(df_1$alpha_dist, 0.98)))
    arrows((1), quantile(df_1$sigma_dist, 0.25), 1, quantile(df_1$sigma_dist, 0.75), angle=90, code=3, length=0.0, lwd=5, lty=1, col=rgb(13, 82, 144, max = 255, alpha = 175))
    arrows((1), quantile(df_1$sigma_dist, 0.025), 1, quantile(df_1$sigma_dist, 0.975), angle=90, code=3, length=0.01, lwd=1.5, lty=1, col=rgb(13, 82, 144, max = 255, alpha = 150))
    mtext(side=1, expression(sigma[d]), line=1.0, at=1, cex=1.2)
    text(x=1, y=0.025, 'D', cex=1.1)
    axis(side=4)
    
    plot(1, median(rho_dist),yaxt='n', xaxt='n', pch=16, col=rgb(13, 82, 144, max = 255, alpha = 255), cex=1.5, ylim=c(0,quantile(df_1$rho_dist, 0.98)), xlim=c(0.5,1.5))
    arrows((1), quantile(df_1$rho_dist, 0.25), 1, quantile(df_1$rho_dist, 0.75), angle=90, code=3, length=0.0, lwd=5, lty=1, col=rgb(13, 82, 144, max = 255, alpha = 175))
    arrows((1), quantile(df_1$rho_dist, 0.025), 1, quantile(df_1$rho_dist, 0.975), angle=90, code=3, length=0.01, lwd=1.5, lty=1, col=rgb(13, 82, 144, max = 255, alpha = 150))
    axis(side=4, at=seq(0,quantile(df_1$rho_dist, 0.98),0.5), labels=seq(0,quantile(df_1$rho_dist, 0.98),0.5)*10)
    mtext(side=1, expression(rho), line=1.0, at=1, cex=1.2)
    mtext(side=4, 'Length Scale (KM)', line=2.5, cex=1)
    text(x=1, y=0.35, 'E', cex=1.1)
    dev.off()
    #plot correlation matrix
    corr_mat <- matrix(unlist(colMedian(df_1[,grep('corr_mat', colnames(df_1))])), nrow=n_pop, ncol=n_pop, byrow=FALSE)[cast_pop,cast_pop]
    colnames(corr_mat) <- unique(dat_list$coho_dat_full$Name)[cast_pop]
    rownames(corr_mat) <- unique(dat_list$coho_dat_full$Name)[cast_pop]
    pdf(file='corr_matrix_plot.pdf')
    corrplot(corr_mat,title =NA, method = "square", outline = T, addgrid.col = "darkgray", order="hclust", mar = c(4,0,4,0), addrect = 1, rect.col = "black", rect.lwd = 1, cl.pos = "b", tl.col = "black", tl.cex = 0.7, cl.cex = 0.7, cl.lim=c(0,1),
             col = colorRampPalette(c("white","white",rgb(13, 82, 144, max = 255, alpha = 255)))(100))
    mtext(side=1, 'Marine survival correlation', at=cast_pop[16], line=0.5, cex=0.85)
    dev.off()
    #Plot posteriors of key model parameters
    pdf('SR_posteriors_AR.pdf')
    n_pop_cast <- length(cast_pop)
    par(mfrow=c(2,1), mar=c(1,0,0,0), oma=c(16,8,6,6))
    mat <- matrix(c(rep(1,n_pop_cast),0,2,0,3,
                    rep(4, n_pop_cast),0,5,0,6,
                    rep(7, n_pop_cast),0,8,0,9),3,n_pop_cast+4, byrow=TRUE)
    layout(mat=mat, widths=c(rep(0.5,n_pop_cast), 0.0, 0.5, 0.25, 0.5, rep(0.5,n_pop_cast), 0.0, 0.5,0.25,0.5), heights=c(1))
    alpha <- df_1[,grep('log_alpha', colnames(df_1))][,cast_pop]
    R_max <- df_1[,grep('log_R_max', colnames(df_1))][,cast_pop]
    sigma_R <- exp(df_1[,grep('log_sigma_R', colnames(df_1))])[,cast_pop]
    plot(1:n_pop_cast, colMedian(alpha), col='skyblue4', ylim=c(min(cred.2(alpha)[1,]),max(cred.2(alpha)[2,])), xlim=c(1,n_pop_cast), xaxt='n', yaxt='n', cex=1)
    axis(side=2, at=c(log(25), log(50), log(100), log(200), log(400), log(800), log(1600), log(3200), log(6400), log(6400*2), log(6400*4), log(6400*8)), labels = c(25, 50, 100, 200, 400,800,1600, 3200, 6400, 6400*2, 6400*4, 6400*8), las=2)
    mtext(side=2, expression(alpha), line=4.2, cex=1)
    #text(x=n_pop_cast, y=log(4650), 'H', srt=90)
    
    mtext(side=2, 'smolts/spawner', line=3.2, cex=0.85)
    arrows(1:n_pop_cast, cred.2(alpha)[1,], 1:n_pop_cast, cred.2(alpha)[2,], angle=90, code=3, length=0.01, lwd=1.5, lty=1, col='skyblue4')
    arrows(1:n_pop_cast, cred.3(alpha)[1,], 1:n_pop_cast, cred.3(alpha)[2,], angle=90, code=3, length=0.0, lwd=3, lty=1, col='skyblue4')
    plot((n_pop_cast+1), median(df_1$mu_alpha), col='skyblue4', ylim=c(min(cred.2(alpha)[1,]),max(cred.2(alpha)[2,])), xaxt='n', yaxt='n',cex=1)
    #text(x=n_pop_cast, y=log(4650), 'I', srt=90)
    
    arrows((n_pop_cast+1), quantile(df_1$mu_alpha, 0.25), (n_pop_cast+1), quantile(df_1$mu_alpha, 0.75), angle=90, code=3, length=0.0, lwd=3, lty=1, col='skyblue4')
    arrows((n_pop_cast+1), quantile(df_1$mu_alpha, 0.025), (n_pop_cast+1), quantile(df_1$mu_alpha, 0.975), angle=90, code=3, length=0.01, lwd=1.5, lty=1, col='skyblue4')
    plot(n_pop_cast+1, median(df_1$sigma_alpha), col='skyblue4', ylim=c(quantile(df_1$sigma_alpha, 0.0), quantile(df_1$sigma_alpha, 0.995)), xaxt='n', yaxt='n',cex=1)
    #text(x=n_pop_cast, y=1.27, 'J', srt=90)
    
    arrows((n_pop_cast+1), quantile(df_1$sigma_alpha, 0.25), (n_pop_cast+1), quantile(df_1$sigma_alpha, 0.75), angle=90, code=3, length=0.0, lwd=3, lty=1, col='skyblue4')
    arrows((n_pop_cast+1), quantile(df_1$sigma_alpha, 0.025), (n_pop_cast+1), quantile(df_1$sigma_alpha, 0.975), angle=90, code=3, length=0.01, lwd=1.5, lty=1, col='skyblue4')
    axis(side=4)
    plot(1:n_pop_cast, colMedian(R_max), col='goldenrod', ylim=c(min(cred.2(R_max)[1,]),max(cred.2(R_max)[2,])), xlim=c(1,n_pop_cast), xaxt='n', yaxt='n', cex=1)
    axis(side=2, at=c(log(7), log(15),log(30),log(60), log(125), log(250), log(500), log(1000), log(2000), log(4000), log(8000)), labels = c(7, 15, 30, 60, 125, 250, 500,1000,2000,4000, 8000), las=2)
    mtext(side=2, expression('R'[max]), line=4.2, cex=1)
    mtext(side=2, 'smolts/KM', line=3.2, cex=0.85)
    #text(x=n_pop_cast, y=log(3800), 'D', srt=90)
    
    arrows(1:n_pop_cast, cred.2(R_max)[1,], 1:n_pop_cast, cred.2(R_max)[2,], angle=90, code=3, length=0.01, lwd=1.5, lty=1, col='goldenrod')
    arrows(1:n_pop_cast, cred.3(R_max)[1,], 1:n_pop_cast, cred.3(R_max)[2,], angle=90, code=3, length=0.0, lwd=3, lty=1, col='goldenrod')
    plot((n_pop_cast+1), median(df_1$mu_R_max), col='goldenrod', ylim=c(min(cred.2(R_max)[1,]),max(cred.2(R_max)[2,])), xaxt='n', yaxt='n',cex=1)
    #text(x=n_pop_cast, y=log(3800), 'E', srt=90)
    arrows((n_pop_cast+1), quantile(df_1$mu_R_max, 0.25), (n_pop_cast+1), quantile(df_1$mu_R_max, 0.75), angle=90, code=3, length=0.0, lwd=3, lty=1, col='goldenrod')
    arrows((n_pop_cast+1), quantile(df_1$mu_R_max, 0.025), (n_pop_cast+1), quantile(df_1$mu_R_max, 0.975), angle=90, code=3, length=0.01, lwd=1.5, lty=1, col='goldenrod')
    plot(n_pop_cast+1, median(df_1$sigma_R_max), col='goldenrod', ylim=c(quantile(df_1$sigma_R_max, 0.0), quantile(df_1$sigma_R_max, 0.995)), xaxt='n', yaxt='n',cex=1)
    #text(x=n_pop_cast, y=1.055, 'F', srt=90)
    arrows((n_pop_cast+1), quantile(df_1$sigma_R_max, 0.25), (n_pop_cast+1), quantile(df_1$sigma_R_max, 0.75), angle=90, code=3, length=0.0, lwd=3, lty=1, col='goldenrod')
    arrows((n_pop_cast+1), quantile(df_1$sigma_R_max, 0.025), (n_pop_cast+1), quantile(df_1$sigma_R_max, 0.975), angle=90, code=3, length=0.01, lwd=1.5, lty=1, col='goldenrod')
    axis(side=4)
    plot(1:n_pop_cast, colMedian(sigma_R), col='tan4', ylim=c(min(cred.2(sigma_R)[1,]),max(cred.2(sigma_R)[2,])), xlim=c(1,n_pop_cast), xaxt='n', yaxt='n',cex=1)
    axis(side=2, las=2)
    mtext(side=2, expression(sigma[R]), line=3.2, cex=1)
    #text(x=n_pop_cast, y=1.375, 'A', srt=90)
    
    arrows(1:n_pop_cast, cred.2(sigma_R)[1,], 1:n_pop_cast, cred.2(sigma_R)[2,], angle=90, code=3, length=0.01, lwd=1.5, lty=1, col='tan4')
    arrows(1:n_pop_cast, cred.3(sigma_R)[1,], 1:n_pop_cast, cred.3(sigma_R)[2,], angle=90, code=3, length=0.0, lwd=3, lty=1, col='tan4')
    axis(side=1, at=c(1:n_pop_cast), labels=unique(coho_dat_full$Name)[cast_pop], las=2, cex.axis=1.1)
    plot((n_pop_cast+1), median(exp(df_1$mu_sigma_R)), col='tan4', ylim=c(min(cred.2(sigma_R)[1,]),max(cred.2(sigma_R)[2,])), xaxt='n', yaxt='n', cex=1)
    mtext(side=1, expression(mu), line=1.0, las=3)
    #text(x=n_pop_cast, y=1.375, 'B', srt=90)
    arrows((n_pop_cast+1), quantile(exp(df_1$mu_sigma_R), 0.25), (n_pop_cast+1), quantile(exp(df_1$mu_sigma_R), 0.75), angle=90, code=3, length=0.0, lwd=3, lty=1, col='tan4')
    arrows((n_pop_cast+1), quantile(exp(df_1$mu_sigma_R), 0.025), (n_pop_cast+1), quantile(exp(df_1$mu_sigma_R), 0.975), angle=90, code=3, length=0.01, lwd=1.5, lty=1, col='tan4')
    axis(side=4, at=c(log(10), log(25), log(50), log(100), log(250), log(500), log(1000)), labels = c(10,25, 50, 100, 250, 500,1000))
    plot(n_pop_cast+1, median(df_1$sigma_sigma_R), col='tan4', ylim=c(quantile(df_1$sigma_sigma_R, 0.005), quantile(df_1$sigma_sigma_R, 0.99)), xaxt='n', yaxt='n', cex=1)
    #text(x=n_pop_cast, y=1.188, 'C', srt=90)
    mtext(side=1, expression(sigma), line=0.8, las=3)
    arrows((n_pop_cast+1), quantile(df_1$sigma_sigma_R, 0.25), (n_pop_cast+1), quantile(df_1$sigma_sigma_R, 0.75), angle=90, code=3, length=0.0, lwd=3, lty=1, col='tan4')
    arrows((n_pop_cast+1), quantile(df_1$sigma_sigma_R, 0.025), (n_pop_cast+1), quantile(df_1$sigma_sigma_R, 0.975), angle=90, code=3, length=0.01, lwd=1.5, lty=1, col='tan4')
    axis(side=4)
    dev.off()
    
    pdf(file='model_fits_AR.pdf')
    #Plot model estimates for smolt abundance
    Smolt_est <- exp(df_1[,grep('log_smolt_est', colnames(df_1))])
    par(mfrow=c(5,2), mar=c(2.5,4,2,1))
    pos <- 1
    J <- 0
    for(j in 1:n_pop){
      if(j %in% cast_pop){
        plot(1:n_year,  colMedian(Smolt_est[,pos:(pos+n_year-1)]), type='o', xlab= 'year', ylab='Smolt abundance', main=unique(coho_dat_full$Name)[j], col='red', lwd=1, pch=16, xaxt='n', ylim=c(min(cred.3(Smolt_est[,pos:(pos+n_year-1)])[1,]), max(cred.3(Smolt_est[,pos:(pos+n_year-1)])[2,])))
        axis(side=1, at=c(1:n_year), labels=c(1:n_year)+yr)
        xx <- c(1:n_year,n_year:1)
        yy <- c(cred.2(Smolt_est[,pos:(pos+n_year-1)])[1,], rev(cred.2(Smolt_est[,pos:(pos+n_year-1)])[2,]))
        polygon(x=xx, y=yy, density = -1, border = NA, col=rgb(238, 64, 0, max = 255, alpha = 90))
        yy <- c(cred.3(Smolt_est[,pos:(pos+n_year-1)])[1,], rev(cred.3(Smolt_est[,pos:(pos+n_year-1)])[2,]))
        polygon(x=xx, y=yy, density = -1, border = NA, col=rgb(238, 64, 0, max = 255, alpha = 150))
        if(j %in% pop_smolt){
          #J <- J + 1
          #points(smolt_true[slice_smolt_start[J]:slice_smolt_end[J]],  smolt_dat[slice_smolt_start[J]:slice_smolt_end[J]], col='black', pch=16)
          points(coho_dat_full$yr[coho_dat_full$pop==j], coho_dat_full$Smolt.Abundance[coho_dat_full$pop==j], col='black', pch=1, lwd=1.25)

        }
      }
      pos <- pos + n_year
    }
    
    err_fun <- function(m, s){
      sig <- -2*log(m) + log(s^2 + m^2)
      sig <- -2*log(m) + log(s^2 + m^2)
      
      return(sqrt(sig))
    }
    
    #Plot model estimates for Escapement abundance
    esc_est <- exp(df_1[,grep('log_adult_est', colnames(df_1))])
    par(mfrow=c(5,2))
    pos <- 1
    J <- 0
    for(j in 1:n_pop){
      if(j %in% cast_pop){
        plot(1:n_year,  colMedian(esc_est[,pos:(pos+n_year-1)]), type='o', xlab= 'year', ylab='escapement abundance', main=unique(coho_dat_full$Name)[j], col='red', lwd=1, pch=16, xaxt='n', ylim=c(min(cred.3(esc_est[,pos:(pos+n_year-1)])[1,]), max(cred.3(esc_est[,pos:(pos+n_year-1)])[2,])))
        axis(side=1, at=c(1:n_year), labels=c(1:n_year)+yr)
        xx <- c(1:n_year,n_year:1)
        yy <- c(cred.2(esc_est[,pos:(pos+n_year-1)])[1,], rev(cred.2(esc_est[,pos:(pos+n_year-1)])[2,]))
        polygon(x=xx, y=yy, density = -1, border = NA, col=rgb(238, 64, 0, max = 255, alpha = 90))
        yy <- c(cred.3(esc_est[,pos:(pos+n_year-1)])[1,], rev(cred.3(esc_est[,pos:(pos+n_year-1)])[2,]))
        polygon(x=xx, y=yy, density = -1, border = NA, col=rgb(238, 64, 0, max = 255, alpha = 150))
        if(j %in% pop_esc){
          #J <- J + 1
          #points(esc_true[slice_esc_start[J]:slice_esc_end[J]],  esc_dat[slice_esc_start[J]:slice_esc_end[J]], col='black', pch=16)
          points(coho_dat_full$yr[coho_dat_full$pop==j], coho_dat_full$Spawners[coho_dat_full$pop==j], col='black', pch=1, lwd=1.25)
        }
      }
      pos <- pos + n_year
    }
    
    #Plot model estimates for Harvest abundance
    harvest_est <- exp(df_1[,grep('log_harvest_est', colnames(df_1))])
    par(mfrow=c(5,2))
    pos <- 1
    J <- 0
    for(j in 1:n_pop){
      if(j %in% cast_pop){
        plot(1:n_year,  colMedian(harvest_est[,pos:(pos+n_year-1)]), type='o', xlab= 'year', ylab='harvest abundance', main=unique(coho_dat_full$Name)[j], col='red', lwd=1, pch=16, xaxt='n', ylim=c(min(cred.3(harvest_est[,pos:(pos+n_year-1)])[1,]), max(cred.3(harvest_est[,pos:(pos+n_year-1)])[2,])))
        axis(side=1, at=c(1:n_year), labels=c(1:n_year)+yr)
        xx <- c(1:n_year,n_year:1)
        yy <- c(cred.2(harvest_est[,pos:(pos+n_year-1)])[1,], rev(cred.2(harvest_est[,pos:(pos+n_year-1)])[2,]))
        polygon(x=xx, y=yy, density = -1, border = NA, col=rgb(238, 64, 0, max = 255, alpha = 90))
        yy <- c(cred.3(harvest_est[,pos:(pos+n_year-1)])[1,], rev(cred.3(harvest_est[,pos:(pos+n_year-1)])[2,]))
        polygon(x=xx, y=yy, density = -1, border = NA, col=rgb(238, 64, 0, max = 255, alpha = 150))
        if(j %in% pop_catch){
         # J <- J + 1
          #points(harvest_true[slice_harvest_start[J]:slice_harvest_end[J]],  harvest_dat[slice_harvest_start[j]:slice_harvest_end[J]], col='black', pch=16)
          points(coho_dat_full$yr[coho_dat_full$pop==j], coho_dat_full$Harvest[coho_dat_full$pop==j], col='black', pch=1, lwd=1.25)
          
        }
      }
      pos <- pos + n_year
    } 
    
    #Plot model estimates for marine survival
    MS_est <-apply(df_1[,grep('logit_smolt_survival_pop', colnames(df_1))],2 , inv.logit)
    hatch_offset <- median(df_1$hatch_offset)
    par(mfrow=c(5,2))
    pos <- 1
    J <- 0
    for(j in 1:n_pop){
      if(j %in% cast_pop){
        plot(1:n_year,  colMedian(MS_est[,pos:(pos+n_year-1)]), type='o', xlab= 'year', ylab='Marine survival', main=unique(coho_dat_full$Name)[j], col='red', lwd=1, pch=16, xaxt='n', ylim=c(min(cred.3(MS_est[,pos:(pos+n_year-1)])[1,]), max(cred.3(MS_est[,pos:(pos+n_year-1)])[2,])))
        axis(side=1, at=c(1:n_year), labels=c(1:n_year)+yr)
        xx <- c(1:n_year,n_year:1)
        yy <- c(cred.2(MS_est[,pos:(pos+n_year-1)])[1,], rev(cred.2(MS_est[,pos:(pos+n_year-1)])[2,]))
        polygon(x=xx, y=yy, density = -1, border = NA, col=rgb(238, 64, 0, max = 255, alpha = 90))
        yy <- c(cred.3(MS_est[,pos:(pos+n_year-1)])[1,],rev(cred.3(MS_est[,pos:(pos+n_year-1)])[2,]))
        polygon(x=xx, y=yy, density = -1, border = NA, col=rgb(238, 64, 0, max = 255, alpha = 150))
        if(j %in% pop_MS){
          #J <- J + 1
          #points(MS_true[slice_MS_start[J]:slice_MS_end[J]],  (MS_dat_x/MS_dat_N)[slice_MS_start[J]:slice_MS_end[J]], col='black', pch=16)
          points(coho_dat_full$yr[coho_dat_full$pop==j], coho_dat_full$Fishery_Plus_Escapement[coho_dat_full$pop==j]/coho_dat_full$Release_No[coho_dat_full$pop==j] , col='black', pch=1, lwd=1.25)
        }
      }
      pos <- pos + n_year
    }   
    
    #Plot adjusted marine survival estimated (with hatchery offset for populations whose marine survival data is from a hatchery)
    MS_est_adj <- df_1[,grep('smolt_survival_adj', colnames(df_1))]
    par(mfrow=c(5,2))
    pos <- 1
    for(j in 1:n_pop_MS){
      if(j %in% cast_pop){
        plot(1:n_year,  (colMedian(MS_est_adj[,pos:(pos+n_year-1)])), type='o', xlab= 'year', ylab='Adjusted marine survival', main=unique(coho_dat_full$Name)[pop_MS[j]], col='red', lwd=1, pch=16, ylim=c(0,0.25), xaxt='n')
        axis(side=1, at=c(1:n_year), labels=c(1:n_year)+yr)
        xx <- c(1:n_year,n_year:1)
        yy <- c((cred.2(MS_est_adj[,pos:(pos+n_year-1)])[1,]), (rev(cred.2(MS_est_adj[,pos:(pos+n_year-1)])[2,])))
        polygon(x=xx, y=yy, density = -1, border = NA, col=rgb(238, 64, 0, max = 255, alpha = 90))
        yy <- c((cred.3(MS_est_adj[,pos:(pos+n_year-1)])[1,]), (rev(cred.3(MS_est_adj[,pos:(pos+n_year-1)])[2,])))
        polygon(x=xx, y=yy, density = -1, border = NA, col=rgb(238, 64, 0, max = 255, alpha = 150))
        #points(MS_true[slice_MS_start[j]:slice_MS_end[j]],  (MS_dat_x/MS_dat_N)[slice_MS_start[j]:slice_MS_end[j]], col='black', pch=16)
        points(coho_dat_full$yr[coho_dat_full$pop==j], coho_dat_full$Fishery_Plus_Escapement[coho_dat_full$pop==j]/coho_dat_full$Release_No[coho_dat_full$pop==j] , col='black', pch=1, lwd=1.25)
      }
      pos <- pos + n_year
    }  
    dev.off()
    
    pdf('state_residuals_AR.pdf')
    #Plot time-series of recruitment deviations with temporal autocorrelation
    r_dev <- df_1[,grep('r_dev', colnames(df_1))]
    par(mfrow=c(4,2))
    pos <- 1
    for(j in 1:n_pop){
      plot(1:n_year,  colMedian(r_dev[,pos:(pos+n_year-1)]), type='o', xlab= 'year', ylab='Recruitment deviations', main=unique(coho_dat_full$Name)[j], col='red', lwd=1, pch=16, ylim=c(min(cred.2(r_dev[,pos:(pos+n_year-1)])[1,]), max(cred.2(r_dev[,pos:(pos+n_year-1)])[2,])), xaxs='i', xaxt='n')
      axis(side=1, at=c(1:n_year), labels=c(1:n_year)+yr)
      xx <- c(1:n_year,n_year:1)
      yy <- c(cred.2(r_dev[,pos:(pos+n_year-1)])[1,], rev(cred.2(r_dev[,pos:(pos+n_year-1)])[2,]))
      polygon(x=xx, y=yy, density = -1, border = NA, col=rgb(238, 64, 0, max = 255, alpha = 90))
      yy <- c(cred.3(r_dev[,pos:(pos+n_year-1)])[1,], rev(cred.3(r_dev[,pos:(pos+n_year-1)])[2,]))
      polygon(x=xx, y=yy, density = -1, border = NA, col=rgb(238, 64, 0, max = 255, alpha = 150))
      abline(h=0, lty=2, col='black')
      acf(colMedian(r_dev[,pos:(pos+n_year-1)]), main='Temporal autocorrelation')
      pos <- pos + n_year
    }
    
    #Plot time-series of marine survival deviations with temporal autocorrelation
    N_year <- n_year -1 
    MS_dev <- df_1[,grep('MS_dev_pop', colnames(df_1))]
    #MS_dev_mean <- df_1[,grep('bias', colnames(df_1))][,(n_pop+2):(n_pop*2+1)]
    
    par(mfrow=c(4,2))
    pos <- 1
    for(j in 1:n_pop){
      plot(2:n_year,  colMedian(MS_dev[,pos:(pos+N_year-1)]), type='o', xlab= 'year', ylab='Marine survival deviations', main=unique(coho_dat_full$Name)[j], col='red', lwd=1, pch=16, ylim=c(min(cred.2(MS_dev[,pos:(pos+N_year-1)])[1,]), max(cred.2(MS_dev[,pos:(pos+N_year-1)])[2,])), xaxs='i', xaxt='n')
      axis(side=1, at=c(2:n_year), labels=c(2:n_year)+yr)
      xx <- c(2:n_year,n_year:2)
      yy <- c(cred.2(MS_dev[,pos:(pos+N_year-1)])[1,], rev(cred.2(MS_dev[,pos:(pos+N_year-1)])[2,]))
      polygon(x=xx, y=yy, density = -1, border = NA, col=rgb(238, 64, 0, max = 255, alpha = 90))
      yy <- c(cred.3(MS_dev[,pos:(pos+N_year-1)])[1,], rev(cred.3(MS_dev[,pos:(pos+N_year-1)])[2,]))
      polygon(x=xx, y=yy, density = -1, border = NA, col=rgb(238, 64, 0, max = 255, alpha = 150))
      abline(h=0, lty=2, col='black')
      #abline(h=colMedian(MS_dev_mean)[j], lty=2, col='blue')
      acf(colMedian(MS_dev[,pos:(pos+N_year-1)]), main='Temporal autocorrelation')
      pos <- pos + N_year
    }
    dev.off()
    
    pdf(file='model_residuals_AR.pdf')
    par(mfrow=c(5,2), mar=c(2.5,4,2,1))
    smolt_resid <- df_1[,grep('smolt_resid', colnames(df_1))]
    esc_resid <- df_1[,grep('esc_resid', colnames(df_1))]
    harvest_resid <- df_1[,grep('harvest_resid', colnames(df_1))]
    MS_resid <- df_1[,grep('MS_resid', colnames(df_1))]
    
    for(j in 1:n_pop_smolt){
      plot(1:n_year + yr, type='n', ylim=c(min(cred.2(smolt_resid[slice_smolt_start[j]:slice_smolt_end[j]])[1,]), max(cred.2(smolt_resid[slice_smolt_start[j]:slice_smolt_end[j]])[2,])), ylab='Smolt residuals', xaxt='n', xlab='year', main=unique(coho_dat_full$Name)[pop_smolt[j]])
      points(smolt_true[slice_smolt_start[j]:slice_smolt_end[j]], colMedian(smolt_resid[slice_smolt_start[j]:slice_smolt_end[j]]), pch=16, col='red')
      arrows(smolt_true[slice_smolt_start[j]:slice_smolt_end[j]], cred.2(smolt_resid[slice_smolt_start[j]:slice_smolt_end[j]])[1,], smolt_true[slice_smolt_start[j]:slice_smolt_end[j]], cred.2(smolt_resid[slice_smolt_start[j]:slice_smolt_end[j]])[2,], angle=90, code=3, length=0.01, lwd=0.5, lty=1, col='red')
      arrows(smolt_true[slice_smolt_start[j]:slice_smolt_end[j]], cred.3(smolt_resid[slice_smolt_start[j]:slice_smolt_end[j]])[1,], smolt_true[slice_smolt_start[j]:slice_smolt_end[j]], cred.3(smolt_resid[slice_smolt_start[j]:slice_smolt_end[j]])[2,], angle=90, code=3, length=0.0, lwd=2, lty=1, col='red')
      abline(h=0, lty=2)
      axis(side=1, at=c(1:n_year), labels=c(1:n_year)+yr)
    }
    par(mfrow=c(5,2), mar=c(2.5,4,2,1))
    for(j in 1:n_pop_esc){
      plot(1:n_year + yr, type='n', ylim=c(min(cred.2(esc_resid[slice_esc_start[j]:slice_esc_end[j]])[1,]), max(cred.2(esc_resid[slice_esc_start[j]:slice_esc_end[j]])[2,])), ylab='Escapement residuals', xaxt='n', xlab='year', main=unique(coho_dat_full$Name)[pop_esc[j]])
      points(esc_true[slice_esc_start[j]:slice_esc_end[j]], colMedian(esc_resid[slice_esc_start[j]:slice_esc_end[j]]), pch=16, col='red')
      arrows(esc_true[slice_esc_start[j]:slice_esc_end[j]], cred.2(esc_resid[slice_esc_start[j]:slice_esc_end[j]])[1,], esc_true[slice_esc_start[j]:slice_esc_end[j]], cred.2(esc_resid[slice_esc_start[j]:slice_esc_end[j]])[2,], angle=90, code=3, length=0.01, lwd=0.5, lty=1, col='red')
      arrows(esc_true[slice_esc_start[j]:slice_esc_end[j]], cred.3(esc_resid[slice_esc_start[j]:slice_esc_end[j]])[1,], esc_true[slice_esc_start[j]:slice_esc_end[j]], cred.3(esc_resid[slice_esc_start[j]:slice_esc_end[j]])[2,], angle=90, code=3, length=0.0, lwd=2, lty=1, col='red')
      abline(h=0, lty=2)
      axis(side=1, at=c(1:n_year), labels=c(1:n_year)+yr)
    }
    par(mfrow=c(5,2), mar=c(2.5,4,2,1))
    for(j in 1:n_pop_catch){
      plot(1:n_year + yr, type='n', ylim=c(min(cred.2(harvest_resid[slice_harvest_start[j]:slice_harvest_end[j]])[1,]), max(cred.2(harvest_resid[slice_harvest_start[j]:slice_harvest_end[j]])[2,])), ylab='Harvest residuals', xaxt='n', xlab='year', main=unique(coho_dat_full$Name)[pop_catch[j]])
      points(harvest_true[slice_harvest_start[j]:slice_harvest_end[j]], colMedian(harvest_resid[slice_harvest_start[j]:slice_harvest_end[j]]), pch=16, col='red')
      arrows(harvest_true[slice_harvest_start[j]:slice_harvest_end[j]], cred.2(harvest_resid[slice_harvest_start[j]:slice_harvest_end[j]])[1,], harvest_true[slice_harvest_start[j]:slice_harvest_end[j]], cred.2(harvest_resid[slice_harvest_start[j]:slice_harvest_end[j]])[2,], angle=90, code=3, length=0.01, lwd=0.5, lty=1, col='red')
      arrows(harvest_true[slice_harvest_start[j]:slice_harvest_end[j]], cred.3(harvest_resid[slice_harvest_start[j]:slice_harvest_end[j]])[1,], harvest_true[slice_harvest_start[j]:slice_harvest_end[j]], cred.3(harvest_resid[slice_harvest_start[j]:slice_harvest_end[j]])[2,], angle=90, code=3, length=0.0, lwd=2, lty=1, col='red')
      abline(h=0, lty=2)
      axis(side=1, at=c(1:n_year), labels=c(1:n_year)+yr)
    }
    par(mfrow=c(5,2), mar=c(2.5,4,2,1))
    for(j in 1:n_pop_MS){
      plot(1:n_year + yr, type='n', ylim=c(min(cred.2(MS_resid[slice_MS_start[j]:slice_MS_end[j]])[1,]), max(cred.2(MS_resid[slice_MS_start[j]:slice_MS_end[j]])[2,])), ylab='MS residuals', xaxt='n', xlab='year', main=unique(coho_dat_full$Name)[pop_MS[j]])
      points(MS_true[slice_MS_start[j]:slice_MS_end[j]], colMedian(MS_resid[slice_MS_start[j]:slice_MS_end[j]]), pch=16, col='red')
      arrows(MS_true[slice_MS_start[j]:slice_MS_end[j]], cred.2(MS_resid[slice_MS_start[j]:slice_MS_end[j]])[1,], MS_true[slice_MS_start[j]:slice_MS_end[j]], cred.2(MS_resid[slice_MS_start[j]:slice_MS_end[j]])[2,], angle=90, code=3, length=0.01, lwd=0.5, lty=1, col='red')
      arrows(MS_true[slice_MS_start[j]:slice_MS_end[j]], cred.3(MS_resid[slice_MS_start[j]:slice_MS_end[j]])[1,], MS_true[slice_MS_start[j]:slice_MS_end[j]], cred.3(MS_resid[slice_MS_start[j]:slice_MS_end[j]])[2,], angle=90, code=3, length=0.0, lwd=2, lty=1, col='red')
      abline(h=0, lty=2)
      axis(side=1, at=c(1:n_year), labels=c(1:n_year)+yr)
    }
    dev.off()
    
    pdf(file='posterior_predictive_check_AR.pdf')
    #Posterior Predictive check for (1) smolt abundance data, (2) escapement abundance data, (3) harvest abundance data
    Smolt_PPC <- df_1[,grep('PPC_smolt', colnames(df_1))]
    Esc_PPC <- df_1[,grep('PPC_esc', colnames(df_1))]
    Harvest_PPC <- df_1[,grep('PPC_catch', colnames(df_1))]
    MS_PPC <- df_1[,grep('PPC_MS', colnames(df_1))]
    n_smolt <- length(smolt_dat)
    n_esc <- length(esc_dat)
    n_harvest <- length(harvest_dat)
    n_MS <- length(MS_dat_x)
    #(1) Smolt adbundance PPC
    par(mfrow=c(5,2), mar=c(2.5,4,2,1))
    PPC_95 <- 0
    PPC_50 <- 0
    pos <- 1
    j <- 1
    for(j in 1:n_pop_smolt){  
      plot(1:n_year,  colMedian(Smolt_PPC[,pos:(pos+n_year-1)]), type='o', xlab= 'year', ylab='smolt PPC', main=unique(coho_dat_full$Name)[pop_smolt[j]], col='red', lwd=1, pch=16, xaxt='n', ylim=c(min(cred.2(Smolt_PPC[,pos:(pos+n_year-1)])[1,]), max(cred.2(Smolt_PPC[,pos:(pos+n_year-1)])[2,])))
      axis(side=1, at=c(1:n_year), labels=c(1:n_year)+yr)
      xx <- c(1:n_year,n_year:1)
      yy <- c(cred.2(Smolt_PPC[,pos:(pos+n_year-1)])[1,], rev(cred.2(Smolt_PPC[,pos:(pos+n_year-1)])[2,]))
      polygon(x=xx, y=yy, density = -1, border = NA, col=rgb(238, 64, 0, max = 255, alpha = 90))
      yy <- c(cred.3(Smolt_PPC[,pos:(pos+n_year-1)])[1,], rev(cred.3(Smolt_PPC[,pos:(pos+n_year-1)])[2,]))
      polygon(x=xx, y=yy, density = -1, border = NA, col=rgb(238, 64, 0, max = 255, alpha = 150))
      points(smolt_true[slice_smolt_start[j]:slice_smolt_end[j]],  smolt_dat[slice_smolt_start[j]:slice_smolt_end[j]], col='black', pch=1, lwd=1.25)
      
      PPC_95 <- PPC_95+length(which(smolt_dat[slice_smolt_start[j]:slice_smolt_end[j]] > cred.2(Smolt_PPC[,(pos:(pos+n_year-1))[smolt_true[slice_smolt_start[j]:slice_smolt_end[j]]]])[1,] &  smolt_dat[slice_smolt_start[j]:slice_smolt_end[j]] < cred.2(Smolt_PPC[,(pos:(pos+n_year-1))[smolt_true[slice_smolt_start[j]:slice_smolt_end[j]]]])[2,]))
      PPC_50 <- PPC_50+length(which(smolt_dat[slice_smolt_start[j]:slice_smolt_end[j]] > cred.3(Smolt_PPC[,(pos:(pos+n_year-1))[smolt_true[slice_smolt_start[j]:slice_smolt_end[j]]]])[1,] &  smolt_dat[slice_smolt_start[j]:slice_smolt_end[j]] < cred.3(Smolt_PPC[,(pos:(pos+n_year-1))[smolt_true[slice_smolt_start[j]:slice_smolt_end[j]]]])[2,]))
      pos <- pos + n_year
    }
    print(PPC_95/length(smolt_true))
    print(PPC_50/length(smolt_true))
    
    #(2) Escapement abundance PPC
    par(mfrow=c(5,2), mar=c(2.5,4,2,1))
    PPC_95 <- 0
    PPC_50 <- 0
    pos <- 1
    j <- 1
    for(j in 1:n_pop_esc){  
      plot(1:n_year,  colMedian(Esc_PPC[,pos:(pos+n_year-1)]), type='o', xlab= 'year', ylab='Esc PPC', main=unique(coho_dat_full$Name)[pop_esc[j]], col='red', lwd=1, pch=16, xaxt='n',ylim=c(min(cred.2(Esc_PPC[,pos:(pos+n_year-1)])[1,]), max(cred.2(Esc_PPC[,pos:(pos+n_year-1)])[2,])))
      axis(side=1, at=c(1:n_year), labels=c(1:n_year)+yr)
      xx <- c(1:n_year,n_year:1)
      yy <- c(cred.2(Esc_PPC[,pos:(pos+n_year-1)])[1,], rev(cred.2(Esc_PPC[,pos:(pos+n_year-1)])[2,]))
      polygon(x=xx, y=yy, density = -1, border = NA, col=rgb(238, 64, 0, max = 255, alpha = 90))
      yy <- c(cred.3(Esc_PPC[,pos:(pos+n_year-1)])[1,], rev(cred.3(Esc_PPC[,pos:(pos+n_year-1)])[2,]))
      polygon(x=xx, y=yy, density = -1, border = NA, col=rgb(238, 64, 0, max = 255, alpha = 150))
      points(esc_true[slice_esc_start[j]:slice_esc_end[j]],  esc_dat[slice_esc_start[j]:slice_esc_end[j]], col='black', pch=1, lwd=1.25)
      
      PPC_95 <- PPC_95+length(which(esc_dat[slice_esc_start[j]:slice_esc_end[j]] > cred.2(Esc_PPC[,(pos:(pos+n_year-1))[esc_true[slice_esc_start[j]:slice_esc_end[j]]]])[1,] &  esc_dat[slice_esc_start[j]:slice_esc_end[j]] < cred.2(Esc_PPC[,(pos:(pos+n_year-1))[esc_true[slice_esc_start[j]:slice_esc_end[j]]]])[2,]))
      PPC_50 <- PPC_50+length(which(esc_dat[slice_esc_start[j]:slice_esc_end[j]] > cred.3(Esc_PPC[,(pos:(pos+n_year-1))[esc_true[slice_esc_start[j]:slice_esc_end[j]]]])[1,] &  esc_dat[slice_esc_start[j]:slice_esc_end[j]] < cred.3(Esc_PPC[,(pos:(pos+n_year-1))[esc_true[slice_esc_start[j]:slice_esc_end[j]]]])[2,]))
      pos <- pos + n_year
    }
    print(PPC_95/length(esc_true))
    print(PPC_50/length(esc_true))
    
    #(3) Harvest adbundance PPC
    par(mfrow=c(5,2), mar=c(2.5,4,2,1))
    PPC_95 <- 0
    PPC_50 <- 0
    pos <- 1
    for(j in 1:n_pop_catch){  
      plot(1:n_year,  colMedian(Harvest_PPC[,pos:(pos+n_year-1)]), type='o', xlab= 'year', ylab='Catch PPC', main=unique(coho_dat_full$Name)[pop_catch[j]], col='red', lwd=1, pch=16, xaxt='n',ylim=c(min(cred.2(Harvest_PPC[,pos:(pos+n_year-1)])[1,]), max(cred.2(Harvest_PPC[,pos:(pos+n_year-1)])[2,])))
      axis(side=1, at=c(1:n_year), labels=c(1:n_year)+yr)
      xx <- c(1:n_year,n_year:1)
      yy <- c(cred.2(Harvest_PPC[,pos:(pos+n_year-1)])[1,], rev(cred.2(Harvest_PPC[,pos:(pos+n_year-1)])[2,]))
      polygon(x=xx, y=yy, density = -1, border = NA, col=rgb(238, 64, 0, max = 255, alpha = 90))
      yy <- c(cred.3(Harvest_PPC[,pos:(pos+n_year-1)])[1,], rev(cred.3(Harvest_PPC[,pos:(pos+n_year-1)])[2,]))
      polygon(x=xx, y=yy, density = -1, border = NA, col=rgb(238, 64, 0, max = 255, alpha = 150))
      points(harvest_true[slice_harvest_start[j]:slice_harvest_end[j]],  harvest_dat[slice_harvest_start[j]:slice_harvest_end[j]], col='black', pch=1, lwd=1.25)
      
      PPC_95 <- PPC_95+length(which(harvest_dat[slice_harvest_start[j]:slice_harvest_end[j]] > cred.2(Harvest_PPC[,(pos:(pos+n_year-1))[harvest_true[slice_harvest_start[j]:slice_harvest_end[j]]]])[1,] &  harvest_dat[slice_harvest_start[j]:slice_harvest_end[j]] < cred.2(Harvest_PPC[,(pos:(pos+n_year-1))[harvest_true[slice_harvest_start[j]:slice_harvest_end[j]]]])[2,]))
      PPC_50 <- PPC_50+length(which(harvest_dat[slice_harvest_start[j]:slice_harvest_end[j]] > cred.3(Harvest_PPC[,(pos:(pos+n_year-1))[harvest_true[slice_harvest_start[j]:slice_harvest_end[j]]]])[1,] &  harvest_dat[slice_harvest_start[j]:slice_harvest_end[j]] < cred.3(Harvest_PPC[,(pos:(pos+n_year-1))[harvest_true[slice_harvest_start[j]:slice_harvest_end[j]]]])[2,]))
      pos <- pos + n_year
    }
    print(PPC_95/length(harvest_true))
    print(PPC_50/length(harvest_true))
    
    #(4) Marine survival PPC
    par(mfrow=c(5,2), mar=c(2.5,4,2,1))
    PPC_95 <- 0
    PPC_50 <- 0
    for(j in 1:n_pop_MS){  
      plot(colMedian(MS_PPC[slice_MS_start[j]:slice_MS_end[j]]), xlab= 'year', ylab='MS PPC', main=unique(coho_dat_full$Name)[pop_MS[j]], col='red', lwd=1, pch=16, ylim=c(min(colMedian(MS_PPC[slice_MS_start[j]:slice_MS_end[j]]))/1000, max(colMedian(MS_PPC[slice_MS_start[j]:slice_MS_end[j]]))*4))
      arrows(1:length(colMedian(MS_PPC[slice_MS_start[j]:slice_MS_end[j]])), cred.3(MS_PPC[slice_MS_start[j]:slice_MS_end[j]])[1,], 1:length(colMedian(MS_PPC[slice_MS_start[j]:slice_MS_end[j]])), cred.3(MS_PPC[slice_MS_start[j]:slice_MS_end[j]])[2,], angle=90, code=3, length=0.0, lwd=2.0, lty=1, col='red')
      arrows(1:length(colMedian(MS_PPC[slice_MS_start[j]:slice_MS_end[j]])), cred.2(MS_PPC[slice_MS_start[j]:slice_MS_end[j]])[1,], 1:length(colMedian(MS_PPC[slice_MS_start[j]:slice_MS_end[j]])), cred.2(MS_PPC[slice_MS_start[j]:slice_MS_end[j]])[2,], angle=90, code=3, length=0.01, lwd=0.75, lty=1, col='red')
      points(MS_dat_x[slice_MS_start[j]:slice_MS_end[j]], col='black', pch=1, lwd=1.25)
      
      PPC_50 <- PPC_50 + length(which(MS_dat_x[slice_MS_start[j]:slice_MS_end[j]] >  cred.3(MS_PPC[slice_MS_start[j]:slice_MS_end[j]])[1,] & MS_dat_x[slice_MS_start[j]:slice_MS_end[j]] < cred.3(MS_PPC[slice_MS_start[j]:slice_MS_end[j]])[2,]))
      PPC_95 <- PPC_95 +length(which(MS_dat_x[slice_MS_start[j]:slice_MS_end[j]] >  cred.2(MS_PPC[slice_MS_start[j]:slice_MS_end[j]])[1,] & MS_dat_x[slice_MS_start[j]:slice_MS_end[j]] < cred.2(MS_PPC[slice_MS_start[j]:slice_MS_end[j]])[2,]))
    }
    print(PPC_95/n_MS)
    print(PPC_50/n_MS)
    
    dev.off()
    
    #Plot the assumed shape of the stock recruit function for each population
    alpha <- exp(df_1[,grep('log_alpha', colnames(df_1))])
    R_max <- exp(df_1[,grep('log_R_max', colnames(df_1))])
    stream_dist <- aggregate(coho_dat_full$KM, by=list('Population'=coho_dat_full$Population), mean, na.rm=TRUE)$x
    pdf('S_R_plots_AR.pdf')
    par(mfrow=c(6,6), mar=c(2,1.3,0.5,0.5), oma=c(2,2,0.6,0.25), mgp=c(3,0.5,0))
    rec_func <- function(S=10000, alpha=60, R_max=0.00017, sigma_R=0.0){
      #calculate recruitment deviation using mean juvenile survival and lognormal bias correction
      epsilon <- exp((rnorm(1, -((sigma_R^2)/2), sigma_R)))
      R <- (S/((1/alpha)+(S/R_max)))*epsilon
      return(R)
    }
    J_esc <- 0
    J_smolt <- 0
    pos <- 1
    for(j in 1:n_pop){
      #if(j %in% pop_smolt){
      if(j %in% cast_pop){
        
        n_spawners = 0:max(ceiling(coho_dat_full$Spawners[coho_dat_full$pop==j]/stream_dist[j]), na.rm=TRUE)
        
        n_recruits <- vector(length=length(n_spawners))
        n_recruits_lower_95 <- vector(length=length(n_spawners))
        n_recruits_upper_95 <- vector(length=length(n_spawners))
        n_recruits_lower_50 <- vector(length=length(n_spawners))
        n_recruits_upper_50 <- vector(length=length(n_spawners))
        
        for(s in 1:length(n_spawners)){
          n_recruits[s] <- rec_func(S=n_spawners[s], alpha=colMedian(alpha)[j], R_max=colMedian(R_max)[j], sigma_R=0.0)
          n_recruits_lower_95[s] <- rec_func(S=n_spawners[s], alpha=cred.2(alpha)[1,j], R_max=cred.2(R_max)[1,j], sigma_R=0.0)
          n_recruits_upper_95[s] <- rec_func(S=n_spawners[s], alpha=cred.2(alpha)[2,j], R_max=cred.2(R_max)[2,j], sigma_R=0.0)
          n_recruits_lower_50[s] <- rec_func(S=n_spawners[s], alpha=cred.3(alpha)[1,j], R_max=cred.3(R_max)[1,j], sigma_R=0.0)
          n_recruits_upper_50[s] <- rec_func(S=n_spawners[s], alpha=cred.3(alpha)[2,j], R_max=cred.3(R_max)[2,j], sigma_R=0.0)
        }

            plot(n_spawners, n_recruits, xaxs='i',type='l', xlim=c(0,max(n_spawners)), cex.main=0.825, cex=0.85, axes=F, lwd=1,
                 ylim=c(0, max(cred.4(Smolt_est[,pos:(pos+n_year-1)]/stream_dist[j])[2,],  1.05*max(n_recruits_upper_95))),
                 col=rgb(4,43,47, max = 255, alpha = 255), xlab=NA, ylab=NA)
            box(col='grey', lwd=0.5)
          axis(side=1, cex.axis=0.85, lwd=0.5)
          axis(side=2, cex.axis=0.85, lwd=0.5)
          mtext(side=3, line=0.1, unique(coho_dat_full$Name)[j], cex=0.6)
          
          xx <- c(n_spawners,rev(n_spawners))
          yy <- c(n_recruits_lower_95, rev(n_recruits_upper_95))
          polygon(x=xx, y=yy, density = -1, border = 'cadetblue4', lwd=0.25, col=rgb(95,158,160, max = 255, alpha = 80))
          yy <- c(n_recruits_lower_50, rev(n_recruits_upper_50))
          polygon(x=xx, y=yy, density = -1, border = 'cadetblue4', lwd=0.25, col=rgb(95,158,160, max = 255, alpha = 140))
          
          points(colMedian(esc_est[,pos:(pos+n_year-1)]/stream_dist[j]),  colMedian(Smolt_est[,pos:(pos+n_year-1)]/stream_dist[j]), pch=16, cex=0.85, col=rgb(7,94,94, max = 255, alpha = 125))
          points(colMedian(esc_est[,pos:(pos+n_year-1)]/stream_dist[j]),  colMedian(Smolt_est[,pos:(pos+n_year-1)]/stream_dist[j]), pch=1, cex=0.85, col=rgb(7,94,94, max = 255, alpha = 125), lwd=0.5)
          #arrows(x0=cred.4(esc_est[,pos:(pos+n_year-1)]/stream_dist[j])[1,], y0=colMedian(Smolt_est[,pos:(pos+n_year-1)]/stream_dist[j]), x1=cred.4(esc_est[,pos:(pos+n_year-1)]/stream_dist[j])[2,], 
                # length=0.01, angle=90, code=3, col=rgb(95,158,160, max = 255, alpha = 175), lwd=0.35)
         # arrows(x0=colMedian(esc_est[,pos:(pos+n_year-1)]/stream_dist[j]), y0=cred.4(Smolt_est[,pos:(pos+n_year-1)]/stream_dist[j])[1,], y1=cred.4(Smolt_est[,pos:(pos+n_year-1)]/stream_dist[j])[2,], 
                 #length=0.01, angle=90, code=3, col=rgb(95,158,160, max = 255, alpha = 175), lwd=0.35)
          
          if(j %in% pop_smolt & j %in% pop_esc){
            x_1 <- cbind(coho_dat_full$Spawners[coho_dat_full$pop==j], coho_dat_full$Smolt.Abundance[coho_dat_full$pop==j])/stream_dist[j]
            x_2 <- na.omit(x_1)
            points(x_2[,1], x_2[,2], pch=1, col='black', cex=0.85, lwd=0.5)
            points(x_2[,1], x_2[,2], pch=16, col=rgb(49,56,56, max = 255, alpha = 175), cex=0.85)
          }
          if(j==cast_pop[33]){
            mtext(side=1, 'Spawners/km',line=1.75, cex=0.9, at=220)
          }
          if(j==cast_pop[19]){
            mtext(side=2, 'Smolts/km', line=1.75, cex=0.9, at=1450)
          }
      }
      pos <- pos + n_year
      }
    dev.off()
    
    #plot posteriors and priors of key model parameters
    pdf(file='priors_posteriors_AR.pdf')
    par(mfrow=c(3,2), mar=c(2,2,2,2), oma=c(2,2,1,1))
    n_int <- length(df_1$mu_alpha)
    
    #Stock-recruit hyperparameters
    mu_alpha <- df_1$mu_alpha
    hist(mu_alpha, col=rgb(95,158,160, max = 255, alpha = 150))
    hist(rnorm(n_int , 4.27, 2), add=T, col=rgb(238, 64, 0, max = 255, alpha = 150))
    
    sigma_alpha <- df_1$sigma_alpha
    hist(sigma_alpha, col=rgb(95,158,160, max = 255, alpha = 150))
    hist(rnorm(n_int , 0.43, 0.25), add=T, col=rgb(238, 64, 0, max = 255, alpha = 150))
    
    mu_R_max <- df_1$mu_R_max
    hist(mu_R_max, col=rgb(95,158,160, max = 255, alpha = 150))
    hist(rnorm(n_int , 7.27, 2), add=T, col=rgb(238, 64, 0, max = 255, alpha = 150))
    
    sigma_R_max <- df_1$sigma_R_max
    hist(sigma_R_max, col=rgb(95,158,160, max = 255, alpha = 150))
    hist(rnorm(n_int , 0.64, 0.25), add=T, col=rgb(238, 64, 0, max = 255, alpha = 150))
    
    mu_sigma_R <- df_1$mu_sigma_R
    hist(mu_sigma_R, col=rgb(95,158,160, max = 255, alpha = 150))
    hist(rnorm(n_int , 0, 5), add=T, col=rgb(238, 64, 0, max = 255, alpha = 150))
    
    sigma_sigma_R <- df_1$sigma_sigma_R
    hist(sigma_sigma_R, col=rgb(95,158,160, max = 255, alpha = 150))
    hist(rnorm(n_int , 0, 5), add=T, col=rgb(238, 64, 0, max = 255, alpha = 150))
    
    #Spatial hyperparameters
    alpha_dist <- df_1$alpha_dist_MS
    hist(alpha_dist, col=rgb(95,158,160, max = 255, alpha = 150))
    hist(rnorm(n_int , 0, 5), add=T, col=rgb(238, 64, 0, max = 255, alpha = 150))
    
    rho_dist <- df_1$rho_dist_MS
    hist(rho_dist, col=rgb(95,158,160, max = 255, alpha = 150), breaks=100)
    hist(rgamma(n_int , 1,0.1), col=rgb(238, 64, 0, max = 255, alpha = 150), add=T)
    
    sigma_dist <- df_1$sigma_dist_MS
    hist(sigma_dist, col=rgb(95,158,160, max = 255, alpha = 150))
    hist(rnorm(n_int , 0,5), col=rgb(238, 64, 0, max = 255, alpha = 150), add=T)
    
    #Hatchery offset terms
    mu_offset <- df_1$mu_offset
    hist(mu_offset, col=rgb(95,158,160, max = 255, alpha = 150))
    hist(rnorm(n_int , 0, 5), add=T, col=rgb(238, 64, 0, max = 255, alpha = 150))
    
    sigma_offset <- df_1$sigma_offset
    hist(sigma_offset, col=rgb(95,158,160, max = 255, alpha = 150))
    hist(rnorm(n_int , 0, 5), add=T, col=rgb(238, 64, 0, max = 255, alpha = 150))
    
    #concentration term for beta-binomial
    k_bb <- df_1$k
    hist(k_bb, col=rgb(95,158,160, max = 255, alpha = 150))
    hist(runif(n_int, 2, 500), add=T, col=rgb(238, 64, 0, max = 255, alpha = 150))
    
    #Initial marine survival terms
    mu_init_surv <- df_1$mu_init_surv
    hist(mu_init_surv, col=rgb(95,158,160, max = 255, alpha = 150))
    hist(rnorm(n_int , 0, 5), add=T, col=rgb(238, 64, 0, max = 255, alpha = 150))
    
    sigma_init_surv <- df_1$sigma_init_surv
    hist(sigma_init_surv, col=rgb(95,158,160, max = 255, alpha = 150))
    hist(rnorm(n_int , 0, 5), add=T, col=rgb(238, 64, 0, max = 255, alpha = 150))
    
    #Mean of the marine survival time series
    mu_mu_surv <- df_1$mu_mu_surv
    hist(mu_mu_surv, col=rgb(95,158,160, max = 255, alpha = 150))
    hist(rnorm(n_int , 0, 5), add=T, col=rgb(238, 64, 0, max = 255, alpha = 150))
    
    sigma_mu_surv <- df_1$sigma_mu_surv
    hist(sigma_mu_surv, col=rgb(95,158,160, max = 255, alpha = 150))
    hist(rnorm(n_int , 0, 5), add=T, col=rgb(238, 64, 0, max = 255, alpha = 150))
    
    #Observation error terms
    sigma_catch <- df_1$sigma_catch 
    hist(sigma_catch , col=rgb(95,158,160, max = 255, alpha = 150))
    hist(rnorm(n_int , 0, 5), add=T, col=rgb(238, 64, 0, max = 255, alpha = 150))
    
    sigma_smolt <- df_1$sigma_smolt 
    hist(sigma_smolt , col=rgb(95,158,160, max = 255, alpha = 150))
    hist(rnorm(n_int , 0, 5), add=T, col=rgb(238, 64, 0, max = 255, alpha = 150))
    
    #Autocorrelation for marine survival time-series
    phi_MS <- df_1[,grep('phi',colnames(df_1))]
    for(i in 1:ncol(phi_MS)){
      hist(phi_MS[,i],col=rgb(95,158,160, max = 255, alpha = 150))
      hist(runif(n_int, -1, 1), add=T, col=rgb(238, 64, 0, max = 255, alpha = 150))
    }
    
    #Autocorrelation for marine survival time-series
    log_adult_init <- df_1[,grep('log_adult_init',colnames(df_1))]
    for(i in 1:ncol(log_adult_init)){
      hist(log_adult_init[,i],col=rgb(95,158,160, max = 255, alpha = 150))
      hist(rnorm(n_int, 0, 10), add=T, col=rgb(238, 64, 0, max = 255, alpha = 150))
    }
    
    dev.off()
  }
  return(mod_fit)
}

mod_fit <- exec_fun(coho_dat_full=dat_list$coho_dat_full, coord=dat_list$coord, exec=FALSE, n_iter=10000, n_thin=1, n_adapt=0.99, n_tree=10.25, n_chain=5, plots=FALSE,
                    sigma_esc = 0.2, 
                    mod_filename='LD_coho_forecast_6_2_4c_fit', cast_pop=c(cast_pop))

fit1.shiny <- as.shinystan(mod_fit)  
launch_shinystan(fit1.shiny)

#Look at posteriors for fitted model object#####
print(mod_fit, pars=c('sigma_R','mu_sigma_R', 'sigma_sigma_R'))
print(mod_fit, pars=c('smolt_survival'))
print(mod_fit, pars=c('log_alpha', 'log_R_max', 'log_adult_init'))
print(mod_fit, pars=c('sigma_smolt', 'sigma_catch'))
print(mod_fit, pars=c('sigma_u_init'))

print(mod_fit, pars=c('u_mat'))
print(mod_fit, pars=c('r_dev'))
print(mod_fit, pars=c('adult_est'))
print(mod_fit, pars=c('smolt_est'))
print(mod_fit, pars=c('harvest_est'))
print(mod_fit, pars=c('sigma_alpha', 'sigma_R_max', 'mu_alpha', 'mu_R_max'))
print(mod_fit, pars=c('Omega'))
print(mod_fit, pars=c('alpha_dist_MS', 'rho_dist_MS', 'sigma_dist_MS'))
print(mod_fit, pars=c('sigma_h', 'rho_h'))

print(mod_fit, pars=c('k'))
print(mod_fit, pars=c('init_surv'))
print(mod_fit, pars=c('u_mat_init'))
print(mod_fit, pars=c('mu_mu_surv'))
print(mod_fit, pars=c('sigma_mu_surv'))
print(mod_fit, pars=c('hatch_offset', 'mu_offset', 'sigma_offset'))
print(mod_fit, pars=c('mu_log_adult_init', 'sigma_log_adult_init', 'z_log_adult_init'))
print(mod_fit, pars=c('log_adult_init'))
print(mod_fit, pars=c('adult_pred'))

df_1 <- as.data.frame(mod_fit)
median(exp(df_1$mu_R_max))
median(exp(df_1$mu_alpha))
median(exp(df_1$mu_sigma_R))

quantile(exp(df_1$mu_sigma_R), c(0.025, 0.975))
quantile(exp(df_1$mu_R_max), c(0.025, 0.975))
quantile(exp(df_1$mu_alpha), c(0.025, 0.975))

unique(dat_list$coho_dat_full$Population)[order(colMedian(exp(df_1[,grep('log_alpha', colnames(df_1))])), decreasing=TRUE)]
unique(dat_list$coho_dat_full$Population)[order(colMedian(exp(df_1[,grep('log_R_max', colnames(df_1))])), decreasing=TRUE)]
unique(dat_list$coho_dat_full$Population)[order(colMedian(exp(df_1[,grep('log_sigma_R', colnames(df_1))])), decreasing=TRUE)]
length(which(colMedian(df_1[,grep('phi', colnames(df_1))]) < 0.25))
length(which(colMedian(df_1[,grep('phi', colnames(df_1))]) > 0.5))
#####

# #Function for executing simple univariate state space time-series models for forecasting -->  RW = random walk, AR= lag-1 AR, MA = lag-1 moving average#####
# exec_fun_RW <- function(coho_dat_full=dat_list$coho_dat_full, exec=TRUE, n_iter=1000, n_thin=1, n_adapt=0.8, n_tree=10, n_chain=4,
#                         sigma_esc = 0.2, form='MA', mod_filename='LD_simple_MA_ind_fit'){
#   #Grab file and save names based on the model being fitted
#   if(form=='RW'){
#     mod = 'LD_coho_forecast_RW.stan'
#   }
#   if(form=='AR'){
#     mod = 'LD_coho_forecast_AR.stan'
#   }  
#   if(form=='MA'){
#     mod = 'LD_coho_forecast_MA_ind.stan'
#   }
#   
#   #Compute number of populations
#   n_pop <- length(unique(coho_dat_full$Population))
#   
#   #Compute populations with escapement data
#   pop_catch <- unique(coho_dat_full$pop[which(is.na(coho_dat_full$Harvest)==FALSE)])
#   pop_esc <- unique(coho_dat_full$pop[which(is.na(coho_dat_full$Spawners)==FALSE)])
#   
#   #Compute number of populations with each data type available
#   n_pop_catch <- length(pop_catch)
#   n_pop_esc <- length(pop_esc)
#   
#   #create vectors to store the number of years of smolt and escapement data available for each population being considered
#   n_year_true_catch <- vector(length=n_pop_catch)
#   n_year_true_esc <- vector(length=n_pop_esc)
#   
#   #Calculate the number of years for which each data type is available for each population
#   #escapement data
#   for (i in 1:n_pop_catch){
#     n_year_true_catch[i] <- length(na.omit(coho_dat_full$Harvest[coho_dat_full$pop==pop_catch[i]]))
#   }
#   for (i in 1:n_pop_esc){
#     n_year_true_esc[i] <- length(na.omit(coho_dat_full$Spawners[coho_dat_full$pop==pop_esc[i]]))
#   }
#   
#   #Compile all the years across populations for which there is escapement and smolt data (not NAs)
#   catch_true <- coho_dat_full$yr[which(is.na(coho_dat_full$Harvest)==FALSE)]
#   esc_true <- coho_dat_full$yr[which(is.na(coho_dat_full$Spawners)==FALSE)]
#   
#   #Compile all the actual observations of smolt and escapement data that are not NA
#   catch_dat <- coho_dat_full$Harvest[which(is.na(coho_dat_full$Harvest)==FALSE)]  
#   esc_dat <- coho_dat_full$Spawners[which(is.na(coho_dat_full$Spawners)==FALSE)]  
#   
#   #Create slice points for cutting up the intact (non-NA) smolt and escapement data by population (stan does not accept ragged data structures
#   slice_catch_start <- vector(length=n_pop_catch)
#   slice_catch_start[1] <- 1
#   slice_catch_end <- vector(length=n_pop_catch)
#   slice_catch_end[n_pop_catch] <- length(catch_true)
#   
#   slice_esc_start <- vector(length=n_pop_esc)
#   slice_esc_start[1] <- 1
#   slice_esc_end <- vector(length=n_pop_esc)
#   slice_esc_end[n_pop_esc] <- length(esc_true)
#   
#   for(p in 2:n_pop_catch){
#     slice_catch_start[p] <- slice_catch_start[p-1] + length(na.omit(coho_dat_full$Harvest[coho_dat_full$pop==pop_catch[p-1]]))
#   }
#   slice_catch_end[1:(n_pop_catch-1)] <- slice_catch_start[2:n_pop_catch]-1
#   
#   for(p in 2:n_pop_esc){
#     slice_esc_start[p] <- slice_esc_start[p-1] + length(na.omit(coho_dat_full$Spawners[coho_dat_full$pop==pop_esc[p-1]]))
#   }
#   slice_esc_end[1:(n_pop_esc-1)] <- slice_esc_start[2:n_pop_esc]-1
#   
#   #test if the cut points are correctly pulling the right populations and plot
#   pdf(file='data_RW.pdf')
#   par(mfrow=c(3,3))
#   for(i in 1:n_pop_catch){
#     print(all(catch_dat[slice_catch_start[i]:slice_catch_end[i]] == na.omit(coho_dat_full$Harvest[coho_dat_full$pop==pop_catch[i]])))
#     plot(coho_dat_full$Harvest[coho_dat_full$pop==pop_catch[i]], ylab='Harvest', xlab='Year', main=paste(unique(coho_dat_full$Name)[pop_catch[i]]), type='o', pch=16, col='black')
#   }
#   for(i in 1:n_pop_esc){
#     print(all(esc_dat[slice_esc_start[i]:slice_esc_end[i]] == na.omit(coho_dat_full$Spawners[coho_dat_full$pop==pop_esc[i]])))
#     plot(coho_dat_full$Spawners[coho_dat_full$pop==pop_esc[i]], ylab='Escapement', xlab='Year', main=paste(unique(coho_dat_full$Name)[pop_esc[i]]), type='o', pch=16, col='black')
#   }
#   dev.off()
# 
#   #Fit stan model
#   if(exec==TRUE){
#     mod_fit_RW <- stan(file = mod, data = list(
#       n_year <- length(unique(coho_dat_full$Calendar.Year))+1,
#       n_pop <- length(unique(coho_dat_full$pop)),
#       pop_esc <- pop_esc,
#       pop_catch <- pop_catch,
#       
#       n_pop_esc <- n_pop_esc,
#       n_pop_catch <- n_pop_catch,
#       
#       esc_true <- esc_true,
#       catch_true <- catch_true,
#       
#       esc_dat <- esc_dat,
#       catch_dat <-catch_dat,
#       
#       n_esc <- length(esc_dat),
#       n_catch <- length(catch_dat),
#       
#       sigma_esc <- sigma_esc,
#       #sigma_catch <- sigma_catch,
#       
#       slice_esc_start <- slice_esc_start,
#       slice_esc_end <-  slice_esc_end,
#     
#       slice_catch_start <- slice_catch_start,
#       slice_catch_end <-  slice_catch_end
#       
#       ),iter = n_iter, chains = n_chain, control=list(adapt_delta=n_adapt, max_treedepth=n_tree), thin=n_thin, seed=666)
#     
#     saveRDS(mod_fit_RW, paste(mod_filename,'.rds',sep=''))
#   }
#   
#   if(exec==FALSE){
#     mod_fit_RW <- readRDS(paste(mod_filename,'.rds',sep=''))
#   }
#   
#   df_1_RW <- as.data.frame(mod_fit_RW)
#   
#   return(mod_fit_RW)
# }
# x <- exec_fun_RW()
# print(x, pars=c('sigma_catch'))
# print(x, pars=c('sigma_proc'))
#####

#Function for executing simple univariate state space time-series models for forecasting -->  RW = random walk, AR= lag-1 AR, MA = lag-1 moving average
exec_fun_RW <- function(coho_dat_full=dat_list$coho_dat_full, exec=TRUE, n_iter=1000, n_thin=1, n_adapt=0.8, n_tree=10, n_chain=4,
                        form='RW', mod_filename='LD_simple_RW_ind_2_fit'){
  #Grab file and save names based on the model being fitted
  if(form=='RW'){
    mod = 'LD_coho_forecast_RW_ind_2.stan'
  }
  if(form=='AR'){
    mod = 'LD_coho_forecast_AR_ind_2.stan'
  }  
  if(form=='MA'){
    mod = 'LD_coho_forecast_MA_ind_2.stan'
  }
  
  #Take sum of total adult returns
  coho_dat_full$Total_return <- coho_dat_full$Spawners + coho_dat_full$Harvest
  
  #Compute number of populations
  n_pop <- length(unique(coho_dat_full$Population))
  
  #Compute populations with escapement data
  pop_tot <- unique(coho_dat_full$pop[which(is.na(coho_dat_full$Total_return)==FALSE)])

  #Compute number of populations with each data type available
  n_pop_tot <- length(pop_tot)

  #create vectors to store the number of years of smolt and escapement data available for each population being considered
  n_year_true_tot <- vector(length=n_pop_tot)

  #Calculate the number of years for which each data type is available for each population
  #escapement data
  for (i in 1:n_pop_tot){
    n_year_true_tot[i] <- length(na.omit(coho_dat_full$Total_return[coho_dat_full$pop==pop_tot[i]]))
  }

  #Compile all the years across populations for which there is escapement and smolt data (not NAs)
  tot_true <- coho_dat_full$yr[which(is.na(coho_dat_full$Total_return)==FALSE)]

  #Compile all the actual observations of smolt and escapement data that are not NA
  tot_dat <- coho_dat_full $Total_return[which(is.na(coho_dat_full$Total_return)==FALSE)]  

  #Create slice points for cutting up the intact (non-NA) smolt and escapement data by population (stan does not accept ragged data structures
  slice_tot_start <- vector(length=n_pop_tot)
  slice_tot_start[1] <- 1
  slice_tot_end <- vector(length=n_pop_tot)
  slice_tot_end[n_pop_tot] <- length(tot_true)
  
  for(p in 2:n_pop_tot){
    slice_tot_start[p] <- slice_tot_start[p-1] + length(na.omit(coho_dat_full$Total_return[coho_dat_full$pop==pop_tot[p-1]]))
  }
  slice_tot_end[1:(n_pop_tot-1)] <- slice_tot_start[2:n_pop_tot]-1
  
  #test if the cut points are correctly pulling the right populations and plot
  pdf(file='data_RW.pdf')
  par(mfrow=c(3,3))
  for(i in 1:n_pop_tot){
    print(all(tot_dat[slice_tot_start[i]:slice_tot_end[i]] == na.omit(coho_dat_full$Total_return[coho_dat_full$pop==pop_tot[i]])))
    plot(coho_dat_full$Total_return[coho_dat_full$pop==pop_tot[i]], ylab='Total_return', xlab='Year', main=paste(unique(coho_dat_full$Name)[pop_tot[i]]), type='o', pch=16, col='black')
  }
  dev.off()
  #####
  #Fit stan model####
  if(exec==TRUE){
    mod_fit_RW <- stan(file = mod, data = list(
      n_year <- length(unique(coho_dat_full$Calendar.Year))+1,
      n_pop <- length(unique(coho_dat_full$pop)),
      pop_tot <- pop_tot,
      
      n_pop_tot <- n_pop_tot,
      
      tot_true <- tot_true,
      
      tot_dat <-tot_dat,
      
      n_tot <- length(tot_dat),
      
      slice_tot_start <- slice_tot_start,
      slice_tot_end <-  slice_tot_end
      
    ),iter = n_iter, chains = n_chain, control=list(adapt_delta=n_adapt, max_treedepth=n_tree), thin=n_thin, seed=666)
    
    saveRDS(mod_fit_RW, paste(mod_filename,'.rds',sep=''))
  }
  
  if(exec==FALSE){
    mod_fit_RW <- readRDS(paste(mod_filename,'.rds',sep=''))
  }
  
  df_1_RW <- as.data.frame(mod_fit_RW)
  
  return(mod_fit_RW)
}
x <- exec_fun_RW()
print(x, pars=c('sigma_obs'))
print(x, pars=c('sigma_proc'))


#Simulate data from known parameters
sim_func <- function(
  coho_dat_full = dat_list$coho_dat_full,
  coord = dat_list$coord,
  sigma_esc = 0.2,
  df_1 = as.data.frame(mod_fit)){
  
  #set the seed
  set.seed(666)
  
  #Compute number of populations
  n_pop <- length(unique(coho_dat_full$Population))
  n_year <- length(unique(coho_dat_full$Calendar.Year))
  
  #Indices for hatchery versus wild populations
  hatchery <- as.vector(na.omit(unique(coho_dat_full$pop[coho_dat_full$Hatchery==1])))
  wild <- (1:n_pop)[-hatchery]
  
  n_hatchery <- length(hatchery)
  n_wild <- length(wild)
  
  #Vectors of releases from CWT data (used to establish reasoable ranges of values for simulating binomial data)
  MS_dat_N <- coho_dat_full$Release_No [which(is.na(coho_dat_full$Release_No)==FALSE)]   
  
  #Create distance matrix
  distance_matrix <- as.matrix(dist(coord, 'euclidean', diag=TRUE, upper=TRUE)/10000)
  
  #Define parameters
  #Marine survival GP
  rho_dist_sim <- median(df_1$rho_dist_MS)
  alpha_dist_sim <- median(df_1$alpha_dist_MS)
  sigma_dist_sim <- median(df_1$sigma_dist_MS)
  
  #Marine survival means and autocorrelation
  mu_mu_surv_sim <- median(df_1$mu_mu_surv)
  sigma_mu_surv_sim <- median(df_1$sigma_mu_surv)
  
  phi_sim <- colMedian(df_1[,grep('phi', colnames(df_1))])
  
  #Stock-recruit hyper parameters
  mu_log_alpha_sim <- median(df_1$`gamma[1,1]`)
  mu_log_R_max_sim <- median(df_1$`gamma[2,1]`)
  sd_log_alpha_sim <- median(df_1$`tau[1]`)
  sd_log_R_max_sim <- median(df_1$`tau[2]`)
  mu_sigma_R_sim <- median(df_1$mu_sigma_R)
  sigma_sigma_R_sim <- median(df_1$sigma_sigma_R)
  
  #Hatchery offset terms
  mu_offset_sim <- median(df_1$mu_offset)
  sigma_offset_sim <- median(df_1$sigma_offset)
  Omega_cor <- median(df_1$`Omega[2,1]`)
  
  #Initial marine survival terms
  #init_surv_sim <- colMedian(df_1[,grep('init_surv', colnames(df_1))])[47:90]
  
  init_surv_sim <- colMedian(df_1[,grep('init_surv', colnames(df_1))])[39:74]
  
  #Average marine survival
  mu_surv_sim <- rnorm(n_pop, mu_mu_surv_sim, sigma_mu_surv_sim)
  
  #Create marine survival variance covariance matrix
  cov_sim <- matrix(nrow=n_pop, ncol=n_pop)
  for(i in 1:n_pop){
    for(j in 1:n_pop){
      if(i==j){
        cov_sim[i,j] <- (alpha_dist_sim^2)*exp(-1*((distance_matrix[i,j]^2)/(2*rho_dist_sim^2))) + sigma_dist_sim^2
      }else{
        cov_sim[i,j] <- (alpha_dist_sim^2)*exp(-1*((distance_matrix[i,j]^2)/(2*rho_dist_sim^2)))
      }
    }
  }
  
  MS_dev_pop_sim <- matrix(nrow=(n_year-1), ncol=n_pop)
  for(y in 1:(n_year-1)){
    MS_dev_pop_sim[y,] <- rmvnorm(1, mean=rep(0, n_pop), sigma=cov_sim)
  }
  
  #PROCESS MODEL
  #1) Population-specific marine survival deviations
  #Marine surival deviations for the coastal populations
  #Population-specific marine survival deviations from the global time-series specified as a Gaussian process
  logit_smolt_survival_pop_sim <- matrix(nrow=n_year, ncol=n_pop)
  
  for(i in 1:n_pop){
    for(y in 1:n_year){
      if(y==1){
        logit_smolt_survival_pop_sim[y,i] = init_surv_sim[i];
      }else{
        logit_smolt_survival_pop_sim[y,i] = mu_surv_sim[i] + phi_sim[i]*(logit_smolt_survival_pop_sim[y-1,i]-mu_surv_sim[i]) + MS_dev_pop_sim[y-1, i];
      }
    }
  }
  
  #Inverse logit transform marine survival deviations #pick up here
  smolt_survival_sim <- inv.logit(logit_smolt_survival_pop_sim)
  smolt_survival_adj_sim <- matrix(nrow=n_year, ncol=n_pop)
  
  #Add hatchery offset terms, create the marine survival estimates
  hatch_offset_sim <- rtnorm(n_hatchery, mu_offset_sim, sigma_offset_sim, upper=0)
  for (y in 1:n_year){
    for(j in 1:n_hatchery){
      smolt_survival_adj_sim[y,hatchery[j]] <-  inv.logit(logit_smolt_survival_pop_sim[y,hatchery[j]] + hatch_offset_sim[j])
    }
  }
  smolt_survival_adj_sim[,wild] <- smolt_survival_sim[,wild]
  
  #3) Stock and recruitment
  #reconstruct covariance matrix from marginal sigmas and correlation matrix
  Omega_cov_sim <- matrix(nrow=2, ncol=2)
  Omega_cov_sim[1,1] <- sd_log_alpha_sim^2
  Omega_cov_sim[2,1] <- Omega_cor*sd_log_alpha_sim*sd_log_R_max_sim
  Omega_cov_sim[1,2] <- Omega_cor*sd_log_alpha_sim*sd_log_R_max_sim
  Omega_cov_sim[2,2] <- sd_log_R_max_sim^2
  
  #Double check using preexisiting software
  #Omega <- matrix(nrow=2, ncol=2)
  #Omega[1,1] <- median(df_1$`Omega[1,1]`)
  #Omega[2,1] <- median(df_1$`Omega[2,1]`)
  #Omega[1,2] <- median(df_1$`Omega[1,2]`)
  #Omega[2,2] <- median(df_1$`Omega[2,2]`)
  #devtools::install_github("easystats/correlation")
  #library("correlation")
  #cor_to_cov(cor=Omega, sd=c(median(df_1$`tau[1]`),median(df_1$`tau[2]`)))
  
  #Generate S-R populations using multivariate normal distribution
  SR_sim <- rmvnorm(n_pop, mean= c(mu_log_alpha_sim, mu_log_R_max_sim), sigma=Omega_cov_sim)
  sigma_R_sim <- exp(rnorm(n_pop, mu_sigma_R_sim, sigma_sigma_R_sim)) 
  alpha_sim <- exp(SR_sim[,1])
  R_max_pre_sim <- exp(SR_sim[,2])
  R_max_sim <- exp(SR_sim[,2])*aggregate(coho_dat_full$KM, by=list('Population'=coho_dat_full$Population), mean)$x
  
  #Generate unobserved spawning populations
  adult_init_sim <- matrix(nrow=2, ncol=n_pop)
  adult_init_sim[1,]  <- colMedian(exp(df_1[,grep('log_adult_init', colnames(df_1))]))[seq(1, n_pop*2, 2)]
  adult_init_sim[2,]  <- colMedian(exp(df_1[,grep('log_adult_init', colnames(df_1))]))[seq(2, n_pop*2, 2)]
  
  #Generate recruitment deviations
  r_dev_sim <- matrix(nrow=n_year, ncol=n_pop)
  for (i in 1:n_pop){
    #r_dev_sim[,i] <- rnorm(n_year, -(sigma_R_sim[i]^2)/2 ,sigma_R_sim[i]); # Recruitment deviations lognormally distributed with bias correction factor
    r_dev_sim[,i] <- rnorm(n_year, 0 ,sigma_R_sim[i]); # Recruitment deviations lognormally distributed with bias correction factor
    
  }
  
  #calculate exploitation rates
  u_mat_init <- matrix(colMedian(df_1[,grep('u_mat_init', colnames(df_1))]), nrow=1, ncol=n_pop, byrow=FALSE)/5
  u_mat_sim <- matrix(ncol=n_pop, nrow=n_year)
  u_mat_sim[1,] <- u_mat_init
  sigma_h_sim <- median(df_1$sigma_h)
  rho_h_sim <- median(df_1$rho_h)
  
  K_h <- matrix(nrow=n_pop, ncol=n_pop)
  for (i in 1:n_pop){
    for(j in 1:n_pop){
      if (i == j){
        K_h[i,j] = sigma_h_sim*sigma_h_sim; 
      } else { 
        K_h[i,j] = rho_h_sim*(sigma_h_sim*sigma_h_sim);
      }
    }
  }
  
  h_dev_sim<- rmvnorm(n_year-1, mean=rep(0, n_pop), sigma=K_h)
  
  
  for(i in 1:n_pop){
    for(y in 2:n_year){
      u_mat_sim[y,i] <- u_mat_sim[y-1,i] + h_dev_sim[y-1,i]
    }
  }
  
  u_mat_sim <- inv.logit(u_mat_sim)
  
  #Calculate smolt recruitment, escapement, and harvest
  harvest_est_sim <-matrix(nrow=n_year, ncol=n_pop)
  smolt_est_sim <- matrix(nrow=n_year, ncol=n_pop)
  adult_est_sim <- matrix(nrow=n_year, ncol=n_pop)
  for(i in 1:n_pop){
    for(y in 1:n_year){
      if (y< 3){ # First two years of the time series where recruitments arise from unobserved spawning events
        harvest_est_sim[y,i] <- adult_init_sim[y, i]*u_mat_sim[y,i]; # harvest calculated as total return times exploitation rate
        
        adult_est_sim[y,i] <- adult_init_sim[y, i]- harvest_est_sim[y,i]; # Escapement calculated as total return minus the harvest
        
        smolt_est_sim[y,i] <-  (adult_est_sim[y, i]/(1/alpha_sim[i]  + adult_est_sim[y, i]/R_max_sim[i]))*exp(r_dev_sim[y,i]); # Calculate the resulting smolt recruitment from the escapement according to the Beverton-Holt function
        
      }else{ #Later years in the time-series where all recruitments arise from observed spawning events
        
        smolt_est_sim[y,i] <- (adult_est_sim[y-2, i]/(1/alpha_sim[i]  + adult_est_sim[y-2, i]/R_max_sim[i]))*exp(r_dev_sim[y,i]); #Smolt recruitment calculated as Beverton-Holt function of spawning abundance 
        
        harvest_est_sim[y,i] <- (smolt_est_sim[y-1,i]*smolt_survival_sim[y-1,i])*u_mat_sim[y,i]; # Harvest is the smolt recruitment, times the marine survival rate of smolts, times the annual explopitation rates
        
        adult_est_sim[y,i] <- smolt_est_sim[y-1,i]*smolt_survival_sim[y-1,i] - harvest_est_sim[y,i]; #Escapement calculated as the smolt recruitment times marine survival minue the harvest
      }
    }
  }
  
  #Generate data with errors from the the true biological states
  smolt_dat_sim <- matrix(nrow=n_year, ncol=n_pop)
  esc_dat_sim <- matrix(nrow=n_year, ncol=n_pop)
  catch_dat_sim <- matrix(nrow=n_year, ncol=n_pop)
  MS_dat_x_sim <- matrix(nrow=n_year, ncol=n_pop)
  MS_dat_N_sim <- matrix(nrow=n_year, ncol=n_pop)
  p_MS_dat_sim <- matrix(nrow=n_year, ncol=n_pop)
  K_sim <- median(df_1$k)
  sigma_smolt_sim <- median(df_1$sigma_smolt)
  sigma_catch_sim <- median(df_1$sigma_catch)
  
  for(i in 1:n_year){
    for(j in 1:n_pop){
      smolt_dat_sim[i,j] <- rlnorm(1, log(smolt_est_sim[i,j]), sigma_smolt_sim)
      esc_dat_sim[i,j] <- rlnorm(1, log(adult_est_sim[i,j]), sigma_esc)
      catch_dat_sim[i,j] <- rlnorm(1, log(harvest_est_sim[i,j]), sigma_catch_sim)
      p_MS_dat_sim[i,j] <- rbeta(1, shape1=smolt_survival_adj_sim[i,j]*(K_sim-2)+1 , shape2=(1-smolt_survival_adj_sim[i,j])*(K_sim-2)+1)
      MS_dat_N_sim[i,j] <- round(rtnorm(1, mean(MS_dat_N), sd=sd(MS_dat_N), lower=1))
      MS_dat_x_sim[i,j] <- rbinom(1, size=MS_dat_N_sim[i,j], prob=p_MS_dat_sim[i,j])
    }
  }
  
  #format data for reading into stan
  coho_dat_sim <- cbind(melt(smolt_dat_sim),melt(esc_dat_sim),melt(catch_dat_sim),melt(MS_dat_x_sim), melt(MS_dat_N_sim))
  coho_dat_sim <- coho_dat_sim[,c(1,2,3,6,9,12,15)]
  colnames(coho_dat_sim) <- c('yr', 'pop', 'Smolt.Abundance', 'Spawners', 'Harvest', 'Fishery_Plus_Escapement', 'Release_No')
  
  #Create output storage list
  dat_sim_list <- list()
  
  #Store the actual data to which the model will be fitted
  dat_sim_list$coho_dat_sim <- coho_dat_sim
  
  #Store the hyperparameters from which the data was generated for later comparison
  dat_sim_list$rho_dist_sim <-rho_dist_sim
  dat_sim_list$alpha_dist_sim <- alpha_dist_sim
  dat_sim_list$sigma_dist_sim <-sigma_dist_sim
  dat_sim_list$init_surv_sim <-init_surv_sim
  dat_sim_list$mu_log_alpha_sim <- mu_log_alpha_sim
  dat_sim_list$mu_log_R_max_sim <- mu_log_R_max_sim
  dat_sim_list$sd_log_alpha_sim <- sd_log_alpha_sim
  dat_sim_list$sd_log_R_max_sim <- sd_log_R_max_sim
  dat_sim_list$mu_sigma_R_sim <- mu_sigma_R_sim
  dat_sim_list$sigma_sigma_R_sim <- sigma_sigma_R_sim
  dat_sim_list$mu_offset_sim <- mu_offset_sim
  dat_sim_list$sigma_offset_sim <- sigma_offset_sim
  dat_sim_list$Omega_cor <- Omega_cor
  dat_sim_list$k_sim <- K_sim
  dat_sim_list$sigma_smolt_sim <- sigma_smolt_sim
  dat_sim_list$sigma_catch_sim <- sigma_catch_sim
  dat_sim_list$sigma_h_sim <- sigma_h_sim
  dat_sim_list$rho_h_sim <- rho_h_sim
  
  #Population specific parameters
  dat_sim_list$sigma_R_sim <-  sigma_R_sim 
  dat_sim_list$alpha_sim  <-  alpha_sim
  dat_sim_list$R_max_sim <-  R_max_sim 
  dat_sim_list$R_max_pre_sim <- R_max_pre_sim
  dat_sim_list$adult_init_sim <- adult_init_sim
  dat_sim_list$r_dev_sim <- r_dev_sim
  dat_sim_list$hatch_offset_sim <- hatch_offset_sim
  
  dat_sim_list$harvest_est_sim <- harvest_est_sim
  dat_sim_list$adult_est_sim <- adult_est_sim
  dat_sim_list$smolt_est_sim <- smolt_est_sim
  
  #marine survival time series
  dat_sim_list$smolt_survival_sim <- smolt_survival_sim
  dat_sim_list$smolt_survival_adj_sim <- smolt_survival_adj_sim
  dat_sim_list$phi_sim <- phi_sim
  dat_sim_list$mu_mu_surv_sim <- mu_mu_surv_sim
  dat_sim_list$sigma_mu_surv_sim<- sigma_mu_surv_sim
  return(dat_sim_list)
}

#Execute simulation function
sim_list <- sim_func(
  coho_dat_full = dat_list$coho_dat_full,
  coord = dat_list$coord,
  sigma_esc = 0.2,
  df_1 = as.data.frame(mod_fit))

#Fit the model to simulated data
exec_func_sim <- function(coho_dat_sim=sim_list$coho_dat_sim, n_iter = 2000, n_chains = 5, n_adapt=0.8, n_tree=10, n_thin=1, 
                          coho_dat_full = dat_list$coho_dat_full,
                          coord = dat_list$coord,
                          sigma_esc=0.2,exec=FALSE, plots=TRUE){
  
  
  #define number of years and populations
  n_year <- length(unique(coho_dat_sim$yr))
  n_pop <- length(unique(coho_dat_sim$pop))
  
  #Define indices derived from the actual data
  #Indices for hatchery versus wild populations
  hatchery <- as.vector(na.omit(unique(coho_dat_full$pop[coho_dat_full$Hatchery==1])))
  wild <- (1:n_pop)[-hatchery]
  
  n_hatchery <- length(hatchery)
  n_wild <- length(wild)
  
  #Vectors of releases from CWT data (used to establish reasoable ranges of values for simulating binomial data)
  MS_dat_N <- coho_dat_full$Release_No [which(is.na(coho_dat_full$Release_No)==FALSE)]   
  
  #Create distance matrix
  distance_matrix <- as.matrix(dist(coord, 'euclidean', diag=TRUE, upper=TRUE)/10000)
  
  stream_dist <- aggregate(coho_dat_full$KM, by=list('Population'=coho_dat_full$Population), mean)$x
  
  #Compute populations with smolt data
  pop_smolt_sim <- unique(coho_dat_sim$pop[which(is.na(coho_dat_sim$Smolt.Abundance)==FALSE)])
  #Compute populations with marine survival data
  pop_MS_sim <- unique(coho_dat_sim$pop[which(is.na(coho_dat_sim$Release_No)==FALSE)])
  #Compute populations with harvest data
  pop_catch_sim <- unique(coho_dat_sim$pop[which(is.na(coho_dat_sim$Harvest)==FALSE)])
  #Compute populations with escapement data
  pop_esc_sim <- unique(coho_dat_sim$pop[which(is.na(coho_dat_sim$Spawners)==FALSE)])
  
  #Compute number of populations with each data type available
  n_pop_smolt_sim <- length(pop_smolt_sim)
  n_pop_MS_sim <- length(pop_MS_sim)
  n_pop_catch_sim <- length(pop_catch_sim)
  n_pop_esc_sim <- length(pop_esc_sim)
  
  #create vectors to store the number of years of smolt and escapement data available for each population being considered
  n_year_true_smolt_sim <- vector(length=n_pop_smolt_sim)
  n_year_true_esc_sim <- vector(length=n_pop_esc_sim)
  n_year_true_harvest_sim <- vector(length=n_pop_catch_sim)
  n_year_true_MS_sim <- vector(length=n_pop_MS_sim)
  
  #Calculate the number of years for which each data type is available for each population
  #Smolt data
  for (i in 1:n_pop_smolt_sim){
    n_year_true_smolt_sim[i] <- length(na.omit(coho_dat_sim$Smolt.Abundance[coho_dat_sim$pop==pop_smolt_sim[i]]))
  }
  #escapement data
  for (i in 1:n_pop_esc_sim){
    n_year_true_esc_sim[i] <- length(na.omit(coho_dat_sim$Spawners[coho_dat_sim$pop==pop_esc_sim[i]]))
  }
  #harvest data
  for (i in 1:n_pop_catch_sim){
    n_year_true_harvest_sim[i] <- length(na.omit(coho_dat_sim$Harvest[coho_dat_sim$pop==pop_catch_sim[i]]))
  }
  #marine survival data
  for (i in 1:n_pop_MS_sim){
    n_year_true_MS_sim[i] <- length(na.omit(coho_dat_sim$Fishery_Plus_Escapement[coho_dat_sim$pop==pop_MS_sim[i]]))
  }
  
  #Compile all the years across populations for which there is escapement and smolt data (not NAs)
  esc_true_sim <- coho_dat_sim$yr[which(is.na(coho_dat_sim$Spawners)==FALSE)]
  smolt_true_sim <- coho_dat_sim$yr[which(is.na(coho_dat_sim$Smolt.Abundance)==FALSE)]
  harvest_true_sim <- coho_dat_sim$yr[which(is.na(coho_dat_sim$Harvest )==FALSE)]
  MS_true_sim <- coho_dat_sim$yr[which(is.na(coho_dat_sim$Fishery_Plus_Escapement)==FALSE)]
  
  #Compile all the actual observations of smolt and escapement data that are not NA
  smolt_dat_sim <- coho_dat_sim$Smolt.Abundance[which(is.na(coho_dat_sim$Smolt.Abundance)==FALSE)]                 
  esc_dat_sim <- coho_dat_sim$Spawners[which(is.na(coho_dat_sim$Spawners)==FALSE)]  
  harvest_dat_sim <- coho_dat_sim$Harvest[which(is.na(coho_dat_sim$Harvest )==FALSE)]      
  MS_dat_x_sim <- coho_dat_sim$Fishery_Plus_Escapement[which(is.na(coho_dat_sim$Fishery_Plus_Escapement)==FALSE)]      
  MS_dat_N_sim <- coho_dat_sim$Release_No [which(is.na(coho_dat_sim$Release_No)==FALSE)]      
  
  #Create slice points for cutting up the intact (non-NA) smolt and escapement data by population (stan does not accept ragged data structures)
  slice_smolt_start_sim <- vector(length=n_pop_smolt_sim)
  slice_smolt_start_sim[1] <- 1
  slice_smolt_end_sim <- vector(length=n_pop_smolt_sim)
  slice_smolt_end_sim[n_pop_smolt_sim] <- length(smolt_true_sim)
  
  slice_esc_start_sim <- vector(length=n_pop_esc_sim)
  slice_esc_start_sim[1] <- 1
  slice_esc_end_sim <- vector(length=n_pop_esc_sim)
  slice_esc_end_sim[n_pop_esc_sim] <- length(esc_true_sim)
  
  slice_harvest_start_sim <- vector(length=n_pop_catch_sim)
  slice_harvest_start_sim[1] <- 1
  slice_harvest_end_sim <- vector(length=n_pop_catch_sim)
  slice_harvest_end_sim[n_pop_catch_sim] <- length(harvest_true_sim)
  
  slice_MS_start_sim <- vector(length=n_pop_MS_sim)
  slice_MS_start_sim[1] <- 1
  slice_MS_end_sim <- vector(length=n_pop_MS_sim)
  slice_MS_end_sim[n_pop_MS_sim] <- length(MS_true_sim)
  
  j <- 1
  for(i in 2:length(smolt_true_sim)){
    if(smolt_true_sim[i] < smolt_true_sim[i-1]){
      j <- j+1
      slice_smolt_start_sim[j] <- i
    }
  }
  slice_smolt_end_sim[1:(n_pop_smolt_sim-1)] <- slice_smolt_start_sim[2:n_pop_smolt_sim]-1
  
  j <- 1
  for(i in 2:length(esc_true_sim)){
    if(esc_true_sim[i] < esc_true_sim[i-1]){
      j <- j+1
      slice_esc_start_sim[j] <- i
    }
  }
  slice_esc_end_sim[1:(n_pop_esc_sim-1)] <- slice_esc_start_sim[2:n_pop_esc_sim]-1
  
  
  j <- 1
  for(i in 2:length(harvest_true_sim)){
    if(harvest_true_sim[i] < harvest_true_sim[i-1]){
      j <- j+1
      slice_harvest_start_sim[j] <- i
    }
  }
  slice_harvest_end_sim[1:(n_pop_catch_sim-1)] <- slice_harvest_start_sim[2:n_pop_catch_sim]-1
  
  j <- 1
  for(i in 2:length(MS_true_sim)){
    if(MS_true_sim[i] < MS_true_sim[i-1]){
      j <- j+1
      slice_MS_start_sim[j] <- i
    }
  }
  slice_MS_end_sim[1:(n_pop_MS_sim-1)] <- slice_MS_start_sim[2:n_pop_MS_sim]-1
  
  #Fit the model to a complete simulated data set
  if(exec==TRUE){
    mod_fit_sim <- stan(file = 'LD_coho_forecast_6_2_4.stan', data = list(
      n_year <- n_year,
      n_pop <- n_pop,
      pop_smolt <- pop_smolt_sim,
      pop_esc <- pop_esc_sim,
      pop_catch <- pop_catch_sim,
      pop_MS <- pop_MS_sim,
      n_pop_smolt <-n_pop_smolt_sim,
      n_pop_esc <- n_pop_esc_sim,
      n_pop_catch <- n_pop_catch_sim,
      n_pop_MS <- n_pop_MS_sim,
      stream_dist <- stream_dist,
      smolt_true <- smolt_true_sim,
      esc_true <- esc_true_sim,
      harvest_true <- harvest_true_sim,
      MS_true <- MS_true_sim,
      smolt_dat <- smolt_dat_sim,
      esc_dat <- esc_dat_sim,
      harvest_dat <- harvest_dat_sim,
      MS_dat_x <- round(MS_dat_x_sim),
      MS_dat_N <- round(MS_dat_N_sim),
      n_smolt <- length(smolt_dat_sim),
      n_esc <- length(esc_dat_sim),
      n_harvest <- length(harvest_dat_sim),
      n_MS <- length(MS_dat_x_sim),
      sigma_esc <- sigma_esc,
      n_hatchery <- length(hatchery),
      hatchery <- hatchery,
      wild <- c(1:n_pop)[-hatchery],
      slice_smolt_start <- slice_smolt_start_sim,
      slice_smolt_end <- slice_smolt_end_sim,
      slice_esc_start <- slice_esc_start_sim,
      slice_esc_end <-  slice_esc_end_sim,
      slice_harvest_start <- slice_harvest_start_sim,
      slice_harvest_end <-  slice_harvest_end_sim,
      slice_MS_start <- slice_MS_start_sim,
      slice_MS_end <-  slice_MS_end_sim,
      u=matrix(1, nrow=1, ncol=n_pop),
      dist=distance_matrix),
      iter =n_iter, chains = n_chains, control=list(adapt_delta=n_adapt, max_treedepth=n_tree), thin=n_thin, seed=666)
    
    saveRDS(mod_fit_sim, 'Coho_forecast_sim_AR.rds')
  }
  if(exec==FALSE){
    mod_fit_sim <- (readRDS('Coho_forecast_sim_AR.rds'))
  }
  
  yr <- min(coho_dat_full$Calendar.Year)-1
  n_year <- length(unique(coho_dat_sim$yr))
  n_pop <- length(levels(coho_dat_full$Population))
  
  if(plots==TRUE){
    df_1_sim <- as.data.frame(mod_fit_sim)
    #Plot the spatial kernal
    pdf(file='spatial_kernel_sim_AR.pdf')
    #fitted model parameters
    rho_dist_fit <- df_1_sim$rho_dist_MS
    alpha_dist_fit <- df_1_sim$alpha_dist_MS
    sigma_dist_fit <- df_1_sim$sigma_dist_MS
    #True parameters
    rho_dist_true <- sim_list$rho_dist_sim
    alpha_dist_true <- sim_list$alpha_dist_sim
    sigma_dist_true <- sim_list$sigma_dist_sim

    
    dist <- seq(0,25,0.1)
    cor_vec_med <- vector(length=length(dist))
    cor_vec_lower_50 <- vector(length=length(dist))
    cor_vec_upper_50 <- vector(length=length(dist))
    cor_vec_lower_95 <- vector(length=length(dist))
    cor_vec_upper_95 <- vector(length=length(dist))
    
    #Puget sound populations
    #Fitted curve
    for(i in 1:length(dist)){
      cor_vec_med[i] <- exp_func(d=dist[i], rho_dist=median(rho_dist_fit), alpha_dist=median(alpha_dist_fit))
      cor_vec_lower_50[i] <- exp_func(d=dist[i], rho_dist=cred.50(rho_dist_fit)[1], alpha_dist=cred.50(alpha_dist_fit)[1])
      cor_vec_upper_50[i] <- exp_func(d=dist[i], rho_dist=cred.50(rho_dist_fit)[2], alpha_dist=cred.50(alpha_dist_fit)[2])
      cor_vec_lower_95[i] <- exp_func(d=dist[i], rho_dist=cred(rho_dist_fit)[1], alpha_dist=cred(alpha_dist_fit)[1])
      cor_vec_upper_95[i] <- exp_func(d=dist[i], rho_dist=cred(rho_dist_fit)[2], alpha_dist=cred(alpha_dist_fit)[2])
    }
    par(mfrow=c(2,1), mar=c(1,0,0,0), oma=c(15,8,15,8))
    mat <- matrix(c(1,1,1,1,1,1,1,1,0,2,0,3,0,4,
                    5,5,5,5,5,5,5,5,0,6,0,7,0,8),2,14, byrow=TRUE)
    layout(mat=mat, widths=rep(c(rep(1, 8), 0.25, 0.5, 0.75, 0.5, 0.75, 0.5),2), heights=rep(1,28))
    
    plot(1:length(dist), cor_vec_med, col='red', type='l', ylim=c(0, 1), xaxt='n', xaxs='i', lwd=2)
    xx <- c(1:(length(dist)),(length(dist)):1)
    yy <- c(cor_vec_lower_50, rev(cor_vec_upper_50))
    polygon(x=xx, y=yy, density = -1, border = NA, col=rgb(238, 64, 0, max = 255, alpha = 150))
    yy <- c(cor_vec_lower_95, rev(cor_vec_upper_95))
    polygon(x=xx, y=yy, density = -1, border = NA, col=rgb(238, 64, 0, max = 255, alpha = 90))
    mtext(side=1, 'Distance (KM)', line=2.5)
    mtext(side=2, 'Marine survival correlation', line=2.5)
    axis(side=1, at=seq(0, length(dist)-1, 10), labels=seq(0, length(dist)-1, 10))
    #True curve
    true_cor_vec <- vector(length=length(dist))
    for(i in 1:length(dist)){
      true_cor_vec[i] <- exp_func(d=dist[i], rho_dist=rho_dist_true, alpha_dist=alpha_dist_true)
    }
    points(1:length(dist), true_cor_vec, col='black', lwd=2, type='l')
    
    #Individual parameters
    plot(1, median(alpha_dist_fit),yaxt='n', xaxt='n', ylim=c(0,quantile(df_1_sim$alpha_dist_MS, 0.98)), xlim=c(0.5,1.5), pch=16, col='red', cex=1.5)
    arrows((1), quantile(df_1_sim$alpha_dist_MS, 0.25), 1, quantile(df_1_sim$alpha_dist_MS, 0.75), angle=90, code=3, length=0.0, lwd=5, lty=1, col=rgb(238, 64, 0, max = 255, alpha = 150))
    arrows((1), quantile(df_1_sim$alpha_dist_MS, 0.025), 1, quantile(df_1_sim$alpha_dist_MS, 0.975), angle=90, code=3, length=0.01, lwd=1.5, lty=1, col=rgb(238, 64, 0, max = 255, alpha = 150))
    mtext(side=1, expression(gamma), line=1.0, at=1)
    points(alpha_dist_true, col='black', pch=16)
    axis(side=4)
    
    plot(1, median(sigma_dist_fit),yaxt='n', xaxt='n', pch=16, col='red', cex=1.5, xlim=c(0.5,1.5))
    arrows((1), quantile(df_1_sim$sigma_dist_MS, 0.25), 1, quantile(df_1_sim$sigma_dist_MS, 0.75), angle=90, code=3, length=0.0, lwd=5, lty=1, col=rgb(238, 64, 0, max = 255, alpha = 150))
    arrows((1), quantile(df_1_sim$sigma_dist_MS, 0.025), 1, quantile(df_1_sim$sigma_dist_MS, 0.975), angle=90, code=3, length=0.01, lwd=1.5, lty=1, col=rgb(238, 64, 0, max = 255, alpha = 150))
    mtext(side=1, expression(sigma[d]), line=1.0, at=1)
    points(sigma_dist_true, col='black', pch=16)
    axis(side=4)
    
    plot(1, median(rho_dist_fit),yaxt='n', xaxt='n', pch=16, col='red', cex=1.5, ylim=c(0,quantile(df_1_sim$rho_dist_MS, 0.98)), xlim=c(0.5,1.5))
    arrows((1), quantile(df_1_sim$rho_dist_MS, 0.25), 1, quantile(df_1_sim$rho_dist_MS, 0.75), angle=90, code=3, length=0.0, lwd=5, lty=1, col=rgb(238, 64, 0, max = 255, alpha = 150))
    arrows((1), quantile(df_1_sim$rho_dist_MS, 0.025), 1, quantile(df_1_sim$rho_dist_MS, 0.975), angle=90, code=3, length=0.01, lwd=1.5, lty=1, col=rgb(238, 64, 0, max = 255, alpha = 150))
    axis(side=4, at=seq(0,quantile(df_1_sim$rho_dist_MS, 0.98),0.5), labels=seq(0,quantile(df_1_sim$rho_dist_MS, 0.98),0.5)*10)
    mtext(side=1, expression(rho), line=1.0, at=1)
    mtext(side=4, 'Length Scale (KM)', line=2.5)
    points(rho_dist_true, col='black', pch=16)
    dev.off()
    
    #Plot posteriors of key model parameters
    pdf('SR_posteriors_sim_AR.pdf')
    par(mfrow=c(2,1), mar=c(1,0,0,0), oma=c(16,6,6,3))
    mat <- matrix(c(rep(1,n_pop),0,2,0,3,
                    rep(4, n_pop),0,5,0,6,
                    rep(7, n_pop),0,8,0,9),3,n_pop+4, byrow=TRUE)
    layout(mat=mat, widths=c(rep(0.5,n_pop), 0.0, 0.5, 0.25, 0.5, rep(0.5,n_pop), 0.0, 0.5,0.25,0.5), heights=c(1))
    alpha <- df_1_sim[,grep('log_alpha', colnames(df_1_sim))]
    R_max <- df_1_sim[,grep('log_R_max', colnames(df_1_sim))]
    sigma_R <- exp(df_1_sim[,grep('log_sigma_R', colnames(df_1_sim))])
    plot(1:n_pop, colMedian(alpha), col='skyblue4', ylim=c(min(cred.2(alpha)[1,]),max(cred.2(alpha)[2,])), xlim=c(1,n_pop), xaxt='n', yaxt='n', cex=1)
    axis(side=2, at=c(log(25), log(50), log(100), log(200), log(400), log(800), log(1600)), labels = c(25, 50, 100, 200, 400,800,1600), las=2)
    mtext(side=2, expression(alpha), line=4.2, cex=1)
    mtext(side=2, 'smolts/spawner', line=3.2, cex=0.85)
    arrows(1:n_pop, cred.2(alpha)[1,], 1:n_pop, cred.2(alpha)[2,], angle=90, code=3, length=0.01, lwd=1.5, lty=1, col='skyblue4')
    arrows(1:n_pop, cred.3(alpha)[1,], 1:n_pop, cred.3(alpha)[2,], angle=90, code=3, length=0.0, lwd=3, lty=1, col='skyblue4')
    points(1:n_pop, log(sim_list$alpha_sim), col='black', pch=16)
    plot((n_pop+1), median(df_1_sim$mu_alpha), col='skyblue4', ylim=c(min(cred.2(alpha)[1,]),max(cred.2(alpha)[2,])), xaxt='n', yaxt='n',cex=1)
    arrows((n_pop+1), quantile(df_1_sim$mu_alpha, 0.25), (n_pop+1), quantile(df_1_sim$mu_alpha, 0.75), angle=90, code=3, length=0.0, lwd=3, lty=1, col='skyblue4')
    arrows((n_pop+1), quantile(df_1_sim$mu_alpha, 0.025), (n_pop+1), quantile(df_1_sim$mu_alpha, 0.975), angle=90, code=3, length=0.01, lwd=1.5, lty=1, col='skyblue4')
    points(n_pop+1, sim_list$mu_log_alpha_sim, col='black', pch=16)
    plot(n_pop+1, median(df_1_sim$sigma_alpha), col='skyblue4', ylim=c(quantile(df_1_sim$sigma_alpha, 0.0), quantile(df_1_sim$sigma_alpha, 0.995)), xaxt='n', yaxt='n',cex=1)
    arrows((n_pop+1), quantile(df_1_sim$sigma_alpha, 0.25), (n_pop+1), quantile(df_1_sim$sigma_alpha, 0.75), angle=90, code=3, length=0.0, lwd=3, lty=1, col='skyblue4')
    arrows((n_pop+1), quantile(df_1_sim$sigma_alpha, 0.025), (n_pop+1), quantile(df_1_sim$sigma_alpha, 0.975), angle=90, code=3, length=0.01, lwd=1.5, lty=1, col='skyblue4')
    points(n_pop+1, sim_list$sd_log_alpha_sim, col='black', pch=16)
    axis(side=4)
    
    plot(1:n_pop, colMedian(R_max), col='goldenrod', ylim=c(min(cred.2(R_max)[1,]),max(cred.2(R_max)[2,])), xlim=c(1,n_pop), xaxt='n', yaxt='n', cex=1)
    axis(side=2, at=c(log(60), log(125), log(250), log(500), log(1000), log(2000), log(4000), log(8000)), labels = c(60, 125, 250, 500,1000,2000,4000, 8000), las=2)
    mtext(side=2, expression('R'[max]), line=4.2, cex=1)
    mtext(side=2, 'smolts/KM', line=3.2, cex=0.85)
    arrows(1:n_pop, cred.2(R_max)[1,], 1:n_pop, cred.2(R_max)[2,], angle=90, code=3, length=0.01, lwd=1.5, lty=1, col='goldenrod')
    arrows(1:n_pop, cred.3(R_max)[1,], 1:n_pop, cred.3(R_max)[2,], angle=90, code=3, length=0.0, lwd=3, lty=1, col='goldenrod')
    points(1:n_pop, log(sim_list$R_max_pre_sim), col='black', pch=16)
    plot((n_pop+1), median(df_1_sim$mu_R_max), col='goldenrod', ylim=c(min(cred.2(R_max)[1,]),max(cred.2(R_max)[2,])), xaxt='n', yaxt='n',cex=1)
    arrows((n_pop+1), quantile(df_1_sim$mu_R_max, 0.25), (n_pop+1), quantile(df_1_sim$mu_R_max, 0.75), angle=90, code=3, length=0.0, lwd=3, lty=1, col='goldenrod')
    arrows((n_pop+1), quantile(df_1_sim$mu_R_max, 0.025), (n_pop+1), quantile(df_1_sim$mu_R_max, 0.975), angle=90, code=3, length=0.01, lwd=1.5, lty=1, col='goldenrod')
    points(n_pop+1, sim_list$mu_log_R_max_sim, col='black', pch=16)
    plot(n_pop+1, median(df_1_sim$sigma_R_max), col='goldenrod', ylim=c(quantile(df_1_sim$sigma_R_max, 0.0), quantile(df_1_sim$sigma_R_max, 0.995)), xaxt='n', yaxt='n',cex=1)
    arrows((n_pop+1), quantile(df_1_sim$sigma_R_max, 0.25), (n_pop+1), quantile(df_1_sim$sigma_R_max, 0.75), angle=90, code=3, length=0.0, lwd=3, lty=1, col='goldenrod')
    arrows((n_pop+1), quantile(df_1_sim$sigma_R_max, 0.025), (n_pop+1), quantile(df_1_sim$sigma_R_max, 0.975), angle=90, code=3, length=0.01, lwd=1.5, lty=1, col='goldenrod')
    points(n_pop+1, sim_list$sd_log_R_max_sim, col='black', pch=16)
    
    axis(side=4)
    plot(1:n_pop, colMedian(sigma_R), col='saddlebrown', ylim=c(min(cred.2(sigma_R)[1,]),max(cred.2(sigma_R)[2,])), xlim=c(1,n_pop), xaxt='n', yaxt='n',cex=1)
    axis(side=2, las=2)
    mtext(side=2, expression(sigma[R]), line=3.2, cex=1)
    arrows(1:n_pop, cred.2(sigma_R)[1,], 1:n_pop, cred.2(sigma_R)[2,], angle=90, code=3, length=0.01, lwd=1.5, lty=1, col='saddlebrown')
    arrows(1:n_pop, cred.3(sigma_R)[1,], 1:n_pop, cred.3(sigma_R)[2,], angle=90, code=3, length=0.0, lwd=3, lty=1, col='saddlebrown')
    points(1:n_pop, sim_list$sigma_R_sim, col='black', pch=16)
    #axis(side=1, at=c(1:n_pop), labels=levels(coho_dat_full$Population), las=2, cex.axis=1.1)
    plot((n_pop+1), median(exp(df_1_sim$mu_sigma_R)), col='saddlebrown', ylim=c(min(cred.2(sigma_R)[1,]),max(cred.2(sigma_R)[2,])), xaxt='n', yaxt='n', cex=1)
    mtext(side=1, expression(mu), line=1.0, las=1)
    arrows((n_pop+1), quantile(exp(df_1_sim$mu_sigma_R), 0.25), (n_pop+1), quantile(exp(df_1_sim$mu_sigma_R), 0.75), angle=90, code=3, length=0.0, lwd=3, lty=1, col='saddlebrown')
    arrows((n_pop+1), quantile(exp(df_1_sim$mu_sigma_R), 0.025), (n_pop+1), quantile(exp(df_1_sim$mu_sigma_R), 0.975), angle=90, code=3, length=0.01, lwd=1.5, lty=1, col='saddlebrown')
    points(n_pop+1, exp(sim_list$mu_sigma_R_sim), col='black', pch=16)
    #axis(side=4, at=c(log(10), log(25), log(50), log(100), log(250), log(500), log(1000)), labels = c(10,25, 50, 100, 250, 500,1000))
    plot(n_pop+1, median(df_1_sim$sigma_sigma_R), col='saddlebrown', ylim=c(quantile(df_1_sim$sigma_sigma_R, 0.005), quantile(df_1_sim$sigma_sigma_R, 0.99)), xaxt='n', yaxt='n', cex=1)
    mtext(side=1, expression(sigma), line=0.8, las=1)
    arrows((n_pop+1), quantile(df_1_sim$sigma_sigma_R, 0.25), (n_pop+1), quantile(df_1_sim$sigma_sigma_R, 0.75), angle=90, code=3, length=0.0, lwd=3, lty=1, col='saddlebrown')
    arrows((n_pop+1), quantile(df_1_sim$sigma_sigma_R, 0.025), (n_pop+1), quantile(df_1_sim$sigma_sigma_R, 0.975), angle=90, code=3, length=0.01, lwd=1.5, lty=1, col='saddlebrown')
    points(n_pop+1, sim_list$sigma_sigma_R_sim, col='black', pch=16)
    axis(side=4)
    dev.off()
    
    pdf(file='model_fits_sim_AR.pdf')
    #Plot model estimates for smolt abundance
    Smolt_est <- exp(df_1_sim[,grep('log_smolt_est', colnames(df_1_sim))])
    par(mfrow=c(5,2), mar=c(2.5,4,2,1))
    pos <- 1
    for(j in 1:n_pop){
      plot(1:n_year,  colMedian(Smolt_est[,pos:(pos+n_year-1)]), type='o', xlab= 'year', ylab='Smolt estimate', col='red', lwd=1, pch=16, xaxt='n', ylim=c(min(cred.3(Smolt_est[,pos:(pos+n_year-1)])[1,]), max(cred.3(Smolt_est[,pos:(pos+n_year-1)])[2,])))
      axis(side=1, at=c(1:n_year), labels=c(1:n_year)+yr)
      xx <- c(1:n_year,n_year:1)
      yy <- c(cred.2(Smolt_est[,pos:(pos+n_year-1)])[1,], rev(cred.2(Smolt_est[,pos:(pos+n_year-1)])[2,]))
      polygon(x=xx, y=yy, density = -1, border = NA, col=rgb(238, 64, 0, max = 255, alpha = 90))
      yy <- c(cred.3(Smolt_est[,pos:(pos+n_year-1)])[1,], rev(cred.3(Smolt_est[,pos:(pos+n_year-1)])[2,]))
      polygon(x=xx, y=yy, density = -1, border = NA, col=rgb(238, 64, 0, max = 255, alpha = 150))
      points(1:n_year,  sim_list$smolt_est_sim[,j], col='black', pch=16)
      pos <- pos + n_year
    }
    
    #Plot model estimates for Escapement abundance
    esc_est <- exp(df_1_sim[,grep('log_adult_est', colnames(df_1_sim))])
    par(mfrow=c(5,2))
    pos <- 1
    for(j in 1:n_pop){
      plot(1:n_year,  colMedian(esc_est[,pos:(pos+n_year-1)]), type='o', xlab= 'year', ylab='escapement estimate', col='red', lwd=1, pch=16, xaxt='n', ylim=c(min(cred.3(esc_est[,pos:(pos+n_year-1)])[1,]), max(cred.3(esc_est[,pos:(pos+n_year-1)])[2,])))
      axis(side=1, at=c(1:n_year), labels=c(1:n_year)+yr)
      xx <- c(1:n_year,n_year:1)
      yy <- c(cred.2(esc_est[,pos:(pos+n_year-1)])[1,], rev(cred.2(esc_est[,pos:(pos+n_year-1)])[2,]))
      polygon(x=xx, y=yy, density = -1, border = NA, col=rgb(238, 64, 0, max = 255, alpha = 90))
      yy <- c(cred.3(esc_est[,pos:(pos+n_year-1)])[1,], rev(cred.3(esc_est[,pos:(pos+n_year-1)])[2,]))
      polygon(x=xx, y=yy, density = -1, border = NA, col=rgb(238, 64, 0, max = 255, alpha = 150))
      points(1:n_year,  sim_list$adult_est_sim[,j], col='black', pch=16)
      pos <- pos + n_year
    }
    
    #Plot model estimates for Harvest abundance
    harvest_est <- exp(df_1_sim[,grep('log_harvest_est', colnames(df_1_sim))])
    par(mfrow=c(5,2))
    pos <- 1
    for(j in 1:n_pop){
      plot(1:n_year,  colMedian(harvest_est[,pos:(pos+n_year-1)]), type='o', xlab= 'year', ylab='harvest estimate', col='red', lwd=1, pch=16, xaxt='n', ylim=c(min(cred.3(harvest_est[,pos:(pos+n_year-1)])[1,]), max(cred.3(harvest_est[,pos:(pos+n_year-1)])[2,])))
      axis(side=1, at=c(1:n_year), labels=c(1:n_year)+yr)
      xx <- c(1:n_year,n_year:1)
      yy <- c(cred.2(harvest_est[,pos:(pos+n_year-1)])[1,], rev(cred.2(harvest_est[,pos:(pos+n_year-1)])[2,]))
      polygon(x=xx, y=yy, density = -1, border = NA, col=rgb(238, 64, 0, max = 255, alpha = 90))
      yy <- c(cred.3(harvest_est[,pos:(pos+n_year-1)])[1,], rev(cred.3(harvest_est[,pos:(pos+n_year-1)])[2,]))
      polygon(x=xx, y=yy, density = -1, border = NA, col=rgb(238, 64, 0, max = 255, alpha = 150))
      points(1:n_year,  sim_list$harvest_est_sim[,j], col='black', pch=16)
      pos <- pos + n_year
    } 
    
    #Plot model estimates for marine survival
    MS_est <-apply(df_1_sim[,grep('logit_smolt_survival_pop', colnames(df_1_sim))],2 , inv.logit)
    par(mfrow=c(5,2))
    pos <- 1
    for(j in 1:n_pop){
      plot(1:n_year,  colMedian(MS_est[,pos:(pos+n_year-1)]), type='o', xlab= 'year', ylab='Marine survival', col='red', lwd=1, pch=16, xaxt='n', ylim=c(min(cred.3(MS_est[,pos:(pos+n_year-1)])[1,]), max(cred.3(MS_est[,pos:(pos+n_year-1)])[2,])))
      axis(side=1, at=c(1:n_year), labels=c(1:n_year)+yr)
      xx <- c(1:n_year,n_year:1)
      yy <- c(cred.2(MS_est[,pos:(pos+n_year-1)])[1,], rev(cred.2(MS_est[,pos:(pos+n_year-1)])[2,]))
      polygon(x=xx, y=yy, density = -1, border = NA, col=rgb(238, 64, 0, max = 255, alpha = 90))
      yy <- c(cred.3(MS_est[,pos:(pos+n_year-1)])[1,],rev(cred.3(MS_est[,pos:(pos+n_year-1)])[2,]))
      polygon(x=xx, y=yy, density = -1, border = NA, col=rgb(238, 64, 0, max = 255, alpha = 150))
      points(1:n_year,  sim_list$smolt_survival_sim[,j], col='black', pch=16)
      pos <- pos + n_year
    }   
    dev.off()
    
    
    #Plot model estimates relative to true simulated values for single/hyperparameters
    pdf(file='sim_test_single_params_AR.pdf')
    par(mfrow=c(4,2), mar=c(1,2,1,1), oma=c(2,2,2,2))
    post_plot_func(par=df_1_sim$k, par_real=sim_list$k_sim, name='k')
    post_plot_func(par=df_1_sim$`Omega[2,1]`, par_real=sim_list$Omega_cor , name='cor')
    post_plot_func(par=df_1_sim$mu_offset, par_real= sim_list$mu_offset_sim, name='mu_offset')
    post_plot_func(par=df_1_sim$sigma_offset, par_real=sim_list$sigma_offset_sim, name='sigma_offset')
    post_plot_func(par=df_1_sim$mu_mu_surv, par_real=sim_list$mu_mu_surv_sim, name='mu_mu_surv')
    post_plot_func(par=df_1_sim$sigma_mu_surv, par_real=sim_list$sigma_mu_surv_sim, name='sigma_mu_surv')
    post_plot_func(par=df_1_sim$rho_h, par_real=sim_list$rho_h_sim, name='rho_h')
    post_plot_func(par=df_1_sim$sigma_h, par_real=sim_list$sigma_h_sim, name='sigma_h')
    post_plot_func(par=df_1_sim$sigma_smolt, par_real=sim_list$sigma_smolt_sim, name='sigma_smolt')
    post_plot_func(par=df_1_sim$sigma_catch, par_real=sim_list$sigma_catch_sim, name='sigma_catch')
    
    phi_sim <- df_1_sim[,grep('phi', colnames(df_1_sim))]
    for(i in 1:ncol(phi_sim)){
      post_plot_func(par=phi_sim[,i], par_real=sim_list$phi_sim[i] , name=paste('phi', i))
    }
    
    dev.off()
  }
  return(mod_fit_sim)
}

mod_fit_sim <-exec_func_sim(coho_dat_sim=sim_list$coho_dat_sim, n_iter = 10000, n_chains = 5, n_adapt=0.8, n_tree=12.25, n_thin=1, 
                            coho_dat_full = dat_list$coho_dat_full,
                            coord = dat_list$coord,
                            sigma_esc=0.2, exec=TRUE, plots=TRUE)


#Evaluate one step ahead forecast accuracy for a set number of years
forecast_fun <- function(coho_dat_full=dat_list$coho_dat_full[dat_list$coho_dat_full$Calendar.Year<2018,],coord=dat_list$coord, n_range=16, exec=FALSE, exec_RW=FALSE,
                         n_iter = 2000, n_chains = 5, n_adapt=0.8, n_tree=10, n_thin=1,sigma_esc=0.2, plots=TRUE, cast_pop=c(cast_pop),
                         form='RW', mode='state'){
  #Total number of years
  n_year <- length(unique(coho_dat_full$Calendar.Year))
  
  #Total number of populations
  n_pop <- length(unique(coho_dat_full$pop))
  
  #Store the predicted values from the spatial model
  pred_list <-replicate(n_range,list())

  #Store the predicted values from the fitted RW model
  pred_list_RW <-replicate(n_range,list())
  
  #Store the observed values being forecasted
  obs_list <-replicate(n_range,list())

  #Store the denomentator (for calculating for forecasts outside the function), numerator, and ASE for the spatial forecast model
  ASE_den_list  <-replicate(n_range,list())
  ASE_num_list  <-replicate(n_range,list())
  ASE_list  <-replicate(n_range,list())
  ME_list  <-replicate(n_range,list())
  PE_list  <-replicate(n_range,list())
  
  #Store the RMSE and MAPE for the spatial model
  RMSE_list <-replicate(n_range,list())
  MAPE_list <-replicate(n_range,list())
  
  #store the ASE, RMSE, and Mape from the random walk (non-fitted) approach
  RW_ASE_list  <- replicate(n_range,list())
  RW_RMSE_list  <- replicate(n_range,list())
  RW_MAPE_list  <- replicate(n_range,list())
  
  #store the ASE, RMSE, and Mape from the fitted random walk model
  RW_fit_ASE_list  <- replicate(n_range,list())
  RW_fit_ME_list  <- replicate(n_range,list())
  RW_fit_PE_list  <- replicate(n_range,list())
  RW_fit_RMSE_list  <- replicate(n_range,list())
  RW_fit_MAPE_list  <- replicate(n_range,list())
  
  #Add column for total adult returns to the full data set
  coho_dat_full$Total_return <- coho_dat_full$Spawners + coho_dat_full$Harvest
  
  #Create data frame containing only the observations subject to forecasting
  coho_dat_forecast_full <- coho_dat_full[coho_dat_full$yr > (n_year - n_range),]
  
  #Populate list containing the observed data points that are being forecasted
  if(mode=='obs'){
    for(i in 1:n_range){
      obs_list[[i]] <- vector(length=n_pop)
      obs_list[[i]] <- coho_dat_forecast_full$Total_return[coho_dat_forecast_full$yr==unique(coho_dat_forecast_full$yr)[i]]
    }
  }
  
  if(mode=='state'){
    for(i in 1:n_range){
      obs_list[[i]] <- vector(length=n_pop)
    }
    mod_fit <- readRDS('LD_coho_forecast_6_2_4c_fit.rds')
    df_1 <- as.data.frame(mod_fit)
    esc_est <- exp(df_1[,grep('log_adult_est', colnames(df_1))])
    harvest_est <- exp(df_1[,grep('log_harvest_est', colnames(df_1))])
    adult_est <- esc_est + harvest_est
    
    pos <- 1
      for(j in 1:n_pop){
        for(i in 1:n_range){
          obs_list[[i]][j] <- colMedian(adult_est[,pos:(pos+n_year-1)])[(n_year-n_range+1):n_year][i]
        }
        pos <- pos + n_year + 1 #the extra +1 is because the full model fit includes data up to 2018 but only using up to 2017 here
        #pos <- pos + n_year #if the year endpoints were to match
      }
  }
  
  #Forecast routine
  for(i in 1:n_range){
    #Alternate counter for reversing order
    I <- rev(1:n_range)[i]
    
    #Subset the data to represent what would be available to the forecaster in a given year
    #First remove all data past the current year
    coho_dat_forecast <- subset(coho_dat_full, coho_dat_full$yr <= (n_year - i ))
    coho_dat_forecast <- droplevels(coho_dat_forecast)
    
    #Add empty vectors to each list element for populating
    #ASE denomenator
    ASE_den_list[[I]] <- vector(length=n_pop)
    #ASE numerator, ASE, RMSE, and MAPE for fitted spatial model
    ASE_num_list[[I]] <- vector(length=n_pop)
    ME_list[[I]] <- vector(length=n_pop)
    PE_list[[I]] <- vector(length=n_pop)
    
    ASE_list[[I]] <- vector(length=n_pop)
    RMSE_list[[I]] <- vector(length=n_pop)
    MAPE_list[[I]] <- vector(length=n_pop)
    
    #ASE, RMSE, and Mape for non-fitted RW approach
    RW_ASE_list[[I]] <- vector(length=n_pop)
    RW_RMSE_list[[I]] <- vector(length=n_pop)
    RW_MAPE_list[[I]] <- vector(length=n_pop)
    
    #ASE, RMSE, and MAPE for fitted RW approach
    RW_fit_ASE_list[[I]] <- vector(length=n_pop)
    RW_fit_ME_list[[I]] <- vector(length=n_pop)
    RW_fit_PE_list[[I]] <- vector(length=n_pop)
    
    RW_fit_RMSE_list[[I]] <- vector(length=n_pop)
    RW_fit_MAPE_list[[I]] <- vector(length=n_pop)
    
    #Calculate the denomenator of the ASE using the training data set (updated with each forecast)
    if(mode=='obs'){
      for(p in 1:n_pop){
        ASE_den_list[[I]][p] <-  mean(abs(diff(na.omit(coho_dat_forecast$Total_return[coho_dat_forecast$pop==p]))))
      }
    }
    if(mode=='state'){
      state_mat <- matrix(colMedian(adult_est), ncol=n_pop, nrow=n_year+1 #+1 because full fit is 33 years, only 32 within func
      ,byrow=FALSE)[1:(n_year-i),] 
      
      for(p in 1:n_pop){
        ASE_den_list[[I]][p] <- mean(abs(diff(state_mat[,p])))
      }
    }
    
    #Forecast assumes that the current year's smolt abundance would not be available in the database yet (but it might be? Probably doesn't matter...)
    #The current year's smolt survival estimate would not be available
    coho_dat_forecast$Release_No[coho_dat_forecast$yr==length(unique(coho_dat_forecast$yr))] <- NA
    coho_dat_forecast$Fishery_Plus_Escapement[coho_dat_forecast$yr==length(unique(coho_dat_forecast$yr))] <- NA
    
    #Optional --> might not have last year's CWT data ready in FRAM yet, remove one more year
    coho_dat_forecast$Release_No[coho_dat_forecast$yr==length(unique(coho_dat_forecast$yr))-1] <- NA
    coho_dat_forecast$Fishery_Plus_Escapement[coho_dat_forecast$yr==length(unique(coho_dat_forecast$yr))-1] <- NA
    
    #Previous Year's harvest and escapement might not be available either
    coho_dat_forecast$Harvest[coho_dat_forecast$yr==length(unique(coho_dat_forecast$yr))] <- NA
    #coho_dat_forecast$Spawners[coho_dat_forecast$yr==length(unique(coho_dat_forecast$yr))] <- NA
    
    coho_dat_forecast <- droplevels(coho_dat_forecast)
    
    #Run the model with the data removed to generate forecasts.
    mod_fit_pred <- exec_fun(coho_dat_full=coho_dat_forecast,
                             coord=coord, exec=exec, n_iter=n_iter, n_thin=n_thin, n_adapt=n_adapt, n_tree=n_tree, n_chain=n_chains, plots=FALSE,
                             sigma_esc = sigma_esc, mod_filename=paste('LD_coho_forecast_6_2_4_fit_yr=',i,'_',n_range, sep=''))
    
    #Run the random walk model
    mod_fit_RW <- exec_fun_RW(coho_dat_full=coho_dat_forecast, exec=exec_RW, n_iter=n_iter, n_thin=n_thin, n_adapt=n_adapt, n_tree=n_tree, n_chain=n_chains,
                              form=form, mod_filename=paste('LD_simple_',form,'_fit_yr=',i,'_',n_range,sep=''))   
    
    #Save the spatial model to data frame and grab the prediction vector
    df_pred <- as.data.frame(mod_fit_pred)
    pred_list[[I]] <- df_pred[,grep('adult_pred', colnames(df_pred))]
    #Save the RW model to a data frame and grab the random walk prediction vector
    df_RW <- as.data.frame(mod_fit_RW)
    pred_list_RW[[I]] <- df_RW[,grep('adult_pred', colnames(df_RW))]
    
    #Now calculate the observed versus predicted (numerator of ASE) and the complete ASE, as well as RMSE and MAPE for all forecasts
    for(p in 1:n_pop){
      #Prediction from Random walk (non-fitted) based on observed return 1 year earlier
      RW_pred <- coho_dat_full$Total_return[coho_dat_full$pop==p & coho_dat_full$yr==(n_year-i)]

      #If previous observation is NA, grab the most recent non-NA value
      if (is.na(RW_pred)==TRUE){
        RW_pred <- na.omit(coho_dat_full$Total_return[coho_dat_full$pop==p & coho_dat_full$yr < (n_year-i)])[length(na.omit(coho_dat_full$Total_return[coho_dat_full$pop==p & coho_dat_full$yr < (n_year-i)]))]
      }
      #Spatial model
      ASE_num_list[[I]][p] <- abs(colMedian(pred_list[[I]])[p] - obs_list[[I]][p])
      ME_list[[I]][p] <- colMedian(pred_list[[I]])[p] - obs_list[[I]][p]
      PE_list[[I]][p] <- (colMedian(pred_list[[I]])[p] - obs_list[[I]][p])/obs_list[[I]][p]
      ASE_list[[I]][p] <- ASE_num_list[[I]][p]/ASE_den_list[[I]][p]
      MAPE_list[[I]][p] <- abs(log(colMedian(pred_list[[I]])[p]/obs_list[[I]][p]))
      RMSE_list[[I]][p] <-  (colMedian(pred_list[[I]])[p] - obs_list[[I]][p])^2
      
      #non-fitted random walk
      RW_ASE_list[[I]][p] <-  abs(RW_pred - obs_list[[I]][p])/ASE_den_list[[I]][p]
      RW_MAPE_list[[I]][p] <- abs(log(RW_pred/obs_list[[I]][p]))
      RW_RMSE_list[[I]][p] <-  (RW_pred - obs_list[[I]][p])^2
      
      #fitted random walk
      RW_fit_ASE_list[[I]][p] <-  abs(colMedian(pred_list_RW[[I]])[p] - obs_list[[I]][p])/ASE_den_list[[I]][p]
      RW_fit_ME_list[[I]][p] <-  colMedian(pred_list_RW[[I]])[p] - obs_list[[I]][p]
      RW_fit_PE_list[[I]][p] <-  (colMedian(pred_list_RW[[I]])[p] - obs_list[[I]][p])/obs_list[[I]][p]
      RW_fit_MAPE_list[[I]][p] <- abs(log(colMedian(pred_list_RW[[I]])[p]/obs_list[[I]][p]))
      RW_fit_RMSE_list[[I]][p] <-  (colMedian(pred_list_RW[[I]])[p] - obs_list[[I]][p])^2
      
    }
  }
  
  if(plots==TRUE){
    pdf(file=paste('forecast_plots_AR',mode,'.pdf'))
    mod_fit_1 <- readRDS(paste('LD_coho_forecast_6_2_4_fit_yr=',n_range,'_',n_range,'.rds',sep=''))
    df_1 <- as.data.frame(mod_fit_1)
    n_year_1 <- length(unique(coho_dat_full$Calendar.Year)) - n_range
    yr <- min(coho_dat_full$Calendar.Year)-1
    n_pop_1 <- length(unique(coho_dat_full$pop))
    #Plot model estimates for smolt abundance
    total_adult_est <- exp(df_1[,grep('log_adult_est', colnames(df_1))]) + exp(df_1[,grep('log_harvest_est', colnames(df_1))])
    par(mfrow=c(9,4), mar=c(1,2.25,0.75,0), oma=c(2,2,0,0.5), mgp=c(3,0.5,0))
    pos <- 1
    J <- 0
    for(j in 1:n_pop_1){
      
      if(j %in% cast_pop){
        pred_pop <- vector(length=n_range)
        pred_pop_lower_50 <- vector(length=n_range)
        pred_pop_upper_50 <- vector(length=n_range)
        pred_pop_lower_95 <- vector(length=n_range)
        pred_pop_upper_95 <- vector(length=n_range)
        for(r in 1:n_range){
          pred_pop[r] <- colMedian(pred_list[[r]])[j]
          pred_pop_lower_50[r] <- cred.3(pred_list[[r]][j])[1]
          pred_pop_upper_50[r] <- cred.3(pred_list[[r]][j])[2]
          pred_pop_lower_95[r] <- cred.2(pred_list[[r]][j])[1]
          pred_pop_upper_95[r] <- cred.2(pred_list[[r]][j])[2]
        }
        
        y_lim <- c(NA, NA)
        y_lim[1] <- min(pred_pop_lower_95, cred.2(total_adult_est[,pos:(pos+n_year_1-1)])[1,], min(coho_dat_full$Total_return[coho_dat_full$pop==j], na.rm=TRUE))
        y_lim[2] <- max(pred_pop_upper_95, cred.2(total_adult_est[,pos:(pos+n_year_1-1)])[2,], max(coho_dat_full$Total_return[coho_dat_full$pop==j], na.rm=TRUE))
        
        plot(1:(n_year_1),  colMedian(total_adult_est[,pos:(pos+n_year_1-1)]), type='o', xlab= NA, ylab=NA, cex=1, axes=F, col='dodgerblue4', lwd=1, xaxs='i', yaxs='i', pch=16, xaxt='n',xlim=c(1,n_year), ylim=c(y_lim))
        points(1:(n_year_1),  colMedian(total_adult_est[,pos:(pos+n_year_1-1)]), col='black', lwd=0.5, pch=1, cex=1)
        box(col='grey', lwd=0.5)
        xx <- c(1:(n_year_1),(n_year_1):1)
        yy <- c(cred.2(total_adult_est[,pos:(pos+n_year_1-1)])[1,], rev(cred.2(total_adult_est[,pos:(pos+n_year_1-1)])[2,]))
        polygon(x=xx, y=yy, density = -1, border = 'dodgerblue4', lwd=0.85, col=rgb(13, 82, 144, max = 255, alpha = 90))
        yy <- c(cred.3(total_adult_est[,pos:(pos+n_year_1-1)])[1,], rev(cred.3(total_adult_est[,pos:(pos+n_year_1-1)])[2,]))
        polygon(x=xx, y=yy, density = -1, border = 'dodgerblue4', lwd=0.85, col=rgb(13, 82, 144, max = 255, alpha = 150))
        points(1:(n_year_1),  colMedian(total_adult_est[,pos:(pos+n_year_1-1)]), type='o', cex=1, col='dodgerblue4', lwd=1, pch=16)
        points(1:(n_year_1),  colMedian(total_adult_est[,pos:(pos+n_year_1-1)]), col='black', lwd=0.5, pch=1, cex=1)
        #abline(v=n_year_1+0.5, col='gold3', lty=2, lwd=1)
        axis(side=1, at=c(1:n_year), labels=c(1:n_year)+yr, cex.axis=0.8, lwd=0.5)
        axis(side=2, cex.axis=0.8, lwd=0.5)
        
        if(mode=='obs'){
          points(coho_dat_full$yr[coho_dat_full$pop==j & coho_dat_full$yr <= n_year_1], coho_dat_full$Total_return[coho_dat_full$pop==j & coho_dat_full$yr <= n_year_1], col='black', pch=16, cex=1)
        }

        xx <- c((n_year_1+1):n_year,rev((n_year_1+1):n_year))
        yy <- c(pred_pop_lower_50, rev(pred_pop_upper_50))
        polygon(x=xx, y=yy, density = -1, border = 'gold3', lwd=0.85, col=rgb(212, 167, 33, max = 255, alpha = 150))
        yy <- c(pred_pop_lower_95, rev(pred_pop_upper_95))
        polygon(x=xx, y=yy, density = -1, border = 'gold3', lwd=0.85, col=rgb(212, 167, 33, max = 255, alpha = 90))
        points((n_year_1+1):n_year, pred_pop, col='gold3', type='o', pch=16, lwd=1.0, cex=1)
        points((n_year_1+1):n_year, pred_pop, col='black', pch=1, lwd=0.5, cex=1)
        
        if(mode=='obs'){
          points((n_year_1+1):n_year, unlist(lapply(obs_list, function(x) x[j])), pch=16, cex=1)
        }
        
        if(mode=='state'){
          points((n_year_1+1):n_year, unlist(lapply(obs_list, function(x) x[j])), type='o', lty=1, lwd=1, pch=16, cex=1, col='dodgerblue4')
        }

        if(j==1){
          if(mode=='obs'){
          legend(x=1.5, y=y_lim[2]*1.07, legend=c('Training', 'Forecast', 'Observed'), col=c('dodgerblue4', 'gold3', 'black'),
                 pch=16, lty=c(1,1,NA), bty='n', y.intersp=0.75, x.intersp=0.25, cex=0.8, pt.cex=1)
          legend(x=1.5, y=y_lim[2]*1.07, legend=c(NA,NA,NA), col=c('black'),
                 pch=1, lwd=0.5, lty=0, bty='n', y.intersp=0.75, x.intersp=0.25, cex=0.8, pt.cex=1)
          }
          if(mode=='state'){
            legend(x=1.5, y=y_lim[2]*1.07, legend=c('State', 'Forecast'), col=c('dodgerblue4', 'gold3'),
                   pch=16, lty=c(1,1,NA), bty='n', y.intersp=0.75, x.intersp=0.25, cex=0.8, pt.cex=1)
            legend(x=1.5, y=y_lim[2]*1.07, legend=c(NA,NA), col=c('black'),
                   pch=1, lwd=0.5, lty=0, bty='n', y.intersp=0.75, x.intersp=0.25, cex=0.8, pt.cex=1)
          }
        }
        mtext(side=3, unique(coho_dat_full$Name)[j], cex=0.6, line=-0.9, at=18)
        
        if(j==cast_pop[17]){
          mtext(side=2, 'Total return', line=1.75, cex=0.9)
        }
        if(j==cast_pop[35]){
          mtext(side=1, 'Year', line=1.75, at=-3, cex=0.9)
        }
      }
      pos <- pos + n_year_1
    }
    dev.off()
  }
  
  master_list <- list()
  #forecasts from the spatial model
  master_list$preds <- pred_list
  
  #ASE denomenator
  master_list$ASE_den <- ASE_den_list
  
  #ASE numerator, ASE, MAPE, and RMSE for spatial model
  master_list$ASE_num <- ASE_num_list
  master_list$ME <- ME_list
  master_list$PE <- PE_list
  master_list$ASE <- ASE_list
  master_list$RMSE <- RMSE_list
  master_list$MAPE <- MAPE_list
  
  #ASE, RMSE, and MAPE for non-fitted random walk approach
  master_list$RW_ASE <- RW_ASE_list
  master_list$RW_RMSE <- RW_RMSE_list
  master_list$RW_MAPE <- RW_MAPE_list
  
  #ASE, RMSE, and MAPE for fitted RW approach
  master_list$RW_fit_ASE <- RW_fit_ASE_list
  master_list$RW_fit_ME <- RW_fit_ME_list
  master_list$RW_fit_PE <- RW_fit_PE_list
  master_list$RW_fit_RMSE <- RW_fit_RMSE_list
  master_list$RW_fit_MAPE <- RW_fit_MAPE_list
  
  #Observed data being forecasted
  master_list$obs <- obs_list
  
  return(master_list)
}

coho_forecast_RW <- forecast_fun(coho_dat_full=dat_list$coho_dat_full[dat_list$coho_dat_full$Calendar.Year<2018,],coord=dat_list$coord, n_range=16, exec=FALSE, exec_RW=FALSE,
                              n_iter = 2000, n_chains = 5, n_adapt=0.8, n_tree=10, n_thin=1, sigma_esc=0.2, plots=TRUE, cast_pop=c(cast_pop),
                              form='RW', mode='state')

coho_forecast_AR <- forecast_fun(coho_dat_full=dat_list$coho_dat_full[dat_list$coho_dat_full$Calendar.Year<2018,],coord=dat_list$coord, n_range=16, exec=FALSE, exec_RW=FALSE,
                                 n_iter = 2000, n_chains = 5, n_adapt=0.8, n_tree=10, n_thin=1, sigma_esc=0.2, plots=TRUE, cast_pop=c(cast_pop),
                                 form='AR', mode='state')

coho_forecast_MA <- forecast_fun(coho_dat_full=dat_list$coho_dat_full[dat_list$coho_dat_full$Calendar.Year<2018,],coord=dat_list$coord, n_range=16, exec=FALSE, exec_RW=FALSE,
                                 n_iter = 2000, n_chains = 5, n_adapt=0.8, n_tree=10, n_thin=1,sigma_esc=0.2, plots=TRUE, cast_pop=c(cast_pop),
                                 form='MA', mode='state')


#Forecast performence statistics
#1) (M)ASE####
#Calculate the mean for each approach across all populations
mean(as.vector(unlist(coho_forecast_RW$ASE)), na.rm=TRUE)
mean(as.vector(unlist(coho_forecast_RW$ME)), na.rm=TRUE)
length(which((as.vector(unlist(coho_forecast_RW$PE))) > 0.5))
length(which((as.vector(unlist(coho_forecast_RW$PE))) < -0.5))

mean(as.vector(unlist(coho_forecast_RW$RW_fit_ASE)), na.rm=TRUE)
mean(as.vector(unlist(coho_forecast_RW$RW_fit_ME)), na.rm=TRUE)
length(which((as.vector(unlist(coho_forecast_RW$RW_fit_PE))) > 0.5))
length(which((as.vector(unlist(coho_forecast_RW$RW_fit_PE))) < -0.5))

mean(as.vector(unlist(coho_forecast_AR$RW_fit_ASE)), na.rm=TRUE)
mean(as.vector(unlist(coho_forecast_AR$RW_fit_ME)), na.rm=TRUE)
length(which((as.vector(unlist(coho_forecast_AR$RW_fit_PE))) > 0.5))
length(which((as.vector(unlist(coho_forecast_AR$RW_fit_PE))) < -0.5))

mean(as.vector(unlist(coho_forecast_MA$RW_fit_ASE)), na.rm=TRUE)
mean(as.vector(unlist(coho_forecast_MA$RW_fit_ME)), na.rm=TRUE)
length(which((as.vector(unlist(coho_forecast_MA$RW_fit_PE))) > 0.5))
length(which((as.vector(unlist(coho_forecast_MA$RW_fit_PE))) < -0.5))

#Calculate forecast accuracy only for the FRAM units included in the MUs of the Published coho forecast
#Includes the Published method
n_range <- 16
n_pop <- 36
ASE_mat_ST <- matrix(ncol=n_pop, nrow=n_range)
ASE_mat_RW <- matrix(ncol=n_pop, nrow=n_range)
ASE_mat_RW_fit <- matrix(ncol=n_pop, nrow=n_range)
ASE_mat_AR <- matrix(ncol=n_pop, nrow=n_range)
ASE_mat_MA <- matrix(ncol=n_pop, nrow=n_range)

for(i in 1:n_pop){
  ASE_mat_ST[,i] <- unlist(lapply(coho_forecast_RW$ASE, function(x) x[i]))
  ASE_mat_RW[,i] <- unlist(lapply(coho_forecast_RW$RW_ASE, function(x) x[i]))
  ASE_mat_RW_fit[,i] <- unlist(lapply(coho_forecast_RW$RW_fit_ASE, function(x) x[i]))
  ASE_mat_AR[,i] <- unlist(lapply(coho_forecast_AR$RW_fit_ASE, function(x) x[i]))
  ASE_mat_MA[,i] <- unlist(lapply(coho_forecast_MA$RW_fit_ASE, function(x) x[i]))
}
cast_pop <- cast_pop
mean(ASE_mat_RW_fit[,cast_pop], na.rm=TRUE)
mean(ASE_mat_AR[,cast_pop], na.rm=TRUE)
mean(ASE_mat_MA[,cast_pop], na.rm=TRUE)
mean(ASE_mat_ST[,cast_pop], na.rm=TRUE)

#For the Published method
cast_pop <- cast_pop
ASE_WDFW <- abs(state_cast[,-c(1,18,19,20)]-matrix(unlist(coho_forecast_RW$obs), nrow=n_pop, ncol=n_range, byrow=FALSE)[cast_pop,])/matrix(unlist(coho_forecast_RW$ASE_den), nrow=n_pop, ncol=n_range, byrow=FALSE)[cast_pop,]
PE_WDFW <- (state_cast[,-c(1,18,19,20)]-matrix(unlist(coho_forecast_RW$obs), nrow=n_pop, ncol=n_range, byrow=FALSE)[cast_pop,])/matrix(unlist(coho_forecast_RW$obs), nrow=n_pop, ncol=n_range, byrow=FALSE)[cast_pop,]
ME_WDFW <- (state_cast[,-c(1,18,19,20)]-matrix(unlist(coho_forecast_RW$obs), nrow=n_pop, ncol=n_range, byrow=FALSE)[cast_pop,])

mean(unlist(ASE_WDFW), na.rm=TRUE)
mean(unlist(ME_WDFW), na.rm=TRUE)

length(which(PE_WDFW > 0.5))
length(which(PE_WDFW<  -0.5))


#Plot relative performence of each method according to MASE
pdf(file='forecast_summary_plot_MASE.pdf')
col.vec <- brewer.pal(n=9, 'YlOrBr')
col.vec_2 <- col2rgb(col.vec)
par(oma=c(8,12,10,0), mar=c(5,1.25,0,0), mfrow=c(1,5))
barplot(apply(ASE_WDFW, 1, mean, na.rm=TRUE), horiz=TRUE, las=2, col='black', names.arg=unique(dat_list$coho_dat_full$Name)[cast_pop], cex.names=1, border=NA, xlim=c(0,2.5))
barplot(apply(ASE_mat_ST[,cast_pop], 2, mean, na.rm=TRUE), add=TRUE, horiz=TRUE, las=2,col=rgb(17,156,240, alpha=150, max=255 ), border=NA)
text(x=1.75, y=1, 'A')
mtext(side=3, 'ST-IPM', line=1.25, cex=0.85)
mtext(side=3, 'vs.', line=0.25, cex=0.85)
mtext(side=3, 'Published', line=-0.75, cex=0.85)

barplot(apply(ASE_mat_RW[,cast_pop], 2, mean, na.rm=TRUE),col=rgb(red=col.vec_2[1,9],green=col.vec_2[2,9],blue=col.vec_2[3,9], alpha=200, max=255 ), horiz=TRUE, las=2, cex.names=0.7, border=NA, xlim=c(0,2.5))
barplot(apply(ASE_mat_ST[,cast_pop], 2, mean, na.rm=TRUE), add=TRUE, horiz=TRUE, las=2,col=rgb(17,156,240, alpha=150, max=255 ), border=NA)
mtext(side=1, 'MASE', line=2.5, cex=0.85, at=2.25)
text(x=1.75, y=1, 'B')
mtext(side=3, 'ST-IPM', line=1.25,cex=0.85)
mtext(side=3, 'vs.', line=0.25, cex=0.85)
mtext(side=3, 'Random walk', line=-0.75, cex=0.85)

#barplot(apply(ASE_mat_RW_fit[,cast_pop], 2, mean, na.rm=TRUE), horiz=TRUE, las=2, col=rgb(red=col.vec_2[1,8],green=col.vec_2[2,8],blue=col.vec_2[3,8], alpha=200, max=255 ), border=NA)
#barplot(apply(ASE_mat_ST[,cast_pop], 2, mean, na.rm=TRUE), add=TRUE, horiz=TRUE, las=2,col=rgb(17,156,240, alpha=150, max=255 ), border=NA)
#mtext(side=1, 'MASE', line=2.5, cex=0.85)

barplot(apply(ASE_mat_AR[,cast_pop], 2, mean, na.rm=TRUE), horiz=TRUE, las=2, col=rgb(red=col.vec_2[1,6],green=col.vec_2[2,6],blue=col.vec_2[3,6], alpha=200, max=255 ), border=NA, xlim=c(0,2.5))
barplot(apply(ASE_mat_ST[,cast_pop], 2, mean, na.rm=TRUE), add=TRUE, horiz=TRUE, las=2,col=rgb(17,156,240, alpha=150, max=255 ), border=NA)
text(x=1.75, y=1, 'C')
mtext(side=3, 'ST-IPM', line=1.25, cex=0.85)
mtext(side=3, 'vs.', line=0.25, cex=0.85)
mtext(side=3, 'Autoregressive', line=-0.75, cex=0.85)

barplot(apply(ASE_mat_MA[,cast_pop], 2, mean, na.rm=TRUE), horiz=TRUE, las=2, col=rgb(red=col.vec_2[1,4],green=col.vec_2[2,4],blue=col.vec_2[3,4], alpha=200, max=255 ), border=NA, xlim=c(0,2.5))
barplot(apply(ASE_mat_ST[,cast_pop], 2, mean, na.rm=TRUE), add=TRUE, horiz=TRUE, las=2,col=rgb(17,156,240, alpha=150, max=255 ), border=NA)
text(x=1.75, y=1, 'D')
mtext(side=3, 'ST-IPM', line=1.25, cex=0.85)
mtext(side=3, 'vs.', line=0.25, cex=0.85)
mtext(side=3, 'Moving average', line=-0.75, cex=0.85)

legend(x=0.85, y=10, pch=15, col=c(col=rgb(17,156,240, alpha=150, max=255 ), 'black', 
                                   rgb(red=col.vec_2[1,8],green=col.vec_2[2,8],blue=col.vec_2[3,8], alpha=200, max=255),
                                   rgb(red=col.vec_2[1,6],green=col.vec_2[2,6],blue=col.vec_2[3,6], alpha=200, max=255),
                                   rgb(red=col.vec_2[1,4],green=col.vec_2[2,4],blue=col.vec_2[3,4], alpha=200, max=255)),
       legend=c('ST-IPM', 'Published', 'RW', 'AR', 'MA'), bty='n', pt.cex=1.5, cex=0.9)

dev.off()
#####

#Published versus ST-IPM only
#Plot relative performence of each method according to MASE
pdf(file='forecast_summary_plot_MASE_only.pdf')
col.vec <- brewer.pal(n=9, 'YlOrBr')
col.vec_2 <- col2rgb(col.vec)
par(oma=c(8,12,10,0), mar=c(5,1.25,0,0), mfrow=c(1,5))
barplot(apply(ASE_WDFW, 1, mean, na.rm=TRUE), horiz=TRUE, las=2, col='black', names.arg=unique(dat_list$coho_dat_full$Name)[cast_pop], cex.names=1, border=NA, xlim=c(0,2.3))
barplot(apply(ASE_mat_ST[,cast_pop], 2, mean, na.rm=TRUE), add=TRUE, horiz=TRUE, las=2,col=rgb(17,156,240, alpha=150, max=255 ), border=NA)
legend(x=0.6, y=4, pch=15, col=c(col=rgb(17,156,240, alpha=150, max=255 ), 'black', 
                                   rgb(red=col.vec_2[1,8],green=col.vec_2[2,8],blue=col.vec_2[3,8], alpha=200, max=255)),
       legend=c('ST-IPM', 'Published'), bty='n', pt.cex=2, cex=1)
mtext(side=1, 'MASE', line=2.7, cex=0.85)
dev.off()

#MASE plot B --> sparklines of ASE over time 
pdf(file='ASE_sparklines.pdf')
par(mfrow=c(length(cast_pop), 1), oma=c(2.5,20,0.1,20), mar=c(0,0,0,0))
for(i in 1:length(cast_pop)){
  plot(ASE_mat_ST[,cast_pop][,i], type='o', cex=0.6, pch=16, xaxt='n', lwd=1.2, col=rgb(17,156,240, alpha=255, max=255 ),
                                                  ylim=c(0, max(max(ASE_mat_ST[,i], na.rm=TRUE),
                                                                max(ASE_WDFW[i,], na.rm=TRUE))))
  points(1:n_range, ASE_WDFW[i,], type='o', cex=0.6, pch=16, xaxt='n', lwd=1.2, col='black')
  mtext(side=4, unique(dat_list$coho_dat_full$Name)[cast_pop][i], las=2, line=0.25, cex=0.6)
}
axis(side=1, at=c(2,6,10,14), labels=c(2,6,10,14) + (2017-n_range))
mtext(side=2, 'ASE', at=62.5, line=2.5, cex=0.85)
dev.off()

#2) RMSE####
#Calculate the mean for each approach across all populations
sqrt(mean(as.vector(unlist(coho_forecast_RW$RMSE)), na.rm=TRUE))
sqrt(mean(as.vector(unlist(coho_forecast_RW$RW_RMSE)), na.rm=TRUE))
sqrt(mean(as.vector(unlist(coho_forecast_RW$RW_fit_RMSE)), na.rm=TRUE))
sqrt(mean(as.vector(unlist(coho_forecast_AR$RW_fit_RMSE)), na.rm=TRUE))
sqrt(mean(as.vector(unlist(coho_forecast_MA$RW_fit_RMSE)), na.rm=TRUE))

#Calculate forecast accuracy only for the FRAM units included in the MUs of the Published coho forecast
#Includes the Published method
n_range <- 16
n_pop <- 36
cast_pop <- cast_pop
RMSE_mat_ST <- matrix(ncol=n_pop, nrow=n_range)
RMSE_mat_RW <- matrix(ncol=n_pop, nrow=n_range)
RMSE_mat_RW_fit <- matrix(ncol=n_pop, nrow=n_range)
RMSE_mat_AR <- matrix(ncol=n_pop, nrow=n_range)
RMSE_mat_MA <- matrix(ncol=n_pop, nrow=n_range)

for(i in 1:n_pop){
  RMSE_mat_ST[,i] <- unlist(lapply(coho_forecast_RW$RMSE, function(x) x[i]))
  RMSE_mat_RW[,i] <- unlist(lapply(coho_forecast_RW$RW_RMSE, function(x) x[i]))
  RMSE_mat_RW_fit[,i] <- unlist(lapply(coho_forecast_RW$RW_fit_RMSE, function(x) x[i]))
  RMSE_mat_AR[,i] <- unlist(lapply(coho_forecast_AR$RW_fit_RMSE, function(x) x[i]))
  RMSE_mat_MA[,i] <- unlist(lapply(coho_forecast_MA$RW_fit_RMSE, function(x) x[i]))
}

sqrt(mean(RMSE_mat_RW_fit[,cast_pop], na.rm=TRUE))
sqrt(mean(RMSE_mat_AR[,cast_pop], na.rm=TRUE))
sqrt(mean(RMSE_mat_MA[,cast_pop], na.rm=TRUE))
sqrt(mean(RMSE_mat_ST[,cast_pop], na.rm=TRUE))

#For the Published method
RMSE_WDFW <- (state_cast[,-c(1,18,19,20)]-matrix(unlist(coho_forecast_RW$obs), nrow=n_pop, ncol=n_range, byrow=FALSE)[cast_pop,])^2
sqrt(mean(unlist(RMSE_WDFW), na.rm=TRUE))

#Create sqrt_mean_function
sqrt_mean <- function(x){
  return(sqrt(mean(x, na.rm=TRUE)))
}

#Plot relative performence of each method according to RMSE
pdf(file='forecast_summary_plot_RMSE.pdf')
col.vec <- brewer.pal(n=9, 'YlOrBr')
col.vec_2 <- col2rgb(col.vec)
par(oma=c(8,12,10,1), mar=c(5,0.25,0,0), mfrow=c(1,5))
barplot(apply(RMSE_WDFW, 1, sqrt_mean), horiz=TRUE, las=2, col='black', names.arg=unique(dat_list$coho_dat_full$Name)[cast_pop], cex.names=1, border=NA)
barplot(apply(RMSE_mat_ST[,cast_pop], 2, sqrt_mean), add=TRUE, horiz=TRUE, las=2,col=rgb(17,156,240, alpha=150, max=255 ), border=NA)
text(115000, 1, 'A')
mtext(side=3, 'ST-IPM', line=1.25, cex=0.85)
mtext(side=3, 'vs.', line=0.25, cex=0.85)
mtext(side=3, 'Published', line=-0.75, cex=0.85)

barplot(apply(RMSE_mat_RW[,cast_pop], 2, sqrt_mean),col=rgb(red=col.vec_2[1,9],green=col.vec_2[2,9],blue=col.vec_2[3,9], alpha=200, max=255 ), horiz=TRUE, las=2, cex.names=0.7, border=NA)
barplot(apply(RMSE_mat_ST[,cast_pop], 2, sqrt_mean), add=TRUE, horiz=TRUE, las=2,col=rgb(17,156,240, alpha=150, max=255 ), border=NA)
mtext(side=1, 'RMSE', line=3.5, cex=0.85, at=90000)
text(55000, 1, 'B')
mtext(side=3, 'ST-IPM', line=1.25, cex=0.85)
mtext(side=3, 'vs.', line=0.25, cex=0.85)
mtext(side=3, 'Random walk', line=-0.75, cex=0.85)

#barplot(apply(RMSE_mat_RW_fit[,cast_pop], 2, sqrt_mean, na.rm=TRUE), horiz=TRUE, las=2, col=rgb(red=col.vec_2[1,8],green=col.vec_2[2,8],blue=col.vec_2[3,8], alpha=200, max=255 ), border=NA)
#barplot(apply(RMSE_mat_ST[,cast_pop], 2, sqrt_mean, na.rm=TRUE), add=TRUE, horiz=TRUE, las=2,col=rgb(17,156,240, alpha=150, max=255 ), border=NA)
#mtext(side=1, 'MRMSE', line=2.5, cex=0.85)

barplot(apply(RMSE_mat_AR[,cast_pop], 2, sqrt_mean), horiz=TRUE, las=2, col=rgb(red=col.vec_2[1,6],green=col.vec_2[2,6],blue=col.vec_2[3,6], alpha=200, max=255 ), border=NA)
barplot(apply(RMSE_mat_ST[,cast_pop], 2, sqrt_mean), add=TRUE, horiz=TRUE, las=2,col=rgb(17,156,240, alpha=150, max=255 ), border=NA)
text(65000, 1, 'C')
mtext(side=3, 'ST-IPM', line=1.25, cex=0.85)
mtext(side=3, 'vs.', line=0.25, cex=0.85)
mtext(side=3, 'Autoregressive', line=-0.75, cex=0.85)

barplot(apply(RMSE_mat_MA[,cast_pop], 2, sqrt_mean), horiz=TRUE, las=2, col=rgb(red=col.vec_2[1,4],green=col.vec_2[2,4],blue=col.vec_2[3,4], alpha=200, max=255 ), border=NA)
barplot(apply(RMSE_mat_ST[,cast_pop], 2, sqrt_mean), add=TRUE, horiz=TRUE, las=2,col=rgb(17,156,240, alpha=150, max=255 ), border=NA)
text(65000, 1, 'D')
mtext(side=3, 'ST-IPM', line=1.25, cex=0.85)
mtext(side=3, 'vs.', line=0.25, cex=0.85)
mtext(side=3, 'Moving average', line=-0.75, cex=0.85)

legend(x=20000, y=34, pch=15, col=c(col=rgb(17,156,240, alpha=150, max=255 ), 'black', 
                                  rgb(red=col.vec_2[1,8],green=col.vec_2[2,8],blue=col.vec_2[3,8], alpha=200, max=255),
                                  rgb(red=col.vec_2[1,6],green=col.vec_2[2,6],blue=col.vec_2[3,6], alpha=200, max=255),
                                  rgb(red=col.vec_2[1,4],green=col.vec_2[2,4],blue=col.vec_2[3,4], alpha=200, max=255)),
       legend=c('ST-IPM', 'Published', 'RW', 'AR', 'MA'), bty='n', pt.cex=1.5, cex=0.85)

dev.off()
#####

#3) LAR/MSAE####
#Calculate the mean for each approach across all populations
exp(median(as.vector(unlist(coho_forecast_RW$MAPE)), na.rm=TRUE)) -1
exp(median(as.vector(unlist(coho_forecast_RW$RW_MAPE)), na.rm=TRUE)) -1
exp(median(as.vector(unlist(coho_forecast_RW$RW_fit_MAPE)), na.rm=TRUE)) -1
exp(median(as.vector(unlist(coho_forecast_AR$RW_fit_MAPE)), na.rm=TRUE)) -1 
exp(median(as.vector(unlist(coho_forecast_MA$RW_fit_MAPE)), na.rm=TRUE)) -1

#Calculate forecast accuracy only for the FRAM units included in the MUs of the Published coho forecast
#Includes the Published method
n_range <- 16
n_pop <- 36
MAPE_mat_ST <- matrix(ncol=n_pop, nrow=n_range)
MAPE_mat_RW <- matrix(ncol=n_pop, nrow=n_range)
MAPE_mat_RW_fit <- matrix(ncol=n_pop, nrow=n_range)
MAPE_mat_AR <- matrix(ncol=n_pop, nrow=n_range)
MAPE_mat_MA <- matrix(ncol=n_pop, nrow=n_range)

for(i in 1:n_pop){
  MAPE_mat_ST[,i] <- unlist(lapply(coho_forecast_RW$MAPE, function(x) x[i]))
  MAPE_mat_RW[,i] <- unlist(lapply(coho_forecast_RW$RW_MAPE, function(x) x[i]))
  MAPE_mat_RW_fit[,i] <- unlist(lapply(coho_forecast_RW$RW_fit_MAPE, function(x) x[i]))
  MAPE_mat_AR[,i] <- unlist(lapply(coho_forecast_AR$RW_fit_MAPE, function(x) x[i]))
  MAPE_mat_MA[,i] <- unlist(lapply(coho_forecast_MA$RW_fit_MAPE, function(x) x[i]))
}

exp(median(MAPE_mat_RW_fit[,cast_pop], na.rm=TRUE)) -1 
exp(median(MAPE_mat_AR[,cast_pop], na.rm=TRUE)) -1 
exp(median(MAPE_mat_MA[,cast_pop], na.rm=TRUE)) -1 
exp(median(MAPE_mat_ST[,cast_pop], na.rm=TRUE)) -1 

#For the Published method
state_cast[state_cast==0] <- 1
MAPE_WDFW <- abs(log(state_cast[,-c(1,18,19,20)]/matrix(unlist(coho_forecast_RW$obs), nrow=n_pop, ncol=n_range, byrow=FALSE)[cast_pop,]))
exp(median(unlist(MAPE_WDFW), na.rm=TRUE))-1

#Create function for calculating LAR from output
lar_fun <- function(x){
  return (exp(median(x, na.rm=TRUE))-1)
}

#Plot relative performence of each method according to MAPE
pdf(file='forecast_summary_plot_MSA.pdf')
col.vec <- brewer.pal(n=9, 'YlOrBr')
col.vec_2 <- col2rgb(col.vec)
par(oma=c(8,12,10,1), mar=c(5,1,0,0), mfrow=c(1,5))
barplot(apply(MAPE_WDFW, 1, lar_fun), horiz=TRUE, las=2, col='black', names.arg=unique(dat_list$coho_dat_full$Name)[cast_pop], cex.names=1, border=NA, xlim=c(0,7))
barplot(apply(MAPE_mat_ST[,cast_pop], 2, lar_fun), add=TRUE, horiz=TRUE, las=2,col=rgb(17,156,240, alpha=150, max=255 ), border=NA)
text(5, 1, 'A')
mtext(side=3, 'ST-IPM', line=1.25, cex=0.85)
mtext(side=3, 'vs.', line=0.25, cex=0.85)
mtext(side=3, 'Published', line=-0.75, cex=0.85)


barplot(apply(MAPE_mat_RW[,cast_pop], 2, lar_fun),col=rgb(red=col.vec_2[1,9],green=col.vec_2[2,9],blue=col.vec_2[3,9], alpha=200, max=255 ), horiz=TRUE, las=2, cex.names=0.7, border=NA, xlim=c(0,7))
barplot(apply(MAPE_mat_ST[,cast_pop], 2, lar_fun), add=TRUE, horiz=TRUE, las=2,col=rgb(17,156,240, alpha=150, max=255 ), border=NA)
mtext(side=1, 'MSA', line=2.5, cex=0.85, at=7.5)
text(5, 1, 'B')
mtext(side=3, 'ST-IPM', line=1.25, cex=0.85)
mtext(side=3, 'vs.', line=0.25, cex=0.85)
mtext(side=3, 'Random walk', line=-0.75, cex=0.85)

#barplot(apply(MAPE_mat_RW_fit[,cast_pop], 2, mean, na.rm=TRUE), horiz=TRUE, las=2, col=rgb(red=col.vec_2[1,8],green=col.vec_2[2,8],blue=col.vec_2[3,8], alpha=200, max=255 ), border=NA)
#barplot(apply(MAPE_mat_ST[,cast_pop], 2, mean, na.rm=TRUE), add=TRUE, horiz=TRUE, las=2,col=rgb(17,156,240, alpha=150, max=255 ), border=NA)
#mtext(side=1, 'MMAPE', line=2.5, cex=0.85)

barplot(apply(MAPE_mat_AR[,cast_pop], 2, lar_fun), horiz=TRUE, las=2, col=rgb(red=col.vec_2[1,6],green=col.vec_2[2,6],blue=col.vec_2[3,6], alpha=200, max=255 ), border=NA, xlim=c(0,7))
barplot(apply(MAPE_mat_ST[,cast_pop], 2, lar_fun), add=TRUE, horiz=TRUE, las=2,col=rgb(17,156,240, alpha=150, max=255 ), border=NA)
text(5, 1, 'C')
mtext(side=3, 'ST-IPM', line=1.25, cex=0.85)
mtext(side=3, 'vs.', line=0.25, cex=0.85)
mtext(side=3, 'Autoregressive', line=-0.75, cex=0.85)

barplot(apply(MAPE_mat_MA[,cast_pop], 2, lar_fun), horiz=TRUE, las=2, col=rgb(red=col.vec_2[1,4],green=col.vec_2[2,4],blue=col.vec_2[3,4], alpha=200, max=255 ), border=NA, xlim=c(0,7))
barplot(apply(MAPE_mat_ST[,cast_pop], 2, lar_fun), add=TRUE, horiz=TRUE, las=2,col=rgb(17,156,240, alpha=150, max=255 ), border=NA)
text(5, 1, 'D')
mtext(side=3, 'ST-IPM', line=1.25, cex=0.85)
mtext(side=3, 'vs.', line=0.25, cex=0.85)
mtext(side=3, 'Published', line=-0.75, cex=0.85)

legend(x=2.5, y=34, pch=15, col=c(col=rgb(17,156,240, alpha=150, max=255 ), 'black', 
                                  rgb(red=col.vec_2[1,8],green=col.vec_2[2,8],blue=col.vec_2[3,8], alpha=200, max=255),
                                  rgb(red=col.vec_2[1,6],green=col.vec_2[2,6],blue=col.vec_2[3,6], alpha=200, max=255),
                                  rgb(red=col.vec_2[1,4],green=col.vec_2[2,4],blue=col.vec_2[3,4], alpha=200, max=255)),
       legend=c('ST-IPM', 'Published', 'RW', 'AR', 'MA'), bty='n', pt.cex=1.5, cex=0.85)

dev.off()
#####

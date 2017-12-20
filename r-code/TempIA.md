Prediction of flight onset of *Cacopsylla melanoneura* and *C. picta*, vectors of apple proliferation disease, using a temperature-based immigration analysis - Appendix A
==========================================================================================================================================================================

Introduction
------------

This document is an appendix to the article: ''Prediction of flight onset of *Cacopsylla melanoneura* and *C. picta*, vectors of apple proliferation disease, using a temperature-based immigration analysis''. The document provides the calculation procedure for a temperature-based immigation analysis based on Tedeschi et al. (2012).

Data
----

##### Loading data and libraries

The dataset "TempIA.RData" is not freely available but access can be requested following [this link](https://berndpanassiti.github.io).

TempIA.RData includes the following R-objects:

-   c.melanoneuraPB

-   c.melanoneuraF1B

-   c.pictaPB

-   c.pictaF1B

-   c.melanoneuraPYT

-   c.melanoneuraF1YT

-   c.pictaPYT

-   c.pictaF1YT

-   modeldata.t

-   sitesInsectVector

-   weatherStationsHourly\_split

-   clusterStations

``` r
rm(list=ls(all=TRUE))
Sys.setenv(TZ="Europe/Rome")
#Sys.getenv("TZ")

source("r-code/00_settings.r")


load(file="data/TempIA.RData")


detachAllPackages() #detach("package:plyr", unload=TRUE) # makes problems with the group_by function
if(!require(ggplot2))
{
  install.packages("ggplot2", repos = mir)
  require(ggplot2)
}
if(!require(plyr))
{
  install.packages("plyr")
  require(plyr)
}
if(!require(dplyr))
{
  install.packages("dplyr")
  require(dplyr)
}

if(!require(stringi))
{
  install.packages("stringi", repos = mir)
  require(stringi)
}

if(!require(scales))
{
  install.packages("scales", repos = mir)
  require(scales)
}
if(!require(MASS))
{
  install.packages("MASS", repos = mir)
  require(MASS)
}
if(!require(qtlcharts))
{
  install.packages("qtlcharts", repos = mir)
  require(qtlcharts)
}

if(!require(lubridate)) # date functions
{
  install.packages("lubridates", repos = mir)
  require(lubridate)
}
```

#### Data preparation

P = parents, re-emigrants, overwintering vectors
F1 = offspring, emigrants
B = vector sampling method: beating tray
YT = vector sampling method: yellow traps

``` r
tempsum <- data.frame()

#pulling out site_id,cluster, year, date and abundance for each insect_id and insect_stage_id
c.melPB <- c.melanoneuraPB[,c(2,3,4,5,6,10,11)]
c.melF1B <- c.melanoneuraF1B[,c(2,3,4,5,6,10,11)]
c.picPB <- c.pictaPB[,c(2,3,4,5,6,10,11)]
c.picF1B <- c.pictaF1B[,c(2,3,4,5,6,10,11)]
c.melPYT <- c.melanoneuraPYT[,c(2,3,4,5,6,10,11)]
c.melF1YT <- c.melanoneuraF1YT[,c(2,3,4,5,6,10,11)]
c.picPYT <- c.pictaPYT[,c(2,3,4,5,6,10,11)]
c.picF1YT <- c.pictaF1YT[,c(2,3,4,5,6,10,11)]

#renaming 'amount'
names(c.melPB)[6] <- 'abundance'
names(c.melF1B)[6] <- 'abundance'
names(c.picPB)[6] <- 'abundance'
names(c.picF1B)[6] <- 'abundance'
names(c.melPYT)[6] <- 'abundance'
names(c.melF1YT)[6] <- 'abundance'
names(c.picPYT)[6] <- 'abundance'
names(c.picF1YT)[6] <- 'abundance'


clusterSites <- unique(modeldata.t[,c(2,3)]) #site_id and cluster id
regionSites <- unique(sitesInsectVector[,c(1,3)])

#Combining beating and yellow trap subsets
c.melPcomb <- merge(c.melPB, c.melPYT, by = c("site_id","region","date","cluster","year"), all = TRUE)
c.melPcomb <- transform(c.melPcomb, abundance = rowSums(c.melPcomb[,c("abundance.x","abundance.y")], na.rm = TRUE))
c.melPcomb[, "presence"] <- apply(c.melPcomb[, c("presence.x","presence.y")], 1, max, na.rm=T)
c.melPcomb <- c.melPcomb[,-which(names(c.melPcomb) %in% c("abundance.x","abundance.y","presence.x","presence.y"))]#dropping uncombined columns
###
c.melF1comb <- merge(c.melF1B, c.melF1YT, by = c("site_id","region","date","cluster","year"), all = TRUE)
c.melF1comb <- transform(c.melF1comb, abundance = rowSums(c.melF1comb[,c("abundance.x","abundance.y")], na.rm = TRUE))
c.melF1comb[, "presence"] <- apply(c.melF1comb[, c("presence.x","presence.y")], 1, max, na.rm=T)
c.melF1comb <- c.melF1comb[,-which(names(c.melF1comb) %in% c("abundance.x","abundance.y","presence.x","presence.y"))]#dropping uncombined columns
###
c.picPcomb <- merge(c.picPB, c.picPYT, by = c("site_id","region","date","cluster","year"), all = TRUE)
c.picPcomb <- transform(c.picPcomb, abundance = rowSums(c.picPcomb[,c("abundance.x","abundance.y")], na.rm = TRUE))
c.picPcomb[, "presence"] <- apply(c.picPcomb[, c("presence.x","presence.y")], 1, max, na.rm=T)
c.picPcomb <- c.picPcomb[,-which(names(c.picPcomb) %in% c("abundance.x","abundance.y","presence.x","presence.y"))]#dropping uncombined columns
###
c.picF1comb <- merge(c.picF1B, c.picF1YT, by = c("site_id","region","date","cluster","year"), all = TRUE)
c.picF1comb <- transform(c.picF1comb, abundance = rowSums(c.picF1comb[,c("abundance.x","abundance.y")], na.rm = TRUE))
c.picF1comb[, "presence"] <- apply(c.picF1comb[, c("presence.x","presence.y")], 1, max, na.rm=T)
c.picF1comb <- c.picF1comb[,-which(names(c.picF1comb) %in% c("abundance.x","abundance.y","presence.x","presence.y"))]#dropping uncombined columns
###

#TEMPSUMLIST CREATION
#listnames <- c('c.melP','c.melF1','c.picP','c.picF1','c.melPYT','c.melF1YT','c.picPYT','c.picF1YT')
#tempsumlist <- list(c.melP,c.melF1,c.picP,c.picF1,c.melPYT,c.melF1YT,c.picPYT,c.picF1YT)
listnames <- c('c.melPcomb','c.melF1comb','c.picPcomb','c.picF1comb')
tempsumlist <- list(c.melPcomb,c.melF1comb,c.picPcomb,c.picF1comb)

names(tempsumlist) <- listnames


############### DECIDE IF SUM BY REGION OR CLUSTER HERE! #########################
#Sum by region
for(i in 1:length(tempsumlist)){
tempsumlist[[i]] <- as.data.frame(tempsumlist[[i]] %>% dplyr::group_by(date,year,region) %>% dplyr::summarise(abundance = sum(abundance,na.rm=T), presence=max(presence))) #combine all sites within same region at same date
}

list2env(tempsumlist,environment())  #overwrite single dataframes in the global environment

#merge all dataframes by columns
tempsum <- join_all(list(tempsum, c.melPcomb,c.melF1comb,c.picPcomb,c.picF1comb), type="full")
```

``` r
# This chunk filters dataset to get only species/region/year combinations with absences before presences!

df <- tempsumlist

SpeciesRegionYear_Absences_before_Presences <- names(which(unlist(
lapply(df, 
    function(j){
                                  sapply(split(j, list(j[,3],j[,2])), function(x){x[1,5]==0})
  }
))==TRUE))


NlistElements <- sapply(df, NROW)


idSpeciesRegionYear=data.frame(rep(names(df),sapply(df,NROW)),paste(rep(names(df),sapply(df,NROW)),unlist(lapply(df, "[", c(3))),unlist(lapply(df, "[", c(2))),sep="."))
colnames(idSpeciesRegionYear) <- c("species","id")

dfNEW <- df

dfNEW$c.melPcomb[["id"]] <- idSpeciesRegionYear[which(idSpeciesRegionYear$species=="c.melPcomb"),2]   #mel remigrants
dfNEW$c.melF1comb[["id"]] <- idSpeciesRegionYear[which(idSpeciesRegionYear$species=="c.melF1comb"),2] #mel emigrants

dfNEW$c.picPcomb[["id"]] <- idSpeciesRegionYear[which(idSpeciesRegionYear$species=="c.picPcomb"),2]   #pic remigrants
dfNEW$c.picF1comb[["id"]] <- idSpeciesRegionYear[which(idSpeciesRegionYear$species=="c.picF1comb"),2] #pic emigrants




# FINAL step - FILTER for absences before presences

dfNEW2 <- lapply(dfNEW, function(j){
                                  j[which(j[,6]%in%SpeciesRegionYear_Absences_before_Presences),]
}
)





# select only records with at least n years of PRESENCES
n=1


SpeciesRegionYear_Presences <- names(which(unlist(
lapply(dfNEW2, 
    function(j){
                                  sapply(split(j, list(j[,3],j[,2])), function(x){sum(x[,5])>0})
  }
))==TRUE))


NrecordsYear <- table(substr(SpeciesRegionYear_Presences,1,nchar(SpeciesRegionYear_Presences)-5))
#NrecordsYear[NrecordsYear>1]

NrecordsYearThreshold <- NrecordsYear[NrecordsYear>=n] # minimum 3 years of records per species and regions



dfNEW3 <- lapply(dfNEW2, function(j){
                                  j[which(gsub("(.*)\\..*","\\1", j[,6])  %in% names(NrecordsYearThreshold)),]
}
) 

tempsumlist <- NULL
tempsumlist <- dfNEW3
```

Temperature-based immigration analysis
--------------------------------------

##### Calculation of thresholds and indices

The code below calculates the following thresholds and indices:

-   **a0** = Date of first occurence for each site and year (+ species, generation and method)

-   aMax = Date of max. abundance for each site (+ species, generation and method)

-   thourly = Temperature for each hour for each weather station

-   thourlyregion = Temperature for each hour for each weather station summarized by region

-   TMax = Max. daily temperature & tmedian from hourly mean temps

-   **T0Max** = max of max daily temperature of 7 days preceding a0 (per observation, insect type and site)

-   T0maxmin and T0maxmean: Minimum of all t0max/t0mean that belong to the same site over several years

-   T0maxyears = Min of t0maxmin by year for id (region || cluster) in order to check for regional differences

-   DD (degree days) = number of hours over t0max per week preceding a0

-   T7n = mean max daily temperatures in 7 days before a0 (using Tmax (max daily temperatures))

-   **Ii** = Immigration Index for the a0; Ii = (T7n - T7th) + DD

1. a0 and aMax

``` r
##################################################################################################
#CALCULATING a0
a0 <-lapply(tempsumlist, function(j){
  sapply(split(j, list(j[,3],j[,2])), function(x){x[,1][which(x[,5]>0)[1]]}) #split by region||cluster||site + year and look for first date where insects are found
})
a0 <- lapply(a0, function(x) x[!is.na(x)]) #removing sites and year where no insect was found thus returning NA
for(i in 1:length(a0)){
  a0[[i]] <- as.data.frame(a0[[i]])
  a0[[i]]$datetime <- as.POSIXct(a0[[i]][,1], origin="1970-01-01")
  a0[[i]]$year <- format(a0[[i]]$datetime, format="%Y")
  #a0[[i]]$t0maxmin <- NA
  #a0[[i]]$t0maxmean <- NA
}

for(i in 1:length(a0)){ #Adding region||cluster||site_id to a0 as id
  for(j in 1:length(a0[[i]][,1])){
    a0[[i]]$region <- rownames(a0[[i]][1])
    a0[[i]]$region <- gsub('.{5}$', '', a0[[i]]$region) #cut off last 5 characters (e.g. Laimburg0001.2016 -> Laimburg0001 || Bozen.2013 <- Bozen for easier filtering of overall site min t0max)
  }
}

##################################################################################################
#CALCULATING aMax
aMax <-lapply(tempsumlist, function(j){
  sapply(split(j, list(j[,3])), function(x){x[,1][which.max(x[,5])]}) #split by region||cluster||site and look for first date where max # of insects are found
})
aMax <- lapply(aMax, function(x) x[!is.na(x)])
for(i in 1:length(aMax)){
  aMax[[i]] <- as.data.frame(aMax[[i]])
  aMax[[i]]$datetime <- as.POSIXct(aMax[[i]][,1], origin="1970-01-01")
  aMax[[i]]$year <- format(aMax[[i]]$datetime, format="%Y")
}

for(i in 1:length(aMax)){ #Adding region||cluster||site_id to aMax as id
  for(j in 1:length(aMax[[i]][,1])){
    aMax[[i]]$region <- rownames(aMax[[i]][1])
  }
}
```

``` r
##################################################################################################
thourly <- weatherStationsHourly_split
thourly <- lapply(thourly, function(x) x[,-c(6:11)]) #removing unneeded columns
for(i in 1:nrow(clusterStations)){
  thourly[[i]]$cluster <- clusterStations$group[i]
}

for(i in 1:length(thourly)){
  thourly[[i]]$datetimenum <- as.POSIXct(thourly[[i]]$datetime, format="%Y-%m-%d")#cut off %H:%M:%S due to problems with numeric transformation
  thourly[[i]]$datetimenum <- as.numeric(thourly[[i]]$datetimenum)
}

#################################################
thourlyregion <- dplyr::bind_rows(thourly) #combine lists to 1 dataframe for dplyr
thourlyregion <- as.data.frame(thourlyregion %>% dplyr::group_by(datetime,datetimenum,region) %>% dplyr::summarise(temphourmean=mean(temphour, na.rm=T)))
thourlyregion <- split(thourlyregion, thourlyregion$region) #combine dataframes back into list
```

2. Tmax and Tmedian

``` r
##################################################################################################
#CALCULATING TMax = Max. daily temperature & tmedian from hourly mean temps
tmax <- list()
tmedian <- list()


for(i in 1:length(weatherStationsHourly_split)){
tmax[[i]] <- aggregate(list(maxtempday=weatherStationsHourly_split[[i]][,4]),list(datetime=cut(as.POSIXct(weatherStationsHourly_split[[i]][,3]), "day")),FUN=max,na.rm=T)
tmax[[i]]$maxtempday[is.infinite(tmax[[i]]$maxtempday)] <- NA
tmax[[i]]$datetime[is.infinite(tmax[[i]]$datetime)] <- NA
tmax[[i]]$site_id <- names(weatherStationsHourly_split)[[i]]
tmax[[i]]$region  <- unique(weatherStationsHourly_split[[i]][1,5]) #repeating the stations region for easier filtering

tmedian[[i]] <- aggregate(list(mediantempday=weatherStationsHourly_split[[i]][,4]),list(datetime=cut(as.POSIXct(weatherStationsHourly_split[[i]][,3]), "day")),FUN=median,na.rm=T)
tmedian[[i]]$mediantempday[is.infinite(tmedian[[i]]$mediantempday)] <- NA
tmedian[[i]]$datetime[is.infinite(tmedian[[i]]$datetime)] <- NA
tmedian[[i]]$site_id <- names(weatherStationsHourly_split)[[i]]
tmedian[[i]]$region  <- unique(weatherStationsHourly_split[[i]][1,5]) #repeating the stations region for easier filtering
}


#naming the stations
names(tmax)    <- names(weatherStationsHourly_split)
names(tmedian) <- names(weatherStationsHourly_split)

#changing the datetime entry to numeric for easier filtering
for(i in 1:length(tmax)){
  tmax[[i]]$datetime <- as.POSIXct(tmax[[i]]$datetime, format="%Y-%m-%d")
  tmax[[i]]$datetime <- as.numeric(tmax[[i]]$datetime)
}

for(i in 1:length(tmedian)){
  tmedian[[i]]$datetime <- as.POSIXct(tmedian[[i]]$datetime, format="%Y-%m-%d")
  tmedian[[i]]$datetime <- as.numeric(tmedian[[i]]$datetime)
}


##################################################################################################
#adding cluster id & region to tmax and tmedian
for(i in 1:nrow(clusterStations)){
  tmax[[i]]$cluster <- clusterStations$group[i]
  tmedian[[i]]$cluster <- clusterStations$group[i]
}
```

3. T0max

``` r
##################################################################################################
#T0Max = max of max daily temperature of 7 days preceding a0 (per observation, insect type and site)
#check a0 and then filter through txmax and tmedian at corresponding weatherstation (from cluster)
#using tmax
for(i in 1:length(a0)){ #loop through all species,generation and method
  for(j in 1:length(a0[[i]][,1])){ #loop through the a0 entries
    tmaxvec <- numeric()
    for(x in 1:length(tmax)){ #loop through ALL weather stations
      if(tmax[[x]]$region[1] == a0[[i]]$region[j] ){ #check which stations cluster id matches the current sites' cluster id
  tmaxvec  <- c(tmaxvec,max(tmax[[x]]$maxtempday[(which(tmax[[x]]$datetime == a0[[i]][,1][j])-7):(which(tmax[[x]]$datetime == a0[[i]][,1][j]))])) #temporary vector saving the max of daily max temperatures of the 7 days before a0 from all weather stations x that share the id with the site j
       print(paste(i,j,x))
      }
    }
    print(tmaxvec)
    a0[[i]]$t0max[j] <- min(tmaxvec, na.rm=T) #min of all weather station's 7 day max
    a0[[i]]$t0mean[j] <- mean(tmaxvec, na.rm=T) #mean of all max weekly temperatures
  }
 a0[[i]]$t0max[is.infinite(a0[[i]]$t0max)] <- NA
 a0[[i]]$t0mean[is.nan(a0[[i]]$t0mean)] <- NA
}
```

4. T0maxmin and T0maxmean

``` r
##################################################################################################
#Minimum of all t0max/t0mean that belong to the same site over several years
t0maxmin <- list()
for(i in 1:length(a0)){
t0maxmin[[i]] <- as.data.frame(a0[[i]] %>% dplyr::group_by(region) %>% dplyr::summarise(t0maxmin = min(t0max,na.rm=T), t0meanmin = min(t0mean,na.rm=T)))
}
names(t0maxmin) <- names(a0)

for(i in 1:length(a0)){
  a0[[i]] <- join_all(list(a0[[i]], t0maxmin[[i]]), by="region", type="left")
}
```

5. T0maxyears

``` r
##################################################################################################
#Min of t0maxmin by year for id (region || cluster) in order to check for regional differences
t0maxyears <- list()
for(i in 1:length(a0)){
t0maxyears[[i]] <- as.data.frame(a0[[i]] %>% dplyr::group_by(region) %>% dplyr::summarise(t0maxyears= min(t0maxmin,na.rm=T)))
}
names(t0maxyears) <- names(a0)

for(i in 1:length(a0)){
  a0[[i]] <- join_all(list(a0[[i]], t0maxyears[[i]]), by="region", type="left")
}
```

6. DD

``` r
##################################################################################################
#Degree days DD = medium daily number of hours over t0max per week preceding a0
for(i in 1:length(a0)){ # 1:4 - Cmel P and F1 and Cpicta P and F1
  for(j in 1:length(a0[[i]][,1])){ # number of total records
    DD <- numeric()
      for(x in 1:length(thourly)){ # x select the different stations (belonging to different regions) in thourly
        if(thourly[[x]]$region[1] == a0[[i]]$region[j]){ # region weather station = region sampling location
          
          df <- thourly[[x]][
                              (which(thourly[[x]]$datetimenum == a0[[i]][,1][j])-(24*7)):
                              (which(thourly[[x]]$datetimenum == a0[[i]][,1][j]))
                            ,c("datetime","temphour")]
          
          df <- df[df$temphour > a0[[i]]$t0maxmin[j],] # select hours during the 7 days before a0 where temp > t0max
          DD <-  c(DD,nrow(df)/7) # number of hours / 7 days
                  
          print(paste(c(i,j,x)))
       }
      }
    print(DD)
    a0[[i]]$DD[j] <- mean(DD, na.rm=T)
  }
}

# *** Explaination loop *** #
# - select vector, e.g. C. melanoneura parental generation (i)
# - loop through survey records (j)
# - select weather station in the same region (x)
```

7. T7n

``` r
############################################################################################
# mean max daily temperatures in 7 days before a0 (using Tmax (max daily temperatures))
for(i in 1:length(a0)){
  for(j in 1:length(a0[[i]][,1])){
    T7n <- numeric()
      for(x in 1:length(tmax)){
       #if(tmax[[x]]$cluster[1] == a0[[i]]$cluster[j]){
        if(thourly[[x]]$region[1] == a0[[i]]$region[j]){
          T7n <-  c(T7n,tmax[[x]]$maxtempday[(which(tmax[[x]]$datetime == a0[[i]][,1][j])-(7)):(which(tmax[[x]]$datetime == a0[[i]][,1][j]))]) #(HINT: 7 days since daily temp data)
       }
      }
    print(T7n)
    a0[[i]]$T7n[j] <- mean(T7n, na.rm=T)
  }
}
```

8. Ii

``` r
##################################################################################################
#Immigration Index for the a0 
#Ii = (T7n - T7th) + DD
#   = mean of max in 7 days before a0 (T7n) - absolute highest temp in week before a0 (t0max) + hours per week above t0max (DD)
for(i in 1:length(a0)){
  for(j in 1:length(a0[[i]][,1])){
    a0[[i]]$Ii[j] <- (a0[[i]]$T7n[j] - a0[[i]]$t0max[j]) + a0[[i]]$DD[j]
  }
}
```

9. Only select regions with at least 2 years of survey data

``` r
a0 <- selectAtleast2Regions(a0) # automatically backups a0 to a0_old
```

10. Species specific selection of temperature and immigration data

``` r
#Index of tempsumlist of dataset you want to look at 
# c.melPcomb:  x=1
# c.melF1comb: x=2
# c.picPcomb:  x=3
# c.picF1comb: x=4

species <- data.frame(species =c("c.melPcomb","c.melF1comb","c.picPcomb","c.picF1comb"),id=1:4)
Species <- species[which(species[,1] %in% colnames(sapply(a0, names))),]$id

Regions <- data.frame(regions =c("Bozen","Burggrafenamt","Eisacktal","Salten-Schlern","Ueberetsch-Unterland","Vinschgau"),id=1:6)


thourlist <- list()

for (s in 1:length(Species)){
  
  
##############################################################################################
  
        x <-Species[s] # for which species is the run? ENTERING DESIRED TEMPSUMLIST INDEX
                                      
##############################################################################################

                                          
r <- Regions[which(Regions$regions %in% unique(a0[[as.character(species[Species[s], 1])]]$region)),2] # ids of regions for the graph



#Calculating max temp per day (tmax) and hours over t0max per day (DD)
for(i in 1:length(thourly)){
  thourly[[i]]$DD <- NA
  thourly[[i]]$tmax <- NA
  for(j in 24:length(thourly[[i]][,3])){
    for(x in x){  #just for one dataset (=x) of tempsumlist
      for(y in 1:length(t0maxyears[[x]][,1])){
      if(thourly[[i]]$region[j] == t0maxyears[[x]]$region[y]){
      thourly[[i]]$DD[(j-24):j] <- sum(as.numeric(thourly[[i]]$temphour[(j-(24)):j])>t0maxyears[[x]]$t0maxyears[y])
      thourly[[i]]$tmax[(j-24):j] <- max(thourly[[i]]$temphour[(j-(24)):j])
        }
      }
    }
  }
}





#Ii aggregated for each id (taking means of all stations within the same id (region || cluster))
thourlyregion <- dplyr::bind_rows(thourly) #combine lists to 1 dataframe for dplyr
thourlyregion <- as.data.frame(thourlyregion %>% dplyr::group_by(datetime,datetimenum,region) %>% dplyr::summarise(temphourmean=mean(temphour, na.rm=T), ddmean=mean(DD,na.rm=T), tmaxmean=mean(tmax, na.rm=T)))

names(t0maxyears[[x]]) <- c("region", "t0maxyears") #rename for use in join_all
thourlyregion <- join_all(list(thourlyregion, t0maxyears[[x]]), by="region", type="left")
thourlyregion$Ii <- NA

thourlyregion <- split(thourlyregion, thourlyregion$region) #combine dataframes back into list


for(i in 1:length(thourlyregion)){
  for(j in 168:length(thourlyregion[[i]][,1])){
  thourlyregion[[i]]$Ii[(j-(24*7)):j] <- (mean(thourlyregion[[i]]$tmaxmean[(j-(24*7)):j]) - thourlyregion[[i]]$t0maxyears[j]) + mean((thourlyregion[[i]]$ddmean[(j-(24*7)):j])>thourlyregion[[i]]$t0maxyears[j])
  }
}





#tmax, DD, Ii region
thourlyregion <- dplyr::bind_rows(thourlyregion) #combine lists to 1 dataframe for dplyr
names(t0maxyears[[x]]) <- c("region", "t0maxyears") #rename for use in join_all
thourlyregion <- join_all(list(thourlyregion, t0maxyears[[x]]), by="region", type="left")
thourlyregion <- split(thourlyregion, thourlyregion$region) #combine dataframes back into list

for(i in 1:length(thourlyregion)){
  thourlyregion[[i]]$DD <- NA
  thourlyregion[[i]]$tmax <- NA
  for(j in 24:length(thourlyregion[[i]][,1])){
    thourlyregion[[i]]$DD[(j-24):j] <- sum(as.numeric(thourlyregion[[i]]$temphourmean[(j-(24)):j])>thourlyregion[[i]]$t0maxyears[j])
      thourlyregion[[i]]$tmax[(j-24):j] <- max(thourlyregion[[i]]$temphourmean[(j-(24)):j])
  }
}

for(i in 1:length(thourlyregion)){
  thourlyregion[[i]]$Ii <- NA
  for(j in 168:length(thourlyregion[[i]][,1])){
  thourlyregion[[i]]$Ii[(j-(24*7)):j] <- (mean(thourlyregion[[i]]$tmax[(j-(24*7)):j]) - thourlyregion[[i]]$t0maxyears[j]) + mean((thourlyregion[[i]]$DD[(j-(24*7)):j])>thourlyregion[[i]]$t0maxyears[j])
  }
}


tmp <- list(a=thourly,b=thourlyregion)
thourlist[[s]] <- tmp
}
```

11. FINAL: Species and region-specific graphics of temperature and immigration index thresholds

``` r
Sys.setlocale("LC_TIME", "English") # to match English date structure
```

    ## [1] ""

``` r
graphtitleSpecies <- c("Cacopsylla melanoneura", "Cacopsylla melanoneura ","Cacopsylla picta","Cacopsylla picta")
graphtitleType <- c("remigrants", "emigrants","remigrants","emigrants")




par(oma=c(5.5,3,2,2),mfrow=c(4,1),mar=c(5,5,4,5)+.1)

for (s in 1:length(Species)){

thourly       <- thourlist[[s]]$a
thourlyregion <- thourlist[[s]]$b
  
##############################################################################################
  
        x <-Species[s] # for which species is the run?
                                      
##############################################################################################

                                          
r <- Regions[which(Regions$regions %in% unique(a0[[as.character(species[Species[s], 1])]]$region)),2] # ids of regions for the graph

mainTitle <- list(as.expression(bquote(~italic(.(graphtitleSpecies[x]))~~.(graphtitleType[x])),cex=cexmlabel2))

for(j in x){ # insect type
for(rr in 1:length(r)){ # region
  i <- r[rr]

  
  
plot(thourlyregion[[i]]$datetime, thourlyregion[[i]]$Ii, main=names(thourlyregion[i]), type="l", ylab="Immigration Index (Ii)",xlab="",
     ylim=c(-15,40),col="darkgrey",xaxt="n",cex.lab=2,cex.main=2, cex.axis=2)




abline(h=0, col="darkgrey") # Immigration index


for(y in c("2013","2014","2015","2016")){
tmp <- subset(thourlyregion[[i]], format(datetime, "%Y") %in% c(y))
abline(v=as.POSIXct(tmp$datetime[which(tmp$Ii > 0)[1]],format="%Y-%m-%d %H:%M%S"), col="darkgrey",lty=2,lwd=4)
abline(v=as.POSIXct(tmp$datetime[which(tmp$temphourmean >= thourlyregion[[i]]$t0maxyears[1])[1]],format="%Y-%m-%d %H:%M%S"), col=adjustcolor("red",0.3),lty=2,lwd=4)
}



points(subset(tempsumlist[[j]], region==names(thourlyregion)[i])$date, subset(tempsumlist[[j]], region==names(thourlyregion)[i])$presence * 40,col = adjustcolor(ifelse(subset(tempsumlist[[j]], region==names(thourlyregion)[i])$presence > 0,'blue','red'),0.5), pch=16,xaxt="n",cex=2)
text(as.POSIXct("2016-06-01"), thourlyregion[[i]]$t0maxyears[1]+5, "T7th", col = "red",cex=3) 
text(as.POSIXct("2013-02-01"), 0+5, "Ii=0", col = "black",cex=3) 
text(as.POSIXct("2014-11-01"), 35, paste("T7th= ",round(thourlyregion[[i]]$t0maxyears[1],digits=2)," °C"),cex=2)


####x axis###
minor <- seq(as.POSIXct("2013-01-01", format = "%Y-%m-%d"), tail(thourlyregion[[i]]$datetime, 1),
             by = "1 months")
labDates <- seq(as.POSIXct("2013-01-15", format = "%Y-%m-%d"),tail(thourlyregion[[i]]$datetime,1), by = "2 months")

yearPos <- c("2013-08-01 CEST","2014-08-01 CET","2015-08-01 CEST","2016-08-01 CEST")

axis.POSIXct(side = 1, thourlyregion[[i]]$datetime, at = labDates, format = "%b", las = 2,cex.axis=2,tick = FALSE)
axis.POSIXct(side = 1, thourlyregion[[i]]$datetime, at = minor, labels = FALSE, tcl = -0.25)
axis.POSIXct(side = 1, thourlyregion[[i]]$datetime, at = yearPos,format = "%Y",las=1, line=4,lty=0,cex.axis=2)


###y temperature###
lines(thourlyregion[[i]]$datetime, thourlyregion[[i]]$temphourmean,yaxt="n",xaxt="n", ylab="",xlab="",col=adjustcolor("red",0.3), type="l", ylim=c(-15,40))
axis(4, cex.axis=2)
abline(h=round(thourlyregion[[i]]$t0maxyears[1],digits=2), col="red")
mtext("Mean Hourly Temperature °C",side=4,line=3,cex=1)
}
#title(mainTitle, outer=T,cex.main=2)
#title(names(tempsumlist)[j], outer=T,cex.main=2)
}

}

##legend###
par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0), mar = c(0, 0, 0, 0), new = TRUE)
plot(0, 0, type = "n", bty = "n", xaxt = "n", yaxt = "n")
legend("bottom", c("Ii","T[°C]","1st Ii>0","1st T>T7th","Presence", "Absence"), col=c("darkgrey","red","darkgrey","red","blue","red"),lty=c(1,1,2,2,NA,NA),pch=c(NA,NA,NA,NA,16,16),bty="n",cex=1.75,horiz = T,pt.cex=2,text.width = 0.21, seg.len=0.75,lwd=2,xpd=2)
```

<img src="TempIA_files/figure-markdown_github/base plots-1.png" style="display: block; margin: auto;" />

References
==========

Tedeschi R., Baldessari M., Mazzoni V., Trona F. and Angeli G. (2012). Population dynamics of *Cacopsylla melanoneura* (Hemiptera: Psyllidae) in northeast Italy and its role in the apple proliferation epidemiology in apple orchards. - Journal of Economic Entomology 105, 322–328.

Session info
------------

| setting  | value                        |
|----------|------------------------------|
| version  | R version 3.3.1 (2016-06-21) |
| system   | x86\_64, darwin13.4.0        |
| ui       | RStudio (1.0.143)            |
| language | (EN)                         |
| collate  | en\_US.UTF-8                 |
| tz       | Europe/Rome                  |

| package    | \*  | version | date       | source         |
|------------|-----|---------|------------|----------------|
| assertthat |     | 0.1     | 2013-12-06 | CRAN (R 3.3.0) |
| backports  |     | 1.0.5   | 2017-01-18 | CRAN (R 3.3.2) |
| colorspace |     | 1.3-2   | 2016-12-14 | CRAN (R 3.3.2) |
| DBI        |     | 0.5-1   | 2016-09-10 | CRAN (R 3.3.0) |
| devtools   |     | 1.12.0  | 2016-06-24 | CRAN (R 3.3.0) |
| digest     |     | 0.6.12  | 2017-01-27 | CRAN (R 3.3.1) |
| dplyr      | \*  | 0.5.0   | 2016-06-24 | CRAN (R 3.3.0) |
| evaluate   |     | 0.10    | 2016-10-11 | cran (@0.10)   |
| ggplot2    | \*  | 2.2.1   | 2016-12-30 | CRAN (R 3.3.2) |
| gridExtra  | \*  | 2.2.1   | 2016-02-29 | CRAN (R 3.3.0) |
| gtable     |     | 0.2.0   | 2016-02-26 | CRAN (R 3.3.0) |
| htmltools  |     | 0.3.5   | 2016-03-21 | CRAN (R 3.3.0) |
| knitr      |     | 1.15.1  | 2016-11-22 | CRAN (R 3.3.2) |
| labeling   |     | 0.3     | 2014-08-23 | CRAN (R 3.3.0) |
| lattice    |     | 0.20-34 | 2016-09-06 | CRAN (R 3.3.0) |
| lazyeval   |     | 0.2.0   | 2016-06-12 | CRAN (R 3.3.0) |
| magrittr   |     | 1.5     | 2014-11-22 | CRAN (R 3.3.0) |
| MASS       | \*  | 7.3-45  | 2016-04-21 | CRAN (R 3.3.1) |
| memoise    |     | 1.0.0   | 2016-01-29 | CRAN (R 3.3.0) |
| munsell    |     | 0.4.3   | 2016-02-13 | CRAN (R 3.3.0) |
| plyr       | \*  | 1.8.4   | 2016-06-08 | CRAN (R 3.3.0) |
| qtlcharts  | \*  | 0.7-8   | 2016-06-28 | CRAN (R 3.3.0) |
| R6         |     | 2.2.0   | 2016-10-05 | cran (@2.2.0)  |
| Rcpp       |     | 0.12.9  | 2017-01-14 | CRAN (R 3.3.2) |
| rmarkdown  |     | 1.3     | 2016-12-21 | CRAN (R 3.3.2) |
| rprojroot  |     | 1.2     | 2017-01-16 | CRAN (R 3.3.2) |
| rsconnect  |     | 0.7     | 2016-12-21 | CRAN (R 3.3.2) |
| scales     | \*  | 0.4.1   | 2016-11-09 | CRAN (R 3.3.2) |
| sp         |     | 1.2-4   | 2016-12-22 | CRAN (R 3.3.2) |
| stringi    | \*  | 1.1.2   | 2016-10-01 | cran (@1.1.2)  |
| stringr    |     | 1.1.0   | 2016-08-19 | cran (@1.1.0)  |
| tibble     |     | 1.2     | 2016-08-26 | CRAN (R 3.3.0) |
| withr      |     | 1.0.2   | 2016-06-20 | CRAN (R 3.3.0) |
| yaml       |     | 2.1.14  | 2016-11-12 | cran (@2.1.14) |

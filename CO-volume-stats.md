colorado-volume-stats
================
Jadyn Park
6/23/2021

The goal of the current script is to combine CO and NU data into a
master spreadsheet. In each CO and NU folder, there are multiple files,
one for each variable/hemisphere(left hem volume, left hem thickness,
right hem volume, right hem thickness, etc.). If there are duplicate
subjects (i.e., subjects that exist in both CO and NU dataset), the
script drops the CO and goes with NU data since the QAs done at NU were
more comprehensive.

To merge the two datasets, they need to share the same subject ID.
However, subject IDs are organized as XXXX\_fs, XXXX\_2\_fs, XXXX\_3\_fs
for CO and XXXX\_MR1, XXXX\_MR2, XXXX\_MR3 for NU. So, before we go
ahead and merge the two, we need to change the \_fs to \_MR.

Because we want to change the subject IDs at once, the plan is to rbind
(i.e., stack the individual files horizontally) and then change the
strings in the subject ID column. BUT area and thickness files have an
extra variable in there, which makes them incompatible with other
files–you can’t rbind volume (35 variables) files with area (36
variables) files.

So, my solution here is to treat these files (volume, curvind, foldind
vs. area, thickness) separately. I’ve also treated left vs. right
hemispheres separately (I don’t know why–just more organized that way,
though it may be more work).

This chunk cleans left hem from CO data (curvind, foldind, gauscurv,
meancurv, thickSD, vol)

``` r
## Set working directory
setwd("~/Downloads/Volume_Stats/")

## load libraries
rm(list = ls())

library(reshape2); library(data.table); library(dplyr); library(ggplot2); library(ggpubr); library(rstatix); library(plyr); library(tidyr)

## List files
co_listfile <- list.files("~/Downloads/Volume_Stats/CO_Prepped_FreeSurferVersion_5.3/", pattern="txt", full.names=T, recursive=T)

## Extract the  files with folder name lh (left hemisphere)
co_listfile.lh <- co_listfile[grep("lh",co_listfile)]

## Inspect the file names
head(co_listfile.lh)

## Combine all txt files in co_listfile.lh and store in dataframe 'co_lh'
for (i in 1:length(co_listfile.lh)){
if(i==1){
  assign(paste0("co_lh"), read.table(co_listfile.lh[i],header = FALSE, sep = "\t"))
}

if(!i==1){

  assign(paste0("Test",i), read.table(co_listfile.lh[i],header = FALSE, sep = "\t"))
  co_lh <- rbind(co_lh,get(paste0("Test",i)))
  rm(list = ls(pattern = "Test"))
}
}

rm(list = ls(pattern = "list.+?"))

## Replace strings in CO data so that it matches NU data
## For example, 1001's baseline and follow-up MRI are set up as 1001_fs and 1001_2_fs in CO data and is 1001_MR1 and 1001_MR2 in NU
## if XXXX_fs=MR1, XXXX_2_fs=MR2, XXXX_3_fs=MR3
co_lh$V1 <- gsub("_fs", "_MR1", co_lh$V1)
co_lh$V1 <- gsub("_2_MR1", "_MR2", co_lh$V1)
co_lh$V1 <- gsub("_3_MR1", "_MR3", co_lh$V1)

## Split co_lh into multiple dataframes
lh.curvind <- co_lh[1:270,]
lh.foldind <- co_lh[271:540,]
lh.gauscurv <- co_lh[541:810,]
lh.meancurv <- co_lh[811:1080,]
lh.thickSD <- co_lh[1081:1350,]
lh.vol <- co_lh[1351:1620,]

## Make first row the header
names(lh.curvind) <- lh.curvind[1,]
lh.curvind <- lh.curvind[-1,]
names(lh.foldind) <- lh.foldind[1,]
lh.foldind <- lh.foldind[-1,]
names(lh.gauscurv) <- lh.gauscurv[1,]
lh.gauscurv <- lh.gauscurv[-1,]
names(lh.meancurv) <- lh.meancurv[1,]
lh.meancurv <- lh.meancurv[-1,]
names(lh.thickSD) <- lh.thickSD[1,]
lh.thickSD <- lh.thickSD[-1,]
names(lh.vol) <- lh.vol[1,]
lh.vol <- lh.vol[-1,]

## Now we have a dataframe associated with each txt file
```

This chunk cleans left hem from CO data (area, thickness)

``` r
## List files
co_listfile <- list.files("~/Downloads/Volume_Stats/co_area/", pattern="txt", full.names=T, recursive=T)

## Extract the  files with folder name lh (left hemisphere)
co_listfile.lh <- co_listfile[grep("lh",co_listfile)]

## Inspect the file names
head(co_listfile.lh)

## Combine all txt files in co_listfile.lh and store in dataframe 'co_lh'
for (i in 1:length(co_listfile.lh)){
if(i==1){
  assign(paste0("co_lh_area"), read.table(co_listfile.lh[i],header = FALSE, sep = "\t"))
}

if(!i==1){

  assign(paste0("Test",i), read.table(co_listfile.lh[i],header = FALSE, sep = "\t"))
  co_lh_area <- rbind(co_lh_area,get(paste0("Test",i)))
  rm(list = ls(pattern = "Test"))
}
}

rm(list = ls(pattern = "list.+?"))

## Replace strings in CO data so that it matches NU data
## For example, 1001's baseline and follow-up MRI are set up as 1001_fs and 1001_2_fs in CO data and is 1001_MR1 and 1001_MR2 in NU
## if XXXX_fs=MR1, XXXX_2_fs=MR2, XXXX_3_fs=MR3
co_lh_area$V1 <- gsub("_fs", "_MR1", co_lh_area$V1)
co_lh_area$V1 <- gsub("_2_MR1", "_MR2", co_lh_area$V1)
co_lh_area$V1 <- gsub("_3_MR1", "_MR3", co_lh_area$V1)

## Split co_lh into multiple dataframes
lh.area <- co_lh_area[1:270,]
lh.thickness <- co_lh_area[271:540,]

## Make first row the header
names(lh.area) <- lh.area[1,]
lh.area <- lh.area[-1,]
names(lh.thickness) <- lh.thickness[1,]
lh.thickness <- lh.thickness[-1,]
```

This chunk cleans left hem from NU data (curvind, thickSD, volume)

``` r
## List files
nu_listfile <- list.files("~/Downloads/Volume_Stats/NU_Prepped_FreeSurfer5.30/", pattern="txt", full.names=T, recursive=T)

## Extract the  files with folder name lh (left hemisphere)
nu_listfile.lh <- nu_listfile[grep("lh",nu_listfile)]

## Inspect the file names
head(nu_listfile.lh)

## Combine all txt files in nu_listfile.lh and store in dataframe 'nu_lh'
for (i in 1:length(nu_listfile.lh)){
if(i==1){
  assign(paste0("nu_lh"), read.table(nu_listfile.lh[i],header = FALSE, sep = "\t"))
}

if(!i==1){

  assign(paste0("Test",i), read.table(nu_listfile.lh[i],header = FALSE, sep = "\t"))
  nu_lh <- rbind(nu_lh,get(paste0("Test",i)))
  rm(list = ls(pattern = "Test"))
}
}

rm(list = ls(pattern = "list.+?"))

## Split co_lh into multiple dataframes
lh.nu.curvind <- nu_lh[1:270,]
lh.nu.thickSD <- nu_lh[271:540,]
lh.nu.vol <- nu_lh[541:810,]

## Make first row the header
names(lh.nu.curvind) <- lh.nu.curvind[1,]
lh.nu.curvind <- lh.nu.curvind[-1,]
names(lh.nu.thickSD) <- lh.nu.thickSD[1,]
lh.nu.thickSD <- lh.nu.thickSD[-1,]
names(lh.nu.vol) <- lh.nu.vol[1,]
lh.nu.vol <- lh.nu.vol[-1,]
```

This chunk cleans left hem from NU data (area, thickness)

``` r
## List files
nu_listfile <- list.files("~/Downloads/Volume_Stats/nu_area/", pattern="txt", full.names=T, recursive=T)

## Extract the  files with folder name lh (left hemisphere)
nu_listfile.lh <- nu_listfile[grep("lh",nu_listfile)]

## Inspect the file names
head(nu_listfile.lh)

## Combine all txt files in co_listfile.lh and store in dataframe 'co_lh'
for (i in 1:length(nu_listfile.lh)){
if(i==1){
  assign(paste0("nu_lh_area"), read.table(nu_listfile.lh[i],header = FALSE, sep = "\t"))
}

if(!i==1){

  assign(paste0("Test",i), read.table(nu_listfile.lh[i],header = FALSE, sep = "\t"))
  nu_lh_area <- rbind(nu_lh_area,get(paste0("Test",i)))
  rm(list = ls(pattern = "Test"))
}
}

rm(list = ls(pattern = "list.+?"))

## Split co_lh into multiple dataframes
lh.nu.area <- nu_lh_area[1:270,]
lh.nu.thickness <- nu_lh_area[271:540,]

## Make first row the header
names(lh.nu.area) <- lh.nu.area[1,]
lh.nu.area <- lh.nu.area[-1,]
names(lh.nu.thickness) <- lh.nu.thickness[1,]
lh.nu.thickness <- lh.nu.thickness[-1,]
```

This chunk binds ALL left hemisphere data together.

``` r
## Area, Curvind, Thickness, ThicknessSD, Vol are in both CO and NU
## Remove subjects in CO datafile
new.area <- lh.area[!(lh.area$lh.aparc.area %in% lh.nu.area$lh.aparc.area),]
new.curvind <- lh.curvind[!(lh.curvind$lh.aparc.curvind %in% lh.nu.curvind$lh.aparc.curvind),]
new.thick <- lh.thickness[!(lh.thickness$lh.aparc.thickness %in% lh.nu.thickness$lh.aparc.thickness),]
new.thickSD <- lh.thickSD[!(lh.thickSD$lh.aparc.thicknessstd %in% lh.nu.thickSD$lh.aparc.thicknessstd),]
new.vol <- lh.vol[!(lh.vol$lh.aparc.volume %in% lh.nu.vol$lh.aparc.volume),]

## Merge NU and CO lh data
lh.nu.area <- rbind(lh.nu.area, new.area)
lh.nu.area <- lh.nu.area[order(lh.nu.area$lh.aparc.area),] ## Reorder so that CO rows are interleaved within NU data

lh.nu.curvind <- rbind(lh.nu.curvind, new.curvind)
lh.nu.curvind <- lh.nu.curvind[order(lh.nu.curvind$lh.aparc.curvind),]

lh.nu.thickness <- rbind(lh.nu.thickness, new.thick)
lh.nu.thickness <- lh.nu.thickness[order(lh.nu.thickness$lh.aparc.thickness),]

lh.nu.thickSD <- rbind(lh.nu.thickSD, new.thickSD)
lh.nu.thickSD <- lh.nu.thickSD[order(lh.nu.thickSD$lh.aparc.thicknessstd),]

lh.nu.vol <- rbind(lh.nu.vol, new.vol)
lh.nu.vol <- lh.nu.vol[order(lh.nu.vol$lh.aparc.volume),]

## Foldind, Gauscurv, Meancurv are only in CO
## Combine area, curvind, thickness, thickSD, vol, foldind, gauscurv, meancurv

## Merging multiple dataframes
## 1. Create a list
## lh.list <- list(lh.nu.area, lh.nu.curvind, lh.nu.thickness, lh.nu.thickSD, lh.nu.vol, lh.foldind, lh.gauscurv, lh.meancurv)

## 2. use 'Reduce' function
## lh.bind <- Reduce(function(x,y,...) merge(x,y,all=T,...), lh.list)

## ^This returns error "vector memory exhausted" so I'll just merge two dataframes at a time instead
lh.bind <- merge(lh.nu.area, lh.nu.curvind, by.x=1, by.y=1, all.x=T, all.y=T)
lh.bind <- merge(lh.bind, lh.nu.thickness, by.x=1, by.y=1, all.x=T, all.y=T)
lh.bind <- merge(lh.bind, lh.nu.thickSD, by.x=1, by.y=1, all.x=T, all.y=T)
lh.bind <- merge(lh.bind, lh.nu.vol, by.x=1, by.y=1, all.x=T, all.y=T)
lh.bind <- merge(lh.bind, lh.foldind, by.x=1, by.y=1, all.x=T, all.y=T)
lh.bind <- merge(lh.bind, lh.gauscurv, by.x=1, by.y=1, all.x=T, all.y=T)
lh.bind <- merge(lh.bind, lh.meancurv, by.x=1, by.y=1, all.x=T, all.y=T)
```

This chunk cleans right hem from CO data (curvind, foldind, gauscurv,
meancurv, thickSD, vol)

``` r
## List files
co_listfile <- list.files("~/Downloads/Volume_Stats/CO_Prepped_FreeSurferVersion_5.3/", pattern="txt", full.names=T, recursive=T)

## Extract the files with folder name rh (right hemisphere)
co_listfile.rh <- co_listfile[grep("rh",co_listfile)]

## Inspect the file names
head(co_listfile.rh)

## Combine all txt files in co_listfile.rh and store in dataframe 'co_rh'
for (i in 1:length(co_listfile.rh)){
if(i==1){
  assign(paste0("co_rh"), read.table(co_listfile.rh[i],header = FALSE, sep = "\t"))
}

if(!i==1){

  assign(paste0("Test",i), read.table(co_listfile.rh[i],header = FALSE, sep = "\t"))
  co_rh <- rbind(co_rh,get(paste0("Test",i)))
  rm(list = ls(pattern = "Test"))
}
}

rm(list = ls(pattern = "list.+?"))

## Replace strings in CO data so that it matches NU data
## For example, 1001's baseline and follow-up MRI are set up as 1001_fs and 1001_2_fs in CO data and is 1001_MR1 and 1001_MR2 in NU
## if XXXX_fs=MR1, XXXX_2_fs=MR2, XXXX_3_fs=MR3
co_rh$V1 <- gsub("_fs", "_MR1", co_rh$V1)
co_rh$V1 <- gsub("_2_MR1", "_MR2", co_rh$V1)
co_rh$V1 <- gsub("_3_MR1", "_MR3", co_rh$V1)

## Split co_rh into multiple dataframes
rh.curvind <- co_rh[1:270,]
rh.foldind <- co_rh[271:540,]
rh.gauscurv <- co_rh[541:810,]
rh.meancurv <- co_rh[811:1080,]
rh.thickSD <- co_rh[1081:1350,]
rh.vol <- co_rh[1351:1620,]

## Make first row the header
names(rh.curvind) <- rh.curvind[1,]
rh.curvind <- rh.curvind[-1,]
names(rh.foldind) <- rh.foldind[1,]
rh.foldind <- rh.foldind[-1,]
names(rh.gauscurv) <- rh.gauscurv[1,]
rh.gauscurv <- rh.gauscurv[-1,]
names(rh.meancurv) <- rh.meancurv[1,]
rh.meancurv <- rh.meancurv[-1,]
names(rh.thickSD) <- rh.thickSD[1,]
rh.thickSD <- rh.thickSD[-1,]
names(rh.vol) <- rh.vol[1,]
rh.vol <- rh.vol[-1,]
```

This chunk cleans right hem from CO data (area, thickness)

``` r
## List files
co_listfile <- list.files("~/Downloads/Volume_Stats/co_area/", pattern="txt", full.names=T, recursive=T)

## Extract the  files with folder name rh (right hemisphere)
co_listfile.rh <- co_listfile[grep("rh",co_listfile)]

## Inspect the file names
head(co_listfile.rh)

## Combine all txt files in co_listfile.rh and store in dataframe 'co_rh'
for (i in 1:length(co_listfile.rh)){
if(i==1){
  assign(paste0("co_rh_area"), read.table(co_listfile.rh[i],header = FALSE, sep = "\t"))
}

if(!i==1){

  assign(paste0("Test",i), read.table(co_listfile.rh[i],header = FALSE, sep = "\t"))
  co_rh_area <- rbind(co_rh_area,get(paste0("Test",i)))
  rm(list = ls(pattern = "Test"))
}
}

rm(list = ls(pattern = "list.+?"))

## Replace strings in CO data so that it matches NU data
## For example, 1001's baseline and follow-up MRI are set up as 1001_fs and 1001_2_fs in CO data and is 1001_MR1 and 1001_MR2 in NU
## if XXXX_fs=MR1, XXXX_2_fs=MR2, XXXX_3_fs=MR3
co_rh_area$V1 <- gsub("_fs", "_MR1", co_rh_area$V1)
co_rh_area$V1 <- gsub("_2_MR1", "_MR2", co_rh_area$V1)
co_rh_area$V1 <- gsub("_3_MR1", "_MR3", co_rh_area$V1)

## Split co_rh into multiple dataframes
rh.area <- co_rh_area[1:270,]
rh.thickness <- co_rh_area[271:540,]

## Make first row the header
names(rh.area) <- rh.area[1,]
rh.area <- rh.area[-1,]
names(rh.thickness) <- rh.thickness[1,]
rh.thickness <- rh.thickness[-1,]
```

This chunk cleans right hem from NU data (curvind, thickSD, volume)

``` r
## List files
nu_listfile <- list.files("~/Downloads/Volume_Stats/NU_Prepped_FreeSurfer5.30/", pattern="txt", full.names=T, recursive=T)

## Extract the  files with folder name rh (right hemisphere)
nu_listfile.rh <- nu_listfile[grep("rh",nu_listfile)]

## Inspect the file names
head(nu_listfile.rh)

## Combine all txt files in nu_listfile.rh and store in dataframe 'nu_rh'
for (i in 1:length(nu_listfile.rh)){
if(i==1){
  assign(paste0("nu_rh"), read.table(nu_listfile.rh[i],header = FALSE, sep = "\t"))
}

if(!i==1){

  assign(paste0("Test",i), read.table(nu_listfile.rh[i],header = FALSE, sep = "\t"))
  nu_rh <- rbind(nu_rh,get(paste0("Test",i)))
  rm(list = ls(pattern = "Test"))
}
}

rm(list = ls(pattern = "list.+?"))

## Split co_rh into multiple dataframes
rh.nu.curvind <- nu_rh[1:270,]
rh.nu.thickSD <- nu_rh[271:540,]
rh.nu.vol <- nu_rh[541:810,]

## Make first row the header
names(rh.nu.curvind) <- rh.nu.curvind[1,]
rh.nu.curvind <- rh.nu.curvind[-1,]
names(rh.nu.thickSD) <- rh.nu.thickSD[1,]
rh.nu.thickSD <- rh.nu.thickSD[-1,]
names(rh.nu.vol) <- rh.nu.vol[1,]
rh.nu.vol <- rh.nu.vol[-1,]
```

This chunk cleans right hem from NU data (area, thickness)

``` r
## List files
nu_listfile <- list.files("~/Downloads/Volume_Stats/nu_area/", pattern="txt", full.names=T, recursive=T)

## Extract the  files with folder name rh (right hemisphere)
nu_listfile.rh <- nu_listfile[grep("rh",nu_listfile)]

## Inspect the file names
head(nu_listfile.rh)

## Combine all txt files in co_listfile.rh and store in dataframe 'co_rh'
for (i in 1:length(nu_listfile.rh)){
if(i==1){
  assign(paste0("nu_rh_area"), read.table(nu_listfile.rh[i],header = FALSE, sep = "\t"))
}

if(!i==1){

  assign(paste0("Test",i), read.table(nu_listfile.rh[i],header = FALSE, sep = "\t"))
  nu_rh_area <- rbind(nu_rh_area,get(paste0("Test",i)))
  rm(list = ls(pattern = "Test"))
}
}

rm(list = ls(pattern = "list.+?"))

## Split co_rh into multiple dataframes
rh.nu.area <- nu_rh_area[1:270,]
rh.nu.thickness <- nu_rh_area[271:540,]

## Make first row the header
names(rh.nu.area) <- rh.nu.area[1,]
rh.nu.area <- rh.nu.area[-1,]
names(rh.nu.thickness) <- rh.nu.thickness[1,]
rh.nu.thickness <- rh.nu.thickness[-1,]
```

This chunk binds ALL right hemisphere data together.

``` r
## Area, Curvind, Thickness, ThicknessSD, Vol are in both CO and NU
## Remove subjects in CO datafile
new.area <- rh.area[!(rh.area$rh.aparc.area %in% rh.nu.area$rh.aparc.area),]
new.curvind <- rh.curvind[!(rh.curvind$rh.aparc.curvind %in% rh.nu.curvind$rh.aparc.curvind),]
new.thick <- rh.thickness[!(rh.thickness$rh.aparc.thickness %in% rh.nu.thickness$rh.aparc.thickness),]
new.thickSD <- rh.thickSD[!(rh.thickSD$rh.aparc.thicknessstd %in% rh.nu.thickSD$rh.aparc.thicknessstd),]
new.vol <- rh.vol[!(rh.vol$rh.aparc.volume %in% rh.nu.vol$rh.aparc.volume),]

## Merge NU and CO lh data
rh.nu.area <- rbind(rh.nu.area, new.area)
rh.nu.area <- rh.nu.area[order(rh.nu.area$rh.aparc.area),] ## Reorder so that CO rows are interleaved within NU data

rh.nu.curvind <- rbind(rh.nu.curvind, new.curvind)
rh.nu.curvind <- rh.nu.curvind[order(rh.nu.curvind$rh.aparc.curvind),]

rh.nu.thickness <- rbind(rh.nu.thickness, new.thick)
rh.nu.thickness <- rh.nu.thickness[order(rh.nu.thickness$rh.aparc.thickness),]

rh.nu.thickSD <- rbind(rh.nu.thickSD, new.thickSD)
rh.nu.thickSD <- rh.nu.thickSD[order(rh.nu.thickSD$rh.aparc.thicknessstd),]

rh.nu.vol <- rbind(rh.nu.vol, new.vol)
rh.nu.vol <- rh.nu.vol[order(rh.nu.vol$rh.aparc.volume),]

## Merging two dataframes at a time
rh.bind <- merge(rh.nu.area, rh.nu.curvind, by.x=1, by.y=1, all.x=T, all.y=T)
rh.bind <- merge(rh.bind, rh.nu.thickness, by.x=1, by.y=1, all.x=T, all.y=T)
rh.bind <- merge(rh.bind, rh.nu.thickSD, by.x=1, by.y=1, all.x=T, all.y=T)
rh.bind <- merge(rh.bind, rh.nu.vol, by.x=1, by.y=1, all.x=T, all.y=T)
rh.bind <- merge(rh.bind, rh.foldind, by.x=1, by.y=1, all.x=T, all.y=T)
rh.bind <- merge(rh.bind, rh.gauscurv, by.x=1, by.y=1, all.x=T, all.y=T)
rh.bind <- merge(rh.bind, rh.meancurv, by.x=1, by.y=1, all.x=T, all.y=T)
```

This chunk imports CO aseg dataset

``` r
## import CO aseg
aseg <- read.delim("~/Downloads/Volume_Stats/CO_Prepped_FreeSurferVersion_5.3/CO_Prep_Volume_Aseg_v5.3.txt")

## replace _fs with _MR
aseg$Measure.volume <- gsub("_fs", "_MR1", aseg$Measure.volume)
aseg$Measure.volume <- gsub("_2_MR1", "_MR2", aseg$Measure.volume)
aseg$Measure.volume <- gsub("_3_MR1", "_MR3", aseg$Measure.volume)

## add suffix _aseg to each column
colnames(aseg) <- paste(colnames(aseg), "aseg", sep="_")
```

This chunk binds lh, rh, and aseg together, creating a master dataset

``` r
master <- merge(lh.bind, rh.bind, by.x=1, by.y=1, all.x=T, all.y=T)
master <- merge(master, aseg, by.x=1, by.y=1, all.x=T, all.y=T)

write.csv(master, "~/Desktop/volume_stats_master.csv")
```

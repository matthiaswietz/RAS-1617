
############################################################################################
    ###  RAS 2016-17 - AMPLICON ANALYSIS  ###
############################################################################################

# This script: load and format dada2-amplicons 
# Load and format metadata

# NOTE: imported ASVs represent four sample sets from RAS 2016-17
# Only samples from moorings F4 and EGC analyzed here
# Moorings HGIV and Fevi34 will be part of a separate study

####################################

# Set working directory; load data
setwd("/AWI_MPI/FRAM/RAS_PPS/2016-17/Rstats")
load("RAS.Rdata")

# Load packages and colors
library(gtools)
library(dplyr)
library(data.table)
library(solartime)
library(oce)
source("RAS_Colors.R")


############################################################################################
  ###  ASV DATA -- LOAD AND FORMATTING  ###
############################################################################################

## ASV IMPORT ##

ASV.bac <- read.table(
  "bac_seqtab.txt",
  h = T,
  sep = "\t",
  check.names = F)

ASV.euk <- read.table(
  "euk_seqtab.txt",
  h = T,
  sep = "\t",
  check.names = F)

############################################################################################

## TAXONOMY ##

TAX.bac <- read.table(
  "bac_tax.txt",
  h = T, 
  sep = "\t", 
  stringsAsFactors = F, 
  row.names = 1)

TAX.euk.raw <- read.table(
  "euk_tax.txt",
  h = T, 
  sep = "\t", 
  stringsAsFactors = F,
  row.names = 1)

####################################

## REFORMAT EUKs

# Delete superphylum/species taxrank
TAX.euk <- TAX.euk.raw %>%
  dplyr::select(-c("Supergroup","Species"))
  
# Export Animalia/Metazoa to new DF 
TAX.meta <- TAX.euk[grep("Metazoa", TAX.euk$Phylum),]
ASV.meta <- merge(ASV.euk, TAX.meta, by = "row.names", all.x=F)
ASV.meta <- ASV.meta[mixedorder(ASV.meta$Row.names),]
rownames(ASV.meta) = ASV.meta$Row.names
ASV.meta <- ASV.meta[, -c(1, 98:103)]

# Remove Animalia/Metazoa from org table
TAX.euk <- TAX.euk[-grep("Metazoa", TAX.euk$Phylum),]

# Reformat to match
ASV.euk <- merge(ASV.euk, TAX.euk, by = "row.names", all.x=F)
ASV.euk <- ASV.euk[mixedorder(ASV.euk$Row.names),]
rownames(ASV.euk) = ASV.euk$Row.names
ASV.euk <- ASV.euk[, -c(1, 98:103)]


############################################################################################
   ###  METADATA  ###
############################################################################################

## Import metadata -- season, nutrients etc ##
ENV <- read.table(
  "metadata.txt", h=T, sep = "\t", 
  stringsAsFactors=F, skipNul=T) %>%
  mutate(date=as.Date(
    date, format = "%Y-%m-%d"))  

# Calculate pressure
ENV$dba <- swPressure(
  ENV$depth, 
  latitude=80, eos="unesco")

# Calculate density
ENV$rho <- swRho(
  ENV$salinity, 
  temperature = ENV$temp, 
  pressure = ENV$dba,
  longitude = NULL, 
  latitude = 80, 
  eos = getOption(
    "oceEOS", default="unesco"))

# Add stratification; convert number for plotting
strat <- read.table(
  "Stratification.txt",
  h = T, sep = "\t",
  check.names = F) %>%
  mutate(date=as.Date(
    date, format = "%Y-%m-%d")) 

strat <- strat %>%
  group_by(date, biome) %>%
  summarize_at(c("strat"), mean) %>%
  mutate(strat=strat*100000) %>%
  ungroup

# Add daylight hours
ENV <- ENV %>%
  mutate(daylight = computeDayLength(
    ENV$date, ENV$lat))

# append to ENV; remove outliers (RAS pushed to >200m depth)
ENV <- ENV %>% 
  inner_join(strat, by = c(
    "date","biome")) %>% 
  filter(!RAS_id %in% c(
    "03_1617_F4_1","04_1617_EGC_2"))

# Generate factors for parameter of interest
ENV$type <- factor(
  ENV$type, levels=c(
    "euk","bac"))
ENV$biome <- factor(
  ENV$biome, levels=c(
    "WSC","EGC"))
ENV$season <- factor(
  ENV$season, levels=c(
    "autumn","winter","spring","summer"))
ENV$month <- factor(
  ENV$month, levels=c(
    "Aug","Sep","Oct","Nov","Dec","Jan",
    "Feb","Mar","Apr","May","Jun","Jul"))
ENV$watermass <- factor(
  ENV$watermass, levels=c(
    "AW","PW","mix"))

# Subset BAC+EUK
ENV.bac <- subset(ENV, type=="bac")
ENV.euk <- subset(ENV, type=="euk") 


############################################################################################
   ###  FORMAT + EXPORT  ###
############################################################################################

# sort  
ASV.bac <- ASV.bac[,mixedsort(names(ASV.bac))]
ASV.euk <- ASV.euk[,mixedsort(names(ASV.euk))]
ENV.bac <- ENV.bac[mixedorder(ENV.bac$clip_id),]
ENV.euk <- ENV.euk[mixedorder(ENV.euk$clip_id),]

ASV.bac <- ASV.bac[,c((match(
  ENV.bac$clip_id, colnames(ASV.bac))))]
ASV.euk <- ASV.euk[,c((
  match(ENV.euk$clip_id, colnames(ASV.euk))))]

# rename 
colnames(ASV.bac) = ENV.bac$RAS_id
colnames(ASV.euk) = ENV.euk$RAS_id

# Add BAC+EUK identifiers
row.names(ASV.bac) <- paste0("bac_", row.names(ASV.bac))
row.names(ASV.euk) <- paste0("euk_", row.names(ASV.euk))
row.names(TAX.bac) <- paste0("bac_", row.names(TAX.bac))
row.names(TAX.euk) <- paste0("euk_", row.names(TAX.euk))

# Reorder columns to match
TAX.bac <- setcolorder(TAX.bac, c(
  "Kingdom","Phylum","Class","Order",
  "Family","Genus"))
TAX.euk <- setcolorder(TAX.euk, c(
  "Kingdom","Phylum","Class","Order",
  "Family","Genus"))

# Verify same roworder of ENV and ASV 
ENV <- ENV[mixedorder(ENV$RAS_id),]
ENV.bac <- ENV.bac[mixedorder(ENV.bac$RAS_id),]
ENV.euk <- ENV.euk[mixedorder(ENV.euk$RAS_id),]
ASV.bac <- ASV.bac[, ENV.bac$RAS_id]
ASV.euk <- ASV.euk[, ENV.euk$RAS_id]

# Add ENV rownames for Phyloseq
row.names(ENV.bac) = ENV.bac$RAS_id
row.names(ENV.euk) = ENV.euk$RAS_id

# Vector for month conversion
month.long = c(
  "2016-08","2016-09","2016-10","2016-11",
  "2016-12","2017-01","2017-02","2017-03",
  "2017-04","2017-05","2017-06","2017-07")


############################################################################################
   ### ANNOTATION COLORS ##
############################################################################################

anno = data.frame(
  "season" = factor(ENV.euk$season),
  "biome" = factor(ENV.euk$biome),
  "month" = factor(ENV.euk$month)) 
rownames(anno) = ENV.euk$RAS_id
anno <- anno[order(anno$biome),]

anno_col = list(
  biome = c(
    "WSC"="indianred2","EGC"="deepskyblue3"),
  season = c(
    "spring"="#ffcc00","summer"="darkorange",
    "autumn"="#aa0000","winter"="#50516b"))


############################################################################################

save.image("RAS.Rdata")

############################################################################################

## continued in Rscript *DataProcessing*
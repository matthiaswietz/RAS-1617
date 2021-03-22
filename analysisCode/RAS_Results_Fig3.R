
############################################################################################
    ###  RAS2016-17 -- AMPLICONS -- FIG. 3 ###
############################################################################################

# set working directory; load data
setwd("/AWI_MPI/FRAM/RAS_PPS/2016-17/Rstats")
load("RAS.Rdata")

# load packages and colors
library(ampvis2)
library(PMCMR)
library(psych)
library(dplyr)
library(tibble)
library(tidyr)
library(ggplot2)
library(ggridges)
library(PNWColors)
library(cowplot)
source("RAS_Colors.R")


############################################################################################
   ###  Fig. 3a --- SEASONAL ENV-PCA  ###
############################################################################################

## PCA -- EGC ##

PCA.EGC.1 <- ENV.euk %>%
  filter(biome=="EGC") %>%
  dplyr::select(-c(
    "salinity","chl_a","CO2","pH","dba",
    "depth","sig","ice_past","strat",
    "full_id","clip_id","mooring",
    "lat","lon","month_full","biome",
    "type","date")) %>%
  drop_na()

dplyr::select(PCA.EGC.1, 
              season, RAS_id, month, 
              watermass, ice_act_conc) -> PCA.EGC.2

PCA.EGC.1  <- PCA.EGC.1  %>% 
  #tibble::column_to_rownames(c("ID_RAS")) %>%
  dplyr::select(-c(
    season, RAS_id, month, 
    watermass, ice_act_conc))

autoplot(
  prcomp(PCA.EGC.1, center=T, scale=T),
  data = PCA.EGC.2, label=F, 
  loadings=T, loadings.colour="gray55",
  loadings.label=T, loadings.label.size=2,
  loadings.label.repel=T, 
  loadings.label.vjust=1.2,
  loadings.label.colour="#b3b3b3") +
  geom_point(aes(
    shape=watermass, 
    color=season, size=ice_act)) +
  scale_shape_manual(values=c(15,17,16)) +
  scale_size(
    range = c(6,11),
    breaks=c(0,40,80)) +
  geom_text(aes(
    label=month), 
    color="white", size=4, vjust=0.5) +
  scale_color_manual(values = c(col.env)) +
  guides(color=F) +
  #scale_x_reverse() + 
  theme_bw() +
  theme(axis.ticks = element_blank(),
        axis.text = element_blank(),
        #panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        legend.key.size = unit(0.05, "cm"),
        legend.position = "bottom")

####################################

PCA.WSC.1 <- ENV.euk %>%
  filter(biome=="WSC") %>%
  dplyr::select(-c(
    salinity, 
    "chl_a","CO2","AW","PW","pH","dba","rho",
    "depth","sig","ice_act","ice_past",
    "full_id","clip_id","mooring",
    "lat","lon","month_full","biome",
    "type","date","watermass",
    "ice_act_conc","strat")) %>%
  drop_na()

dplyr::select(PCA.WSC.1, 
              season, month, RAS_id) -> PCA.WSC.2

PCA.WSC.1  <- PCA.WSC.1  %>% 
  dplyr::select(-c(season, month, RAS_id))

autoplot(
  prcomp(PCA.WSC.1, center=T, scale=T),
  data = PCA.WSC.2, label=F, 
  #frame=T, frame.colour = 'season', 
  loadings=T, loadings.colour="#b3b3b3",
  loadings.label=T, loadings.label.size=2,
  loadings.label.vjust=1.5,
  loadings.label.colour="black") +
  geom_point(aes(
    color=season, 
    fill=season), 
    size=9, shape=22) +
  geom_text(aes(
    label=month), 
    color="white", size=4, vjust=0.5) +
  scale_shape_manual(values=c(15,17,16)) +
  scale_color_manual(values=col.env) +
  scale_fill_manual(values=col.env) +
  #scale_y_reverse() +
  theme_bw() +
  theme(axis.ticks = element_blank(),
        axis.text = element_blank(),
        #panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        legend.key.size = unit(0.05, "cm"),
        legend.position = "bottom")


############################################################################################
    ###  Fig. 3b -- SEASONAL ORDINATION  ### 
############################################################################################

NMDS.plot1 <- amp_ordinate(
  amp_subset_samples(
     ampvis.bac, biome=="WSC"),
  type = "NMDS",
  distmeasure = "jsd",
  transform = "hellinger", 
  sample_color_by = "season", sample_colorframe = T,
  sample_point_size = 2.44) +
  scale_color_manual(values=col.env) +
  scale_fill_manual(values=col.env) +
  theme_bw() + 
  theme(axis.text = element_blank(),
        axis.title = element_blank(),
        axis.ticks = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = "none")

NMDS.plot2 <- amp_ordinate(
  amp_subset_samples(
     ampvis.euk, biome=="WSC"),
  type = "NMDS",
  distmeasure = "jsd",
  transform = "hellinger", 
  sample_color_by = "season", sample_colorframe = T,
  sample_point_size = 2.44) +
  scale_color_manual(values=col.env) +
  scale_fill_manual(values=col.env) +
  theme_bw() + 
  theme(axis.text = element_blank(),
        axis.title = element_blank(),
        axis.ticks = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = "none")

####################################

NMDS.plot3 <- amp_ordinate(
  amp_subset_samples(
     ampvis.bac, biome=="EGC"),
  type = "NMDS",
  distmeasure = "jsd",
  transform = "hellinger", 
  sample_color_by = "season", sample_colorframe = T,
  sample_point_size = 2.44) +
  scale_color_manual(values=col.env) +
  scale_fill_manual(values=col.env) +
  theme_bw() + 
  theme(axis.text = element_blank(),
        axis.title = element_blank(),
        axis.ticks = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = "none")

NMDS.plot4 <- amp_ordinate(
  amp_subset_samples(
     ampvis.euk, biome=="EGC"),
  type = "NMDS",  
  transform = "hellinger", 
  distmeasure = "jsd",
  sample_color_by = "season", sample_colorframe = T,
  sample_point_size = 2.44) +
  scale_color_manual(values=col.env) +
  scale_fill_manual(values=col.env) +
  theme_bw() + 
  theme(axis.text = element_blank(),
        axis.title = element_blank(),
        axis.ticks = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = "none")

####################################

plot_grid(
  NMDS.plot1, 
  NMDS.plot3, 
  NMDS.plot2, 
  NMDS.plot4, 
  labels = c(
    "PRO-WSC","PRO-EGC", 
    "EUK-WSC","EUK-EGC"),
  label_x = 0.09,
  label_y = 0.97,
  label_size = 8,
  align = "v",
  ncol = 2)


############################################################################################
   ###  Fig. 3b -- SEASONAL DISTANCES  ### 
############################################################################################

# Calculate bacs & euks per site
dist1 <- as.data.frame(as.matrix(dist.bac.WSC)) %>% 
  rownames_to_column() %>% 
  pivot_longer(-rowname) %>%
  dplyr::rename(RAS_id = rowname) %>%
  left_join(dplyr::select(
    ENV.bac, RAS_id, season), by=c("RAS_id")) %>%
  dplyr::rename(RAS_id = name, RAS_id.2 = RAS_id) %>%
  left_join(dplyr::select(
    ENV.bac, RAS_id, season), by=c("RAS_id")) %>%
  mutate(dist = paste(season.x, season.y, sep="_")) %>%
  group_by(dist) %>%
  summarize_at(c("value"), mean) %>%
  ungroup() %>%
  separate(dist, c("season1", "season2"), sep="_") %>%
  distinct(value, .keep_all=T) %>%  
  mutate(biome="WSC", type="bac")

dist2 <- as.data.frame(as.matrix(dist.bac.EGC)) %>% 
  rownames_to_column() %>% 
  pivot_longer(-rowname) %>%
  dplyr::rename(RAS_id = rowname) %>%
  left_join(dplyr::select(
    ENV.bac, RAS_id, season), by=c("RAS_id")) %>%
  dplyr::rename(RAS_id = name, RAS_id.2 = RAS_id) %>%
  left_join(dplyr::select(
    ENV.bac, RAS_id, season), by=c("RAS_id")) %>%
  mutate(dist = paste(season.x, season.y, sep="_")) %>%
  group_by(dist) %>%
  summarize_at(c("value"), mean) %>%
  ungroup() %>%
  separate(dist, c("season1", "season2"), sep="_") %>%
  distinct(value, .keep_all=T) %>% 
  mutate(biome="EGC", type="bac")

dist3 <- as.data.frame(as.matrix(dist.euk.WSC)) %>% 
  rownames_to_column() %>% 
  pivot_longer(-rowname) %>%
  dplyr::rename(RAS_id = rowname) %>%
  left_join(dplyr::select(
    ENV.euk, RAS_id, season), by=c("RAS_id")) %>%
  dplyr::rename(RAS_id = name, RAS_id.2 = RAS_id) %>%
  left_join(dplyr::select(
    ENV.euk, RAS_id, season), by=c("RAS_id")) %>%
  mutate(dist = paste(season.x, season.y, sep="_")) %>%
  group_by(dist) %>%
  summarize_at(c("value"), mean) %>%
  ungroup() %>%
  separate(dist, c("season1", "season2"), sep="_") %>%
  distinct(value, .keep_all=T) %>%  
  mutate(biome="WSC", type="euk")

dist4 <- as.data.frame(as.matrix(dist.euk.EGC)) %>% 
  rownames_to_column() %>% 
  pivot_longer(-rowname) %>%
  dplyr::rename(RAS_id = rowname) %>%
  left_join(dplyr::select(
    ENV.euk, RAS_id, season), by=c("RAS_id")) %>%
  dplyr::rename(RAS_id = name, RAS_id.2 = RAS_id) %>%
  left_join(dplyr::select(
    ENV.euk, RAS_id, season), by=c("RAS_id")) %>%
  mutate(dist = paste(season.x, season.y, sep="_")) %>%
  group_by(dist) %>%
  summarize_at(c("value"), mean) %>%
  ungroup() %>%
  separate(dist, c("season1", "season2"), sep="_") %>%
  distinct(value, .keep_all=T) %>%
  mutate(biome="EGC", type="euk")

# Merge all + plot
dist.ssn <- rbind(
  dist1, dist2, dist3, dist4) %>%
  reshape2::melt() %>%
  mutate(biome=factor(
    biome, levels = c("WSC","EGC"))) %>%
  mutate(type=factor(
    type, levels = c("euk","bac"))) 

ggplot(dist.ssn, 
  aes(x=season1, y=season2, fill=value)) +
geom_tile() +
geom_text(
  aes(label=round(value, digits=2)), 
  color="white", size=4) +
facet_grid(type~biome) +
scale_fill_gradientn(
  colors = rev(pnw_palette("Lake", 5))) +
#coord_flip() +
theme_bw() +
theme(axis.ticks = element_blank(),
      panel.border = element_blank(),
      axis.title = element_blank())

# remove temp-data
rm(dist1, dist2, dist3, dist4)


###################################################################################

sessionInfo()
save.image("RAS.Rdata")

###################################################################################

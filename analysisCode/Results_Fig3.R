
############################################################################################
   ###  RAS2016-17 - AMPLICON ANALYSIS  ###  
   ###  Fig. 3a --- SEASONAL ENV-PCA 
############################################################################################

library(ampvis2)
library(PMCMR)
library(psych)
library(dplyr)
library(tibble)
library(tidyr)
library(ggplot2)
library(ggfortify)
library(PNWColors)
library(cowplot)
source("Colors.R")


################################
   
## PCA -- EGC ##

PCA.EGC.1 <- ENV.euk %>%
  filter(site=="EGC") %>%
  dplyr::select(-c(
    "salinity","chl_a","CO2","pH","dba",
    "depth","sig","ice_past","strat",
    "clip_id","mooring","lat","lon",
    "site","locus_tag","date")) %>%
  drop_na()

dplyr::select(
  PCA.EGC.1, 
  season, RAS_id, month, 
  watermass, ice_cover) -> PCA.EGC.2

PCA.EGC.1  <- PCA.EGC.1  %>% 
  #tibble::column_to_rownames(c("ID_RAS")) %>%
  dplyr::select(-c(
    season, RAS_id, month, 
    watermass, ice_cover))

autoplot(
  prcomp(PCA.EGC.1, center=T, scale=T),
  data = PCA.EGC.2, label=F, 
  loadings=T, loadings.colour="gray55",
  loadings.label=T, loadings.label.size=2,
  loadings.label.repel=T, 
  loadings.label.vjust=1.2,
  loadings.label.colour="black") +
geom_point(aes(
  shape=watermass, 
  color=season, 
  size=ice_conc)) +
scale_shape_manual(values=c(15,17,16)) +
scale_size(
  range = c(6,11),
  breaks=c(0,40,80)) +
geom_text(aes(
  label=month), 
  color="white", 
  size=4, vjust=0.5) +
scale_color_manual(values = c(col.env)) +
guides(color="none") +
theme_bw() +
theme(axis.ticks = element_blank(),
      axis.text = element_blank(),
      panel.border = element_blank(),
      legend.key.size = unit(0.05, "cm"),
      legend.position = "bottom")

####################################

PCA.WSC.1 <- ENV.euk %>%
  filter(site=="WSC") %>%
  dplyr::select(-c(
    salinity, 
    "chl_a","CO2","AW","PW","pH",
    "dba","rho","depth","sig",
    "ice_conc","ice_past",
    "clip_id","mooring","lat","lon",
    "site","locus_tag","date",
    "watermass","ice_cover","strat")) %>%
  drop_na()

dplyr::select(
  PCA.WSC.1, 
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
     ampvis.bac, site=="WSC"),
  type = "NMDS",
  distmeasure = "bray",
  transform = "hellinger", 
  sample_color_by = "season", sample_colorframe = T,
  sample_point_size = 3) +
  scale_color_manual(values=col.env) +
  scale_fill_manual(values=col.env) +
  theme_bw() + 
  theme(axis.text = element_blank(),
        axis.title = element_blank(),
        axis.ticks = element_blank(),
        panel.border = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = "none")

NMDS.plot2 <- amp_ordinate(
  amp_subset_samples(
     ampvis.euk, site=="WSC"),
  type = "NMDS",
  distmeasure = "bray",
  transform = "hellinger", 
  sample_color_by = "season", sample_colorframe = T,
  sample_point_size = 3) +
  scale_color_manual(values=col.env) +
  scale_fill_manual(values=col.env) +
  theme_bw() + 
  theme(axis.text = element_blank(),
        axis.title = element_blank(),
        axis.ticks = element_blank(),
        panel.border = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = "none")

####################################

NMDS.plot3 <- amp_ordinate(
  amp_subset_samples(
     ampvis.bac, site=="EGC"),
  type = "NMDS",
  distmeasure = "bray",
  transform = "hellinger", 
  sample_color_by = "season", 
  sample_colorframe = T,
  sample_point_size = 3) +
  scale_color_manual(values=col.env) +
  scale_fill_manual(values=col.env) +
  theme_bw() + 
  theme(axis.text = element_blank(),
        axis.title = element_blank(),
        axis.ticks = element_blank(),
        panel.border = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = "none")

NMDS.plot4 <- amp_ordinate(
  amp_subset_samples(
     ampvis.euk, site=="EGC"),
  type = "NMDS",  
  transform = "hellinger", 
  distmeasure = "bray",
  sample_color_by = "season", 
  sample_colorframe = T,
  sample_point_size = 3) +
  scale_color_manual(values=col.env) +
  scale_fill_manual(values=col.env) +
  theme_bw() + 
  theme(axis.text = element_blank(),
        axis.title = element_blank(),
        axis.ticks = element_blank(),
        panel.border = element_blank(),
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
  mutate(site="WSC", locus_tag="16S")

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
  mutate(site="EGC", locus_tag="16S")

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
  mutate(site="WSC", locus_tag="18S")

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
  mutate(site="EGC", locus_tag="18S")

# Merge all + plot
dist.ssn <- rbind(
   dist1, dist2, dist3, dist4) %>%
 reshape2::melt() %>%
  mutate(site=factor(
    site, levels = c("WSC","EGC"))) %>%
  mutate(type=factor(
    locus_tag, levels = c("18S","16S"))) %>%
  mutate(across(everything(),~gsub("autumn","AU", .))) %>%
  mutate(across(everything(),~gsub("winter","WI", .))) %>%
  mutate(across(everything(),~gsub("spring","SP", .))) %>%
  mutate(across(everything(),~gsub("summer","SU", .))) %>%
  mutate(across(value, as.numeric))

ggplot(dist.ssn, 
  aes(x=season1, y=season2, fill=value)) +
geom_tile() +
geom_text(
  aes(label=round(value, digits=2)), 
  color="white", size=4.4) +
facet_grid(locus_tag~site) +
scale_fill_gradientn(
  colors = rev(pnw_palette("Lake", 5))) +
theme_bw() +
theme(axis.ticks = element_blank(),
      panel.grid = element_blank(),
      panel.border = element_blank(),
      axis.title = element_blank())

# remove temp-data
rm(dist1, dist2, dist3, dist4)


###################################################################################

sessionInfo()
save.image("RAS.Rdata")

###################################################################################

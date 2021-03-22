
############################################################################################
    ###  RAS 2016-17 - AMPLICONS -- Fig.1B  ###
############################################################################################

# set working directory; load data
setwd("/AWI_MPI/FRAM/RAS_PPS/2016-17/Rstats")
load("RAS.Rdata")

# load packages and colors
library(ampvis2)
library(dplyr)
library(gtools)
library(tibble)
library(tidyr)
library(ggplot2)
library(ggfortify)
library(cowplot)
library(PNWColors)
source("RAS_Colors.R")


############################################################################################
   ###  ANNUAL ENV COMPARISON - Fig. 1B  ###
############################################################################################

ggplot(data = subset(
  ENV.melt, variable %in% c(
    "temp","O2_conc","daylight",
    "AW","NO3_NO2","ice_act"))) +
  aes(date, value, 
      group=variable, color=variable) +
  geom_line(size=1.2, na.rm=T) + 
 geom_point(size=1.2, na.rm=T) + 
  scale_color_manual(
    values = pnw_palette("Cascades",7)) +
  scale_x_date() +
  scale_y_continuous(expand = c(0.2,0.1)) +
  facet_grid(variable~biome, scales="free") +
  theme_bw() +
  theme(panel.grid.minor = element_blank(),
        axis.title = element_blank(),
        axis.text.x = element_text(
          angle=90, hjust=0.5, vjust=0.5),
        axis.text = element_text(size=6),
        legend.position = "none")



###################################################################################

sessionInfo()
save.image("RAS.Rdata")

###################################################################################
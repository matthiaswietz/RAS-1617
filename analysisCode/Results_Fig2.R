
############################################################################################
   ###  RAS2016-17 -- AMPLICON ANALYSIS  ###
   ###  Fig. 2 -- ANNUAL OVERVIEW
############################################################################################

library(phyloseq)
library(ampvis2)
library(dplyr)
library(ggplot2)
library(viridis)
library(cowplot)
source("Colors.R")

###########################################
   
# Subset taxa
year.bac <- sum.order %>%
  filter(Order %in% c(
    "Flavobacteriales","Rhodobacterales",
    "Rhodospirillales","Planctomycetales",
    "Gammaproteobacteria_uc","Arctic97B-4",
    "Nitrosopumilales")) %>%
  group_by(Order, month, locus_tag, site) %>%
  summarize_at(c("Abundance"), mean) %>% 
  ungroup %>%
  dplyr::rename(tax=Order) 

year.euk <- sum.class %>%
  filter(Class %in% c(
    "Bacillariophyta","Syndiniales",
    "RAD-C","Prymnesiophyceae",
    "Pirsonia_Clade","MAST-3")) %>%
  group_by(Class, month, locus_tag, site) %>%
  summarize_at(c("Abundance"), mean) %>% 
  ungroup %>%
  dplyr::rename(tax=Class) 

# Bind and fix date
yearplot <- rbind(year.bac, year.euk)
yearplot$date <- month.long[yearplot$month]
yearplot$date <- paste0(yearplot$date, "-01")
yearplot <- yearplot %>% 
  mutate(date = as.Date(
    date, format = "%Y-%m-%d"))

# Order appearance in plot
yearplot$tax <- factor(
  yearplot$tax, levels=c(
    "Flavobacteriales","Rhodospirillales",
    "Rhodobacterales","Nitrosopumilales",
    "Gammaproteobacteria_uc","Arctic97B-4",
    "Planctomycetales","Syndiniales","RAD-C",
    "Bacillariophyta","Prymnesiophyceae",
    "MAST-3","Pirsonia_Clade"))
yearplot$month <- factor(
  yearplot$month, levels=c(
    "Aug","Sep","Oct","Nov","Dec","Jan",
    "Feb","Mar","Apr","May","Jun","Jul"))

plot_year1 <- ggplot(data=subset(
  yearplot, locus_tag=="18S"), 
  aes(x=date, y=Abundance)) + 
geom_area(aes(fill=tax), 
  position=position_dodge(0.2), 
  alpha=0.98) +
scale_x_date(
  expand = c(0.04,0.04), 
  breaks=NULL) +
scale_y_continuous(
  expand = c(0,4),
  breaks = c(0,10,20,30,40,50,60)) +
scale_fill_manual(values = col.taxon) +
facet_grid(~site) +
theme_bw() +
theme(axis.title = element_blank(),
      axis.text.x = element_blank(), 
      legend.position = "none", 
      strip.background = element_blank(), 
      strip.text = element_blank(),
      panel.border = element_blank(),
      panel.grid.minor = element_blank(),
      plot.margin = unit(c(0,0,0,0),"cm"))

plot_year2 <- ggplot(data=subset(
  yearplot, locus_tag=="16S"), 
  aes(x=date, y=Abundance)) + 
geom_area(aes(fill=tax), 
  position=position_dodge(0.2), 
  alpha=0.98) +
scale_x_date(
  expand = c(0.04,0.04), 
  breaks=NULL) +
scale_y_continuous(
  expand = c(0,4),
  breaks = c(0,10,20,30,40,50,60)) +
scale_fill_manual(values = col.taxon) +
facet_grid(~site) +  
theme_bw() +
theme(axis.title = element_blank(),
      axis.text.x = element_blank(), 
      legend.position = "none", 
      strip.background = element_blank(), 
      strip.text = element_blank(),
      panel.border = element_blank(),
      panel.grid.minor = element_blank(),
      plot.margin = unit(c(0,0,0,0),"cm"))


############################################################################################
   ###  DIVERSITY + ENV ###
############################################################################################

plot_dist <- ggplot(div.avg, 
  aes(x=month, y=dist, 
  group=locus_tag, colour=locus_tag)) +
geom_line(size=1.5) +
geom_point(size=2.6) +
scale_x_discrete(
  expand = c(0.04,0.04)) +
scale_y_continuous(
  limits=c(0.35,1)) +
scale_color_manual(values=col.env) +
facet_grid(~site) +
theme_bw() + 
theme(legend.position = "none", 
      axis.title = element_blank(), 
      axis.text.x = element_blank(), 
      axis.text.y = element_text(size=7),
      strip.background = element_blank(), 
      strip.text = element_blank(),
      panel.border = element_blank(),
      panel.grid.minor = element_blank(),
      plot.margin = unit(c(0,0,0,0),"cm"))

plot_div <- ggplot(div.avg, 
    aes(x=month, y=simpson, 
    group=locus_tag, colour=locus_tag)) +
  geom_line(size=1.5) +
  geom_point(size=2.6) +
  scale_x_discrete(
    expand = c(0.04,0.04)) +
  scale_y_continuous(
    limits=c(0,130),
    breaks = c(40,80,120)) +
  scale_color_manual(values=col.env) +
  facet_grid(~site) +
  theme_bw() + 
  theme(legend.position = "none", 
        axis.title = element_blank(), 
        axis.text.x = element_blank(), 
        axis.text.y = element_text(size=7),
        strip.background = element_blank(), 
        strip.text = element_blank(),
        panel.border = element_blank(),
        panel.grid.minor = element_blank(),
        plot.margin = unit(c(0,0,0,0),"cm"))

plot_light <- avg.month %>%
  reshape2::melt() %>%
 filter(variable=="daylight") %>%
  ggplot() + 
  aes(x=month, y=variable, fill=value) + 
  geom_tile() +
  scale_fill_viridis(
    option="B", 
    begin=0, end=0.92,
    breaks=c(0,6,12,18,24)) +
  scale_x_discrete(
    expand = c(0.04,0.04)) +
  facet_wrap(~site) +
  theme_bw() +
  theme(legend.position = "top", 
        axis.title = element_blank(), 
        axis.text.x = element_blank(),
        axis.ticks = element_blank(),
        axis.text.y = element_text(size=7),
        strip.background = element_blank(), 
        strip.text = element_blank(),
        panel.border = element_blank(),
        panel.grid.minor = element_blank(),
        plot.margin = unit(c(0,0,0,0),"cm"))

############################################################################################

## PLOTTING and STATS ##

plot_grid(
  plot_light,
  plot_year1, 
  plot_year2,
  plot_dist,
  plot_div,
  ncol = 1,
  label_x = 0.1,
  label_size = 8,
  rel_heights = c(
    1.8, 5, 5.2, 2.2, 2.2, 2.2),
  align = "v", 
  axis = "tblr")

####################################

## Stats - daylight:
# WSC-bac 0.5** // WSC-euk: 0.8*****
# EGC-bac 0.4* // EGC-euk 0.4*

## Stats - temp:
# WSC-bac 0.8******** // WSC-euk: 0.6***
# EGC-bac 0.7**** // EGC-euk n.sig.

with(subset(div.all, 
  site==c("EGC") & locus_tag==c("16S")),
cor.test(
   dist, temp, 
   use = "complete.obs",  
   exact = F,
   method = c("spearman")))

with(subset(div.all, 
   site==c("EGC") & locus_tag==c("18S")),
cor.test(
  dist, temp, 
  use = "complete.obs",  
  exact = F,
  method = c("spearman")))

####################################

# remove temp-data
rm(year.bac, year.euk)


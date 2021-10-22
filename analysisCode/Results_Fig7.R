
############################################################################################
   ###  RAS 2016-17 -- AMPLICON ANALYSIS  ###
   ###  FIG. 7 -- SUMMER 
############################################################################################

library(ggridges)

## WSC ##

summer1 <- sum.genus %>% filter(
  site=="WSC" &
    month %in% c("Apr","May","Jun","Jul") & 
    date != "2017-07-22" &
    Genus %in% c(
      "OCS116_clade_uc","NS10",
      "NS2b","Aurantivirga","Formosa",
      "Polaribacter","Amylibacter",
      "Cyclobacteriaceae_uc",
      "Phaeocystis","Grammonema",
      "Chytriodinium","Chaetoceros",
      "Micromonas","Woloszynskia",
      "Thalassiosira","Odontella")) %>%
  mutate(date = as.Date(
    date, format = "%Y-%m-%d")) %>%
  ungroup %>%
  mutate(Genus = factor(Genus, levels=c(
    "Cyclobacteriaceae_uc","Polaribacter",
    "NS10","OCS116_clade_uc","Aurantivirga",
    "NS2b","Formosa","Amylibacter",
    "Thalassiosira","Odontella",
    "Micromonas","Chytriodinium",
    "Woloszynskia","Chaetoceros",
    "Phaeocystis","Grammonema"))) %>%
  select_if(~all(!is.na(.))) %>%  
  distinct() 

####################################

plot_summer01 <- ggplot(data=subset(
  summer1, locus_tag=="18S"), aes(
    x=date, y=Genus, 
    height=Abundance, fill=Genus)) +
  geom_ridgeline(
    alpha=0.94, 
    size=0.1,
    scale=0.25, 
    show.legend=F) +
  scale_x_continuous(
    expand = c(0.02,0.02)) +
  scale_y_discrete(
    expand = c(0.02,0.02)) +
  scale_fill_manual(
    values = rev(pnw_palette("Sailboat",8))) +
  theme_classic() +
  theme(axis.title=element_blank(),
        axis.text.x=element_blank(),
        axis.text.y=element_text(size=10),
        panel.grid = element_blank(),
        axis.ticks = element_blank())
  
plot_summer02 <- ggplot(data=subset(
  summer1, locus_tag=="16S"), aes(
    x=date, y=Genus, 
    height=Abundance, 
    fill=Genus)) +
geom_ridgeline(
  alpha=0.94, 
  size=0.1,
  scale=1.15, 
  show.legend=F) +
scale_x_continuous(
  expand = c(0.02,0.02)) +
scale_y_discrete(
  expand = c(0.02,0.02)) +
scale_fill_manual(
    values = pnw_palette("Cascades",8)) +
theme_classic() +
theme(axis.title=element_blank(),
        axis.text.x=element_blank(),
        axis.text.y=element_text(size=10),
        panel.grid = element_blank(),
        axis.ticks = element_blank())

####################################

plot_summer03 <- ENV.bac %>%
  filter(site=="WSC" & month %in% c(
    "Apr","May","Jun","Jul")) %>%
  ggplot() +
  geom_point(aes(
    x=date, y=chl_a), 
    shape=16, size=2.6, color="gray22") +
  scale_y_continuous(
    na.value = 0, 
    limits = c(-0.3,3.2), 
    breaks = c(0,1,2,3),
    position="right") +
  theme_classic() + 
  theme(axis.title.x = element_blank(), 
        axis.title.y = element_text(size=8), 
        axis.text.x = element_blank(),
        axis.text.y = element_text(size=10),
        axis.ticks = element_blank(),
        panel.grid.major = element_line(size=0.4)) 

####################################

plot_summer04 <- ENV.bac %>%
  filter(site=="WSC" & month %in% c(
    "Apr","May","Jun","Jul")) %>%
  ggplot() +
  geom_point(aes(
    x=date, y=NO3_NO2), 
    shape=16, size=2.6, color="gray22") +
  scale_y_continuous(
    limits = c(1,14), 
    breaks = c(1,4,8,12),
    position="right") +
  theme_classic() + 
  theme(axis.title.x = element_blank(), 
        axis.title.y = element_text(size=8), 
        axis.text.x = element_blank(),
        axis.text.y = element_text(size=10),
        axis.ticks = element_blank(),
        panel.grid.major = element_line(size=0.4)) 

####################################

plot_summer05 <- ENV.bac %>%
  filter(site=="WSC" & month %in% c(
    "Apr","May","Jun","Jul")) %>%
  ggplot() +
  geom_point(aes(
    x=date, y=PO4), 
    shape=16, size=2.6, color="gray22") +
  scale_y_continuous(
    limits = c(0,1), 
    breaks = c(0.25,0.5,0.75),
    position="right") +
  theme_classic() + 
  theme(axis.title.x = element_blank(), 
        axis.title.y = element_text(size=8), 
        axis.text.x = element_blank(),
        axis.text.y = element_text(size=10),
        axis.ticks = element_blank(),
        panel.grid.major = element_line(size=0.4)) 

####################################

plot_summer06 <- ENV.bac %>%
  filter(site=="WSC" & month %in% c(
    "Apr7","May","Jun","Jul")) %>%
  ggplot() +
  geom_point(aes(
    x=date, y=O2_conc), 
    shape=16, size=2.6, color="gray22") +
  scale_y_continuous(
    limits = c(300,340), 
    breaks = c(310,320,330),
    position="right") +
  theme_classic() + 
  theme(axis.title.x = element_blank(), 
        axis.title.y = element_text(size=8), 
        axis.text.x = element_blank(),
        axis.text.y = element_text(size=10),
        axis.ticks = element_blank(),
        panel.grid.major = element_line(size=0.4)) 

##################################################################################
   
### EGC ###

summer2 <- sum.genus %>% filter(
  site=="EGC" & month %in% c(
  "Apr","May","Jun","Jul") & 
    !RAS_id %in% c("07_2016_EGC_1","07_2016_EGC_2") &
  Genus %in% c(
    "Colwellia","Aurantivirga",
    "Marinimicrobia_SAR406_uc",
     "Cryomorphaceae_uc","Lentimonas",
     "OCS116_clade_uc","Luteolibacter",
     "Marinomonas","Fragilariopsis",
     "Chaetoceros","Leegaardiella",
     "Phaeocystis","Chromidina",
    "Chytriodinium","Woloszynskia",
    "Gyrodinium")) %>%
  mutate(date = as.Date(
    date, format = "%Y-%m-%d")) %>%
  ungroup %>%
  mutate(Genus = factor(Genus, levels=c(
    "Lentimonas","Marinomonas",
    "Cryomorphaceae_uc","OCS116_clade_uc",
    "Aurantivirga","Luteolibacter",
    "Marinimicrobia_SAR406_uc",
    "Colwellia","Leegaardiella",
    "Fragilariopsis","Chromidina",
    "Chytriodinium","Woloszynskia",
    "Chaetoceros","Phaeocystis",
    "Gyrodinium"))) %>%
  select_if(~all(!is.na(.))) %>%  
  distinct() 

plot_summer07 <- ggplot(data=subset(
  summer2, locus_tag=="18S"), aes(
    x=date, y=Genus, 
    height=Abundance, fill=Genus)) +
  geom_ridgeline(
    alpha=0.94, 
    size=0.1,
    scale=0.28, 
    show.legend=F) +
  scale_x_continuous(
    expand = c(0.02,0.02)) +
  scale_y_discrete(
    expand = c(0.02,0.02),
    position="right") +
  #scale_fill_manual(values=col.taxon) +
  scale_fill_manual(
    values = pnw_palette("Shuksan2",8)) +
  theme_classic() +
  theme(axis.title=element_blank(),
        axis.text.x=element_blank(),
        axis.text.y=element_text(size=10),
        panel.grid = element_blank(),
        axis.ticks = element_blank())

plot_summer08 <- ggplot(data=subset(
  summer2, locus_tag=="16S"), aes(
    x=date, y=Genus, 
    height=Abundance, fill=Genus)) +
  geom_ridgeline(
    alpha=0.94, 
    size=0.1,
    scale=0.38, 
    show.legend=F) +
  scale_x_continuous(
    expand = c(0.02,0.02)) +
  scale_y_discrete(
    expand = c(0.02,0.02),
    position="right") +
  #scale_fill_manual(values=col.taxon) +
  scale_fill_manual(
    values = pnw_palette("Anemone",8)) +
  theme_classic() +
  theme(axis.title=element_blank(),
        axis.text.x=element_blank(),
        axis.text.y=element_text(size=10),
        panel.grid = element_blank(),
        axis.ticks = element_blank())

plot_summer09 <- ENV.euk %>%
  filter(site=="EGC" & month %in% c(
    "Apr","May","Jun","Jul") & 
      !RAS_id %in% c("07_2016_EGC_1","07_2016_EGC_2")) %>%
  ggplot() +
  geom_point(aes(
    x=date, y=chl_a), 
    shape=16, size=2.6, 
    color="gray22") +
  scale_y_continuous(
    limits = c(-0.3,3.2), 
    breaks = c(0,1,2,3)) +
  theme_classic() + 
  theme(axis.title.x = element_blank(), 
        axis.title.y = element_text(size=8), 
        axis.text.x = element_blank(),
        axis.text.y = element_text(size=10),
        axis.ticks = element_blank(),
        panel.grid.major = element_line(size=0.4)) 

####################################

plot_summer10 <- ENV.euk %>%
  filter(site=="EGC" & month %in% c(
    "Apr","May","Jun","Jul") & 
    !RAS_id %in% c("07_2016_EGC_1","07_2016_EGC_2")) %>%
  ggplot() +
  geom_point(aes(
    x=date, y=NO3_NO2), 
    shape=16, size=2.6, color="gray22") +
  scale_y_continuous(
    limits = c(1,14), 
    breaks = c(1,4,8,12)) +
  theme_classic() + 
  theme(axis.title.x = element_blank(), 
        axis.title.y = element_text(size=8), 
        axis.text.x = element_blank(),
        axis.text.y = element_text(size=10),
        axis.ticks = element_blank(),
        panel.grid.major = element_line(size=0.4)) 

####################################

plot_summer11 <- ENV.euk %>%
  filter(site=="EGC" & month %in% c(
    "Apr","May","Jun","Jul") & 
      !RAS_id %in% c("07_2016_EGC_1","07_2016_EGC_2")) %>%
  ggplot() +
  geom_point(aes(
    x=date, y=PO4), 
    shape=16, size=2.6, color="gray22") +
  scale_y_continuous(
    limits = c(0,1), 
    breaks = c(0.25,0.5,0.75)) +
  theme_classic() +
  theme(axis.title.x = element_blank(), 
        axis.title.y = element_text(size=8), 
        axis.text.x = element_blank(),
        axis.text.y = element_text(size=10),
        axis.ticks = element_blank(),
        panel.grid.major = element_line(size=0.4)) 

####################################

plot_summer12 <- ENV.euk %>%
  filter(site=="EGC" & month %in% c(
    "Apr","May","Jun","Jul") & 
      !RAS_id %in% c("07_2016_EGC_1","07_2016_EGC_2")) %>%
  ggplot() +
  scale_y_continuous(
    limits = c(300,340), 
    breaks = c(310,320,330)) +
  geom_point(aes(
    x=date, y=O2_conc), 
    shape=16, size=2.6, color="gray22") +
  theme_classic() + 
  theme(axis.title.x = element_blank(), 
        axis.title.y = element_text(size=8), 
        axis.text.x = element_blank(),
        axis.text.y = element_text(size=10),
        axis.ticks = element_blank(),
        panel.grid.major = element_line(size=0.4)) 

####################################

# export size 7.4 x 9
plot_grid(
  plot_summer01, 
  plot_summer07, 
  plot_summer02, 
  plot_summer08,
  plot_summer03, 
  plot_summer09,
  plot_summer04,
  plot_summer10, 
  plot_summer05,
  plot_summer11, 
  plot_summer06,
  plot_summer12, 
  ncol = 2,
  rel_heights = c(
    2.5, 2.5, 1, 1, 1, 1),
  align = "v",
  axis = "tblr")

# remove temp-data
rm(summer1, summer2)
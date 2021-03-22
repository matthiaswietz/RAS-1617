############################################################################################
   ###  RAS2016-17 -- AMPLICONS -- FIG. 6 ###
############################################################################################

summer1 <- sum.genus %>%
  filter(biome=="WSC" &
           month %in% c("Apr","May","Jun","Jul") & 
           date != "2017-07-22" &
           Genus %in% c(
             "NS2b_marine_group","NS10_marine_group",
             "Aurantivirga","Formosa",
             "Polaribacter","Cyclobacteriaceae_uc",
             "Amylibacter","Dadabacteriales_uc",
             "NS4_marine_group","Phaeocystis",
             "Grammonema","Chytriodinium",
             "Micromonas","Woloszynskia","Gyrodinium",
             "Chaetoceros","Thalassiosira","Odontella",
             "Ephelota","Pelagophyceae_XXX",
             "Pseudo-nitzschia","OCS116_clade_uc")) %>%
  mutate(date = as.Date(
    date, format = "%Y-%m-%d")) %>%
  ungroup %>%
  mutate(Genus = factor(Genus, levels=c(
    "Cyclobacteriaceae_uc","Dadabacteriales_uc",
    "NS4_marine_group","NS10_marine_group",
    "Polaribacter","OCS116_clade_uc",
    "NS2b_marine_group","Formosa",
    "Amylibacter","Aurantivirga",
    "Woloszynskia","Odontella",
    "Ephelota","Chytriodinium",
    "Micromonas","Pelagophyceae_XXX",
    "Chaetoceros","Gyrodinium",
    "Grammonema","Pseudo-nitzschia",
    "Thalassiosira","Phaeocystis"))) %>%
  select_if(~all(!is.na(.))) %>%  
  distinct() 

####################################

plot_summer01 <- ggplot(data=subset(
  summer1, type=="euk"), aes(
    x=date, y=Genus, 
    height=Abundance, fill=Genus)) +
  geom_ridgeline(
    alpha=0.97, 
    scale=0.3, 
    show.legend=F) +
  scale_fill_manual(values=col.taxon) +
  theme_bw() +
  theme(axis.title=element_blank(),
        axis.text.x=element_blank(),
        panel.grid = element_blank(),
        axis.ticks = element_blank())

plot_summer02 <- ggplot(data=subset(
  summer1, type=="bac"), aes(
    x=date, y=Genus, 
    height=Abundance, fill=Genus)) +
  geom_ridgeline(
    alpha=0.97, 
    scale=1.1, 
    show.legend=F) +
  scale_fill_manual(values=col.taxon) +
  theme_bw() +
  theme(axis.title=element_blank(),
        axis.text.x=element_blank(),
        panel.grid = element_blank(),
        axis.ticks = element_blank())

####################################

plot_summer03 <- ENV.bac %>%
  filter(biome=="WSC" & month_full %in% c(
    "Apr_17","May_17","Jun_17","Jul_17")) %>%
  ggplot() +
  geom_point(aes(
    x=date, y=chl_a), 
    shape=16, size=2.6, color="gray22") +
  scale_y_continuous(
    na.value = 0, 
    limits = c(-0.3,3.2), 
    breaks = c(0,1,2,3)) +
  theme_bw() + 
  theme(axis.title.x = element_blank(), 
        axis.title.y = element_text(size=6), 
        axis.text.x = element_blank(),
        axis.text.y = element_text(size=6),
        axis.ticks.x = element_blank(),
        panel.grid.minor = element_blank()) 

####################################

plot_summer04 <- ENV.bac %>%
  filter(biome=="WSC" & month_full %in% c(
    "Apr_17","May_17","Jun_17","Jul_17")) %>%
  ggplot() +
  geom_point(aes(
    x=date, y=NO3_NO2), 
    shape=16, size=2.6, color="gray22") +
  scale_y_continuous(
    limits = c(1,14), 
    breaks = c(1,4,8,12)) +
  theme_bw() + 
  theme(axis.title.x = element_blank(), 
        axis.title.y = element_text(size=6), 
        axis.text.x = element_blank(),
        axis.text.y = element_text(size=6),
        axis.ticks.x = element_blank(),
        panel.grid.minor = element_blank()) 

####################################

plot_summer05 <- ENV.bac %>%
  filter(biome=="WSC" & month_full %in% c(
    "Apr_17","May_17","Jun_17","Jul_17")) %>%
  ggplot() +
  geom_point(aes(
    x=date, y=PO4), 
    shape=16, size=2.6, color="gray22") +
  scale_y_continuous(
    limits = c(0,1), 
    breaks = c(0.25,0.5,0.75)) +
  theme_bw() + 
  theme(axis.title.x = element_blank(), 
        axis.title.y = element_text(size=6), 
        axis.text.x = element_blank(),
        axis.text.y = element_text(size=6),
        axis.ticks.x = element_blank(),
        panel.grid.minor = element_blank()) 

####################################

plot_summer06 <- ENV.bac %>%
  filter(biome=="WSC" & month_full %in% c(
    "Apr_17","May_17","Jun_17","Jul_17")) %>%
  ggplot() +
  geom_point(
    aes(x=date, y=rho), 
    shape=16, size=2.6, color="gray22") +
  scale_y_continuous(
    limits = c(1027.6,1028.5),
    breaks = c(1027.7,1028.0,1028.3)) +
  theme_bw() + 
  theme(axis.title.x=element_blank(), 
        axis.title.y=element_text(size=6), 
        axis.text.x=element_blank(),
        axis.text.y=element_text(size=6),
        axis.ticks.x=element_blank(),
        panel.grid.minor=element_blank()) 

####################################

plot_summer07 <- ENV.bac %>%
  filter(biome=="WSC" & month_full %in% c(
    "Apr_17","May_17","Jun_17","Jul_17")) %>%
  ggplot() +
  geom_point(aes(
    x=date, y=O2_conc), 
    shape=16, size=2.6, color="gray22") +
  scale_y_continuous(
    limits = c(300,340), 
    breaks = c(310,320,330)) +
  theme_bw() + 
  theme(axis.title.x = element_blank(), 
        axis.title.y = element_text(size=6), 
        axis.text.x = element_blank(),
        axis.text.y = element_text(size=6),
        axis.ticks.x = element_blank(),
        panel.grid.minor = element_blank())


##################################################################################
###  SUMMER -- EGC  ###
###################################################################################

summer2 <- sum.genus %>%
  filter(biome=="EGC" & month_full %in% c(
    "Apr_17","May_17","Jun_17","Jul_17") & 
      Genus %in% c(
        "Colwellia","Aurantivirga",
        "Marinimicrobia_SAR406_uc",
        "Cryomorphaceae_uc","Lentimonas",
        "Novel-clade-2_X","Marinomonas",
        "OCS116_clade_uc","Luteolibacter",
        "Dadabacteriales_uc","Chaetoceros",
        "Fragilariopsis","TAGIRI1-lineage_X",
        "Leegaardiella","Nitrincolaceae_uc",
        "Stramenopiles_XXXXX","Phaeocystis",
        "Chromidina","Heterocapsa","Chytriodinium",
        "Woloszynskia","Gyrodinium")) %>%
  mutate(date = as.Date(
    date, format = "%Y-%m-%d")) %>%
  ungroup %>%
  mutate(Genus = factor(Genus, levels=c(
    "Lentimonas","OCS116_clade_uc",
    "Cryomorphaceae_uc","Aurantivirga",
    "Luteolibacter","Nitrincolaceae_uc",
    "Marinomonas","Dadabacteriales_uc",
    "Marinimicrobia_SAR406_uc","Colwellia",
    "Chytriodinium","Novel-clade-2_X",
    "Heterocapsa","Leegaardiella",
    "TAGIRI1-lineage_X","Stramenopiles_XXXXX",
    "Woloszynskia","Chaetoceros",
    "Fragilariopsis","Phaeocystis",
    "Chromidina","Gyrodinium"))) %>%
  select_if(~all(!is.na(.))) %>%  
  distinct() 

####################################

plot_summer08 <- ggplot(
  data=subset(summer2, type=="euk"), 
  aes(x=date, y=Genus, 
      height=Abundance, fill=Genus)) +
  geom_ridgeline(
    alpha=0.97, 
    scale=0.36, 
    show.legend=F) +
  scale_fill_manual(values=col.taxon) +
  theme_bw() +
  theme(axis.title=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks = element_blank(),
        panel.grid = element_blank())

plot_summer09 <- ggplot(data=subset(
  summer2, type=="bac"), aes(
    x=date, y=Genus, 
    height=Abundance, fill=Genus)) +
  geom_ridgeline(
    alpha=0.97, 
    scale=0.42, 
    show.legend=F) +
  scale_fill_manual(values=col.taxon) +
  theme_bw() +
  theme(axis.title=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks = element_blank(),
        panel.grid = element_blank())

####################################

plot_summer10 <- ENV.euk %>%
  filter(biome=="EGC" & month_full %in% c(
    "Apr_17","May_17","Jun_17","Jul_17")) %>%
  ggplot() +
  geom_point(aes(
    x=date, y=chl_a), 
    shape=16, size=2.6, color="gray22") +
  scale_y_continuous(
    limits = c(-0.3,3.2), 
    breaks = c(0,1,2,3)) +
  theme_bw() + 
  theme(axis.title.x=element_blank(), 
        axis.title.y=element_text(size=6), 
        axis.text.x=element_blank(),
        axis.text.y=element_text(size=6),
        axis.ticks.x=element_blank(),
        panel.grid.minor=element_blank()) 

####################################

plot_summer11 <- ENV.euk %>%
  filter(biome=="EGC" & month_full %in% c(
    "Apr_17","May_17","Jun_17","Jul_17")) %>%
  ggplot() +
  geom_point(aes(
    x=date, y=NO3_NO2), 
    shape=16, size=2.6, color="gray22") +
  scale_y_continuous(
    limits = c(1,14), 
    breaks = c(1,4,8,12)) +
  theme_bw() + 
  theme(axis.title.x = element_blank(), 
        axis.title.y = element_text(size=6), 
        axis.text.x = element_blank(),
        axis.text.y = element_text(size=6),
        axis.ticks.x = element_blank(),
        panel.grid.minor = element_blank()) 

####################################

plot_summer12 <- ENV.euk %>%
  filter(biome=="EGC" & month_full %in% c(
    "Apr_17","May_17","Jun_17","Jul_17")) %>%
  ggplot() +
  geom_point(aes(
    x=date, y=PO4), 
    shape=16, size=2.6, color="gray22") +
  scale_y_continuous(
    limits = c(0,1), 
    breaks = c(0.25,0.5,0.75)) +
  theme_bw() + 
  theme(axis.title.x = element_blank(), 
        axis.title.y = element_text(size=6), 
        axis.text.x = element_blank(),
        axis.text.y = element_text(size=6),
        axis.ticks.x = element_blank(),
        panel.grid.minor = element_blank()) 

####################################

plot_summer13 <- ENV.euk %>%
  filter(biome=="EGC" & month_full %in% c(
    "Apr_17","May_17","Jun_17","Jul_17")) %>%
  ggplot() +
  geom_point(
    aes(x=date, y=rho), 
    shape=16, size=2.6, color="gray22") +
  scale_y_continuous(
    limits = c(1027.6,1028.5),
    breaks = c(1027.7,1028.0,1028.3)) +
  theme_bw() + 
  theme(axis.title.x=element_blank(), 
        axis.title.y=element_text(size=6), 
        axis.text.x=element_blank(),
        axis.text.y=element_text(size=6),
        axis.ticks.x=element_blank(),
        panel.grid.minor=element_blank()) 

####################################

plot_summer14 <- ENV.euk %>%
  filter(biome=="EGC" & month_full %in% c(
    "Apr_17","May_17","Jun_17","Jul_17")) %>%
  ggplot() +
  scale_y_continuous(
    limits = c(300,340), 
    breaks = c(310,320,330)) +
  geom_point(aes(
    x=date, y=O2_conc), 
    shape=16, size=2.6, color="gray22") +
  theme_bw() + 
  theme(axis.title.x=element_blank(), 
        axis.title.y=element_text(size=6), 
        axis.text.x = element_blank(),
        axis.text.y = element_text(size=6),
        axis.ticks.x=element_blank(),
        panel.grid.minor=element_blank())

####################################

plot_grid(
  plot_summer01, 
  plot_summer08, 
  plot_summer02, 
  plot_summer09,
  plot_summer03, 
  plot_summer10,
  plot_summer04,
  plot_summer11, 
  plot_summer05,
  plot_summer12, 
  plot_summer06,
  plot_summer13, 
  plot_summer07,
  plot_summer14, 
  ncol = 2,
  rel_heights = c(
    2.6, 2.24, 0.8, 0.8, 0.8, 0.8, 0.8, 0.8),
  align = "v",
  axis = "tblr")

# remove temp-data
rm(summer1, summer2)
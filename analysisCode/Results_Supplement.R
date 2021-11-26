
############################################################################################
   ### RAS 2016-17 -- AMPLICON ANALYSIS ###
   ### Supplementary Figures
############################################################################################

###  Suppl. Fig 1a --  site ENV-PCA  ###

PCA.1 <- ENV.bac %>%
  dplyr::select(-c(
    salinity, chl_a, CO2, pH, dba,
    depth,sig, ice_past, strat,
    clip_id, date, mooring,
    lat, lon, locus_tag)) %>%
  drop_na()

dplyr::select(PCA.1, 
  season, watermass, RAS_id, 
  month, site) -> PCA.2 

PCA.1 <- PCA.1 %>% 
  dplyr::select(-c(
    season, watermass, month, 
    RAS_id, site, ice_cover))

autoplot(
  prcomp(PCA.1, center=T, scale=T),
  data = PCA.2, label=F, 
  loadings=T, loadings.colour="gray44",
  loadings.label=T, loadings.label.size=2.5,
  loadings.label.vjust=1.3,
  loadings.label.colour="black") +
geom_point(aes(
  shape = watermass, 
  color = site,
  size=AW), 
  size=4) +
scale_color_manual(values=col.env) +
scale_shape_manual(values=c(15,17,16)) +
theme_bw() +
theme(
  panel.grid.minor = element_blank(),
  axis.ticks = element_blank(),
  axis.text = element_blank(),
  legend.position = "none")


############################################################################################
   ###  Suppl. Fig 1b  ###
############################################################################################

ENV %>%
  reshape2::melt(id.vars=c(
    "site","date","watermass")) %>%
  filter(variable %in% c(
    "PO4","NO3_NO2","temp","ice_conc")) %>%
  mutate(watermass=factor(
    watermass, levels = c("AW","mix","PW"))) %>%
  mutate(across(value, as.numeric)) %>%
ggplot() +
geom_boxplot(aes(
  x=watermass, y=value, fill=watermass), 
  lwd=0.3, outlier.size=0.5) +
facet_grid(variable~site, scales="free") +
scale_fill_manual(values = c(
  "AW"="#bf5dab",
  "PW"="#80b3ff",
  "mix"="#b0b0b0")) +
theme_bw() +
theme(axis.ticks = element_blank(),
      panel.grid.minor = element_blank(),
      axis.title = element_blank(),
      axis.text.x = element_blank(),
      legend.position = "right")


############################################################################################
   ###  Suppl Fig. 2a ###
############################################################################################

amp_ordinate(ampvis.bac, 
  filter_species = 0.1, 
  type = "NMDS",
  transform = "hellinger", 
  distmeasure = "bray", 
  sample_color_by = "site", 
  sample_colorframe = T,
  sample_label_size = 1.5) +
scale_colour_manual(values=col.env) +
scale_fill_manual(values=col.env) +
geom_point(aes(size=AW)) +
scale_size(range = c(1, 5)) +
theme_bw() +
theme(legend.position = "none",
      axis.text = element_blank(),
      axis.ticks = element_blank(),
      axis.title = element_text(size=9),
      panel.grid.minor = element_blank())

amp_ordinate(ampvis.euk, 
  filter_species = 0.1, 
  type = "NMDS",
  transform = "hellinger", 
  distmeasure = "bray", 
  sample_color_by = "site", 
  sample_colorframe = T,
  sample_label_size = 1.5) +
scale_colour_manual(values=col.env) +
scale_fill_manual(values=col.env) +
geom_point(aes(size=AW)) +
scale_size(range = c(1, 5)) +
theme_bw() +
theme(legend.position = "none",
      axis.text = element_blank(),
      axis.ticks = element_blank(),
      axis.title = element_text(size=9),
      panel.grid.minor = element_blank())


###################################################################################
   ###  Suppl. Fig. 3 ###
###################################################################################

pheatmap(
  as.matrix(dist.euk), 
  cellwidth = 7, 
  cellheight = 7, 
  show_rownames = T, 
  show_colnames = T, 
  fontsize_col=8,
  na_col = "white",
  fontsize=7, 
  border_color = "gray98",
  cluster_rows = T, 
  cluster_cols = T, 
  labels_row = anno$month,
  labels_col = anno$month,
  clustering_method = "complete",
  col = two.colors(
    n=33, start='#defcf6', 
    end="#643e78", 
    middle="#99ffdd"),
  annotation = anno,
  annotation_colors = anno_col,
  angle_col = c("90"))

###################################################################################

dist1 <- as.data.frame(as.matrix(dist.bac)) %>% 
  rownames_to_column() %>% 
  pivot_longer(-rowname) %>%
  dplyr::rename(RAS_id = rowname) %>%
  left_join(dplyr::select(
    ENV.bac, RAS_id, AW), by=c("RAS_id")) %>%
  dplyr::rename(ID1 = name, ID2 = RAS_id) %>%
  filter(!grepl("[EGC]", ID1)) %>% 
  filter(!grepl("[F4]", ID2)) %>%
  distinct(value, .keep_all=T)   

dist2 <- as.data.frame(as.matrix(dist.euk)) %>% 
  rownames_to_column() %>% 
  pivot_longer(-rowname) %>%
  dplyr::rename(RAS_id = rowname) %>%
  left_join(dplyr::select(
    ENV.euk, RAS_id, AW), by=c("RAS_id")) %>%
  dplyr::rename(ID1 = name, ID2 = RAS_id) %>%
  filter(!grepl("[EGC]", ID1)) %>% 
  filter(!grepl("[F4]", ID2)) %>%
  distinct(value, .keep_all=T)  

# AW: BAC -0.4 (e-16) // EUK -0.2 (e-5)
# ice_act: BAC 0.3 (e-16) // EUK 0.2 (e-5) [ice_past lower]
cor.test(
  dist2$value,
  dist2$AW, 
  use = "complete.obs", 
  method = "spearman", 
  exact = F)

# remove temp-data
rm(dist1, dist2)


############################################################################################
   ###  Suppl. Fig. 5  ###
############################################################################################

## FAMILIES -- full overview ##

plot_fam1 <- amp_heatmap(
  ampvis.euk, 
  tax_aggregate = "Family",
  tax_show = c(
    "RAD-C uc","RAD-B-Group-IV",
    "Dino-II-7","Dino-I-1",
    "Opalinopsidae","Pelagomonadaceae",
    "Heterocapsaceae","Pelagophyceae uc",
    "Pirsonia XX","MAST-3I uc",
    "Coscinodiscophyceae","Raphid-pennate",
    "Suessiaceae","Tovelliaceae",
    "Chytriodiniaceae","Ceratiaceae",
    "Strombidiidae_M","Gymnodiniaceae",
    "Mamiellaceae","Araphid-pennate",
    "Phaeocystaceae","Mediophyceae"),
  order_y_by = c(
    "Pirsonia XX","MAST-3I uc",
    "RAD-C uc","Dino-II-7","Dino-I-1",
    "RAD-B-Group-IV","Opalinopsidae",
    "Pelagomonadaceae","Heterocapsaceae",
    "Pelagophyceae uc","Coscinodiscophyceae",
    "Raphid-pennate","Suessiaceae",
    "Tovelliaceae","Chytriodiniaceae",
    "Ceratiaceae","Strombidiidae_M",
    "Gymnodiniaceae","Mamiellaceae",
    "Araphid-pennate","Phaeocystaceae",
    "Mediophyceae"),
  plot_values = F, 
  plot_colorscale = "log10",
  plot_legendbreaks = c(1, 5, 10, 20),
  min_abundance = 0.7, 
  max_abundance = 25,
  normalise = T,
  group_by = "month",
  facet_by = "site",
  color_vector = c( 
    "gray97","seashell2",
    "lightsteelblue2","skyblue4",
    "#1f2252", "#0a0b2b")) +
geom_text(aes(
  label = round(Abundance),
  color = ifelse(
    Abundance < 4.4, "black","white"))) +
scale_color_identity() +
theme(axis.text.x = element_blank(),
      axis.text.y = element_text(
        size = 10, face = "italic"),
      strip.background = element_blank(), 
      strip.text = element_blank(),
      legend.position = "none",
      axis.ticks = element_blank())

plot_fam2 <- amp_heatmap(
  ampvis.bac, 
  tax_aggregate = "Family",
  tax_show = c(
    "Colwelliaceae","Marinimicrobia_SAR406 uc",
    "Woeseiaceae","Pirellulaceae",
    "Nitrospinaceae","Arctic97B-4 uc",
    "Dadabacteriales uc","Thioglobaceae",
    "Nitrosopumilaceae","Gammaproteobacteria uc",
    "Defluviicoccales uc","Magnetospiraceae",
    "AEGEAN-169","Nitrincolaceae",
    "SAR116","OCS116",
    "Rubritaleaceae","Rhodobacteraceae",
    "Cryomorphaceae","SAR11 Clade I",
    "Flavobacteriaceae","SAR86 uc"),
  order_y_by = c(
    "Marinimicrobia_SAR406 uc","Colwelliaceae",
    "Woeseiaceae","Gammaproteobacteria uc",
    "Pirellulaceae","Nitrospinaceae",
    "Rubritaleaceae","Arctic97B-4 uc",
    "Dadabacteriales uc","Nitrosopumilaceae",
    "Defluviicoccales uc","Magnetospiraceae",
    "AEGEAN-169","SAR116",
    "Rhodobacteraceae","OCS116",
    "Cryomorphaceae","Flavobacteriaceae",
    "Nitrincolaceae","Thioglobaceae",
    "SAR86 uc","SAR11 Clade I"),
  plot_values = F, 
  plot_colorscale = "log10",
  plot_legendbreaks = c(1, 5, 10, 20),
  min_abundance = 0.7, 
  max_abundance = 25,
  normalise = T,
  group_by = "month",
  facet_by = "site",
  color_vector = c( 
    "gray97","seashell2",
    "lightsteelblue2","skyblue4",
    "#1f2252", "#0a0b2b")) +
geom_text(aes(
  label = round(Abundance),
  color = ifelse(
    Abundance < 4.4,"black","white"))) +
scale_color_identity() +
theme(axis.text.x = element_blank(),
      axis.text.y = element_text(
        size = 10, face = "italic"),
      strip.background = element_blank(), 
      strip.text = element_blank(),
      legend.position = "none",
      axis.ticks = element_blank())

# export size 11 x 9.5
plot_grid(
  plot_fam1, 
  plot_fam2, 
  ncol = 1,
  rel_heights = c(
    1, 1),
  align = "v",
  axis = "tb") 

####################################

## DIATOMS -- details ##

plot_dia1 <- amp_heatmap(
  amp_subset_taxa(
    ampvis.euk, c("Bacillariophyta"), 
    normalise = T),
  tax_aggregate = "Class",
  normalise = F,
  plot_values = F,
  round = 0,
  plot_values_size = 4,
  min_abundance = 1,
  max_abundance = 60,
  group_by ="season",
  facet_by = "site",
  plot_colorscale = "sqrt",
  color_vector = c( 
    "gray97","seashell2",
    "lightsteelblue2","skyblue4",
    "#1f2252", "#0a0b2b")) +
geom_text(aes(
  label = round(Abundance),
  color = "white")) +
scale_color_identity() +
theme(axis.text.x = element_blank(),
      axis.text.y = element_text(size=10),
      strip.background = element_blank(), 
      strip.text = element_blank(),
      legend.position = "none",
      axis.ticks.x = element_blank())

## Major genera 
plot_dia2 <- amp_heatmap(
  ampvis.euk, 
  tax_aggregate = "Genus",
  tax_show = c(
    "Chaetoceros","Fragilariopsis",
    "Thalassiosira","Pseudo-nitzschia",
    "Phaeocystis","Rhizosolenia",
    "Grammonema","Proboscia","Corethron"),
  normalise = T,
  plot_values = F,
  round = 0,
  plot_values_size = 4,
  min_abundance = 1,
  max_abundance = 40,
  group_by ="season",
  facet_by = "site",
  plot_colorscale = "log10",
  order_y_by = c(
    "Grammonema","Chaetoceros",
    "Phaeocystis","Thalassiosira",
    "Fragilariopsis","Pseudo-nitzschia",
    "Corethron","Proboscia","Rhizosolenia"),
  color_vector = c( 
    "gray97","seashell2",
    "lightsteelblue2","skyblue4",
    "#1f2252", "#0a0b2b")) +
geom_text(aes(
  label = round(Abundance),
  color = ifelse(
    Abundance < 4.4,"black","white"))) +
scale_color_identity() +
theme(axis.text.y = element_text(
        size = 10,face = "italic"),
      axis.text.x = element_text(size = 10),
      axis.ticks = element_blank(),
      legend.position = "none",
      strip.background = element_blank(),
      strip.text = element_blank()) 

# export size 3.5 x 6
plot_grid(
  plot_dia1, 
  plot_dia2, 
  ncol = 1,
  rel_heights = c(
    0.4, 2),
  align = "v",
  axis = "tb") 


############################################################################################
   ###  Suppl. Fig. 6 -- TRANSITIONS  ###
############################################################################################

plot_trans1 <- ENV %>%
  reshape2::melt(id.vars=c(
    "site","season","RAS_id")) %>%
  filter(site=="WSC" & variable %in% c(
   "daylight") & season !="summer") %>%
mutate(across(value, as.numeric)) %>%
ggplot() +
geom_tile(aes(
  x=RAS_id, y=variable, fill=value)) +
scale_fill_viridis(
  option="B", 
  begin=0, end=0.92,
  breaks=c(0,6,12,18,24)) +
theme_bw() +
theme(axis.title = element_blank(),
      axis.text = element_blank(),
      axis.ticks.x = element_blank(),
      panel.grid.minor = element_blank(),
      legend.position = "none")

plot_trans2 <- sum.genus %>% filter(
  site=="WSC" & Genus %in% c(
    "Chaetoceros","Phaeocystis",
    "Chytridiomycetes uc","SAR92",
    "Labyrinthulaceae uc","Aurantivirga",
    "Dino-I-1 uc","Dino-II-7 uc") & 
    season !="summer") %>%
  ungroup %>%
  mutate(Group = case_when(
    Genus %in% c(
      "Dino-I-1 uc",
      "Dino-II-7 uc")~"Dinos",
    Genus %in% c(
     "Phaeocystis",
     "Chaetoceros")~"Photos",
    Genus %in% c(
     "SAR92",
     "Aurantivirga")~"Bacs",
    Genus %in% c(
     "Chytridiomycetes uc",
     "Labyrinthulaceae uc")~"Fungi")) %>%
  mutate(Group = factor(Group, levels = c(
    "Fungi","Dinos","Photos","Bacs"))) %>%
ggplot() + 
geom_line(aes(
  factor(date), Abundance, 
  color=Genus, group=Genus), size=1.4) +
geom_point(aes(
  factor(date), Abundance, 
  group=1, color=Genus), size=2.5) +
scale_x_discrete(
  expand = c(0.02,0.02))+
scale_y_continuous(
  expand = c(0.1,0.1)) +
scale_color_manual(values=col.taxon) +
facet_wrap(
  ~Group, ncol=1, scales="free",
  strip.position="right") +
theme_bw() +
theme(axis.title = element_blank(),
      axis.text.x = element_blank(),
      panel.grid.minor = element_blank(),
      axis.ticks.x = element_blank(),
      legend.position = "none")

####################################

plot_grid(
  plot_trans1, 
  plot_trans2, 
  ncol = 1,
  rel_heights = c(0.3, 3),
  align = "v", 
  axis = "tblr")


############################################################################################
   ###  Suppl Fig 7a -- EGC WINTER LOOP  ###
############################################################################################

sum.genus %>% 
  filter(site=="EGC" & Genus %in% c(
    "Bacillaria","Naviculales",
    "Polarella","Flavobacterium",
    "Chrysophyceae C uc")) %>% 
ggplot() + 
geom_bar(aes(
  x=date, y=Abundance, 
  fill=Genus, group=Genus),
  stat="identity", position="stack") +
geom_line(aes(
  x=date, y=ice_conc/3.4, 
  group=1), size=1.2, color="darkturquoise") +  
geom_point(aes(
  x=date, y=ice_conc/3.4), 
  color="darkturquoise", size=2) +  
scale_y_continuous(
  breaks = c(0,10,20,30), 
  limits = c(0,30),
  expand = c(0.01,0.01),
  sec.axis = sec_axis(~.*3.4)) +
scale_fill_manual(
  values=pnw_palette("Mushroom",6)) +
theme_bw() +
theme(axis.title = element_blank(),
      axis.text.x = element_blank(),
      panel.grid.minor = element_blank(),
      legend.position = "right")


############################################################################################
  ###  Suppl Fig 7b -- EGC taxa  ###
############################################################################################

sum.genus %>%
  filter(site=="EGC" & Genus %in% c(
    "Nitrospina","Magnetospiraceae uc",
    "Aurantivirga","RAD-B-Group-IV uc",
    "Phaeocystis","OCS116 uc","Colwellia",
    "Thalassiosira","Novel-clade-2 uc",
    "Marinimicrobia_SAR406 uc")) %>%
  filter(watermass !="mix") %>%
  filter(ice_cover !="medium") %>%
  group_by(
    Genus, site, locus_tag,
    watermass, ice_cover) %>%
  summarize_at(c("Abundance"), mean) %>%
  mutate(fraction = case_when(
    watermass=="PW" & ice_cover=="high" ~ Abundance * 1, 
    watermass=="AW" & ice_cover=="low" ~ Abundance * -1,
    TRUE ~ NA_real_)) %>%
  mutate(ID = case_when(
    watermass=="PW" ~ "PW",
    watermass=="AW" ~ "AW")) %>%
  mutate(ID = factor(ID, levels = c(
    "low","high","PW","AW"))) %>%
  drop_na() %>%
ggplot() +
geom_bar(aes(
    x=Genus, y=fraction, fill=ID),
    stat="identity", position="stack", 
    lwd=0.6, color="white") +
scale_x_discrete(limits=c(
  "Colwellia","Marinimicrobia_SAR406 uc",
  "RAD-B-Group-IV uc","Magnetospiraceae uc",
  "Nitrospina","OCS116 uc",
  "Aurantivirga","Novel-clade-2 uc",
  "Thalassiosira","Phaeocystis")) +
scale_fill_manual(values=c(
  "PW"="deepskyblue3", 
  "AW"="indianred1")) +
scale_y_continuous(
  limits = c(-7.5, 7.5)) +
geom_hline(
  yintercept=0, size=0.7,
  linetype="dashed", 
  color="gray22") + 
theme_bw() + 
theme(axis.title = element_blank(),
      axis.text.x = element_text(angle=45, h=1),
      panel.grid.minor = element_blank(),
      panel.grid.major = element_line(size=0.3),
      legend.position = "none")


############################################################################################
   ###  Suppl Fig 7 -- Summer details  ###
############################################################################################

summer1 <- sum.genus %>% filter(
  site=="WSC" &
    month %in% c("Apr","May","Jun","Jul") & 
    date != "2017-07-22" &
    Genus %in% c(
      "Cyclobacteriaceae uc",
      "NS10","Polaribacter",
      "OCS116 uc","Aurantivirga",
      "NS2b","Formosa","Amylibacter",
      "Odontella","Micromonas",
      "Chytriodinium","Woloszynskia",
      "Chaetoceros","Thalassiosira",
      "Phaeocystis","Grammonema")) %>%
  mutate(date = as.Date(
    date, format = "%Y-%m-%d")) %>%
  ungroup %>%
  mutate(Genus = factor(Genus, levels=c(
    "Amylibacter","Formosa",
    "NS2b","Aurantivirga",
    "OCS116 uc","NS10", 
    "Polaribacter","Cyclobacteriaceae uc",
    "Grammonema","Phaeocystis",
    "Chaetoceros","Woloszynskia",
    "Chytriodinium","Micromonas",
    "Odontella","Thalassiosira"))) %>%
  select_if(~all(!is.na(.))) %>%  
  distinct() 

####################################

plot_summer11 <- ggplot(data=subset(
  summer1, locus_tag=="16S"), aes(
    x=date, y=Abundance, 
    fill=Genus, group=Genus)) +
geom_area(size=1) +
scale_fill_manual(
  values = rev(pnw_palette("Cascades",8))) +
facet_wrap(
  ~Genus, ncol=1, 
  scale="free", 
  strip.position="right")+
scale_x_date(
  expand=c(0,0)) +
scale_y_continuous(
  expand=c(0.02,0.02),
  n.breaks = 3) +
theme_classic() +
theme(axis.title=element_blank(),
        #strip.placement = "outside",
        axis.text.x=element_blank(),
        strip.background = element_blank(),
        panel.grid = element_blank(),
        axis.ticks = element_blank(),
        legend.position = "none")

plot_summer12 <- ggplot(data=subset(
  summer1, locus_tag=="18S"), aes(
    x=date, y=Abundance, 
    fill=Genus, group=Genus)) +
geom_area(size=1) +
scale_fill_manual(
  values = pnw_palette("Sailboat",8)) +
scale_y_continuous(
  expand=c(0.02,0.02),
  n.breaks = 3) +
facet_wrap(
  ~Genus, ncol=1, 
  scale="free", 
  strip.position="right")+
theme_classic() +
theme(axis.title=element_blank(),
        axis.text.x=element_blank(),
        strip.background = element_blank(),
        panel.grid = element_blank(),
        axis.ticks = element_blank(),
        legend.position = "none")

#############

summer2 <- sum.genus %>%
  filter(site=="EGC" & month %in% c(
    "Apr","May","Jun","Jul") & !RAS_id %in% c(
    "07_2016_EGC_1","07_2016_EGC_2") & Genus %in% c(
    "Colwellia","Marinimicrobia_SAR406 uc",
    "Luteolibacter","Aurantivirga",
    "OCS116 uc","Cryomorphaceae uc",
    "Marinomonas","Lentimonas",
    "Gyrodinium","Phaeocystis",
    "Chaetoceros","Woloszynskia",
    "Chytriodinium","Chromidina",
    "Fragilariopsis","Leegaardiella")) %>%
  mutate(date = as.Date(
    date, format = "%Y-%m-%d")) %>%
  ungroup %>%
  mutate(Genus = factor(Genus, levels=c(
    "Colwellia","Marinimicrobia_SAR406 uc",
    "Luteolibacter","Aurantivirga",
    "OCS116 uc","Cryomorphaceae uc",
    "Marinomonas","Lentimonas",
    "Gyrodinium","Phaeocystis",
    "Chaetoceros","Woloszynskia",
    "Chytriodinium","Chromidina",
    "Fragilariopsis","Leegaardiella"))) %>%
  select_if(~all(!is.na(.))) %>%  
  distinct() 

####################################

plot_summer13 <- ggplot(data=subset(
  summer2, locus_tag=="16S"), aes(
    x=date, y=Abundance, 
    fill=Genus, group=Genus)) +
geom_area(size=1) +
scale_fill_manual(
  values = rev(pnw_palette("Anemone",8))) +
facet_wrap(
  ~Genus, ncol=1, 
  scale="free", 
  strip.position = "right")+
scale_x_date(
  expand=c(0,0)) +
scale_y_continuous(
  expand=c(0.02,0.02),
  n.breaks = 3) +
theme_classic() +
theme(axis.title=element_blank(),
      axis.text.x=element_blank(),
      strip.background = element_blank(),
      panel.grid = element_blank(),
      axis.ticks = element_blank(),
      legend.position = "none")

plot_summer14 <- ggplot(data=subset(
  summer2, locus_tag=="18S"), aes(
    x=date, y=Abundance, 
    fill=Genus, group=Genus)) +
geom_area(size=1) +
scale_fill_manual(
  values = pnw_palette("Shuksan2",8)) +
scale_y_continuous(
  expand=c(0.02,0.02),
  n.breaks = 3) +
facet_wrap(
  ~Genus, ncol=1, 
  scale="free", 
  strip.position = "right")+
theme_classic() +
theme(axis.title=element_blank(),
      axis.text.x=element_blank(),
      strip.background = element_blank(),
      panel.grid = element_blank(),
      axis.ticks = element_blank(),
      legend.position = "none")

####################################

# export size 12 x 6
plot_grid(
  plot_summer12, 
  plot_summer14, 
  plot_summer11, 
  plot_summer13,
  ncol = 2,
  align = "v",
  axis = "tblr")


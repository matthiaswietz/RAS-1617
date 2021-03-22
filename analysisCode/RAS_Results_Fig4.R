
############################################################################################
   ###  RAS 2016-17 -- AMPLICONS -- Fig.4  ###
############################################################################################

## CLASSES

plot_top1 <- amp_heatmap(
  ampvis.euk, 
  tax_aggregate = "Class",
  tax_show = c(
    "Bacillariophyta",
    "Dinophyceae","Syndiniales",
    "Prymnesiophyceae","RAD-C"),
  order_y_by = c(
    "RAD-C","Syndiniales","Dinophyceae",
    "Prymnesiophyceae","Bacillariophyta"),
  plot_values = T, 
  plot_colorscale = "log10",
  plot_legendbreaks = c(1, 5, 10, 20),
  normalise = T,
  group_by = "season",
  facet_by = "biome") +
  geom_tile(aes(Group, Display), fill="white")+
  geom_point(aes(size=Abundance, color=season)) +
  geom_text(aes(label=round(Abundance, digits=0))) +
  scale_size(range = c(0.01, 13))   +
  scale_color_manual(values = col.env) +
  theme(
    axis.text.x = element_blank(),
    axis.text.y = element_text(size=8),
    legend.position = "none",
    axis.ticks.x = element_blank()) 

plot_top2 <- amp_heatmap(
  ampvis.bac, 
  tax_aggregate = "Class",
  tax_show = c(
    "Alphaproteobacteria","Bacteroidia",
    "Gammaproteobacteria","Verrucomicrobiae",
    "Nitrososphaeria"),
  order_y_by = c(
    "Verrucomicrobiae","Nitrososphaeria",
    "Bacteroidia","Gammaproteobacteria",
    "Alphaproteobacteria"),
  plot_values = T, 
  normalise = T,
  group_by = "season",
  facet_by = "biome") +
  geom_tile(aes(Group, Display), fill="white")+
  geom_point(aes(size=Abundance, color=season)) +
  geom_text(aes(label=round(Abundance, digits=0))) +
  scale_size(range = c(0.01, 13)) +
  scale_color_manual(values = col.env) +
  #theme_classic() +
  theme(
    axis.text.x = element_blank(),
    axis.text.y = element_text(size=8),
    legend.position = "none",
    axis.ticks.x = element_blank()) 

####################################

## FAMILIES

plot_top3 <- amp_heatmap(
  ampvis.euk, 
  tax_aggregate = "Family",
  tax_show = c(
   "RAD-C_XX","RAD-B-Group-IV",
   "Phaeocystaceae","Mamiellaceae",
   "Gymnodiniaceae","Raphid-pennate",
   "Araphid-pennate","Pelagophyceae_XX",
   "Dino-Group-I-Clade-1","Dino-Group-II-Clade-7",
   "Strombidiidae_M","Suessiaceae",
   "Opalinopsidae","Chytriodiniaceae",
   "Radial-centric-basal-Coscinodiscophyceae",
   "Tovelliaceae","Polar-centric-Mediophyceae",
   "MAST-3I_X","Ceratiaceae",
   "Heterocapsaceae","Pelagomonadaceae",
   "Pirsonia_Clade_XX"),
 order_y_by = c(
   "Pirsonia_Clade_XX","MAST-3I_X",
   "RAD-B-Group-IV","RAD-C_XX",
   "Dino-Group-II-Clade-7",
   "Dino-Group-I-Clade-1",
   "Opalinopsidae","Ceratiaceae",
   "Pelagomonadaceae","Heterocapsaceae",
   "Radial-centric-basal-Coscinodiscophyceae",
   "Suessiaceae","Pelagophyceae_XX",
   "Tovelliaceae","Chytriodiniaceae",
   "Strombidiidae_M","Gymnodiniaceae",
   "Mamiellaceae","Araphid-pennate",
   "Phaeocystaceae","Polar-centric-Mediophyceae",
   "Raphid-pennate"),
 plot_values = F, 
 plot_colorscale = "log10",
 plot_legendbreaks = c(1, 5, 10, 20),
 min_abundance = 0.6, 
 max_abundance = 32,
 normalise = T,
 group_by = "season",
 facet_by = "biome",
 color_vector = c( 
   "white","seashell2",
   "lightsteelblue2","skyblue3",
   "#1f2252","#0a0b2b")) +
 theme(
   axis.text.x = element_blank(),
   axis.text.y = element_text(size=8),
   strip.background = element_blank(), 
   strip.text = element_blank(),
   legend.position = "bottom",
   axis.ticks.x = element_blank()) 

plot_top4 <- amp_heatmap(
  ampvis.bac, 
  tax_aggregate = "Family",
  tax_show = c(
    "Colwelliaceae","Marinimicrobia_SAR406_uc",
    "Woeseiaceae","Pirellulaceae",
    "Nitrospinaceae","Arctic97B-4_uc",
    "Dadabacteriales_uc","Thioglobaceae",
    "Nitrosopumilaceae","Gammaproteobacteria_uc",
    "Defluviicoccales_uc","Magnetospiraceae",
    "AEGEAN-169","Nitrincolaceae",
    "SAR116_clade","OCS116_clade",
    "Rubritaleaceae","Rhodobacteraceae",
    "Cryomorphaceae","Flavobacteriaceae",
    "SAR11_Clade_I","SAR86_clade_uc"),
  order_y_by = c(
    "Marinimicrobia_SAR406_uc","Colwelliaceae",
    "Woeseiaceae","Gammaproteobacteria_uc",
    "Pirellulaceae","Nitrospinaceae",
    "Rubritaleaceae","Arctic97B-4_uc",
    "Dadabacteriales_uc","Nitrosopumilaceae",
    "Defluviicoccales_uc","Magnetospiraceae",
    "AEGEAN-169","SAR116_clade",
    "Rhodobacteraceae","OCS116_clade",
    "Cryomorphaceae","Flavobacteriaceae",
    "Nitrincolaceae","Thioglobaceae",
    "SAR86_clade_uc","SAR11_Clade_I"),
  plot_values = F, 
  plot_colorscale = "log10",
  plot_legendbreaks = c(1, 5, 10, 20),
  min_abundance = 0.6, 
  max_abundance = 32,
  normalise = T,
  group_by = "season",
  facet_by = "biome",
  color_vector = c(
    "white","seashell2",
    "lightsteelblue2","skyblue3",
    "#1f2252","#0a0b2b")) +
  theme(
    axis.text.x = element_blank(),
    axis.text.y = element_text(size=8),
    strip.background = element_blank(), 
    strip.text = element_blank(),
    legend.position = "bottom",
    axis.ticks.x = element_blank())

##############

plot_grid(
  plot_top1,
  plot_top2,
  plot_top3,
  plot_top4,
  ncol = 2,
  rel_heights = c(0.9, 1.8),
  align = "v",
  axis = "tblr")


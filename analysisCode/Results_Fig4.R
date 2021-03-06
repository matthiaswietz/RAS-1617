
############################################################################################
  ###  RAS2016-17 -- AMPLICON ANALYSIS  ###
  ###  Fig. 4 -- TOP FAMILIES
############################################################################################

plot_top1 <- amp_heatmap(
  ampvis.euk, 
  tax_aggregate = "Family",
  tax_show = c(
    "RAD-C uc","Dino-II-7","Dino-I-1",
    "Opalinopsidae","Pelagomonadaceae",
    "Coscinodiscophyceae","Raphid-pennate",
    "Suessiaceae","Tovelliaceae",
    "Chytriodiniaceae","Ceratiaceae",
    "Strombidiidae_M","Gymnodiniaceae",
    "Mamiellaceae","Araphid-pennate",
    "Phaeocystaceae","Mediophyceae"),
 order_y_by = c(
   "RAD-C uc","Dino-II-7","Dino-I-1",
   "Opalinopsidae","Pelagomonadaceae",
   "Coscinodiscophyceae","Raphid-pennate",
   "Suessiaceae","Tovelliaceae",
   "Chytriodiniaceae","Ceratiaceae",
   "Strombidiidae_M","Gymnodiniaceae",
   "Mamiellaceae","Araphid-pennate",
   "Phaeocystaceae","Mediophyceae"),
 plot_values = F, 
 plot_colorscale = "log10",
 plot_legendbreaks = c(1, 5, 10, 20),
 min_abundance = 0.7, 
 max_abundance = 30,
 normalise = T,
 group_by = "season",
 facet_by = "site",
 color_vector = c( 
   "gray97","gray76",
   "lightsteelblue2","skyblue4",
   "#4d034f","#6b0101")) +
 theme(
   axis.text.x = element_blank(),
   axis.text.y = element_text(size=8),
   strip.background = element_blank(), 
   strip.text = element_blank(),
   legend.position = "bottom",
   axis.ticks.x = element_blank()) 

plot_top2 <- amp_heatmap(
  ampvis.bac, 
  tax_aggregate = "Family",
  tax_show = c(
    "Colwelliaceae","Marinimicrobia_SAR406 uc",
    "Woeseiaceae","Pirellulaceae",
    "Nitrospinaceae","Arctic97B-4 uc",
    "Nitrosopumilaceae","Gammaproteobacteria uc",
    "Defluviicoccales uc","Magnetospiraceae",
    "AEGEAN-169","Nitrincolaceae",
    "SAR116","OCS116",
    "Rubritaleaceae","Rhodobacteraceae",
    "Flavobacteriaceae"),
  order_y_by = c(
    "Marinimicrobia_SAR406 uc","Colwelliaceae",
    "Woeseiaceae","Gammaproteobacteria uc",
    "Pirellulaceae","Nitrospinaceae",
    "Rubritaleaceae","Arctic97B-4 uc",
    "Nitrosopumilaceae","Defluviicoccales uc",
    "Magnetospiraceae","AEGEAN-169",
    "SAR116","Nitrincolaceae",
    "Rhodobacteraceae","OCS116",
    "Flavobacteriaceae"),
  plot_values = F, 
  plot_colorscale = "log10",
  plot_legendbreaks = c(1, 5, 10, 20),
  min_abundance = 0.7, 
  max_abundance = 30,
  normalise = T,
  group_by = "season",
  facet_by = "site",
  color_vector = c( 
    "gray97","gray76",
    "lightsteelblue2","skyblue4",
    "#4d034f","#6b0101")) +
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
  ncol = 2,
  rel_heights = c(0.9, 1.8),
  align = "v",
  axis = "tblr")


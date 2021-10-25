############################################################################################
   ###  RAS 2016-17 -- AMPLICONS ANALYSIS
   ###  Fig. 6 -- Autumn/winter 
############################################################################################

# NO3 + stratification
plot_night1 <- avg.month %>%
  filter(season %in% c(
    "autumn","winter")) %>%
  dplyr::select(c(
    "month","site",
    "NO3_NO2","SiO4","strat")) %>%
  distinct() %>%
ggplot() +
geom_line(aes(
  x=month, y=NO3_NO2, group=1), 
  size=1, color="gray22", linetype="solid") +  
geom_point(aes(
  x=month, y=NO3_NO2, group=1), 
  shape=15, size=2.4, color="gray22") +  
geom_line(aes(
  x=month, y=SiO4, group=1), 
  size=1, color="gray22", linetype="solid") +  
geom_point(aes(
  x=month, y=SiO4, group=1), 
  shape=17, size=2.4, color="gray22") +  
geom_line(aes(
  x=month, y=strat, group=1), 
  size=1.5, color="dodgerblue", linetype="dotted") +  
geom_point(aes(
  x=month, y=strat, group=1), 
  size=2.4, color="dodgerblue") + 
scale_y_continuous(
  breaks = c(4,8,12), 
  limits = c(0,15),
  expand = c(0.01,0.01)) +
facet_wrap(~site, scales="free") +
theme_bw() +
theme(panel.grid.minor = element_blank(),
      axis.title = element_blank(),
      axis.text = element_blank(),
      legend.position = "none")

# Taxa as barplot
plot_night2 <- avg.month %>%
  filter(Genus %in% c(
    "Proboscia","Fragilariopsis",
    "Dadabacteriales_uc",
    "Pseudo-nitzschia","Rhizosolenia",
    "Corethron","Tripos","RAD-C_uc",
    "Dino-I-1_uc","Bacillaria",
    "Dino-II-7_uc","Naviculales",
    "SAR116_clade_uc","Cand_Puniceispirillum",
    "Amylibacter","Planktomarina","Polarella",
    "Ascidiaceihabitans","Luteolibacter",
    "Flavobacterium","Chrysophyceae_C_uc",
    "Cand_Nitrosopumilus","Magnetospira",
    "Magnetospiraceae_uc","Chaetoceros",
    "Gammaproteobacteria_uc","Colwellia",
    "Arctic97B-4_uc","Marinimicrobia_SAR406_uc") & 
      season %in% c(
        "autumn","winter"))  %>% 
  mutate(Genus2 = case_when(
    Genus %in% c(
      "Dino-I-1_uc",
      "Dino-II-7_uc")~'Syndiniales',
    Genus %in% c(
      "Chrysophyceae_C_X",
      "Polarella","Bacillaria",
      "Naviculales") ~ "icetaxa",
    Genus %in% c(
      "Cand_Puniceispirillum",
      "SAR116_clade_uc")~"SAR116",
    Genus %in% c(
      "Magnetospiraceae_uc",
      "Magnetospira")~"Magnetospiraceae",
    TRUE ~ as.character(.$Genus))) %>%
  filter(Abundance > 1.45) %>%
  select_if(~all(!is.na(.))) %>%  
  distinct() %>%
ggplot() + 
geom_bar(aes(
  x=month, y=Abundance, 
  group=Genus2, fill=Genus2), 
  size=1.3, stat="identity") +
scale_fill_manual(values=col.gradient) +
scale_y_continuous(
  expand = c(0.01,0.01)) +
facet_grid(locus_tag~site, scales="free") +
theme_bw() +
theme(panel.grid.minor = element_blank(),
      axis.title = element_blank(),
      axis.text.x = element_blank(),
      axis.ticks = element_blank(),
      strip.background = element_blank(), 
      strip.text = element_blank(),
      legend.position = "right",
      legend.text = element_text(size=7),
      legend.title = element_blank())

# PW + pH
plot_night3 <- avg.month %>%
  filter(season %in% c(
    "autumn","winter")) %>%
  dplyr::select(c(
    "month","site",
    "PW","pH")) %>%
  distinct() %>%
ggplot() +
geom_line(aes(
  month, PW, group=1), 
  size=1, colour="gray11", 
  linetype="solid") + 
geom_point(aes(
  month, PW, 
  group=PW, colour=pH, size=pH)) + 
scale_colour_gradientn(
  colours = c("gray11","cyan"),
  breaks=c(8.16,8.12,8.08)) +
scale_y_continuous(
  expand = c(0.12,0.08)) +
scale_size(
  breaks=c(8.16,8.12,8.08),
  range=c(1,5)) +
facet_grid(~site, scales="free") +
theme_bw() +
theme(panel.grid.minor = element_blank(),
      axis.title = element_blank(),
      axis.text.x = element_blank(),
      axis.ticks = element_blank(),
      strip.background = element_blank(), 
      strip.text = element_blank(),
      legend.position = "bottom",
      legend.direction = "horizontal",
      legend.text = element_text(size=7))

####################################

plot_grid(
  plot_night1, 
  plot_night2, 
  plot_night3,
  ncol=1,
  rel_heights = c(
    0.2, 0.3, 0.22),
  align="v",
  axis = "tblr")

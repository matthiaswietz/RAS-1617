
############################################################################################
   ###  RAS2016-17 - AMPLICON ANALYSIS  ###
   ###  Fig. 1B  -- ANNUAL ENV COMPARISON
############################################################################################

ENV %>%
  reshape2::melt(id.vars=c(
    "site","date")) %>%
  filter(variable %in% c(
    "temp","O2_conc","daylight",
    "AW","NO3_NO2","ice_act")) %>%
  mutate(across(value, as.numeric)) %>%
ggplot() +
aes(date, value, 
    group=variable, color=variable) +
geom_line(size=1.8, na.rm=T) + 
geom_point(size=1.2, na.rm=T) + 
scale_color_manual(
  values = pnw_palette("Starfish",6)) +
scale_x_date(
  expand = c(0.02,0.02)) +
  scale_y_continuous(
  expand = c(0.2,0.1)) +
facet_grid(variable~site, scales="free") +
theme_bw() +
theme(panel.grid.minor = element_blank(),
      axis.title = element_blank(),
      axis.ticks = element_blank(),
      axis.text.x = element_text(
        angle=90, hjust=0.5, vjust=0.5),
      axis.text = element_text(size=10),
      legend.position = "none")

# export size 7.5 x 9
# aesthetic corrections (e.g. connect dots) in Inkscape


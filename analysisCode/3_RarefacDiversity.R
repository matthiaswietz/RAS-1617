
############################################################################################
    ###  RAS2016-17 - AMPLICON ANALYSIS  ###
    ###  Calculate rarefaction and diversity indices
############################################################################################

library(iNEXT)
library(olsrr)
library(ggplot2)
library(cowplot)


############################################################################################
   ###  RAREFACTION AND COVERAGE  ###
############################################################################################

iNEXT.bac <- otu_table(
  ASV.bac, taxa_are_rows=F)
iNEXT.euk <- otu_table(
  ASV.euk, taxa_are_rows=F)

iNEXT.bac <- iNEXT(
  as.data.frame(otu_table(iNEXT.bac)), q=c(0),
  datatype="abundance", conf = 0.95, nboot = 100)
iNEXT.euk <- iNEXT(
  as.data.frame(otu_table(iNEXT.euk)), q=c(0),
  datatype="abundance", conf = 0.95, nboot = 100)

###################################

rarefac.bac <- fortify(iNEXT.bac, type=1)
rarefac.point.bac <- rarefac.bac[which(
  rarefac.bac$method == "observed"),]
rarefac.line.bac <- rarefac.bac[which(
  rarefac.bac$method != "observed"),]
rarefac.line.bac$method <- factor(rarefac.line.bac$method,
    c("interpolated","extrapolated"),
    c("interpolation","extrapolation"))

rarefac.euk <- fortify(iNEXT.euk, type=1)
rarefac.point.euk <- rarefac.euk[which(
  rarefac.euk$method == "observed"),]
rarefac.line.euk <- rarefac.euk[which(
  rarefac.euk$method != "observed"),]
rarefac.line.euk$method <- factor(rarefac.line.euk$method,
    c("interpolated","extrapolated"),
    c("interpolation","extrapolation"))

###################################

cover.bac <- fortify(iNEXT.bac, type=2)
cover.euk <- fortify(iNEXT.euk, type=2)

cover.point.bac <- cover.bac [which(
  cover.bac$method == "observed"),]
cover.line.bac <- cover.bac [which(
  cover.bac$method != "observed"),]
cover.line.bac$method <- factor(cover.line.bac$method,
    c("interpolated","extrapolated"),
    c("interpolation","extrapolation"))

cover.point.euk <- cover.euk [which(
  cover.euk$method == "observed"),]
cover.line.euk <- cover.euk [which(
  cover.euk$method != "observed"),]
cover.line.euk$method <- factor(cover.line.euk$method,
    c("interpolated","extrapolated"),
    c("interpolation","extrapolation"))

coverage.bac <- ggplot(cover.bac, 
  aes(x=x, y=y, colour=site))+ 
  geom_line(aes(linetype = method), 
  lwd = 0.5, data = cover.line.bac) +
  scale_colour_discrete(guide = F) +
  scale_x_continuous(
    limits = c(0,1e+5)) +
  scale_y_continuous(
    breaks = seq(0.9,1,0.05), 
    limits = c(0.9,1)) +
  labs(x="Sample size", 
       y="Sample coverage") +
  theme_bw(base_size = 12) + 
  theme(legend.position="bottom")

coverage.euk <- ggplot(cover.euk, 
  aes(x=x, y=y, colour=site))+ 
  geom_line(aes(linetype = method), 
  lwd = 0.5, data = cover.line.euk) +
  scale_colour_discrete(guide = F) +
  scale_x_continuous(
    limits = c(0,1e+5)) +
  scale_y_continuous(
    breaks = seq(0.9,1,0.05), 
    limits = c(0.9,1)) +
  labs(x="Sample size", 
       y="Sample coverage") +
  theme_bw(base_size = 12) + 
  theme(legend.position="bottom")

###################################

richness.bac <- ggplot(
  rarefac.bac, aes(x=x, y=y, colour=site)) +
  geom_line(
    aes(linetype = method), 
    lwd = 0.5, data = rarefac.line.bac) +
scale_colour_discrete(guide="none") +
scale_x_continuous(limits = c(0,1e+5)) +
labs(x="Sample size", 
     y="Species richness") +
theme_bw() + 
theme(legend.position="bottom")

richness.euk <- ggplot(
  rarefac.euk, aes(x=x, y=y, colour=site)) +
geom_line(
    aes(linetype = method), 
    lwd = 0.5, data = rarefac.line.euk)+
scale_colour_discrete(guide="none") +
scale_x_continuous(limits = c(0,1e+5)) +
labs(x="Sample size", 
     y="Species richness") +
theme_bw() + 
  theme(legend.position="bottom")


############################################################################################
   ###  ALPHA-DIVERSITY SUMMARIES ###
############################################################################################

# Combine BAC data
richness <- iNEXT.bac$AsyEst[
  iNEXT.bac$AsyEst$Diversity=="Species richness",] %>%
  arrange(Site)
simpson <- iNEXT.bac$AsyEst[
  iNEXT.bac$AsyEst$Diversity=="Simpson diversity",] %>%
  arrange(Site)

# Calculate Shannon Index & evenness
shannon <- estimate_richness(
  pseq.bac.abs, measures=c("Shannon")) %>%
  rownames_to_column("RAS_id") %>%
  arrange(RAS_id)

# compile, and calculate evenness
div.bac <- data.frame(
  RAS_id = iNEXT.bac$DataInfo$site,
  locus_tag = ENV.bac$locus_tag,
  richness = richness$Observed,
  simpson = simpson$Observed,
  shannon = shannon$Shannon) %>%
  mutate(evenness = shannon/log(richness)) 

###################################

# Combine EUK data
richness <- iNEXT.euk$AsyEst[
  iNEXT.euk$AsyEst$Diversity=="Species richness",] %>%
  arrange(Site)
simpson <- iNEXT.euk$AsyEst[
  iNEXT.euk$AsyEst$Diversity=="Simpson diversity",] %>%
  arrange(Site)

# Calculate Shannon Index & evenness
shannon <- estimate_richness(
  pseq.euk.abs, measures=c("Shannon")) %>%
  rownames_to_column("RAS_id") %>%
  arrange(RAS_id)

# compile, and calculate evenness
div.euk <- data.frame(
  RAS_id = iNEXT.euk$DataInfo$site,
  locus_tag = ENV.euk$locus_tag,
  richness = richness$Observed,
  simpson = simpson$Observed,
  shannon = shannon$Shannon) %>%
  mutate(evenness = shannon/log(richness)) 

###################################

# merge BAC/EUK data
div.all <- rbind(
  div.bac, div.euk)

# Plot curves
plot_grid(
  richness.bac, 
  coverage.bac, 
  richness.euk, 
  coverage.euk,  
labels = c("BAC","BAC","EUK","EUK"),
label_x = 0.1,
label_size = 9,
ncol = 2)


############################################################################################
   ###  DISTANCES -- YEAR-ROUND  ### 
############################################################################################

## Distances for PERMANOVA + BetaDiv ##

dist.bac <- merge_phyloseq(
  hel.WSC.bac, hel.EGC.bac) %>%
  phyloseq::distance(
    dist.bac, method = "jsd")
dist.bac.WSC <- phyloseq::distance(
  hel.WSC.bac, method = "jsd")
dist.bac.EGC <- phyloseq::distance(
  hel.EGC.bac, method = "jsd")

dist.euk <- merge_phyloseq(
  hel.WSC.euk, hel.EGC.euk) %>%
  phyloseq::distance(
    dist.euk, method = "jsd")
dist.euk.WSC <- phyloseq::distance(
  hel.WSC.euk, method = "jsd")
dist.euk.EGC <- phyloseq::distance(
  hel.EGC.euk, method = "jsd")

###################################

d.bac.WSC <- as.matrix(dist.bac.WSC) 
d.bac.WSC <- as.data.frame(1 - d.bac.WSC) %>%
  dplyr::select(matches("08_2016_F4_1")) %>%
  rownames_to_column("RAS_id") %>%
  rename("dist"="08_2016_F4_1") 

d.euk.WSC <- as.matrix(dist.euk.WSC) 
d.euk.WSC <- as.data.frame(1 - d.euk.WSC) %>%
  dplyr::select(matches("08_2016_F4_1")) %>%
  rownames_to_column("RAS_id") %>%
  rename("dist"="08_2016_F4_1") 

################### 

d.bac.EGC <- as.matrix(dist.bac.EGC) 
d.bac.EGC <- as.data.frame(1 - d.bac.EGC) %>%
  dplyr::select(matches("07_2016_EGC_1")) %>%
  rownames_to_column("RAS_id") %>%
  rename("dist"="07_2016_EGC_1") 

d.euk.EGC <- as.matrix(dist.euk.EGC) 
d.euk.EGC <- as.data.frame(1 - d.euk.EGC) %>%
  dplyr::select(matches("07_2016_EGC_1")) %>%
  rownames_to_column("RAS_id") %>%
  rename("dist"="07_2016_EGC_1") 


############################################################################################
   ###  MERGE and PLOT  ###
############################################################################################

div.bac <- rbind(
  d.bac.WSC, d.bac.EGC) %>%
  left_join(div.bac, by="RAS_id") %>%
  mutate(locus_tag="16S")

div.euk <- rbind(
  d.euk.WSC, d.euk.EGC) %>%
  left_join(div.euk, by="RAS_id") %>%
  mutate(locus_tag="18S")

# Combine everything
div.all <- rbind(div.bac, div.euk) %>%
  inner_join(ENV, by = c("RAS_id","locus_tag")) 

# Calculate averages
div.avg <- div.all %>%
  group_by(
    month, locus_tag, site) %>%
  summarize_at(c(
    "simpson","richness",
    "shannon","evenness",
    "dist","daylight"), mean) %>% 
  mutate_if(is.numeric, round, 2) %>%
  ungroup 


###################################

# remove temp-data
rm(list = ls(
  pattern = "cover.point.*|cover.line.*|
   rarefac.point.*|rarefac.line.*"))


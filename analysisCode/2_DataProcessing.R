
############################################################################################
     ###  RAS2016-17 - AMPLICON ANALYSIS  ###
############################################################################################

# This script: Prepare phyloseq & ampvis objects
# Data transformations
# Calculate taxon summaries

library(phyloseq)
library(ampvis2)
library(vegan)
library(dplyr)


############################################################################################
   ###  PHYLOSEQ  ###
############################################################################################

## BACTERIA
asv = otu_table(ASV.bac, taxa_are_rows=T)
tax = tax_table(TAX.bac)
rownames(tax) <- rownames(asv)

pseq.bac.abs <- phyloseq(
  otu_table(asv, taxa_are_rows = F), 
  sample_data(ENV.bac), 
  tax_table(tax)) %>%
  filter_taxa(function(x) 
    sum(x>3) > (0.03*length(x)), T)

# Fix tax-IDs; remove Cutibacteria contamination
colnames(pseq.bac.abs@tax_table) <- c(
  "Kingdom","Phylum","Class","Order",
  "Family","Genus","Species")
pseq.bac.abs <- subset_taxa(pseq.bac.abs, 
  Order != "Propionibacteriales")

####################################

## EUKARYOTES
asv = otu_table(ASV.euk, taxa_are_rows=T)
tax = tax_table(TAX.euk)
rownames(tax) <- rownames(asv)

pseq.euk.abs <- phyloseq(
  otu_table(asv, taxa_are_rows=F), 
  sample_data(ENV.euk), 
  tax_table(tax)) %>%
  filter_taxa(function(x) 
    sum(x>3) > (0.03*length(x)), T)

# Fix tax-IDs
colnames(pseq.euk.abs@tax_table) <- c(
 "Kingdom","Phylum","Class",
 "Order","Family","Genus","Species")

####################################

## transform to rel. abundance

pseq.bac.rel = transform_sample_counts(
  pseq.bac.abs, function(x) x / sum(x) * 100) 
pseq.euk.rel = transform_sample_counts(
  pseq.euk.abs, function(x) x / sum(x) * 100) 


############################################################################################
   ###  AMPVIS LOAD  ###
############################################################################################

ampvis.bac <- data.frame(
  OTU = rownames(
    phyloseq::otu_table(pseq.bac.abs)@.Data),
  phyloseq::otu_table(pseq.bac.abs)@.Data,
  phyloseq::tax_table(pseq.bac.abs)@.Data,
  check.names = F)
bac.env <- data.frame(
  phyloseq::sample_data(pseq.bac.abs), 
  check.names = F) %>%
  mutate(Sample=rownames(.), .before = mooring)
ampvis.bac <- amp_load(ampvis.bac, bac.env)

ampvis.euk <- data.frame(
  OTU = rownames(
    phyloseq::otu_table(pseq.euk.abs)@.Data),
  phyloseq::otu_table(pseq.euk.abs)@.Data,
  phyloseq::tax_table(pseq.euk.abs)@.Data,
  check.names = F)
euk.env <- data.frame(
  phyloseq::sample_data(pseq.euk.abs), 
  check.names = F) %>%
  mutate(Sample=rownames(.), .before = mooring)
ampvis.euk <- amp_load(ampvis.euk, euk.env)


############################################################################################
   ###  Hellinger-transform + NMDS  ###
############################################################################################

hel.WSC.bac = subset_samples(
  pseq.bac.rel, site %in% c("WSC")) 
otu_table(hel.WSC.bac) = otu_table(decostand(
  otu_table(hel.WSC.bac), method="hellinger"), 
  taxa_are_rows=T)

hel.WSC.euk = subset_samples(
  pseq.euk.rel, site %in% c("WSC")) 
otu_table(hel.WSC.euk) = otu_table(decostand(
  otu_table(hel.WSC.euk), method="hellinger"), 
  taxa_are_rows=T)

hel.EGC.bac = subset_samples(
  pseq.bac.rel, site %in% c("EGC")) 
otu_table(hel.EGC.bac) = otu_table(decostand(
  otu_table(hel.EGC.bac), method="hellinger"), 
  taxa_are_rows=T)

hel.EGC.euk = subset_samples(
  pseq.euk.rel, site %in% c("EGC")) 
otu_table(hel.EGC.euk) = otu_table(decostand(
  otu_table(hel.EGC.euk), method="hellinger"), 
  taxa_are_rows=T)


############################################################################################
   ###  SUM ASVs BY TAXRANK  ###
############################################################################################

## Taxranks Phylum >> Genus ##

sum.phylum <- psmelt(pseq.bac.rel) %>% as_tibble %>% 
  rbind(psmelt(pseq.euk.rel)) %>%
  group_by(Phylum, locus_tag, site, date) %>%
  summarize_at(c("Abundance"), sum) %>%
  left_join(ENV, by = c(
    "date","locus_tag","site"), na.rm=F) %>%
  dplyr::select(-c(
    "clip_id","mooring",
    "lat","lon"))

sum.class <- psmelt(pseq.bac.rel) %>% as_tibble %>% 
  rbind(psmelt(pseq.euk.rel)) %>%
  group_by(Class, locus_tag, site, date) %>%
  summarize_at(c("Abundance"), sum) %>%
  left_join(ENV, by = c(
    "date","locus_tag","site"), na.rm=F) %>%
  dplyr::select(-c(
    "clip_id","mooring",
    "lat","lon"))

sum.order <- psmelt(pseq.bac.rel) %>% as_tibble %>% 
  rbind(psmelt(pseq.euk.rel)) %>%
  group_by(Order, locus_tag, site, date) %>%
  summarize_at(c("Abundance"), sum) %>%
  left_join(ENV, by = c(
    "date","locus_tag","site"), na.rm=F) %>%
  dplyr::select(-c(
    "clip_id","mooring",
    "lat","lon"))

sum.family <- psmelt(pseq.bac.rel) %>% as_tibble %>% 
  rbind(psmelt(pseq.euk.rel)) %>%
  group_by(Family, locus_tag, site, date) %>%
  summarize_at(c("Abundance"), sum) %>%
  left_join(ENV, by = c(
    "date","locus_tag","site"), na.rm=F) %>%
  dplyr::select(-c(
   "clip_id","mooring",
   "lat","lon"))

sum.genus <- psmelt(pseq.bac.rel) %>% as_tibble %>% 
  rbind(psmelt(pseq.euk.rel)) %>%
  group_by(Genus, locus_tag, site, date) %>%
  summarize_at(c("Abundance"), sum) %>%
  left_join(ENV, by = c(
    "date","locus_tag","site"), na.rm=F) %>%
  dplyr::select(-c(
   "clip_id","mooring",
   "lat","lon"))

####################################

# Combine everything
sum.all <- rbindlist(list(
  sum.phylum, 
  sum.class,
  sum.order,
  sum.family, 
  sum.genus), 
  fill=T) %>%
  mutate(date=as.Date(
    date, format = "%Y-%m-%d"))   

avg.month <- sum.all %>%
  group_by(
    Phylum, Class, Order, Family, Genus, 
    locus_tag, site, month, season) %>%
  summarize_at(c(
    "Abundance","NO3_NO2","NO2","PO4","SiO4",
    "temp","chl_a","O2_conc","O2_sat","pH",
    "daylight","rho","ice_conc","ice_past",
    "CO2","PW","AW","strat"), 
    mean, na.rm=T) %>% 
  mutate_if(is.numeric, round, 2) %>% 
  ungroup 


############################################################################################

# remove temporary datasets
rm(asv, tax, bac.env, euk.env)

############################################################################################

## continued in Rscript *RarefacDiversity*
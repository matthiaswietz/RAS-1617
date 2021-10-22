
############################################################################################
   ###  RAS2016-17 - AMPLICON ANALYSIS  ###  
   ###  Fig 5 -- PLS  
############################################################################################

library(mixOmics)

pls.WSC.tax <- sum.family %>%
  filter(site=="WSC") %>% 
  pivot_wider(
    data = .,
    id_cols = RAS_id, 
    names_from = Family, 
    values_from = c("Abundance")) 

pls.WSC.env <- ENV.bac %>%
  filter(site=="WSC") %>%
  dplyr::select_if(., is.numeric) %>%
  dplyr::select(-c(
    "sig","depth","dba","ice_conc",
    "ice_past","AW","strat","rho",
    "lat","lon")) 

# order and reformat
pls.WSC.tax <- pls.WSC.tax[
  mixedorder(pls.WSC.tax$RAS_id),] %>%
  column_to_rownames("RAS_id") %>%
  drop_na()

# subset to most abundant families
pls.WSC.tax <- pls.WSC.tax[
  ,colMeans(pls.WSC.tax) > 1] %>%
  as.matrix()

####################################

pls.EGC.tax <- sum.family %>%
  filter(site=="EGC") %>% 
  pivot_wider(
    data = .,
    id_cols = RAS_id, 
    names_from = Family, 
    values_from = c("Abundance")) 

pls.EGC.env <- ENV.bac %>%
  filter(site=="EGC") %>%
  dplyr::select_if(., is.numeric) %>%
  dplyr::select(-c(
    "sig","depth","dba","AW",
    "strat","rho","lat","lon")) 

# order and reformat
pls.EGC.tax <- pls.EGC.tax[
  mixedorder(pls.EGC.tax$RAS_id),] %>%
  column_to_rownames("RAS_id") 

#subset
pls.EGC.tax <- pls.EGC.tax[
  ,colMeans(pls.EGC.tax) > 1]

####################################

pls.WSC <- pls(
  pls.WSC.tax,
  pls.WSC.env,
  ncomp = 2,
  scale = T,
  mode = c("regression"),
  tol = 1e-06,
  max.iter = 100,
  near.zero.var = F,
  logratio = "none",
  multilevel = NULL,
  all.outputs = T)

pls.EGC <- pls(
  pls.EGC.tax,
  pls.EGC.env,
  ncomp = 2,
  scale = T,
  mode = c("regression"),
  tol = 1e-06,
  max.iter = 100,
  near.zero.var = F,
  logratio = "none",
  multilevel = NULL,
  all.outputs = T)

####################################

cim(pls.WSC,
    color = col.heat,
    row.names = T, col.names = T,
    row.cex = 0.8, col.cex = 0.8,
    cutoff = 0.5,
    cluster = "both",
    dist.method = c(
      "euclidean","euclidean"),
    clust.method = c(
      "complete","complete"),
    cut.tree = c(0,0),
    symkey = T,
    keysize = c(1, 1),
    keysize.label = 1,
    margins = c(18,26),
    center = T,
    transpose = F,
    mapping = "XY",
    save = "pdf",
    name.save = "Fig5_PLS-WSC")

cim(pls.EGC,
    color = col.heat,
    row.names = T, col.names = T,
    row.cex = 0.8, col.cex = 0.8,
    cutoff = 0.5,
    cluster = "both",
    dist.method = c(
      "euclidean","euclidean"),
    clust.method = c(
      "complete","complete"),
    cut.tree = c(0,0),
    symkey = T,
    keysize = c(1, 1),
    keysize.label = 1,
    margins = c(18,26),
    center = T,
    transpose = F,
    mapping = "XY",
    save = "pdf",
    name.save = "Fig5_PLS-EGC")

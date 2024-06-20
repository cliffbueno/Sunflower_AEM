# Sunflower Microbiome Analysis
# By Cliff Bueno de Mesquita, Kane Lab, Spring 2024
# 16S and ITS sequencing for 10 genotypes from 15 locations
# Take clean input files from Corinne Walsh and perform downstream analyses

# Document outline (click view document outline in top right):
# 1. Setup
# 2. Sclerotinia
# 3. Map
# 4. Alpha
# 5. Beta
# 6. Taxa
# 7. Venn
# 8. Site subsets
# 9. Genotype subsets
# 10. NCBI



#### 1. Setup ####
# Libraries
library(plyr) # Data manipulation
library(tidyverse) # Data manipulation
library(mctoolsr) # Microbial analyses
library(RColorBrewer) # Colors
library(vegan) # Multivariate analyses
library(indicspecies) # Indicator species
library(car) # Stats
library(FSA) # SE
library(magrittr) # Set names
library(PMCMRplus) # Stats
library(readxl) # Excel
library(writexl) # Excel
library(plotly) # Interactive plots
library(ggmap) # Maps
library(ggsn) # Maps
library(multcomp) # Tukey HSD and significance letters
library(emmeans) # Tukey HSD and significance letters
library(scales) # View colors
library(cowplot) # Multipanels
library(qvalue) # q values for indicator species
library(reshape2) # melt
library(gridExtra) # graphs
library(grid) # graphs
library(cowplot) # graphs
library(ggpubr) # graphs
library(ggExtra) # graphs
library(ggh4x) # graphs
library(dendextend) # graphs
library(corrplot) # correlation plots
library(pheatmap) # heatmaps
library(zCompositions) # CLR
library(compositions) # Aitchison
library(mobr) # rarefaction curves
library(plotly) # interactive graphs
library(patchwork) # insets
library(ggbiplot) # ordinations
library(ggfortify) # ordinations
library(ggrepel) # repel text
library(usmap) # for map
library(microseq) # for fastas
library(biomformat) # for exporting .biom
library(ggpmisc) # for graphs
library(geos); packageVersion("geos") # for geospatial
library(sp) # for raster
library(fields) # for geographic distance
library(gdm) # generalized dissimilarity modeling
library(ggdendro) # dendrograms
library(TSdist) # for gdm
library(gdm) # for gdm
library(pairwiseAdonis) # pairwise permanova
library(FUNGuildR) # funguild in R
library(phyloseq) # for UniFrac
library(ape) # For trees
library(phytools) # For trees
library(bestglm) # Find best combo of predictor variables
library(rsq) # Partial R2 values
library(MicEco) # Venn diagrams
library(SpiecEasi) # Networks
library(igraph) # View networks
library(rnetcarto) # Networks
library(gamlss) # Zero inflated regression
library(NBZIMM) # Neg. binomial, 0 inflated mixed models
library(ggspatial) # map annotation

# Functions
find_hull <- function(df) df[chull(df$Axis01, df$Axis02),]

`%notin%` <- Negate(`%in%`)

save_pheatmap_pdf <- function(x, filename, width = 4, height = 7) {
  stopifnot(!missing(x))
  stopifnot(!missing(filename))
  pdf(filename, width=width, height=height)
  grid::grid.newpage()
  grid::grid.draw(x$gtable)
  dev.off()
}

save_pheatmap_png <- function(x, filename, width = 4, height = 7) {
  stopifnot(!missing(x))
  stopifnot(!missing(filename))
  png(filename, width=width, height=height, units = "in", res = 300)
  grid::grid.newpage()
  grid::grid.draw(x$gtable)
  dev.off()
}

dist_geo <- function(lat_a, lon_a, lat_b, lon_b) { 
  if(anyNA(c(lat_a, lon_a, lat_b, lon_b))) return(NA) 
  round(distm(c(lon_a, lat_a), c(lon_b, lat_b), fun = distHaversine)/1000,2) 
}

MyDiamond <- function(coords, v=NULL, params) {
  vertex.color <- params("vertex", "color")
  if (length(vertex.color) != 1 && !is.null(v)) {
    vertex.color <- vertex.color[v]
  }
  vertex.frame.color <- params("vertex", "frame.color")
  if (length(vertex.frame.color) != 1 && !is.null(v)) {
    vertex.frame.color <- vertex.frame.color[v]
  }
  vertex.size <- 1/200 * params("vertex", "size")
  if (length(vertex.size) != 1 && !is.null(v)) {
    vertex.size <- vertex.size[v]
  }
  
  symbols(x=coords[,1], y=coords[,2], bg=vertex.color, fg=vertex.frame.color,
          stars=cbind(vertex.size, vertex.size, vertex.size, vertex.size),
          add=TRUE, inches=FALSE)
}
add_shape("diamond", clip=shapes("circle")$clip,
          plot=MyDiamond, parameters=list(vertex.frame.color="white",
                                          vertex.frame.width=1))

mytriangle <- function(coords, v=NULL, params) {
  vertex.color <- params("vertex", "color")
  if (length(vertex.color) != 1 && !is.null(v)) {
    vertex.color <- vertex.color[v]
  }
  vertex.size <- 1/200 * params("vertex", "size")
  if (length(vertex.size) != 1 && !is.null(v)) {
    vertex.size <- vertex.size[v]
  }
  
  symbols(x=coords[,1], y=coords[,2], bg=vertex.color,
          stars=cbind(vertex.size, vertex.size, vertex.size),
          add=TRUE, inches=FALSE)
}
add_shape("triangle", clip=shapes("circle")$clip,
          plot=mytriangle)

# Repository path
setwd("~/Documents/GitHub/SunflowerGxE/")

# Plotting functions
source("code/cliffplot_taxa_bars.R")
source("code/plot_multipatt_asv.R")

# Effect size (from Jack Darcy, Science Advances 2018)
source("code/effectSize.R")



#### _Data ####
# Prepare data. Do once. After, can just load in the _Start Here section

# Corinne got CMI and GDD, but we should get basic temp and precip as well
# Use CHELSA 2.1 https://chelsa-climate.org/downloads/
# Download BIO1 (mean annual T) and BIO12 (mean annual P) for 1981-2010
# Also get BIO4 (T seasonality) and BIO15 (P seasonality) for 1981-2010
# These 4 variables were used by Cliff in the chimpanzee gut microbiome paper
# Couldn't get raster working on PC so do on Innes server. Made input_filt_wClim.RDS.



#### __16S ####
# Check different input files from Corinne
input_16S <- readRDS("data/input_filt.RDS")
nrow(input_16S$data_loaded) # 27928 ASVs
nrow(input_16S$map_loaded) # 587 samples
ncol(input_16S$map_loaded) # 27 columns

input_16S <- readRDS("data/input_filt_rar9k.RDS")
nrow(input_16S$data_loaded) # 27913 ASVs
nrow(input_16S$map_loaded) # 575 samples
ncol(input_16S$map_loaded) # 27 columns

input_16S <- readRDS("data/input_filt_relab.RDS")
nrow(input_16S$data_loaded) # 27928 ASVs
nrow(input_16S$map_loaded) # 587 samples
ncol(input_16S$map_loaded) # 27 columns

input_16S <- readRDS("data/input.filt2C.RDS")
nrow(input_16S$data_loaded) # 10960 ASVs, removed ASVs present in fewer than 30 samples
nrow(input_16S$map_loaded) # 587 samples
ncol(input_16S$map_loaded) # 31 columns

input_16S_wCMI_GDD <- readRDS("data/input.filtC.rds")
nrow(input_16S_wCMI_GDD$data_loaded) # 27928 ASVs
nrow(input_16S_wCMI_GDD$map_loaded) # 587 samples
ncol(input_16S_wCMI_GDD$map_loaded) # 31 columns

# Use file with climate data added by Cliff and all ASVs
# Need to also add CMI and GDD from Corinne
input_16S <- readRDS("data/input_filt_wClim.RDS")
nrow(input_16S$data_loaded) # 27928 ASVs
nrow(input_16S$map_loaded) # 587 samples
ncol(input_16S$map_loaded) # 31 columns
sum(input_16S$map_loaded$sample != input_16S_wCMI_GDD$map_loaded$sample)
input_16S$map_loaded$CMI <- input_16S_wCMI_GDD$map_loaded$CMI
input_16S$map_loaded$GDD <- input_16S_wCMI_GDD$map_loaded$GDD
input_16S$map_loaded$CMI_scaled <- input_16S_wCMI_GDD$map_loaded$cmi.scl
input_16S$map_loaded$GDD_scaled <- input_16S_wCMI_GDD$map_loaded$gdd.scl

# Check that Chloroplast, Mitochondria, Eukaryota, Domain NA are all filtered
input_filt_16S <- input_16S
input_filt_16S <- filter_taxa_from_input(input_filt_16S,
                                     taxa_to_remove = "Chloroplast") # None found
input_filt_16S <- filter_taxa_from_input(input_filt_16S,
                                     taxa_to_remove = "Mitochondria") # None found
input_filt_16S <- filter_taxa_from_input(input_filt_16S,
                                     taxa_to_remove = "Eukaryota") # None found
input_filt_16S <- filter_taxa_from_input(input_filt_16S,
                                     taxa_to_remove = "NA",
                                     at_spec_level = 1) # None found

# Good. Now check for singletons and doubletons
singdoub_16S <- data.frame("count" = rowSums(input_filt_16S$data_loaded)) %>%
  filter(count < 3) %>%
  mutate(ASV = rownames(.)) # 0 observations

sort(colSums(input_filt_16S$data_loaded))
mean(colSums(input_filt_16S$data_loaded)) # 23783.96
se(colSums(input_filt_16S$data_loaded)) # 532.2505

# Check rarefaction curve
rarecurve(t(input_filt_16S$data_loaded), step = 1000)
rarecurve(t(input_filt_16S$data_loaded), step = 1000, label = FALSE)
rc_16S <- rarecurve(t(input_filt_16S$data_loaded), step = 1000, label = FALSE, tidy = TRUE)
ggplot(rc_16S, aes(Sample, Species, group = Site)) +
  geom_line() +
  geom_vline(xintercept = 6027, colour = "red") +
  theme_bw()

# Add unrarefied ASV richness and Shannon diversity
input_filt_16S$map_loaded$rich <- specnumber(input_filt_16S$data_loaded, 
                                             MARGIN = 2)
input_filt_16S$map_loaded$shannon <- vegan::diversity(input_filt_16S$data_loaded, 
                                                      index = "shannon", 
                                                      MARGIN = 2)

# Change Sclerotinia variable name
input_filt_16S$map_loaded$Sclerotinia <- input_filt_16S$map_loaded$Average.of.percent.incidence

# Rarefy. Use for taxonomic analysis. Use unrarefied data from compositional analysis.
# Rarefy at 6027 (lowest). Use this.
set.seed(530)
input_filt_rare_16S <- single_rarefy(input_filt_16S, 6027) # n = 587 still
sort(colSums(input_filt_rare_16S$data_loaded))
table(input_filt_rare_16S$map_loaded$Site)
# Only 3 sites without 40 samples

# Rarefy at 10269? Better but drops 23 samples.
# set.seed(530)
# input_filt_rare_16S <- single_rarefy(input_filt_16S, 10269) # n = 564
# sort(colSums(input_filt_rare_16S$data_loaded))
# table(input_filt_rare_16S$map_loaded$Site)
# Now 8 sites without 40 samples.

# Add rarefied richness and shannon
input_filt_rare_16S$map_loaded$rich <- specnumber(input_filt_rare_16S$data_loaded, 
                                                  MARGIN = 2)
input_filt_rare_16S$map_loaded$shannon <- vegan::diversity(input_filt_rare_16S$data_loaded, 
                                                           index = "shannon", 
                                                           MARGIN = 2)



#### __ITS ####
input_ITS <- readRDS("data/input_filt_ITS_CW.RDS")
nrow(input_ITS$data_loaded) # 5625 ASVs
nrow(input_ITS$map_loaded) # 587 samples
ncol(input_ITS$map_loaded) # 27 columns

# Add climate data
sum(input_ITS$map_loaded$sample != input_16S_wCMI_GDD$map_loaded$sample) # Match
sum(input_ITS$map_loaded$sample != input_16S_wCMI_GDD$map_loaded$sample) # Match
input_ITS$map_loaded$CMI <- input_16S_wCMI_GDD$map_loaded$CMI
input_ITS$map_loaded$GDD <- input_16S_wCMI_GDD$map_loaded$GDD
input_ITS$map_loaded$CMI_scaled <- input_16S_wCMI_GDD$map_loaded$cmi.scl
input_ITS$map_loaded$GDD_scaled <- input_16S_wCMI_GDD$map_loaded$gdd.scl
input_ITS$map_loaded$Temp <- input_16S$map_loaded$Temp
input_ITS$map_loaded$P <- input_16S$map_loaded$P
input_ITS$map_loaded$Tseas <- input_16S$map_loaded$Tseas
input_ITS$map_loaded$Pseas <- input_16S$map_loaded$Pseas

# Check that Bacteria, Archaea, Domain NA are all filtered.
input_filt_ITS <- input_ITS
input_filt_ITS <- filter_taxa_from_input(input_filt_ITS,
                                         taxa_to_remove = "Bacteria") # None found
input_filt_ITS <- filter_taxa_from_input(input_filt_ITS,
                                         taxa_to_remove = "Archaea") # None found
input_filt_ITS <- filter_taxa_from_input(input_filt_ITS,
                                         taxa_to_remove = "NA",
                                         at_spec_level = 1) # None found

# Important for ITS: filter out all non-fungi! keep k__Fungi
# (ITS includes DNA from plants and animals)
# Check other "kingdoms" here:
tax_sum_kingdom <- summarize_taxonomy(input = input_filt_ITS, level = 1, report_higher_tax = F)
# k__fungi is >99% of reads but still, for accuracy, we need to filter.
input_filt_ITS <- filter_taxa_from_input(input_filt_ITS,
                                         taxa_to_keep = "k__Fungi",
                                         at_spec_level = 1) # 304 taxa removed

# Good. Now check for singletons and doubletons
singdoub_ITS <- data.frame("count" = rowSums(input_filt_ITS$data_loaded)) %>%
  filter(count < 3) %>%
  mutate(ASV = rownames(.)) # 0 observations

sort(colSums(input_filt_ITS$data_loaded))
mean(colSums(input_filt_ITS$data_loaded)) # 23836.49
se(colSums(input_filt_ITS$data_loaded)) # 394.1155

# Check rarefaction curve
rarecurve(t(input_filt_ITS$data_loaded), step = 1000)
rarecurve(t(input_filt_ITS$data_loaded), step = 1000, label = FALSE)
rc_ITS <- rarecurve(t(input_filt_ITS$data_loaded), step = 1000, label = FALSE, tidy = TRUE)
ggplot(rc_ITS, aes(Sample, Species, group = Site)) +
  geom_line() +
  geom_vline(xintercept = 2581, colour = "red") +
  geom_vline(xintercept = 4809, colour = "red") +
  geom_vline(xintercept = 10160, colour = "red") +
  theme_bw()
# 2581 is good for some samples but not for the most rich samples
# 4809 is better and would only drop 1 sample
# 10160 good for all samples

# Add unrarefied ASV richness and Shannon diversity
input_filt_ITS$map_loaded$rich <- specnumber(input_filt_ITS$data_loaded, 
                                             MARGIN = 2)
input_filt_ITS$map_loaded$shannon <- vegan::diversity(input_filt_ITS$data_loaded, 
                                                      index = "shannon", 
                                                      MARGIN = 2)

# Change Sclerotinia variable name
input_filt_ITS$map_loaded$Sclerotinia <- input_filt_ITS$map_loaded$Average.of.percent.incidence

# Rarefy. Use for taxonomic analysis. Use unrarefied data from compositional analysis.
# Rarefy at 2581 (lowest)?
# set.seed(530)
# input_filt_rare_ITS <- single_rarefy(input_filt_ITS, 2581) # n = 587 still
# sort(colSums(input_filt_rare_ITS$data_loaded))
# table(input_filt_rare_ITS$map_loaded$Site)

# Also rarefy at 4809. Drops one sample from Kirkmeyer. Use this.
set.seed(530)
input_filt_rare_ITS <- single_rarefy(input_filt_ITS, 4809) # n = 586
sort(colSums(input_filt_rare_ITS$data_loaded))
table(input_filt_rare_ITS$map_loaded$Site)

# Rarefy at 10160?
# set.seed(530)
# input_filt_rare_ITS <- single_rarefy(input_filt_ITS, 10160) # n = 555
# sort(colSums(input_filt_rare_ITS$data_loaded))
# table(input_filt_rare_ITS$map_loaded$Site)

# Add rarefied richness and shannon
input_filt_rare_ITS$map_loaded$rich <- specnumber(input_filt_rare_ITS$data_loaded, 
                                                  MARGIN = 2)
input_filt_rare_ITS$map_loaded$shannon <- vegan::diversity(input_filt_rare_ITS$data_loaded, 
                                                           index = "shannon", 
                                                           MARGIN = 2)



#### __Function ####
# FAPROTAX annotation (16S)
faprotax_input <- input_filt_rare_16S$data_loaded %>%
  mutate(taxonomy = paste(input_filt_rare_16S$taxonomy_loaded$taxonomy1,
                          input_filt_rare_16S$taxonomy_loaded$taxonomy2,
                          input_filt_rare_16S$taxonomy_loaded$taxonomy3,
                          input_filt_rare_16S$taxonomy_loaded$taxonomy4,
                          input_filt_rare_16S$taxonomy_loaded$taxonomy5,
                          input_filt_rare_16S$taxonomy_loaded$taxonomy6,
                          sep = "; ")) %>%
  dplyr::select(taxonomy, everything())
#write.table(faprotax_input, "data/faprotax_input.txt", sep = "\t", row.names = F)
# Note: Had to cut and paste "taxonomy" header in order for it to be found properly.

# Run in terminal
# ./Desktop/FAPROTAX_1.2.7/collapse_table.py -i ~/Documents/GitHub/SunflowerGxE/data/faprotax_input.txt -o ~/Desktop/faprotax_output.tsv -r ~/Desktop/faprotax_report.txt -g ./Desktop/FAPROTAX_1.2.7/FAPROTAX.txt -d "taxonomy" --group_leftovers_as "Unassigned" -v
# 5373 out of 27893 records (19.2629 %) were assigned to at least one group



# PiCRUST2 annotation (16S)
# Need repset.fasta file and ASV table as .biom
asv_biom <- make_biom(input_filt_rare_16S$data_loaded)
#write_biom(asv_biom, "data/ASV_table.biom")

# Need to make repset with just ASV_ID
gxe_repset <- readFasta("data/repset_16S_filt.fasta") %>%
  separate(Header, into = c("Header", "taxonomy"), sep = " ", remove = T) %>%
  dplyr::select(Header, Sequence, -taxonomy)
writeFasta(gxe_repset, "data/repset_16S_filt_ASVid.fasta")

# Run the default script in terminal on Innes Server
# cd ~/Documents/GitHub/SunflowerGxE
# picrust2_pipeline.py -s data/repset_16S_filt_ASVid.fasta -i data/ASV_table.biom -o ./picrust2_out_pipeline -p 10



# FUNGuild annotation (ITS)
funguild_input <- input_filt_rare_ITS$data_loaded %>%
  rownames_to_column(var = "ASV_ID") %>%
  dplyr::select(ASV_ID, everything()) %>%
  mutate(taxonomy = paste(input_filt_rare_ITS$taxonomy_loaded$taxonomy1,
                          input_filt_rare_ITS$taxonomy_loaded$taxonomy2,
                          input_filt_rare_ITS$taxonomy_loaded$taxonomy3,
                          input_filt_rare_ITS$taxonomy_loaded$taxonomy4,
                          input_filt_rare_ITS$taxonomy_loaded$taxonomy5,
                          input_filt_rare_ITS$taxonomy_loaded$taxonomy6,
                          input_filt_rare_ITS$taxonomy_loaded$taxonomy8,
                          sep = ";"))
fung_guilds <- funguild_assign(funguild_input,
                               tax_col = "taxonomy")
#saveRDS(fung_guilds, "data/fung_guilds.rds")
fung_guilds <- readRDS("data/fung_guilds.rds")


#### __Soil ####
# Got soil data later on from Brent. Add to input data here. And calculate "soil distance"
# There are 2 replicates from each site. Join by Site_Rep
# 15 variables measured at 0-6 cm depth
# Cadmium had a couple below detection - make numeric and it will set those to NA
# For now replace those with 0 since they were so low
soil <- read_excel("data/soil tests with Cd included.xlsx", sheet = 2) %>%
  mutate(Cd_ppm = as.numeric(Cd_ppm)) %>%
  replace_na(list(Cd_ppm = 0))
names(soil)
clay_sort <- soil %>%
  group_by(Site) %>%
  summarise(Clay = mean(Clay_perc, na.rm = T)) %>%
  arrange(Clay)
soil_long <- soil %>%
  dplyr::select(Site, 7:21) %>%
  pivot_longer(cols = c(2:16),
               names_to = "variable") %>%
  mutate(variable = factor(variable,
                           levels = c("Clay_perc", "Silt_perc", "Sand_perc",
                                      "OM_perc", "NO3_ppm", "P_ppm",
                                      "pH", "Salts_mmhos_cm", "CEC_meq_g", 
                                      "K_ppm", "Ca_ppm", "Mg_ppm",
                                      "Na_ppm", "S_ppm", "Cd_ppm"))) %>%
  mutate(Site = factor(Site,
                       levels = clay_sort$Site))
facet_names <- c("Clay_perc" = "Clay~('%')",
                 "Silt_perc" = "Silt~('%')",
                 "Sand_perc" = "Sand~('%')",
                 "OM_perc" = "OM~('%')",
                 "NO3_ppm" = "NO[3]^'-'~(ppm)",
                 "P_ppm" = "P~(ppm)",
                 "pH" = "pH",
                 "Salts_mmhos_cm" = "Salts~(mmhos/cm)",
                 "CEC_meq_g" = "CEC~(meq/g)",
                 "K_ppm" = "K~(ppm)",
                 "Ca_ppm" = "Ca~(ppm)",
                 "Mg_ppm" = "Mg~(ppm)",
                 "Na_ppm" = "Na~(ppm)",
                 "S_ppm" = "S~(ppm)",
                 "Cd_ppm" = "Cd~(ppm)")
pdf("InitialFigs/Soil.pdf", width = 7, height = 7)
ggplot(soil_long, aes(x = Site, y = value)) +
  geom_point(size = 2, alpha = 0.75, pch = 16) +
  labs(x = "Site",
       y = NULL) +
  scale_y_continuous(expand = c(0.1, 0.1)) +
  facet_wrap(~ variable, ncol = 3, scales = "free_y",
             labeller = as_labeller(facet_names, default = label_parsed)) +
  theme_bw() +
  theme(axis.text.x = element_text(size = 8, angle = 45, hjust = 1))
dev.off()

texture <- soil %>%
  group_by(Site) %>%
  slice_head(n = 1) %>%
  dplyr::select(Site, Texture)

# By Site (originally did this)
soil_toMerge <- soil %>%
  dplyr::select(-Ref_No, -Lab_No) %>%
  group_by(Site) %>%
  summarise_if(is.numeric, mean) %>%
  left_join(., texture, by = "Site")

# Check
sum(soil_toMerge$Site %notin% input_filt_16S$map_loaded$Site)

# By Site/Replicate (use this)
soil_toMerge <- soil %>%
  separate(SampleID, into = c("Code", "Rep"), sep = "-") %>%
  mutate(Site_Rep = paste(Site, Rep, sep = "_")) %>%
  dplyr::select(Site_Rep, everything(), -Ref_No, -Lab_No, -Code, -Rep, -Site)



#### __Phenotype ####
# Summarize chlorophyl for Table
chloro <- input_filt_rare_16S$map_loaded %>%
  group_by(Entry) %>%
  summarise(Chlorophyll = mean(chlorophyll, na.rm = TRUE))

#### ___GRIN ####
# Get other phenotype data from USDA GRIN Website
# Download "Observations" spreadsheet for each
# Note HA 466 = PI 667183
# Read in, manipulate, and merge the spreadsheets
pt <- list()
f <- list.files("data/phenotype/")
setwd("data/phenotype/")
for (i in 1:length(f)) {
  pt[[i]] <- read_excel(f[i], skip = 1)
  pt[[i]]$Genotype <- levels(as.factor(input_filt_rare_16S$map_loaded$Entry))[i]
}
View(pt[[1]])
pt_merged <- pt[[1]]
for (i in 2:length(pt)) {
  pt_merged <- rbind(pt_merged, pt[[i]])
}

# Messy dataset. Some items with multiple values. Some just on 1-9 scale.
# Try getting just continuous variables. Then we can group and average by genotype.
# Also need to get variables that were recorded for all 10 genotypes
pt_merged_cont <- pt_merged %>%
  mutate(Value = as.numeric(Value)) %>%
  drop_na(Value) %>%
  group_by(Description, Genotype) %>%
  summarise(Value = mean(Value),
            SampleSize = n())

# Count number of rows and get variables with 9 or 10 genotypes
freq <- as.data.frame(table(pt_merged_cont$Description)) %>%
  filter(Freq >= 9)

# Filter
pt_merged_cont_filt <- pt_merged_cont %>%
  filter(Description %in% freq$Var1)

# Widen
pt_wide <- pt_merged_cont_filt %>%
  dplyr::select(-SampleSize) %>%
  pivot_wider(names_from = Description, values_from = Value) %>%
  set_names(c("Entry",
              "DaysToFlower", "StemLengthHighRange", "StemLengthLowRange", "PerSeed_14_64", 
              "PerSeed_20_64", "PerNoBranch", "GRIN_Sclerotinia", "SeedWeight")) %>%
  mutate(GRIN_Sclerotinia = 100 - GRIN_Sclerotinia) # Convert % resistance to % incidence!

sclero_GRIN <- pt_wide %>%
  dplyr::select(Entry, GRIN_Sclerotinia)
# Note those data are from "SUNFLOWER.ND.SCLEROTINIA.08-09". 

# Write for Table S2
write_xlsx(pt_wide, "~/Desktop/Sunflower/phenotype.xlsx", format_headers = F)

# Reset working directory
setwd("../../")



#### ____Sclerotinia ####
# Plot % Sclerotinia by genotype (note: only 1 value repeated for all samples of genotype)
# Have our own data from Carrington 2017. Also have data on GRIN (SD 08/09). Compare.
ggplot(input_filt_16S$map_loaded, aes(reorder(Entry, Sclerotinia, mean), Sclerotinia)) +
  geom_point(size = 5) +
  labs(x = "Genotype",
       y = "% Sclerotinia incidence") +
  theme_bw()

ggplot(pt_wide, aes(reorder(Entry, GRIN_Sclerotinia, mean), GRIN_Sclerotinia)) +
  geom_point(size = 5) +
  labs(x = "Genotype",
       y = "% Sclerotinia incidence") +
  theme_bw()

# Points are duplicated, remake
sclero_df <- input_filt_16S$map_loaded %>%
  group_by(Entry) %>%
  summarise(Sclerotinia = mean(Sclerotinia))

ggplot(sclero_df, aes(reorder(Entry, Sclerotinia, mean), Sclerotinia)) +
  geom_point(size = 5) +
  labs(x = "Genotype",
       y = "% Sclerotinia incidence") +
  theme_bw()

ggplot(sclero_df, aes(Entry, Sclerotinia)) +
  geom_point(size = 5) +
  labs(x = "Genotype",
       y = "% Sclerotinia incidence") +
  coord_flip() +
  theme_bw()

# Plot with genetic distance
gen.dist <- as.dist(readRDS("data/gen.dist.mat.rds"))
hc <- hclust(gen.dist, method = "ward.D2")
plot(hc)
ggdendrogram(hc, rotate = TRUE, size = 2)
dhc <- as.dendrogram(hc)
ddata <- dendro_data(dhc, type = "rectangle")
ggplot(segment(ddata)) + 
  geom_segment(aes(x = x, y = y, xend = xend, yend = yend)) + 
  geom_text(data = ddata$labels, 
            aes(x = x, y = y-1, label = label), size = 3, vjust = 0,
            inherit.aes = F) +
  coord_flip() + 
  scale_y_reverse(expand = c(0.1, 0)) +
  theme_dendro()

g_dend <- ggplot(segment(ddata)) + 
  geom_segment(aes(x = x, y = y, xend = xend, yend = yend)) + 
  coord_flip() + 
  scale_y_reverse(expand = c(0,0)) +
  scale_x_continuous(expand = c(0.06,0.06)) +
  theme_dendro()
g_dend

sclero_df <- sclero_df %>%
  mutate(Entry = gsub(" ", "", Entry)) %>%
  mutate(Entry = factor(Entry,
                        levels = c(ddata$labels$label)))

g_sclero <- ggplot(sclero_df, aes(Entry, Sclerotinia)) +
  geom_point(size = 5) +
  labs(x = "Genotype",
       y = "% Sclerotinia incidence") +
  coord_flip() +
  theme_bw() +
  theme(axis.title.y = element_blank())
g_sclero

pdf("InitialFigs/SunflowerGenotypes_Sclero.pdf", width = 7, height = 5)
plot_grid(g_dend, g_sclero, align = "h", rel_widths = c(0.3, 0.7))
dev.off()

# Add GRIN
sclero_df <- input_filt_16S$map_loaded %>%
  group_by(Entry) %>%
  summarise(Sclerotinia = mean(Sclerotinia)) %>%
  left_join(., sclero_GRIN, by = "Entry") %>%
  pivot_longer(cols = c(Sclerotinia, GRIN_Sclerotinia)) %>%
  mutate(Entry = gsub(" ", "", Entry)) %>%
  mutate(Entry = factor(Entry,
                        levels = c(ddata$labels$label)))

g_sclero <- ggplot(sclero_df, aes(Entry, value, colour = name)) +
  geom_point(size = 5) +
  labs(x = "Genotype",
       y = "% Sclerotinia incidence",
       colour = "Dataset") +
  scale_colour_manual(values = c("blue", "red")) +
  coord_flip() +
  theme_bw() +
  theme(axis.title.y = element_blank())
g_sclero

pdf("InitialFigs/SunflowerGenotypes_Sclero_Both.pdf", width = 7, height = 5)
plot_grid(g_dend, g_sclero, align = "h", rel_widths = c(0.3, 0.7))
dev.off()

# Compare
sclero_df <- input_filt_16S$map_loaded %>%
  group_by(Entry) %>%
  summarise(Sclerotinia = mean(Sclerotinia)) %>%
  left_join(., sclero_GRIN, by = "Entry")
summary(lm(Sclerotinia ~ GRIN_Sclerotinia, data = sclero_df))
pdf("InitialFigs/CompareSclerotinia.pdf", width = 7, height = 5)
ggplot(sclero_df, aes(Sclerotinia, GRIN_Sclerotinia)) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed") +
  geom_smooth(method = "lm", linetype = "dashed") +
  geom_point(size = 5) +
  labs(x = "% Sclerotinia incidence (Carrington)",
       y = "% Sclerotinia incidence (GRIN)",
       colour = "Dataset") +
  xlim(0,100) +
  ylim(0,100) +
  theme_bw()
dev.off() 



#### _Merge ####
# Merge to the mctoolsr input data
# Make sure rownames are retained in map_loaded, certain tidyverse manipulations remove them
# Merge plant data by Entry. Merge soil data by Site_Rep.
rownames(input_filt_16S$map_loaded)
rownames(input_filt_rare_16S$map_loaded)
rownames(input_filt_ITS$map_loaded)
rownames(input_filt_rare_ITS$map_loaded)

input_filt_16S$map_loaded$Site_Rep <- paste(input_filt_16S$map_loaded$Site, 
                                            input_filt_16S$map_loaded$Rep,
                                            sep = "_")
input_filt_rare_16S$map_loaded$Site_Rep <- paste(input_filt_rare_16S$map_loaded$Site, 
                                                 input_filt_rare_16S$map_loaded$Rep,
                                                 sep = "_")
input_filt_ITS$map_loaded$Site_Rep <- paste(input_filt_ITS$map_loaded$Site, 
                                            input_filt_ITS$map_loaded$Rep,
                                            sep = "_")
input_filt_rare_ITS$map_loaded$Site_Rep <- paste(input_filt_rare_ITS$map_loaded$Site, 
                                                 input_filt_rare_ITS$map_loaded$Rep,
                                                 sep = "_")

input_filt_16S$map_loaded <- input_filt_16S$map_loaded %>%
  mutate(sampleID = rownames(.)) %>%
  left_join(., soil_toMerge, by = "Site_Rep") %>%
  left_join(., pt_wide, by = "Entry") %>%
  column_to_rownames(var = "sampleID") %>%
  mutate(sampleID = rownames(.))

input_filt_rare_16S$map_loaded <- input_filt_rare_16S$map_loaded %>%
  mutate(sampleID = rownames(.)) %>%
  left_join(., soil_toMerge, by = "Site_Rep") %>%
  left_join(., pt_wide, by = "Entry") %>%
  column_to_rownames(var = "sampleID") %>%
  mutate(sampleID = rownames(.))

input_filt_ITS$map_loaded <- input_filt_ITS$map_loaded %>%
  mutate(sampleID = rownames(.)) %>%
  left_join(., soil_toMerge, by = "Site_Rep") %>%
  left_join(., pt_wide, by = "Entry") %>%
  column_to_rownames(var = "sampleID") %>%
  mutate(sampleID = rownames(.))

input_filt_rare_ITS$map_loaded <- input_filt_rare_ITS$map_loaded %>%
  mutate(sampleID = rownames(.)) %>%
  left_join(., soil_toMerge, by = "Site_Rep") %>%
  left_join(., pt_wide, by = "Entry") %>%
  column_to_rownames(var = "sampleID") %>%
  mutate(sampleID = rownames(.))

rownames(input_filt_16S$map_loaded)
rownames(input_filt_rare_16S$map_loaded)
rownames(input_filt_ITS$map_loaded)
rownames(input_filt_rare_ITS$map_loaded)

# Save
# saveRDS(input_filt_16S, "data/input_filt_16S.rds")
# saveRDS(input_filt_rare_16S, "data/input_filt_rare_16S.rds")
# saveRDS(input_filt_ITS, "data/input_filt_ITS.rds")
# saveRDS(input_filt_rare_ITS, "data/input_filt_rare_ITS.rds")



#### _Start Here ####
# Import mctoolsr inputs
input_filt_16S <- readRDS("data/input_filt_16S.rds")
input_filt_rare_16S <- readRDS("data/input_filt_rare_16S.rds")
nrow(input_filt_16S$data_loaded) # 27928
nrow(input_filt_rare_16S$data_loaded) # 27893

input_filt_ITS <- readRDS("data/input_filt_ITS.rds")
input_filt_rare_ITS <- readRDS("data/input_filt_rare_ITS.rds")
nrow(input_filt_ITS$data_loaded) # 5321
nrow(input_filt_rare_ITS$data_loaded) # 5294

# Change Kirkmeyer to Brighton
input_filt_16S$map_loaded$Site <- gsub("Kirkmeyer", "Brighton", input_filt_16S$map_loaded$Site)
input_filt_rare_16S$map_loaded$Site <- gsub("Kirkmeyer", "Brighton", input_filt_rare_16S$map_loaded$Site)
input_filt_ITS$map_loaded$Site <- gsub("Kirkmeyer", "Brighton", input_filt_ITS$map_loaded$Site)
input_filt_rare_ITS$map_loaded$Site <- gsub("Kirkmeyer", "Brighton", input_filt_rare_ITS$map_loaded$Site)

# Make factors
input_filt_rare_16S$map_loaded$Site <- as.factor(input_filt_rare_16S$map_loaded$Site)
input_filt_rare_ITS$map_loaded$Site <- as.factor(input_filt_rare_ITS$map_loaded$Site)
input_filt_rare_16S$map_loaded$Entry <- as.factor(input_filt_rare_16S$map_loaded$Entry)
input_filt_rare_ITS$map_loaded$Entry <- as.factor(input_filt_rare_ITS$map_loaded$Entry)

# Import distance matrices (these were originally made in the "_Distance" section)
# Unifrace and genetic distance took a while to create so this is more efficient
unifrac.distance <- readRDS("data/unifrac.distance.rds")
geography.distance.mat <- readRDS("data/geography.distance.mat.rds")
geography.distance.mat_ITS <- readRDS("data/geography.distance.mat_ITS.rds")
climate.distance <- readRDS("data/climate.distance.rds")
climate.distance_ITS <- readRDS("data/climate.distance_ITS.rds")
genetic.distance <- readRDS("data/genetic.distance.rds")
genetic.distance_ITS <- readRDS("data/genetic.distance_ITS.rds")
soil.distance <- readRDS("data/soil.distance.rds")
soil.distance_ITS <- readRDS("data/soil.distance_ITS.rds")
plant.distance <- readRDS("data/plant.distance.rds")
plant.distance_ITS <- readRDS("data/plant.distance_ITS.rds")

# Import OTU-ASV ID mapping file
otu_asv_map_16S <- readRDS("data/OTU_ID_16S.rds") %>%
  dplyr::select(ASV_ID, OTU_ID)
otu_asv_map_ITS <- readRDS("data/OTU_ID_ITS.rds") %>%
  dplyr::select(ASV_ID, OTU_ID)



#### 2. Map ####
# Make a site map. Good for supplementary figure.
# Note: Per email from Brent 5/28/24, Brighton coords are: 39.926920,-104.660118
coords <- input_filt_16S$map_loaded %>%
  dplyr::select(Site, LONG, LAT) %>%
  group_by(Site) %>%
  summarise(Latitude = mean(LAT),
            Longitude = mean(LONG))
coords$Latitude[1] <- 39.926920
coords$Longitude[2] <- -104.660118
test_data <- data.frame(lon = coords$Longitude, lat = coords$Latitude)
transformed_data <- usmap_transform(test_data)
coords <- coords %>%
  left_join(., transformed_data, by = c("Latitude" = "lat"))
pdf("InitialFigs/Map.pdf", width = 8, height = 6)
map <- plot_usmap(exclude = c("AK", "HI"),
           color = "white",
           fill = "grey80",
           size = 0.3) +
  geom_point(data = coords, 
             aes(x = x, y = y),
             fill = "red",
             color = "black",
             size = 4,
             shape = 21) +
  geom_text_repel(data = coords,
                  aes(x = x, y = y, label = Site),
                  size = 4) +
  annotation_scale(location = "bl") +
  annotation_north_arrow(location = "br") +
  theme(legend.position = "none")
map
dev.off()
#write.csv(coords, "~/Desktop/Sunflower/GIS/coords.csv")

# Compare BIOCLIM and CHELSA
chelsa <- input_filt_rare_16S$map_loaded %>%
  group_by(Site) %>%
  summarise(MAT = mean(Temp),
            MAP = mean(P),
            Tseas = mean(Tseas),
            Pseas = mean(Pseas))
bioclim <- read.csv("~/Desktop/Sunflower/info.csv")
sum(chelsa$Site != bioclim$Site)
plot(chelsa$MAT, bioclim$MeanT) # Good
plot(chelsa$MAP, bioclim$MeanP) # Good
plot(chelsa$Tseas, bioclim$Tseas) # Good
plot(chelsa$Pseas, bioclim$Pseas) # Bad



#### 3. Alpha ####
#### _16S ####
# Get descriptive info
min(input_filt_rare_16S$map_loaded$rich) # 569
max(input_filt_rare_16S$map_loaded$rich) # 2098
mean(input_filt_rare_16S$map_loaded$rich) # 1500
se(input_filt_rare_16S$map_loaded$rich) # 12

# Test and plot
leveneTest(input_filt_rare_16S$map_loaded$rich ~ input_filt_rare_16S$map_loaded$Site) # Not homogeneous
leveneTest(input_filt_rare_16S$map_loaded$rich ~ input_filt_rare_16S$map_loaded$Entry) # Homogeneous
m <- aov(rich ~ Site * Entry, data = input_filt_rare_16S$map_loaded)
summary(m)
Anova(m, type = "III", singular.ok = TRUE) # Site, Interaction
hist(m$residuals)
shapiro.test(m$residuals)
plot(m$fitted.values, m$residuals)
m <- aov(rich ~ Entry, data = input_filt_rare_16S$map_loaded)
summary(m)
m <- aov(rich ~ Site, data = input_filt_rare_16S$map_loaded)
summary(m)
TukeyHSD(m)
t <- emmeans(object = m, specs = "Site") %>%
  cld(object = ., adjust = "Tukey", Letters = letters, alpha = 0.05) %>%
  mutate(name = "rich",
         y = max(input_filt_rare_16S$map_loaded$rich)+
           (max(input_filt_rare_16S$map_loaded$rich)-
              min(input_filt_rare_16S$map_loaded$rich))/10,
         #y = 5650,
         Dataset = "16S")

leveneTest(input_filt_rare_16S$map_loaded$shannon ~ input_filt_rare_16S$map_loaded$Site) # Not homogenous
leveneTest(input_filt_rare_16S$map_loaded$shannon ~ input_filt_rare_16S$map_loaded$Entry) # Homogeneous
m1 <- aov(shannon ~ Site * Entry, data = input_filt_rare_16S$map_loaded)
summary(m1)
Anova(m1, type = "III", singular.ok = TRUE) # Site, Interaction
hist(m1$residuals)
shapiro.test(m1$residuals)
plot(m1$fitted.values, m$residuals)
m1 <- aov(shannon ~ Entry, data = input_filt_rare_16S$map_loaded)
summary(m1)
m1 <- aov(shannon ~ Site, data = input_filt_rare_16S$map_loaded)
summary(m1)
TukeyHSD(m1)
t1 <- emmeans(object = m1, specs = "Site") %>%
  cld(object = ., adjust = "Tukey", Letters = letters, alpha = 0.05) %>%
  mutate(name = "shannon",
         y = max(input_filt_rare_16S$map_loaded$shannon)+
           (max(input_filt_rare_16S$map_loaded$shannon)-
              min(input_filt_rare_16S$map_loaded$shannon))/10,
         #y = 5650,
         Dataset = "16S")

label_df_16S <- rbind(t, t1)
alpha_long_16S <- input_filt_rare_16S$map_loaded %>%
  pivot_longer(cols = c("rich", "shannon"))
facet_names <- c("rich" = "ASV Richness",
                 "shannon" = "Shannon Diversity")

pdf("InitialFigs/Alpha_16S.pdf", width = 8, height = 6)
g1 <- ggplot(alpha_long_16S, aes(reorder(Site, value, mean), value)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(size = 2, alpha = 0.75, width = 0.2, pch = 16) +
  geom_text(data = label_df_16S, aes(Site, y, label = str_trim(.group)), 
            size = 3, color = "black") +
  labs(x = "Site", y = "ASV richness") +
  facet_wrap(~ name, scales = "free_y", labeller = as_labeller(facet_names)) +
  theme_bw() +
  theme(legend.position = "right",
        axis.title.y = element_blank(),
        axis.title.x = element_text(size = 12),
        axis.text.y = element_text(size = 10),
        axis.text.x = element_text(size = 10, angle = 45, hjust = 1),
        strip.text = element_text(size = 14),
        strip.background = element_rect(fill = "white"))
g1
dev.off()

# Add insets
m <- aov(rich ~ Site * Entry, data = input_filt_rare_16S$map_loaded)
eta_sq_m <- eta_sq(m) %>%
  as.data.frame() %>%
  rownames_to_column(var = "group") %>%
  set_names(c("group", "value")) %>%
  mutate(group = c("S***", "G", "I", "R")) %>%
  arrange(desc(group)) %>%
  mutate(ypos = cumsum(value) - 0.5*value)
g2 <- ggplot(eta_sq_m, aes(x = "", y = value, fill = group)) +
  geom_bar(stat = "identity", width = 1, color = NA) +
  coord_polar("y", start = 0) +
  geom_text_repel(data = eta_sq_m,
                  aes(y = ypos, label = group), color = "white", size = 3,
                  box.padding = 0.1) +
  scale_fill_manual(values = c("#5DC863FF", "grey30", "grey70", "#21908CFF")) +
  theme_void() + 
  theme(legend.position = "none")
g2

m1 <- aov(shannon ~ Site * Entry, data = input_filt_rare_16S$map_loaded)
eta_sq_m1 <- eta_sq(m1) %>%
  as.data.frame() %>%
  rownames_to_column(var = "group") %>%
  set_names(c("group", "value")) %>%
  mutate(group = c("S***", "G", "I", "R")) %>%
  arrange(desc(group)) %>%
  mutate(ypos = cumsum(value) - 0.5*value)
g3 <- ggplot(eta_sq_m1, aes(x = "", y = value, fill = group)) +
  geom_bar(stat = "identity", width = 1, color = NA) +
  coord_polar("y", start = 0) +
  geom_text_repel(data = eta_sq_m1,
                  aes(y = ypos, label = group), color = "white", size = 3,
                  box.padding = 0.1) +
  scale_fill_manual(values = c("#5DC863FF", "grey30", "grey70", "#21908CFF")) +
  theme_void() + 
  theme(legend.position = "none")
g3

plot.with.inset <-
  ggdraw() +
  draw_plot(g1) +
  draw_plot(g2, x = 0.31, y = 0.25, width = 0.2, height = 0.2) +
  draw_plot(g3, x = 0.8, y = 0.25, width = 0.2, height = 0.2)
plot.with.inset
pdf("InitialFigs/Alpha_wInset_16S.pdf", width = 8, height = 6)
plot.with.inset
dev.off()



#### _ITS ####
# Get descriptive info
min(input_filt_rare_ITS$map_loaded$rich) # 43
max(input_filt_rare_ITS$map_loaded$rich) # 212
mean(input_filt_rare_ITS$map_loaded$rich) # 124
se(input_filt_rare_ITS$map_loaded$rich) # 1

# Test and plot
leveneTest(input_filt_rare_ITS$map_loaded$rich ~ input_filt_rare_ITS$map_loaded$Site) # Not homogeneous
leveneTest(input_filt_rare_ITS$map_loaded$rich ~ input_filt_rare_ITS$map_loaded$Entry) # Homogeneous
m <- aov(rich ~ Site * Entry, data = input_filt_rare_ITS$map_loaded)
summary(m)
Anova(m, type = "III", singular.ok = TRUE) # Site
hist(m$residuals)
shapiro.test(m$residuals)
plot(m$fitted.values, m$residuals)
m <- aov(rich ~ Entry, data = input_filt_rare_ITS$map_loaded)
summary(m)
m <- aov(rich ~ Site, data = input_filt_rare_ITS$map_loaded)
summary(m)
TukeyHSD(m)
t <- emmeans(object = m, specs = "Site") %>%
  cld(object = ., adjust = "Tukey", Letters = letters, alpha = 0.05) %>%
  mutate(name = "rich",
         y = max(input_filt_rare_ITS$map_loaded$rich)+
           (max(input_filt_rare_ITS$map_loaded$rich)-
              min(input_filt_rare_ITS$map_loaded$rich))/10,
         #y = 5650,
         Dataset = "ITS")

leveneTest(input_filt_rare_ITS$map_loaded$shannon ~ input_filt_rare_ITS$map_loaded$Site) # Not homogenous
leveneTest(input_filt_rare_ITS$map_loaded$shannon ~ input_filt_rare_ITS$map_loaded$Entry) # Homogeneous
m1 <- aov(shannon ~ Site * Entry, data = input_filt_rare_ITS$map_loaded)
summary(m1)
Anova(m1, type = "III", singular.ok = TRUE) # Site, Entry, Marginal Interaction
hist(m1$residuals)
shapiro.test(m1$residuals)
plot(m1$fitted.values, m$residuals)
m1 <- aov(shannon ~ Entry, data = input_filt_rare_ITS$map_loaded)
summary(m1)
m1 <- aov(shannon ~ Site, data = input_filt_rare_ITS$map_loaded)
summary(m1)
TukeyHSD(m1)
t1 <- emmeans(object = m1, specs = "Site") %>%
  cld(object = ., adjust = "Tukey", Letters = letters, alpha = 0.05) %>%
  mutate(name = "shannon",
         # y = max(input_filt_rare_ITS$map_loaded$shannon)+
         #   (max(input_filt_rare_ITS$map_loaded$shannon)-
         #      min(input_filt_rare_ITS$map_loaded$shannon))/10,
         y = 5,
         Dataset = "ITS")

label_df_ITS <- rbind(t, t1)
alpha_long_ITS <- input_filt_rare_ITS$map_loaded %>%
  pivot_longer(cols = c("rich", "shannon"))
facet_names <- c("rich" = "ASV Richness",
                 "shannon" = "Shannon Diversity")

pdf("InitialFigs/Alpha_ITS.pdf", width = 8, height = 6)
g4 <- ggplot(alpha_long_ITS, aes(reorder(Site, value, mean), value)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(size = 2, alpha = 0.75, width = 0.2, pch = 16) +
  geom_text(data = label_df_ITS, aes(Site, y, label = str_trim(.group)), 
            size = 3, color = "black") +
  labs(x = "Site", y = "ASV richness") +
  facet_wrap(~ name, scales = "free_y", labeller = as_labeller(facet_names)) +
  theme_bw() +
  theme(legend.position = "right",
        axis.title.y = element_blank(),
        axis.title.x = element_text(size = 12),
        axis.text.y = element_text(size = 10),
        axis.text.x = element_text(size = 10, angle = 45, hjust = 1),
        strip.text = element_text(size = 14),
        strip.background = element_rect(fill = "white"))
g4
dev.off()

# Add insets
m <- aov(rich ~ Site * Entry, data = input_filt_rare_ITS$map_loaded)
eta_sq_m <- eta_sq(m) %>%
  as.data.frame() %>%
  rownames_to_column(var = "group") %>%
  set_names(c("group", "value")) %>%
  mutate(group = c("S***", "G", "I", "R")) %>%
  arrange(desc(group)) %>%
  mutate(ypos = cumsum(value) - 0.5*value)
g5 <- ggplot(eta_sq_m, aes(x = "", y = value, fill = group)) +
  geom_bar(stat = "identity", width = 1, color = NA) +
  coord_polar("y", start = 0) +
  geom_text_repel(data = eta_sq_m,
                  aes(y = ypos, label = group), color = "white", size = 3,
                  box.padding = 0.1) +
  scale_fill_manual(values = c("#5DC863FF", "grey30", "grey70", "#21908CFF")) +
  theme_void() + 
  theme(legend.position = "none")
g5

m1 <- aov(shannon ~ Site * Entry, data = input_filt_rare_ITS$map_loaded)
eta_sq_m1 <- eta_sq(m1) %>%
  as.data.frame() %>%
  rownames_to_column(var = "group") %>%
  set_names(c("group", "value")) %>%
  mutate(group = c("S***", "G", "I", "R")) %>%
  arrange(desc(group)) %>%
  mutate(ypos = cumsum(value) - 0.5*value)
g6 <- ggplot(eta_sq_m1, aes(x = "", y = value, fill = group)) +
  geom_bar(stat = "identity", width = 1, color = NA) +
  coord_polar("y", start = 0) +
  geom_text_repel(data = eta_sq_m1,
                  aes(y = ypos, label = group), color = "white", size = 3,
                  box.padding = 0.1) +
  scale_fill_manual(values = c("#5DC863FF", "grey30", "grey70", "#21908CFF")) +
  theme_void() + 
  theme(legend.position = "none")
g6

plot.with.inset <-
  ggdraw() +
  draw_plot(g4) +
  draw_plot(g5, x = 0.02, y = 0.7, width = 0.2, height = 0.2) +
  draw_plot(g6, x = 0.52, y = 0.7, width = 0.2, height = 0.2)
plot.with.inset
pdf("InitialFigs/Alpha_wInset_ITS.pdf", width = 8, height = 6)
plot.with.inset
dev.off()

# Combined
alpha_long_16S <- alpha_long_16S %>%
  mutate(Dataset = "16S") %>%
  dplyr::select(Site, value, name, Dataset)
alpha_long_ITS <- alpha_long_ITS %>%
  mutate(Dataset = "ITS") %>%
  dplyr::select(Site, value, name, Dataset)
alpha_long <- rbind(alpha_long_16S, alpha_long_ITS)
label_df_long <- rbind(label_df_16S, label_df_ITS)
facet_names <- c("rich" = "ASV Richness",
                 "shannon" = "Shannon Diversity",
                 "16S" = "a) Archaea/Bacteria",
                 "ITS" = "b) Fungi")
g7 <- ggplot(alpha_long, aes(reorder(Site, value, mean), value)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(size = 2, alpha = 0.75, width = 0.2, pch = 16) +
  geom_text(data = label_df_long, aes(Site, y, label = str_trim(.group)), 
            size = 3, color = "black") +
  labs(x = "Site", y = "ASV richness") +
  facet_grid2(name ~ Dataset, scales = "free", independent = "y",
              labeller = as_labeller(facet_names)) +
  theme_bw() +
  theme(legend.position = "right",
        axis.title.y = element_blank(),
        axis.title.x = element_text(size = 12),
        axis.text.y = element_text(size = 10),
        axis.text.x = element_text(size = 10, angle = 45, hjust = 1),
        strip.text = element_text(size = 14),
        strip.background = element_rect(fill = "white"))
g7
plot.with.inset <-
  ggdraw() +
  draw_plot(g7) +
  draw_plot(g2, x = 0.35, y = 0.58, width = 0.15, height = 0.15) +
  draw_plot(g3, x = 0.35, y = 0.21, width = 0.15, height = 0.15) +
  draw_plot(g5, x = 0.51, y = 0.76, width = 0.15, height = 0.15) +
  draw_plot(g6, x = 0.51, y = 0.4, width = 0.15, height = 0.15)
plot.with.inset
pdf("FinalFigs/Figure1.pdf", width = 8, height = 6)
plot.with.inset
dev.off()



#### _Predictors ####
# Which combination of variables best predicts ASV richness?
# Remove chlorophyll and Cd - NAs present
env_16S <- input_filt_rare_16S$map_loaded %>%
  dplyr::select(LAT, LONG, Sclerotinia, CMI, GDD, Temp, P, Tseas, Pseas,
                pH, OM_perc, Salts_mmhos_cm, NO3_ppm, P_ppm, K_ppm, Ca_ppm, Mg_ppm,
                Na_ppm, S_ppm, CEC_meq_g, Sand_perc, Silt_perc, Clay_perc)
n_na <- as.data.frame(matrix(nrow = ncol(env_16S), ncol = 2, NA)) %>%
  set_names(c("Variable", "NA_count"))
for (i in 1:ncol(env_16S)) {
  n_na$Variable[i] <- names(env_16S)[i]
  n_na$NA_count[i] <- sum(is.na(env_16S[,i]))
}
n_na
nrow(env_16S) # 587
sum(is.na(env_16S))
M_16S <- cor(env_16S)
corrplot(M_16S, method = "number", type = "lower")
# Remove highly correlated variables. Only accepts 15 variables.
env_16S <- env_16S %>%
  dplyr::select(LAT, LONG, Sclerotinia, CMI, GDD, P,
                pH, OM_perc, NO3_ppm, P_ppm, K_ppm, Ca_ppm, Mg_ppm,
                Na_ppm, Clay_perc)
y <- input_filt_rare_16S$map_loaded$rich
X <- env_16S
Xy <- as.data.frame(cbind(X,y))
bestmodel <- bestglm(Xy, IC = "AIC", RequireFullEnumerationQ = TRUE, TopModels=10)
bestmodel
bestmodel$BestModels
nullmodel <- lm(rich ~ 1, data = input_filt_rare_16S$map_loaded)
bestmodel <- lm(rich ~ LAT + CMI + GDD + P + pH + OM_perc + NO3_ppm + P_ppm + 
                  Ca_ppm + Mg_ppm + Clay_perc, 
                data = input_filt_rare_16S$map_loaded)
summary(bestmodel)
plot(bestmodel)
AIC(nullmodel)
AIC(bestmodel)
AIC(bestmodel) - AIC(nullmodel)
anova(nullmodel, bestmodel)
# Partial R2
rp_16S <- as.data.frame(rsq.partial(bestmodel, adj = TRUE)) %>%
  dplyr::select(-adjustment)
pdf("InitialFigs/AlphaPredictors_16S.pdf", width = 6, height = 4)
ggplot(rp_16S, aes(reorder(variable, partial.rsq, mean), partial.rsq)) +
  geom_bar(stat = "identity") +
  labs(x = "Predictor",
       y = expression("Partial "*R^2*"")) +
  theme_bw() +
  theme(axis.title = element_text(size = 12, face = "bold"),
        axis.text.y = element_text(size = 10),
        axis.text.x = element_text(size = 10, angle = 45, hjust = 1))
dev.off()



env_ITS <- input_filt_rare_ITS$map_loaded %>%
  dplyr::select(LAT, LONG, Sclerotinia, CMI, GDD, Temp, P, Tseas, Pseas,
                pH, OM_perc, Salts_mmhos_cm, NO3_ppm, P_ppm, K_ppm, Ca_ppm, Mg_ppm,
                Na_ppm, S_ppm, CEC_meq_g, Sand_perc, Silt_perc, Clay_perc)
n_na <- as.data.frame(matrix(nrow = ncol(env_ITS), ncol = 2, NA)) %>%
  set_names(c("Variable", "NA_count"))
for (i in 1:ncol(env_ITS)) {
  n_na$Variable[i] <- names(env_ITS)[i]
  n_na$NA_count[i] <- sum(is.na(env_ITS[,i]))
}
n_na
nrow(env_ITS) # 586
sum(is.na(env_ITS))
M_ITS <- cor(env_ITS)
corrplot(M_ITS, method = "number", type = "lower")
# Remove highly correlated variables. Only accepts 15 variables.
env_ITS <- env_ITS %>%
  dplyr::select(LAT, LONG, Sclerotinia, CMI, GDD, P,
                pH, OM_perc, NO3_ppm, P_ppm, K_ppm, Ca_ppm, Mg_ppm,
                Na_ppm, Clay_perc)
y <- input_filt_rare_ITS$map_loaded$rich
X <- env_ITS
Xy <- as.data.frame(cbind(X,y))
bestmodel <- bestglm(Xy, IC = "AIC", RequireFullEnumerationQ = TRUE, TopModels=10)
bestmodel
bestmodel$BestModels
nullmodel <- lm(rich ~ 1, data = input_filt_rare_ITS$map_loaded)
bestmodel <- lm(rich ~ LAT + LONG + Sclerotinia + CMI + P + pH + OM_perc + 
                  NO3_ppm + P_ppm + Ca_ppm + Clay_perc, 
                data = input_filt_rare_ITS$map_loaded)
summary(bestmodel)
ggplot(input_filt_rare_ITS$map_loaded, aes(P_ppm, rich)) +
  geom_point() +
  geom_smooth(method = "lm")
plot(bestmodel)
AIC(nullmodel)
AIC(bestmodel)
AIC(bestmodel) - AIC(nullmodel)
anova(nullmodel, bestmodel)
# Partial R2
rp_ITS <- as.data.frame(rsq.partial(bestmodel, adj = TRUE)) %>%
  dplyr::select(-adjustment)
pdf("InitialFigs/AlphaPredictors_ITS.pdf", width = 6, height = 4)
ggplot(rp_ITS, aes(reorder(variable, partial.rsq, mean), partial.rsq)) +
  geom_bar(stat = "identity") +
  labs(x = "Predictor",
       y = expression("Partial "*R^2*"")) +
  theme_bw() +
  theme(axis.title = element_text(size = 12, face = "bold"),
        axis.text.y = element_text(size = 10),
        axis.text.x = element_text(size = 10, angle = 45, hjust = 1))
dev.off()

# Combine
rp_16S <- rp_16S %>%
  mutate(Dataset = "16S")
rp_ITS <- rp_ITS %>%
  mutate(Dataset = "ITS")
rp <- rbind(rp_16S, rp_ITS)
facet_names <- c("16S" = "a) Archaea/Bacteria",
                 "ITS" = "b) Fungi")
pdf("FinalFigs/FigureS2.pdf", width = 6, height = 4)
ggplot(rp, aes(reorder(variable, partial.rsq, mean), partial.rsq)) +
  geom_bar(stat = "identity") +
  labs(x = "Predictor",
       y = expression("Partial "*R^2*"")) +
  facet_wrap(~ Dataset, ncol = 1, labeller = as_labeller(facet_names)) +
  theme_bw() +
  theme(axis.title = element_text(size = 12),
        axis.text.y = element_text(size = 10),
        axis.text.x = element_text(size = 10, angle = 45, hjust = 1),
        strip.text = element_text(size = 14),
        strip.background = element_rect(fill = "white"))
dev.off()



#### 4. Beta ####

#### _16S ####

#### __Aitchison ####
# CLR transformation. Do on Innes server as it takes a lot of memory.
# Ran with and without default 80% filter
#asv_czm_16S <- cmultRepl(t(input_filt_16S$data_loaded), label = 0, method = "CZM")
asv_czm_16S <- readRDS("data/asv_czm_16S.rds")
asv_clr_16S <- clr(asv_czm_16S)
aclr_16S <- compositions::dist(asv_clr_16S)

# Filter dropped samples
# Columns and rows containing more than 80% zeros/unobserved values were deleted
dim(t(input_filt_16S$data_loaded)) # 587 samples, 27928 ASVs
dim(asv_czm_16S) # 554 samples, 2454 ASVs. Note: Only highly prevalent ASVs retained!
input_filt_clr_16S <- filter_data(input_filt_16S,
                                  filter_cat = "BoxCell",
                                  keep_vals = rownames(asv_czm_16S))

set.seed(1150)
ado1 <- adonis2(aclr_16S ~ input_filt_clr_16S$map_loaded$Site * input_filt_clr_16S$map_loaded$Entry)
ado1
# Site, Int significant. 0.57, 0.001***; 0.01, 0.075; 0.10, 0.019*
set.seed(1150)
ado2 <- adonis2(aclr_16S ~ input_filt_clr_16S$map_loaded$Entry * input_filt_clr_16S$map_loaded$Site)
ado2
# Both and interaction significant. Small effect of order. 0.57, 0.001***; 0.01, 0.047*; 0.10, 0.019*
anova(betadisper(aclr_16S, input_filt_clr_16S$map_loaded$Site)) # Dispersion not homogeneous
anova(betadisper(aclr_16S, input_filt_clr_16S$map_loaded$Entry)) # Dispersion homogeneous

# PCA with vectors
d.pcx_16S <- prcomp(aclr_16S)
env_clr <- input_filt_clr_16S$map_loaded %>%
  dplyr::select(LAT, LONG, # Geography
                CMI, GDD, Temp, P, Tseas, Pseas, # Climate
                pH, OM_perc, Salts_mmhos_cm, NO3_ppm, P_ppm, K_ppm, Ca_ppm, Mg_ppm, # Soil
                Na_ppm, S_ppm, CEC_meq_g, Sand_perc, Silt_perc, Clay_perc,
                chlorophyll, Sclerotinia, GRIN_Sclerotinia, DaysToFlower, StemLengthHighRange, # Plant
                StemLengthLowRange, PerSeed_14_64, PerSeed_20_64, PerNoBranch, SeedWeight)
set.seed(100)
ef_clr_16S <- envfit(d.pcx_16S, env_clr, permutations = 999, na.rm = TRUE)
ef_clr_16S # All significant but plant
biplot(d.pcx_16S)
plot(ef_clr_16S, p.max = 0.05, cex = 0.5)
ordiplot(d.pcx_16S)
plot(ef_clr_16S, p.max = 0.05, cex = 0.5)
multiplier_clr_16S <- ordiArrowMul(ef_clr_16S)
vec.df_clr_16S <- as.data.frame(ef_clr_16S$vectors$arrows*sqrt(ef_clr_16S$vectors$r)) %>%
  mutate(PC1 = PC1 * multiplier_clr_16S,
         PC2 = PC2 * multiplier_clr_16S) %>%
  mutate(variables = rownames(.)) %>%
  filter(ef_clr_16S$vectors$pvals < 0.05) %>%
  mutate(shortnames = c("LAT", "LONG", "GDD", "Temp", "Precip", "Tseas", "Pseas",
                        "pH", "OM", "Salts", "NO3", "PO4", "K", "Ca", "Mg", "Na", "S",
                        "CEC", "Sand", "Silt", "Clay"))
pcoaA1_clr_16S <- paste("PC1: ", round((eigenvals(d.pcx_16S)/sum(eigenvals(d.pcx_16S)))[1]*100, 1), "%")
pcoaA2_clr_16S <- paste("PC2: ", round((eigenvals(d.pcx_16S)/sum(eigenvals(d.pcx_16S)))[2]*100, 1), "%")
input_filt_clr_16S$map_loaded$Axis01 <- vegan::scores(d.pcx_16S)[,1]
input_filt_clr_16S$map_loaded$Axis02 <- vegan::scores(d.pcx_16S)[,2]
# input_filt_clr_16S$map_loaded$Axis01 <- input_filt_clr_16S$map_loaded$Axis01/sqrt(sum((input_filt_clr_16S$map_loaded$Axis01 - mean(input_filt_clr_16S$map_loaded$Axis01))^2))
# input_filt_clr_16S$map_loaded$Axis02 <- input_filt_clr_16S$map_loaded$Axis02/sqrt(sum((input_filt_clr_16S$map_loaded$Axis02 - mean(input_filt_clr_16S$map_loaded$Axis02))^2))
micro.hulls_clr_16S <- ddply(input_filt_clr_16S$map_loaded, c("Site"), find_hull)
g1 <- ggplot(input_filt_clr_16S$map_loaded, aes(Axis01, Axis02)) +
  geom_polygon(data = micro.hulls_clr_16S, 
               aes(colour = Site),
               alpha = 0, show.legend = F) +
  geom_point(size = 3, pch = 16, alpha = 0.75, aes(colour = Site)) +
  geom_segment(data = vec.df_clr_16S,
               aes(x = 0, xend = PC1, y = 0, yend = PC2),
               arrow = arrow(length = unit(0.35, "cm")),
               colour = "gray", alpha = 0.6,
               inherit.aes = FALSE) +
  geom_text_repel(data = vec.df_clr_16S,
                  aes(x = PC1, y = PC2, label = shortnames),
                  size = 3, color = "black") +
  labs(x = pcoaA1_clr_16S, 
       y = pcoaA2_clr_16S,
       colour = "Site") +
  scale_colour_viridis_d() +
  theme_bw() +  
  theme(legend.position = "right",
        axis.title = element_text(face = "bold", size = 12), 
        axis.text = element_text(size = 10))
g1
pdf("InitialFigs/Beta_Aitch_16S.pdf", width = 7, height = 5)
g1
dev.off()

# Add inset
eta_sq_m2 <- eta_sq_adonis(ado1) %>%
  as.data.frame() %>%
  rownames_to_column(var = "group") %>%
  filter(group != "Total") %>%
  set_names(c("group", "value")) %>%
  mutate(group = c("Site***", "Gen.", "Int.*", "Resid.")) %>%
  arrange(desc(group)) %>%
  mutate(ypos = cumsum(value) - 0.5*value)
g2 <- ggplot(eta_sq_m2, aes(x = "", y = value, fill = group)) +
  geom_bar(stat = "identity", width = 1, color = NA) +
  coord_polar("y", start = 0) +
  geom_text_repel(data = eta_sq_m2,
                  aes(y = ypos, label = group), color = "white", size = 3,
                  box.padding = 0.1) +
  scale_fill_manual(values = c("#5DC863FF", "grey30", "grey70", "#21908CFF")) +
  theme_void() + 
  theme(legend.position = "none")
g2

plot.with.inset2 <-
  ggdraw() +
  draw_plot(g1) +
  draw_plot(g2, x = 0.4, y = 0.1, width = 0.25, height = 0.25)
plot.with.inset2
pdf("InitialFigs/Beta_Aitch_wInset_16S.pdf", width = 7, height = 5)
plot.with.inset2
dev.off()



# Rerun with all samples and taxa
asv_czm_16S <- readRDS("data/asv_czm_keep_16S.rds")
asv_clr_16S <- clr(asv_czm_16S)
aclr_16S <- compositions::dist(asv_clr_16S)

set.seed(1150)
ado3 <- adonis2(aclr_16S ~ input_filt_16S$map_loaded$Site * input_filt_16S$map_loaded$Entry)
ado3
# Site and interaction significant. 0.48, 0.001***; 0.01, 0.132; 0.12, 0.004**
set.seed(1150)
ado4 <- adonis2(aclr_16S ~ input_filt_16S$map_loaded$Entry * input_filt_16S$map_loaded$Site)
ado4
# Site and interaction significant. No effect of order. 0.48, 0.001***; 0.01, 0.105; 0.12, 0.004**
anova(betadisper(aclr_16S, input_filt_16S$map_loaded$Site)) # Dispersion not homogeneous
anova(betadisper(aclr_16S, input_filt_16S$map_loaded$Entry)) # Dispersion homogeneous

# PCA with vectors
d.pcx_16S <- prcomp(aclr_16S)
env_clr <- input_filt_16S$map_loaded %>%
  dplyr::select(LAT, LONG, # Geography
                CMI, GDD, Temp, P, Tseas, Pseas, # Climate
                pH, OM_perc, Salts_mmhos_cm, NO3_ppm, P_ppm, K_ppm, Ca_ppm, Mg_ppm, # Soil
                Na_ppm, S_ppm, CEC_meq_g, Sand_perc, Silt_perc, Clay_perc,
                chlorophyll, Sclerotinia, GRIN_Sclerotinia, DaysToFlower, StemLengthHighRange, # Plant
                StemLengthLowRange, PerSeed_14_64, PerSeed_20_64, PerNoBranch, SeedWeight)
set.seed(100)
ef_clr_16S <- envfit(d.pcx_16S, env_clr, permutations = 999, na.rm = TRUE)
ef_clr_16S # Same but chlorophyll sig and silt not
biplot(d.pcx_16S)
plot(ef_clr_16S, p.max = 0.05, cex = 0.5)
ordiplot(d.pcx_16S)
plot(ef_clr_16S, p.max = 0.05, cex = 0.5)
multiplier_16S <- ordiArrowMul(ef_clr_16S)
vec.df_16S <- as.data.frame(ef_clr_16S$vectors$arrows*sqrt(ef_clr_16S$vectors$r)) %>%
  mutate(PC1 = PC1 * multiplier_16S,
         PC2 = PC2 * multiplier_16S) %>%
  mutate(variables = rownames(.)) %>%
  filter(ef_clr_16S$vectors$pvals < 0.05) %>%
  mutate(shortnames = c("LAT", "LONG", "CMI", "GDD", "Temp", "Precip", "Tseas", "Pseas",
                        "pH", "OM", "Salts", "NO3", "PO4", "K", "Ca", "Mg", "Na", "S",
                        "CEC", "Clay", "chlorophyll"))
pcoaA1_16S <- paste("PC1: ", round((eigenvals(d.pcx_16S)/sum(eigenvals(d.pcx_16S)))[1]*100, 1), "%")
pcoaA2_16S <- paste("PC2: ", round((eigenvals(d.pcx_16S)/sum(eigenvals(d.pcx_16S)))[2]*100, 1), "%")
input_filt_16S$map_loaded$Axis01 <- vegan::scores(d.pcx_16S)[,1]
input_filt_16S$map_loaded$Axis02 <- vegan::scores(d.pcx_16S)[,2]
# input_filt_16S$map_loaded$Axis01 <- input_filt_16S$map_loaded$Axis01/sqrt(sum((input_filt_16S$map_loaded$Axis01 - mean(input_filt_16S$map_loaded$Axis01))^2))
# input_filt_16S$map_loaded$Axis02 <- input_filt_16S$map_loaded$Axis02/sqrt(sum((input_filt_16S$map_loaded$Axis02 - mean(input_filt_16S$map_loaded$Axis02))^2))
micro.hulls_16S <- ddply(input_filt_16S$map_loaded, c("Site"), find_hull)
g1 <- ggplot(input_filt_16S$map_loaded, aes(Axis01, Axis02)) +
  geom_polygon(data = micro.hulls_16S, 
               aes(colour = Site),
               alpha = 0, show.legend = F) +
  geom_point(size = 3, pch = 16, alpha = 0.75, aes(colour = Site)) +
  geom_segment(data = vec.df_16S,
               aes(x = 0, xend = PC1, y = 0, yend = PC2),
               arrow = arrow(length = unit(0.35, "cm")),
               colour = "gray", alpha = 0.6,
               inherit.aes = FALSE) +
  geom_text_repel(data = vec.df_16S,
                  aes(x = PC1, y = PC2, label = shortnames),
                  size = 3, color = "black", max.overlaps = 20) +
  labs(x = pcoaA1_16S, 
       y = pcoaA2_16S,
       colour = "Site") +
  scale_colour_viridis_d() +
  theme_bw() +  
  theme(legend.position = "right",
        axis.title = element_text(face = "bold", size = 12), 
        axis.text = element_text(size = 10))
g1

pdf("InitialFigs/Beta_Aitch_All_16S.pdf", width = 7, height = 5)
g1
dev.off()



#### __Bray ####
bc_16S <- calc_dm(input_filt_rare_16S$data_loaded)
set.seed(1150)
ado5 <- adonis2(bc_16S ~ input_filt_rare_16S$map_loaded$Site * input_filt_rare_16S$map_loaded$Entry)
ado5
# Both and interaction significant. 0.56, 0.001; 0.01, 0.042; 0.10, 0.001
set.seed(1150)
ado6 <- adonis2(bc_16S ~ input_filt_rare_16S$map_loaded$Entry * input_filt_rare_16S$map_loaded$Site)
ado6
# Both and interaction significant. Small effect of order. 0.56, 0.001; 0.01, 0.027; 0.10, 0.001
anova(betadisper(bc_16S, input_filt_rare_16S$map_loaded$Site)) # Dispersion not homogeneous
anova(betadisper(bc_16S, input_filt_rare_16S$map_loaded$Entry)) # Dispersion homogeneous

# Test effect of genotype while stratifying permutations by Site
set.seed(1150)
adonis2(bc_16S ~ input_filt_rare_16S$map_loaded$Entry, 
        strata = input_filt_rare_16S$map_loaded$Site)

# Pairwise by Site
set.seed(1150)
pwa_16S <- as.data.frame(pairwiseAdonis::pairwise.adonis(bc_16S, 
                                                         input_filt_rare_16S$map_loaded$Site))

# PCoA with vectors
pcoa_16S <- cmdscale(bc_16S, k = nrow(input_filt_rare_16S$map_loaded) - 1, eig = T)
env_16S <- input_filt_rare_16S$map_loaded %>%
  dplyr::select(LAT, LONG, # Geography
                CMI, GDD, Temp, P, Tseas, Pseas, # Climate
                pH, OM_perc, Salts_mmhos_cm, NO3_ppm, P_ppm, K_ppm, Ca_ppm, Mg_ppm, # Soil
                Na_ppm, S_ppm, CEC_meq_g, Sand_perc, Silt_perc, Clay_perc,
                chlorophyll, Sclerotinia, GRIN_Sclerotinia, DaysToFlower, StemLengthHighRange, # Plant
                StemLengthLowRange, PerSeed_14_64, PerSeed_20_64, PerNoBranch, SeedWeight)
set.seed(100)
ef_16S <- envfit(pcoa_16S, env_16S, permutations = 999, na.rm = TRUE)
ef_16S # All sig but plant variables (except chlorophyll)
ordiplot(pcoa_16S)
plot(ef_16S, p.max = 0.05, cex = 0.5)
multiplier_16S <- ordiArrowMul(ef_16S)
vec.df_16S <- as.data.frame(ef_16S$vectors$arrows*sqrt(ef_16S$vectors$r)) %>%
  mutate(Dim1 = Dim1 * multiplier_16S,
         Dim2 = Dim2 * multiplier_16S) %>%
  mutate(variables = rownames(.)) %>%
  filter(ef_16S$vectors$pvals < 0.05) %>%
  mutate(shortnames = c("LAT", "LONG", "CMI", "GDD", "Temp", "Precip", "Tseas", "Pseas",
                        "pH", "OM", "Salts", "NO3", "PO4", "K", "Ca", "Mg", "Na", "S",
                        "CEC", "Sand", "Silt", "Clay", "Chlorophyll"))
pcoaA1_16S <- paste("PC1: ", round((eigenvals(pcoa_16S)/sum(eigenvals(pcoa_16S)))[1]*100, 1), "%")
pcoaA2_16S <- paste("PC2: ", round((eigenvals(pcoa_16S)/sum(eigenvals(pcoa_16S)))[2]*100, 1), "%")
input_filt_rare_16S$map_loaded$Axis01 <- vegan::scores(pcoa_16S)[,1]
input_filt_rare_16S$map_loaded$Axis02 <- vegan::scores(pcoa_16S)[,2]
micro.hulls_16S <- ddply(input_filt_rare_16S$map_loaded, c("Site"), find_hull)
g8 <- ggplot(input_filt_rare_16S$map_loaded, aes(Axis01, Axis02)) +
  geom_polygon(data = micro.hulls_16S, 
               aes(colour = Site),
               alpha = 0, show.legend = F) +
  geom_point(size = 3, pch = 16, alpha = 0.75, aes(colour = Site)) +
  geom_segment(data = vec.df_16S,
               aes(x = 0, xend = Dim1, y = 0, yend = Dim2),
               arrow = arrow(length = unit(0.35, "cm")),
               colour = "gray", alpha = 0.6,
               inherit.aes = FALSE) + 
  geom_text_repel(data = vec.df_16S,
                  aes(x = Dim1, y = Dim2, label = shortnames),
                  size = 3, color = "black") +
  labs(x = pcoaA1_16S, 
       y = pcoaA2_16S,
       colour = "Site") +
  scale_colour_viridis_d() +
  theme_bw() +  
  theme(legend.position = "right",
        axis.title = element_text(face = "bold", size = 12), 
        axis.text = element_text(size = 10))
g8
pdf("InitialFigs/Beta_Bray_16S.pdf", width = 7, height = 5)
g8
dev.off()

leg <- get_legend((g8))

# Add inset
eta_sq_m3 <- eta_sq_adonis(ado5) %>%
  as.data.frame() %>%
  rownames_to_column(var = "group") %>%
  filter(group != "Total") %>%
  set_names(c("group", "value")) %>%
  mutate(group = c("S***", "G*", "I***", "R")) %>%
  arrange(desc(group)) %>%
  mutate(ypos = cumsum(value) - 0.5*value)
g9 <- ggplot(eta_sq_m3, aes(x = "", y = value, fill = group)) +
  geom_bar(stat = "identity", width = 1, color = NA) +
  coord_polar("y", start = 0) +
  geom_text_repel(data = eta_sq_m3,
                  aes(y = ypos, label = group), color = "white", size = 3,
                  box.padding = 0.1) +
  scale_fill_manual(values = c("#5DC863FF", "grey30", "grey70", "#21908CFF")) +
  theme_void() + 
  theme(legend.position = "none")
g9

plot.with.inset3 <-
  ggdraw() +
  draw_plot(g4) +
  draw_plot(g5, x = 0.54, y = 0.08, width = 0.25, height = 0.25)
plot.with.inset3
pdf("InitialFigs/Beta_Bray_wInset_16S.pdf", width = 7, height = 5)
plot.with.inset3
dev.off()

# Make ordinations with site centroids and labels
# Also flip axes to better match map
pd <- betadisper(bc_16S, input_filt_rare_16S$map_loaded$Site)
centroids_16S <- as.data.frame(pd$centroids) %>%
  dplyr::select(PCoA1, PCoA2) %>%
  rownames_to_column(var = "Site")
pcoa_16S <- cmdscale(bc_16S, k = nrow(input_filt_rare_16S$map_loaded) - 1, eig = T)
pcoaA1 <- paste("PC1: ", round((eigenvals(pcoa_16S)/sum(eigenvals(pcoa_16S)))[1]*100, 1), "%")
pcoaA2 <- paste("PC2: ", round((eigenvals(pcoa_16S)/sum(eigenvals(pcoa_16S)))[2]*100, 1), "%")
input_filt_rare_16S$map_loaded$Axis01 <- vegan::scores(pcoa_16S)[,1]
input_filt_rare_16S$map_loaded$Axis02 <- vegan::scores(pcoa_16S)[,2]
micro.hulls <- ddply(input_filt_rare_16S$map_loaded, c("Site"), find_hull)
pdf("InitialFigs/Beta_Bray_Centroids_16S.pdf", width = 7, height = 5)
ggplot(input_filt_rare_16S$map_loaded, aes(-Axis01, -Axis02)) +
  stat_ellipse(mapping = aes(x = -Axis01, y = -Axis02, color = Site, fill = Site), 
               geom = "polygon", alpha = 0.2, linewidth = NA) +
  geom_point(size = 1.5, pch = 16, alpha = 0.375, aes(colour = Site)) +
  geom_point(data = centroids_16S, aes(-PCoA1, -PCoA2, fill = Site),
             size = 5, shape = 24, stroke = 1) +
  geom_text_repel(data = centroids_16S, aes(-PCoA1, -PCoA2, label = Site),
                  box.padding = 0.7) +
  labs(x = pcoaA1, 
       y = pcoaA2,
       colour = "Site") +
  scale_colour_viridis_d() +
  scale_fill_viridis_d() +
  theme_bw() +  
  theme(legend.position = "none",
        axis.title = element_text(face = "bold", size = 12), 
        axis.text = element_text(size = 10),
        plot.title = element_text(vjust = -1))
dev.off()

# Add vectors
pdf("InitialFigs/Beta_Bray_Centroids_wVec_16S.pdf", width = 7, height = 5)
g10 <- ggplot(input_filt_rare_16S$map_loaded, aes(-Axis01, -Axis02)) +
  stat_ellipse(mapping = aes(x = -Axis01, y = -Axis02, color = Site, fill = Site), 
               geom = "polygon", alpha = 0.2, linewidth = NA) +
  geom_point(size = 1.5, pch = 16, alpha = 0.375, aes(colour = Site)) +
  geom_point(data = centroids_16S, aes(-PCoA1, -PCoA2, fill = Site),
             size = 5, shape = 24, stroke = 1) +
  geom_segment(data = vec.df_16S,
               aes(x = 0, xend = -Dim1, y = 0, yend = -Dim2),
               arrow = arrow(length = unit(0.35, "cm")),
               colour = "gray", alpha = 0.25,
               inherit.aes = FALSE) + 
  geom_text_repel(data = vec.df_16S,
                  aes(x = -Dim1, y = -Dim2, label = shortnames),
                  size = 2, color = "red", box.padding = 0.05) +
  # geom_text(data = subset(vec.df_16S, shortnames %in% c("LAT")),
  #           aes(x = -Dim1, y = -Dim2, label = shortnames),
  #           size = 2, color = "red", nudge_y = 0.03) +
  # geom_text(data = subset(vec.df_16S, shortnames %in% c("Precip", "CEC")),
  #           aes(x = -Dim1, y = -Dim2, label = shortnames),
  #           size = 2, color = "red", nudge_y = 0.01) +
  # geom_text(data = subset(vec.df_16S, shortnames %in% c("P04")),
  #           aes(x = -Dim1, y = -Dim2, label = shortnames),
  #           size = 2, color = "red", nudge_y = -0.01) +
  geom_text_repel(data = centroids_16S, aes(-PCoA1, -PCoA2, label = Site),
                  box.padding = 0.7, size = 3) +
  labs(x = pcoaA1, 
       y = pcoaA2,
       colour = "Site") +
  scale_colour_viridis_d() +
  scale_fill_viridis_d() +
  theme_bw() +  
  theme(legend.position = "none",
        axis.title = element_text(face = "bold", size = 12), 
        axis.text = element_text(size = 10),
        plot.title = element_text(vjust = -1),
        panel.grid = element_blank())
g10
dev.off()

# Check ordination at all taxonomic levels. Need to see if it's consistent with ASV level.
tax_sum_phyla <- summarize_taxonomy(input = input_filt_rare_16S, 
                                    level = 2, 
                                    report_higher_tax = T)
phy_bc <- calc_dm(tax_sum_phyla)
pcoa_16S <- cmdscale(phy_bc, k = nrow(input_filt_rare_16S$map_loaded) - 1, eig = T)
pcoaA1 <- paste("PC1: ", round((eigenvals(pcoa_16S)/sum(eigenvals(pcoa_16S)))[1]*100, 1), "%")
pcoaA2 <- paste("PC2: ", round((eigenvals(pcoa_16S)/sum(eigenvals(pcoa_16S)))[2]*100, 1), "%")
input_filt_rare_16S$map_loaded$Axis01 <- vegan::scores(pcoa_16S)[,1]
input_filt_rare_16S$map_loaded$Axis02 <- vegan::scores(pcoa_16S)[,2]
micro.hulls <- ddply(input_filt_rare_16S$map_loaded, c("Site"), find_hull)
g_p <- ggplot(input_filt_rare_16S$map_loaded, aes(Axis01, Axis02)) +
  geom_polygon(data = micro.hulls, 
               aes(colour = Site),
               alpha = 0, show.legend = F) +
  geom_point(size = 2, pch = 16, alpha = 0.75, aes(colour = Site)) +
  labs(x = pcoaA1, 
       y = pcoaA2,
       colour = "Site") +
  scale_colour_viridis_d() +
  ggtitle("Phylum") +
  theme_bw() +  
  theme(legend.position = "none",
        axis.title = element_text(face = "bold", size = 12), 
        axis.text = element_text(size = 10),
        plot.title = element_text(vjust = -1))
g_p

tax_sum_classes <- summarize_taxonomy(input = input_filt_rare_16S, 
                                      level = 3, 
                                      report_higher_tax = T)
cla_bc <- calc_dm(tax_sum_classes)
pcoa_16S <- cmdscale(cla_bc, k = nrow(input_filt_rare_16S$map_loaded) - 1, eig = T)
pcoaA1 <- paste("PC1: ", round((eigenvals(pcoa_16S)/sum(eigenvals(pcoa_16S)))[1]*100, 1), "%")
pcoaA2 <- paste("PC2: ", round((eigenvals(pcoa_16S)/sum(eigenvals(pcoa_16S)))[2]*100, 1), "%")
input_filt_rare_16S$map_loaded$Axis01 <- vegan::scores(pcoa_16S)[,1]
input_filt_rare_16S$map_loaded$Axis02 <- vegan::scores(pcoa_16S)[,2]
micro.hulls <- ddply(input_filt_rare_16S$map_loaded, c("Site"), find_hull)
g_c <- ggplot(input_filt_rare_16S$map_loaded, aes(Axis01, -Axis02)) +
  geom_polygon(data = micro.hulls, 
               aes(colour = Site),
               alpha = 0, show.legend = F) +
  geom_point(size = 2, pch = 16, alpha = 0.75, aes(colour = Site)) +
  labs(x = pcoaA1, 
       y = pcoaA2,
       colour = "Site") +
  scale_colour_viridis_d() +
  ggtitle("Class") +
  theme_bw() +  
  theme(legend.position = "none",
        axis.title = element_text(face = "bold", size = 12), 
        axis.text = element_text(size = 10),
        plot.title = element_text(vjust = -1))
g_c

tax_sum_orders <- summarize_taxonomy(input = input_filt_rare_16S, 
                                     level = 4, 
                                     report_higher_tax = T)
ord_bc <- calc_dm(tax_sum_orders)
pcoa_16S <- cmdscale(ord_bc, k = nrow(input_filt_rare_16S$map_loaded) - 1, eig = T)
pcoaA1 <- paste("PC1: ", round((eigenvals(pcoa_16S)/sum(eigenvals(pcoa_16S)))[1]*100, 1), "%")
pcoaA2 <- paste("PC2: ", round((eigenvals(pcoa_16S)/sum(eigenvals(pcoa_16S)))[2]*100, 1), "%")
input_filt_rare_16S$map_loaded$Axis01 <- vegan::scores(pcoa_16S)[,1]
input_filt_rare_16S$map_loaded$Axis02 <- vegan::scores(pcoa_16S)[,2]
micro.hulls <- ddply(input_filt_rare_16S$map_loaded, c("Site"), find_hull)
g_o <- ggplot(input_filt_rare_16S$map_loaded, aes(Axis01, Axis02)) +
  geom_polygon(data = micro.hulls, 
               aes(colour = Site),
               alpha = 0, show.legend = F) +
  geom_point(size = 2, pch = 16, alpha = 0.75, aes(colour = Site)) +
  labs(x = pcoaA1, 
       y = pcoaA2,
       colour = "Site") +
  scale_colour_viridis_d() +
  ggtitle("Order") +
  theme_bw() +  
  theme(legend.position = "none",
        axis.title = element_text(face = "bold", size = 12), 
        axis.text = element_text(size = 10),
        plot.title = element_text(vjust = -1))
g_o

tax_sum_families <- summarize_taxonomy(input = input_filt_rare_16S, 
                                       level = 5, 
                                       report_higher_tax = T)
fam_bc <- calc_dm(tax_sum_families)
pcoa_16S <- cmdscale(fam_bc, k = nrow(input_filt_rare_16S$map_loaded) - 1, eig = T)
pcoaA1 <- paste("PC1: ", round((eigenvals(pcoa_16S)/sum(eigenvals(pcoa_16S)))[1]*100, 1), "%")
pcoaA2 <- paste("PC2: ", round((eigenvals(pcoa_16S)/sum(eigenvals(pcoa_16S)))[2]*100, 1), "%")
input_filt_rare_16S$map_loaded$Axis01 <- vegan::scores(pcoa_16S)[,1]
input_filt_rare_16S$map_loaded$Axis02 <- vegan::scores(pcoa_16S)[,2]
micro.hulls <- ddply(input_filt_rare_16S$map_loaded, c("Site"), find_hull)
g_f <- ggplot(input_filt_rare_16S$map_loaded, aes(Axis01, -Axis02)) +
  geom_polygon(data = micro.hulls, 
               aes(colour = Site),
               alpha = 0, show.legend = F) +
  geom_point(size = 2, pch = 16, alpha = 0.75, aes(colour = Site)) +
  labs(x = pcoaA1, 
       y = pcoaA2,
       colour = "Site") +
  scale_colour_viridis_d() +
  ggtitle("Family") +
  theme_bw() +  
  theme(legend.position = "none",
        axis.title = element_text(face = "bold", size = 12), 
        axis.text = element_text(size = 10),
        plot.title = element_text(vjust = -1))
g_f

tax_sum_genera <- summarize_taxonomy(input = input_filt_rare_16S, 
                                     level = 6, 
                                     report_higher_tax = T)
gen_bc <- calc_dm(tax_sum_genera)
pcoa_16S <- cmdscale(gen_bc, k = nrow(input_filt_rare_16S$map_loaded) - 1, eig = T)
pcoaA1 <- paste("PC1: ", round((eigenvals(pcoa_16S)/sum(eigenvals(pcoa_16S)))[1]*100, 1), "%")
pcoaA2 <- paste("PC2: ", round((eigenvals(pcoa_16S)/sum(eigenvals(pcoa_16S)))[2]*100, 1), "%")
input_filt_rare_16S$map_loaded$Axis01 <- vegan::scores(pcoa_16S)[,1]
input_filt_rare_16S$map_loaded$Axis02 <- vegan::scores(pcoa_16S)[,2]
micro.hulls <- ddply(input_filt_rare_16S$map_loaded, c("Site"), find_hull)
g_g <- ggplot(input_filt_rare_16S$map_loaded, aes(Axis01, Axis02)) +
  geom_polygon(data = micro.hulls, 
               aes(colour = Site),
               alpha = 0, show.legend = F) +
  geom_point(size = 2, pch = 16, alpha = 0.75, aes(colour = Site)) +
  labs(x = pcoaA1, 
       y = pcoaA2,
       colour = "Site") +
  scale_colour_viridis_d() +
  ggtitle("Genus") +
  theme_bw() +  
  theme(legend.position = "none",
        axis.title = element_text(face = "bold", size = 12), 
        axis.text = element_text(size = 10),
        plot.title = element_text(vjust = -1))
g_g

pcoa_16S <- cmdscale(bc_16S, k = nrow(input_filt_rare_16S$map_loaded) - 1, eig = T)
pcoaA1 <- paste("PC1: ", round((eigenvals(pcoa_16S)/sum(eigenvals(pcoa_16S)))[1]*100, 1), "%")
pcoaA2 <- paste("PC2: ", round((eigenvals(pcoa_16S)/sum(eigenvals(pcoa_16S)))[2]*100, 1), "%")
input_filt_rare_16S$map_loaded$Axis01 <- vegan::scores(pcoa_16S)[,1]
input_filt_rare_16S$map_loaded$Axis02 <- vegan::scores(pcoa_16S)[,2]
micro.hulls <- ddply(input_filt_rare_16S$map_loaded, c("Site"), find_hull)
g_asv <- ggplot(input_filt_rare_16S$map_loaded, aes(Axis01, Axis02)) +
  geom_polygon(data = micro.hulls, 
               aes(colour = Site),
               alpha = 0, show.legend = F) +
  geom_point(size = 2, pch = 16, alpha = 0.75, aes(colour = Site)) +
  labs(x = pcoaA1, 
       y = pcoaA2,
       colour = "Site") +
  scale_colour_viridis_d() +
  ggtitle("ASV") +
  theme_bw() +  
  theme(legend.position = "none",
        axis.title = element_text(face = "bold", size = 12), 
        axis.text = element_text(size = 10),
        plot.title = element_text(vjust = -1))
g_asv

ggplotly(g_g) # Ralls and Kirkmeyer cluster. Kind of following lat/long.
mp1 <- plot_grid(g_p, g_c, g_o, g_f, g_g, g_asv, ncol = 3)
plot_grid(mp1, leg, rel_widths = c(0.85, 0.15))

pdf("InitialFigs/Beta_Bray_AllLevels_16S.pdf", width = 9, height = 6)
plot_grid(mp1, leg, rel_widths = c(0.83, 0.17))
dev.off()

# Check (plot, test) Archaea and Bacteria individually
arc <- filter_taxa_from_input(input_filt_rare_16S,
                              taxa_to_keep = "Archaea",
                              at_spec_level = 1)
nrow(arc$data_loaded) # 479 Archaea
arc_bc <- calc_dm(arc$data_loaded)
set.seed(1000)
adonis2(arc_bc ~ arc$map_loaded$Site * arc$map_loaded$Entry) # Both sig and int
set.seed(1000)
adonis2(arc_bc ~ arc$map_loaded$Entry * arc$map_loaded$Site) # No effect of order
anova(betadisper(arc_bc, arc$map_loaded$Site)) # Dispersion not homogeneous
anova(betadisper(arc_bc, arc$map_loaded$Entry)) # Dispersion homogeneous
arc_pcoa <- cmdscale(arc_bc, k = nrow(arc$map_loaded) - 1, eig = T)
pcoaA1 <- paste("PC1: ", round((eigenvals(arc_pcoa)/sum(eigenvals(arc_pcoa)))[1]*100, 1), "%")
pcoaA2 <- paste("PC2: ", round((eigenvals(arc_pcoa)/sum(eigenvals(arc_pcoa)))[2]*100, 1), "%")
arc$map_loaded$Axis01 <- vegan::scores(arc_pcoa)[,1]
arc$map_loaded$Axis02 <- vegan::scores(arc_pcoa)[,2]
micro.hulls <- ddply(arc$map_loaded, c("Site"), find_hull)
pdf("InitialFigs/Beta_Bray_Arc.pdf", width = 8, height = 6)
ggplot(arc$map_loaded, aes(Axis01, Axis02)) +
  geom_polygon(data = micro.hulls, 
               aes(colour = Site),
               alpha = 0, show.legend = F) +
  geom_point(size = 3, pch = 16, alpha = 0.75, aes(colour = Site)) +
  labs(x = pcoaA1, 
       y = pcoaA2,
       colour = "Site") +
  scale_colour_viridis_d() +
  ggtitle("Archaea") +
  theme_bw() +  
  theme(legend.position = "right",
        axis.title = element_text(face = "bold", size = 12), 
        axis.text = element_text(size = 10),
        plot.title = element_text(vjust = -1))
dev.off()

bac <- filter_taxa_from_input(input_filt_rare_16S,
                              taxa_to_keep = "Bacteria",
                              at_spec_level = 1)
nrow(bac$data_loaded) # 27414 Bacteria
bac_bc <- calc_dm(bac$data_loaded)
set.seed(1000)
adonis2(bac_bc ~ bac$map_loaded$Site * bac$map_loaded$Entry) # Both sig and int
set.seed(1000)
adonis2(bac_bc ~ bac$map_loaded$Entry * bac$map_loaded$Site) # No effect of order
anova(betadisper(bac_bc, bac$map_loaded$Site)) # Dispersion not homogeneous
anova(betadisper(bac_bc, bac$map_loaded$Entry)) # Dispersion homogeneous
bac_pcoa <- cmdscale(bac_bc, k = nrow(bac$map_loaded) - 1, eig = T)
pcoaA1 <- paste("PC1: ", round((eigenvals(bac_pcoa)/sum(eigenvals(bac_pcoa)))[1]*100, 1), "%")
pcoaA2 <- paste("PC2: ", round((eigenvals(bac_pcoa)/sum(eigenvals(bac_pcoa)))[2]*100, 1), "%")
bac$map_loaded$Axis01 <- vegan::scores(bac_pcoa)[,1]
bac$map_loaded$Axis02 <- vegan::scores(bac_pcoa)[,2]
micro.hulls <- ddply(bac$map_loaded, c("Site"), find_hull)
pdf("InitialFigs/Beta_Bray_Bac.pdf", width = 8, height = 6)
ggplot(bac$map_loaded, aes(Axis01, Axis02)) +
  geom_polygon(data = micro.hulls, 
               aes(colour = Site),
               alpha = 0, show.legend = F) +
  geom_point(size = 3, pch = 16, alpha = 0.75, aes(colour = Site)) +
  labs(x = pcoaA1, 
       y = pcoaA2,
       colour = "Site") +
  scale_colour_viridis_d() +
  ggtitle("Bacteria") +
  theme_bw() +  
  theme(legend.position = "right",
        axis.title = element_text(face = "bold", size = 12), 
        axis.text = element_text(size = 10),
        plot.title = element_text(vjust = -1))
dev.off()



#### __Jaccard ####
ja_16S <- calc_dm(input_filt_rare_16S$data_loaded, method = "jaccard")
set.seed(1150)
ado7 <- adonis2(ja_16S ~ input_filt_rare_16S$map_loaded$Site * input_filt_rare_16S$map_loaded$Entry)
ado7 # Site and interaction significant.
set.seed(1150)
ado8 <- adonis2(ja_16S ~ input_filt_rare_16S$map_loaded$Entry * input_filt_rare_16S$map_loaded$Site)
ado8
# Site and interaction significant. Small effect of order (genotype significant 0.045*)
anova(betadisper(ja_16S, input_filt_rare_16S$map_loaded$Site)) # Dispersion not homogeneous
anova(betadisper(ja_16S, input_filt_rare_16S$map_loaded$Entry)) # Dispersion homogeneous

# PCoA with vectors
pcoa_16S <- cmdscale(ja_16S, k = nrow(input_filt_rare_16S$map_loaded) - 1, eig = T)
set.seed(100)
ef_16S <- envfit(pcoa_16S, env_16S, permutations = 999, na.rm = TRUE)
ef_16S # All sig but plant variables (except chlorophyll)
ordiplot(pcoa_16S)
plot(ef_16S, p.max = 0.05, cex = 0.5)
multiplier_16S <- ordiArrowMul(ef_16S)
vec.df_16S <- as.data.frame(ef_16S$vectors$arrows*sqrt(ef_16S$vectors$r)) %>%
  mutate(Dim1 = Dim1 * multiplier_16S,
         Dim2 = Dim2 * multiplier_16S) %>%
  mutate(variables = rownames(.)) %>%
  filter(ef_16S$vectors$pvals < 0.05) %>%
  mutate(shortnames = c("LAT", "LONG", "CMI", "GDD", "Temp", "Precip", "Tseas", "Pseas",
                        "pH", "OM", "Salts", "NO3", "PO4", "K", "Ca", "Mg", "Na", "S",
                        "CEC", "Sand", "Silt", "Clay", "Chlorophyll"))
pcoaA1_16S <- paste("PC1: ", round((eigenvals(pcoa_16S)/sum(eigenvals(pcoa_16S)))[1]*100, 1), "%")
pcoaA2_16S <- paste("PC2: ", round((eigenvals(pcoa_16S)/sum(eigenvals(pcoa_16S)))[2]*100, 1), "%")
input_filt_rare_16S$map_loaded$Axis01 <- vegan::scores(pcoa_16S)[,1]
input_filt_rare_16S$map_loaded$Axis02 <- vegan::scores(pcoa_16S)[,2]
micro.hulls_16S <- ddply(input_filt_rare_16S$map_loaded, c("Site"), find_hull)
g4 <- ggplot(input_filt_rare_16S$map_loaded, aes(Axis01, Axis02)) +
  geom_polygon(data = micro.hulls_16S, 
               aes(colour = Site),
               alpha = 0, show.legend = F) +
  geom_point(size = 3, pch = 16, alpha = 0.75, aes(colour = Site)) +
  geom_segment(data = vec.df_16S,
               aes(x = 0, xend = Dim1, y = 0, yend = Dim2),
               arrow = arrow(length = unit(0.35, "cm")),
               colour = "gray", alpha = 0.6,
               inherit.aes = FALSE) + 
  geom_text_repel(data = vec.df_16S,
                  aes(x = Dim1, y = Dim2, label = shortnames),
                  size = 3, color = "black") +
  labs(x = pcoaA1_16S, 
       y = pcoaA2_16S,
       colour = "Site") +
  scale_colour_viridis_d() +
  theme_bw() +  
  theme(legend.position = "right",
        axis.title = element_text(face = "bold", size = 12), 
        axis.text = element_text(size = 10))
g4
pdf("InitialFigs/Beta_Jaccard_16S.pdf", width = 7, height = 5)
g4
dev.off()



#### __UniFrac ####
un_16S <- unifrac.distance
set.seed(1150)
ado9 <- adonis2(un_16S ~ input_filt_rare_16S$map_loaded$Site * input_filt_rare_16S$map_loaded$Entry)
ado9
# Site and interaction significant.
set.seed(1150)
ado10 <- adonis2(un_16S ~ input_filt_rare_16S$map_loaded$Entry * input_filt_rare_16S$map_loaded$Site)
ado10
# Site and interaction significant. No effect of order.
anova(betadisper(un_16S, input_filt_rare_16S$map_loaded$Site)) # Dispersion not homogeneous
anova(betadisper(un_16S, input_filt_rare_16S$map_loaded$Entry)) # Dispersion homogeneous

# PCoA with vectors
pcoa_16S <- cmdscale(un_16S, k = nrow(input_filt_rare_16S$map_loaded) - 1, eig = T)
set.seed(100)
ef_16S <- envfit(pcoa_16S, env_16S, permutations = 999, na.rm = TRUE)
ef_16S # All sig but plant variables (except chlorophyll)
ordiplot(pcoa_16S)
plot(ef_16S, p.max = 0.05, cex = 0.5)
multiplier_16S <- ordiArrowMul(ef_16S)
vec.df_16S <- as.data.frame(ef_16S$vectors$arrows*sqrt(ef_16S$vectors$r)) %>%
  mutate(Dim1 = Dim1 * multiplier_16S,
         Dim2 = Dim2 * multiplier_16S) %>%
  mutate(variables = rownames(.)) %>%
  filter(ef_16S$vectors$pvals < 0.05) %>%
  mutate(shortnames = c("LAT", "LONG", "CMI", "GDD", "Temp", "Precip", "Tseas", "Pseas",
                        "pH", "OM", "Salts", "NO3", "PO4", "K", "Ca", "Mg", "Na", "S",
                        "CEC", "Sand", "Silt", "Clay", "Chlorophyll"))
pcoaA1_16S <- paste("PC1: ", round((eigenvals(pcoa_16S)/sum(eigenvals(pcoa_16S)))[1]*100, 1), "%")
pcoaA2_16S <- paste("PC2: ", round((eigenvals(pcoa_16S)/sum(eigenvals(pcoa_16S)))[2]*100, 1), "%")
input_filt_rare_16S$map_loaded$Axis01 <- vegan::scores(pcoa_16S)[,1]
input_filt_rare_16S$map_loaded$Axis02 <- vegan::scores(pcoa_16S)[,2]
micro.hulls_16S <- ddply(input_filt_rare_16S$map_loaded, c("Site"), find_hull)
g4 <- ggplot(input_filt_rare_16S$map_loaded, aes(Axis01, Axis02)) +
  geom_polygon(data = micro.hulls_16S, 
               aes(colour = Site),
               alpha = 0, show.legend = F) +
  geom_point(size = 3, pch = 16, alpha = 0.75, aes(colour = Site)) +
  geom_segment(data = vec.df_16S,
               aes(x = 0, xend = Dim1, y = 0, yend = Dim2),
               arrow = arrow(length = unit(0.35, "cm")),
               colour = "gray", alpha = 0.6,
               inherit.aes = FALSE) + 
  geom_text_repel(data = vec.df_16S,
                  aes(x = Dim1, y = Dim2, label = shortnames),
                  size = 3, color = "black") +
  labs(x = pcoaA1_16S, 
       y = pcoaA2_16S,
       colour = "Site") +
  scale_colour_viridis_d() +
  theme_bw() +  
  theme(legend.position = "right",
        axis.title = element_text(face = "bold", size = 12), 
        axis.text = element_text(size = 10))
g4
pdf("InitialFigs/Beta_UniFrac_16S.pdf", width = 7, height = 5)
g4
dev.off()



#### __Compare ####
# Compare the 4 distance matrices
set.seed(402)
mantel(bc_16S, aclr_16S) # 0.50, 0.001
set.seed(402)
mantel(bc_16S, unifrac.distance) # 0.77, 0.001
set.seed(402)
mantel(bc_16S, ja_16S) # 0.9965, 0.001



#### _ITS ####

#### __Aitchison ####
# CLR transformation. Can do on PC since not as many ASVs as 16S
# Ran with default 80% filter
asv_czm_ITS <- cmultRepl(t(input_filt_ITS$data_loaded), label = 0, method = "CZM")
asv_clr_ITS <- clr(asv_czm_ITS)
aclr_ITS <- compositions::dist(asv_clr_ITS)

# Filter dropped samples
# Columns and rows containing more than 80% zeros/unobserved values were deleted
dim(t(input_filt_ITS$data_loaded)) # 587 samples, 5321 ASVs
dim(asv_czm_ITS) # 569 samples, 138 ASVs. Note: Only highly prevalent ASVs retained!
input_filt_clr_ITS <- filter_data(input_filt_ITS,
                                  filter_cat = "BoxCell",
                                  keep_vals = rownames(asv_czm_ITS))

set.seed(1150)
ado11 <- adonis2(aclr_ITS ~ input_filt_clr_ITS$map_loaded$Site * input_filt_clr_ITS$map_loaded$Entry)
ado11
# Site and interaction significant. 0.64, 0.001***; 0.01, 0.098; 0.08, 0.013*
set.seed(1150)
ado12 <- adonis2(aclr_ITS ~ input_filt_clr_ITS$map_loaded$Entry * input_filt_clr_ITS$map_loaded$Site)
# Both and interaction significant. Small effect of order. 0.64, 0.001***; 0.01, 0.021*; 0.08, 0.013*
anova(betadisper(aclr_ITS, input_filt_clr_ITS$map_loaded$Site)) # Dispersion not homogeneous
anova(betadisper(aclr_ITS, input_filt_clr_ITS$map_loaded$Entry)) # Dispersion homogeneous

# PCA with vectors
d.pcx_ITS <- prcomp(aclr_ITS)
env_clr <- input_filt_clr_ITS$map_loaded %>%
  dplyr::select(LAT, LONG, # Geography
                CMI, GDD, Temp, P, Tseas, Pseas, # Climate
                pH, OM_perc, Salts_mmhos_cm, NO3_ppm, P_ppm, K_ppm, Ca_ppm, Mg_ppm, # Soil
                Na_ppm, S_ppm, CEC_meq_g, Sand_perc, Silt_perc, Clay_perc,
                chlorophyll, Sclerotinia, GRIN_Sclerotinia, DaysToFlower, StemLengthHighRange, # Plant
                StemLengthLowRange, PerSeed_14_64, PerSeed_20_64, PerNoBranch, SeedWeight)
set.seed(100)
ef_clr_ITS <- envfit(d.pcx_ITS, env_clr, permutations = 999, na.rm = TRUE)
ef_clr_ITS # All significant but plant (except chlorophyll)
ordiplot(d.pcx_ITS)
plot(ef_clr_ITS, p.max = 0.05, cex = 0.5)
multiplier_clr_ITS <- ordiArrowMul(ef_clr_ITS)
vec.df_clr_ITS <- as.data.frame(ef_clr_ITS$vectors$arrows*sqrt(ef_clr_ITS$vectors$r)) %>%
  mutate(PC1 = PC1 * multiplier_clr_ITS,
         PC2 = PC2 * multiplier_clr_ITS) %>%
  mutate(variables = rownames(.)) %>%
  filter(ef_clr_ITS$vectors$pvals < 0.05) %>%
  mutate(shortnames = c("LAT", "LONG", "CMI", "GDD", "Temp", "Precip", "Tseas", "Pseas",
                        "pH", "OM", "Salts", "NO3", "PO4", "K", "Ca", "Mg", "Na", "S",
                        "CEC", "Sand", "Silt", "Clay", "Chlorophyll"))
pcoaA1_clr_ITS <- paste("PC1: ", round((eigenvals(d.pcx_ITS)/sum(eigenvals(d.pcx_ITS)))[1]*100, 1), "%")
pcoaA2_clr_ITS <- paste("PC2: ", round((eigenvals(d.pcx_ITS)/sum(eigenvals(d.pcx_ITS)))[2]*100, 1), "%")
input_filt_clr_ITS$map_loaded$Axis01 <- vegan::scores(d.pcx_ITS)[,1]
input_filt_clr_ITS$map_loaded$Axis02 <- vegan::scores(d.pcx_ITS)[,2]
# input_filt_clr_ITS$map_loaded$Axis01 <- input_filt_clr_ITS$map_loaded$Axis01/sqrt(sum((input_filt_clr_ITS$map_loaded$Axis01 - mean(input_filt_clr_ITS$map_loaded$Axis01))^2))
# input_filt_clr_ITS$map_loaded$Axis02 <- input_filt_clr_ITS$map_loaded$Axis02/sqrt(sum((input_filt_clr_ITS$map_loaded$Axis02 - mean(input_filt_clr_ITS$map_loaded$Axis02))^2))
micro.hulls_clr_ITS <- ddply(input_filt_clr_ITS$map_loaded, c("Site"), find_hull)
g6 <- ggplot(input_filt_clr_ITS$map_loaded, aes(Axis01, Axis02)) +
  geom_polygon(data = micro.hulls_clr_ITS, 
               aes(colour = Site),
               alpha = 0, show.legend = F) +
  geom_point(size = 3, pch = 16, alpha = 0.75, aes(colour = Site)) +
  geom_segment(data = vec.df_clr_ITS,
               aes(x = 0, xend = PC1, y = 0, yend = PC2),
               arrow = arrow(length = unit(0.35, "cm")),
               colour = "gray", alpha = 0.6,
               inherit.aes = FALSE) +
  geom_text_repel(data = vec.df_clr_ITS,
                  aes(x = PC1, y = PC2, label = shortnames),
                  size = 3, color = "black") +
  labs(x = pcoaA1_clr_ITS, 
       y = pcoaA2_clr_ITS,
       colour = "Site") +
  scale_colour_viridis_d() +
  theme_bw() +  
  theme(legend.position = "right",
        axis.title = element_text(face = "bold", size = 12), 
        axis.text = element_text(size = 10))
g6
pdf("InitialFigs/Beta_Aitch_ITS.pdf", width = 7, height = 5)
g6
dev.off()

# Add inset
eta_sq_m4 <- eta_sq_adonis(ado11) %>%
  as.data.frame() %>%
  rownames_to_column(var = "group") %>%
  filter(group != "Total") %>%
  set_names(c("group", "value")) %>%
  mutate(group = c("Site***", "Gen.", "Int.*", "Resid.")) %>%
  arrange(desc(group)) %>%
  mutate(ypos = cumsum(value) - 0.5*value)
g7 <- ggplot(eta_sq_m4, aes(x = "", y = value, fill = group)) +
  geom_bar(stat = "identity", width = 1, color = NA) +
  coord_polar("y", start = 0) +
  geom_text_repel(data = eta_sq_m4,
                  aes(y = ypos, label = group), color = "white", size = 3,
                  box.padding = 0.1) +
  scale_fill_manual(values = c("#5DC863FF", "grey30", "grey70", "#21908CFF")) +
  theme_void() + 
  theme(legend.position = "none")
g7

plot.with.inset4 <-
  ggdraw() +
  draw_plot(g6) +
  draw_plot(g7, x = 0.54, y = 0.08, width = 0.25, height = 0.25)
plot.with.inset4
pdf("InitialFigs/Beta_Aitch_wInset_ITS.pdf", width = 7, height = 5)
plot.with.inset4
dev.off()



# Rerun with all samples and taxa
set.seed(1150)
asv_czm_ITS <- cmultRepl(t(input_filt_ITS$data_loaded), label = 0, method = "CZM",
                         z.warning = 1)
asv_clr_ITS <- clr(asv_czm_ITS)
aclr_ITS <- compositions::dist(asv_clr_ITS)

set.seed(1150)
ado13 <- adonis2(aclr_ITS ~ input_filt_ITS$map_loaded$Site * input_filt_ITS$map_loaded$Entry)
ado13 # Site and interaction significant.
set.seed(1150)
ado14 <- adonis2(aclr_ITS ~ input_filt_ITS$map_loaded$Entry * input_filt_ITS$map_loaded$Site)
ado14 # Both and interaction significant. Small effect of order.
anova(betadisper(aclr_ITS, input_filt_ITS$map_loaded$Site)) # Dispersion not homogeneous
anova(betadisper(aclr_ITS, input_filt_ITS$map_loaded$Entry)) # Dispersion homogeneous

# PCA with vectors
d.pcx_ITS <- prcomp(aclr_ITS)
env_clr <- input_filt_ITS$map_loaded %>%
  dplyr::select(LAT, LONG, # Geography
                CMI, GDD, Temp, P, Tseas, Pseas, # Climate
                pH, OM_perc, Salts_mmhos_cm, NO3_ppm, P_ppm, K_ppm, Ca_ppm, Mg_ppm, # Soil
                Na_ppm, S_ppm, CEC_meq_g, Sand_perc, Silt_perc, Clay_perc,
                chlorophyll, Sclerotinia, GRIN_Sclerotinia, DaysToFlower, StemLengthHighRange, # Plant
                StemLengthLowRange, PerSeed_14_64, PerSeed_20_64, PerNoBranch, SeedWeight)
set.seed(100)
ef_clr_ITS <- envfit(d.pcx_ITS, env_clr, permutations = 999, na.rm = TRUE)
ef_clr_ITS # Same
ordiplot(d.pcx_ITS)
plot(ef_clr_ITS, p.max = 0.05, cex = 0.5)
multiplier_ITS <- ordiArrowMul(ef_clr_ITS)
vec.df_ITS <- as.data.frame(ef_clr_ITS$vectors$arrows*sqrt(ef_clr_ITS$vectors$r)) %>%
  mutate(PC1 = PC1 * multiplier_ITS,
         PC2 = PC2 * multiplier_ITS) %>%
  mutate(variables = rownames(.)) %>%
  filter(ef_clr_ITS$vectors$pvals < 0.05) %>%
  mutate(shortnames = c("LAT", "LONG", "CMI", "GDD", "Temp", "Precip", "Tseas", "Pseas",
                        "pH", "OM", "Salts", "NO3", "PO4", "K", "Ca", "Mg", "Na", "S",
                        "CEC", "Sand", "Silt", "Clay", "Chlorophyll"))
pcoaA1_ITS <- paste("PC1: ", round((eigenvals(d.pcx_ITS)/sum(eigenvals(d.pcx_ITS)))[1]*100, 1), "%")
pcoaA2_ITS <- paste("PC2: ", round((eigenvals(d.pcx_ITS)/sum(eigenvals(d.pcx_ITS)))[2]*100, 1), "%")
input_filt_ITS$map_loaded$Axis01 <- vegan::scores(d.pcx_ITS)[,1]
input_filt_ITS$map_loaded$Axis02 <- vegan::scores(d.pcx_ITS)[,2]
# input_filt_ITS$map_loaded$Axis01 <- input_filt_ITS$map_loaded$Axis01/sqrt(sum((input_filt_ITS$map_loaded$Axis01 - mean(input_filt_ITS$map_loaded$Axis01))^2))
# input_filt_ITS$map_loaded$Axis02 <- input_filt_ITS$map_loaded$Axis02/sqrt(sum((input_filt_ITS$map_loaded$Axis02 - mean(input_filt_ITS$map_loaded$Axis02))^2))
micro.hulls_ITS <- ddply(input_filt_ITS$map_loaded, c("Site"), find_hull)
g1 <- ggplot(input_filt_ITS$map_loaded, aes(Axis01, Axis02)) +
  geom_polygon(data = micro.hulls_ITS, 
               aes(colour = Site),
               alpha = 0, show.legend = F) +
  geom_point(size = 3, pch = 16, alpha = 0.75, aes(colour = Site)) +
  geom_segment(data = vec.df_ITS,
               aes(x = 0, xend = PC1, y = 0, yend = PC2),
               arrow = arrow(length = unit(0.35, "cm")),
               colour = "gray", alpha = 0.6,
               inherit.aes = FALSE) +
  geom_text_repel(data = vec.df_ITS,
                  aes(x = PC1, y = PC2, label = shortnames),
                  size = 3, color = "black") +
  labs(x = pcoaA1_ITS, 
       y = pcoaA2_ITS,
       colour = "Site") +
  scale_colour_viridis_d() +
  theme_bw() +  
  theme(legend.position = "right",
        axis.title = element_text(face = "bold", size = 12), 
        axis.text = element_text(size = 10))
g1

pdf("InitialFigs/Beta_Aitch_All_ITS.pdf", width = 7, height = 5)
g1
dev.off()



#### __Bray ####
bc_ITS <- calc_dm(input_filt_rare_ITS$data_loaded)
set.seed(1150)
ado15 <- adonis2(bc_ITS ~ input_filt_rare_ITS$map_loaded$Site * input_filt_rare_ITS$map_loaded$Entry)
ado15
# Site and interaction significant. 0.62, 0.001; 0.01, 0.072; 0.09, 0.001
set.seed(1150)
ado16 <- adonis2(bc_ITS ~ input_filt_rare_ITS$map_loaded$Entry * input_filt_rare_ITS$map_loaded$Site)
ado16
# Both and interaction significant. Small effect of order. 0.62, 0.001; 0.01, 0.024; 0.09, 0.001
anova(betadisper(bc_ITS, input_filt_rare_ITS$map_loaded$Site)) # Dispersion not homogeneous
anova(betadisper(bc_ITS, input_filt_rare_ITS$map_loaded$Entry)) # Dispersion homogeneous

# Test effect of genotype while stratifying permutations by Site
set.seed(1150)
adonis2(bc_ITS ~ input_filt_rare_ITS$map_loaded$Entry, 
        strata = input_filt_rare_ITS$map_loaded$Site)

# Pairwise by Site
set.seed(1150)
pwa_ITS <- as.data.frame(pairwiseAdonis::pairwise.adonis(bc_ITS, 
                                                         input_filt_rare_ITS$map_loaded$Site))

# PCoA with vectors
pcoa_ITS <- cmdscale(bc_ITS, k = nrow(input_filt_rare_ITS$map_loaded) - 1, eig = T)
env_ITS <- input_filt_rare_ITS$map_loaded %>%
  dplyr::select(LAT, LONG, # Geography
                CMI, GDD, Temp, P, Tseas, Pseas, # Climate
                pH, OM_perc, Salts_mmhos_cm, NO3_ppm, P_ppm, K_ppm, Ca_ppm, Mg_ppm, # Soil
                Na_ppm, S_ppm, CEC_meq_g, Sand_perc, Silt_perc, Clay_perc,
                chlorophyll, Sclerotinia, GRIN_Sclerotinia, DaysToFlower, StemLengthHighRange, # Plant
                StemLengthLowRange, PerSeed_14_64, PerSeed_20_64, PerNoBranch, SeedWeight)
set.seed(100)
ef_ITS <- envfit(pcoa_ITS, env_ITS, permutations = 999, na.rm = TRUE)
ef_ITS # All significant but plant (except chlorophyll)
ordiplot(pcoa_ITS)
plot(ef_ITS, p.max = 0.05, cex = 0.5)
multiplier_ITS <- ordiArrowMul(ef_ITS)
vec.df_ITS <- as.data.frame(ef_ITS$vectors$arrows*sqrt(ef_ITS$vectors$r)) %>%
  mutate(Dim1 = Dim1 * multiplier_ITS,
         Dim2 = Dim2 * multiplier_ITS) %>%
  mutate(variables = rownames(.)) %>%
  filter(ef_ITS$vectors$pvals < 0.05) %>%
  mutate(shortnames = c("LAT", "LONG", "CMI", "GDD", "Temp", "Precip", "Tseas", "Pseas",
                        "pH", "OM", "Salts", "NO3", "PO4", "K", "Ca", "Mg", "Na", "S",
                        "CEC", "Sand", "Silt", "Clay", "Chlorophyll"))
pcoaA1_ITS <- paste("PC1: ", round((eigenvals(pcoa_ITS)/sum(eigenvals(pcoa_ITS)))[1]*100, 1), "%")
pcoaA2_ITS <- paste("PC2: ", round((eigenvals(pcoa_ITS)/sum(eigenvals(pcoa_ITS)))[2]*100, 1), "%")
input_filt_rare_ITS$map_loaded$Axis01 <- vegan::scores(pcoa_ITS)[,1]
input_filt_rare_ITS$map_loaded$Axis02 <- vegan::scores(pcoa_ITS)[,2]
micro.hulls_ITS <- ddply(input_filt_rare_ITS$map_loaded, c("Site"), find_hull)
g11 <- ggplot(input_filt_rare_ITS$map_loaded, aes(Axis01, Axis02)) +
  geom_polygon(data = micro.hulls_ITS, 
               aes(colour = Site),
               alpha = 0, show.legend = F) +
  geom_point(size = 3, pch = 16, alpha = 0.75, aes(colour = Site)) +
  geom_segment(data = vec.df_ITS,
               aes(x = 0, xend = Dim1, y = 0, yend = Dim2),
               arrow = arrow(length = unit(0.35, "cm")),
               colour = "gray", alpha = 0.6,
               inherit.aes = FALSE) + 
  geom_text_repel(data = vec.df_ITS,
                  aes(x = Dim1, y = Dim2, label = shortnames),
                  size = 3, color = "black") +
  labs(x = pcoaA1_ITS, 
       y = pcoaA2_ITS,
       colour = "Site") +
  scale_colour_viridis_d() +
  theme_bw() +  
  theme(legend.position = "right",
        axis.title = element_text(face = "bold", size = 12), 
        axis.text = element_text(size = 10))
g11
pdf("InitialFigs/Beta_Bray_ITS.pdf", width = 7, height = 5)
g11
dev.off()

leg <- get_legend((g11))

# Add inset
eta_sq_m5 <- eta_sq_adonis(ado15) %>%
  as.data.frame() %>%
  rownames_to_column(var = "group") %>%
  filter(group != "Total") %>%
  set_names(c("group", "value")) %>%
  mutate(group = c("S***", "G*", "I***", "R")) %>%
  arrange(desc(group)) %>%
  mutate(ypos = cumsum(value) - 0.5*value)
g12 <- ggplot(eta_sq_m5, aes(x = "", y = value, fill = group)) +
  geom_bar(stat = "identity", width = 1, color = NA) +
  coord_polar("y", start = 0) +
  geom_text_repel(data = eta_sq_m5,
                  aes(y = ypos, label = group), color = "white", size = 3,
                  box.padding = 0.1) +
  scale_fill_manual(values = c("#5DC863FF", "grey30", "grey70", "#21908CFF")) +
  theme_void() + 
  theme(legend.position = "none")
g12

plot.with.inset5 <-
  ggdraw() +
  draw_plot(g8) +
  draw_plot(g9, x = 0.54, y = 0.75, width = 0.25, height = 0.25)
plot.with.inset5
pdf("InitialFigs/Beta_Bray_wInset_ITS.pdf", width = 7, height = 5)
plot.with.inset5
dev.off()

# Make ordinations with site centroids and labels
# Also flip axes to better match map
pd <- betadisper(bc_ITS, input_filt_rare_ITS$map_loaded$Site)
centroids <- as.data.frame(pd$centroids) %>%
  dplyr::select(PCoA1, PCoA2) %>%
  rownames_to_column(var = "Site")
pcoa_ITS <- cmdscale(bc_ITS, k = nrow(input_filt_rare_ITS$map_loaded) - 1, eig = T)
pcoaA1 <- paste("PC1: ", round((eigenvals(pcoa_ITS)/sum(eigenvals(pcoa_ITS)))[1]*100, 1), "%")
pcoaA2 <- paste("PC2: ", round((eigenvals(pcoa_ITS)/sum(eigenvals(pcoa_ITS)))[2]*100, 1), "%")
input_filt_rare_ITS$map_loaded$Axis01 <- vegan::scores(pcoa_ITS)[,1]
input_filt_rare_ITS$map_loaded$Axis02 <- vegan::scores(pcoa_ITS)[,2]
micro.hulls <- ddply(input_filt_rare_ITS$map_loaded, c("Site"), find_hull)
pdf("InitialFigs/Beta_Bray_Centroids_ITS.pdf", width = 7, height = 5)
ggplot(input_filt_rare_ITS$map_loaded, aes(-Axis01, -Axis02)) +
  stat_ellipse(mapping = aes(x = -Axis01, y = -Axis02, color = Site, fill = Site), 
               geom = "polygon", alpha = 0.2, linewidth = NA) +
  geom_point(size = 1.5, pch = 16, alpha = 0.375, aes(colour = Site)) +
  geom_point(data = centroids, aes(-PCoA1, -PCoA2, fill = Site),
             size = 5, shape = 24, stroke = 1) +
  geom_text_repel(data = centroids, aes(-PCoA1, -PCoA2, label = Site),
                  box.padding = 0.7) +
  labs(x = pcoaA1, 
       y = pcoaA2,
       colour = "Site") +
  scale_colour_viridis_d() +
  scale_fill_viridis_d() +
  theme_bw() +  
  theme(legend.position = "none",
        axis.title = element_text(face = "bold", size = 12), 
        axis.text = element_text(size = 10),
        plot.title = element_text(vjust = -1))
dev.off()

# Add vectors
pdf("InitialFigs/Beta_Bray_Centroids_wVec_ITS.pdf", width = 7, height = 5)
g13 <- ggplot(input_filt_rare_ITS$map_loaded, aes(-Axis01, -Axis02)) +
  stat_ellipse(mapping = aes(x = -Axis01, y = -Axis02, color = Site, fill = Site), 
               geom = "polygon", alpha = 0.2, linewidth = NA) +
  geom_point(size = 1.5, pch = 16, alpha = 0.375, aes(colour = Site)) +
  geom_point(data = centroids, aes(-PCoA1, -PCoA2, fill = Site),
             size = 5, shape = 24, stroke = 1) +
  geom_segment(data = vec.df_ITS,
               aes(x = 0, xend = -Dim1, y = 0, yend = -Dim2),
               arrow = arrow(length = unit(0.35, "cm")),
               colour = "gray", alpha = 0.25,
               inherit.aes = FALSE) + 
  geom_text_repel(data = vec.df_ITS,
                  aes(x = -Dim1, y = -Dim2, label = shortnames),
                  size = 2, color = "red", box.padding = 0.1) +
  # geom_text(data = subset(vec.df_ITS, shortnames %in% c("CMI")),
  #           aes(x = -Dim1, y = -Dim2, label = shortnames),
  #           size = 2, color = "red", nudge_y = -0.01, nudge_x = 0.03) +
  # geom_text(data = subset(vec.df_ITS, shortnames %in% c("Mg")),
  #           aes(x = -Dim1, y = -Dim2, label = shortnames),
  #           size = 2, color = "red", nudge_y = 0.005) +
  # geom_text(data = subset(vec.df_ITS, shortnames %in% c("OM")),
  #           aes(x = -Dim1, y = -Dim2, label = shortnames),
  #           size = 2, color = "red", nudge_y = -0.01, nudge_x = 0.01) +
  # geom_text(data = subset(vec.df_ITS, shortnames %in% c("P04")),
  #           aes(x = -Dim1, y = -Dim2, label = shortnames),
  #           size = 2, color = "red", nudge_y = -0.02) +
  geom_text_repel(data = centroids, aes(-PCoA1, -PCoA2, label = Site),
                  box.padding = 0.7, size = 3) +
  labs(x = pcoaA1_ITS, 
       y = pcoaA2_ITS,
       colour = "Site") +
  scale_colour_viridis_d() +
  scale_fill_viridis_d() +
  theme_bw() +  
  theme(legend.position = "none",
        axis.title = element_text(face = "bold", size = 12), 
        axis.text = element_text(size = 10),
        plot.title = element_text(vjust = -1),
        panel.grid = element_blank())
g13
dev.off()

# Check ordination at all taxonomic levels. Need to see if it's consistent with ASV level.
tax_sum_phyla <- summarize_taxonomy(input = input_filt_rare_ITS, 
                                    level = 2, 
                                    report_higher_tax = T)
phy_bc <- calc_dm(tax_sum_phyla)
pcoa_ITS <- cmdscale(phy_bc, k = nrow(input_filt_rare_ITS$map_loaded) - 1, eig = T)
pcoaA1 <- paste("PC1: ", round((eigenvals(pcoa_ITS)/sum(eigenvals(pcoa_ITS)))[1]*100, 1), "%")
pcoaA2 <- paste("PC2: ", round((eigenvals(pcoa_ITS)/sum(eigenvals(pcoa_ITS)))[2]*100, 1), "%")
input_filt_rare_ITS$map_loaded$Axis01 <- vegan::scores(pcoa_ITS)[,1]
input_filt_rare_ITS$map_loaded$Axis02 <- vegan::scores(pcoa_ITS)[,2]
micro.hulls <- ddply(input_filt_rare_ITS$map_loaded, c("Site"), find_hull)
g_p <- ggplot(input_filt_rare_ITS$map_loaded, aes(Axis01, Axis02)) +
  geom_polygon(data = micro.hulls, 
               aes(colour = Site),
               alpha = 0, show.legend = F) +
  geom_point(size = 2, pch = 16, alpha = 0.75, aes(colour = Site)) +
  labs(x = pcoaA1, 
       y = pcoaA2,
       colour = "Site") +
  scale_colour_viridis_d() +
  ggtitle("Phylum") +
  theme_bw() +  
  theme(legend.position = "none",
        axis.title = element_text(face = "bold", size = 12), 
        axis.text = element_text(size = 10),
        plot.title = element_text(vjust = -1))
g_p

tax_sum_classes <- summarize_taxonomy(input = input_filt_rare_ITS, 
                                      level = 3, 
                                      report_higher_tax = T)
cla_bc <- calc_dm(tax_sum_classes)
pcoa_ITS <- cmdscale(cla_bc, k = nrow(input_filt_rare_ITS$map_loaded) - 1, eig = T)
pcoaA1 <- paste("PC1: ", round((eigenvals(pcoa_ITS)/sum(eigenvals(pcoa_ITS)))[1]*100, 1), "%")
pcoaA2 <- paste("PC2: ", round((eigenvals(pcoa_ITS)/sum(eigenvals(pcoa_ITS)))[2]*100, 1), "%")
input_filt_rare_ITS$map_loaded$Axis01 <- vegan::scores(pcoa_ITS)[,1]
input_filt_rare_ITS$map_loaded$Axis02 <- vegan::scores(pcoa_ITS)[,2]
micro.hulls <- ddply(input_filt_rare_ITS$map_loaded, c("Site"), find_hull)
g_c <- ggplot(input_filt_rare_ITS$map_loaded, aes(Axis01, Axis02)) +
  geom_polygon(data = micro.hulls, 
               aes(colour = Site),
               alpha = 0, show.legend = F) +
  geom_point(size = 2, pch = 16, alpha = 0.75, aes(colour = Site)) +
  labs(x = pcoaA1, 
       y = pcoaA2,
       colour = "Site") +
  scale_colour_viridis_d() +
  ggtitle("Class") +
  theme_bw() +  
  theme(legend.position = "none",
        axis.title = element_text(face = "bold", size = 12), 
        axis.text = element_text(size = 10),
        plot.title = element_text(vjust = -1))
g_c

tax_sum_orders <- summarize_taxonomy(input = input_filt_rare_ITS, 
                                     level = 4, 
                                     report_higher_tax = T)
ord_bc <- calc_dm(tax_sum_orders)
pcoa_ITS <- cmdscale(ord_bc, k = nrow(input_filt_rare_ITS$map_loaded) - 1, eig = T)
pcoaA1 <- paste("PC1: ", round((eigenvals(pcoa_ITS)/sum(eigenvals(pcoa_ITS)))[1]*100, 1), "%")
pcoaA2 <- paste("PC2: ", round((eigenvals(pcoa_ITS)/sum(eigenvals(pcoa_ITS)))[2]*100, 1), "%")
input_filt_rare_ITS$map_loaded$Axis01 <- vegan::scores(pcoa_ITS)[,1]
input_filt_rare_ITS$map_loaded$Axis02 <- vegan::scores(pcoa_ITS)[,2]
micro.hulls <- ddply(input_filt_rare_ITS$map_loaded, c("Site"), find_hull)
g_o <- ggplot(input_filt_rare_ITS$map_loaded, aes(Axis01, -Axis02)) +
  geom_polygon(data = micro.hulls, 
               aes(colour = Site),
               alpha = 0, show.legend = F) +
  geom_point(size = 2, pch = 16, alpha = 0.75, aes(colour = Site)) +
  labs(x = pcoaA1, 
       y = pcoaA2,
       colour = "Site") +
  scale_colour_viridis_d() +
  ggtitle("Order") +
  theme_bw() +  
  theme(legend.position = "none",
        axis.title = element_text(face = "bold", size = 12), 
        axis.text = element_text(size = 10),
        plot.title = element_text(vjust = -1))
g_o

tax_sum_families <- summarize_taxonomy(input = input_filt_rare_ITS, 
                                       level = 5, 
                                       report_higher_tax = T)
fam_bc <- calc_dm(tax_sum_families)
pcoa_ITS <- cmdscale(fam_bc, k = nrow(input_filt_rare_ITS$map_loaded) - 1, eig = T)
pcoaA1 <- paste("PC1: ", round((eigenvals(pcoa_ITS)/sum(eigenvals(pcoa_ITS)))[1]*100, 1), "%")
pcoaA2 <- paste("PC2: ", round((eigenvals(pcoa_ITS)/sum(eigenvals(pcoa_ITS)))[2]*100, 1), "%")
input_filt_rare_ITS$map_loaded$Axis01 <- vegan::scores(pcoa_ITS)[,1]
input_filt_rare_ITS$map_loaded$Axis02 <- vegan::scores(pcoa_ITS)[,2]
micro.hulls <- ddply(input_filt_rare_ITS$map_loaded, c("Site"), find_hull)
g_f <- ggplot(input_filt_rare_ITS$map_loaded, aes(Axis01, -Axis02)) +
  geom_polygon(data = micro.hulls, 
               aes(colour = Site),
               alpha = 0, show.legend = F) +
  geom_point(size = 2, pch = 16, alpha = 0.75, aes(colour = Site)) +
  labs(x = pcoaA1, 
       y = pcoaA2,
       colour = "Site") +
  scale_colour_viridis_d() +
  ggtitle("Family") +
  theme_bw() +  
  theme(legend.position = "none",
        axis.title = element_text(face = "bold", size = 12), 
        axis.text = element_text(size = 10),
        plot.title = element_text(vjust = -1))
g_f

tax_sum_genera <- summarize_taxonomy(input = input_filt_rare_ITS, 
                                     level = 6, 
                                     report_higher_tax = T)
gen_bc <- calc_dm(tax_sum_genera)
pcoa_ITS <- cmdscale(gen_bc, k = nrow(input_filt_rare_ITS$map_loaded) - 1, eig = T)
pcoaA1 <- paste("PC1: ", round((eigenvals(pcoa_ITS)/sum(eigenvals(pcoa_ITS)))[1]*100, 1), "%")
pcoaA2 <- paste("PC2: ", round((eigenvals(pcoa_ITS)/sum(eigenvals(pcoa_ITS)))[2]*100, 1), "%")
input_filt_rare_ITS$map_loaded$Axis01 <- vegan::scores(pcoa_ITS)[,1]
input_filt_rare_ITS$map_loaded$Axis02 <- vegan::scores(pcoa_ITS)[,2]
micro.hulls <- ddply(input_filt_rare_ITS$map_loaded, c("Site"), find_hull)
g_g <- ggplot(input_filt_rare_ITS$map_loaded, aes(Axis01, Axis02)) +
  geom_polygon(data = micro.hulls, 
               aes(colour = Site),
               alpha = 0, show.legend = F) +
  geom_point(size = 2, pch = 16, alpha = 0.75, aes(colour = Site)) +
  labs(x = pcoaA1, 
       y = pcoaA2,
       colour = "Site") +
  scale_colour_viridis_d() +
  ggtitle("Genus") +
  theme_bw() +  
  theme(legend.position = "none",
        axis.title = element_text(face = "bold", size = 12), 
        axis.text = element_text(size = 10),
        plot.title = element_text(vjust = -1))
g_g

pcoa_ITS <- cmdscale(bc_ITS, k = nrow(input_filt_rare_ITS$map_loaded) - 1, eig = T)
pcoaA1 <- paste("PC1: ", round((eigenvals(pcoa_ITS)/sum(eigenvals(pcoa_ITS)))[1]*100, 1), "%")
pcoaA2 <- paste("PC2: ", round((eigenvals(pcoa_ITS)/sum(eigenvals(pcoa_ITS)))[2]*100, 1), "%")
input_filt_rare_ITS$map_loaded$Axis01 <- vegan::scores(pcoa_ITS)[,1]
input_filt_rare_ITS$map_loaded$Axis02 <- vegan::scores(pcoa_ITS)[,2]
micro.hulls <- ddply(input_filt_rare_ITS$map_loaded, c("Site"), find_hull)
g_asv <- ggplot(input_filt_rare_ITS$map_loaded, aes(Axis01, Axis02)) +
  geom_polygon(data = micro.hulls, 
               aes(colour = Site),
               alpha = 0, show.legend = F) +
  geom_point(size = 2, pch = 16, alpha = 0.75, aes(colour = Site)) +
  labs(x = pcoaA1, 
       y = pcoaA2,
       colour = "Site") +
  scale_colour_viridis_d() +
  ggtitle("ASV") +
  theme_bw() +  
  theme(legend.position = "none",
        axis.title = element_text(face = "bold", size = 12), 
        axis.text = element_text(size = 10),
        plot.title = element_text(vjust = -1))
g_asv

ggplotly(g_asv) # Ralls and Kirkmeyer dont't cluster like 16S.
# Kirkmeyer and Lindsborg-dryland cluster.
mp1 <- plot_grid(g_p, g_c, g_o, g_f, g_g, g_asv, ncol = 3)
plot_grid(mp1, leg, rel_widths = c(0.85, 0.15))

pdf("InitialFigs/Beta_Bray_AllLevels_ITS.pdf", width = 9, height = 6)
plot_grid(mp1, leg, rel_widths = c(0.83, 0.17))
dev.off()



#### __Jaccard ####
ja_ITS <- calc_dm(input_filt_rare_ITS$data_loaded, method = "jaccard")
set.seed(1150)
ado17 <- adonis2(ja_ITS ~ input_filt_rare_ITS$map_loaded$Site * input_filt_rare_ITS$map_loaded$Entry)
ado17 # Site and interaction significant.
set.seed(1150)
ado18 <- adonis2(ja_ITS ~ input_filt_rare_ITS$map_loaded$Entry * input_filt_rare_ITS$map_loaded$Site)
ado18 # Site and interaction significant. Small effect of order (genotype significant 0.033*)
anova(betadisper(ja_ITS, input_filt_rare_ITS$map_loaded$Site)) # Dispersion not homogeneous
anova(betadisper(ja_ITS, input_filt_rare_ITS$map_loaded$Entry)) # Dispersion homogeneous

# PCoA with vectors
pcoa_ITS <- cmdscale(ja_ITS, k = nrow(input_filt_rare_ITS$map_loaded) - 1, eig = T)
set.seed(100)
ef_ITS <- envfit(pcoa_ITS, env_ITS, permutations = 999, na.rm = TRUE)
ef_ITS # All sig but plant variables (except chlorophyll)
ordiplot(pcoa_ITS)
plot(ef_ITS, p.max = 0.05, cex = 0.5)
multiplier_ITS <- ordiArrowMul(ef_ITS)
vec.df_ITS <- as.data.frame(ef_ITS$vectors$arrows*sqrt(ef_ITS$vectors$r)) %>%
  mutate(Dim1 = Dim1 * multiplier_ITS,
         Dim2 = Dim2 * multiplier_ITS) %>%
  mutate(variables = rownames(.)) %>%
  filter(ef_ITS$vectors$pvals < 0.05) %>%
  mutate(shortnames = c("LAT", "LONG", "CMI", "GDD", "Temp", "Precip", "Tseas", "Pseas",
                        "pH", "OM", "Salts", "NO3", "PO4", "K", "Ca", "Mg", "Na", "S",
                        "CEC", "Sand", "Silt", "Clay", "Chlorophyll"))
pcoaA1_ITS <- paste("PC1: ", round((eigenvals(pcoa_ITS)/sum(eigenvals(pcoa_ITS)))[1]*100, 1), "%")
pcoaA2_ITS <- paste("PC2: ", round((eigenvals(pcoa_ITS)/sum(eigenvals(pcoa_ITS)))[2]*100, 1), "%")
input_filt_rare_ITS$map_loaded$Axis01 <- vegan::scores(pcoa_ITS)[,1]
input_filt_rare_ITS$map_loaded$Axis02 <- vegan::scores(pcoa_ITS)[,2]
micro.hulls_ITS <- ddply(input_filt_rare_ITS$map_loaded, c("Site"), find_hull)
g4 <- ggplot(input_filt_rare_ITS$map_loaded, aes(Axis01, Axis02)) +
  geom_polygon(data = micro.hulls_ITS, 
               aes(colour = Site),
               alpha = 0, show.legend = F) +
  geom_point(size = 3, pch = 16, alpha = 0.75, aes(colour = Site)) +
  geom_segment(data = vec.df_ITS,
               aes(x = 0, xend = Dim1, y = 0, yend = Dim2),
               arrow = arrow(length = unit(0.35, "cm")),
               colour = "gray", alpha = 0.6,
               inherit.aes = FALSE) + 
  geom_text_repel(data = vec.df_ITS,
                  aes(x = Dim1, y = Dim2, label = shortnames),
                  size = 3, color = "black") +
  labs(x = pcoaA1_ITS, 
       y = pcoaA2_ITS,
       colour = "Site") +
  scale_colour_viridis_d() +
  theme_bw() +  
  theme(legend.position = "right",
        axis.title = element_text(face = "bold", size = 12), 
        axis.text = element_text(size = 10))
g4
pdf("InitialFigs/Beta_Jaccard_ITS.pdf", width = 7, height = 5)
g4
dev.off()



#### __Compare ####
# Compare the 4 distance matrices
# Need to drop the 1 sample from the CLR data
input_filt_drop <- filter_data(input_filt_rare_ITS,
                               filter_cat = "sampleID",
                               keep_vals = input_filt_rare_ITS$map_loaded$sampleID)
asv_czm_ITS_drop <- cmultRepl(t(input_filt_drop$data_loaded), label = 0, method = "CZM",
                              z.warning = 1)
asv_clr_ITS_drop <- clr(asv_czm_ITS_drop)
aclr_ITS_drop <- compositions::dist(asv_clr_ITS_drop)
set.seed(402)
mantel(bc_ITS, aclr_ITS_drop) # 0.89, 0.001
set.seed(402)
mantel(bc_ITS, ja_ITS) # 0.95, 0.001



#### _Combine ####
# Make multipanel figure 2
set.seed(1)
g10 <- ggplot(input_filt_rare_16S$map_loaded, aes(-Axis01, -Axis02)) +
  stat_ellipse(mapping = aes(x = -Axis01, y = -Axis02, color = Site, fill = Site), 
               geom = "polygon", alpha = 0.2, linewidth = NA) +
  geom_point(size = 1.5, pch = 16, alpha = 0.375, aes(colour = Site)) +
  geom_point(data = centroids_16S, aes(-PCoA1, -PCoA2, fill = Site),
             size = 5, shape = 24, stroke = 1) +
  geom_segment(data = vec.df_16S,
               aes(x = 0, xend = -Dim1, y = 0, yend = -Dim2),
               arrow = arrow(length = unit(0.35, "cm")),
               colour = "gray", alpha = 0.25,
               inherit.aes = FALSE) + 
  geom_text(data = subset(vec.df_16S,
                          shortnames %notin% c("Salts", "S")),
                  aes(x = -Dim1, y = -Dim2, label = shortnames),
                  size = 2, color = "red") +
  # geom_text_repel(data = subset(vec.df_16S,
  #                               shortnames %notin% c("LAT", "Precip", "P04")),
  #                 aes(x = -Dim1, y = -Dim2, label = shortnames),
  #                 size = 2, color = "red", box.padding = 0.05) +
  geom_text(data = subset(vec.df_16S, shortnames %in% c("Salts")),
            aes(x = -Dim1, y = -Dim2, label = shortnames),
            size = 2, color = "red", nudge_x = 0.03) +
  geom_text(data = subset(vec.df_16S, shortnames %in% c("S")),
            aes(x = -Dim1, y = -Dim2, label = shortnames),
            size = 2, color = "red", nudge_y = 0.01) +
  geom_text_repel(data = centroids_16S, aes(-PCoA1, -PCoA2, label = Site),
                  box.padding = 0.7, size = 3) +
  labs(x = pcoaA1, 
       y = pcoaA2,
       colour = "Site") +
  scale_colour_viridis_d() +
  scale_fill_viridis_d() +
  theme_bw() +  
  theme(legend.position = "none",
        axis.title = element_text(face = "bold", size = 12), 
        axis.text = element_text(size = 10),
        plot.title = element_text(vjust = -1),
        panel.grid = element_blank())
g10

set.seed(1)
g13 <- ggplot(input_filt_rare_ITS$map_loaded, aes(-Axis01, -Axis02)) +
  stat_ellipse(mapping = aes(x = -Axis01, y = -Axis02, color = Site, fill = Site), 
               geom = "polygon", alpha = 0.2, linewidth = NA) +
  geom_point(size = 1.5, pch = 16, alpha = 0.375, aes(colour = Site)) +
  geom_point(data = centroids, aes(-PCoA1, -PCoA2, fill = Site),
             size = 5, shape = 24, stroke = 1) +
  geom_segment(data = vec.df_ITS,
               aes(x = 0, xend = -Dim1, y = 0, yend = -Dim2),
               arrow = arrow(length = unit(0.35, "cm")),
               colour = "gray", alpha = 0.25,
               inherit.aes = FALSE) + 
  geom_text(data = subset(vec.df_ITS,
                          shortnames %notin% c("CMI", "OM", "Mg", "PO4", "GDD")),
            aes(x = -Dim1, y = -Dim2, label = shortnames),
            size = 2, color = "red") +
  # geom_text_repel(data = subset(vec.df_ITS,
  #                               shortnames %notin% c("Tseas", "CMI", "P04")),
  #                 aes(x = -Dim1, y = -Dim2, label = shortnames),
  #                 size = 2, color = "red", box.padding = 0.05) +
  geom_text(data = subset(vec.df_ITS, shortnames %in% c("CMI", "OM")),
            aes(x = -Dim1, y = -Dim2, label = shortnames),
            size = 2, color = "red", nudge_x = 0.03) +
  geom_text(data = subset(vec.df_ITS, shortnames %in% c("MG")),
            aes(x = -Dim1, y = -Dim2, label = shortnames),
            size = 2, color = "red", nudge_y = -0.01) +
  geom_text(data = subset(vec.df_ITS, shortnames %in% c("PO4")),
            aes(x = -Dim1, y = -Dim2, label = shortnames),
            size = 2, color = "red", nudge_y = -0.02) +
  geom_text(data = subset(vec.df_ITS, shortnames %in% c("GDD")),
            aes(x = -Dim1, y = -Dim2, label = shortnames),
            size = 2, color = "red", nudge_y = 0.02) +
  geom_text_repel(data = centroids, aes(-PCoA1, -PCoA2, label = Site),
                box.padding = 0.7, size = 3) +
  labs(x = pcoaA1_ITS, 
       y = pcoaA2_ITS,
       colour = "Site") +
  scale_colour_viridis_d() +
  scale_fill_viridis_d() +
  theme_bw() +  
  theme(legend.position = "none",
        axis.title = element_text(face = "bold", size = 12), 
        axis.text = element_text(size = 10),
        plot.title = element_text(vjust = -1),
        panel.grid = element_blank())
g13

g10 <- g10 + ggtitle("a) Archaea/Bacteria") +
  theme(plot.title = element_text(hjust = 0.5))
g13 <- g13 + ggtitle("b) Fungi") +
  theme(plot.title = element_text(hjust = 0.5))
fig2 <- plot_grid(g10, g13, g9, g12, ncol = 2, rel_heights = c(0.8, 0.2))
fig2

map <- plot_usmap(exclude = c("AK", "HI"),
                  color = "white",
                  fill = "grey80",
                  linewidth = 0.1) +
  geom_point(data = coords, 
             aes(x = x, y = y, fill = Site),
             color = "black",
             size = 2,
             shape = 24) +
  geom_text_repel(data = coords,
                  aes(x = x, y = y, label = Site),
                  size = 2, max.overlaps = 15, min.segment.length = 0) +
  scale_fill_viridis_d() +
  annotation_scale(location = "bl", height = unit(0.25, "cm")) +
  annotation_north_arrow(location = "br", height = unit(0.5, "cm"), width = unit(0.5, "cm")) +
  theme(legend.position = "none",
        plot.margin = margin(0,0,0,0, unit = "pt"))
map

upper <- plot_grid(g10, g13, ncol = 2)
lower <- plot_grid(g9, map, g12, ncol = 3, rel_widths = c(0.2, 0.6, 0.2))

pdf("InitialFigs/Figure2.pdf", width = 9, height = 6)
set.seed(1)
plot_grid(upper, lower, ncol = 1, rel_heights = c(0.7, 0.3))
dev.off()
# In Inkscape, make some final adjustments (map size, pie letters)



#### _RDA ####
# RDA is very slow with many variables. Run on server.
env_16S <- input_filt_rare_16S$map_loaded %>%
  dplyr::select(LAT, LONG, # Geography
                CMI, GDD, Temp, P, Tseas, Pseas, # Climate
                pH, OM_perc, Salts_mmhos_cm, NO3_ppm, P_ppm, K_ppm, Ca_ppm, Mg_ppm, # Soil
                Na_ppm, S_ppm, CEC_meq_g, Sand_perc, Silt_perc, Clay_perc,
                chlorophyll, Sclerotinia, GRIN_Sclerotinia, DaysToFlower, StemLengthHighRange, # Plant
                StemLengthLowRange, PerSeed_14_64, PerSeed_20_64, PerNoBranch, SeedWeight)
nrow(env_16S) # 587
sum(is.na(env_16S))
env_nona_16S <- env_16S %>%
  drop_na()
nrow(env_nona_16S) # 387
comm_nona_16S <- as.data.frame(t(input_filt_rare_16S$data_loaded)) %>%
  filter(rownames(.) %in% rownames(env_nona_16S))
mod0_16S <- rda(comm_nona_16S ~ 1, env_nona_16S)  # Model with intercept only
mod1_16S <- rda(comm_nona_16S ~ ., env_nona_16S)  # Model with all explanatory variables
mod_16S <- ordistep(mod0_16S, scope = formula(mod1_16S))
mod_16S
mod_16S$anova # LONG, CMI, Pseas, P, Tseas, T, LAT, GDD, Chlorophyll

env_ITS <- input_filt_rare_ITS$map_loaded %>%
  dplyr::select(LAT, LONG, # Geography
                CMI, GDD, Temp, P, Tseas, Pseas, # Climate
                pH, OM_perc, Salts_mmhos_cm, NO3_ppm, P_ppm, K_ppm, Ca_ppm, Mg_ppm, # Soil
                Na_ppm, S_ppm, CEC_meq_g, Sand_perc, Silt_perc, Clay_perc,
                chlorophyll, Sclerotinia, GRIN_Sclerotinia, DaysToFlower, StemLengthHighRange, # Plant
                StemLengthLowRange, PerSeed_14_64, PerSeed_20_64, PerNoBranch, SeedWeight)
nrow(env_ITS) # 586
sum(is.na(env_ITS))
env_nona_ITS <- env_ITS %>%
  drop_na()
nrow(env_nona_ITS) # 586
comm_nona_ITS <- as.data.frame(t(input_filt_rare_ITS$data_loaded)) %>%
  filter(rownames(.) %in% rownames(env_nona_ITS))
mod0_ITS <- rda(comm_nona_ITS ~ 1, env_nona_ITS)  # Model with intercept only
mod1_ITS <- rda(comm_nona_ITS ~ ., env_nona_ITS)  # Model with all explanatory variables
mod_ITS <- ordistep(mod0_ITS, scope = formula(mod1_ITS))
mod_ITS
mod_ITS$anova # CMI, GDD, Tseas, LONG, P, Pseas, Temp, LAT, chlorophyll



#### _Varpart ####
# varpart can take 2, 3, or 4 explanatory matrices
# Partition variation into geography, climate, soil, plant
mod_16S <- varpart(comm_nona_16S, 
                   ~ LAT + LONG, 
                   ~ CMI + GDD + Temp + P + Tseas + Pseas,
                   ~ pH + OM_perc + Salts_mmhos_cm + NO3_ppm + P_ppm + K_ppm + Ca_ppm + Mg_ppm +
                     Na_ppm + S_ppm + CEC_meq_g + Sand_perc + Silt_perc + Clay_perc, 
                   ~ chlorophyll + Sclerotinia + GRIN_Sclerotinia + DaysToFlower + StemLengthHighRange +
                     StemLengthLowRange + PerSeed_14_64 + PerSeed_20_64 + PerNoBranch + SeedWeight,
                   data = env_nona_16S)
mod_16S
summary(mod_16S)
pdf("InitialFigs/Varpart_16S.pdf", width = 7, height = 5)
plot(mod_16S, bg = 2:5, Xnames = c('Geog.', 'Clim.', 'Soil', 'Plant'))
dev.off()

mod_ITS <- varpart(comm_nona_ITS, 
                   ~ LAT + LONG, 
                   ~ CMI + GDD + Temp + P + Tseas + Pseas,
                   ~ pH + OM_perc + Salts_mmhos_cm + NO3_ppm + P_ppm + K_ppm + Ca_ppm + Mg_ppm +
                     Na_ppm + S_ppm + CEC_meq_g + Sand_perc + Silt_perc + Clay_perc, 
                   ~ chlorophyll + Sclerotinia + GRIN_Sclerotinia + DaysToFlower + StemLengthHighRange +
                     StemLengthLowRange + PerSeed_14_64 + PerSeed_20_64 + PerNoBranch + SeedWeight,
                   data = env_nona_ITS)
mod_ITS
summary(mod_ITS)
pdf("InitialFigs/Varpart_ITS.pdf", width = 7, height = 5)
plot(mod_ITS, bg = 2:5, Xnames = c('Geog.', 'Clim.', 'Soil', 'Plant'))
dev.off()

plot_grid(plot(mod_16S, bg = 2:5, Xnames = c('Geog.', 'Clim.', 'Soil', 'Plant')),
          plot(mod_ITS, bg = 2:5, Xnames = c('Geog.', 'Clim.', 'Soil', 'Plant')))



#### _Distance ####
# Make distance matrices

# Phylogenetic (16S only)
# Import phylogenetic tree. Use Phyloseq to calculate weighted UniFrac distance
otu <- otu_table(input_filt_rare_16S$data_loaded, taxa_are_rows = T)
tax <- tax_table(as.matrix(input_filt_rare_16S$taxonomy_loaded))
map <- sample_data(input_filt_rare_16S$map_loaded)
tree <- read_tree("data/rep_phylo_16S.tre")
is.rooted(tree) # FALSE. Tree is not rooted
tree <- midpoint.root(tree)
is.rooted(tree) # TRUE. Tree is now rooted at midpoint
input.phy <- phyloseq(otu, tax, map, tree)
unifrac.distance <- distance(input.phy, 
                             method = "wunifrac", 
                             type = "samples")
# Save (took a bit to calculate)
#saveRDS(unifrac.distance, "data/unifrac.distance.rds")
class(unifrac.distance)
dim(unifrac.distance)

# Geography
geography <- dplyr::select(input_filt_rare_16S$map_loaded, LONG, LAT)
geography.distance.mat <- rdist.earth(geography, miles = FALSE, R = 6371)
dim(geography.distance.mat)
#saveRDS(geography.distance.mat, "data/geography.distance.mat.rds")

geography_ITS <- dplyr::select(input_filt_rare_ITS$map_loaded, LONG, LAT)
geography.distance.mat_ITS <- rdist.earth(geography_ITS, miles = FALSE, R = 6371)
dim(geography.distance.mat_ITS)
#saveRDS(geography.distance.mat_ITS, "data/geography.distance.mat_ITS.rds")

# Climate
climate <- dplyr::select(input_filt_rare_16S$map_loaded, 
                         CMI, GDD, Temp, P, Tseas, Pseas)
climate <- as.data.frame(scale(climate))
climate.distance <- dist(climate, method = "euclidean")
dim(climate.distance)
#saveRDS(climate.distance, "data/climate.distance.rds")

climate_ITS <- dplyr::select(input_filt_rare_ITS$map_loaded, 
                             CMI, GDD, Temp, P, Tseas, Pseas)
climate_ITS <- as.data.frame(scale(climate_ITS))
climate.distance_ITS <- dist(climate_ITS, method = "euclidean")
dim(climate.distance_ITS)
#saveRDS(climate.distance_ITS, "data/climate.distance_ITS.rds")

# Plant
# Have 10 x 10 matrix, need to expand this to encompass the whole data set
gen_dist_mat <- readRDS("data/gen.dist.mat.rds")
gen_dist_mat[upper.tri(gen_dist_mat, diag = TRUE)] <- NA
gen_dist_df <- as.data.frame(gen_dist_mat)
gen_dist_df$sampleID <- rownames(gen_dist_df)
gen_dist_df_long <- reshape2::melt(gen_dist_df, id.vars = "sampleID")
levels(as.factor(gen_dist_df_long$sampleID))
gen_dist_df_long <- gen_dist_df_long %>%
  mutate(XY = paste(sampleID, variable, sep = "_"),
         YX = paste(variable, sampleID, sep = "_"))
# Note, the length of this dataframe should now equal (n*(n-1))/2 if NA omitted.
nrow(gen_dist_df_long)
nrow(gen_dist_df_long) == (10*(10-1))/2
for (i in 1:nrow(gen_dist_df_long)) {
  ifelse(gen_dist_df_long$sampleID[i] == gen_dist_df_long$variable[i],
         gen_dist_df_long$comparison[i] <- "within",
         gen_dist_df_long$comparison[i] <- "between")
}
gen_dist_df_long <- gen_dist_df_long %>%
  mutate(value = ifelse(comparison == "within",
                        0,
                        value))
levels(as.factor(gen_dist_df_long$sampleID))
levels(as.factor(gen_dist_df_long$variable))

# Put one of the other distance matrices into long format and then replace values
# Convert from dist object to matrix object
prep_mat <- as.matrix(climate.distance)
prep_mat[upper.tri(prep_mat, diag = TRUE)] <- NA
prep_df <- as.data.frame(prep_mat)
prep_df$sampleID <- rownames(prep_df)
prep_df_long <- reshape2::melt(prep_df, id.vars = "sampleID")
# Note, the length of this dataframe should now equal (n*(n-1))/2 if NA omitted
nrow(prep_df_long)
nrow(prep_df_long) == (587*(587-1))/2
levels(as.factor(prep_df_long$sampleID)) # 587
levels(as.factor(prep_df_long$variable)) # 587
input_filt_rare_16S$map_loaded$sampleID <- rownames(input_filt_rare_16S$map_loaded)
genotype_cols <- dplyr::select(input_filt_rare_16S$map_loaded, sampleID, Entry)
prep_df_long <- inner_join(prep_df_long, genotype_cols, 
                           by = c("sampleID" = "sampleID"))
prep_df_long <- inner_join(prep_df_long, genotype_cols, 
                           by = c("variable" = "sampleID"))
for (i in 1:nrow(prep_df_long)) {
  ifelse(prep_df_long$Entry.x[i] == prep_df_long$Entry.y[i],
         prep_df_long$comparison[i] <- "within",
         prep_df_long$comparison[i] <- "between")
}

# Sub out spaces to enable exact matching
table(prep_df_long$comparison)
prep_df_long <- prep_df_long %>%
  mutate(Entry.x = gsub(" ", "", Entry.x),
         Entry.y = gsub(" ", "", Entry.y)) %>%
  mutate(XY = paste(Entry.x, Entry.y, sep = "_"),
         YX = paste(Entry.y, Entry.x, sep = "_"))

# Check matches
sum(prep_df_long$XY %notin% gen_dist_df_long$XY)
sum(prep_df_long$XY %notin% gen_dist_df_long$YX)
sum(gen_dist_df_long$XY %notin% prep_df_long$XY)
sum(gen_dist_df_long$XY %notin% prep_df_long$YX)

# Need 4 left_joins to get everything!
prep_df_long <- prep_df_long %>%
  left_join(., gen_dist_df_long, by = "XY") %>%
  left_join(., gen_dist_df_long, by = c("XY" = "YX")) %>%
  left_join(., gen_dist_df_long, by = c("YX.x" = "XY")) %>%
  left_join(., gen_dist_df_long, by = c("YX.x" = "YX")) %>%
  mutate(gendist = ifelse(comparison.x == "within",
                          0,
                          value.y)) %>%
  mutate(gendist = ifelse(is.na(gendist) == TRUE,
                          value.x.x,
                          gendist)) %>%
  mutate(gendist = ifelse(is.na(gendist) == TRUE,
                          value.y.y,
                          gendist)) %>%
  mutate(gendist = ifelse(is.na(gendist) == TRUE,
                          value,
                          gendist)) %>%
  dplyr::select(sampleID.x, variable.x, gendist) %>%
  set_names(c("sampleID", "variable", "value"))
sum(is.na(prep_df_long$value)) # 0

# Now get back to wide format - make distance matrix
genetic.distance <- prep_df_long %>%
  pivot_wider(names_from = variable, values_from = value) %>%
  column_to_rownames(var = "sampleID")

#saveRDS(genetic.distance, "data/genetic.distance.rds")

# Check all dimensions. Should be 587 x 587
dim(geography.distance.mat)
dim(climate.distance)
dim(genetic.distance)
dim(bc_16S)

# Need to remake for ITS because 1 sample dropped!
# Copy code from above and rerun with climate.distance_ITS
prep_mat <- as.matrix(climate.distance_ITS)
prep_mat[upper.tri(prep_mat, diag = TRUE)] <- NA
prep_df <- as.data.frame(prep_mat)
prep_df$sampleID <- rownames(prep_df)
prep_df_long <- reshape2::melt(prep_df, id.vars = "sampleID")
input_filt_rare_ITS$map_loaded$sampleID <- rownames(input_filt_rare_ITS$map_loaded)
genotype_cols <- dplyr::select(input_filt_rare_ITS$map_loaded, sampleID, Entry)
prep_df_long <- inner_join(prep_df_long, genotype_cols, 
                           by = c("sampleID" = "sampleID"))
prep_df_long <- inner_join(prep_df_long, genotype_cols, 
                           by = c("variable" = "sampleID"))
for (i in 1:nrow(prep_df_long)) {
  ifelse(prep_df_long$Entry.x[i] == prep_df_long$Entry.y[i],
         prep_df_long$comparison[i] <- "within",
         prep_df_long$comparison[i] <- "between")
}
prep_df_long <- prep_df_long %>%
  mutate(Entry.x = gsub(" ", "", Entry.x),
         Entry.y = gsub(" ", "", Entry.y)) %>%
  mutate(XY = paste(Entry.x, Entry.y, sep = "_"),
         YX = paste(Entry.y, Entry.x, sep = "_")) %>%
  left_join(., gen_dist_df_long, by = "XY") %>%
  left_join(., gen_dist_df_long, by = c("XY" = "YX")) %>%
  left_join(., gen_dist_df_long, by = c("YX.x" = "XY")) %>%
  left_join(., gen_dist_df_long, by = c("YX.x" = "YX")) %>%
  mutate(gendist = ifelse(comparison.x == "within",
                          0,
                          value.y)) %>%
  mutate(gendist = ifelse(is.na(gendist) == TRUE,
                          value.x.x,
                          gendist)) %>%
  mutate(gendist = ifelse(is.na(gendist) == TRUE,
                          value.y.y,
                          gendist)) %>%
  mutate(gendist = ifelse(is.na(gendist) == TRUE,
                          value,
                          gendist)) %>%
  dplyr::select(sampleID.x, variable.x, gendist) %>%
  set_names(c("sampleID", "variable", "value"))
sum(is.na(prep_df_long$value)) # 0

# Now get back to wide format - make distance matrix
genetic.distance_ITS <- prep_df_long %>%
  pivot_wider(names_from = variable, values_from = value) %>%
  column_to_rownames(var = "sampleID")
#saveRDS(genetic.distance_ITS, "data/genetic.distance_ITS.rds")

# Soil distance (note just 1 number per site, just like climate)
soil <- dplyr::select(input_filt_rare_16S$map_loaded,
                      pH, OM_perc, Salts_mmhos_cm, NO3_ppm, P_ppm, K_ppm, Ca_ppm, Mg_ppm, 
                      Na_ppm, S_ppm, CEC_meq_g, Sand_perc, Silt_perc, Clay_perc)
soil <- as.data.frame(scale(soil))
soil.distance <- dist(soil, method = "euclidean")
dim(soil.distance)
#saveRDS(soil.distance, "data/soil.distance.rds")

soil_ITS <- dplyr::select(input_filt_rare_ITS$map_loaded,
                          pH, OM_perc, Salts_mmhos_cm, NO3_ppm, P_ppm, K_ppm, Ca_ppm, Mg_ppm, 
                          Na_ppm, S_ppm, CEC_meq_g, Sand_perc, Silt_perc, Clay_perc)
soil_ITS <- as.data.frame(scale(soil_ITS))
soil.distance_ITS <- dist(soil_ITS, method = "euclidean")
dim(soil.distance_ITS)
#saveRDS(soil.distance_ITS, "data/soil.distance_ITS.rds")

# Plant phenotype distance (note only chlorophyll has more data than just 1 per Entry)
plant <- dplyr::select(input_filt_rare_16S$map_loaded, 
                       chlorophyll, Sclerotinia, GRIN_Sclerotinia, DaysToFlower, StemLengthHighRange,
                       StemLengthLowRange, PerSeed_14_64, PerSeed_20_64, PerNoBranch, SeedWeight)
plant <- as.data.frame(scale(plant))
plant.distance <- dist(plant, method = "euclidean")
dim(plant.distance)
#saveRDS(plant.distance, "data/plant.distance.rds")

plant_ITS <- dplyr::select(input_filt_rare_ITS$map_loaded, 
                           chlorophyll, Sclerotinia, GRIN_Sclerotinia, DaysToFlower, StemLengthHighRange,
                           StemLengthLowRange, PerSeed_14_64, PerSeed_20_64, PerNoBranch, SeedWeight)
plant_ITS <- as.data.frame(scale(plant_ITS))
plant.distance_ITS <- dist(plant_ITS, method = "euclidean")
dim(plant.distance_ITS)
#saveRDS(plant.distance_ITS, "data/plant.distance_ITS.rds")



#### __Mantel ####
# Predictors and predictors
set.seed(402)
mantel(geography.distance.mat, climate.distance) # r = 0.96, p = 0.001
set.seed(402)
mantel(geography.distance.mat, soil.distance) # r = 0.36, p = 0.001
set.seed(402)
mantel(geography.distance.mat, plant.distance) # r = 0.07, p = 0.001
set.seed(402)
mantel(geography.distance.mat, genetic.distance) # r = -0.009, p = 0.904



# Predictors and community composition
set.seed(402)
mantel(bc_16S, geography.distance.mat) # r = 0.45, p = 0.001
mc_16S <- mantel.correlog(bc_16S, geography.distance.mat)
plot(mc_16S)
set.seed(402)
mantel(bc_16S, climate.distance) # r = 0.48, p = 0.001
set.seed(402)
mantel(bc_16S, soil.distance) # r = 0.69, p = 0.001
set.seed(402)
mantel(bc_16S, plant.distance) # r = 0.02, p = 0.001
set.seed(402)
mantel(bc_16S, genetic.distance) # r = -0.009, p = 0.935
set.seed(402)
mantel(bc_16S, genetic.distance, method = "spearman") # r = 0.002, p = 0.483

set.seed(402)
mantel(unifrac.distance, geography.distance.mat) # r = 0.36, p = 0.001
mc_16S <- mantel.correlog(unifrac.distance, geography.distance.mat)
plot(mc_16S)
set.seed(402)
mantel(unifrac.distance, climate.distance) # r = 0.37, p = 0.001
set.seed(402)
mantel(unifrac.distance, soil.distance) # r = 0.54, p = 0.001
set.seed(402)
mantel(unifrac.distance, plant.distance) # r = 0.01, p = 0.139
set.seed(402)
mantel(unifrac.distance, genetic.distance) # r = -0.004, p = 0.669
set.seed(402)
mantel(unifrac.distance, genetic.distance, method = "spearman") # r = 0.002, p = 0.436

set.seed(402)
mantel(bc_ITS, geography.distance.mat_ITS) # r = 0.54, p = 0.001
mc_ITS <- mantel.correlog(bc_ITS, geography.distance.mat_ITS)
plot(mc_ITS)
set.seed(402)
mantel(bc_ITS, climate.distance_ITS) # r = 0.59, p = 0.001
set.seed(402)
mantel(bc_ITS, soil.distance_ITS) # r = 0.63, p = 0.001
set.seed(402)
mantel(bc_ITS, plant.distance_ITS) # r = 0.05, p = 0.001
set.seed(402)
mantel(bc_ITS, genetic.distance_ITS) # r = -0.01, p = 0.964
set.seed(402)
mantel(bc_ITS, genetic.distance_ITS, method = "spearman") # r = -0.01, p = 0.699

# Partial test - test effects of climate and genetics while controlling for geography
set.seed(402)
mantel.partial(bc_16S, climate.distance, geography.distance.mat) # r = 0.19, p = 0.001
set.seed(402)
mantel.partial(bc_16S, soil.distance, geography.distance.mat) # r = 0.64, p = 0.001
set.seed(402)
mantel.partial(bc_16S, plant.distance, geography.distance.mat) # r = -0.01, p = 0.896
set.seed(402)
mantel.partial(bc_16S, genetic.distance, geography.distance.mat) # r = -0.005, p = 0.817

set.seed(402)
mantel.partial(unifrac.distance, climate.distance, geography.distance.mat) # r = 0.09, p = 0.001
set.seed(402)
mantel.partial(unifrac.distance, soil.distance, geography.distance.mat) # r = 0.47, p = 0.001
set.seed(402)
mantel.partial(unifrac.distance, plant.distance, geography.distance.mat) # r = -0.01, p = 0.896
set.seed(402)
mantel.partial(unifrac.distance, genetic.distance, geography.distance.mat) # r = -0.0004, p = 0.504

set.seed(402)
mantel.partial(bc_ITS, climate.distance_ITS, geography.distance.mat_ITS) # r = 0.33, p = 0.001
set.seed(402)
mantel.partial(bc_ITS, soil.distance_ITS, geography.distance.mat_ITS) # r = 0.55, p = 0.001
set.seed(402)
mantel.partial(bc_ITS, plant.distance_ITS, geography.distance.mat_ITS) # r = 0.02, p = 0.041
set.seed(402)
mantel.partial(bc_ITS, genetic.distance_ITS, geography.distance.mat_ITS) # r = -0.01, p = 0.913

# Quick Plot
png("InitialFigs/Geography_16S.png", width = 7, height = 5, units = "in", res = 300)
qplot(as.dist(geography.distance.mat), bc_16S, geom = c("point","smooth"), alpha = I(0.01)) +
  labs(x = "Geographic Distance (km)",
       y = "Bray-Curtis Dissimilarity") +
  ylim(0, 1) +
  theme_bw() +
  theme(axis.title = element_text(face = "bold", size = 12), 
        axis.text = element_text(size = 10),
        plot.margin = unit(c(0.1,0.1,0.1,0.15),"cm"))
dev.off()

png("InitialFigs/Geography_ITS.png", width = 7, height = 5, units = "in", res = 300)
qplot(as.dist(geography.distance.mat_ITS), bc_ITS, geom = c("point","smooth"), alpha = I(0.01)) +
  labs(x = "Geographic Distance (km)",
       y = "Bray-Curtis Dissimilarity") +
  ylim(0, 1) +
  theme_bw() +
  theme(axis.title = element_text(face = "bold", size = 12), 
        axis.text = element_text(size = 10),
        plot.margin = unit(c(0.1,0.1,0.1,0.15),"cm"))
dev.off()

png("InitialFigs/ClimateDistance_16S.png", width = 7, height = 5, units = "in", res = 300)
qplot(as.dist(climate.distance), bc_16S, geom = c("point","smooth"), alpha = I(0.01)) +
  labs(x = "Climate distance",
       y = "Bray-Curtis Dissimilarity") +
  ylim(0, 1) +
  theme_bw() +
  theme(axis.title = element_text(face = "bold", size = 12), 
        axis.text = element_text(size = 10),
        plot.margin = unit(c(0.1,0.1,0.1,0.15),"cm"))
dev.off()

png("InitialFigs/ClimateDistance_ITS.png", width = 7, height = 5, units = "in", res = 300)
qplot(as.dist(climate.distance_ITS), bc_ITS, geom = c("point","smooth"), alpha = I(0.01)) +
  labs(x = "Climate distance",
       y = "Bray-Curtis Dissimilarity") +
  ylim(0, 1) +
  theme_bw() +
  theme(axis.title = element_text(face = "bold", size = 12), 
        axis.text = element_text(size = 10),
        plot.margin = unit(c(0.1,0.1,0.1,0.15),"cm"))
dev.off()

png("InitialFigs/SoilDistance_16S.png", width = 7, height = 5, units = "in", res = 300)
qplot(as.dist(soil.distance), bc_16S, geom = c("point","smooth"), alpha = I(0.01)) +
  labs(x = "Soil distance",
       y = "Bray-Curtis Dissimilarity") +
  ylim(0, 1) +
  theme_bw() +
  theme(axis.title = element_text(face = "bold", size = 12), 
        axis.text = element_text(size = 10),
        plot.margin = unit(c(0.1,0.1,0.1,0.15),"cm"))
dev.off()

png("InitialFigs/SoilDistance_ITS.png", width = 7, height = 5, units = "in", res = 300)
qplot(as.dist(soil.distance_ITS), bc_ITS, geom = c("point","smooth"), alpha = I(0.01)) +
  labs(x = "Soil distance",
       y = "Bray-Curtis Dissimilarity") +
  ylim(0, 1) +
  theme_bw() +
  theme(axis.title = element_text(face = "bold", size = 12), 
        axis.text = element_text(size = 10),
        plot.margin = unit(c(0.1,0.1,0.1,0.15),"cm"))
dev.off()

png("InitialFigs/PlantDistance_16S.png", width = 7, height = 5, units = "in", res = 300)
qplot(as.dist(plant.distance), bc_16S, geom = c("point","smooth"), alpha = I(0.01)) +
  labs(x = "Plant distance",
       y = "Bray-Curtis Dissimilarity") +
  ylim(0, 1) +
  theme_bw() +
  theme(axis.title = element_text(face = "bold", size = 12), 
        axis.text = element_text(size = 10),
        plot.margin = unit(c(0.1,0.1,0.1,0.15),"cm"))
dev.off()

png("InitialFigs/PlantDistance_ITS.png", width = 7, height = 5, units = "in", res = 300)
qplot(as.dist(plant.distance_ITS), bc_ITS, geom = c("point","smooth"), alpha = I(0.01)) +
  labs(x = "Plant distance",
       y = "Bray-Curtis Dissimilarity") +
  ylim(0, 1) +
  theme_bw() +
  theme(axis.title = element_text(face = "bold", size = 12), 
        axis.text = element_text(size = 10),
        plot.margin = unit(c(0.1,0.1,0.1,0.15),"cm"))
dev.off()

png("InitialFigs/GeneticDistance_16S.png", width = 7, height = 5, units = "in", res = 300)
qplot(as.dist(genetic.distance), bc_16S, geom = c("point","smooth"), alpha = I(0.01)) +
  labs(x = "Plant genetic distance",
       y = "Bray-Curtis Dissimilarity") +
  ylim(0, 1) +
  theme_bw() +
  theme(axis.title = element_text(face = "bold", size = 12), 
        axis.text = element_text(size = 10),
        plot.margin = unit(c(0.1,0.1,0.1,0.15),"cm"))
dev.off()

png("InitialFigs/GeneticDistance_ITS.png", width = 7, height = 5, units = "in", res = 300)
qplot(as.dist(genetic.distance_ITS), bc_ITS, geom = c("point","smooth"), alpha = I(0.01)) +
  labs(x = "Plant genetic distance",
       y = "Bray-Curtis Dissimilarity") +
  ylim(0, 1) +
  theme_bw() +
  theme(axis.title = element_text(face = "bold", size = 12), 
        axis.text = element_text(size = 10),
        plot.margin = unit(c(0.1,0.1,0.1,0.15),"cm"))
dev.off()



#### __GDM ####
# Code adapted from Bueno de Mesquita et al. 2021, mSystems (Chimpanzee Gut Microbiomes)

# Distance matrices prepared in Mantel section
# Prepped here for GDM

# 16S
sum(rownames(as.matrix(bc_16S)) != input_filt_rare_16S$map_loaded$sampleID)
gdm.sampID <- rownames(as.matrix(bc_16S))
bacteria.distance.v.mat.gdm <- cbind(gdm.sampID, as.matrix(bc_16S))
geography.distance.v.mat.gdm <- cbind(gdm.sampID, geography.distance.mat)
climate.distance.v.mat.gdm <- cbind(gdm.sampID, as.matrix(climate.distance))
soil.distance.v.mat.gdm <- cbind(gdm.sampID, as.matrix(soil.distance))
plant.distance.v.mat.gdm <- cbind(gdm.sampID, as.matrix(plant.distance))
genetic.distance.mat.gdm <- cbind(gdm.sampID, genetic.distance)
dim(bacteria.distance.v.mat.gdm)
dim(geography.distance.v.mat.gdm)
dim(climate.distance.v.mat.gdm)
dim(soil.distance.v.mat.gdm)
dim(plant.distance.v.mat.gdm)
dim(genetic.distance.mat.gdm)
input_filt_rare_16S$map_loaded$gdm.sampID <- input_filt_rare_16S$map_loaded$sampleID
gdm.data <- dplyr::select(input_filt_rare_16S$map_loaded, gdm.sampID, LAT, LONG)
gdm.bray <- formatsitepair(bioData = bacteria.distance.v.mat.gdm,
                           bioFormat = 3,
                           siteColumn = "gdm.sampID",
                           XColumn = "LONG",
                           YColumn = "LAT",
                           predData = gdm.data,
                           distPreds = list(climate.distance.v.mat.gdm,
                                            soil.distance.v.mat.gdm,
                                            plant.distance.v.mat.gdm,
                                            genetic.distance.mat.gdm))
gdm.1 <- gdm(gdm.bray, geo = TRUE)
summary(gdm.1)
# Variable Importance
gdm.varImp(gdm.bray, geo = TRUE)
# Retained GEOG, 2 (soil), 4 (genetic)
# Percent deviance explained by final model = 55.973
# Final set of predictors returned: 
# Geographic
# matrix_1
# matrix_2
# matrix_4
# $`Model assessment`
# All predictors
# Model deviance                   5923.872
# Percent deviance explained         55.973
# Model p-value                       0.000
# Fitted permutations                50.000
# 
# $`Predictor Importance`
# All predictors
# Geographic          2.324
# matrix_1            0.294
# matrix_2           36.989
# matrix_4            0.004
# 
# $`Predictor p-values`
# All predictors
# Geographic              0
# matrix_1                0
# matrix_2                0
# matrix_4                0
# 
# $`Model Convergence`
# All predictors
# Geographic             50
# matrix_1               50
# matrix_2               50
# matrix_4               50
str(gdm.1)
length(gdm.1$predictors)
plot(gdm.1, plot.layout = c(2,3))
gdm.1.splineDat <- isplineExtract(gdm.1)
str(gdm.1.splineDat)
par(mfrow = c(1,5))
plot(gdm.1.splineDat$x[,"Geographic"], gdm.1.splineDat$y[,"Geographic"], lwd=3,
     type="l", xlab="Geographic distance", ylab="Partial ecological distance")
plot(gdm.1.splineDat$x[,"matrix_1"], gdm.1.splineDat$y[,"matrix_1"], lwd=3,
     type="l", xlab="Climate distance", ylab="Partial ecological distance")
plot(gdm.1.splineDat$x[,"matrix_2"], gdm.1.splineDat$y[,"matrix_2"], lwd=3,
     type="l", xlab="Soil distance", ylab="Partial ecological distance")
plot(gdm.1.splineDat$x[,"matrix_3"], gdm.1.splineDat$y[,"matrix_3"], lwd=3,
     type="l", xlab="Plant distance", ylab="Partial ecological distance")
plot(gdm.1.splineDat$x[,"matrix_4"], gdm.1.splineDat$y[,"matrix_4"], lwd=3,
     type="l", xlab="Genetic distance", ylab="Partial ecological distance")
max(gdm.1.splineDat$y[,"Geographic"])
max(gdm.1.splineDat$y[,"matrix_1"])
max(gdm.1.splineDat$y[,"matrix_2"])
max(gdm.1.splineDat$y[,"matrix_3"])
max(gdm.1.splineDat$y[,"matrix_4"])
gdm.1.pred <- predict(gdm.1, gdm.bray)
head(gdm.1.pred)
par(mfrow = c(1,1))
plot(gdm.bray$distance, gdm.1.pred, xlab="Observed dissimilarity",
     ylab="Predicted dissimilarity", xlim=c(0,1), ylim=c(0,1), pch=20, col=rgb(0,0,1,0.5))
lines(c(-1,2), c(-1,2))



# UniFrac
bacteria.distance.v.mat.gdm <- cbind(gdm.sampID, as.matrix(un_16S))
gdm.bray <- formatsitepair(bioData = bacteria.distance.v.mat.gdm,
                           bioFormat = 3,
                           siteColumn = "gdm.sampID",
                           XColumn = "LONG",
                           YColumn = "LAT",
                           predData = gdm.data,
                           distPreds = list(climate.distance.v.mat.gdm,
                                            soil.distance.v.mat.gdm,
                                            plant.distance.v.mat.gdm,
                                            genetic.distance.mat.gdm))
gdm.1 <- gdm(gdm.bray, geo = TRUE)
summary(gdm.1)
# Variable Importance
gdm.varImp(gdm.bray, geo = TRUE)
# Retained GEOG, 2 (soil), 4 (genetic)
# Percent deviance explained by final model = 36.149
# Final set of predictors returned: 
# Geographic
# matrix_1
# matrix_2
# matrix_4
# $`Model assessment`
# All predictors
# Model deviance                    849.328
# Percent deviance explained         36.149
# Model p-value                       0.000
# Fitted permutations                50.000
# 
# $`Predictor Importance`
# All predictors
# Geographic          4.974
# matrix_1            0.341
# matrix_2           36.478
# matrix_4            0.008
# 
# $`Predictor p-values`
# All predictors
# Geographic              0
# matrix_1                0
# matrix_2                0
# matrix_4                0
# 
# $`Model Convergence`
# All predictors
# Geographic             50
# matrix_1               50
# matrix_2               50
# matrix_4               50
str(gdm.1)
length(gdm.1$predictors)
plot(gdm.1, plot.layout = c(2,3))
gdm.1.splineDat <- isplineExtract(gdm.1)
str(gdm.1.splineDat)
par(mfrow = c(1,5))
plot(gdm.1.splineDat$x[,"Geographic"], gdm.1.splineDat$y[,"Geographic"], lwd=3,
     type="l", xlab="Geographic distance", ylab="Partial ecological distance")
plot(gdm.1.splineDat$x[,"matrix_1"], gdm.1.splineDat$y[,"matrix_1"], lwd=3,
     type="l", xlab="Climate distance", ylab="Partial ecological distance")
plot(gdm.1.splineDat$x[,"matrix_2"], gdm.1.splineDat$y[,"matrix_2"], lwd=3,
     type="l", xlab="Soil distance", ylab="Partial ecological distance")
plot(gdm.1.splineDat$x[,"matrix_3"], gdm.1.splineDat$y[,"matrix_3"], lwd=3,
     type="l", xlab="Plant distance", ylab="Partial ecological distance")
plot(gdm.1.splineDat$x[,"matrix_4"], gdm.1.splineDat$y[,"matrix_4"], lwd=3,
     type="l", xlab="Genetic distance", ylab="Partial ecological distance")
max(gdm.1.splineDat$y[,"Geographic"])
max(gdm.1.splineDat$y[,"matrix_1"])
max(gdm.1.splineDat$y[,"matrix_2"])
max(gdm.1.splineDat$y[,"matrix_3"])
max(gdm.1.splineDat$y[,"matrix_4"])
gdm.1.pred <- predict(gdm.1, gdm.bray)
head(gdm.1.pred)
par(mfrow = c(1,1))
plot(gdm.bray$distance, gdm.1.pred, xlab="Observed dissimilarity",
     ylab="Predicted dissimilarity", xlim=c(0,1), ylim=c(0,1), pch=20, col=rgb(0,0,1,0.5))
lines(c(-1,2), c(-1,2))



# ITS
sum(rownames(as.matrix(bc_ITS)) != input_filt_rare_ITS$map_loaded$sampleID)
gdm.sampID <- rownames(as.matrix(bc_ITS))
fungi.distance.v.mat.gdm <- cbind(gdm.sampID, as.matrix(bc_ITS))
geography.distance.v.mat.gdm <- cbind(gdm.sampID, geography.distance.mat_ITS)
climate.distance.v.mat.gdm <- cbind(gdm.sampID, as.matrix(climate.distance_ITS))
soil.distance.v.mat.gdm <- cbind(gdm.sampID, as.matrix(soil.distance_ITS))
plant.distance.v.mat.gdm <- cbind(gdm.sampID, as.matrix(plant.distance_ITS))
genetic.distance.mat.gdm <- cbind(gdm.sampID, genetic.distance_ITS)
dim(fungi.distance.v.mat.gdm)
dim(geography.distance.v.mat.gdm)
dim(climate.distance.v.mat.gdm)
dim(soil.distance.v.mat.gdm)
dim(plant.distance.v.mat.gdm)
dim(genetic.distance.mat.gdm)
input_filt_rare_ITS$map_loaded$gdm.sampID <- input_filt_rare_ITS$map_loaded$sampleID
gdm.data <- dplyr::select(input_filt_rare_ITS$map_loaded, gdm.sampID, LAT, LONG)
gdm.bray <- formatsitepair(bioData = fungi.distance.v.mat.gdm,
                           bioFormat = 3,
                           siteColumn = "gdm.sampID",
                           XColumn = "LONG",
                           YColumn = "LAT",
                           predData = gdm.data,
                           distPreds = list(climate.distance.v.mat.gdm,
                                            soil.distance.v.mat.gdm,
                                            plant.distance.v.mat.gdm,
                                            genetic.distance.mat.gdm))
gdm.2 <- gdm(gdm.bray, geo = TRUE)
summary(gdm.2)
# Variable Importance
gdm.varImp(gdm.bray, geo = TRUE)
# Percent deviance explained by the full model =  57.632
# Fitting GDMs to the permuted site-pair tables...
# Assessing importance of geographic distance...
# Assessing importance of matrix_1...
# Assessing importance of matrix_2...
# Assessing importance of matrix_3...
# All remaining predictors are significant, ceasing assessment.
# Percent deviance explained by final model = 57.632
# Final set of predictors returned: 
# Geographic
# matrix_1
# matrix_2
# matrix_3
# $`Model assessment`
# All predictors
# Model deviance                   5107.440
# Percent deviance explained         57.632
# Model p-value                       0.000
# Fitted permutations                50.000
# 
# $`Predictor Importance`
# All predictors
# Geographic          0.305
# matrix_1            2.982
# matrix_2           19.520
# matrix_3            0.023
# 
# $`Predictor p-values`
# All predictors
# Geographic              0
# matrix_1                0
# matrix_2                0
# matrix_3                0
# 
# $`Model Convergence`
# All predictors
# Geographic             50
# matrix_1               50
# matrix_2               50
# matrix_3               50
str(gdm.2)
length(gdm.2$predictors)
plot(gdm.2, plot.layout = c(2,3))
gdm.2.splineDat <- isplineExtract(gdm.2)
str(gdm.2.splineDat)
par(mfrow = c(1,5))
plot(gdm.2.splineDat$x[,"Geographic"], gdm.2.splineDat$y[,"Geographic"], lwd=3,
     type="l", xlab="Geographic distance", ylab="Partial ecological distance")
plot(gdm.2.splineDat$x[,"matrix_1"], gdm.2.splineDat$y[,"matrix_1"], lwd=3,
     type="l", xlab="Climate distance", ylab="Partial ecological distance")
plot(gdm.2.splineDat$x[,"matrix_2"], gdm.2.splineDat$y[,"matrix_2"], lwd=3,
     type="l", xlab="Soil distance", ylab="Partial ecological distance")
plot(gdm.2.splineDat$x[,"matrix_3"], gdm.2.splineDat$y[,"matrix_3"], lwd=3,
     type="l", xlab="Plant distance", ylab="Partial ecological distance")
plot(gdm.2.splineDat$x[,"matrix_4"], gdm.2.splineDat$y[,"matrix_4"], lwd=3,
     type="l", xlab="Genetic distance", ylab="Partial ecological distance")
max(gdm.2.splineDat$y[,"Geographic"])
max(gdm.2.splineDat$y[,"matrix_1"])
max(gdm.2.splineDat$y[,"matrix_2"])
max(gdm.2.splineDat$y[,"matrix_3"])
max(gdm.2.splineDat$y[,"matrix_4"])
gdm.2.pred <- predict(gdm.2, gdm.bray)
head(gdm.2.pred)
par(mfrow = c(1,1))
plot(gdm.bray$distance, gdm.2.pred, xlab="Observed dissimilarity",
     ylab="Predicted dissimilarity", xlim=c(0,1), ylim=c(0,1), pch=20, col=rgb(0,0,1,0.5))
lines(c(-1,2), c(-1,2))



#### __Predictors ####
pcoa_clim <- cmdscale(climate.distance, k = nrow(input_filt_rare_16S$map_loaded) - 1, eig = T)
pcoaA1 <- paste("PC1: ", round((eigenvals(pcoa_clim)/sum(eigenvals(pcoa_clim)))[1]*100, 1), "%")
pcoaA2 <- paste("PC2: ", round((eigenvals(pcoa_clim)/sum(eigenvals(pcoa_clim)))[2]*100, 1), "%")
input_filt_rare_16S$map_loaded$Axis01 <- vegan::scores(pcoa_clim)[,1]
input_filt_rare_16S$map_loaded$Axis02 <- vegan::scores(pcoa_clim)[,2]
micro.hulls <- ddply(input_filt_rare_16S$map_loaded, c("Site"), find_hull)
ggplot(input_filt_rare_16S$map_loaded, aes(Axis01, Axis02)) +
  geom_point(size = 3, pch = 16, alpha = 0.75, aes(colour = Site)) +
  labs(x = pcoaA1, 
       y = pcoaA2,
       colour = "Site") +
  scale_colour_viridis_d() +
  theme_bw() +  
  theme(legend.position = "none",
        axis.title = element_text(face = "bold", size = 12), 
        axis.text = element_text(size = 10),
        plot.title = element_text(vjust = -1))

pcoa_soil <- cmdscale(soil.distance, k = nrow(input_filt_rare_16S$map_loaded) - 1, eig = T)
pcoaA1 <- paste("PC1: ", round((eigenvals(pcoa_soil)/sum(eigenvals(pcoa_soil)))[1]*100, 1), "%")
pcoaA2 <- paste("PC2: ", round((eigenvals(pcoa_soil)/sum(eigenvals(pcoa_soil)))[2]*100, 1), "%")
input_filt_rare_16S$map_loaded$Axis01 <- vegan::scores(pcoa_soil)[,1]
input_filt_rare_16S$map_loaded$Axis02 <- vegan::scores(pcoa_soil)[,2]
micro.hulls <- ddply(input_filt_rare_16S$map_loaded, c("Site"), find_hull)
ggplot(input_filt_rare_16S$map_loaded, aes(Axis01, Axis02)) +
  geom_point(size = 3, pch = 16, alpha = 0.75, aes(colour = Site)) +
  labs(x = pcoaA1, 
       y = pcoaA2,
       colour = "Site") +
  scale_colour_viridis_d() +
  theme_bw() +  
  theme(legend.position = "none",
        axis.title = element_text(face = "bold", size = 12), 
        axis.text = element_text(size = 10),
        plot.title = element_text(vjust = -1))

pcoa_plant <- cmdscale(plant.distance, k = nrow(input_filt_rare_16S$map_loaded) - 1, eig = T)
pcoaA1 <- paste("PC1: ", round((eigenvals(pcoa_plant)/sum(eigenvals(pcoa_plant)))[1]*100, 1), "%")
pcoaA2 <- paste("PC2: ", round((eigenvals(pcoa_plant)/sum(eigenvals(pcoa_plant)))[2]*100, 1), "%")
input_filt_rare_16S$map_loaded$Axis01 <- vegan::scores(pcoa_plant)[,1]
input_filt_rare_16S$map_loaded$Axis02 <- vegan::scores(pcoa_plant)[,2]
micro.hulls <- ddply(input_filt_rare_16S$map_loaded, c("Site"), find_hull)
ggplot(input_filt_rare_16S$map_loaded, aes(Axis01, Axis02)) +
  geom_point(size = 3, pch = 16, alpha = 0.75, aes(colour = Site)) +
  labs(x = pcoaA1, 
       y = pcoaA2,
       colour = "Site") +
  scale_colour_viridis_d() +
  theme_bw() +  
  theme(legend.position = "none",
        axis.title = element_text(face = "bold", size = 12), 
        axis.text = element_text(size = 10),
        plot.title = element_text(vjust = -1))



#### 5. Taxa ####
# Check number of ASVs with assigned genera
nrow(input_filt_rare_16S$taxonomy_loaded)
sum(input_filt_rare_16S$taxonomy_loaded$taxonomy6 == "NA")
17619/27893 # 63%

nrow(input_filt_rare_ITS$taxonomy_loaded)
sum(input_filt_rare_ITS$taxonomy_loaded$taxonomy6 == "NA")
1263/5294 # 24%


#### _Indicators ####
# Multipatt
# Run on Innes server to speed up. Also just consider single groups.
# Load results - ASV level
mp_site_16S <- readRDS("data/mp_site_16S.rds")
mp_site_ITS <- readRDS("data/mp_site_ITS.rds")
mp_Entry_16S <- readRDS("data/mp_Entry_16S.rds")
mp_Entry_ITS <- readRDS("data/mp_Entry_ITS.rds")
input_filt_rare_16S$taxonomy_loaded$taxonomy8 <- input_filt_rare_16S$taxonomy_loaded$taxonomy7
input_filt_rare_ITS$taxonomy_loaded$taxonomy9 <- input_filt_rare_ITS$taxonomy_loaded$taxonomy7
input_filt_rare_ITS$taxonomy_loaded$taxonomy7 <- input_filt_rare_ITS$taxonomy_loaded$taxonomy8
input_filt_rare_ITS$taxonomy_loaded$taxonomy8 <- input_filt_rare_ITS$taxonomy_loaded$taxonomy9
input_filt_rare_ITS$taxonomy_loaded$taxonomy9 <- NULL
input_filt_rare_16S$map_loaded$Site <- as.factor(input_filt_rare_16S$map_loaded$Site)
input_filt_rare_ITS$map_loaded$Site <- as.factor(input_filt_rare_ITS$map_loaded$Site)
input_filt_rare_16S$map_loaded$Entry <- as.factor(input_filt_rare_16S$map_loaded$Entry)
input_filt_rare_ITS$map_loaded$Entry <- as.factor(input_filt_rare_ITS$map_loaded$Entry)

# Adjust ITS taxonomy
View(input_filt_rare_ITS$taxonomy_loaded)
input_filt_rare_ITS$taxonomy_loaded <- input_filt_rare_ITS$taxonomy_loaded %>%
  mutate(taxonomy1 = gsub("k__", "", taxonomy1),
         taxonomy2 = gsub("p__", "", taxonomy2),
         taxonomy3 = gsub("c__", "", taxonomy3),
         taxonomy4 = gsub("o__", "", taxonomy4),
         taxonomy5 = gsub("f__", "", taxonomy5),
         taxonomy6 = gsub("g__", "", taxonomy6),
         taxonomy7 = gsub("s__", "", taxonomy7))

# Add OTU_IDs. Make it taxonomy9. Make sure ASVs are row names.
View(input_filt_rare_16S$taxonomy_loaded)
View(input_filt_rare_ITS$taxonomy_loaded)
input_filt_rare_16S$taxonomy_loaded <- input_filt_rare_16S$taxonomy_loaded %>%
  left_join(., otu_asv_map_16S, by = c("taxonomy8" = "ASV_ID"))
rownames(input_filt_rare_16S$taxonomy_loaded) <- input_filt_rare_16S$taxonomy_loaded$taxonomy8
input_filt_rare_ITS$taxonomy_loaded <- input_filt_rare_ITS$taxonomy_loaded %>%
  left_join(., otu_asv_map_ITS, by = c("taxonomy8" = "ASV_ID"))
rownames(input_filt_rare_ITS$taxonomy_loaded) <- input_filt_rare_ITS$taxonomy_loaded$taxonomy8



#### __Site ####
# Indicators of site with all 10 genotypes as input

### 16S 
tax_sum_asv_16S <- summarize_taxonomy(input_filt_rare_16S, 
                                      level = 7, 
                                      report_higher_tax = F,
                                      relative = T)
# (faster way)
input_filt_rare_16S_rel <- convert_to_relative_abundances(input_filt_rare_16S)
tax_sum_asv_16S <- input_filt_rare_16S_rel$data_loaded

set.seed(425)
mp_site_16S <- multipatt(t(input_filt_rare_16S$data_loaded), 
                         input_filt_rare_16S$map_loaded$Site, 
                         func = "r.g", 
                         duleg = TRUE,
                         control = how(nperm=999))
summary(mp_site_16S) # Number of species associated to 1 group: 19642 (p.fdr < 0.01)

# Plot significant and high r 
png("InitialFigs/Multipatt_Site_16S_r90.png", width = 6, height = 8, units = "in", res = 300)
plot_multipatt_asv(mp_obj = mp_site_16S, 
                   input = input_filt_rare_16S,
                   tax_sum = tax_sum_asv_16S,
                   group = "Site",
                   filter = FALSE,
                   #filter_vals = "Field",
                   abund = "% Rel. Abund.",
                   qcut = 0.05,
                   rcut = 0.9)
dev.off()

# Get list of top
mp_site_16S_sig_top <- mp_site_16S$sign %>%
  rename("s.Brighton" = "s.Kirkmeyer") %>%
  mutate(p.fdr = p.adjust(p.value, "fdr")) %>%
  filter(p.fdr < 0.05) %>%
  filter(stat > 0.9) %>%
  mutate("ASV_ID" = rownames(.)) %>%
  left_join(., input_filt_rare_16S$taxonomy_loaded, by = c("ASV_ID" = "taxonomy8")) %>%
  set_names(c(substr(names(.)[1:15], 3, nchar(.))), names(.)[16:ncol(.)]) %>%
  arrange(across(names(.)[1:15], desc)) 
  
# Plot number of phyla that are indicators
mp_site_16S_sig <- mp_site_16S$sign %>%
  rename("s.Brighton" = "s.Kirkmeyer") %>%
  mutate(p.fdr = p.adjust(p.value, "fdr")) %>%
  filter(p.fdr < 0.01) %>%
  mutate("ASV_ID" = rownames(.)) %>%
  left_join(., input_filt_rare_16S$taxonomy_loaded, by = c("ASV_ID" = "taxonomy8")) %>%
  set_names(c(substr(names(.)[1:15], 3, nchar(.))), names(.)[16:ncol(.)])
top_taxa_16S <- as.data.frame(table(mp_site_16S_sig$taxonomy2)) %>%
  arrange(desc(Freq)) %>%
  slice_head(n = 12) %>%
  mutate(Var1 = as.character(Var1))
mp_site_16S_sig_toPlot <- mp_site_16S_sig %>%
  pivot_longer(cols = names(.)[1:15]) %>%
  filter(value > 0) %>%
  group_by(name, taxonomy2) %>%
  summarise(n_taxa = n()) %>%
  drop_na(taxonomy2)
mp_site_16S_sig_toPlot$taxonomy2[mp_site_16S_sig_toPlot$taxonomy2 %notin% top_taxa_16S$Var1] = 'Other'
mp_site_16S_sig_toPlot <- mp_site_16S_sig_toPlot %>%
  group_by(name, taxonomy2) %>%
  summarise(n_taxa = sum(n_taxa)) %>%
  mutate(taxonomy2 = factor(taxonomy2,
                            levels = c("Other", rev(top_taxa_16S$Var1))))
n_per_site_16S <- mp_site_16S_sig_toPlot %>%
  group_by(name) %>%
  summarise(site_n = sum(n_taxa))
pdf("InitialFigs/Multipatt_Site_16S_nTax.pdf", width = 7, height = 5)
ggplot(mp_site_16S_sig_toPlot, aes(name, n_taxa, fill = taxonomy2, group = taxonomy2)) +
  geom_bar(stat = "identity", colour = NA, linewidth = 0.25) +
  labs(x = "Site",
       y = "Number of indicator ASVs",
       fill = "Phylum") +
  scale_fill_manual(values = c("grey90", brewer.pal(12, "Paired")[12:1])) +
  scale_y_continuous(expand = c(0.01, 0.01)) +  
  theme_classic() +
  theme(axis.title.y = element_text(size = 12),
        axis.title.x = element_blank(),
        axis.text.y = element_text(size = 10),
        axis.text.x = element_text(size = 10, angle = 45, hjust = 1))
dev.off()

fig5a <- ggplot(mp_site_16S_sig_toPlot, aes(name, n_taxa, fill = taxonomy2, group = taxonomy2)) +
  geom_bar(stat = "identity", colour = NA, linewidth = 0.25) +
  labs(x = "Site",
       y = "# of indicator ASVs",
       fill = "Phylum") +
  scale_fill_manual(values = c("grey90", brewer.pal(12, "Paired")[12:1])) +
  scale_y_continuous(expand = c(0.01, 0.01)) +  
  theme_classic() +
  theme(axis.title.y = element_text(size = 12),
        axis.title.x = element_blank(),
        axis.text.y = element_text(size = 10),
        axis.text.x = element_blank(),
        legend.margin = margin(0, 0, 0, 0, unit = "pt"),
        legend.box.margin = margin(0, 0, 0, 0, unit = "pt"),
        legend.key.size = unit(0.3, "cm"),
        legend.text = element_text(size = 6))
fig5a
leg5a <- get_legend(fig5a)
fig5a <- fig5a + theme(legend.position = "none")

# Rerun at OTU level
# Make a new input file with OTU_ID as taxonomy 8
input_filt_rare_16S_indic <- input_filt_rare_16S
input_filt_rare_16S_indic$taxonomy_loaded$taxonomy8 <- input_filt_rare_16S_indic$taxonomy_loaded$OTU_ID
set.seed(425)
mp_site_16S <- multipatt(t(tax_sum_OTU_16S), 
                         input_filt_rare_16S$map_loaded$Site, 
                         func = "r.g", 
                         duleg = TRUE,
                         control = how(nperm=999))
summary(mp_site_16S) # Number of species associated to 1 group: 7959 (p.fdr < 0.01)

# Plot significant and high r 
png("InitialFigs/Multipatt_Site_16S_OTU_r90.png", width = 6, height = 8, units = "in", res = 300)
plot_multipatt_asv(mp_obj = mp_site_16S, 
                   input = input_filt_rare_16S_indic,
                   tax_sum = tax_sum_OTU_16S,
                   group = "Site",
                   filter = FALSE,
                   #filter_vals = "Field",
                   abund = "% Rel. Abund.",
                   qcut = 0.05,
                   rcut = 0.9)
dev.off()

# Get list of top
mp_site_16S_sig_top <- mp_site_16S$sign %>%
  mutate(p.fdr = p.adjust(p.value, "fdr")) %>%
  filter(p.fdr < 0.05) %>%
  filter(stat > 0.9) %>%
  mutate("ASV_ID" = rownames(.)) %>%
  left_join(., otu_tax_16S, by = c("ASV_ID" = "OTU_ID")) %>%
  set_names(c(substr(names(.)[1:15], 3, nchar(.))), names(.)[16:ncol(.)]) %>%
  arrange(across(names(.)[1:15], desc)) 

# Plot number of phyla that are indicators
mp_site_16S_sig <- mp_site_16S$sign %>%
  mutate(p.fdr = p.adjust(p.value, "fdr")) %>%
  filter(p.fdr < 0.01) %>%
  mutate("ASV_ID" = rownames(.)) %>%
  left_join(., otu_tax_16S, by = c("ASV_ID" = "OTU_ID")) %>%
  set_names(c(substr(names(.)[1:15], 3, nchar(.))), names(.)[16:ncol(.)])
top_taxa_16S <- as.data.frame(table(mp_site_16S_sig$taxonomy2)) %>%
  arrange(desc(Freq)) %>%
  slice_head(n = 12) %>%
  mutate(Var1 = as.character(Var1))
mp_site_16S_sig_toPlot <- mp_site_16S_sig %>%
  pivot_longer(cols = names(.)[1:15]) %>%
  filter(value > 0) %>%
  group_by(name, taxonomy2) %>%
  summarise(n_taxa = n()) %>%
  drop_na(taxonomy2)
mp_site_16S_sig_toPlot$taxonomy2[mp_site_16S_sig_toPlot$taxonomy2 %notin% top_taxa_16S$Var1] = 'Other'
mp_site_16S_sig_toPlot <- mp_site_16S_sig_toPlot %>%
  group_by(name, taxonomy2) %>%
  summarise(n_taxa = sum(n_taxa)) %>%
  mutate(taxonomy2 = factor(taxonomy2,
                            levels = c("Other", rev(top_taxa_16S$Var1))))
pdf("InitialFigs/Multipatt_Site_16S_OTU_nTax.pdf", width = 7, height = 5)
ggplot(mp_site_16S_sig_toPlot, aes(name, n_taxa, fill = taxonomy2, group = taxonomy2)) +
  geom_bar(stat = "identity", colour = NA, linewidth = 0.25) +
  labs(x = "Site",
       y = "Number of indicator OTUs",
       fill = "Phylum") +
  scale_fill_manual(values = c("grey90", brewer.pal(12, "Paired")[12:1])) +
  scale_y_continuous(expand = c(0.01, 0.01)) +  
  theme_classic() +
  theme(axis.title.y = element_text(size = 12),
        axis.title.x = element_blank(),
        axis.text.y = element_text(size = 10),
        axis.text.x = element_text(size = 10, angle = 45, hjust = 1))
dev.off()



### ITS
tax_sum_asv_ITS <- summarize_taxonomy(input_filt_rare_ITS, 
                                      level = 8, 
                                      report_higher_tax = F, 
                                      relative = T)
set.seed(425)
mp_site_ITS <- multipatt(t(input_filt_rare_ITS$data_loaded), 
                         input_filt_rare_ITS$map_loaded$Site, 
                         func = "r.g", 
                         duleg = TRUE,
                         control = how(nperm=999))
summary(mp_site_ITS) # Number of species associated to 1 group: 2089 (p.fdr < 0.05)

# Plot significant and high r 
png("InitialFigs/Multipatt_Site_ITS_r65.png", width = 6, height = 8, units = "in", res = 300)
plot_multipatt_asv(mp_obj = mp_site_ITS, 
                   input = input_filt_rare_ITS,
                   tax_sum = tax_sum_asv_ITS,
                   group = "Site",
                   filter = FALSE,
                   #filter_vals = "Field",
                   abund = "% Rel. Abund.",
                   qcut = 0.05,
                   rcut = 0.65)
dev.off()

# Get list of top
mp_site_ITS_sig_top <- mp_site_ITS$sign %>%
  rename("s.Brighton" = "s.Kirkmeyer") %>%
  mutate(p.fdr = p.adjust(p.value, "fdr")) %>%
  filter(p.fdr < 0.05) %>%
  filter(stat > 0.65) %>%
  mutate("ASV_ID" = rownames(.)) %>%
  left_join(., input_filt_rare_16S$taxonomy_loaded, by = c("ASV_ID" = "taxonomy8")) %>%
  set_names(c(substr(names(.)[1:15], 3, nchar(.))), names(.)[16:ncol(.)]) %>%
  arrange(across(names(.)[1:15], desc)) 

# Plot number of classes that are indicators
mp_site_ITS_sig <- mp_site_ITS$sign %>%
  rename("s.Brighton" = "s.Kirkmeyer") %>%
  mutate(p.fdr = p.adjust(p.value, "fdr")) %>%
  filter(p.fdr < 0.05) %>%
  mutate("ASV_ID" = rownames(.)) %>%
  left_join(., input_filt_rare_ITS$taxonomy_loaded, by = c("ASV_ID" = "taxonomy8")) %>%
  set_names(c(substr(names(.)[1:15], 3, nchar(.))), names(.)[16:ncol(.)])
top_taxa_ITS <- as.data.frame(table(mp_site_ITS_sig$taxonomy3)) %>%
  arrange(desc(Freq)) %>%
  slice_head(n = 12) %>%
  mutate(Var1 = as.character(Var1))
mp_site_ITS_sig_toPlot <- mp_site_ITS_sig %>%
  pivot_longer(cols = names(.)[1:15]) %>%
  filter(value > 0) %>%
  group_by(name, taxonomy3) %>%
  summarise(n_taxa = n()) %>%
  drop_na(taxonomy3)
mp_site_ITS_sig_toPlot$taxonomy3[mp_site_ITS_sig_toPlot$taxonomy3 %notin% top_taxa_ITS$Var1] = 'Other'
mp_site_ITS_sig_toPlot <- mp_site_ITS_sig_toPlot %>%
  group_by(name, taxonomy3) %>%
  summarise(n_taxa = sum(n_taxa)) %>%
  mutate(taxonomy3 = factor(taxonomy3,
                            levels = c("Other", "NA",
                                       rev(subset(top_taxa_ITS$Var1, top_taxa_ITS$Var1 != "NA")))))
n_per_site_ITS <- mp_site_ITS_sig_toPlot %>%
  group_by(name) %>%
  summarise(site_n = sum(n_taxa))
pdf("InitialFigs/Multipatt_Site_ITS_nTax.pdf", width = 7, height = 5)
ggplot(mp_site_ITS_sig_toPlot, aes(name, n_taxa, fill = taxonomy3, group = taxonomy3)) +
  geom_bar(stat = "identity", colour = NA, linewidth = 0.25) +
  labs(x = "Site",
       y = "Number of indicator ASVs",
       fill = "Class") +
  scale_fill_manual(values = c("grey90", "grey75", brewer.pal(11, "Paired")[11:1])) +
  scale_y_continuous(expand = c(0.01, 0.01)) +  
  theme_classic() +
  theme(axis.title.y = element_text(size = 12),
        axis.title.x = element_blank(),
        axis.text.y = element_text(size = 10),
        axis.text.x = element_text(size = 10, angle = 45, hjust = 1))
dev.off()

fig5b <- ggplot(mp_site_ITS_sig_toPlot, aes(name, n_taxa, fill = taxonomy3, group = taxonomy3)) +
  geom_bar(stat = "identity", colour = NA, linewidth = 0.25) +
  labs(x = "Site",
       y = "# of indicator ASVs",
       fill = "Class") +
  scale_fill_manual(values = c("grey90", "grey75", brewer.pal(11, "Paired")[11:1])) +
  scale_y_continuous(expand = c(0.01, 0.01)) +  
  theme_classic() +
  theme(axis.title.y = element_text(size = 12),
        axis.title.x = element_blank(),
        axis.text.y = element_text(size = 10),
        axis.text.x = element_text(size = 10, angle = 45, hjust = 1),
        legend.margin = margin(0, 0, 0, 0, unit = "pt"),
        legend.box.margin = margin(0, 0, 0, 0, unit = "pt"),
        legend.key.size = unit(0.3, "cm"),
        legend.text = element_text(size = 6))
fig5b
leg5b <- get_legend(fig5b)
fig5b <- fig5b + theme(legend.position = "none")

left <- plot_grid(fig5a, fig5b, ncol = 1, align = "v", rel_heights = c(0.45, 0.55),
                  labels = "auto")
right <- plot_grid(leg5a, leg5b, NULL, ncol = 1, align = "v", rel_heights = c(0.4, 0.5, 0.1))
pdf("FinalFigs/Figure5.pdf", width = 7, height = 6)
plot_grid(left, right, rel_widths = c(0.8, 0.2))
dev.off()



# Rerun at OTU level
# Make a new input file with OTU_ID as taxonomy 8
input_filt_rare_ITS_indic <- input_filt_rare_ITS
input_filt_rare_ITS_indic$taxonomy_loaded$taxonomy8 <- input_filt_rare_ITS_indic$taxonomy_loaded$OTU_ID
set.seed(425)
mp_site_ITS <- multipatt(t(tax_sum_OTU_ITS), 
                         input_filt_rare_ITS$map_loaded$Site, 
                         func = "r.g", 
                         duleg = TRUE,
                         control = how(nperm=999))
summary(mp_site_ITS) # Number of species associated to 1 group: 1295 (p.fdr < 0.05)

# Plot significant and high r 
png("InitialFigs/Multipatt_Site_ITS_OTU_r65.png", width = 6, height = 8, units = "in", res = 300)
plot_multipatt_asv(mp_obj = mp_site_ITS, 
                   input = input_filt_rare_ITS_indic,
                   tax_sum = tax_sum_OTU_ITS,
                   group = "Site",
                   filter = FALSE,
                   #filter_vals = "Field",
                   abund = "% Rel. Abund.",
                   qcut = 0.05,
                   rcut = 0.65)
dev.off()

# Get list of top
mp_site_ITS_sig_top <- mp_site_ITS$sign %>%
  mutate(p.fdr = p.adjust(p.value, "fdr")) %>%
  filter(p.fdr < 0.05) %>%
  filter(stat > 0.65) %>%
  mutate("ASV_ID" = rownames(.)) %>%
  left_join(., otu_tax_ITS, by = c("ASV_ID" = "OTU_ID")) %>%
  set_names(c(substr(names(.)[1:15], 3, nchar(.))), names(.)[16:ncol(.)]) %>%
  arrange(across(names(.)[1:15], desc)) 

# Plot number of classes that are indicators
mp_site_ITS_sig <- mp_site_ITS$sign %>%
  mutate(p.fdr = p.adjust(p.value, "fdr")) %>%
  filter(p.fdr < 0.05) %>%
  mutate("ASV_ID" = rownames(.)) %>%
  left_join(., otu_tax_ITS, by = c("ASV_ID" = "OTU_ID")) %>%
  set_names(c(substr(names(.)[1:15], 3, nchar(.))), names(.)[16:ncol(.)])
top_taxa_ITS <- as.data.frame(table(mp_site_ITS_sig$taxonomy3)) %>%
  arrange(desc(Freq)) %>%
  slice_head(n = 12) %>%
  mutate(Var1 = as.character(Var1))
mp_site_ITS_sig_toPlot <- mp_site_ITS_sig %>%
  pivot_longer(cols = names(.)[1:15]) %>%
  filter(value > 0) %>%
  group_by(name, taxonomy3) %>%
  summarise(n_taxa = n()) %>%
  drop_na(taxonomy3)
mp_site_ITS_sig_toPlot$taxonomy3[mp_site_ITS_sig_toPlot$taxonomy3 %notin% top_taxa_ITS$Var1] = 'Other'
mp_site_ITS_sig_toPlot <- mp_site_ITS_sig_toPlot %>%
  group_by(name, taxonomy3) %>%
  summarise(n_taxa = sum(n_taxa)) %>%
  mutate(taxonomy3 = factor(taxonomy3,
                            levels = c("Other", "NA",
                                       rev(subset(top_taxa_ITS$Var1, top_taxa_ITS$Var1 != "NA")))))
pdf("InitialFigs/Multipatt_Site_ITS_OTU_nTax.pdf", width = 7, height = 5)
ggplot(mp_site_ITS_sig_toPlot, aes(name, n_taxa, fill = taxonomy3, group = taxonomy3)) +
  geom_bar(stat = "identity", colour = NA, linewidth = 0.25) +
  labs(x = "Site",
       y = "Number of indicator OTUs",
       fill = "Class") +
  scale_fill_manual(values = c("grey90", "grey75", brewer.pal(11, "Paired")[11:1])) +
  scale_y_continuous(expand = c(0.01, 0.01)) +  
  theme_classic() +
  theme(axis.title.y = element_text(size = 12),
        axis.title.x = element_blank(),
        axis.text.y = element_text(size = 10),
        axis.text.x = element_text(size = 10, angle = 45, hjust = 1))
dev.off()



#### __Genotype ####
# Indicators of genotype with all 15 sites as input

### 16S
set.seed(425)
mp_Entry_16S <- multipatt(t(input_filt_rare_16S$data_loaded), 
                          input_filt_rare_16S$map_loaded$Entry, 
                          func = "r.g", 
                          duleg = TRUE,
                          control = how(nperm=999))
summary(mp_Entry_16S) # Number of species associated to 1 group: 830 (raw p < 0.05)

# Plot significant and high r
png("InitialFigs/Multipatt_Entry_16S_r20.png", width = 6, height = 8, units = "in", res = 300)
plot_multipatt_asv(mp_obj = mp_Entry_16S, 
                   input = input_filt_rare_16S,
                   tax_sum = tax_sum_asv_16S,
                   group = "Entry",
                   filter = FALSE,
                   #filter_vals = "Field",
                   abund = "% Rel. Abund.",
                   qcut = 1.1,
                   rcut = 0.2)
dev.off()

# Get list of top
mp_site_16S_sig_top <- mp_Entry_16S$sign %>%
  mutate(p.fdr = p.adjust(p.value, "fdr")) %>%
  filter(p.value < 0.05) %>%
  filter(stat > 0.2) %>%
  mutate("ASV_ID" = rownames(.)) %>%
  left_join(., input_filt_rare_16S$taxonomy_loaded, by = c("ASV_ID" = "taxonomy8")) %>%
  set_names(c(substr(names(.)[1:10], 3, nchar(.))), names(.)[11:ncol(.)]) %>%
  arrange(across(names(.)[1:10], desc)) 

# Plot number of classes that are indicators
mp_Entry_16S_sig <- mp_Entry_16S$sign %>%
  mutate(p.fdr = p.adjust(p.value, "fdr")) %>%
  filter(p.value < 0.05) %>%
  mutate("ASV_ID" = rownames(.)) %>%
  left_join(., input_filt_rare_16S$taxonomy_loaded, by = c("ASV_ID" = "taxonomy8")) %>%
  set_names(c(substr(names(.)[1:10], 3, nchar(.))), names(.)[11:ncol(.)])
top_taxa_16S <- as.data.frame(table(mp_Entry_16S_sig$taxonomy2)) %>%
  arrange(desc(Freq)) %>%
  slice_head(n = 12) %>%
  mutate(Var1 = as.character(Var1))
mp_Entry_16S_sig_toPlot <- mp_Entry_16S_sig %>%
  pivot_longer(cols = names(.)[1:10]) %>%
  filter(value > 0) %>%
  group_by(name, taxonomy2) %>%
  summarise(n_taxa = n()) %>%
  drop_na(taxonomy2)
mp_Entry_16S_sig_toPlot$taxonomy2[mp_Entry_16S_sig_toPlot$taxonomy2 %notin% top_taxa_16S$Var1] = 'Other'
mp_Entry_16S_sig_toPlot <- mp_Entry_16S_sig_toPlot %>%
  group_by(name, taxonomy2) %>%
  summarise(n_taxa = sum(n_taxa)) %>%
  mutate(taxonomy2 = factor(taxonomy2,
                            levels = c("Other", 
                                       rev(top_taxa_16S$Var1))))
n_per_geno_16S <- mp_Entry_16S_sig_toPlot %>%
  group_by(name) %>%
  summarise(site_n = sum(n_taxa))
pdf("InitialFigs/Multipatt_Entry_16S_nTax.pdf", width = 7, height = 5)
ggplot(mp_Entry_16S_sig_toPlot, aes(name, n_taxa, fill = taxonomy2, group = taxonomy2)) +
  geom_bar(stat = "identity", colour = NA, linewidth = 0.25) +
  labs(x = "Entry",
       y = "Number of indicator ASVs",
       fill = "Phylum") +
  scale_fill_manual(values = c("grey90", brewer.pal(12, "Paired")[12:1])) +
  scale_y_continuous(expand = c(0.01, 0.01)) +  
  theme_classic() +
  theme(axis.title.y = element_text(size = 12),
        axis.title.x = element_blank(),
        axis.text.y = element_text(size = 10),
        axis.text.x = element_text(size = 10, angle = 45, hjust = 1))
dev.off()

figs8a <- ggplot(mp_Entry_16S_sig_toPlot, aes(name, n_taxa, fill = taxonomy2, group = taxonomy2)) +
  geom_bar(stat = "identity", colour = NA, linewidth = 0.25) +
  labs(x = "Entry",
       y = "# of indicator ASVs",
       fill = "Phylum") +
  scale_fill_manual(values = c("grey90", brewer.pal(12, "Paired")[12:1])) +
  scale_y_continuous(expand = c(0.01, 0.01)) +  
  theme_classic() +
  theme(axis.title.y = element_text(size = 12),
        axis.title.x = element_blank(),
        axis.text.y = element_text(size = 10),
        axis.text.x = element_blank(),
        legend.margin = margin(0, 0, 0, 0, unit = "pt"),
        legend.box.margin = margin(0, 0, 0, 0, unit = "pt"),
        legend.key.size = unit(0.3, "cm"),
        legend.text = element_text(size = 6))
figs8a
legs8a <- get_legend(figs8a)
figs8a <- figs8a + theme(legend.position = "none")



# Redo at OTU level
set.seed(425)
mp_Entry_16S <- multipatt(t(tax_sum_OTU_16S), 
                          input_filt_rare_16S$map_loaded$Entry, 
                          func = "r.g", 
                          duleg = TRUE,
                          control = how(nperm=999))
summary(mp_Entry_16S) # Number of species associated to 1 group: 302 (raw p < 0.05)

# Plot significant and high r
png("InitialFigs/Multipatt_Entry_16S_OTU_r20.png", width = 6, height = 8, units = "in", res = 300)
plot_multipatt_asv(mp_obj = mp_Entry_16S, 
                   input = input_filt_rare_16S_indic,
                   tax_sum = tax_sum_OTU_16S,
                   group = "Entry",
                   filter = FALSE,
                   #filter_vals = "Field",
                   abund = "% Rel. Abund.",
                   qcut = 1.1,
                   rcut = 0.2)
dev.off()

# Get list of top
mp_site_16S_sig_top <- mp_site_16S$sign %>%
  mutate(p.fdr = p.adjust(p.value, "fdr")) %>%
  filter(p.value < 0.05) %>%
  filter(stat > 0.2) %>%
  mutate("ASV_ID" = rownames(.)) %>%
  left_join(., otu_tax_16S, by = c("ASV_ID" = "OTU_ID")) %>%
  set_names(c(substr(names(.)[1:10], 3, nchar(.))), names(.)[11:ncol(.)]) %>%
  arrange(across(names(.)[1:10], desc)) 

# Plot number of classes that are indicators
mp_Entry_16S_sig <- mp_Entry_16S$sign %>%
  mutate(p.fdr = p.adjust(p.value, "fdr")) %>%
  filter(p.value < 0.05) %>%
  mutate("ASV_ID" = rownames(.)) %>%
  left_join(., otu_tax_16S, by = c("ASV_ID" = "OTU_ID")) %>%
  set_names(c(substr(names(.)[1:10], 3, nchar(.))), names(.)[11:ncol(.)])
top_taxa_16S <- as.data.frame(table(mp_Entry_16S_sig$taxonomy2)) %>%
  arrange(desc(Freq)) %>%
  slice_head(n = 12) %>%
  mutate(Var1 = as.character(Var1))
mp_Entry_16S_sig_toPlot <- mp_Entry_16S_sig %>%
  pivot_longer(cols = names(.)[1:10]) %>%
  filter(value > 0) %>%
  group_by(name, taxonomy2) %>%
  summarise(n_taxa = n()) %>%
  drop_na(taxonomy2)
mp_Entry_16S_sig_toPlot$taxonomy2[mp_Entry_16S_sig_toPlot$taxonomy2 %notin% top_taxa_16S$Var1] = 'Other'
mp_Entry_16S_sig_toPlot <- mp_Entry_16S_sig_toPlot %>%
  group_by(name, taxonomy2) %>%
  summarise(n_taxa = sum(n_taxa)) %>%
  mutate(taxonomy2 = factor(taxonomy2,
                            levels = c("Other", 
                                       rev(top_taxa_16S$Var1))))
pdf("InitialFigs/Multipatt_Entry_16S_OTU_nTax.pdf", width = 7, height = 5)
ggplot(mp_Entry_16S_sig_toPlot, aes(name, n_taxa, fill = taxonomy2, group = taxonomy2)) +
  geom_bar(stat = "identity", colour = NA, linewidth = 0.25) +
  labs(x = "Entry",
       y = "Number of indicator OTUs",
       fill = "Phylum") +
  scale_fill_manual(values = c("grey90", brewer.pal(12, "Paired")[12:1])) +
  scale_y_continuous(expand = c(0.01, 0.01)) +  
  theme_classic() +
  theme(axis.title.y = element_text(size = 12),
        axis.title.x = element_blank(),
        axis.text.y = element_text(size = 10),
        axis.text.x = element_text(size = 10, angle = 45, hjust = 1))
dev.off()



### ITS
set.seed(425)
mp_Entry_ITS <- multipatt(t(input_filt_rare_ITS$data_loaded), 
                          input_filt_rare_ITS$map_loaded$Entry, 
                          func = "r.g", 
                          duleg = TRUE,
                          control = how(nperm=999))
summary(mp_Entry_ITS) # Number of species associated to 1 group: 141 (raw p)

# Plot significant and high r
png("InitialFigs/Multipatt_Entry_ITS_r17.png", width = 6, height = 8, units = "in", res = 300)
plot_multipatt_asv(mp_obj = mp_Entry_ITS, 
                   input = input_filt_rare_ITS,
                   tax_sum = tax_sum_asv_ITS,
                   group = "Entry",
                   filter = FALSE,
                   #filter_vals = "Field",
                   abund = "% Rel. Abund.",
                   qcut = 1.1,
                   rcut = 0.17)
dev.off()

# Get list of top
mp_Entry_ITS_sig_top <- mp_Entry_ITS$sign %>%
  mutate(p.fdr = p.adjust(p.value, "fdr")) %>%
  filter(p.value <= 0.05) %>%
  filter(stat > 0.17) %>%
  mutate("ASV_ID" = rownames(.)) %>%
  left_join(., input_filt_rare_16S$taxonomy_loaded, by = c("ASV_ID" = "taxonomy8")) %>%
  set_names(c(substr(names(.)[1:10], 3, nchar(.))), names(.)[11:ncol(.)]) %>%
  arrange(across(names(.)[1:10], desc)) 

# Plot number of classes that are indicators
mp_Entry_ITS_sig <- mp_Entry_ITS$sign %>%
  mutate(p.fdr = p.adjust(p.value, "fdr")) %>%
  filter(p.value <= 0.05) %>%
  mutate("ASV_ID" = rownames(.)) %>%
  left_join(., input_filt_rare_ITS$taxonomy_loaded, by = c("ASV_ID" = "taxonomy8")) %>%
  set_names(c(substr(names(.)[1:10], 3, nchar(.))), names(.)[11:ncol(.)])
top_taxa_ITS <- as.data.frame(table(mp_Entry_ITS_sig$taxonomy3)) %>%
  arrange(desc(Freq)) %>%
  slice_head(n = 12) %>%
  mutate(Var1 = as.character(Var1))
mp_Entry_ITS_sig_toPlot <- mp_Entry_ITS_sig %>%
  pivot_longer(cols = names(.)[1:10]) %>%
  filter(value > 0) %>%
  group_by(name, taxonomy3) %>%
  summarise(n_taxa = n()) %>%
  drop_na(taxonomy3)
mp_Entry_ITS_sig_toPlot$taxonomy3[mp_Entry_ITS_sig_toPlot$taxonomy3 %notin% top_taxa_ITS$Var1] = 'Other'
mp_Entry_ITS_sig_toPlot <- mp_Entry_ITS_sig_toPlot %>%
  group_by(name, taxonomy3) %>%
  summarise(n_taxa = sum(n_taxa)) %>%
  mutate(taxonomy3 = factor(taxonomy3,
                            levels = c("Other", "NA",
                                       rev(subset(top_taxa_ITS$Var1, top_taxa_ITS$Var1 != "NA")))))
n_per_geno_ITS <- mp_Entry_ITS_sig_toPlot %>%
  group_by(name) %>%
  summarise(site_n = sum(n_taxa))
pdf("InitialFigs/Multipatt_Entry_ITS_nTax.pdf", width = 7, height = 5)
ggplot(mp_Entry_ITS_sig_toPlot, aes(name, n_taxa, fill = taxonomy3, group = taxonomy3)) +
  geom_bar(stat = "identity", colour = NA, linewidth = 0.25) +
  labs(x = "Entry",
       y = "Number of indicator ASVs",
       fill = "Class") +
  scale_fill_manual(values = c("grey90", "grey75", brewer.pal(11, "Paired")[11:1])) +
  scale_y_continuous(expand = c(0.01, 0.01)) +  
  theme_classic() +
  theme(axis.title.y = element_text(size = 12),
        axis.title.x = element_blank(),
        axis.text.y = element_text(size = 10),
        axis.text.x = element_text(size = 10, angle = 45, hjust = 1))
dev.off()

figs8b <- ggplot(mp_Entry_ITS_sig_toPlot, aes(name, n_taxa, fill = taxonomy3, group = taxonomy3)) +
  geom_bar(stat = "identity", colour = NA, linewidth = 0.25) +
  labs(x = "Entry",
       y = "# of indicator ASVs",
       fill = "Class") +
  scale_fill_manual(values = c("grey90", "grey75", brewer.pal(11, "Paired")[11:1])) +
  scale_y_continuous(expand = c(0.01, 0.01)) +  
  theme_classic() +
  theme(axis.title.y = element_text(size = 12),
        axis.title.x = element_blank(),
        axis.text.y = element_text(size = 10),
        axis.text.x = element_text(size = 10, angle = 45, hjust = 1),
        legend.margin = margin(0, 0, 0, 0, unit = "pt"),
        legend.box.margin = margin(0, 0, 0, 0, unit = "pt"),
        legend.key.size = unit(0.3, "cm"),
        legend.text = element_text(size = 6))
figs8b
legs8b <- get_legend(figs8b)
figs8b <- figs8b + theme(legend.position = "none")

left <- plot_grid(figs8a, figs8b, ncol = 1, align = "v", rel_heights = c(0.45, 0.55),
                  labels = "auto")
right <- plot_grid(legs8a, legs8b, NULL, ncol = 1, align = "v", rel_heights = c(0.4, 0.5, 0.1))
pdf("FinalFigs/FigureS8.pdf", width = 7, height = 6)
plot_grid(left, right, rel_widths = c(0.8, 0.2))
dev.off()


# Redo at OTU level
set.seed(425)
mp_Entry_ITS <- multipatt(t(tax_sum_OTU_ITS), 
                          input_filt_rare_ITS_indic$map_loaded$Entry, 
                          func = "r.g", 
                          duleg = TRUE,
                          control = how(nperm=999))
summary(mp_Entry_ITS) # Number of species associated to 1 group: 73 (raw p)

# Plot significant and high r
png("InitialFigs/Multipatt_Entry_ITS_OTU_r17.png", width = 6, height = 8, units = "in", res = 300)
plot_multipatt_asv(mp_obj = mp_Entry_ITS, 
                   input = input_filt_rare_ITS_indic,
                   tax_sum = tax_sum_OTU_ITS,
                   group = "Entry",
                   filter = FALSE,
                   #filter_vals = "Field",
                   abund = "% Rel. Abund.",
                   qcut = 1.1,
                   rcut = 0.17)
dev.off()

# Get list of top
mp_Entry_ITS_sig_top <- mp_Entry_ITS$sign %>%
  mutate(p.fdr = p.adjust(p.value, "fdr")) %>%
  filter(p.value < 0.05) %>%
  filter(stat > 0.17) %>%
  mutate("ASV_ID" = rownames(.)) %>%
  left_join(., otu_tax_ITS, by = c("ASV_ID" = "OTU_ID")) %>%
  set_names(c(substr(names(.)[1:10], 3, nchar(.))), names(.)[11:ncol(.)]) %>%
  arrange(across(names(.)[1:10], desc)) 

# Plot number of classes that are indicators
mp_Entry_ITS_sig <- mp_Entry_ITS$sign %>%
  mutate(p.fdr = p.adjust(p.value, "fdr")) %>%
  filter(p.value < 0.05) %>%
  mutate("ASV_ID" = rownames(.)) %>%
  left_join(., otu_tax_ITS, by = c("ASV_ID" = "OTU_ID")) %>%
  set_names(c(substr(names(.)[1:10], 3, nchar(.))), names(.)[11:ncol(.)])
top_taxa_ITS <- as.data.frame(table(mp_Entry_ITS_sig$taxonomy3)) %>%
  arrange(desc(Freq)) %>%
  slice_head(n = 12) %>%
  mutate(Var1 = as.character(Var1))
mp_Entry_ITS_sig_toPlot <- mp_Entry_ITS_sig %>%
  pivot_longer(cols = names(.)[1:10]) %>%
  filter(value > 0) %>%
  group_by(name, taxonomy3) %>%
  summarise(n_taxa = n()) %>%
  drop_na(taxonomy3)
mp_Entry_ITS_sig_toPlot$taxonomy3[mp_Entry_ITS_sig_toPlot$taxonomy3 %notin% top_taxa_ITS$Var1] = 'Other'
mp_Entry_ITS_sig_toPlot <- mp_Entry_ITS_sig_toPlot %>%
  group_by(name, taxonomy3) %>%
  summarise(n_taxa = sum(n_taxa)) %>%
  mutate(taxonomy3 = factor(taxonomy3,
                            levels = c("Other", "NA",
                                       rev(subset(top_taxa_ITS$Var1, top_taxa_ITS$Var1 != "NA")))))
pdf("InitialFigs/Multipatt_Entry_ITS_OTU_nTax.pdf", width = 7, height = 5)
ggplot(mp_Entry_ITS_sig_toPlot, aes(name, n_taxa, fill = taxonomy3, group = taxonomy3)) +
  geom_bar(stat = "identity", colour = NA, linewidth = 0.25) +
  labs(x = "Entry",
       y = "Number of indicator OTUs",
       fill = "Class") +
  scale_fill_manual(values = c("grey90", "grey75", brewer.pal(11, "Paired")[11:1])) +
  scale_y_continuous(expand = c(0.01, 0.01)) +  
  theme_classic() +
  theme(axis.title.y = element_text(size = 12),
        axis.title.x = element_blank(),
        axis.text.y = element_text(size = 10),
        axis.text.x = element_text(size = 10, angle = 45, hjust = 1))
dev.off()



#### _Barplots ####
# Need to make an object of site labels
facet_names <- c("Brookings" = "Brookings",
                 "Burlington" = "Burlington",
                 "Carrington" = "Carrington",
                 "Casselton" = "Casselton",
                 "Grandin" = "Grandin",
                 "Brighton" = "Brighton",
                 "Land Institute" = "Land\nInstitute",
                 "Lindsborg-dryland" = "Lindsborg\ndryland",
                 "Lindsborg-irrigated" = "Lindsborg\nirrigated",
                 "Mandan" = "Mandan",
                 "McLaughlin" = "McLaughlin",
                 "Mentor" = "Mentor",
                 "Pierre" = "Pierre",
                 "Ralls" = "Ralls",
                 "Velva" = "Velva")

#### __16S ####
# For 16S, get functional guild information with FAPROTAX
fap <- read.delim("data/faprotax_output.tsv") %>%
  column_to_rownames(var = "group") %>%
  t() %>%
  as.data.frame() %>%
  mutate(sampleID = rownames(.))
names(fap)

# Quick look
cliffplot_taxa_bars(input = input_filt_rare_16S, level = 1, variable = "Site")
cliffplot_taxa_bars(input = input_filt_rare_16S, level = 2, variable = "Site")
cliffplot_taxa_bars(input = input_filt_rare_16S, level = 3, variable = "Site")
cliffplot_taxa_bars(input = input_filt_rare_16S, level = 4, variable = "Site")
cliffplot_taxa_bars(input = input_filt_rare_16S, level = 5, variable = "Site")
cliffplot_taxa_bars(input = input_filt_rare_16S, level = 6, variable = "Site")
cliffplot_taxa_bars(input = input_filt_rare_16S, level = 7, variable = "Site")

# For each level, summarize, sort by abundance, plot all samples

# Phylum
tax_sum_phyla_16S_587 <- summarize_taxonomy(input = input_filt_rare_16S, 
                                            level = 2, 
                                            report_higher_tax = F)
tax_sum_phyla_16S <- summarize_taxonomy(input = input_filt_rare_16S, 
                                        level = 2, 
                                        report_higher_tax = F)
bars_phyla_16S <- plot_taxa_bars(tax_sum_phyla_16S,
                                 input_filt_rare_16S$map_loaded,
                                 "sampleID",
                                 num_taxa = 12,
                                 data_only = TRUE) %>%
  mutate(taxon = fct_rev(taxon)) %>%
  left_join(., input_filt_rare_16S$map_loaded, by = c("group_by" = "sampleID"))
top_phyla_16S <- bars_phyla_16S %>%
  group_by(taxon) %>%
  summarise(mean = mean(mean_value)) %>%
  filter(taxon != "Other") %>%
  arrange(-mean) %>%
  mutate(taxon = as.character(taxon))
bars_phyla_16S <- bars_phyla_16S %>%
  mutate(taxon = factor(taxon,
                        levels = c("Other", rev(top_phyla_16S$taxon))))
pdf("InitialFigs/Taxa_Phyla_16S.pdf", width = 8, height = 6)
ggplot(bars_phyla_16S, aes(group_by, mean_value, fill = taxon)) +
  geom_bar(stat = "identity", colour = NA, size = 0.25) +
  labs(x = "Sample", y = "Relative abundance", fill = "Phylum") +
  scale_fill_manual(values = c("grey90", brewer.pal(12, "Paired")[12:1])) +
  scale_y_continuous(expand = c(0.01, 0.01)) + 
  facet_grid(~ Site, space = "free", scales = "free_x", 
             labeller = as_labeller(facet_names)) +
  theme_classic() +
  theme(axis.title.y = element_text(face = "bold", size = 12),
        axis.title.x = element_blank(),
        axis.text.y = element_text(size = 10),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.line.x = element_blank(),
        axis.line.y = element_blank(),
        strip.text = element_text(size = 4),
        strip.background = element_rect(linewidth = 0.2),
        legend.position = "right",
        legend.margin = margin(0, 0, 0, 0, unit = "pt"),
        legend.box.margin = margin(0, 0, 0, 0, unit = "pt"),
        legend.key.size = unit(0.3, "cm"),
        panel.spacing.x = unit(c(0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 
                                 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2), 
                               "lines"))
dev.off()

# Guild
tax_sum_guilds_16S <- read.delim("data/faprotax_output.tsv") %>%
  column_to_rownames(var = "group") %>%
  set_names(names(tax_sum_phyla_16S))
bars_guilds_16S <- plot_taxa_bars(tax_sum_guilds_16S,
                                 input_filt_rare_16S$map_loaded,
                                 "sampleID",
                                 num_taxa = 12,
                                 data_only = TRUE) %>%
  mutate(taxon = fct_rev(taxon)) %>%
  left_join(., input_filt_rare_16S$map_loaded, by = c("group_by" = "sampleID"))
top_guilds_16S <- bars_guilds_16S %>%
  group_by(taxon) %>%
  summarise(mean = mean(mean_value)) %>%
  filter(taxon != "Other") %>%
  filter(taxon != "Unassigned") %>%
  arrange(-mean) %>%
  mutate(taxon = as.character(taxon))
bars_guilds_16S <- bars_guilds_16S %>%
  mutate(taxon = factor(taxon,
                        levels = c("Unassigned", "Other", rev(top_guilds_16S$taxon))))
pdf("InitialFigs/Taxa_Guilds_16S.pdf", width = 8, height = 6)
ggplot(bars_guilds_16S, aes(group_by, mean_value, fill = taxon)) +
  geom_bar(stat = "identity", colour = NA, size = 0.25) +
  labs(x = "Sample", y = "Abundance (# of rarefied reads)", fill = "Guild") +
  scale_fill_manual(values = c("grey75", "grey90", brewer.pal(12, "Paired")[12:1])) +
  scale_y_continuous(expand = c(0.01, 0.01)) + 
  facet_grid(~ Site, space = "free", scales = "free_x", 
             labeller = as_labeller(facet_names)) +
  theme_classic() +
  theme(axis.title.y = element_text(face = "bold", size = 12),
        axis.title.x = element_blank(),
        axis.text.y = element_text(size = 10),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.line.x = element_blank(),
        axis.line.y = element_blank(),
        strip.text = element_text(size = 4),
        strip.background = element_rect(linewidth = 0.2),
        legend.position = "right",
        legend.margin = margin(0, 0, 0, 0, unit = "pt"),
        legend.box.margin = margin(0, 0, 0, 0, unit = "pt"),
        legend.key.size = unit(0.3, "cm"),
        panel.spacing.x = unit(c(0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 
                                 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2), 
                               "lines"))
dev.off()

# Class
tax_sum_class_16S <- summarize_taxonomy(input = input_filt_rare_16S, 
                                        level = 3, 
                                        report_higher_tax = F)
bars_class_16S <- plot_taxa_bars(tax_sum_class_16S,
                                 input_filt_rare_16S$map_loaded,
                                 "sampleID",
                                 num_taxa = 12,
                                 data_only = TRUE) %>%
  mutate(taxon = fct_rev(taxon)) %>%
  left_join(., input_filt_rare_16S$map_loaded, by = c("group_by" = "sampleID"))
top_class_16S <- bars_class_16S %>%
  group_by(taxon) %>%
  summarise(mean = mean(mean_value)) %>%
  filter(taxon != "Other") %>%
  arrange(-mean) %>%
  mutate(taxon = as.character(taxon))
bars_class_16S <- bars_class_16S %>%
  mutate(taxon = factor(taxon,
                        levels = c("Other", rev(top_class_16S$taxon))))
pdf("InitialFigs/Taxa_Class_16S.pdf", width = 8, height = 6)
ggplot(bars_class_16S, aes(group_by, mean_value, fill = taxon)) +
  geom_bar(stat = "identity", colour = NA, size = 0.25) +
  labs(x = "Sample", y = "Relative abundance", fill = "Class") +
  scale_fill_manual(values = c("grey90", brewer.pal(12, "Paired")[12:1])) +
  scale_y_continuous(expand = c(0.01, 0.01)) + 
  facet_grid(~ Site, space = "free", scales = "free_x", 
             labeller = as_labeller(facet_names)) +
  theme_classic() +
  theme(axis.title.y = element_text(face = "bold", size = 12),
        axis.title.x = element_blank(),
        axis.text.y = element_text(size = 10),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.line.x = element_blank(),
        axis.line.y = element_blank(),
        strip.text = element_text(size = 4),
        strip.background = element_rect(linewidth = 0.2),
        legend.position = "right",
        legend.margin = margin(0, 0, 0, 0, unit = "pt"),
        legend.box.margin = margin(0, 0, 0, 0, unit = "pt"),
        legend.key.size = unit(0.3, "cm"),
        panel.spacing.x = unit(c(0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 
                                 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2), 
                               "lines"))
dev.off()

# Order
tax_sum_order_16S <- summarize_taxonomy(input = input_filt_rare_16S, 
                                        level = 4, 
                                        report_higher_tax = F)
bars_order_16S <- plot_taxa_bars(tax_sum_order_16S,
                                 input_filt_rare_16S$map_loaded,
                                 "sampleID",
                                 num_taxa = 12,
                                 data_only = TRUE) %>%
  mutate(taxon = fct_rev(taxon)) %>%
  left_join(., input_filt_rare_16S$map_loaded, by = c("group_by" = "sampleID"))
top_order_16S <- bars_order_16S %>%
  group_by(taxon) %>%
  summarise(mean = mean(mean_value)) %>%
  filter(taxon != "Other") %>%
  filter(taxon != "NA") %>%
  arrange(-mean) %>%
  mutate(taxon = as.character(taxon))
bars_order_16S <- bars_order_16S %>%
  mutate(taxon = factor(taxon,
                        levels = c("NA", "Other", rev(top_order_16S$taxon))))
pdf("InitialFigs/Taxa_Order_16S.pdf", width = 8, height = 6)
ggplot(bars_order_16S, aes(group_by, mean_value, fill = taxon)) +
  geom_bar(stat = "identity", colour = NA, size = 0.25) +
  labs(x = "Sample", y = "Relative abundance", fill = "Order") +
  scale_fill_manual(values = c("grey75", "grey90", brewer.pal(11, "Paired")[11:1])) +
  scale_y_continuous(expand = c(0.01, 0.01)) + 
  facet_grid(~ Site, space = "free", scales = "free_x", 
             labeller = as_labeller(facet_names)) +
  theme_classic() +
  theme(axis.title.y = element_text(face = "bold", size = 12),
        axis.title.x = element_blank(),
        axis.text.y = element_text(size = 10),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.line.x = element_blank(),
        axis.line.y = element_blank(),
        strip.text = element_text(size = 4),
        strip.background = element_rect(linewidth = 0.2),
        legend.position = "right",
        legend.margin = margin(0, 0, 0, 0, unit = "pt"),
        legend.box.margin = margin(0, 0, 0, 0, unit = "pt"),
        legend.key.size = unit(0.3, "cm"),
        panel.spacing.x = unit(c(0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 
                                 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2), 
                               "lines"))
dev.off()

# Family
tax_sum_family_16S <- summarize_taxonomy(input = input_filt_rare_16S, 
                                        level = 5, 
                                        report_higher_tax = F)
bars_family_16S <- plot_taxa_bars(tax_sum_family_16S,
                                 input_filt_rare_16S$map_loaded,
                                 "sampleID",
                                 num_taxa = 12,
                                 data_only = TRUE) %>%
  mutate(taxon = fct_rev(taxon)) %>%
  left_join(., input_filt_rare_16S$map_loaded, by = c("group_by" = "sampleID"))
top_family_16S <- bars_family_16S %>%
  group_by(taxon) %>%
  summarise(mean = mean(mean_value)) %>%
  filter(taxon != "Other") %>%
  filter(taxon != "NA") %>%
  arrange(-mean) %>%
  mutate(taxon = as.character(taxon))
bars_family_16S <- bars_family_16S %>%
  mutate(taxon = factor(taxon,
                        levels = c("NA", "Other", rev(top_family_16S$taxon))))
pdf("InitialFigs/Taxa_Family_16S.pdf", width = 8, height = 6)
ggplot(bars_family_16S, aes(group_by, mean_value, fill = taxon)) +
  geom_bar(stat = "identity", colour = NA, size = 0.25) +
  labs(x = "Sample", y = "Relative abundance", fill = "Family") +
  scale_fill_manual(values = c("grey75", "grey90", brewer.pal(11, "Paired")[11:1])) +
  scale_y_continuous(expand = c(0.01, 0.01)) + 
  facet_grid(~ Site, space = "free", scales = "free_x", 
             labeller = as_labeller(facet_names)) +
  theme_classic() +
  theme(axis.title.y = element_text(face = "bold", size = 12),
        axis.title.x = element_blank(),
        axis.text.y = element_text(size = 10),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.line.x = element_blank(),
        axis.line.y = element_blank(),
        strip.text = element_text(size = 4),
        strip.background = element_rect(linewidth = 0.2),
        legend.position = "right",
        legend.margin = margin(0, 0, 0, 0, unit = "pt"),
        legend.box.margin = margin(0, 0, 0, 0, unit = "pt"),
        legend.key.size = unit(0.3, "cm"),
        panel.spacing.x = unit(c(0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 
                                 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2), 
                               "lines"))
dev.off()

# Genus
tax_sum_genus_16S <- summarize_taxonomy(input = input_filt_rare_16S, 
                                        level = 6, 
                                        report_higher_tax = F)
bars_genus_16S <- plot_taxa_bars(tax_sum_genus_16S,
                                 input_filt_rare_16S$map_loaded,
                                 "sampleID",
                                 num_taxa = 12,
                                 data_only = TRUE) %>%
  mutate(taxon = fct_rev(taxon)) %>%
  left_join(., input_filt_rare_16S$map_loaded, by = c("group_by" = "sampleID"))
top_genus_16S <- bars_genus_16S %>%
  group_by(taxon) %>%
  summarise(mean = mean(mean_value)) %>%
  filter(taxon != "Other") %>%
  filter(taxon != "NA") %>%
  arrange(-mean) %>%
  mutate(taxon = as.character(taxon))
bars_genus_16S <- bars_genus_16S %>%
  mutate(taxon = factor(taxon,
                        levels = c("NA", "Other", rev(top_genus_16S$taxon))))
pdf("InitialFigs/Taxa_Genus_16S.pdf", width = 8, height = 6)
ggplot(bars_genus_16S, aes(group_by, mean_value, fill = taxon)) +
  geom_bar(stat = "identity", colour = NA, size = 0.25) +
  labs(x = "Sample", y = "Relative abundance", fill = "Genus") +
  scale_fill_manual(values = c("grey75", "grey90", brewer.pal(11, "Paired")[11:1])) +
  scale_y_continuous(expand = c(0.01, 0.01)) + 
  facet_grid(~ Site, space = "free", scales = "free_x", 
             labeller = as_labeller(facet_names)) +
  theme_classic() +
  theme(axis.title.y = element_text(face = "bold", size = 12),
        axis.title.x = element_blank(),
        axis.text.y = element_text(size = 10),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.line.x = element_blank(),
        axis.line.y = element_blank(),
        strip.text = element_text(size = 4),
        strip.background = element_rect(linewidth = 0.2),
        legend.position = "right",
        legend.margin = margin(0, 0, 0, 0, unit = "pt"),
        legend.box.margin = margin(0, 0, 0, 0, unit = "pt"),
        legend.key.size = unit(0.3, "cm"),
        panel.spacing.x = unit(c(0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 
                                 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2), 
                               "lines"))
dev.off()

# Need prevalence of certain genera for discussion
gen_int_16S <- c("Pseudomonas", "Bacillus", "Paenibacillus", "Variovorax", "Sphingobacterium")
gen_int_16S_df <- tax_sum_genus_16S %>%
  filter(rownames(.) %in% gen_int_16S) %>%
  t() %>%
  as.data.frame() %>%
  mutate(Site = input_filt_rare_16S$map_loaded$Site)
prev_5_all <- gen_int_16S_df %>%
  pivot_longer(cols = c(1:5), names_to = "Genus") %>%
  group_by(Genus) %>%
  summarise(Genus_n = sum(value > 0)) %>%
  mutate(Gen_per = Genus_n/587*100)
prev_5 <- gen_int_16S_df %>%
  pivot_longer(cols = c(1:5), names_to = "Genus") %>%
  group_by(Site, Genus) %>%
  summarise(site_n = n(),
            Genus_n = sum(value > 0)) %>%
  mutate(Gen_per = Genus_n/site_n*100)
prev_5_site <- prev_5 %>%
  group_by(Site) %>%
  summarise(n_genus = sum(Genus_n > 0)) %>%
  mutate(per_genus = n_genus/5*100)

# OTU
tax_sum_OTU_16S <- summarize_taxonomy(input = input_filt_rare_16S, 
                                      level = 9, 
                                      report_higher_tax = F)
bars_OTU_16S <- plot_taxa_bars(tax_sum_OTU_16S,
                               input_filt_rare_16S$map_loaded,
                               "sampleID",
                               num_taxa = 12,
                               data_only = TRUE) %>%
  mutate(taxon = fct_rev(taxon)) %>%
  left_join(., input_filt_rare_16S$map_loaded, by = c("group_by" = "sampleID"))
top_OTU_16S <- bars_OTU_16S %>%
  group_by(taxon) %>%
  summarise(mean = mean(mean_value)) %>%
  filter(taxon != "Other") %>%
  arrange(-mean) %>%
  mutate(taxon = as.character(taxon))
bars_OTU_16S <- bars_OTU_16S %>%
  mutate(taxon = factor(taxon,
                        levels = c("Other", rev(top_OTU_16S$taxon))))
pdf("InitialFigs/Taxa_OTU_16S.pdf", width = 8, height = 6)
ggplot(bars_OTU_16S, aes(group_by, mean_value, fill = taxon)) +
  geom_bar(stat = "identity", colour = NA, linewidth = 0.25) +
  labs(x = "Sample", y = "Relative abundance", fill = "OTU") +
  scale_fill_manual(values = c("grey90", brewer.pal(12, "Paired")[12:1])) +
  scale_y_continuous(expand = c(0.01, 0.01)) + 
  facet_grid(~ Site, space = "free", scales = "free_x", 
             labeller = as_labeller(facet_names)) +
  theme_classic() +
  theme(axis.title.y = element_text(face = "bold", size = 12),
        axis.title.x = element_blank(),
        axis.text.y = element_text(size = 10),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.line.x = element_blank(),
        axis.line.y = element_blank(),
        strip.text = element_text(size = 4),
        strip.background = element_rect(linewidth = 0.2),
        legend.position = "right",
        legend.margin = margin(0, 0, 0, 0, unit = "pt"),
        legend.box.margin = margin(0, 0, 0, 0, unit = "pt"),
        legend.key.size = unit(0.3, "cm"),
        panel.spacing.x = unit(c(0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 
                                 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2), 
                               "lines"))
dev.off()

# ASV (for this one, do rank abundance)
bars_asv_16S <- plot_taxa_bars(tax_sum_asv_16S,
                               input_filt_rare_16S$map_loaded,
                               "sampleID",
                               num_taxa = 12,
                               data_only = TRUE) %>%
  mutate(taxon = fct_rev(taxon)) %>%
  left_join(., input_filt_rare_16S$map_loaded, by = c("group_by" = "sampleID"))
top_asv_16S <- bars_asv_16S %>%
  group_by(taxon) %>%
  summarise(mean = mean(mean_value)) %>%
  filter(taxon != "Other") %>%
  filter(taxon != "NA") %>%
  arrange(-mean) %>%
  mutate(taxon = as.character(taxon)) %>%
  left_join(., input_filt_rare_16S$taxonomy_loaded, by = c("taxon" = "taxonomy7")) %>%
  mutate(taxonomy = paste(taxonomy5, taxonomy6, sep = "; ")) %>%
  mutate(taxonomy = ifelse(taxonomy == "NA; NA",
                           taxonomy4,
                           taxonomy))
top_asv_16S_taxa <- top_asv_16S %>%
  mutate(mean = ifelse(mean > 0.004,
                       mean - mean/1.2,
                       mean))
pdf("InitialFigs/Taxa_ASV_16S.pdf", width = 5, height = 7)
ggplot(top_asv_16S, aes(reorder(taxon, mean, median), mean, fill = taxonomy6)) +
  geom_bar(stat = "identity", alpha = 0.75) +
  geom_text(data = top_asv_16S_taxa,
            aes(taxon, mean, label = taxonomy),
            size = 3, angle = 90, inherit.aes = F, hjust = 0) +
  scale_x_discrete(limits = rev) +
  scale_fill_brewer(palette = "Paired") +
  labs(x = "ASV_ID",
       y = "Mean relative abundance",
       fill = "Genus") +
  theme_bw() +
  theme(legend.position = "none",
        axis.title = element_text(size = 12),
        axis.text.y = element_text(size = 10),
        axis.text.x = element_text(size = 10, angle = 45, hjust = 1))
dev.off()

# Prevalence
ncol(tax_sum_asv_16S)
prev_16S <- data.frame("ASV_ID" = rownames(tax_sum_asv_16S),
                       "Absent" = rowSums(tax_sum_asv_16S==0)) %>%
  mutate(Present = ncol(tax_sum_asv_16S) - Absent) %>%
  mutate(Present_Perc = Present/ncol(tax_sum_asv_16S)*100)
top_prev_16S <- prev_16S %>%
  filter(Present_Perc > 80) %>%
  left_join(., input_filt_rare_16S$taxonomy_loaded, by = c("ASV_ID" = "taxonomy7")) %>%
  mutate(taxonomy = paste(taxonomy2, taxonomy3, taxonomy4, taxonomy5, taxonomy6, 
                          sep = "; "))
pdf("InitialFigs/Taxa_ASV_prev_16S.pdf", width = 5, height = 7)
ggplot(top_prev_16S, aes(reorder(ASV_ID, Present_Perc, median), Present_Perc, fill = taxonomy3)) +
  geom_bar(stat = "identity", alpha = 0.75) +
  geom_text(data = top_prev_16S,
            aes(ASV_ID, Present_Perc-80, label = taxonomy),
            size = 2, angle = 90, inherit.aes = F, hjust = 0) +
  scale_x_discrete(limits = rev) +
  scale_fill_brewer(palette = "Paired") +
  labs(x = "ASV_ID",
       y = "% Prevalence",
       fill = "Family, Genus, Species") +
  theme_bw() +
  theme(legend.position = "none",
        axis.title = element_text(size = 12),
        axis.text.y = element_text(size = 10),
        axis.text.x = element_text(size = 10, angle = 45, hjust = 1))
dev.off()



#### __ITS ####
# For ITS, get functional guild information from FUNGuild
# Load
fung_guilds <- readRDS("data/fung_guilds.rds")

# How many ASVs assigned guilds?
sum(is.na(fung_guilds$guild)) # 1100 unassigned
sum(is.na(fung_guilds$guild) == FALSE) # 4194 assigned

# One tricky thins is that confidence in the assignment varies
table(fung_guilds$confidenceRanking)
# 522 Highly Probable
# 1133 Possible
# 2539 Probable

# Another tricky thing is that assignments often include multiple guilds
guildcounts <- as.data.frame(table(fung_guilds$guild))

# Need to Set NA to "Unassigned" multiple guild to "Multiple" and single guilds to that guild.
# Also curate a new column to merge redundant guilds even if they are multiple
# E.g., dung saprotroph and wood saprotroph can be just "Saprotroph"
# Also, the notes column specifies if it is a dark septate endophyte so pull that info too
fung_guilds <- fung_guilds %>%
  mutate(guild2 = guild) %>%
  replace_na(list(guild2 = "Unassigned")) %>%
  mutate(guild2 = gsub("\\|", "", guild2)) %>%
  mutate(guild2 = ifelse(grepl("septate", notes, ignore.case = TRUE) == TRUE,
                         "Dark Septate Endophyte",
                         guild2)) %>%
  mutate(guild3 = ifelse(grepl("-", guild2) == TRUE,
                         "Multiple",
                         guild2)) %>%
  mutate(guild2 = as.factor(guild2)) %>%
  mutate(guild4 = recode_factor(guild2,
                                `Animal Parasite-Fungal Parasite` = "Parasite",
                                `Dung Saprotroph` = "Saprotroph",
                                `Dung Saprotroph-Plant Saprotroph` = "Saprotroph",
                                `Dung Saprotroph-Plant Saprotroph-Soil Saprotroph` = "Saprotroph",
                                `Dung Saprotroph-Plant Saprotroph-Wood Saprotroph` = "Saprotroph",
                                `Dung Saprotroph-Soil Saprotroph` = "Saprotroph",
                                `Dung Saprotroph-Soil Saprotroph-Wood Saprotroph` = "Saprotroph",
                                `Dung Saprotroph-Undefined Saprotroph` = "Saprotroph",
                                `Fungal Parasite` = "Parasite",
                                `Lichen Parasite` = "Parasite",
                                `Litter Saprotroph` = "Saprotroph",
                                `Litter Saprotroph-Soil Saprotroph-Wood Saprotroph` = "Saprotroph",
                                `Plant Saprotroph` = "Saprotroph",
                                `Plant Saprotroph-Undefined Saprotroph` = "Saprotroph",
                                `Plant Saprotroph-Undefined Saprotroph-Wood Saprotroph` = "Saprotroph",
                                `Plant Saprotroph-Wood Saprotroph` = "Saprotroph",
                                `Pollen Saprotroph` = "Saprotroph",
                                `Soil Saprotroph-Undefined Saprotroph`= "Saprotroph",
                                `Undefined Saprotroph` = "Saprotroph",
                                `Undefined Saprotroph-Wood Saprotroph` = "Saprotroph",
                                `Wood Saprotroph` = "Saprotroph")) %>%
  mutate(guild4 = as.character(guild4)) %>%
  mutate(guild_curated = ifelse(grepl("-", guild4) == TRUE,
                                "Multiple",
                                guild4))
guildcounts <- as.data.frame(table(fung_guilds$guild2))
formcounts <- as.data.frame(table(fung_guilds$growthForm))
guildcounts3 <- as.data.frame(table(fung_guilds$guild3))
guildcounts_curated <- as.data.frame(table(fung_guilds$guild_curated))

# Add guild_curated to taxonomy table
# (First pass was guild3)
sum(fung_guilds$ASV_ID != input_filt_rare_ITS$taxonomy_loaded$taxonomy8)
input_filt_rare_ITS$taxonomy_loaded$Guild <- fung_guilds$guild_curated

# Quick look
cliffplot_taxa_bars(input = input_filt_rare_ITS, level = 1, variable = "Site")
cliffplot_taxa_bars(input = input_filt_rare_ITS, level = 2, variable = "Site")
cliffplot_taxa_bars(input = input_filt_rare_ITS, level = 3, variable = "Site")
cliffplot_taxa_bars(input = input_filt_rare_ITS, level = 4, variable = "Site")
cliffplot_taxa_bars(input = input_filt_rare_ITS, level = 5, variable = "Site")
cliffplot_taxa_bars(input = input_filt_rare_ITS, level = 6, variable = "Site")
cliffplot_taxa_bars(input = input_filt_rare_ITS, level = 9, variable = "Site")

# For each level, summarize, sort by abundance, plot all samples

# Phylum
tax_sum_phyla_ITS <- summarize_taxonomy(input = input_filt_rare_ITS, 
                                        level = 2, 
                                        report_higher_tax = F)
bars_phyla_ITS <- plot_taxa_bars(tax_sum_phyla_ITS,
                                 input_filt_rare_ITS$map_loaded,
                                 "sampleID",
                                 num_taxa = 12,
                                 data_only = TRUE) %>%
  mutate(taxon = fct_rev(taxon)) %>%
  left_join(., input_filt_rare_ITS$map_loaded, by = c("group_by" = "sampleID"))
top_phyla_ITS <- bars_phyla_ITS %>%
  group_by(taxon) %>%
  summarise(mean = mean(mean_value)) %>%
  filter(taxon != "Other") %>%
  arrange(-mean) %>%
  mutate(taxon = as.character(taxon))
bars_phyla_ITS <- bars_phyla_ITS %>%
  mutate(taxon = factor(taxon,
                        levels = c("Other", rev(top_phyla_ITS$taxon))))
pdf("InitialFigs/Taxa_Phyla_ITS.pdf", width = 8, height = 6)
ggplot(bars_phyla_ITS, aes(group_by, mean_value, fill = taxon)) +
  geom_bar(stat = "identity", colour = NA, size = 0.25) +
  labs(x = "Sample", y = "Relative abundance", fill = "Phylum") +
  scale_fill_manual(values = c("grey90", brewer.pal(12, "Paired")[12:1])) +
  scale_y_continuous(expand = c(0.01, 0.01)) + 
  facet_grid(~ Site, space = "free", scales = "free_x", 
             labeller = as_labeller(facet_names)) +
  theme_classic() +
  theme(axis.title.y = element_text(face = "bold", size = 12),
        axis.title.x = element_blank(),
        axis.text.y = element_text(size = 10),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.line.x = element_blank(),
        axis.line.y = element_blank(),
        strip.text = element_text(size = 4),
        strip.background = element_rect(linewidth = 0.2),
        legend.position = "right",
        legend.margin = margin(0, 0, 0, 0, unit = "pt"),
        legend.box.margin = margin(0, 0, 0, 0, unit = "pt"),
        legend.key.size = unit(0.3, "cm"),
        panel.spacing.x = unit(c(0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 
                                 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2), 
                               "lines"))
dev.off()

# Guild
tax_sum_guilds_ITS <- summarize_taxonomy(input = input_filt_rare_ITS, 
                                         level = 10, 
                                         report_higher_tax = F)
bars_guilds_ITS <- plot_taxa_bars(tax_sum_guilds_ITS,
                                  input_filt_rare_ITS$map_loaded,
                                  "sampleID",
                                  num_taxa = 13, # There are only 13 curated guilds anyways
                                  data_only = TRUE) %>%
  mutate(taxon = fct_rev(taxon)) %>%
  left_join(., input_filt_rare_ITS$map_loaded, by = c("group_by" = "sampleID"))
top_guilds_ITS <- bars_guilds_ITS %>%
  group_by(taxon) %>%
  summarise(mean = mean(mean_value)) %>%
  #filter(taxon != "Other") %>%
  filter(taxon != "Unassigned") %>%
  filter(taxon != "Multiple") %>%
  arrange(-mean) %>%
  mutate(taxon = as.character(taxon))
bars_guilds_ITS <- bars_guilds_ITS %>%
  mutate(taxon = factor(taxon,
                        levels = c("Unassigned", "Multiple", rev(top_guilds_ITS$taxon))))
pdf("InitialFigs/Taxa_Guilds_ITS_Curated.pdf", width = 8, height = 6)
ggplot(bars_guilds_ITS, aes(group_by, mean_value, fill = taxon)) +
  geom_bar(stat = "identity", colour = NA, size = 0.25) +
  labs(x = "Sample", y = "Relative abundance", fill = "Guild") +
  scale_fill_manual(values = c("grey75", "grey90", brewer.pal(11, "Paired")[11:1])) +
  scale_y_continuous(expand = c(0.01, 0.01)) + 
  facet_grid(~ Site, space = "free", scales = "free_x", 
             labeller = as_labeller(facet_names)) +
  theme_classic() +
  theme(axis.title.y = element_text(face = "bold", size = 12),
        axis.title.x = element_blank(),
        axis.text.y = element_text(size = 10),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.line.x = element_blank(),
        axis.line.y = element_blank(),
        strip.text = element_text(size = 4),
        strip.background = element_rect(linewidth = 0.2),
        legend.position = "right",
        legend.margin = margin(0, 0, 0, 0, unit = "pt"),
        legend.box.margin = margin(0, 0, 0, 0, unit = "pt"),
        legend.key.size = unit(0.3, "cm"),
        panel.spacing.x = unit(c(0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 
                                 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2), 
                               "lines"))
dev.off()

# Class
tax_sum_class_ITS <- summarize_taxonomy(input = input_filt_rare_ITS, 
                                        level = 3, 
                                        report_higher_tax = F)
bars_class_ITS <- plot_taxa_bars(tax_sum_class_ITS,
                                 input_filt_rare_ITS$map_loaded,
                                 "sampleID",
                                 num_taxa = 12,
                                 data_only = TRUE) %>%
  mutate(taxon = fct_rev(taxon)) %>%
  left_join(., input_filt_rare_ITS$map_loaded, by = c("group_by" = "sampleID"))
top_class_ITS <- bars_class_ITS %>%
  group_by(taxon) %>%
  summarise(mean = mean(mean_value)) %>%
  filter(taxon != "Other") %>%
  filter(taxon != "NA") %>%
  arrange(-mean) %>%
  mutate(taxon = as.character(taxon))
bars_class_ITS <- bars_class_ITS %>%
  mutate(taxon = factor(taxon,
                        levels = c("NA", "Other", rev(top_class_ITS$taxon))))
pdf("InitialFigs/Taxa_Class_ITS.pdf", width = 8, height = 6)
ggplot(bars_class_ITS, aes(group_by, mean_value, fill = taxon)) +
  geom_bar(stat = "identity", colour = NA, size = 0.25) +
  labs(x = "Sample", y = "Relative abundance", fill = "Class") +
  scale_fill_manual(values = c("grey75", "grey90", brewer.pal(11, "Paired")[11:1])) +
  scale_y_continuous(expand = c(0.01, 0.01)) + 
  facet_grid(~ Site, space = "free", scales = "free_x", 
             labeller = as_labeller(facet_names)) +
  theme_classic() +
  theme(axis.title.y = element_text(face = "bold", size = 12),
        axis.title.x = element_blank(),
        axis.text.y = element_text(size = 10),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.line.x = element_blank(),
        axis.line.y = element_blank(),
        strip.text = element_text(size = 4),
        strip.background = element_rect(linewidth = 0.2),
        legend.position = "right",
        legend.margin = margin(0, 0, 0, 0, unit = "pt"),
        legend.box.margin = margin(0, 0, 0, 0, unit = "pt"),
        legend.key.size = unit(0.3, "cm"),
        panel.spacing.x = unit(c(0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 
                                 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2), 
                               "lines"))
dev.off()

# Order
tax_sum_order_ITS <- summarize_taxonomy(input = input_filt_rare_ITS, 
                                        level = 4, 
                                        report_higher_tax = F)
bars_order_ITS <- plot_taxa_bars(tax_sum_order_ITS,
                                 input_filt_rare_ITS$map_loaded,
                                 "sampleID",
                                 num_taxa = 12,
                                 data_only = TRUE) %>%
  mutate(taxon = fct_rev(taxon)) %>%
  left_join(., input_filt_rare_ITS$map_loaded, by = c("group_by" = "sampleID"))
top_order_ITS <- bars_order_ITS %>%
  group_by(taxon) %>%
  summarise(mean = mean(mean_value)) %>%
  filter(taxon != "Other") %>%
  filter(taxon != "NA") %>%
  arrange(-mean) %>%
  mutate(taxon = as.character(taxon))
bars_order_ITS <- bars_order_ITS %>%
  mutate(taxon = factor(taxon,
                        levels = c("NA", "Other", rev(top_order_ITS$taxon))))
pdf("InitialFigs/Taxa_Order_ITS.pdf", width = 8, height = 6)
ggplot(bars_order_ITS, aes(group_by, mean_value, fill = taxon)) +
  geom_bar(stat = "identity", colour = NA, size = 0.25) +
  labs(x = "Sample", y = "Relative abundance", fill = "Order") +
  scale_fill_manual(values = c("grey75", "grey90", brewer.pal(11, "Paired")[11:1])) +
  scale_y_continuous(expand = c(0.01, 0.01)) + 
  facet_grid(~ Site, space = "free", scales = "free_x", 
             labeller = as_labeller(facet_names)) +
  theme_classic() +
  theme(axis.title.y = element_text(face = "bold", size = 12),
        axis.title.x = element_blank(),
        axis.text.y = element_text(size = 10),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.line.x = element_blank(),
        axis.line.y = element_blank(),
        strip.text = element_text(size = 4),
        strip.background = element_rect(linewidth = 0.2),
        legend.position = "right",
        legend.margin = margin(0, 0, 0, 0, unit = "pt"),
        legend.box.margin = margin(0, 0, 0, 0, unit = "pt"),
        legend.key.size = unit(0.3, "cm"),
        panel.spacing.x = unit(c(0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 
                                 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2), 
                               "lines"))
dev.off()

# Family
tax_sum_family_ITS <- summarize_taxonomy(input = input_filt_rare_ITS, 
                                         level = 5, 
                                         report_higher_tax = F)
bars_family_ITS <- plot_taxa_bars(tax_sum_family_ITS,
                                  input_filt_rare_ITS$map_loaded,
                                  "sampleID",
                                  num_taxa = 12,
                                  data_only = TRUE) %>%
  mutate(taxon = fct_rev(taxon)) %>%
  left_join(., input_filt_rare_ITS$map_loaded, by = c("group_by" = "sampleID"))
top_family_ITS <- bars_family_ITS %>%
  group_by(taxon) %>%
  summarise(mean = mean(mean_value)) %>%
  filter(taxon != "Other") %>%
  arrange(-mean) %>%
  mutate(taxon = as.character(taxon))
bars_family_ITS <- bars_family_ITS %>%
  mutate(taxon = factor(taxon,
                        levels = c("Other", rev(top_family_ITS$taxon))))
pdf("InitialFigs/Taxa_Family_ITS.pdf", width = 8, height = 6)
ggplot(bars_family_ITS, aes(group_by, mean_value, fill = taxon)) +
  geom_bar(stat = "identity", colour = NA, size = 0.25) +
  labs(x = "Sample", y = "Relative abundance", fill = "Family") +
  scale_fill_manual(values = c("grey90", brewer.pal(12, "Paired")[12:1])) +
  scale_y_continuous(expand = c(0.01, 0.01)) + 
  facet_grid(~ Site, space = "free", scales = "free_x", 
             labeller = as_labeller(facet_names)) +
  theme_classic() +
  theme(axis.title.y = element_text(face = "bold", size = 12),
        axis.title.x = element_blank(),
        axis.text.y = element_text(size = 10),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.line.x = element_blank(),
        axis.line.y = element_blank(),
        strip.text = element_text(size = 4),
        strip.background = element_rect(linewidth = 0.2),
        legend.position = "right",
        legend.margin = margin(0, 0, 0, 0, unit = "pt"),
        legend.box.margin = margin(0, 0, 0, 0, unit = "pt"),
        legend.key.size = unit(0.3, "cm"),
        panel.spacing.x = unit(c(0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 
                                 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2), 
                               "lines"))
dev.off()

# Genus
tax_sum_genus_ITS <- summarize_taxonomy(input = input_filt_rare_ITS, 
                                        level = 6, 
                                        report_higher_tax = F)
bars_genus_ITS <- plot_taxa_bars(tax_sum_genus_ITS,
                                 input_filt_rare_ITS$map_loaded,
                                 "sampleID",
                                 num_taxa = 12,
                                 data_only = TRUE) %>%
  mutate(taxon = fct_rev(taxon)) %>%
  left_join(., input_filt_rare_ITS$map_loaded, by = c("group_by" = "sampleID"))
top_genus_ITS <- bars_genus_ITS %>%
  group_by(taxon) %>%
  summarise(mean = mean(mean_value)) %>%
  filter(taxon != "Other") %>%
  filter(taxon != "NA") %>%
  arrange(-mean) %>%
  mutate(taxon = as.character(taxon))
bars_genus_ITS <- bars_genus_ITS %>%
  mutate(taxon = factor(taxon,
                        levels = c("NA", "Other", rev(top_genus_ITS$taxon))))
pdf("InitialFigs/Taxa_Genus_ITS.pdf", width = 8, height = 6)
ggplot(bars_genus_ITS, aes(group_by, mean_value, fill = taxon)) +
  geom_bar(stat = "identity", colour = NA, size = 0.25) +
  labs(x = "Sample", y = "Relative abundance", fill = "Genus") +
  scale_fill_manual(values = c("grey75", "grey90", brewer.pal(11, "Paired")[11:1])) +
  scale_y_continuous(expand = c(0.01, 0.01)) + 
  facet_grid(~ Site, space = "free", scales = "free_x", 
             labeller = as_labeller(facet_names)) +
  theme_classic() +
  theme(axis.title.y = element_text(face = "bold", size = 12),
        axis.title.x = element_blank(),
        axis.text.y = element_text(size = 10),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.line.x = element_blank(),
        axis.line.y = element_blank(),
        strip.text = element_text(size = 4),
        strip.background = element_rect(linewidth = 0.2),
        legend.position = "right",
        legend.margin = margin(0, 0, 0, 0, unit = "pt"),
        legend.box.margin = margin(0, 0, 0, 0, unit = "pt"),
        legend.key.size = unit(0.3, "cm"),
        panel.spacing.x = unit(c(0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 
                                 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2), 
                               "lines"))
dev.off()

# Need prevalence of certain genera for discussion
gen_int_ITS <- c("Mortierella", "Podila", "Gloeoporus")
gen_int_ITS_df <- tax_sum_genus_ITS %>%
  filter(rownames(.) %in% gen_int_ITS) %>%
  t() %>%
  as.data.frame() %>%
  mutate(Site = input_filt_rare_ITS$map_loaded$Site)
prev_3_all <- gen_int_ITS_df %>%
  pivot_longer(cols = c(1:3), names_to = "Genus") %>%
  group_by(Genus) %>%
  summarise(Genus_n = sum(value > 0)) %>%
  mutate(Gen_per = Genus_n/586*100)
prev_3 <- gen_int_ITS_df %>%
  pivot_longer(cols = c(1:3), names_to = "Genus") %>%
  group_by(Site, Genus) %>%
  summarise(site_n = n(),
            Genus_n = sum(value > 0)) %>%
  mutate(Gen_per = Genus_n/site_n*100)
prev_3_site <- prev_3 %>%
  group_by(Site) %>%
  summarise(n_genus = sum(Genus_n > 0)) %>%
  mutate(per_genus = n_genus/3*100)

# OTU
tax_sum_OTU_ITS <- summarize_taxonomy(input = input_filt_rare_ITS, 
                                      level = 9, 
                                      report_higher_tax = F)
bars_OTU_ITS <- plot_taxa_bars(tax_sum_OTU_ITS,
                               input_filt_rare_ITS$map_loaded,
                               "sampleID",
                               num_taxa = 12,
                               data_only = TRUE) %>%
  mutate(taxon = fct_rev(taxon)) %>%
  left_join(., input_filt_rare_ITS$map_loaded, by = c("group_by" = "sampleID"))
top_OTU_ITS <- bars_OTU_ITS %>%
  group_by(taxon) %>%
  summarise(mean = mean(mean_value)) %>%
  filter(taxon != "Other") %>%
  arrange(-mean) %>%
  mutate(taxon = as.character(taxon))
bars_OTU_ITS <- bars_OTU_ITS %>%
  mutate(taxon = factor(taxon,
                        levels = c("Other", rev(top_OTU_ITS$taxon))))
pdf("InitialFigs/Taxa_OTU_ITS.pdf", width = 8, height = 6)
ggplot(bars_OTU_ITS, aes(group_by, mean_value, fill = taxon)) +
  geom_bar(stat = "identity", colour = NA, linewidth = 0.25) +
  labs(x = "Sample", y = "Relative abundance", fill = "OTU") +
  scale_fill_manual(values = c("grey90", brewer.pal(12, "Paired")[12:1])) +
  scale_y_continuous(expand = c(0.01, 0.01)) + 
  facet_grid(~ Site, space = "free", scales = "free_x", 
             labeller = as_labeller(facet_names)) +
  theme_classic() +
  theme(axis.title.y = element_text(face = "bold", size = 12),
        axis.title.x = element_blank(),
        axis.text.y = element_text(size = 10),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.line.x = element_blank(),
        axis.line.y = element_blank(),
        strip.text = element_text(size = 4),
        strip.background = element_rect(linewidth = 0.2),
        legend.position = "right",
        legend.margin = margin(0, 0, 0, 0, unit = "pt"),
        legend.box.margin = margin(0, 0, 0, 0, unit = "pt"),
        legend.key.size = unit(0.3, "cm"),
        panel.spacing.x = unit(c(0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 
                                 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2), 
                               "lines"))
dev.off()

# ASV (for this one, do rank abundance)
bars_asv_ITS <- plot_taxa_bars(tax_sum_asv_ITS,
                               input_filt_rare_ITS$map_loaded,
                               "sampleID",
                               num_taxa = 12,
                               data_only = TRUE) %>%
  mutate(taxon = fct_rev(taxon)) %>%
  left_join(., input_filt_rare_ITS$map_loaded, by = c("group_by" = "sampleID"))
top_asv_ITS <- bars_asv_ITS %>%
  group_by(taxon) %>%
  summarise(mean = mean(mean_value)) %>%
  filter(taxon != "Other") %>%
  filter(taxon != "NA") %>%
  arrange(-mean) %>%
  mutate(taxon = as.character(taxon)) %>%
  left_join(., input_filt_rare_ITS$taxonomy_loaded, by = c("taxon" = "taxonomy8")) %>%
  mutate(taxonomy = paste(taxonomy5, taxonomy6, taxonomy7, sep = "; "))
top_asv_ITS_taxa <- top_asv_ITS %>%
  mutate(mean = ifelse(mean > 0.05,
                       mean - mean/1.5,
                       mean))
pdf("InitialFigs/Taxa_ASV_ITS.pdf", width = 5, height = 7)
ggplot(top_asv_ITS, aes(reorder(taxon, mean, median), mean, fill = taxonomy)) +
  geom_bar(stat = "identity", alpha = 0.75) +
  geom_text(data = top_asv_ITS_taxa, 
            aes(taxon, mean, label = taxonomy), 
            size = 3, angle = 90, inherit.aes = F, hjust = 0) +
  scale_x_discrete(limits = rev) +
  scale_fill_brewer(palette = "Paired") +
  labs(x = "ASV_ID",
       y = "Mean relative abundance",
       fill = "Family, Genus, Species") +
  theme_bw() +
  theme(legend.position = "none",
        axis.title = element_text(size = 12),
        axis.text.y = element_text(size = 10),
        axis.text.x = element_text(size = 10, angle = 45, hjust = 1))
dev.off()

# Prevalence
ncol(tax_sum_asv_ITS)
prev_ITS <- data.frame("ASV_ID" = rownames(tax_sum_asv_ITS),
                       "Absent" = rowSums(tax_sum_asv_ITS==0)) %>%
  mutate(Present = ncol(tax_sum_asv_ITS) - Absent) %>%
  mutate(Present_Perc = Present/ncol(tax_sum_asv_ITS)*100)
top_prev_ITS <- prev_ITS %>%
  filter(Present_Perc > 70) %>%
  left_join(., input_filt_rare_ITS$taxonomy_loaded, by = c("ASV_ID" = "taxonomy8")) %>%
  mutate(taxonomy = paste(taxonomy5, taxonomy6, taxonomy7, sep = "; ")) %>%
  mutate(taxonomy = ifelse(taxonomy == "NA; NA; NA",
                           taxonomy3,
                           taxonomy))
pdf("InitialFigs/Taxa_ASV_prev_ITS.pdf", width = 5, height = 7)
ggplot(top_prev_ITS, aes(reorder(ASV_ID, Present_Perc, median), Present_Perc, fill = taxonomy)) +
  geom_bar(stat = "identity", alpha = 0.75) +
  geom_text(data = top_prev_ITS,
            aes(ASV_ID, Present_Perc-50, label = taxonomy),
            size = 3, angle = 90, inherit.aes = F, hjust = 0) +
  scale_x_discrete(limits = rev) +
  scale_fill_brewer(palette = "Paired") +
  labs(x = "ASV_ID",
       y = "% Prevalence",
       fill = "Family, Genus, Species") +
  theme_bw() +
  theme(legend.position = "none",
        axis.title = element_text(size = 12),
        axis.text.y = element_text(size = 10),
        axis.text.x = element_text(size = 10, angle = 45, hjust = 1))
dev.off()



#### __Combined ####
# Need to combine 16S and ITS for paper
# Phyla (Figure 4), Genera (Figure S4), and Guilds (Figure S5)
# In order to align plots vertically, need to remove 1 sample from 16S that's not in ITS
input_filt_rare_16S_586 <- filter_data(input_filt_rare_16S,
                                       filter_cat = "sampleID",
                                       keep_vals = input_filt_rare_ITS$map_loaded$sampleID)
sum(input_filt_rare_16S_586$map_loaded$sampleID != input_filt_rare_ITS$map_loaded$sampleID) # Good

# Phylum
tax_sum_phyla_16S <- summarize_taxonomy(input = input_filt_rare_16S_586, 
                                        level = 2, 
                                        report_higher_tax = F)
bars_phyla_16S <- plot_taxa_bars(tax_sum_phyla_16S,
                                 input_filt_rare_16S_586$map_loaded,
                                 "sampleID",
                                 num_taxa = 12,
                                 data_only = TRUE) %>%
  mutate(taxon = fct_rev(taxon)) %>%
  left_join(., input_filt_rare_16S_586$map_loaded, by = c("group_by" = "sampleID"))
top_phyla_16S <- bars_phyla_16S %>%
  group_by(taxon) %>%
  summarise(mean = mean(mean_value)) %>%
  filter(taxon != "Other") %>%
  arrange(-mean) %>%
  mutate(taxon = as.character(taxon))
bars_phyla_16S <- bars_phyla_16S %>%
  mutate(taxon = factor(taxon,
                        levels = c("Other", rev(top_phyla_16S$taxon))))
fig4a <- ggplot(bars_phyla_16S, aes(group_by, mean_value, fill = taxon)) +
  geom_bar(stat = "identity", colour = NA, linewidth = 0.25) +
  labs(x = "Sample", y = "Relative abundance", fill = "Phylum") +
  scale_fill_manual(values = c("grey90", brewer.pal(12, "Paired")[12:1])) +
  scale_y_continuous(expand = c(0.01, 0.01)) + 
  facet_grid(~ Site, space = "free", scales = "free_x", 
             labeller = as_labeller(facet_names)) +
  theme_classic() +
  theme(axis.title.y = element_text(face = "bold", size = 12),
        axis.title.x = element_blank(),
        axis.text.y = element_text(size = 10),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.line.x = element_blank(),
        axis.line.y = element_blank(),
        strip.text = element_text(size = 4),
        strip.background = element_rect(linewidth = 0.2),
        legend.position = "right",
        legend.margin = margin(0, 0, 0, 0, unit = "pt"),
        legend.box.margin = margin(0, 0, 0, 0, unit = "pt"),
        legend.key.size = unit(0.3, "cm"),
        panel.spacing.x = unit(c(0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 
                                 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2), 
                               "lines"))
fig4a
leg4a <- get_legend(fig4a)
fig4a <- fig4a + theme(legend.position = "none")

tax_sum_phyla_ITS <- summarize_taxonomy(input = input_filt_rare_ITS, 
                                        level = 2, 
                                        report_higher_tax = F)
bars_phyla_ITS <- plot_taxa_bars(tax_sum_phyla_ITS,
                                 input_filt_rare_ITS$map_loaded,
                                 "sampleID",
                                 num_taxa = 12,
                                 data_only = TRUE) %>%
  mutate(taxon = fct_rev(taxon)) %>%
  left_join(., input_filt_rare_ITS$map_loaded, by = c("group_by" = "sampleID"))
top_phyla_ITS <- bars_phyla_ITS %>%
  group_by(taxon) %>%
  summarise(mean = mean(mean_value)) %>%
  filter(taxon != "Other") %>%
  arrange(-mean) %>%
  mutate(taxon = as.character(taxon))
bars_phyla_ITS <- bars_phyla_ITS %>%
  mutate(taxon = factor(taxon,
                        levels = c("Other", rev(top_phyla_ITS$taxon))))
fig4b <- ggplot(bars_phyla_ITS, aes(group_by, mean_value, fill = taxon)) +
  geom_bar(stat = "identity", colour = NA, linewidth = 0.25) +
  labs(x = "Sample", y = "Relative abundance", fill = "Phylum") +
  scale_fill_manual(values = c("grey90", brewer.pal(12, "Paired")[12:1])) +
  scale_y_continuous(expand = c(0.01, 0.01)) + 
  facet_grid(~ Site, space = "free", scales = "free_x", 
             labeller = as_labeller(facet_names)) +
  theme_classic() +
  theme(axis.title.y = element_text(face = "bold", size = 12),
        axis.title.x = element_blank(),
        axis.text.y = element_text(size = 10),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.line.x = element_blank(),
        axis.line.y = element_blank(),
        strip.text = element_blank(),
        strip.background = element_blank(),
        legend.position = "right",
        legend.margin = margin(0, 0, 0, 0, unit = "pt"),
        legend.box.margin = margin(0, 0, 0, 0, unit = "pt"),
        legend.key.size = unit(0.3, "cm"),
        panel.spacing.x = unit(c(0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 
                                 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2), 
                               "lines"))
fig4b
leg4b <- get_legend(fig4b)
fig4b <- fig4b + theme(legend.position = "none")

left <- plot_grid(fig4a, fig4b, ncol = 1, align = "v", labels = "auto")
left
right <- plot_grid(leg4a, leg4b, ncol = 1, align = "v")
right
pdf("FinalFigs/Figure4.pdf", width = 8, height = 6)
plot_grid(left, right, ncol = 2, rel_widths = c(0.8, 0.2))
dev.off()



# Genus
tax_sum_genera_16S <- summarize_taxonomy(input = input_filt_rare_16S_586, 
                                        level = 6, 
                                        report_higher_tax = F)
bars_genera_16S <- plot_taxa_bars(tax_sum_genera_16S,
                                 input_filt_rare_16S_586$map_loaded,
                                 "sampleID",
                                 num_taxa = 12,
                                 data_only = TRUE) %>%
  mutate(taxon = fct_rev(taxon)) %>%
  left_join(., input_filt_rare_16S_586$map_loaded, by = c("group_by" = "sampleID"))
top_genera_16S <- bars_genera_16S %>%
  group_by(taxon) %>%
  summarise(mean = mean(mean_value)) %>%
  filter(taxon != "Other") %>%
  filter(taxon != "NA") %>%
  arrange(-mean) %>%
  mutate(taxon = as.character(taxon))
bars_genera_16S <- bars_genera_16S %>%
  mutate(taxon = factor(taxon,
                        levels = c("NA", "Other", rev(top_genera_16S$taxon))))
figs4a <- ggplot(bars_genera_16S, aes(group_by, mean_value, fill = taxon)) +
  geom_bar(stat = "identity", colour = NA, linewidth = 0.25) +
  labs(x = "Sample", y = "Relative abundance", fill = "Genus") +
  scale_fill_manual(values = c("grey75", "grey90", brewer.pal(11, "Paired")[11:1])) +
  scale_y_continuous(expand = c(0.01, 0.01)) + 
  facet_grid(~ Site, space = "free", scales = "free_x", 
             labeller = as_labeller(facet_names)) +
  theme_classic() +
  theme(axis.title.y = element_text(face = "bold", size = 12),
        axis.title.x = element_blank(),
        axis.text.y = element_text(size = 10),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.line.x = element_blank(),
        axis.line.y = element_blank(),
        strip.text = element_text(size = 4),
        strip.background = element_rect(linewidth = 0.2),
        legend.position = "right",
        legend.margin = margin(0, 0, 0, 0, unit = "pt"),
        legend.box.margin = margin(0, 0, 0, 0, unit = "pt"),
        legend.key.size = unit(0.3, "cm"),
        legend.text = element_text(size = 5),
        panel.spacing.x = unit(c(0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 
                                 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2), 
                               "lines"))
figs4a
legs4a <- get_legend(figs4a)
figs4a <- figs4a + theme(legend.position = "none")

tax_sum_genera_ITS <- summarize_taxonomy(input = input_filt_rare_ITS, 
                                        level = 6, 
                                        report_higher_tax = F)
bars_genera_ITS <- plot_taxa_bars(tax_sum_genera_ITS,
                                 input_filt_rare_ITS$map_loaded,
                                 "sampleID",
                                 num_taxa = 12,
                                 data_only = TRUE) %>%
  mutate(taxon = fct_rev(taxon)) %>%
  left_join(., input_filt_rare_ITS$map_loaded, by = c("group_by" = "sampleID"))
top_genera_ITS <- bars_genera_ITS %>%
  group_by(taxon) %>%
  summarise(mean = mean(mean_value)) %>%
  filter(taxon != "Other") %>%
  filter(taxon != "NA") %>%
  arrange(-mean) %>%
  mutate(taxon = as.character(taxon))
bars_genera_ITS <- bars_genera_ITS %>%
  mutate(taxon = factor(taxon,
                        levels = c("NA", "Other", rev(top_genera_ITS$taxon))))
figs4b <- ggplot(bars_genera_ITS, aes(group_by, mean_value, fill = taxon)) +
  geom_bar(stat = "identity", colour = NA, size = 0.25) +
  labs(x = "Sample", y = "Relative abundance", fill = "Genus") +
  scale_fill_manual(values = c("grey75", "grey90", brewer.pal(11, "Paired")[11:1])) +
  scale_y_continuous(expand = c(0.01, 0.01)) + 
  facet_grid(~ Site, space = "free", scales = "free_x", 
             labeller = as_labeller(facet_names)) +
  theme_classic() +
  theme(axis.title.y = element_text(face = "bold", size = 12),
        axis.title.x = element_blank(),
        axis.text.y = element_text(size = 10),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.line.x = element_blank(),
        axis.line.y = element_blank(),
        strip.text = element_blank(),
        strip.background = element_blank(),
        legend.position = "right",
        legend.margin = margin(0, 0, 0, 0, unit = "pt"),
        legend.box.margin = margin(0, 0, 0, 0, unit = "pt"),
        legend.key.size = unit(0.3, "cm"),
        legend.text = element_text(size = 5),
        panel.spacing.x = unit(c(0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 
                                 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2), 
                               "lines"))
figs4b
legs4b <- get_legend(figs4b)
figs4b <- figs4b + theme(legend.position = "none")

left <- plot_grid(figs4a, figs4b, ncol = 1, align = "v", labels = "auto")
left
right <- plot_grid(legs4a, legs4b, ncol = 1, align = "v")
right
pdf("FinalFigs/FigureS4.pdf", width = 8.5, height = 6)
plot_grid(left, right, ncol = 2, rel_widths = c(0.8, 0.2))
dev.off()



# Guild
tax_sum_guilds_16S <- read.delim("data/faprotax_output.tsv") %>%
  column_to_rownames(var = "group") %>%
  set_names(names(tax_sum_phyla_16S_587)) %>%
  dplyr::select(names(tax_sum_phyla_16S))
tax_sum_guilds_16S_t <- as.data.frame(t(tax_sum_guilds_16S))
plot(tax_sum_guilds_16S_t$nitrification, tax_sum_guilds_16S_t$aerobic_ammonia_oxidation)
bars_guilds_16S <- plot_taxa_bars(tax_sum_guilds_16S,
                                  input_filt_rare_16S$map_loaded,
                                  "sampleID",
                                  num_taxa = 12,
                                  data_only = TRUE) %>%
  mutate(taxon = fct_rev(taxon)) %>%
  left_join(., input_filt_rare_16S$map_loaded, by = c("group_by" = "sampleID"))
top_guilds_16S <- bars_guilds_16S %>%
  group_by(taxon) %>%
  summarise(mean = mean(mean_value)) %>%
  filter(taxon != "Other") %>%
  filter(taxon != "Unassigned") %>%
  arrange(-mean) %>%
  mutate(taxon = as.character(taxon))
bars_guilds_16S <- bars_guilds_16S %>%
  mutate(taxon = factor(taxon,
                        levels = c("Unassigned", "Other", rev(top_guilds_16S$taxon))))
figs5a <- ggplot(bars_guilds_16S, aes(group_by, mean_value, fill = taxon)) +
  geom_bar(stat = "identity", colour = NA, size = 0.25) +
  labs(x = "Sample", y = "# of rarefied reads", fill = "Guild") +
  scale_fill_manual(values = c("grey75", "grey90", brewer.pal(12, "Paired")[12:1])) +
  scale_y_continuous(expand = c(0.01, 0.01)) + 
  facet_grid(~ Site, space = "free", scales = "free_x", 
             labeller = as_labeller(facet_names)) +
  theme_classic() +
  theme(axis.title.y = element_text(face = "bold", size = 12),
        axis.title.x = element_blank(),
        axis.text.y = element_text(size = 10),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.line.x = element_blank(),
        axis.line.y = element_blank(),
        strip.text = element_text(size = 4),
        strip.background = element_rect(linewidth = 0.2),
        legend.position = "right",
        legend.margin = margin(0, 0, 0, 0, unit = "pt"),
        legend.box.margin = margin(0, 0, 0, 0, unit = "pt"),
        legend.key.size = unit(0.3, "cm"),
        legend.text = element_text(size = 5),
        panel.spacing.x = unit(c(0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 
                                 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2), 
                               "lines"))
figs5a
legs5a <- get_legend(figs5a)
figs5a <- figs5a + theme(legend.position = "none")

tax_sum_guilds_ITS <- summarize_taxonomy(input = input_filt_rare_ITS, 
                                         level = 10, 
                                         report_higher_tax = F)
bars_guilds_ITS <- plot_taxa_bars(tax_sum_guilds_ITS,
                                  input_filt_rare_ITS$map_loaded,
                                  "sampleID",
                                  num_taxa = 13, # There are only 13 curated guilds anyways
                                  data_only = TRUE) %>%
  mutate(taxon = fct_rev(taxon)) %>%
  left_join(., input_filt_rare_ITS$map_loaded, by = c("group_by" = "sampleID"))
top_guilds_ITS <- bars_guilds_ITS %>%
  group_by(taxon) %>%
  summarise(mean = mean(mean_value)) %>%
  #filter(taxon != "Other") %>%
  filter(taxon != "Unassigned") %>%
  filter(taxon != "Multiple") %>%
  arrange(-mean) %>%
  mutate(taxon = as.character(taxon))
bars_guilds_ITS <- bars_guilds_ITS %>%
  mutate(taxon = factor(taxon,
                        levels = c("Unassigned", "Multiple", rev(top_guilds_ITS$taxon))))
figs5b <- ggplot(bars_guilds_ITS, aes(group_by, mean_value, fill = taxon)) +
  geom_bar(stat = "identity", colour = NA, size = 0.25) +
  labs(x = "Sample", y = "Relative abundance", fill = "Guild") +
  scale_fill_manual(values = c("grey75", "grey90", brewer.pal(11, "Paired")[11:1])) +
  scale_y_continuous(expand = c(0.01, 0.01)) + 
  facet_grid(~ Site, space = "free", scales = "free_x", 
             labeller = as_labeller(facet_names)) +
  theme_classic() +
  theme(axis.title.y = element_text(face = "bold", size = 12),
        axis.title.x = element_blank(),
        axis.text.y = element_text(size = 10),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.line.x = element_blank(),
        axis.line.y = element_blank(),
        strip.text = element_blank(),
        strip.background = element_blank(),
        legend.position = "right",
        legend.margin = margin(0, 0, 0, 0, unit = "pt"),
        legend.box.margin = margin(0, 0, 0, 0, unit = "pt"),
        legend.key.size = unit(0.3, "cm"),
        legend.text = element_text(size = 5),
        panel.spacing.x = unit(c(0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 
                                 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2), 
                               "lines"))
figs5b
legs5b <- get_legend(figs5b)
figs5b <- figs5b + theme(legend.position = "none")

left <- plot_grid(figs5a, figs5b, ncol = 1, align = "v", labels = "auto")
left
right <- plot_grid(legs5a, legs5b, ncol = 1, align = "v")
right
pdf("FinalFigs/FigureS5.pdf", width = 8, height = 6)
plot_grid(left, right, ncol = 2, rel_widths = c(0.8, 0.2))
dev.off()



# Prevalence
ncol(tax_sum_asv_16S)
prev_16S <- data.frame("ASV_ID" = rownames(tax_sum_asv_16S),
                       "Absent" = rowSums(tax_sum_asv_16S==0)) %>%
  mutate(Present = ncol(tax_sum_asv_16S) - Absent) %>%
  mutate(Present_Perc = Present/ncol(tax_sum_asv_16S)*100)
top_prev_16S <- prev_16S %>%
  #filter(Present_Perc > 80) %>%
  arrange(desc(Present_Perc)) %>%
  slice_head(n = 20) %>%
  left_join(., input_filt_rare_16S$taxonomy_loaded, by = c("ASV_ID" = "taxonomy7")) %>%
  mutate(taxonomy = paste(taxonomy2, taxonomy5, taxonomy6, sep = "; "))
figs7a <- ggplot(top_prev_16S, aes(reorder(ASV_ID, Present_Perc, median), Present_Perc, 
                                   fill = taxonomy2)) +
  geom_bar(stat = "identity", alpha = 0.3) +
  geom_text(data = top_prev_16S,
            aes(ASV_ID, 1, label = taxonomy),
            size = 2.5, angle = 90, inherit.aes = F, hjust = 0) +
  scale_x_discrete(limits = rev) +
  scale_y_continuous(expand = c(0.01, 0.01),
                     limits = c(0, 100)) +
  scale_fill_brewer(palette = "Paired") +
  labs(x = "ASV_ID",
       y = "% Prevalence") +
  theme_bw() +
  theme(legend.position = "none",
        axis.title.y = element_text(size = 12, face = "bold"),
        axis.title.x = element_blank(),
        axis.text.y = element_text(size = 10),
        axis.text.x = element_text(size = 8, angle = 45, hjust = 1),
        panel.grid = element_blank())
figs7a

ncol(tax_sum_asv_ITS)
prev_ITS <- data.frame("ASV_ID" = rownames(tax_sum_asv_ITS),
                       "Absent" = rowSums(tax_sum_asv_ITS==0)) %>%
  mutate(Present = ncol(tax_sum_asv_ITS) - Absent) %>%
  mutate(Present_Perc = Present/ncol(tax_sum_asv_ITS)*100)
top_prev_ITS <- prev_ITS %>%
  #filter(Present_Perc > 70) %>%
  arrange(desc(Present_Perc)) %>%
  slice_head(n = 20) %>%
  left_join(., input_filt_rare_ITS$taxonomy_loaded, by = c("ASV_ID" = "taxonomy8")) %>%
  mutate(taxonomy = paste(taxonomy5, taxonomy6, taxonomy7, sep = "; ")) %>%
  mutate(taxonomy = ifelse(taxonomy == "NA; NA; NA",
                           taxonomy3,
                           taxonomy))
figs7b <- ggplot(top_prev_ITS, aes(reorder(ASV_ID, Present_Perc, median), Present_Perc, 
                                   fill = taxonomy3)) +
  geom_bar(stat = "identity", alpha = 0.3) +
  geom_text(data = top_prev_ITS,
            aes(ASV_ID, 1, label = taxonomy),
            size = 2.5, angle = 90, inherit.aes = F, hjust = 0) +
  scale_x_discrete(limits = rev) +
  scale_y_continuous(expand = c(0.01, 0.01),
                     limits = c(0, 100)) +
  scale_fill_brewer(palette = "Paired") +
  labs(x = "ASV_ID",
       y = "% Prevalence") +
  theme_bw() +
  theme(legend.position = "none",
        axis.title.y = element_blank(),
        axis.title.x = element_blank(),
        axis.text.y = element_blank(),
        axis.text.x = element_text(size = 8, angle = 45, hjust = 1),
        panel.grid = element_blank())
figs7b

pdf("FinalFigs/FigureS7.pdf", width = 8, height = 6)
plot_grid(figs7a, figs7b, ncol = 2, rel_widths = c(0.53, 0.47),
          labels = "auto", label_x = c(0.9, 0.9), label_y = c(0.99, 0.99), align = "h")
dev.off()


#### _Of Interest ####
# Test and plot

### Sclerotinia
# Note: Taxonomy is Fungi, Ascomycota, Leotiomycetes, Helotiales, Sclerotiniaceae
tax_sum_genus_ITS_t <- as.data.frame(t(tax_sum_genus_ITS))
input_filt_rare_ITS$map_loaded$SclerotiniaAbund <- tax_sum_genus_ITS_t$Sclerotinia
leveneTest(input_filt_rare_ITS$map_loaded$SclerotiniaAbund ~ input_filt_rare_ITS$map_loaded$Site) # Homogeneous
leveneTest(input_filt_rare_ITS$map_loaded$SclerotiniaAbund ~ input_filt_rare_ITS$map_loaded$Entry) # Homogeneous
m <- aov(SclerotiniaAbund ~ Site * Entry, data = input_filt_rare_ITS$map_loaded)
summary(m) # All
Anova(m, type = "III", singular.ok = TRUE) # Interaction
hist(m$residuals)
shapiro.test(m$residuals)
plot(m$fitted.values, m$residuals)
m <- aov(SclerotiniaAbund ~ Entry, data = input_filt_rare_ITS$map_loaded)
summary(m)
m <- aov(SclerotiniaAbund ~ Site, data = input_filt_rare_ITS$map_loaded)
summary(m)
TukeyHSD(m)
t <- emmeans(object = m, specs = "Site") %>%
  cld(object = ., adjust = "Tukey", Letters = letters, alpha = 0.05) %>%
  mutate(name = "SclerotiniaAbund",
         y = max(input_filt_rare_ITS$map_loaded$SclerotiniaAbund)+
           (max(input_filt_rare_ITS$map_loaded$SclerotiniaAbund)-
              min(input_filt_rare_ITS$map_loaded$SclerotiniaAbund))/10)

pdf("InitialFigs/Taxa_Sclerotinia.pdf", width = 7, height = 5)
ggplot(input_filt_rare_ITS$map_loaded, 
       aes(reorder(Site, SclerotiniaAbund, mean), SclerotiniaAbund)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(size = 2, alpha = 0.75, width = 0.2, pch = 16) +
  geom_text(data = t, aes(Site, y, label = str_trim(.group)), 
            size = 3, color = "black") +
  labs(x = "Site", y = "Sclerotinia relative abundance") +
  theme_bw() +
  theme(legend.position = "right",
        axis.title = element_text(size = 12),
        axis.text.y = element_text(size = 10),
        axis.text.x = element_text(size = 10, angle = 45, hjust = 1),
        strip.text = element_text(size = 14),
        strip.background = element_rect(fill = "white"))
dev.off()

# Label points to check high ones
ggplot(input_filt_rare_ITS$map_loaded, 
       aes(reorder(Site, SclerotiniaAbund, mean), SclerotiniaAbund)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(size = 2, alpha = 0.75, width = 0.2, pch = 16) +
  geom_text(data = t, aes(Site, y, label = str_trim(.group)), 
            size = 3, color = "black") +
  geom_text(aes(Site, SclerotiniaAbund, label = Entry)) +
  labs(x = "Site", y = "Sclerotinia relative abundance") +
  theme_bw() +
  theme(legend.position = "right",
        axis.title = element_text(size = 12),
        axis.text.y = element_text(size = 10),
        axis.text.x = element_text(size = 10, angle = 45, hjust = 1),
        strip.text = element_text(size = 14),
        strip.background = element_rect(fill = "white"))

sum(input_filt_rare_ITS$map_loaded$SclerotiniaAbund > 0) # Present in 21 samples, 7 sites.



### ASV1 (Podila)
input_filt_rare_ITS$map_loaded$ASV1 <- tax_sum_asv_ITS_t$ASV_1
leveneTest(input_filt_rare_ITS$map_loaded$ASV1 ~ input_filt_rare_ITS$map_loaded$Site) # Not Homogeneous
leveneTest(input_filt_rare_ITS$map_loaded$ASV1 ~ input_filt_rare_ITS$map_loaded$Entry) # Homogeneous
m <- aov(ASV1 ~ Site * Entry, data = input_filt_rare_ITS$map_loaded)
summary(m)
Anova(m, type = "III", singular.ok = TRUE) # Site
hist(m$residuals)
shapiro.test(m$residuals)
plot(m$fitted.values, m$residuals)
m <- aov(ASV1 ~ Entry, data = input_filt_rare_ITS$map_loaded)
summary(m)
m <- aov(ASV1 ~ Site, data = input_filt_rare_ITS$map_loaded)
summary(m)
TukeyHSD(m)
t <- emmeans(object = m, specs = "Site") %>%
  cld(object = ., adjust = "Tukey", Letters = letters, alpha = 0.05) %>%
  mutate(name = "ASV1",
         y = max(input_filt_rare_ITS$map_loaded$ASV1)+
           (max(input_filt_rare_ITS$map_loaded$ASV1)-
              min(input_filt_rare_ITS$map_loaded$ASV1))/10)

pdf("InitialFigs/Taxa_Podila.pdf", width = 7, height = 5)
ggplot(input_filt_rare_ITS$map_loaded, 
       aes(reorder(Site, ASV1, mean), ASV1)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(size = 2, alpha = 0.75, width = 0.2, pch = 16) +
  geom_text(data = t, aes(Site, y, label = str_trim(.group)), 
            size = 3, color = "black") +
  labs(x = "Site", y = "Podila minutissima ASV_1 relative abundance") +
  theme_bw() +
  theme(legend.position = "right",
        axis.title = element_text(size = 12),
        axis.text.y = element_text(size = 10),
        axis.text.x = element_text(size = 10, angle = 45, hjust = 1),
        strip.text = element_text(size = 14),
        strip.background = element_rect(fill = "white"))
dev.off()

sum(input_filt_rare_ITS$map_loaded$ASV1 > 0) # Present in 358 samples, 7 sites.



### Naganishia
input_filt_rare_ITS$map_loaded$Naganishia <- tax_sum_genus_ITS_t$Naganishia
leveneTest(input_filt_rare_ITS$map_loaded$Naganishia ~ input_filt_rare_ITS$map_loaded$Site) # Not Homogeneous
leveneTest(input_filt_rare_ITS$map_loaded$Naganishia ~ input_filt_rare_ITS$map_loaded$Entry) # Homogeneous
m <- aov(Naganishia ~ Site * Entry, data = input_filt_rare_ITS$map_loaded)
summary(m) # All
Anova(m, type = "III", singular.ok = TRUE) # Interaction
hist(m$residuals)
shapiro.test(m$residuals)
plot(m$fitted.values, m$residuals)
m <- aov(Naganishia ~ Entry, data = input_filt_rare_ITS$map_loaded)
summary(m)
m <- aov(Naganishia ~ Site, data = input_filt_rare_ITS$map_loaded)
summary(m)
TukeyHSD(m)
t <- emmeans(object = m, specs = "Site") %>%
  cld(object = ., adjust = "Tukey", Letters = letters, alpha = 0.05) %>%
  mutate(name = "Naganishia",
         y = max(input_filt_rare_ITS$map_loaded$Naganishia)+
           (max(input_filt_rare_ITS$map_loaded$Naganishia)-
              min(input_filt_rare_ITS$map_loaded$Naganishia))/10)
pdf("InitialFigs/Taxa_Naganishia.pdf", width = 7, height = 5)
ggplot(input_filt_rare_ITS$map_loaded, 
       aes(reorder(Site, Naganishia, mean), Naganishia)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(size = 2, alpha = 0.75, width = 0.2, pch = 16) +
  geom_text(data = t, aes(Site, y, label = str_trim(.group)), 
            size = 3, color = "black") +
  labs(x = "Site", y = "Naganishia relative abundance") +
  theme_bw() +
  theme(legend.position = "right",
        axis.title = element_text(size = 12),
        axis.text.y = element_text(size = 10),
        axis.text.x = element_text(size = 10, angle = 45, hjust = 1),
        strip.text = element_text(size = 14),
        strip.background = element_rect(fill = "white"))
dev.off()
sum(input_filt_rare_ITS$map_loaded$Naganishia > 0) # Present in 64 samples
sum_nag <- input_filt_rare_ITS$map_loaded %>%
  group_by(Site) %>%
  summarise(Naganishia = mean(Naganishia)) # Present in 7 sites.



### Saprotroph
tax_sum_guilds_ITS_t <- as.data.frame(t(tax_sum_guilds_ITS))
input_filt_rare_ITS$map_loaded$Saprotroph <- tax_sum_guilds_ITS_t$Saprotroph
leveneTest(input_filt_rare_ITS$map_loaded$Saprotroph ~ input_filt_rare_ITS$map_loaded$Site) # NH
leveneTest(input_filt_rare_ITS$map_loaded$Saprotroph ~ input_filt_rare_ITS$map_loaded$Entry) # H
m <- aov(Saprotroph ~ Site * Entry, data = input_filt_rare_ITS$map_loaded)
summary(m) # Site***, Entry., Int*
Anova(m, type = "III", singular.ok = TRUE) # Site***, Int*
hist(m$residuals)
shapiro.test(m$residuals)
plot(m$fitted.values, m$residuals)
m <- aov(Saprotroph ~ Entry, data = input_filt_rare_ITS$map_loaded)
summary(m)
m <- aov(Saprotroph ~ Site, data = input_filt_rare_ITS$map_loaded)
summary(m)
TukeyHSD(m)
t <- emmeans(object = m, specs = "Site") %>%
  cld(object = ., adjust = "Tukey", Letters = letters, alpha = 0.05) %>%
  mutate(name = "Saprotroph",
         y = max(input_filt_rare_ITS$map_loaded$Saprotroph)+
           (max(input_filt_rare_ITS$map_loaded$Saprotroph)-
              min(input_filt_rare_ITS$map_loaded$Saprotroph))/10)
pdf("InitialFigs/Taxa_Guilds_Saprotroph.pdf", width = 7, height = 5)
ggplot(input_filt_rare_ITS$map_loaded, 
       aes(reorder(Site, Saprotroph, mean), Saprotroph)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(size = 2, alpha = 0.75, width = 0.2, pch = 16) +
  geom_text(data = t, aes(Site, y, label = str_trim(.group)), 
            size = 3, color = "black") +
  labs(x = "Site", y = "Saprotroph relative abundance") +
  theme_bw() +
  theme(legend.position = "right",
        axis.title = element_text(size = 12),
        axis.text.y = element_text(size = 10),
        axis.text.x = element_text(size = 10, angle = 45, hjust = 1),
        strip.text = element_text(size = 14),
        strip.background = element_rect(fill = "white"))
dev.off()
sum(input_filt_rare_ITS$map_loaded$Saprotroph > 0) # Present in all samples
sum_sap <- input_filt_rare_ITS$map_loaded %>%
  group_by(Site) %>%
  summarise(Saprotroph = mean(Saprotroph)) 
sum(sum_sap$Saprotroph > 0) # Present in 15 sites.



### Parasite
input_filt_rare_ITS$map_loaded$Parasite <- tax_sum_guilds_ITS_t$Parasite
leveneTest(input_filt_rare_ITS$map_loaded$Parasite ~ input_filt_rare_ITS$map_loaded$Site) # NH
leveneTest(input_filt_rare_ITS$map_loaded$Parasite ~ input_filt_rare_ITS$map_loaded$Entry) # H
m <- aov(Parasite ~ Site * Entry, data = input_filt_rare_ITS$map_loaded)
summary(m) # Site***
Anova(m, type = "III", singular.ok = TRUE) # Site.
hist(m$residuals)
shapiro.test(m$residuals)
plot(m$fitted.values, m$residuals)
m <- aov(Parasite ~ Entry, data = input_filt_rare_ITS$map_loaded)
summary(m)
m <- aov(Parasite ~ Site, data = input_filt_rare_ITS$map_loaded)
summary(m)
TukeyHSD(m)
t <- emmeans(object = m, specs = "Site") %>%
  cld(object = ., adjust = "Tukey", Letters = letters, alpha = 0.05) %>%
  mutate(name = "Parasite",
         y = max(input_filt_rare_ITS$map_loaded$Parasite)+
           (max(input_filt_rare_ITS$map_loaded$Parasite)-
              min(input_filt_rare_ITS$map_loaded$Parasite))/10)
pdf("InitialFigs/Taxa__Guilds_Parasite.pdf", width = 7, height = 5)
ggplot(input_filt_rare_ITS$map_loaded, 
       aes(reorder(Site, Parasite, mean), Parasite)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(size = 2, alpha = 0.75, width = 0.2, pch = 16) +
  geom_text(data = t, aes(Site, y, label = str_trim(.group)), 
            size = 3, color = "black") +
  labs(x = "Site", y = "Parasite relative abundance") +
  theme_bw() +
  theme(legend.position = "right",
        axis.title = element_text(size = 12),
        axis.text.y = element_text(size = 10),
        axis.text.x = element_text(size = 10, angle = 45, hjust = 1),
        strip.text = element_text(size = 14),
        strip.background = element_rect(fill = "white"))
dev.off()
sum(input_filt_rare_ITS$map_loaded$Parasite > 0) # Present in 476 samples
sum_par <- input_filt_rare_ITS$map_loaded %>%
  group_by(Site) %>%
  summarise(Parasite = mean(Parasite)) 
sum(sum_par$Parasite > 0) # Present in 15 sites.



### Parasite
input_filt_rare_ITS$map_loaded$Parasite <- tax_sum_guilds_ITS_t$Parasite
leveneTest(input_filt_rare_ITS$map_loaded$Parasite ~ input_filt_rare_ITS$map_loaded$Site) # NH
leveneTest(input_filt_rare_ITS$map_loaded$Parasite ~ input_filt_rare_ITS$map_loaded$Entry) # H
m <- aov(Parasite ~ Site * Entry, data = input_filt_rare_ITS$map_loaded)
summary(m) # Site***
Anova(m, type = "III", singular.ok = TRUE) # Site.
hist(m$residuals)
shapiro.test(m$residuals)
plot(m$fitted.values, m$residuals)
m <- aov(Parasite ~ Entry, data = input_filt_rare_ITS$map_loaded)
summary(m)
m <- aov(Parasite ~ Site, data = input_filt_rare_ITS$map_loaded)
summary(m)
TukeyHSD(m)
t <- emmeans(object = m, specs = "Site") %>%
  cld(object = ., adjust = "Tukey", Letters = letters, alpha = 0.05) %>%
  mutate(name = "Parasite",
         y = max(input_filt_rare_ITS$map_loaded$Parasite)+
           (max(input_filt_rare_ITS$map_loaded$Parasite)-
              min(input_filt_rare_ITS$map_loaded$Parasite))/10)
pdf("InitialFigs/Taxa_Guilds_Parasite.pdf", width = 7, height = 5)
ggplot(input_filt_rare_ITS$map_loaded, 
       aes(reorder(Site, Parasite, mean), Parasite)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(size = 2, alpha = 0.75, width = 0.2, pch = 16) +
  geom_text(data = t, aes(Site, y, label = str_trim(.group)), 
            size = 3, color = "black") +
  labs(x = "Site", y = "Parasite relative abundance") +
  theme_bw() +
  theme(legend.position = "right",
        axis.title = element_text(size = 12),
        axis.text.y = element_text(size = 10),
        axis.text.x = element_text(size = 10, angle = 45, hjust = 1),
        strip.text = element_text(size = 14),
        strip.background = element_rect(fill = "white"))
dev.off()
sum(input_filt_rare_ITS$map_loaded$Parasite > 0) # Present in 476 samples
sum_par <- input_filt_rare_ITS$map_loaded %>%
  group_by(Site) %>%
  summarise(Parasite = mean(Parasite)) 
sum(sum_par$Parasite > 0) # Present in 15 sites.



### AMF
input_filt_rare_ITS$map_loaded$AMF <- tax_sum_guilds_ITS_t$`Arbuscular Mycorrhizal`
leveneTest(input_filt_rare_ITS$map_loaded$AMF ~ input_filt_rare_ITS$map_loaded$Site) # NH
leveneTest(input_filt_rare_ITS$map_loaded$AMF ~ input_filt_rare_ITS$map_loaded$Entry) # H
m <- aov(AMF ~ Site * Entry, data = input_filt_rare_ITS$map_loaded)
summary(m) # Site***
Anova(m, type = "III", singular.ok = TRUE) # Site.
hist(m$residuals)
shapiro.test(m$residuals)
plot(m$fitted.values, m$residuals)
m <- aov(AMF ~ Entry, data = input_filt_rare_ITS$map_loaded)
summary(m)
m <- aov(AMF ~ Site, data = input_filt_rare_ITS$map_loaded)
summary(m)
TukeyHSD(m)
t <- emmeans(object = m, specs = "Site") %>%
  cld(object = ., adjust = "Tukey", Letters = letters, alpha = 0.05) %>%
  mutate(name = "AMF",
         y = max(input_filt_rare_ITS$map_loaded$AMF)+
           (max(input_filt_rare_ITS$map_loaded$AMF)-
              min(input_filt_rare_ITS$map_loaded$AMF))/10)
pdf("InitialFigs/Taxa_Guilds_AMF.pdf", width = 7, height = 5)
ggplot(input_filt_rare_ITS$map_loaded, 
       aes(reorder(Site, AMF, mean), AMF)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(size = 2, alpha = 0.75, width = 0.2, pch = 16) +
  geom_text(data = t, aes(Site, y, label = str_trim(.group)), 
            size = 3, color = "black") +
  labs(x = "Site", y = "AMF relative abundance") +
  theme_bw() +
  theme(legend.position = "right",
        axis.title = element_text(size = 12),
        axis.text.y = element_text(size = 10),
        axis.text.x = element_text(size = 10, angle = 45, hjust = 1),
        strip.text = element_text(size = 14),
        strip.background = element_rect(fill = "white"))
dev.off()
sum(input_filt_rare_ITS$map_loaded$AMF > 0) # Present in 476 samples
sum_amf <- input_filt_rare_ITS$map_loaded %>%
  group_by(Site) %>%
  summarise(AMF = mean(AMF)) 
sum(sum_amf$AMF > 0) # Present in 15 sites.



### Pathogen
input_filt_rare_ITS$map_loaded$Pathogen <- tax_sum_guilds_ITS_t$`Plant Pathogen`
leveneTest(input_filt_rare_ITS$map_loaded$Pathogen ~ input_filt_rare_ITS$map_loaded$Site) # NH
leveneTest(input_filt_rare_ITS$map_loaded$Pathogen ~ input_filt_rare_ITS$map_loaded$Entry) # H
m <- aov(Pathogen ~ Site * Entry, data = input_filt_rare_ITS$map_loaded)
summary(m) # Site***
Anova(m, type = "III", singular.ok = TRUE) # Site.
hist(m$residuals)
shapiro.test(m$residuals)
plot(m$fitted.values, m$residuals)
m <- aov(Pathogen ~ Entry, data = input_filt_rare_ITS$map_loaded)
summary(m)
m <- aov(Pathogen ~ Site, data = input_filt_rare_ITS$map_loaded)
summary(m)
TukeyHSD(m)
t <- emmeans(object = m, specs = "Site") %>%
  cld(object = ., adjust = "Tukey", Letters = letters, alpha = 0.05) %>%
  mutate(name = "Pathogen",
         y = max(input_filt_rare_ITS$map_loaded$Pathogen)+
           (max(input_filt_rare_ITS$map_loaded$Pathogen)-
              min(input_filt_rare_ITS$map_loaded$Pathogen))/10)
pdf("InitialFigs/Taxa_Guilds_Pathogen.pdf", width = 7, height = 5)
ggplot(input_filt_rare_ITS$map_loaded, 
       aes(reorder(Site, Pathogen, mean), Pathogen)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(size = 2, alpha = 0.75, width = 0.2, pch = 16) +
  geom_text(data = t, aes(Site, y, label = str_trim(.group)), 
            size = 3, color = "black") +
  labs(x = "Site", y = "Plant Pathogen relative abundance") +
  theme_bw() +
  theme(legend.position = "right",
        axis.title = element_text(size = 12),
        axis.text.y = element_text(size = 10),
        axis.text.x = element_text(size = 10, angle = 45, hjust = 1),
        strip.text = element_text(size = 14),
        strip.background = element_rect(fill = "white"))
dev.off()
sum(input_filt_rare_ITS$map_loaded$Pathogen > 0) # Present in 515 samples
sum_pat <- input_filt_rare_ITS$map_loaded %>%
  group_by(Site) %>%
  summarise(Pathogen = mean(Pathogen)) 
sum(sum_pat$Pathogen > 0) # Present in 15 sites.



### Epiphyte
input_filt_rare_ITS$map_loaded$Epiphyte <- tax_sum_guilds_ITS_t$Epiphyte
leveneTest(input_filt_rare_ITS$map_loaded$Epiphyte ~ input_filt_rare_ITS$map_loaded$Site) # NH
leveneTest(input_filt_rare_ITS$map_loaded$Epiphyte ~ input_filt_rare_ITS$map_loaded$Entry) # H
m <- aov(Epiphyte ~ Site * Entry, data = input_filt_rare_ITS$map_loaded)
summary(m) # Site***
Anova(m, type = "III", singular.ok = TRUE) # Site***
hist(m$residuals)
shapiro.test(m$residuals)
plot(m$fitted.values, m$residuals)
m <- aov(Epiphyte ~ Entry, data = input_filt_rare_ITS$map_loaded)
summary(m)
m <- aov(Epiphyte ~ Site, data = input_filt_rare_ITS$map_loaded)
summary(m)
TukeyHSD(m)
t <- emmeans(object = m, specs = "Site") %>%
  cld(object = ., adjust = "Tukey", Letters = letters, alpha = 0.05) %>%
  mutate(name = "Epiphyte",
         y = max(input_filt_rare_ITS$map_loaded$Epiphyte)+
           (max(input_filt_rare_ITS$map_loaded$Epiphyte)-
              min(input_filt_rare_ITS$map_loaded$Epiphyte))/10)
pdf("InitialFigs/Taxa_Guilds_Epiphyte.pdf", width = 7, height = 5)
ggplot(input_filt_rare_ITS$map_loaded, 
       aes(reorder(Site, Epiphyte, mean), Epiphyte)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(size = 2, alpha = 0.75, width = 0.2, pch = 16) +
  geom_text(data = t, aes(Site, y, label = str_trim(.group)), 
            size = 3, color = "black") +
  labs(x = "Site", y = "Epiphyte relative abundance") +
  theme_bw() +
  theme(legend.position = "right",
        axis.title = element_text(size = 12),
        axis.text.y = element_text(size = 10),
        axis.text.x = element_text(size = 10, angle = 45, hjust = 1),
        strip.text = element_text(size = 14),
        strip.background = element_rect(fill = "white"))
dev.off()
sum(input_filt_rare_ITS$map_loaded$Epiphyte > 0) # Present in 288 samples
sum_epi <- input_filt_rare_ITS$map_loaded %>%
  group_by(Site) %>%
  summarise(Epiphyte = mean(Epiphyte)) 
sum(sum_epi$Epiphyte > 0) # Present in 12 sites.



### Ectomycorrhizal
input_filt_rare_ITS$map_loaded$Ectomycorrhizal <- tax_sum_guilds_ITS_t$Ectomycorrhizal
leveneTest(input_filt_rare_ITS$map_loaded$Ectomycorrhizal ~ input_filt_rare_ITS$map_loaded$Site) # NH
leveneTest(input_filt_rare_ITS$map_loaded$Ectomycorrhizal ~ input_filt_rare_ITS$map_loaded$Entry) # H
m <- aov(Ectomycorrhizal ~ Site * Entry, data = input_filt_rare_ITS$map_loaded)
summary(m) # Site***
Anova(m, type = "III", singular.ok = TRUE) # N.S.
hist(m$residuals)
shapiro.test(m$residuals)
plot(m$fitted.values, m$residuals)
m <- aov(Ectomycorrhizal ~ Entry, data = input_filt_rare_ITS$map_loaded)
summary(m)
m <- aov(Ectomycorrhizal ~ Site, data = input_filt_rare_ITS$map_loaded)
summary(m)
TukeyHSD(m)
t <- emmeans(object = m, specs = "Site") %>%
  cld(object = ., adjust = "Tukey", Letters = letters, alpha = 0.05) %>%
  mutate(name = "Ectomycorrhizal",
         y = max(input_filt_rare_ITS$map_loaded$Ectomycorrhizal)+
           (max(input_filt_rare_ITS$map_loaded$Ectomycorrhizal)-
              min(input_filt_rare_ITS$map_loaded$Ectomycorrhizal))/10)
pdf("InitialFigs/Taxa_Guilds_Ectomycorrhizal.pdf", width = 7, height = 5)
ggplot(input_filt_rare_ITS$map_loaded, 
       aes(reorder(Site, Ectomycorrhizal, mean), Ectomycorrhizal)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(size = 2, alpha = 0.75, width = 0.2, pch = 16) +
  geom_text(data = t, aes(Site, y, label = str_trim(.group)), 
            size = 3, color = "black") +
  labs(x = "Site", y = "Ectomycorrhizal relative abundance") +
  theme_bw() +
  theme(legend.position = "right",
        axis.title = element_text(size = 12),
        axis.text.y = element_text(size = 10),
        axis.text.x = element_text(size = 10, angle = 45, hjust = 1),
        strip.text = element_text(size = 14),
        strip.background = element_rect(fill = "white"))
dev.off()
sum(input_filt_rare_ITS$map_loaded$Ectomycorrhizal > 0) # Present in 149 samples
sum_ect <- input_filt_rare_ITS$map_loaded %>%
  group_by(Site) %>%
  summarise(Ectomycorrhizal = mean(Ectomycorrhizal)) 
sum(sum_ect$Ectomycorrhizal > 0) # Present in 15 sites.



### DSE
input_filt_rare_ITS$map_loaded$DSE <- tax_sum_guilds_ITS_t$`Dark Septate Endophyte`
leveneTest(input_filt_rare_ITS$map_loaded$DSE ~ input_filt_rare_ITS$map_loaded$Site) # NH
leveneTest(input_filt_rare_ITS$map_loaded$DSE ~ input_filt_rare_ITS$map_loaded$Entry) # H
m <- aov(DSE ~ Site * Entry, data = input_filt_rare_ITS$map_loaded)
summary(m) # Site*
Anova(m, type = "III", singular.ok = TRUE) # N.S.
hist(m$residuals)
shapiro.test(m$residuals)
plot(m$fitted.values, m$residuals)
m <- aov(DSE ~ Entry, data = input_filt_rare_ITS$map_loaded)
summary(m)
m <- aov(DSE ~ Site, data = input_filt_rare_ITS$map_loaded)
summary(m)
TukeyHSD(m)
t <- emmeans(object = m, specs = "Site") %>%
  cld(object = ., adjust = "Tukey", Letters = letters, alpha = 0.05) %>%
  mutate(name = "DSE",
         y = max(input_filt_rare_ITS$map_loaded$DSE)+
           (max(input_filt_rare_ITS$map_loaded$DSE)-
              min(input_filt_rare_ITS$map_loaded$DSE))/10)
pdf("InitialFigs/Taxa_Guilds_DSE.pdf", width = 7, height = 5)
ggplot(input_filt_rare_ITS$map_loaded, 
       aes(reorder(Site, DSE, mean), DSE)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(size = 2, alpha = 0.75, width = 0.2, pch = 16) +
  geom_text(data = t, aes(Site, y, label = str_trim(.group)), 
            size = 3, color = "black") +
  labs(x = "Site", y = "DSE relative abundance") +
  theme_bw() +
  theme(legend.position = "right",
        axis.title = element_text(size = 12),
        axis.text.y = element_text(size = 10),
        axis.text.x = element_text(size = 10, angle = 45, hjust = 1),
        strip.text = element_text(size = 14),
        strip.background = element_rect(fill = "white"))
dev.off()
sum(input_filt_rare_ITS$map_loaded$DSE > 0) # Present in 79 samples
sum_dse <- input_filt_rare_ITS$map_loaded %>%
  group_by(Site) %>%
  summarise(DSE = mean(DSE)) 
sum(sum_dse$DSE > 0) # Present in 15 sites.



### Endophyte
input_filt_rare_ITS$map_loaded$Endophyte <- tax_sum_guilds_ITS_t$Endophyte
leveneTest(input_filt_rare_ITS$map_loaded$Endophyte ~ input_filt_rare_ITS$map_loaded$Site) # NH
leveneTest(input_filt_rare_ITS$map_loaded$Endophyte ~ input_filt_rare_ITS$map_loaded$Entry) # H
m <- aov(Endophyte ~ Site * Entry, data = input_filt_rare_ITS$map_loaded)
summary(m) # Site*
Anova(m, type = "III", singular.ok = TRUE) # Site.
hist(m$residuals)
shapiro.test(m$residuals)
plot(m$fitted.values, m$residuals)
m <- aov(Endophyte ~ Entry, data = input_filt_rare_ITS$map_loaded)
summary(m)
m <- aov(Endophyte ~ Site, data = input_filt_rare_ITS$map_loaded)
summary(m)
TukeyHSD(m)
t <- emmeans(object = m, specs = "Site") %>%
  cld(object = ., adjust = "Tukey", Letters = letters, alpha = 0.05) %>%
  mutate(name = "Endophyte",
         y = max(input_filt_rare_ITS$map_loaded$Endophyte)+
           (max(input_filt_rare_ITS$map_loaded$Endophyte)-
              min(input_filt_rare_ITS$map_loaded$Endophyte))/10)
pdf("InitialFigs/Taxa_Guilds_Endophyte.pdf", width = 7, height = 5)
ggplot(input_filt_rare_ITS$map_loaded, 
       aes(reorder(Site, Endophyte, mean), Endophyte)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(size = 2, alpha = 0.75, width = 0.2, height = 0, pch = 16) +
  geom_text(data = t, aes(Site, y, label = str_trim(.group)), 
            size = 3, color = "black") +
  labs(x = "Site", y = "Endophyte relative abundance") +
  theme_bw() +
  theme(legend.position = "right",
        axis.title = element_text(size = 12),
        axis.text.y = element_text(size = 10),
        axis.text.x = element_text(size = 10, angle = 45, hjust = 1),
        strip.text = element_text(size = 14),
        strip.background = element_rect(fill = "white"))
dev.off()
sum(input_filt_rare_ITS$map_loaded$Endophyte > 0) # Present in 127 samples
sum_Endophyte <- input_filt_rare_ITS$map_loaded %>%
  group_by(Site) %>%
  summarise(Endophyte = mean(Endophyte)) 
sum(sum_Endophyte$Endophyte > 0) # Present in 13 sites.



### Lichenized
input_filt_rare_ITS$map_loaded$Lichenized <- tax_sum_guilds_ITS_t$Lichenized
leveneTest(input_filt_rare_ITS$map_loaded$Lichenized ~ input_filt_rare_ITS$map_loaded$Site) # H
leveneTest(input_filt_rare_ITS$map_loaded$Lichenized ~ input_filt_rare_ITS$map_loaded$Entry) # H
m <- aov(Lichenized ~ Site * Entry, data = input_filt_rare_ITS$map_loaded)
summary(m) # N.S.
Anova(m, type = "III", singular.ok = TRUE) # N.S.
hist(m$residuals)
shapiro.test(m$residuals)
plot(m$fitted.values, m$residuals)
m <- aov(Lichenized ~ Entry, data = input_filt_rare_ITS$map_loaded)
summary(m)
m <- aov(Lichenized ~ Site, data = input_filt_rare_ITS$map_loaded)
summary(m)
TukeyHSD(m)
t <- emmeans(object = m, specs = "Site") %>%
  cld(object = ., adjust = "Tukey", Letters = letters, alpha = 0.05) %>%
  mutate(name = "Lichenized",
         y = max(input_filt_rare_ITS$map_loaded$Lichenized)+
           (max(input_filt_rare_ITS$map_loaded$Lichenized)-
              min(input_filt_rare_ITS$map_loaded$Lichenized))/10)
pdf("InitialFigs/Taxa_Guilds_Lichenized.pdf", width = 7, height = 5)
ggplot(input_filt_rare_ITS$map_loaded, 
       aes(reorder(Site, Lichenized, mean), Lichenized)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(size = 2, alpha = 0.75, width = 0.2, height = 0, pch = 16) +
  geom_text(data = t, aes(Site, y, label = str_trim(.group)), 
            size = 3, color = "black") +
  labs(x = "Site", y = "Lichenized relative abundance") +
  theme_bw() +
  theme(legend.position = "right",
        axis.title = element_text(size = 12),
        axis.text.y = element_text(size = 10),
        axis.text.x = element_text(size = 10, angle = 45, hjust = 1),
        strip.text = element_text(size = 14),
        strip.background = element_rect(fill = "white"))
dev.off()
sum(input_filt_rare_ITS$map_loaded$Lichenized > 0) # Present in 41 samples
sum_Lichenized <- input_filt_rare_ITS$map_loaded %>%
  group_by(Site) %>%
  summarise(Lichenized = mean(Lichenized)) 
sum(sum_Lichenized$Lichenized > 0) # Present in 9 sites.



#### _42 OTU ####
# Check the 42 OTUs that were positively associated with Sclerotinia resistance
# From Pogoda et al. 2023
# Test and plot
# Used negative binomial
# # Other options:
# # Negative binomial regression
# nb1 <- glm.nb(OTU_1191 ~ Site*Entry, data = df42, link = log)
# summary(nb1)
# Anova(nb1)
# nb1
# null <- glm.nb(OTU_1191 ~ 1, data = df42, link = log)
# anova(null, nb1)
# 
# # Zero-inflated Beta regression
# b1 <- gamlss(df42$OTU_1191 ~ df42$Site * df42$Entry,  family = BEZI, trace = F)
# summary(b1)
# plot(b1)
# emmeans(b1, "Environment", type = "response")
# pairs(.Last.value)



#### Import
# Note: Pogoda seqs are 200 bp. GxE seqs are longer, usually ~253.
pogoda_42 <- read.delim("data/OTUs_42_Sclero.txt", header = F)
pogoda_repset <- readFasta("data/Pogoda_repset.fasta") %>%
  filter(Header %in% pogoda_42$V1)
gxe_repset <- readFasta("data/repset_16S_filt.fasta")
asv_42 <- readFasta("data/repset_16S_filt.fasta") %>%
  filter(grepl(paste(pogoda_repset$Sequence, collapse = "|"), Sequence)) %>%
  separate(Header, into = c("ASV_ID", "taxonomy"), sep = " ", remove = T) %>%
  left_join(., otu_asv_map_16S, by = "ASV_ID")
# 48 selected. Multiple ASV sequences contain the same OTU. How many unique OTU?
length(unique(asv_42$OTU_ID)) # 37

# Which OTUs were not found? Get those and BLAST them. There might be a close (>97%) but not exact match
# It looks like only 3 OTUs were not found.
gxe_repset$Sequence200 <- substr(gxe_repset$Sequence, start = 1, stop = 200)
asv_42$Sequence200 <- substr(asv_42$Sequence, start = 1, stop = 200)
pogoda_repset_toBLAST <- pogoda_repset %>%
  left_join(., asv_42, by = c("Sequence" = "Sequence200"))
length(unique(pogoda_repset_toBLAST$Header)) # 42 Pogoda OTUs
length(unique(pogoda_repset_toBLAST$OTU_ID)) # 38 GxE OTUs
asv_42 <- asv_42 %>%
  left_join(., pogoda_repset_toBLAST, by = "ASV_ID")
length(unique(asv_42$Header)) # 39
length(unique(asv_42$ASV_ID)) # 48
# OTU 697 = 99.50% ID to ASV_2266
# OTU 991 = 99.00% ID to ASV_2024
# OTU 1022 = 99.00% ID to ASV_1640. ASV_1640 was already in an OTU that was kept. So we gained 2.
asv_IDs_42 <- c(asv_42$ASV_ID, "ASV_2266", "ASV_2024", "ASV_1640")
# Rerun the above code but adjust the filter to include these
# We need to filter by asv_42$OTU_ID
# Note: All 42 OTU sequences were found but they are contained in 39 OTUs in the GxE dataset
# This is because the GxE dataset had longer reads
# The 97% similarity threshold could be met over 253 bp when it wasn't over only 200 bp

asv_42 <- readFasta("data/repset_16S_filt.fasta") %>%
  separate(Header, into = c("ASV_ID", "taxonomy"), sep = " ", remove = T) %>%
  filter(ASV_ID %in% asv_IDs_42) %>%
  left_join(., otu_asv_map_16S, by = "ASV_ID")
length(unique(asv_42$OTU_ID)) # 39

ts42 <- tax_sum_OTU_16S %>%
  filter(rownames(.) %in% asv_42$OTU_ID) %>%
  t() %>%
  as.data.frame() %>%
  rownames_to_column(var = "sampleID")

df42 <- input_filt_rare_16S$map_loaded %>%
  dplyr::select(sampleID, Site, Entry) %>%
  left_join(., ts42, by = "sampleID")

top4_otus <- c("OTU_2454", "OTU_883", "OTU_2972", "OTU_9410")

# Site prevalence: Plot %OTU sample prev by site, also %OTU prev by site
prev_42 <- df42 %>%
  pivot_longer(cols = c(4:42), names_to = "OTU") %>%
  group_by(Site, OTU) %>%
  summarise(site_n = n(),
            otu_n = sum(value > 0)) %>%
  mutate(otu_per = otu_n/site_n*100)

prev_42_site <- prev_42 %>%
  group_by(Site) %>%
  summarise(n_otu = sum(otu_n > 0)) %>%
  mutate(per_otu = n_otu/39*100)

strip_site <- strip_themed(background_x = elem_list_rect(fill = c(rep("yellow", 4), 
                                                                  rep("white", 35))))
pdf("InitialFigs/Taxa_42_SitePrevalence.pdf", width = 7, height = 7)
ggplot(prev_42, aes(Site, otu_per)) +
  geom_bar(stat = "identity") +
  facet_wrap2(~ OTU, ncol = 4, strip = strip_site) +
  labs(x = "Site",
       y = "% Prevalence") +
  theme_bw() +
  theme(legend.position = "none",
        axis.text.y = element_text(size = 6),
        axis.text.x = element_text(size = 6, angle = 45, hjust = 1),
        strip.text = element_text(size = 6),
        panel.spacing.x = unit(c(0.5, 0.5, 0.5), "cm"))
dev.off()

pdf("InitialFigs/Taxa_42_OTUSitePrevalence.pdf", width = 7, height = 5)
figs9 <- ggplot(prev_42_site, aes(Site, per_otu)) +
  geom_bar(stat = "identity") +
  scale_y_continuous(limits = c(0, 100),
                     breaks = c(0, 25, 50, 75, 100)) +
  labs(x = "Site",
       y = "% 39 OTU Prevalence") +
  theme_bw() +
  theme(legend.position = "none",
        axis.text.y = element_text(size = 8),
        axis.text.x = element_text(size = 8, angle = 45, hjust = 1))
figs9
dev.off()

pdf("FinalFigs/FigureS9.pdf", width = 7, height = 5)
figs9
dev.off()

            

#### Test
results42 <- as.data.frame(matrix(NA, nrow = ncol(df42), ncol = 12)) %>%
  set_names(c("OTU", "LeveneS", "LeveneG", "ShapiroAOV", "ShapiroNB1", "ShapiroNB2",
              "SitePaov", "GenotypePaov", "IntPaov", "SitePnb", "GenotypePnb",
              "GenotypePnb_Car"))
for (i in 4:ncol(df42)) {
  # OTU name
  results42$OTU[i] <- names(df42)[i]
  
  # Levene Test
  l1 <- leveneTest(df42[,i] ~ df42$Site)
  l2 <- leveneTest(df42[,i] ~ df42$Entry)
  results42$LeveneS[i] <- l1$`Pr(>F)`[1]
  results42$LeveneG[i] <- l2$`Pr(>F)`[1]
  
  # Make data frame with the predictors and each OTU
  df <- df42 %>%
    dplyr::select(Site, Entry, i) %>%
    set_names(c("Site", "Entry", "OTU"))
  df_car <- df %>%
    filter(Site == "Carrington")
  
  # Models
  m <- aov(OTU ~ Site * Entry, data = df)
  nb1 <- glmm.nb(OTU ~ Site, random = ~ 1 | Entry, data = df)
  nb2 <- glmm.nb(OTU ~ Entry, random = ~ 1 | Site, data = df)
  nb3 <- MASS::glm.nb(OTU ~ Entry, data = df_car) # Not working but at least estimates theta
  nb4 <- glm(OTU ~ Entry, data = df_car, family = negative.binomial(nb3$theta))
  c1 <- Anova(m, type = "III", singular.ok = TRUE)
  c2 <- Anova(nb1)
  c3 <- Anova(nb2)
  c4 <- Anova(nb4)
  
  results42$SitePaov[i] <- c1$`Pr(>F)`[1]
  results42$GenotypePaov[i] <- c1$`Pr(>F)`[2]
  results42$IntPaov[i] <- c1$`Pr(>F)`[3]
  results42$SitePnb[i] <- c2$`Pr(>Chisq)`[1]
  results42$GenotypePnb[i] <- c3$`Pr(>Chisq)`[1]
  results42$GenotypePnb_Car[i] <- c4$`Pr(>Chisq)`[1]
  
  # Shapiro Test
  s1 <- shapiro.test(m$residuals)
  s2 <- shapiro.test(nb1$residuals)
  s3 <- shapiro.test(nb2$residuals)
  results42$ShapiroAOV[i] <- s1$p.value
  results42$ShapiroNB1[i] <- s2$p.value
  results42$ShapiroNB2[i] <- s3$p.value
}

results42 <- results42 %>%
  drop_na(OTU) %>%
  mutate(SitePFDRaov = p.adjust(SitePaov, method = "fdr"),
         GenotypePFDRaov = p.adjust(GenotypePaov, method = "fdr"),
         IntPFDRaov = p.adjust(IntPaov, method = "fdr"),
         SitePFDRnb = p.adjust(SitePnb, method = "fdr"),
         GenotypePFDRnb = p.adjust(GenotypePnb, method = "fdr"),
         GenotypePFDRnb_Car = p.adjust(GenotypePnb_Car, method = "fdr")) %>%
  mutate(SymbolSite = ifelse(SitePnb <= 0.001,
                             "***",
                             ifelse(SitePnb > 0.001 & SitePnb <= 0.01,
                                    "**",
                                    ifelse(SitePnb > 0.01 & SitePnb <= 0.05,
                                           "*",
                                           "N.S.")))) %>%
  mutate(SymbolGeno = ifelse(GenotypePnb <= 0.001,
                             "***",
                             ifelse(GenotypePnb > 0.001 & GenotypePnb <= 0.01,
                                    "**",
                                    ifelse(GenotypePnb > 0.01 & GenotypePnb <= 0.05,
                                           "*",
                                           "N.S.")))) %>%
  mutate(SymbolGeno_Car = ifelse(GenotypePnb_Car <= 0.001,
                             "***",
                             ifelse(GenotypePnb_Car > 0.001 & GenotypePnb_Car <= 0.01,
                                    "**",
                                    ifelse(GenotypePnb_Car > 0.01 & GenotypePnb_Car <= 0.05,
                                           "*",
                                           "N.S."))))
sigSite <- results42 %>%
  filter(SitePFDRaov < 0.05) # 32 p, 32 pFDR
sigGeno <- results42 %>%
  filter(GenotypePaov < 0.05) # 0 p, 0 pFDR
sigInt <- results42 %>%
  filter(IntPFDRaov < 0.05) # 13 p, 9 pFDR
sigSiteNB <- results42 %>%
  filter(SitePFDRnb < 0.05) # 39 p, 39 pFDR
sigGenoNB <- results42 %>%
  filter(GenotypePnb < 0.05) # 11 p, 5 pFDR
sigGenoNB_Car <- results42 %>%
  filter(GenotypePnb_Car < 0.05) # 12 p, 7 pFDR



#### Plot
## Genotype
sigGeno_otus <- sigGenoNB$OTU[sigGenoNB$OTU %notin% top4_otus]
otherGeno_otus <- unique(asv_42$OTU_ID)[unique(asv_42$OTU_ID) %notin% c(top4_otus, sigGeno_otus)]
df42_long_geno <- df42 %>%
  pivot_longer(cols = c(4:ncol(.))) %>%
  mutate(top4 = ifelse(name %in% top4_otus,
                       "top4",
                       "not")) %>%
  left_join(., results42, by = c("name" = "OTU")) %>%
  left_join(., otu_tax_16S, by = c("name" = "OTU_ID")) %>%
  mutate(name = factor(name,
                       levels = c(top4_otus, sigGeno_otus, otherGeno_otus)))
sigLab_geno <- df42_long_geno %>%
  group_by(name) %>%
  summarise(maxy = max(value)) %>%
  mutate(y = maxy + maxy/2) %>%
  left_join(., results42, by = c("name" = "OTU")) %>%
  left_join(., otu_tax_16S, by = c("name" = "OTU_ID")) %>%
  mutate(name = factor(name,
                       levels = c(top4_otus, sigGeno_otus, otherGeno_otus))) %>%
  arrange(name = desc(levels(df42_long_geno$name)))
strip_geno <- strip_themed(background_x = elem_list_rect(fill = c(rep("yellow", 4), 
                                                                  rep("lightblue", 11),
                                                                  rep("white", 24))))
facet_labs_geno <- paste(sigLab_geno$taxonomy2, sigLab_geno$name, sep = ": ")
names(facet_labs_geno) <- sigLab_geno$name
pdf("InitialFigs/Taxa_42_Genotype.pdf", width = 7, height = 7)
figs10 <- ggplot(df42_long_geno, aes(Entry, value*100)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(size = 0.25, width = 0.1) +
  geom_text(data = subset(sigLab_geno, SymbolGeno != "N.S."),
            aes(x = 5.5, y = y*100, label = SymbolGeno),
            inherit.aes = F, size = 8) +
  geom_text(data = subset(sigLab_geno, SymbolGeno == "N.S."), 
            aes(x = 5.5, y = (y*100)+(y*100/10), label = SymbolGeno),
            inherit.aes = F, size = 2.5) +
  scale_colour_viridis_d() +
  facet_wrap2(~ name, ncol = 4, scales = "free_y", strip = strip_geno,
              labeller = labeller(name = facet_labs_geno)) +
  scale_y_continuous(expand = expand_scale(mult = c(0.1, 0.6), 
                                           add = c(0, 0))) +
  labs(x = "Genotype",
       y = "% Relative abundance") +
  theme_bw() +
  theme(legend.position = "none",
        axis.text.y = element_text(size = 6),
        axis.text.x = element_text(size = 8, angle = 45, hjust = 1),
        strip.text = element_text(size = 6))
figs10
dev.off()

pdf("FinalFigs/FigureS10.pdf", width = 7, height = 7)
figs10
dev.off()


## Site
sigSite_otus <- sigSiteNB$OTU[sigSiteNB$OTU %notin% top4_otus]
otherSite_otus <- unique(asv_42$OTU_ID)[unique(asv_42$OTU_ID) %notin% c(top4_otus)]
df42_long_site <- df42 %>%
  pivot_longer(cols = c(4:ncol(.))) %>%
  mutate(top4 = ifelse(name %in% top4_otus,
                       "top4",
                       "not")) %>%
  left_join(., results42, by = c("name" = "OTU")) %>%
  left_join(., otu_tax_16S, by = c("name" = "OTU_ID")) %>%
  mutate(name = factor(name,
                       levels = c(top4_otus, otherSite_otus)))
sigLab_site <- df42_long_site %>%
  group_by(name) %>%
  summarise(maxy = max(value)) %>%
  mutate(y = maxy + maxy/10) %>%
  left_join(., results42, by = c("name" = "OTU")) %>%
  left_join(., otu_tax_16S, by = c("name" = "OTU_ID")) %>%
  mutate(name = factor(name,
                       levels = c(top4_otus, otherSite_otus))) %>%
  arrange(name = desc(levels(df42_long_site$name)))
strip_site <- strip_themed(background_x = elem_list_rect(fill = c(rep("yellow", 4), 
                                                                  rep("lightblue", 35))))
facet_labs_site <- paste(sigLab_site$taxonomy2, sigLab_site$name, sep = ": ")
names(facet_labs_site) <- sigLab_site$name
pdf("InitialFigs/Taxa_42_Site.pdf", width = 7, height = 7)
fig6 <- ggplot(df42_long_site, aes(Site, value*100)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(size = 0.25, width = 0.1) +
  geom_text(data = sigLab_site,
            aes(x = 8, y = y*100, label = SymbolSite),
            inherit.aes = F, size = 8) +
  scale_colour_viridis_d() +
  facet_wrap2(~ name, ncol = 4, scales = "free_y", strip = strip_site,
              labeller = labeller(name = facet_labs_site)) +
  scale_y_continuous(expand = expand_scale(mult = c(0.1, 0.6), 
                                           add = c(0, 0))) +
  labs(x = "Site",
       y = "% Relative abundance") +
  theme_bw() +
  theme(legend.position = "none",
        axis.text.y = element_text(size = 6),
        axis.text.x = element_text(size = 6, angle = 45, hjust = 1),
        strip.text = element_text(size = 6))
fig6
dev.off()
pdf("FinalFigs/Figure6.pdf", width = 7, height = 7)
fig6
dev.off()


## Carrington
# Check Carrington - this is where these OTUs were originally found.
# These OTUs were affected by genotype when there were 95 genotypes, and correlated with Sclero resistance
# What about in this dataset? Affected by 10 genotypes?
sigGeno_Car_otus <- sigGenoNB_Car$OTU[sigGenoNB_Car$OTU %notin% top4_otus]
otherGeno_Car_otus <- unique(asv_42$OTU_ID)[unique(asv_42$OTU_ID) %notin% c(top4_otus, sigGeno_Car_otus)]
df42_long_geno_Car <- df42 %>%
  pivot_longer(cols = c(4:ncol(.))) %>%
  mutate(top4 = ifelse(name %in% top4_otus,
                       "top4",
                       "not")) %>%
  left_join(., results42, by = c("name" = "OTU")) %>%
  left_join(., otu_tax_16S, by = c("name" = "OTU_ID")) %>%
  mutate(name = factor(name,
                       levels = c(top4_otus, sigGeno_Car_otus[sigGeno_Car_otus %notin% top4_otus], otherGeno_Car_otus)))
sigLab_geno_Car <- df42_long_geno_Car %>%
  group_by(name) %>%
  summarise(maxy = max(value)) %>%
  mutate(y = maxy + maxy/10) %>%
  left_join(., results42, by = c("name" = "OTU")) %>%
  left_join(., otu_tax_16S, by = c("name" = "OTU_ID")) %>%
  mutate(name = factor(name,
      levels = c(top4_otus, sigGeno_Car_otus[sigGeno_Car_otus %notin% top4_otus], otherGeno_Car_otus))) %>%
  arrange(name = desc(levels(df42_long_geno_Car$name)))
strip_geno_Car <- strip_themed(background_x = elem_list_rect(fill = c(rep("yellow", 4), 
                                                                      rep("lightblue", 11),
                                                                      rep("white", 24))))
facet_labs_geno_Car <- paste(sigLab_geno_Car$taxonomy2, sigLab_geno_Car$name, sep = ": ")
names(facet_labs_geno_Car) <- sigLab_geno_Car$name
pdf("InitialFigs/Taxa_42_Genotype_Carrington.pdf", width = 7, height = 7)
ggplot(df42_long_geno_Car, aes(Entry, value*100)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(size = 0.25, width = 0.1) +
  geom_text(data = subset(sigLab_geno_Car, SymbolGeno_Car != "N.S."),
            aes(x = 5.5, y = y*100, label = SymbolGeno_Car),
            inherit.aes = F, size = 8) +
  geom_text(data = subset(sigLab_geno_Car, SymbolGeno_Car == "N.S."), 
            aes(x = 5.5, y = y*100, label = SymbolGeno_Car),
            inherit.aes = F, size = 3) +
  scale_colour_viridis_d() +
  facet_wrap2(~ name, ncol = 4, scales = "free_y", strip = strip_geno,
              labeller = labeller(name = facet_labs_geno_Car)) +
  scale_y_continuous(expand = expand_scale(mult = c(0.1, 0.5), 
                                           add = c(0, 0))) +
  labs(x = "Genotype",
       y = "% Relative abundance") +
  theme_bw() +
  theme(legend.position = "none",
        axis.text.y = element_text(size = 6),
        axis.text.x = element_text(size = 8, angle = 45, hjust = 1),
        strip.text = element_text(size = 6))
dev.off()

# Shared?
sigGeno_otus
sum(sigGeno_Car_otus %in% sigGeno_otus)



#### _6 OTU ####
# 13 strains were shown in the laboratory to inhibit Sclerotinia
# These correspond to 6 unique OTUs
# We have the Pogoda IDs - match to this dataset
inhibitors <- c("OTU_10441", "OTU_12372", "OTU_978", "OTU_4850", "OTU_275", "OTU_1672")
# Any in the 42 OTUs?
sum(inhibitors %in% pogoda_42$V1) # 0
# NCBI Blast with my repset to get the ASV ID
in_ASVs <- c("ASV_13525", "ASV_10535", "ASV_3116", "ASV_15020", "ASV_806", "ASV_1866")

# Get OTU IDs for those 6
input_6 <- filter_taxa_from_input(input_filt_rare_16S,
                                  taxa_IDs_to_keep = in_ASVs)
View(input_6$taxonomy_loaded)
facet_names_6 <- c("OTU_1757" = "Pseudomonas\nsp.",
                   "OTU_1988" = "Rhodococcus\nsp.",
                   "OTU_2406" = "Bacillus\nsp.",
                   "OTU_3404" = "Variovorax\nparadoxus",
                   "OTU_3523" = "Paenibacillus\npolymyxa",
                   "OTU_9352" = "Sphingobacterium\nkitahiroshimense")
ts6 <- tax_sum_OTU_16S %>%
  filter(rownames(.) %in% input_6$taxonomy_loaded$OTU_ID) %>%
  t() %>%
  as.data.frame() %>%
  rownames_to_column(var = "sampleID")
df6 <- input_filt_rare_16S$map_loaded %>%
  dplyr::select(sampleID, Site, Rep, Entry) %>%
  left_join(., ts6, by = "sampleID")
prev_all <- df6 %>%
  pivot_longer(cols = c(4:9), names_to = "OTU") %>%
  group_by(OTU) %>%
  summarise(otu_n = sum(value > 0)) %>%
  mutate(otu_per = otu_n/587*100)
prev_6 <- df6 %>%
  pivot_longer(cols = c(4:9), names_to = "OTU") %>%
  group_by(Site, OTU) %>%
  summarise(site_n = n(),
            otu_n = sum(value > 0)) %>%
  mutate(otu_per = otu_n/site_n*100)
pdf("InitialFigs/Taxa_6_SitePrevalence.pdf", width = 7, height = 4)
ggplot(prev_6, aes(Site, otu_per)) +
  geom_bar(stat = "identity") +
  facet_wrap(~ OTU, ncol = 6, labeller = as_labeller(facet_names_6)) +
  labs(x = "Site",
       y = "% Prevalence") +
  scale_y_continuous(expand = c(0.01, 0.01)) +
  theme_bw() +
  theme(legend.position = "none",
        axis.text.y = element_text(size = 6),
        axis.text.x = element_text(size = 4, angle = 45, hjust = 1),
        strip.text = element_text(size = 6, face = "italic"),
        strip.background = element_rect(linewidth = 0.2, fill = "white"))
dev.off()

results6 <- as.data.frame(matrix(NA, nrow = ncol(df6), ncol = 13)) %>%
  set_names(c("OTU", "Heritability", "LeveneS", "LeveneG", "ShapiroAOV", "ShapiroNB1", "ShapiroNB2",
              "SitePaov", "GenotypePaov", "IntPaov", "SitePnb", "GenotypePnb",
              "GenotypePnb_Car"))
for (i in 5:ncol(df6)) {
  # OTU name
  results6$OTU[i] <- names(df6)[i]
  
  # Levene Test
  l1 <- leveneTest(df6[,i] ~ df6$Site)
  l2 <- leveneTest(df6[,i] ~ df6$Entry)
  results6$LeveneS[i] <- l1$`Pr(>F)`[1]
  results6$LeveneG[i] <- l2$`Pr(>F)`[1]
  
  # Make data frame with the predictors and each OTU
  df <- df6 %>%
    dplyr::select(Site, Rep, Entry, i) %>%
    set_names(c("Site", "Rep", "Entry", "OTU"))
  df_car <- df %>%
    filter(Site == "Carrington")
  
  # Models
  m <- aov(OTU ~ Site * Entry, data = df)
  nb1 <- glmm.nb(OTU ~ Site, random = ~ 1 | Entry, data = df)
  nb2 <- glmm.nb(OTU ~ Entry, random = ~ 1 | Site, data = df)
  nb2.1 <- glm.nb(OTU ~ Site + Rep + Entry, data = df)
  #nb2.2 <- glm.nb(OTU ~ Entry, data = df)
  nb3 <- MASS::glm.nb(OTU ~ Entry, data = df_car) # Not working but at least estimates theta
  nb4 <- glm(OTU ~ Entry, data = df_car, family = negative.binomial(nb3$theta))
  c1 <- Anova(m, type = "III", singular.ok = TRUE)
  c2 <- Anova(nb1, type = "II", test.statistic = "F")
  c3 <- car::Anova(nb2, type = "II", test.statistic = "F")
  c4 <- Anova(nb4)
  
  results6$SitePaov[i] <- c1$`Pr(>F)`[1]
  results6$GenotypePaov[i] <- c1$`Pr(>F)`[2]
  results6$IntPaov[i] <- c1$`Pr(>F)`[3]
  results6$SitePnb[i] <- c2$`Pr(>Chisq)`[1]
  results6$GenotypePnb[i] <- c3$`Pr(>Chisq)`[1]
  results6$GenotypePnb_Car[i] <- c4$`Pr(>Chisq)`[1]
  
  eta <- eta_sq_glm(nb2.1) %>%
    as.data.frame() %>%
    rownames_to_column(var = "variable") %>%
    set_names(c("variable", "eta"))
  
  results6$Heritability[i] <- eta$eta[3]
  
  # Shapiro Test
  s1 <- shapiro.test(m$residuals)
  s2 <- shapiro.test(nb1$residuals)
  s3 <- shapiro.test(nb2$residuals)
  results6$ShapiroAOV[i] <- s1$p.value
  results6$ShapiroNB1[i] <- s2$p.value
  results6$ShapiroNB2[i] <- s3$p.value
}
# Rerun manually for i = 9 and i = 10 to get around error
i = 9
i = 10

# Process
results6 <- results6 %>%
  drop_na(OTU) %>%
  mutate(SitePFDRaov = p.adjust(SitePaov, method = "fdr"),
         GenotypePFDRaov = p.adjust(GenotypePaov, method = "fdr"),
         IntPFDRaov = p.adjust(IntPaov, method = "fdr"),
         SitePFDRnb = p.adjust(SitePnb, method = "fdr"),
         GenotypePFDRnb = p.adjust(GenotypePnb, method = "fdr"),
         GenotypePFDRnb_Car = p.adjust(GenotypePnb_Car, method = "fdr")) %>%
  mutate(SymbolSite = ifelse(SitePnb <= 0.001,
                             "***",
                             ifelse(SitePnb > 0.001 & SitePnb <= 0.01,
                                    "**",
                                    ifelse(SitePnb > 0.01 & SitePnb <= 0.05,
                                           "*",
                                           "N.S.")))) %>%
  mutate(SymbolGeno = ifelse(GenotypePnb <= 0.001,
                             "***",
                             ifelse(GenotypePnb > 0.001 & GenotypePnb <= 0.01,
                                    "**",
                                    ifelse(GenotypePnb > 0.01 & GenotypePnb <= 0.05,
                                           "*",
                                           "N.S.")))) %>%
  mutate(SymbolGeno_Car = ifelse(GenotypePnb_Car <= 0.001,
                                 "***",
                                 ifelse(GenotypePnb_Car > 0.001 & GenotypePnb_Car <= 0.01,
                                        "**",
                                        ifelse(GenotypePnb_Car > 0.01 & GenotypePnb_Car <= 0.05,
                                               "*",
                                               "N.S."))))
sigSite <- results6 %>%
  filter(SitePFDRaov < 0.05) # 1 p, 1 pFDR
sigGeno <- results6 %>%
  filter(GenotypePaov < 0.05) # 0 p, 0 pFDR
sigInt <- results6 %>%
  filter(IntPFDRaov < 0.05) # 0 p, 0 pFDR
sigSiteNB <- results6 %>%
  filter(SitePFDRnb < 0.05) # 6 p, 6 pFDR
sigGenoNB <- results6 %>%
  filter(GenotypePnb < 0.05) # 3 p, 2 pFDR
sigGenoNB_Car <- results6 %>%
  filter(GenotypePnb_Car < 0.05) # 1 p, 0 pFDR

## Genotype
sigGeno_otus <- sigGenoNB$OTU
otherGeno_otus <- results6$OTU[results6$OTU %notin% sigGenoNB$OTU]
df6_long_geno <- df6 %>%
  pivot_longer(cols = c(4:ncol(.))) %>%
  mutate(top4 = ifelse(name %in% top4_otus,
                       "top4",
                       "not")) %>%
  left_join(., results6, by = c("name" = "OTU")) %>%
  left_join(., otu_tax_16S, by = c("name" = "OTU_ID")) %>%
  mutate(name = factor(name,
                       levels = c(sigGeno_otus, otherGeno_otus))) %>%
  mutate(name = as.character(name))
sigLab_geno <- df6_long_geno %>%
  group_by(name) %>%
  summarise(maxy = max(value)) %>%
  mutate(y = maxy + maxy/10) %>%
  left_join(., results6, by = c("name" = "OTU")) %>%
  left_join(., otu_tax_16S, by = c("name" = "OTU_ID")) %>%
  mutate(name = factor(name,
                       levels = c(sigGeno_otus, otherGeno_otus))) %>%
  arrange(name = desc(levels(df6_long_geno$name))) %>%
  mutate(name = as.character(name))
pdf("InitialFigs/Taxa_6_Genotype.pdf", width = 7, height = 4)
ggplot(df6_long_geno, aes(Entry, value*100)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(size = 0.25, width = 0.15, height = 0) +
  geom_text(data = subset(sigLab_geno, SymbolGeno != "N.S."),
            aes(x = 5.5, y = y*100, label = SymbolGeno),
            inherit.aes = F, size = 8) +
  geom_text(data = subset(sigLab_geno, SymbolGeno == "N.S."), 
            aes(x = 5.5, y = y*100, label = SymbolGeno),
            inherit.aes = F, size = 3) +
  facet_wrap(~ name, ncol = 6, scales = "free_y", labeller = as_labeller(facet_names_6)) +
  scale_y_continuous(expand = expand_scale(mult = c(0.01, 0.15), 
                                           add = c(0, 0))) +
  labs(x = "Genotype",
       y = "% Relative abundance") +
  theme_bw() +
  theme(legend.position = "none",
        axis.text.y = element_text(size = 6),
        axis.text.x = element_text(size = 4, angle = 45, hjust = 1),
        strip.text = element_text(size = 6, face = "italic"),
        strip.background = element_rect(linewidth = 0.2, fill = "white"))
dev.off()



## Site
sigSite_otus <- sigSiteNB$OTU
df6_long_site <- df6 %>%
  pivot_longer(cols = c(4:ncol(.))) %>%
  mutate(top4 = ifelse(name %in% top4_otus,
                       "top4",
                       "not")) %>%
  left_join(., results6, by = c("name" = "OTU")) %>%
  left_join(., otu_tax_16S, by = c("name" = "OTU_ID")) %>%
  mutate(name = factor(name,
                       levels = c(sigSite_otus))) %>%
  mutate(name = as.character(name))
sigLab_site <- df6_long_site %>%
  group_by(name) %>%
  summarise(maxy = max(value)) %>%
  mutate(y = maxy + maxy/10) %>%
  left_join(., results6, by = c("name" = "OTU")) %>%
  left_join(., otu_tax_16S, by = c("name" = "OTU_ID")) %>%
  mutate(name = factor(name,
                       levels = c(sigSite_otus))) %>%
  arrange(name = desc(levels(df6_long_site$name))) %>%
  mutate(name = as.character(name))
pdf("InitialFigs/Taxa_6_Site.pdf", width = 7, height = 4)
ggplot(df6_long_site, aes(Site, value*100)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(size = 0.25, width = 0.15, height = 0) +
  geom_text(data = sigLab_site,
            aes(x = 7.5, y = y*100, label = SymbolSite),
            inherit.aes = F, size = 8) +
  facet_wrap(~ name, ncol = 6, scales = "free_y", labeller = as_labeller(facet_names_6)) +
  scale_y_continuous(expand = expand_scale(mult = c(0.01, 0.15), 
                                           add = c(0, 0))) +
  labs(x = "Site",
       y = "% Relative abundance") +
  theme_bw() +
  theme(legend.position = "none",
        axis.text.y = element_text(size = 6),
        axis.text.x = element_text(size = 4, angle = 45, hjust = 1),
        strip.text = element_text(size = 6, face = "italic"),
        strip.background = element_rect(linewidth = 0.2, fill = "white"))
dev.off()



#### _PiCRUST2 ####
# Check abundances of genes of interest
# Schmidt et al. 2019: nifH, amoA [bacterial and archaeal], nirK, nrfA, and nosZ
# nifH = K02588
# amoA = K10944
# nirK = K00368
# nirS = K15864
# nosZ = K00376
# nrfA = K03385
# narG K00370 nitrate reductase
# napA K02567 nitrate reductase
# norB K04561 nitric oxide reductase

# Others:
# Kechid et al. 2013
# These play an essential role in plant growth promotion by the rhizospheric bacterium STM196
# NRT2.5 = K02575 nitrate transporter2.5
# NRT2.6 = K02575 high affinity nitrate transporter 2.6

# Duan et al. 2013
# Genes potentially involved in plant growth promotion such as indole-3-acetic acid (IAA) biosynthesis, trehalose production, siderophore production, acetoin synthesis, and phosphate solubilization were determined.
# acdS = K01505 K ACC deaminase
  # acdR = No KO; ACC deaminase upstream regulation
  # pvdF = No KO; K pyoverdine synthesis
# pvdYII = K03896 hydroxyornithine acetylase (siderophore)
# tms1 = K00466 IAA biosynthesis (tryptophan 2'-monooxygenase)
# tms2 = K21801 IAA biosynthesis	(indoleacetamide hydrolase)
# K01721, K20807 IAA biosynthesis (nitrile hydratase)
  # AMI1 = No KO; indole-3-acetamide amidohydrolase
# iaaH = K21801 indoleacetamide hydrolase
# K01652, K01653, K11258 acetolactate synthase
  # No KO; zinc-containing alcohol dehydrogenase
# ubiC = K03181 chorismate lyase encoded 
  # treS = no KO trehalose synthase
# treT = K13057 trehalose synthase
# treY = K06044 maltooligosyltrehalose synthase
# treZ = K01236 maltooligosyltrehalose trehalohydrolase
# OtsA = K00697 other trehalose pathway
# OtsB = K01087 other trehalose pathway

# P solubilization (Rawat et al. 2021)
# pqqABCD = K06135, K06136, K06137, K06138 pyrroloquinoline quinine synthesis
# gcd = K00117 glucose dehydrogenase

# Import PiCRUST2 KO annotation table
ko <- read.delim("data/pred_metagenome_unstrat.tsv") %>%
  column_to_rownames(var = "function.")

# Format
ko_t <- ko %>%
  t() %>%
  as.data.frame() %>%
  mutate(rn = substr(rownames(.), start = 2, stop = nchar(rownames(.)))) %>%
  mutate(rn = gsub("\\.", "\\-", rn)) %>%
  rownames_to_column(var = "rowname") %>%
  dplyr::select(-rowname) %>%
  column_to_rownames(var = "rn")

# Check order and reorder if necessary with match
sum(rownames(input_filt_rare_16S$map_loaded) != rownames(ko_t)) # Good

# Make df
meta <- input_filt_rare_16S$map_loaded
meta <- meta %>%
  dplyr::select(Site, Entry) %>%
  mutate(amoA = ko_t$K10944,
         nifH = ko_t$K02588,
         narG = ko_t$K00370,
         napA = ko_t$K02567,
         nirK = ko_t$K00368,
         nirS = ko_t$K15864,
         nrfA = ko_t$K03385,
         norB = ko_t$K04561,
         nosZ = ko_t$K00376,
         acdS = ko_t$K01505,
         pvdYII = ko_t$K03896,
         tms1 = ko_t$K00466,
         #tms2 = ko_t$K21801, # Not found
         nhyd1 = ko_t$K01721,
         #nhyd2 = ko_t$K20807, # Not found
         #iaaH = ko_t$K21801, # Not found
         acsyn1 = ko_t$K01652, 
         acsyn2 = ko_t$K01653, 
         acsyn3 = ko_t$K11258,
         ubiC = ko_t$K03181,
         treT = ko_t$K13057,
         #treY = ko_t$K0604, # Not found
         treZ = ko_t$K01236,
         OtsA = ko_t$K00697,
         OtsB = ko_t$K01087,
         pqqA = ko_t$K06135,
         pqqB = ko_t$K06136,
         pqqC = ko_t$K06137,
         pqqD = ko_t$K06138,
         gcd = ko_t$K00117) %>%
  dplyr::select(-Site, -Entry, Site, Entry)

# Note: Plot 1-9 "Nitrogen Cycling Genes. 10-24 as "Other PGP Genes", 25-29 as "P-solubilization Genes"

# Stats and plot - nifH demo
leveneTest(meta$nifH ~ meta$Site) # Not Homogeneous
leveneTest(meta$nifH ~ meta$Entry) # Homogeneous
m <- aov(nifH ~ Site * Entry, data = meta)
summary(m) # All
Anova(m, type = "III", singular.ok = TRUE) # Site***, Int*
hist(m$residuals)
shapiro.test(m$residuals)
plot(m$fitted.values, m$residuals)
m <- aov(nifH ~ Entry, data = meta)
summary(m)
m <- aov(nifH ~ Site, data = meta)
summary(m)
TukeyHSD(m)
t <- emmeans(object = m, specs = "Site") %>%
  cld(object = ., adjust = "Tukey", Letters = letters, alpha = 0.05) %>%
  mutate(name = "nifH",
         y = max(meta$nifH)+(max(meta$nifH)-min(meta$nifH))/10)
pdf("InitialFigs/Function_nifH.pdf", width = 7, height = 5)
ggplot(meta, aes(reorder(Site, nifH, mean), nifH)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(size = 2, alpha = 0.75, width = 0.2, pch = 16) +
  geom_text(data = t, aes(Site, y, label = str_trim(.group)), 
            size = 3, color = "black") +
  labs(x = "Site", y = "nifH abundance") +
  theme_bw() +
  theme(legend.position = "right",
        axis.title = element_text(size = 12),
        axis.text.y = element_text(size = 10),
        axis.text.x = element_text(size = 10, angle = 45, hjust = 1),
        strip.text = element_text(size = 14),
        strip.background = element_rect(fill = "white"))
dev.off()

# Nitrogen Cycling
# Stats - for loop
m <- list()
t <- list()
for (i in 1:9) {
  df <- meta %>%
    dplyr::select(Site, Entry, i) %>%
    set_names(c("Site", "Entry", "Gene"))
  m[[i]] <- aov(Gene ~ Site, data = df)
  t[[i]] <- emmeans(object = m[[i]], specs = "Site") %>%
    cld(object = ., adjust = "Tukey", Letters = letters, alpha = 0.05) %>%
    mutate(name = names(meta)[i],
           y = max(meta[,i]) + max(meta[,i])/10)
}
t_merge <- t[[1]]
for (i in 2:length(t)) {
  t_merge <- rbind(t_merge, t[[i]])
}
t_merge <- t_merge %>%
  mutate(name = factor(name,
                       levels = c("amoA", "nifH", "narG", "napA", "nirK", "nirS", 
                                  "nrfA", "norB", "nosZ"))) %>%
  mutate(func = recode_factor(name,
                              amoA = "NH3 ox.", 
                              nifH = "N fix.", 
                              narG = "NO3 red.", 
                              napA = "NO3 red.", 
                              nirK = "NO2 red.", 
                              nirS = "NO2 red.", 
                              nrfA = "NO2 red.", 
                              norB = "NO red.", 
                              nosZ = "N2O red.")) %>%
  mutate(func = factor(func,
                       levels = c("N fix.", "NH3 ox.", 
                                  "NO3 red.", "NO2 red.", "NO red.", "N2O red.")))
func <- t_merge %>%
  group_by(name) %>%
  slice_head(n = 1) %>%
  dplyr::select(name, func)

# Plot - multipanel
gene_plot <- meta %>%
  pivot_longer(cols = c("amoA", "nifH", "narG", "napA", "nirK", "nirS", 
                        "nrfA", "norB", "nosZ")) %>%
  mutate(name = factor(name,
                       levels = c("amoA", "nifH", "narG", "napA", "nirK", "nirS", 
                                  "nrfA", "norB", "nosZ"))) %>%
  dplyr::select(Site, Entry, name, value) %>%
  left_join(., func, by = "name") %>%
  mutate(func = factor(func,
                       levels = c("N fix.", "NH3 ox.", 
                                  "NO3 red.", "NO2 red.", "NO red.", "N2O red."))) %>%
  mutate(pc = "Nitrogen Cycling Genes")
pdf("InitialFigs/Function_Ncyc.pdf", width = 6, height = 8)
ggplot(gene_plot, aes(Site, value)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(size = 1, alpha = 0.5, width = 0.2, pch = 16) +
  geom_text(data = t_merge, aes(Site, y, label = str_trim(.group)), 
            size = 3, color = "black") +
  scale_y_continuous(expand = c(0.1, 0.1)) +
  labs(x = NULL,
       y = "PICRUSt2 abundance") +
  facet_nested(func + name ~ pc, scale = "free") +
  theme_bw() +
  theme(legend.position = "none",
        axis.title.y = element_text(size = 12),
        axis.text.x = element_text(size = 10, angle = 45, hjust = 1),
        axis.text.y = element_text(size = 8),
        strip.text.x = element_text(size = 12),
        strip.text.y = element_text(size = 6.5, face = "italic"),
        strip.background = element_rect(fill = "white"))
dev.off()

# Other PGP Genes
# Stats - for loop
meta_PGP <- meta %>%
  dplyr::select(10:21, Site, Entry)
m <- list()
t <- list()
for (i in 1:12) {
  df <- meta_PGP %>%
    dplyr::select(Site, Entry, i) %>%
    set_names(c("Site", "Entry", "Gene"))
  m[[i]] <- aov(Gene ~ Site, data = df)
  t[[i]] <- emmeans(object = m[[i]], specs = "Site") %>%
    cld(object = ., adjust = "Tukey", Letters = letters, alpha = 0.05) %>%
    mutate(name = names(meta_PGP)[i],
           y = max(meta_PGP[,i]) + max(meta_PGP[,i])/10)
}
t_merge <- t[[1]]
for (i in 2:length(t)) {
  t_merge <- rbind(t_merge, t[[i]])
}
t_merge <- t_merge %>%
  mutate(name = factor(name,
                       levels = names(meta_PGP)[1:12])) %>%
  mutate(func = recode_factor(name,
                              acdS = "ACC", 
                              pvdYII = "Sidero.", 
                              tms1 = "IAA", 
                              nhyd1 = "IAA", 
                              acsyn1 = "Acetoin", 
                              acsyn2 = "Acetoin", 
                              acsyn3 = "Acetoin", 
                              ubiC = "Antimic.", 
                              treT = "Trehalose",
                              treZ = "Trehalose",
                              OtsA = "Trehalose",
                              OtsB = "Trehalose")) %>%
  mutate(func = factor(func,
                       levels = c("ACC", "Sidero.", 
                                  "IAA", "Acetoin", "Antimic.", "Trehalose")))
func <- t_merge %>%
  group_by(name) %>%
  slice_head(n = 1) %>%
  dplyr::select(name, func)

# Plot - multipanel
gene_plot <- meta_PGP %>%
  pivot_longer(cols = names(meta_PGP)[1:12]) %>%
  mutate(name = factor(name,
                       levels = names(meta_PGP)[1:12])) %>%
  dplyr::select(Site, Entry, name, value) %>%
  left_join(., func, by = "name") %>%
  mutate(func = factor(func,
                       levels = c("ACC", "Sidero.", 
                                  "IAA", "Acetoin", "Antimic.", "Trehalose"))) %>%
  mutate(pc = "PGP Genes")
pdf("InitialFigs/Function_PGP.pdf", width = 6, height = 8)
ggplot(gene_plot, aes(Site, value)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(size = 1, alpha = 0.5, width = 0.2, pch = 16) +
  geom_text(data = t_merge, aes(Site, y, label = str_trim(.group)), 
            size = 3, color = "black") +
  scale_y_continuous(expand = c(0.2, 0.2)) +
  labs(x = NULL,
       y = "PICRUSt2 abundance") +
  facet_nested(func + name ~ pc, scale = "free") +
  theme_bw() +
  theme(legend.position = "none",
        axis.title.y = element_text(size = 12),
        axis.text.x = element_text(size = 10, angle = 45, hjust = 1),
        axis.text.y = element_text(size = 6),
        strip.text.x = element_text(size = 12),
        strip.text.y = element_text(size = 6.5, face = "italic"),
        strip.background = element_rect(fill = "white"))
dev.off()



# P solubilization
# Stats - for loop
meta_Psol <- meta %>%
  dplyr::select(22:26, Site, Entry)
m <- list()
t <- list()
for (i in 1:5) {
  df <- meta_Psol %>%
    dplyr::select(Site, Entry, i) %>%
    set_names(c("Site", "Entry", "Gene"))
  m[[i]] <- aov(Gene ~ Site, data = df)
  t[[i]] <- emmeans(object = m[[i]], specs = "Site") %>%
    cld(object = ., adjust = "Tukey", Letters = letters, alpha = 0.05) %>%
    mutate(name = names(meta_Psol)[i],
           y = max(meta_Psol[,i]) + max(meta_Psol[,i])/10)
}
t_merge <- t[[1]]
for (i in 2:length(t)) {
  t_merge <- rbind(t_merge, t[[i]])
}
t_merge <- t_merge %>%
  mutate(name = factor(name,
                       levels = names(meta_Psol)[1:5])) %>%
  mutate(func = recode_factor(name,
                              pqqA = "Pyrroloquinoline quinine", 
                              pqqB = "Pyrroloquinoline quinine", 
                              pqqC = "Pyrroloquinoline quinine", 
                              pqqD = "Pyrroloquinoline quinine", 
                              gcd = "Glucose dehydrogenase")) %>%
  mutate(func = factor(func,
                       levels = c("Pyrroloquinoline quinine", "Glucose dehydrogenase")))
func <- t_merge %>%
  group_by(name) %>%
  slice_head(n = 1) %>%
  dplyr::select(name, func)

# Plot - multipanel
gene_plot <- meta_Psol %>%
  pivot_longer(cols = names(meta_Psol)[1:5]) %>%
  mutate(name = factor(name,
                       levels = names(meta_Psol)[1:5])) %>%
  dplyr::select(Site, Entry, name, value) %>%
  left_join(., func, by = "name") %>%
  mutate(func = factor(func,
                       levels = c("Pyrroloquinoline quinine", "Glucose dehydrogenase"))) %>%
  mutate(pc = "P Solubilization Genes")
pdf("InitialFigs/Function_Psol.pdf", width = 6, height = 8)
ggplot(gene_plot, aes(Site, value)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(size = 1, alpha = 0.5, width = 0.2, pch = 16) +
  geom_text(data = t_merge, aes(Site, y, label = str_trim(.group)), 
            size = 3, color = "black") +
  scale_y_continuous(expand = c(0.1, 0.1)) +
  labs(x = NULL,
       y = "PICRUSt2 abundance") +
  facet_nested(func + name ~ pc, scale = "free") +
  theme_bw() +
  theme(legend.position = "none",
        axis.title.y = element_text(size = 12),
        axis.text.x = element_text(size = 10, angle = 45, hjust = 1),
        axis.text.y = element_text(size = 8),
        strip.text.x = element_text(size = 12),
        strip.text.y = element_text(size = 6.5, face = "italic"),
        strip.background = element_rect(fill = "white"))
dev.off()



#### 6. Networks ####
# Make networks for 16S and ITS using SPIEC-EASI
# Identify key taxa with betweenness, degree, participation metrics.
# Then merge the datasets and look for cross domain interactions
# Look at top 100 taxa

# Basic network plot with phyloseq.
top_16S <- as.data.frame(rowMeans(input_filt_rare_16S$data_loaded)) %>%
  set_names("MeanCount") %>%
  arrange(desc(MeanCount)) %>%
  slice_head(n = 100)
input_abund_16S <- filter_taxa_from_input(input_filt_rare_16S,
                                          taxa_IDs_to_keep = rownames(top_16S),
                                          at_spec_level = 7)
nrow(input_abund_16S$data_loaded) # 100
otu_16S <- otu_table(input_abund_16S$data_loaded, taxa_are_rows = T)
tax_16S <- tax_table(as.matrix(input_abund_16S$taxonomy_loaded))
map_16S <- sample_data(input_abund_16S$map_loaded)
input_phy_16S <- phyloseq(otu_16S, tax_16S, map_16S)
se_mb_16S <- spiec.easi(input_phy_16S, 
                        method = 'mb', 
                        lambda.min.ratio = 1e-2,
                        nlambda = 20, 
                        pulsar.params = list(rep.num = 50))
net_16S <- adj2igraph(getRefit(se_mb_16S),  
                      vertex.attr=list(name = taxa_names(input_phy_16S)))
pdf("InitialFigs/Network_16S.pdf", width = 7, height = 5)
plot_network(net_16S, 
             input_phy_16S, 
             type = 'taxa', 
             color = "taxonomy2",
             shape = "taxonomy1",
             point_size = 3,
             label = NULL) +
  labs(shape = "Domain",
       color = "Phylum") +
  guides(shape = guide_legend(order = 1),
         color = guide_legend(override.aes = list(shape = c(16,16,16,16,
                                                            17,16,16,16,
                                                            16,16,17,16)))) +
  ggtitle("16S Network (top 100)") +
  theme(plot.title = element_text(size = 14, face = "bold", hjust = 0.5))
dev.off()

top_ITS <- as.data.frame(rowMeans(input_filt_rare_ITS$data_loaded)) %>%
  set_names("MeanCount") %>%
  arrange(desc(MeanCount)) %>%
  slice_head(n = 100)
input_abund_ITS <- filter_taxa_from_input(input_filt_rare_ITS,
                                          taxa_IDs_to_keep = rownames(top_ITS),
                                          at_spec_level = 7)
nrow(input_abund_ITS$data_loaded) # 100
otu_ITS <- otu_table(input_abund_ITS$data_loaded, taxa_are_rows = T)
tax_ITS <- tax_table(as.matrix(input_abund_ITS$taxonomy_loaded))
map_ITS <- sample_data(input_abund_ITS$map_loaded)
input_phy_ITS <- phyloseq(otu_ITS, tax_ITS, map_ITS)
se_mb_ITS <- spiec.easi(input_phy_ITS, 
                        method = 'mb', 
                        lambda.min.ratio = 1e-2,
                        nlambda = 20, 
                        pulsar.params = list(rep.num = 50))
net_ITS <- adj2igraph(getRefit(se_mb_ITS),  
                      vertex.attr=list(name = taxa_names(input_phy_ITS)))
pdf("InitialFigs/Network_ITS.pdf", width = 7, height = 5)
plot_network(net_ITS, 
             input_phy_ITS, 
             type = 'taxa', 
             color = "taxonomy3",
             point_size = 3,
             label = NULL) +
  labs(color = "Class") +
  ggtitle("ITS Network (top 100)") +
  theme(plot.title = element_text(size = 14, face = "bold", hjust = 0.5))
dev.off()

# Need to create unique prok and euk ESV (ASV) IDs
# Unique ASV IDs are needed for rownames of $data_loaded and $taxonomy_loaded
# Also filter 16S data to ITS data
input_abund_16S <- filter_data(input_abund_16S,
                               filter_cat = "sampleID",
                               keep_vals = input_filt_rare_ITS$map_loaded$sampleID)
input_abund_16S_toMerge <- input_abund_16S
num_16S <- seq(1:nrow(input_abund_16S_toMerge$data_loaded))
ASVlabel_16S <- rep("ASV_Prok_", nrow(input_abund_16S_toMerge$data_loaded))
IDs_16S <- paste(ASVlabel_16S, num_16S, sep = "")
rownames(input_abund_16S_toMerge$data_loaded) <- IDs_16S
rownames(input_abund_16S_toMerge$taxonomy_loaded) <- IDs_16S
input_abund_16S_toMerge$taxonomy_loaded$taxonomy9 <- IDs_16S

input_abund_ITS_toMerge <- input_abund_ITS
num_ITS <- seq(1:nrow(input_abund_ITS_toMerge$data_loaded))
ASVlabel_ITS <- rep("ASV_Euk_", nrow(input_abund_ITS_toMerge$data_loaded))
IDs_ITS <- paste(ASVlabel_ITS, num_ITS, sep = "")
rownames(input_abund_ITS_toMerge$data_loaded) <- IDs_ITS
rownames(input_abund_ITS_toMerge$taxonomy_loaded) <- IDs_ITS
input_abund_ITS_toMerge$taxonomy_loaded$taxonomy9 <- IDs_ITS

input_filt_abund_combined <- input_abund_16S_toMerge
input_filt_abund_combined$data_loaded <- rbind(input_filt_abund_combined$data_loaded,
                                               input_abund_ITS_toMerge$data_loaded)
input_filt_abund_combined$taxonomy_loaded <- rbind(input_filt_abund_combined$taxonomy_loaded,
                                                   input_abund_ITS_toMerge$taxonomy_loaded)
nrow(input_filt_abund_combined$data_loaded) # 200 total taxa

# Convert mctoolsr to phyloseq
names(input_filt_abund_combined$taxonomy_loaded) <- c("Domain", "Phylum", "Class", "Order",
                                                      "Family", "Genus", "Species", "ASV_ID", "ASV_ID2")
otu <- otu_table(input_filt_abund_combined$data_loaded, taxa_are_rows = T)
tax <- tax_table(as.matrix(input_filt_abund_combined$taxonomy_loaded))
map <- sample_data(input_filt_abund_combined$map_loaded)
input.phy <- phyloseq(otu, tax, map)
se.mb2 <- spiec.easi(input.phy, 
                     method='mb', 
                     lambda.min.ratio=1e-2,
                     nlambda=20,
                     pulsar.params=list(rep.num=50))
net <- adj2igraph(getRefit(se.mb2),  
                  vertex.attr=list(name=taxa_names(input.phy)))

# Phyloseq plot
plot_network(net, 
             input.phy, 
             type = 'taxa', 
             color = "Phylum",
             shape = "Domain",
             point_size = 3,
             label = NULL)
plot_network(net, 
             input.phy, 
             type = 'taxa', 
             color = "Phylum",
             shape = "Domain",
             point_size = 3,
             label = NULL,
             layout.method = layout.circle)

# Stats
E(net) # 683
V(net) # 200
transitivity(net) # Average clustering coefficient. 0.177
deg <- degree(net, mode="all")
mean(deg) # 6.83

# Add taxonomic information and color
# Add phylum, which is in the taxonomy_loaded table
# Check phylum numbers
table(input_filt_abund_combined$taxonomy_loaded$Phylum) # 17

# Check match
sum(V(net)$name != rownames(input_filt_abund_combined$taxonomy_loaded)) # Check match

# Make factor
input_filt_abund_combined$taxonomy_loaded$Phylum <- as.factor(input_filt_abund_combined$taxonomy_loaded$Phylum)

# Confer to network
V(net)$phylum = input_filt_abund_combined$taxonomy_loaded$Phylum

# Check levels
levels(input_filt_abund_combined$taxonomy_loaded$Phylum) # There are 17 phyla

# Set n to number of levels
n <- length(levels(input_filt_abund_combined$taxonomy_loaded$Phylum))

# Save taxonomy and colors in tax
tax <- input_filt_abund_combined$taxonomy_loaded

# Get colors for n levels
colrs <- colorRampPalette(brewer.pal(12, "Paired"))(n) # expanded Rcolorbrewer paired palette
tax$color <- tax$Phylum
levels(tax$Phylum)
tax$color <- recode_factor(tax$color,
                           "Acidobacteriota" = colrs[1],
                           "Actinobacteriota" = colrs[2],
                           "Ascomycota" = colrs[3],
                           "Bacteroidota" = colrs[4],
                           "Basidiomycota" = colrs[5],
                           "Chloroflexi" = colrs[6],
                           "Chytridiomycota" = colrs[7],
                           "Crenarchaeota" = colrs[8],
                           "Entotheonellaeota" = colrs[9],
                           "Firmicutes" = colrs[10],
                           "Mortierellomycota" = colrs[11],
                           "Nitrospirota" = colrs[12],
                           "Olpidiomycota" = colrs[13],
                           "Planctomycetota" = colrs[14],
                           "Proteobacteria" = colrs[15],
                           "Thermoplasmatota" = colrs[16],
                           "Verrucomicrobiota" = colrs[17])
V(net)$color <- as.character(tax$color)

# Layout in circle
par(mar = c(0,8,0,0), xpd = TRUE) # bottom, left, top, right
plot(net,
     vertex.color = V(net)$color, 
     vertex.size = deg*2, 
     vertex.shape = "circle", 
     vertex.frame.color = "black",
     vertex.label = NA, 
     #edge.color = ifelse(cor.matrix$r > 0, "#619CFF","#F8766D"),
     edge.curved = 0.2,
     edge.width = 0.2,
     layout = layout_in_circle(net, order = order(V(net)$phylum)))
legend(x = -1.8, y = 0.5, levels(input_filt_abund_combined$taxonomy_loaded$Phylum), 
       pch = 21, col = "black", pt.bg = colrs, pt.cex = 2, cex = 0.8, bty = "n", ncol = 1)

# Refine taxonomy, shapes, color edged, other stylistic things

# Check Modules
rng_adj <- igraph::get.adjacency(net, sparse = FALSE)
netcarto(rng_adj)
edge.betweenness.community(net) # 28 groups
fastgreedy.community(net) # 10 groups
walktrap.community(net) # 22 groups
spinglass.community(net) # Can't work with unconnected graph
leading.eigenvector.community(net) # 8 groups
label.propagation.community(net) # 4 groups
cluster_louvain(net) # 11 groups

# Color edge by weight
se.mb2$lambda
bm <- symBeta(getOptBeta(se.mb2), mode = "maxabs")
diag(bm) <- 0
weights <- Matrix::summary(t(bm))[,3]
E(net)$weight <- weights
E(net)[weight > 0]$color <-"blue" # positive blue
E(net)[weight < 0]$color <-"red"  # negative red
edge_attr(net)

# Vertex shape by role
adj <- as.matrix(as_adjacency_matrix(net))
role <- netcarto(adj)
role <- role[[1]]
v.names <- data.frame(name = V(net)$name)
v.names.ni <- v.names %>%
  filter(name %notin% role$name) %>%
  mutate(module = NA,
         connectivity = NA,
         participation = NA,
         role = NA)
role <- rbind(role, v.names.ni) %>%
  mutate(role = as.factor(role)) %>%
  droplevels() %>%
  mutate(shape = recode_factor(role,
                               #"Peripheral Hub" = "square",
                               "Connector Hub" = "square",
                               #"Kinless Hub" = "square",
                               "Connector" = "diamond",
                               #"Kinless" = "triangle",
                               "Peripheral" = "circle",
                               "Ultra peripheral" = "circle")) %>%
  mutate(shape = as.character(shape)) %>%
  mutate(shape = replace_na(shape, "circle"))
vertex_attr(net)
V(net)$shape <- role$shape
hubcon <- role %>%
  filter(role == "Peripheral Hub" | role == "Connector Hub") %>%
  left_join(., input_filt_abund_combined$taxonomy_loaded, by = c("name" = "ASV_ID2"))

# Info for labels
length(V(net))
length(E(net))
round(mean(deg), 1)
round(transitivity(net), 3)

#### _Color, Shape ####
# Color edges, shape by role

# Plot Circle
par(mar = c(2,8,0,0), xpd = TRUE) # bottom, left, top, right
plot(net,
     vertex.color = V(net)$color,
     vertex.size = deg,
     vertex.shape = V(net)$shape,
     vertex.frame.color = "black",
     vertex.label = NA, 
     edge.curved = 0.2,
     edge.width = 0.5,
     layout = layout_in_circle(net, order = order(V(net)$phylum)))
legend(x = -1.7, y = 1.05, levels(input_filt_abund_combined$taxonomy_loaded$Phylum), 
       pch = 21, col = "black", pt.bg = colrs, pt.cex = 2, cex = 0.8, bty = "n", ncol = 1,
       title = "Taxonomy", title.adj = 0)
legend(x = -1.7, y = -0.3, c("Peripheral", "Connector", "Hub"), pch = c(21, 23, 22), 
       col = "black", pt.cex = 2, cex = 0.8, bty = "n", ncol = 1, 
       title = "Role", title.adj = 0)
legend(x = -1.7, y = -0.65, c("Positive", "Negative"), lwd = 2,
       col = c("blue", "red"), cex = 0.8, bty = "n", ncol = 1, 
       title = "Association", title.adj = 0)
title(main = "Sunflower Rhizosphere Network (SPIEC-EASI)", adj = 0.5, line = -1.5, cex.main = 1.5)
mtext("Nodes = 200", side = 1, line = -2, cex = 0.8)
mtext("Edges = 683", side = 1, line = -1, cex = 0.8)
mtext("Mean Degrees = 6.8", side = 1, line = 0, cex = 0.8)
mtext("Clustering Coefficient = 0.177", side = 1, line = 1, cex = 0.8)

# Update taxonomy
# Need to label Archaea, Bacteria, Fungi (A_, B_, F_)
View(input_filt_abund_combined$taxonomy_loaded)
input_filt_abund_combined$taxonomy_loaded <- input_filt_abund_combined$taxonomy_loaded %>%
  mutate(Phylum = as.character(Phylum)) %>%
  mutate(taxonomy = paste(Domain, Phylum, sep = "_")) %>%
  mutate(taxonomy = gsub("Archaea", "A", taxonomy)) %>%
  mutate(taxonomy = gsub("Bacteria", "B", taxonomy)) %>%
  mutate(taxonomy = gsub("Fungi", "F", taxonomy))

# Run as above but with taxonomy instead of Phylum
sum(V(net)$name != rownames(input_filt_abund_combined$taxonomy_loaded)) # Check match
input_filt_abund_combined$taxonomy_loaded$taxonomy <- as.factor(input_filt_abund_combined$taxonomy_loaded$taxonomy)
V(net)$taxonomy = input_filt_abund_combined$taxonomy_loaded$taxonomy
levels(input_filt_abund_combined$taxonomy_loaded$taxonomy) # There are 17 phyla
n <- length(levels(input_filt_abund_combined$taxonomy_loaded$taxonomy))
tax <- input_filt_abund_combined$taxonomy_loaded
colrs <- colorRampPalette(brewer.pal(12, "Paired"))(n) # expanded Rcolorbrewer paired palette
tax$color <- tax$taxonomy
levels(tax$taxonomy)
tax$color <- recode_factor(tax$color,
                           "A_Crenarchaeota" = colrs[1],
                           "A_Thermoplasmatota" = colrs[2],
                           "B_Acidobacteriota" = colrs[3],
                           "B_Actinobacteriota" = colrs[4],
                           "B_Bacteroidota" = colrs[5],
                           "B_Chloroflexi" = colrs[6],
                           "B_Entotheonellaeota" = colrs[7],
                           "B_Firmicutes" = colrs[8],
                           "B_Nitrospirota" = colrs[9],
                           "B_Planctomycetota" = colrs[10],
                           "B_Proteobacteria" = colrs[11],
                           "B_Verrucomicrobiota" = colrs[12],
                           "F_Ascomycota" = colrs[13],
                           "F_Basidiomycota" = colrs[14],
                           "F_Chytridiomycota" = colrs[15],
                           "F_Mortierellomycota" = colrs[16],
                           "F_Olpidiomycota" = colrs[17])
V(net)$color <- as.character(tax$color)

# And add ASV ID as vertex attribute
V(net)$ASV_ID2 <- input_filt_abund_combined$taxonomy_loaded$ASV_ID2

#### _Good Taxonomy ####
# Replot, Circle, shape by Role
par(mar = c(3,8,1,0), xpd = TRUE) # bottom, left, top, right
plot(net,
     vertex.color = V(net)$color,
     vertex.size = deg,
     vertex.shape = V(net)$shape,
     vertex.frame.color = "black",
     vertex.label = NA, 
     edge.curved = 0.2,
     edge.width = 0.5,
     layout = layout_in_circle(net, order = order(V(net)$taxonomy)))
legend(x = -1.8, y = 1.1, levels(input_filt_abund_combined$taxonomy_loaded$taxonomy), 
       pch = 21, col = "black", pt.bg = colrs, pt.cex = 2, cex = 0.8, bty = "n", ncol = 1,
       title = "Taxonomy", title.adj = 0)
legend(x = -1.8, y = -0.3, c("Peripheral", "Connector", "Hub"), pch = c(21, 23, 22), 
       col = "black", pt.cex = 2, cex = 0.8, bty = "n", ncol = 1, 
       title = "Role", title.adj = 0)
legend(x = -1.8, y = -0.7, c("Positive", "Negative"), lwd = 2,
       col = c("blue", "red"), cex = 0.8, bty = "n", ncol = 1, 
       title = "Association", title.adj = 0)
title(main = "Sunflower Rhizosphere Network (SPIEC-EASI)", adj = 0.5, line = -1, cex.main = 1.5)
mtext("Nodes = 200", side = 1, line = -1.5, cex = 0.8)
mtext("Edges = 683", side = 1, line = -0.5, cex = 0.8)
mtext("Mean Degrees = 6.8", side = 1, line = 0.5, cex = 0.8)
mtext("Clustering Coefficient = 0.177", side = 1, line = 1.5, cex = 0.8)

# Replot, Circle, shape by Domain. Save this one.
sum(V(net)$name != rownames(input_filt_abund_combined$taxonomy_loaded)) # Check match
V(net)$Domain = input_filt_abund_combined$taxonomy_loaded$Domain
V(net)$Domain <- gsub("Archaea", "square", V(net)$Domain)
V(net)$Domain <- gsub("Bacteria", "circle", V(net)$Domain)
V(net)$Domain <- gsub("Fungi", "triangle", V(net)$Domain)

pdf("InitialFigs/Network_Combined.pdf", width = 7, height = 6)
par(mar = c(3,8,1,0), xpd = TRUE) # bottom, left, top, right
plot(net,
     vertex.color = V(net)$color,
     vertex.size = deg,
     vertex.shape = V(net)$Domain,
     vertex.frame.color = "black",
     vertex.label = NA, 
     edge.curved = 0.2,
     edge.width = 0.5,
     layout = layout_in_circle(net, order = order(V(net)$taxonomy)))
legend(x = -1.8, y = 1.1, levels(input_filt_abund_combined$taxonomy_loaded$taxonomy), 
       pch = c(rep(22, 2), rep(21, 10), rep(24, 15)), col = "black", pt.bg = colrs, pt.cex = 1.5, 
       cex = 0.8, bty = "n", ncol = 1,
       title = "Taxonomy", title.cex = 1.1, title.adj = 0, y.intersp = 1.1)
legend(x = -1.8, y = -0.4, c("Archaea", "Bacteria", "Fungi"), pch = c(22, 21, 24), 
       col = "black", pt.cex = 1.5, cex = 0.8, bty = "n", ncol = 1, 
       title = "Domain", title.cex = 1.1, title.adj = 0, y.intersp = 1.1)
legend(x = -1.8, y = -0.8, c("Positive", "Negative"), lwd = 2,
       col = c("blue", "red"), cex = 0.8, bty = "n", ncol = 1, 
       title = "Association", title.cex = 1.1, title.adj = 0)
title(main = "Sunflower Rhizosphere Network (SPIEC-EASI)", adj = 0.5, line = -1, cex.main = 1.2)
mtext("Nodes = 200", side = 1, line = -1.5, cex = 0.8)
mtext("Edges = 683", side = 1, line = -0.5, cex = 0.8)
mtext("Mean Degrees = 6.8", side = 1, line = 0.5, cex = 0.8)
mtext("Clustering Coefficient = 0.177", side = 1, line = 1.5, cex = 0.8)
dev.off()




#### _Betweenness ####
# Plot degree versus betweenness for each network
se <- se.mb2
net <- adj2igraph(getRefit(se),  vertex.attr=list(name=taxa_names(input.phy)))
bw <- data.frame("Degree" = degree(net, mode="all"),
                 "Betweenness" = betweenness(net)) %>%
  mutate("ASV" = rownames(.)) %>%
  left_join(., input_filt_abund_combined$taxonomy_loaded, by = c("ASV" = "ASV_ID2")) %>%
  mutate(Phylum = as.factor(Phylum))

nb.cols <- length(levels(bw$Phylum))
mycolors <- colorRampPalette(brewer.pal(12, "Paired"))(nb.cols)
ggplot(bw, aes(Degree, Betweenness, colour = taxonomy, shape = Domain)) +
  geom_jitter(size = 2, alpha = 1, width = 0.15) +
  labs(x = "Degree",
       y = "Betweenness",
       colour = "Phylum") +
  scale_color_manual(values = mycolors) +
  scale_x_continuous(limits = c(-0.5, 18.5),
                     breaks = c(0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18)) +
  theme_bw() +
  theme(legend.key.size = unit(0.5, "cm"),
        legend.spacing.y = unit(0, "cm"))

# Check ASVs, add labels
ggplot(bw, aes(Degree, Betweenness, colour = taxonomy, shape = Domain)) +
  geom_jitter(size = 2, alpha = 1, width = 0.15) +
  geom_text(data = bw,
            aes(x = Degree, y = Betweenness, label = ASV), 
            size = 3, inherit.aes = F) +
  labs(x = "Degree",
       y = "Betweenness",
       colour = "Phylum") +
  scale_color_manual(values = mycolors) +
  scale_x_continuous(limits = c(-0.5, 18.5),
                     breaks = c(0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18)) +
  theme_bw() +
  theme(legend.key.size = unit(0.5, "cm"),
        legend.spacing.y = unit(0, "cm"))

bw_asvs <- bw %>%
  filter(Betweenness > 600) %>%
  mutate("HighTax" = ifelse(Genus != "NA",
                            Genus,
                            ifelse(Family != "NA",
                                   Family,
                                   ifelse(Order != "NA",
                                          Order,
                                          ifelse(Class != "NA",
                                                 Class,
                                                 ifelse(Phylum != "NA",
                                                        Phylum,
                                                        Domain)))))) %>%
  mutate("HightTax_ASV" = paste(HighTax, ASV_ID, sep = "_"))

pdf("InitialFigs/Network_Betweenness.pdf", width = 7, height = 5)
ggplot(bw, aes(Degree, Betweenness, colour = taxonomy, shape = Domain)) +
  geom_jitter(size = 2, alpha = 1, width = 0.15) +
  geom_text_repel(data = bw_asvs,
                  #min.segment.length = 0,
                  aes(x = Degree, y = Betweenness, label = HightTax_ASV), 
                  size = 2, inherit.aes = F,
                  position = position_jitter(width = 0.15)) +
  labs(x = "Degree",
       y = "Betweenness",
       colour = "Phylum") +
  scale_color_manual(values = mycolors) +
  scale_x_continuous(limits = c(-0.5, 18.5),
                     breaks = c(0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18)) +
  theme_bw() +
  theme(legend.key.size = unit(0.5, "cm"),
        legend.spacing.y = unit(0, "cm"))
dev.off()



#### _Participation ####
# Plot participation coefficient and within-module degree (Barnes et al.)
# Dashed lines at 0.61 and 2.2
# z-score (within module degree) is "connectivity"
# participation coefficient P is "participation
roles <- role %>%
  left_join(., input_filt_abund_combined$taxonomy_loaded, by = c("name" = "ASV_ID2")) %>%
  mutate(Phylum = as.factor(Phylum),
         role = as.factor(role)) %>%
  filter(is.na(module) == F) %>%
  droplevels()
nb.cols <- length(levels(roles$Phylum))
mycolors <- colorRampPalette(brewer.pal(12, "Paired"))(nb.cols)
ggplot(roles, aes(participation, connectivity, colour = taxonomy, shape = role)) +
  geom_vline(xintercept = 0.62, linetype = "dashed", colour = "gray", linewidth = 0.25) +
  geom_hline(yintercept = 2.5, linetype = "dashed", colour = "gray", linewidth = 0.25) +
  geom_point(size = 2, alpha = 0.75) +
  labs(x = "Participation coefficient (P)",
       y = "Within module degree (z)",
       colour = "Phylum",
       shape = "Network Role") +
  scale_color_manual(values = mycolors) +
  theme_bw() +
  theme(legend.key.size = unit(0.5, "cm"),
        legend.spacing.y = unit(0, "cm"),
        panel.grid = element_blank())

# Check ASVs, add labels
ggplot(roles, aes(participation, connectivity, colour = taxonomy, shape = role)) +
  geom_vline(xintercept = 0.62, linetype = "dashed", colour = "gray", linewidth = 0.25) +
  geom_hline(yintercept = 2.5, linetype = "dashed", colour = "gray", linewidth = 0.25) +
  geom_point(size = 2, alpha = 0.75) +
  geom_text(data = roles,
            aes(x = participation, y = connectivity, label = name), 
            size = 3, inherit.aes = F) +
  labs(x = "Participation coefficient (P)",
       y = "Within module degree (z)",
       colour = "Phylum",
       shape = "Network Role") +
  scale_color_manual(values = mycolors) +
  theme_bw() +
  theme(legend.key.size = unit(0.5, "cm"),
        legend.spacing.y = unit(0, "cm"),
        panel.grid = element_blank())

pz_asvs <- roles %>%
  filter(connectivity > 2.5 | participation > 0.62) %>%
  mutate("HighTax" = ifelse(Genus != "NA",
                            Genus,
                            ifelse(Family != "NA",
                                   Family,
                                   ifelse(Order != "NA",
                                          Order,
                                          Class))))

pdf("InitialFigs/Network_Participation.pdf", width = 7, height = 5)
ggplot(roles, aes(participation, connectivity, colour = taxonomy, shape = role)) +
  geom_vline(xintercept = 0.62, linetype = "dashed", colour = "gray", linewidth = 0.25) +
  geom_hline(yintercept = 2.5, linetype = "dashed", colour = "gray", linewidth = 0.25) +
  geom_point(size = 2, alpha = 0.75) +
  geom_text_repel(data = pz_asvs,
                  min.segment.length = 0,
                  max.overlaps = 20,
                  aes(x = participation, y = connectivity, label = HighTax), 
                  size = 2, inherit.aes = F) +
  labs(x = "Participation coefficient (P)",
       y = "Within module degree (z)",
       colour = "Phylum",
       shape = "Network Role") +
  scale_color_manual(values = mycolors) +
  theme_bw() +
  theme(legend.key.size = unit(0.4, "cm"),
        legend.margin = margin(0,0,0,0),
        panel.grid = element_blank())
dev.off()

keystone <- roles %>%
  filter(participation > 0.62 | connectivity > 2.5)

#write_xlsx(keystone, format_headers = F, "keystone_taxa.xlsx")



#### 7. Venn ####
# Any taxa in all samples?
# 1 bacterium, 0 fungi

# Taxa in all sites?
# 343 16S, 28 fungi

# Taxa in all genotypes?
# 9330 16S, 606 fungi

# How many ASVs are in all sites or only in 1 site? 
# How many ASVs are in all genotypes or only in 1 genotype?
# Need to do Venn Diagram with 15 or 10 levels!
# Not possible with mctoolsr - Try ps_venn from MicEco. 
# Input is phyloseq object
# Output is list of taxa in single sites and all combinations
# -Can then filter out single sites and all sites
# Can make similar plot to indicator species with number in all or unique, colored by taxonomy

# Do with ASVs and OTUs
# Do with original data and then with OTUs present in only 1 sample filtered out
# Already calculated prevalence in the Taxa Barplots section
prev_16S
prev_ITS

# ASVs present in only 1 sample
prev_16S_n1 <- prev_16S %>%
  filter(Present == 1) # n = 902
prev_ITS_n1 <- prev_ITS %>%
  filter(Present == 1) # n = 1885

# OTUs present in only 1 sample
prev_16S_OTU <- data.frame("OTU_ID" = rownames(tax_sum_OTU_16S),
                           "Absent" = rowSums(tax_sum_OTU_16S==0)) %>%
  mutate(Present = ncol(tax_sum_OTU_16S) - Absent) %>%
  mutate(Present_Perc = Present/ncol(tax_sum_OTU_16S)*100)
prev_ITS_OTU <- data.frame("OTU_ID" = rownames(tax_sum_OTU_ITS),
                           "Absent" = rowSums(tax_sum_OTU_ITS==0)) %>%
  mutate(Present = ncol(tax_sum_OTU_ITS) - Absent) %>%
  mutate(Present_Perc = Present/ncol(tax_sum_OTU_ITS)*100)

prev_16S_OTU_n1 <- prev_16S_OTU %>%
  filter(Present == 1) # n = 180
prev_ITS_OTU_n1 <- prev_ITS_OTU %>%
  filter(Present == 1) # n = 661

# Filter
input_filt_rare_16S_n2 <- filter_taxa_from_input(input_filt_rare_16S,
                                                 taxa_IDs_to_remove = prev_16S_n1$ASV_ID)
input_filt_rare_ITS_n2 <- filter_taxa_from_input(input_filt_rare_ITS,
                                                 taxa_IDs_to_remove = prev_ITS_n1$ASV_ID)
nrow(input_filt_rare_16S_n2$data_loaded) # 26991
nrow(input_filt_rare_ITS_n2$data_loaded) # 3409



#### _Site ####
# 16S ASV
# 16S ASV 2 sample
# 16S OTU
# 16S OTU 2 sample
# ITS ASV
# ITS ASV 2 sample
# ITS OTU
# ITS OTU 2 sample

# For paper use ASV 2 sample (Figure S9)

#### __16S ####
# ASV level, n = 1 included
otu_16S <- otu_table(input_filt_rare_16S$data_loaded, taxa_are_rows = T)
tax_16S <- tax_table(as.matrix(input_filt_rare_16S$taxonomy_loaded))
map_16S <- sample_data(input_filt_rare_16S$map_loaded)
input_phy_16S <- phyloseq(otu_16S, tax_16S, map_16S)
venn_site_16S <- ps_venn(ps = input_phy_16S,
                         group = "Site",
                         fraction = 0.000000000000000000001,
                         weight = FALSE,
                         relative = FALSE,
                         plot = FALSE)
#saveRDS(venn_site_16S, "data/venn_site_16S.rds") # Slow, so save.
length(venn_site_16S)
length(venn_site_16S[[16]]) # 343 ASV present once at all sites
names(venn_site_16S)[1:16]
venn_site_16S_sub <- venn_site_16S[1:16]
site_names <- c(levels(input_filt_rare_16S$map_loaded$Site), "All")
for (i in 1:length(venn_site_16S_sub)) {
  venn_site_16S_sub[[i]] <- as.data.frame(venn_site_16S_sub[[i]]) %>%
    set_names(c("ASV_ID")) %>%
    mutate(Site = site_names[i])
}
venn_site_16S_toPlot <- venn_site_16S_sub[[1]]
for (i in 2:length(venn_site_16S_sub)) {
  venn_site_16S_toPlot <- rbind(venn_site_16S_toPlot, venn_site_16S_sub[[i]])
}
venn_site_16S_toPlot <- venn_site_16S_toPlot %>%
  left_join(., input_filt_rare_16S$taxonomy_loaded, by = c("ASV_ID" = "taxonomy7"))
top_taxa_16S <- as.data.frame(table(venn_site_16S_toPlot$taxonomy2)) %>%
  arrange(desc(Freq)) %>%
  slice_head(n = 12) %>%
  mutate(Var1 = as.character(Var1))
venn_site_16S_toPlot <- venn_site_16S_toPlot %>%
  group_by(Site, taxonomy2) %>%
  summarise(n_taxa = n()) %>%
  drop_na(taxonomy2)
venn_site_16S_toPlot$taxonomy2[venn_site_16S_toPlot$taxonomy2 %notin% top_taxa_16S$Var1] = 'Other'
venn_site_16S_toPlot <- venn_site_16S_toPlot %>%
  group_by(Site, taxonomy2) %>%
  summarise(n_taxa = sum(n_taxa)) %>%
  mutate(taxonomy2 = factor(taxonomy2,
                            levels = c("Other", rev(top_taxa_16S$Var1))))
pdf("InitialFigs/Venn_Site_16S.pdf", width = 7, height = 5)
ggplot(venn_site_16S_toPlot, aes(Site, n_taxa, fill = taxonomy2, group = taxonomy2)) +
  geom_bar(stat = "identity", colour = NA, linewidth = 0.25) +
  labs(x = "Site",
       y = "Number of ASVs",
       fill = "Phylum") +
  scale_fill_manual(values = c("grey75", brewer.pal(12, "Paired")[12:1])) +
  scale_y_continuous(expand = c(0.01, 0.01)) +  
  theme_classic() +
  theme(axis.title.y = element_text(size = 12),
        axis.title.x = element_blank(),
        axis.text.y = element_text(size = 10),
        axis.text.x = element_text(size = 10, angle = 45, hjust = 1))
dev.off()



#ASV level, n = 2 cutoff
otu_16S <- otu_table(input_filt_rare_16S_n2$data_loaded, taxa_are_rows = T)
tax_16S <- tax_table(as.matrix(input_filt_rare_16S_n2$taxonomy_loaded))
map_16S <- sample_data(input_filt_rare_16S_n2$map_loaded)
input_phy_16S <- phyloseq(otu_16S, tax_16S, map_16S)
venn_site_16S <- ps_venn(ps = input_phy_16S,
                         group = "Site",
                         fraction = 0.000000000000000000001,
                         weight = FALSE,
                         relative = FALSE,
                         plot = FALSE)
#saveRDS(venn_site_16S, "data/venn_site_16S_n2.rds") # Slow, so save.
venn_site_16S <- readRDS("data/venn_site_16S_n2.rds")
length(venn_site_16S)
length(venn_site_16S[[16]]) # 343 ASV present once at all sites
names(venn_site_16S)[1:16]
names(venn_site_16S)[6] <- "Brighton"
venn_site_16S_sub <- venn_site_16S[1:16]
site_names <- c(names(venn_site_16S)[1:15], "All")
for (i in 1:length(venn_site_16S_sub)) {
  venn_site_16S_sub[[i]] <- as.data.frame(venn_site_16S_sub[[i]]) %>%
    set_names(c("ASV_ID")) %>%
    mutate(Site = site_names[i])
}
venn_site_16S_toPlot <- venn_site_16S_sub[[1]]
for (i in 2:length(venn_site_16S_sub)) {
  venn_site_16S_toPlot <- rbind(venn_site_16S_toPlot, venn_site_16S_sub[[i]])
}
venn_site_16S_toPlot <- venn_site_16S_toPlot %>%
  left_join(., input_filt_rare_16S$taxonomy_loaded, by = c("ASV_ID" = "taxonomy7"))
top_taxa_16S <- as.data.frame(table(venn_site_16S_toPlot$taxonomy2)) %>%
  arrange(desc(Freq)) %>%
  slice_head(n = 12) %>%
  mutate(Var1 = as.character(Var1))
venn_site_16S_toPlot <- venn_site_16S_toPlot %>%
  group_by(Site, taxonomy2) %>%
  summarise(n_taxa = n()) %>%
  drop_na(taxonomy2)
venn_site_16S_toPlot$taxonomy2[venn_site_16S_toPlot$taxonomy2 %notin% top_taxa_16S$Var1] = 'Other'
venn_site_16S_toPlot <- venn_site_16S_toPlot %>%
  group_by(Site, taxonomy2) %>%
  summarise(n_taxa = sum(n_taxa)) %>%
  mutate(taxonomy2 = factor(taxonomy2,
                            levels = c("Other", rev(top_taxa_16S$Var1))))
n_per_site_16S_venn <- venn_site_16S_toPlot %>%
  group_by(Site) %>%
  summarise(site_n = sum(n_taxa))
pdf("InitialFigs/Venn_Site_16S_n2.pdf", width = 7, height = 5)
ggplot(venn_site_16S_toPlot, aes(Site, n_taxa, fill = taxonomy2, group = taxonomy2)) +
  geom_bar(stat = "identity", colour = NA, linewidth = 0.25) +
  labs(x = "Site",
       y = "Number of ASVs",
       fill = "Phylum") +
  scale_fill_manual(values = c("grey75", brewer.pal(12, "Paired")[12:1])) +
  scale_y_continuous(expand = c(0.01, 0.01)) +  
  theme_classic() +
  theme(axis.title.y = element_text(size = 12),
        axis.title.x = element_blank(),
        axis.text.y = element_text(size = 10),
        axis.text.x = element_text(size = 10, angle = 45, hjust = 1))
dev.off()

figs8a <- ggplot(venn_site_16S_toPlot, aes(Site, n_taxa, fill = taxonomy2, group = taxonomy2)) +
  geom_bar(stat = "identity", colour = NA, linewidth = 0.25) +
  labs(x = "Site",
       y = "# of ASVs",
       fill = "Phylum") +
  scale_fill_manual(values = c("grey75", brewer.pal(12, "Paired")[12:1])) +
  scale_y_continuous(expand = c(0.01, 0.01)) +  
  theme_classic() +
  theme(axis.title.y = element_text(size = 12),
        axis.title.x = element_blank(),
        axis.text.y = element_text(size = 10),
        axis.text.x = element_blank(),
        legend.margin = margin(0, 0, 0, 0, unit = "pt"),
        legend.box.margin = margin(0, 0, 0, 0, unit = "pt"),
        legend.key.size = unit(0.3, "cm"),
        legend.text = element_text(size = 6))
figs8a
legs8a <- get_legend(figs8a)
figs8a <- figs8a + theme(legend.position = "none")



# OTU level, n = 1 included
otu_tax_16S <- input_filt_rare_16S$taxonomy_loaded %>%
  group_by(OTU_ID) %>%
  slice_head(n = 1) %>%
  as.data.frame()
rownames(otu_tax_16S) <- otu_tax_16S$OTU_ID
otu_16S <- otu_table(tax_sum_OTU_16S, taxa_are_rows = T)
tax_16S <- tax_table(as.matrix(otu_tax_16S))
map_16S <- sample_data(input_filt_rare_16S$map_loaded)
input_phy_16S <- phyloseq(otu_16S, tax_16S, map_16S)
venn_site_16S <- ps_venn(ps = input_phy_16S,
                         group = "Site",
                         fraction = 0.000000000000000000001,
                         weight = FALSE,
                         relative = FALSE,
                         plot = FALSE)
#saveRDS(venn_site_16S, "data/venn_site_OTU_16S.rds") # Slow, so save.
length(venn_site_16S)
length(venn_site_16S[[16]]) # 580 OTU present once at all sites
names(venn_site_16S)[1:16]
venn_site_16S_sub <- venn_site_16S[1:16]
site_names <- c(levels(input_filt_rare_16S$map_loaded$Site), "All")
for (i in 1:length(venn_site_16S_sub)) {
  venn_site_16S_sub[[i]] <- as.data.frame(venn_site_16S_sub[[i]]) %>%
    set_names(c("ASV_ID")) %>%
    mutate(Site = site_names[i])
}
venn_site_16S_toPlot <- venn_site_16S_sub[[1]]
for (i in 2:length(venn_site_16S_sub)) {
  venn_site_16S_toPlot <- rbind(venn_site_16S_toPlot, venn_site_16S_sub[[i]])
}
venn_site_16S_toPlot <- venn_site_16S_toPlot %>%
  left_join(., otu_tax_16S, by = c("ASV_ID" = "OTU_ID"))
top_taxa_16S <- as.data.frame(table(venn_site_16S_toPlot$taxonomy2)) %>%
  arrange(desc(Freq)) %>%
  slice_head(n = 12) %>%
  mutate(Var1 = as.character(Var1))
venn_site_16S_toPlot <- venn_site_16S_toPlot %>%
  group_by(Site, taxonomy2) %>%
  summarise(n_taxa = n()) %>%
  drop_na(taxonomy2)
venn_site_16S_toPlot$taxonomy2[venn_site_16S_toPlot$taxonomy2 %notin% top_taxa_16S$Var1] = 'Other'
venn_site_16S_toPlot <- venn_site_16S_toPlot %>%
  group_by(Site, taxonomy2) %>%
  summarise(n_taxa = sum(n_taxa)) %>%
  mutate(taxonomy2 = factor(taxonomy2,
                            levels = c("Other", rev(top_taxa_16S$Var1))))
pdf("InitialFigs/Venn_Site_OTU_16S.pdf", width = 7, height = 5)
ggplot(venn_site_16S_toPlot, aes(Site, n_taxa, fill = taxonomy2, group = taxonomy2)) +
  geom_bar(stat = "identity", colour = NA, linewidth = 0.25) +
  labs(x = "Site",
       y = "Number of OTUs",
       fill = "Phylum") +
  scale_fill_manual(values = c("grey75", brewer.pal(12, "Paired")[12:1])) +
  scale_y_continuous(expand = c(0.01, 0.01)) +  
  theme_classic() +
  theme(axis.title.y = element_text(size = 12),
        axis.title.x = element_blank(),
        axis.text.y = element_text(size = 10),
        axis.text.x = element_text(size = 10, angle = 45, hjust = 1))
dev.off()



# OTU level, n = 2 cutoff
otu_tax_16S_n2 <- input_filt_rare_16S_n2$taxonomy_loaded %>%
  group_by(OTU_ID) %>%
  slice_head(n = 1) %>%
  as.data.frame()
rownames(otu_tax_16S_n2) <- otu_tax_16S_n2$OTU_ID
tax_sum_OTU_16S_n2 <- summarize_taxonomy(input_filt_rare_16S_n2, level = 9, report_higher_tax = F)
otu_16S <- otu_table(tax_sum_OTU_16S_n2, taxa_are_rows = T)
tax_16S <- tax_table(as.matrix(otu_tax_16S_n2))
map_16S <- sample_data(input_filt_rare_16S_n2$map_loaded)
input_phy_16S <- phyloseq(otu_16S, tax_16S, map_16S)
venn_site_16S <- ps_venn(ps = input_phy_16S,
                         group = "Site",
                         fraction = 0.000000000000000000001,
                         weight = FALSE,
                         relative = FALSE,
                         plot = FALSE)
#saveRDS(venn_site_16S, "data/venn_site_OTU_16S_n2.rds") # Slow, so save.
length(venn_site_16S)
length(venn_site_16S[[16]]) # 580 OTU present once at all sites
names(venn_site_16S)[1:16]
venn_site_16S_sub <- venn_site_16S[1:16]
site_names <- c(levels(input_filt_rare_16S_n2$map_loaded$Site), "All")
for (i in 1:length(venn_site_16S_sub)) {
  venn_site_16S_sub[[i]] <- as.data.frame(venn_site_16S_sub[[i]]) %>%
    set_names(c("ASV_ID")) %>%
    mutate(Site = site_names[i])
}
venn_site_16S_toPlot <- venn_site_16S_sub[[1]]
for (i in 2:length(venn_site_16S_sub)) {
  venn_site_16S_toPlot <- rbind(venn_site_16S_toPlot, venn_site_16S_sub[[i]])
}
venn_site_16S_toPlot <- venn_site_16S_toPlot %>%
  left_join(., otu_tax_16S_n2, by = c("ASV_ID" = "OTU_ID"))
top_taxa_16S <- as.data.frame(table(venn_site_16S_toPlot$taxonomy2)) %>%
  arrange(desc(Freq)) %>%
  slice_head(n = 12) %>%
  mutate(Var1 = as.character(Var1))
venn_site_16S_toPlot <- venn_site_16S_toPlot %>%
  group_by(Site, taxonomy2) %>%
  summarise(n_taxa = n()) %>%
  drop_na(taxonomy2)
venn_site_16S_toPlot$taxonomy2[venn_site_16S_toPlot$taxonomy2 %notin% top_taxa_16S$Var1] = 'Other'
venn_site_16S_toPlot <- venn_site_16S_toPlot %>%
  group_by(Site, taxonomy2) %>%
  summarise(n_taxa = sum(n_taxa)) %>%
  mutate(taxonomy2 = factor(taxonomy2,
                            levels = c("Other", rev(top_taxa_16S$Var1))))
pdf("InitialFigs/Venn_Site_OTU_16S_n2.pdf", width = 7, height = 5)
ggplot(venn_site_16S_toPlot, aes(Site, n_taxa, fill = taxonomy2, group = taxonomy2)) +
  geom_bar(stat = "identity", colour = NA, linewidth = 0.25) +
  labs(x = "Site",
       y = "Number of OTUs",
       fill = "Phylum") +
  scale_fill_manual(values = c("grey75", brewer.pal(12, "Paired")[12:1])) +
  scale_y_continuous(expand = c(0.01, 0.01)) +  
  theme_classic() +
  theme(axis.title.y = element_text(size = 12),
        axis.title.x = element_blank(),
        axis.text.y = element_text(size = 10),
        axis.text.x = element_text(size = 10, angle = 45, hjust = 1))
dev.off()



#### __ITS ####
# ASV Level, n = 1 included
otu_ITS <- otu_table(input_filt_rare_ITS$data_loaded, taxa_are_rows = T)
tax_ITS <- tax_table(as.matrix(input_filt_rare_ITS$taxonomy_loaded))
map_ITS <- sample_data(input_filt_rare_ITS$map_loaded)
input_phy_ITS <- phyloseq(otu_ITS, tax_ITS, map_ITS)
venn_site_ITS <- ps_venn(ps = input_phy_ITS,
                         group = "Site",
                         fraction = 0.000000000000000000001,
                         weight = FALSE,
                         relative = FALSE,
                         plot = FALSE)
length(venn_site_ITS)
length(venn_site_ITS[[16]]) # 28 ASV present once at all sites
names(venn_site_ITS)[1:16]
venn_site_ITS_sub <- venn_site_ITS[1:16]
site_names <- c(levels(input_filt_rare_ITS_n2$map_loaded$Site), "All")
for (i in 1:length(venn_site_ITS_sub)) {
  venn_site_ITS_sub[[i]] <- as.data.frame(venn_site_ITS_sub[[i]]) %>%
    set_names(c("ASV_ID")) %>%
    mutate(Site = site_names[i])
}
venn_site_ITS_toPlot <- venn_site_ITS_sub[[1]]
for (i in 2:length(venn_site_ITS_sub)) {
  venn_site_ITS_toPlot <- rbind(venn_site_ITS_toPlot, venn_site_ITS_sub[[i]])
}
venn_site_ITS_toPlot <- venn_site_ITS_toPlot %>%
  left_join(., input_filt_rare_ITS_n2$taxonomy_loaded, by = c("ASV_ID" = "taxonomy8"))
top_taxa_ITS <- as.data.frame(table(venn_site_ITS_toPlot$taxonomy3)) %>%
  arrange(desc(Freq)) %>%
  slice_head(n = 12) %>%
  mutate(Var1 = as.character(Var1))
venn_site_ITS_toPlot <- venn_site_ITS_toPlot %>%
  group_by(Site, taxonomy3) %>%
  summarise(n_taxa = n()) %>%
  drop_na(taxonomy3)
venn_site_ITS_toPlot$taxonomy3[venn_site_ITS_toPlot$taxonomy3 %notin% top_taxa_ITS$Var1] = 'Other'
venn_site_ITS_toPlot <- venn_site_ITS_toPlot %>%
  group_by(Site, taxonomy3) %>%
  summarise(n_taxa = sum(n_taxa)) %>%
  mutate(taxonomy3 = factor(taxonomy3,
                            levels = c("Other", "NA",
                                       rev(subset(top_taxa_ITS$Var1, top_taxa_ITS$Var1 != "NA")))))
pdf("InitialFigs/Venn_Site_ITS.pdf", width = 7, height = 5)
ggplot(venn_site_ITS_toPlot, aes(Site, n_taxa, fill = taxonomy3, group = taxonomy3)) +
  geom_bar(stat = "identity", colour = NA, linewidth = 0.25) +
  labs(x = "Site",
       y = "Number of ASVs",
       fill = "Class") +
  scale_fill_manual(values = c("grey90", "grey75", brewer.pal(11, "Paired")[11:1])) +
  scale_y_continuous(expand = c(0.01, 0.01)) +  
  theme_classic() +
  theme(axis.title.y = element_text(size = 12),
        axis.title.x = element_blank(),
        axis.text.y = element_text(size = 10),
        axis.text.x = element_text(size = 10, angle = 45, hjust = 1))
dev.off()



# ASV level, n = 2 cutoff
otu_ITS <- otu_table(input_filt_rare_ITS_n2$data_loaded, taxa_are_rows = T)
tax_ITS <- tax_table(as.matrix(input_filt_rare_ITS_n2$taxonomy_loaded))
map_ITS <- sample_data(input_filt_rare_ITS_n2$map_loaded)
input_phy_ITS <- phyloseq(otu_ITS, tax_ITS, map_ITS)
venn_site_ITS <- ps_venn(ps = input_phy_ITS,
                         group = "Site",
                         fraction = 0.000000000000000000001,
                         weight = FALSE,
                         relative = FALSE,
                         plot = FALSE)
length(venn_site_ITS)
length(venn_site_ITS[[16]]) # 28 ASV present once at all sites
names(venn_site_ITS)[1:16]
venn_site_ITS_sub <- venn_site_ITS[1:16]
site_names <- c(levels(as.factor(input_filt_rare_ITS_n2$map_loaded$Site)), "All")
for (i in 1:length(venn_site_ITS_sub)) {
  venn_site_ITS_sub[[i]] <- as.data.frame(venn_site_ITS_sub[[i]]) %>%
    set_names(c("ASV_ID")) %>%
    mutate(Site = site_names[i])
}
venn_site_ITS_toPlot <- venn_site_ITS_sub[[1]]
for (i in 2:length(venn_site_ITS_sub)) {
  venn_site_ITS_toPlot <- rbind(venn_site_ITS_toPlot, venn_site_ITS_sub[[i]])
}
venn_site_ITS_toPlot <- venn_site_ITS_toPlot %>%
  left_join(., input_filt_rare_ITS_n2$taxonomy_loaded, by = c("ASV_ID" = "taxonomy8"))
top_taxa_ITS <- as.data.frame(table(venn_site_ITS_toPlot$taxonomy3)) %>%
  arrange(desc(Freq)) %>%
  slice_head(n = 12) %>%
  mutate(Var1 = as.character(Var1))
venn_site_ITS_toPlot <- venn_site_ITS_toPlot %>%
  group_by(Site, taxonomy3) %>%
  summarise(n_taxa = n()) %>%
  drop_na(taxonomy3)
venn_site_ITS_toPlot$taxonomy3[venn_site_ITS_toPlot$taxonomy3 %notin% top_taxa_ITS$Var1] = 'Other'
venn_site_ITS_toPlot <- venn_site_ITS_toPlot %>%
  group_by(Site, taxonomy3) %>%
  summarise(n_taxa = sum(n_taxa)) %>%
  mutate(taxonomy3 = factor(taxonomy3,
                            levels = c("Other", "NA",
                                       rev(subset(top_taxa_ITS$Var1, top_taxa_ITS$Var1 != "NA")))))
n_per_site_ITS_venn <- venn_site_ITS_toPlot %>%
  group_by(Site) %>%
  summarise(site_n = sum(n_taxa))
pdf("InitialFigs/Venn_Site_ITS_n2.pdf", width = 7, height = 5)
ggplot(venn_site_ITS_toPlot, aes(Site, n_taxa, fill = taxonomy3, group = taxonomy3)) +
  geom_bar(stat = "identity", colour = NA, linewidth = 0.25) +
  labs(x = "Site",
       y = "Number of ASVs",
       fill = "Class") +
  scale_fill_manual(values = c("grey90", "grey75", brewer.pal(11, "Paired")[11:1])) +
  scale_y_continuous(expand = c(0.01, 0.01)) +  
  theme_classic() +
  theme(axis.title.y = element_text(size = 12),
        axis.title.x = element_blank(),
        axis.text.y = element_text(size = 10),
        axis.text.x = element_text(size = 10, angle = 45, hjust = 1))
dev.off()

figs8b <- ggplot(venn_site_ITS_toPlot, aes(Site, n_taxa, fill = taxonomy3, group = taxonomy3)) +
  geom_bar(stat = "identity", colour = NA, linewidth = 0.25) +
  labs(x = "Site",
       y = "# of ASVs",
       fill = "Class") +
  scale_fill_manual(values = c("grey90", "grey75", brewer.pal(11, "Paired")[11:1])) +
  scale_y_continuous(expand = c(0.01, 0.01)) +  
  theme_classic() +
  theme(axis.title.y = element_text(size = 12),
        axis.title.x = element_blank(),
        axis.text.y = element_text(size = 10),
        axis.text.x = element_text(size = 10, angle = 45, hjust = 1),
        legend.margin = margin(0, 0, 0, 0, unit = "pt"),
        legend.box.margin = margin(0, 0, 0, 0, unit = "pt"),
        legend.key.size = unit(0.3, "cm"),
        legend.text = element_text(size = 6))
figs8b
legs8b <- get_legend(figs8b)
figs8b <- figs8b + theme(legend.position = "none")

left <- plot_grid(figs8a, figs8b, ncol = 1, align = "v", rel_heights = c(0.45, 0.55),
                  labels = "auto")
right <- plot_grid(legs8a, legs8b, NULL, ncol = 1, align = "v", rel_heights = c(0.4, 0.5, 0.1))
pdf("FinalFigs/FigureS8.pdf", width = 7, height = 6)
plot_grid(left, right, rel_widths = c(0.8, 0.2))
dev.off()


# OTU level, n = 1 included
otu_tax_ITS <- input_filt_rare_ITS$taxonomy_loaded %>%
  group_by(OTU_ID) %>%
  slice_head(n = 1) %>%
  as.data.frame()
rownames(otu_tax_ITS) <- otu_tax_ITS$OTU_ID
otu_ITS <- otu_table(tax_sum_OTU_ITS, taxa_are_rows = T)
tax_ITS <- tax_table(as.matrix(otu_tax_ITS))
map_ITS <- sample_data(input_filt_rare_ITS$map_loaded)
input_phy_ITS <- phyloseq(otu_ITS, tax_ITS, map_ITS)
venn_site_ITS <- ps_venn(ps = input_phy_ITS,
                         group = "Site",
                         fraction = 0.000000000000000000001,
                         weight = FALSE,
                         relative = FALSE,
                         plot = FALSE)
#saveRDS(venn_site_ITS, "data/venn_site_OTU_ITS.rds") # Slow, so save.
length(venn_site_ITS)
length(venn_site_ITS[[16]]) # 32 OTU present once at all sites
names(venn_site_ITS)[1:16]
venn_site_ITS_sub <- venn_site_ITS[1:16]
site_names <- c(levels(input_filt_rare_ITS$map_loaded$Site), "All")
for (i in 1:length(venn_site_ITS_sub)) {
  venn_site_ITS_sub[[i]] <- as.data.frame(venn_site_ITS_sub[[i]]) %>%
    set_names(c("ASV_ID")) %>%
    mutate(Site = site_names[i])
}
venn_site_ITS_toPlot <- venn_site_ITS_sub[[1]]
for (i in 2:length(venn_site_ITS_sub)) {
  venn_site_ITS_toPlot <- rbind(venn_site_ITS_toPlot, venn_site_ITS_sub[[i]])
}
venn_site_ITS_toPlot <- venn_site_ITS_toPlot %>%
  left_join(., otu_tax_ITS, by = c("ASV_ID" = "OTU_ID"))
top_taxa_ITS <- as.data.frame(table(venn_site_ITS_toPlot$taxonomy3)) %>%
  arrange(desc(Freq)) %>%
  slice_head(n = 12) %>%
  mutate(Var1 = as.character(Var1))
venn_site_ITS_toPlot <- venn_site_ITS_toPlot %>%
  group_by(Site, taxonomy3) %>%
  summarise(n_taxa = n()) %>%
  drop_na(taxonomy3)
venn_site_ITS_toPlot$taxonomy3[venn_site_ITS_toPlot$taxonomy3 %notin% top_taxa_ITS$Var1] = 'Other'
venn_site_ITS_toPlot <- venn_site_ITS_toPlot %>%
  group_by(Site, taxonomy3) %>%
  summarise(n_taxa = sum(n_taxa)) %>%
  mutate(taxonomy3 = factor(taxonomy3,
                            levels = c("Other", "NA",
                                       rev(subset(top_taxa_ITS$Var1, top_taxa_ITS$Var1 != "NA")))))
pdf("InitialFigs/Venn_Site_OTU_ITS.pdf", width = 7, height = 5)
ggplot(venn_site_ITS_toPlot, aes(Site, n_taxa, fill = taxonomy3, group = taxonomy3)) +
  geom_bar(stat = "identity", colour = NA, linewidth = 0.25) +
  labs(x = "Site",
       y = "Number of OTUs",
       fill = "Class") +
  scale_fill_manual(values = c("grey90", "grey75", brewer.pal(11, "Paired")[11:1])) +
  scale_y_continuous(expand = c(0.01, 0.01)) +  
  theme_classic() +
  theme(axis.title.y = element_text(size = 12),
        axis.title.x = element_blank(),
        axis.text.y = element_text(size = 10),
        axis.text.x = element_text(size = 10, angle = 45, hjust = 1))
dev.off()


# OTU level with n = 2 cutoff
otu_tax_ITS_n2 <- input_filt_rare_ITS_n2$taxonomy_loaded %>%
  group_by(OTU_ID) %>%
  slice_head(n = 1) %>%
  as.data.frame()
rownames(otu_tax_ITS_n2) <- otu_tax_ITS_n2$OTU_ID
otu_ITS <- otu_table(tax_sum_OTU_ITS_n2, taxa_are_rows = T)
tax_ITS <- tax_table(as.matrix(otu_tax_ITS_n2))
map_ITS <- sample_data(input_filt_rare_ITS_n2$map_loaded)
input_phy_ITS <- phyloseq(otu_ITS, tax_ITS, map_ITS)
venn_site_ITS <- ps_venn(ps = input_phy_ITS,
                         group = "Site",
                         fraction = 0.000000000000000000001,
                         weight = FALSE,
                         relative = FALSE,
                         plot = FALSE)
#saveRDS(venn_site_ITS, "data/venn_site_OTU_ITS.rds") # Slow, so save.
length(venn_site_ITS)
length(venn_site_ITS[[16]]) # 32 OTU present once at all sites
names(venn_site_ITS)[1:16]
venn_site_ITS_sub <- venn_site_ITS[1:16]
site_names <- c(levels(input_filt_rare_ITS_n2$map_loaded$Site), "All")
for (i in 1:length(venn_site_ITS_sub)) {
  venn_site_ITS_sub[[i]] <- as.data.frame(venn_site_ITS_sub[[i]]) %>%
    set_names(c("ASV_ID")) %>%
    mutate(Site = site_names[i])
}
venn_site_ITS_toPlot <- venn_site_ITS_sub[[1]]
for (i in 2:length(venn_site_ITS_sub)) {
  venn_site_ITS_toPlot <- rbind(venn_site_ITS_toPlot, venn_site_ITS_sub[[i]])
}
venn_site_ITS_toPlot <- venn_site_ITS_toPlot %>%
  left_join(., otu_tax_ITS_n2, by = c("ASV_ID" = "OTU_ID"))
top_taxa_ITS <- as.data.frame(table(venn_site_ITS_toPlot$taxonomy3)) %>%
  arrange(desc(Freq)) %>%
  slice_head(n = 12) %>%
  mutate(Var1 = as.character(Var1))
venn_site_ITS_toPlot <- venn_site_ITS_toPlot %>%
  group_by(Site, taxonomy3) %>%
  summarise(n_taxa = n()) %>%
  drop_na(taxonomy3)
venn_site_ITS_toPlot$taxonomy3[venn_site_ITS_toPlot$taxonomy3 %notin% top_taxa_ITS$Var1] = 'Other'
venn_site_ITS_toPlot <- venn_site_ITS_toPlot %>%
  group_by(Site, taxonomy3) %>%
  summarise(n_taxa = sum(n_taxa)) %>%
  mutate(taxonomy3 = factor(taxonomy3,
                            levels = c("Other", "NA",
                                       rev(subset(top_taxa_ITS$Var1, top_taxa_ITS$Var1 != "NA")))))
pdf("InitialFigs/Venn_Site_OTU_ITS_n2.pdf", width = 7, height = 5)
ggplot(venn_site_ITS_toPlot, aes(Site, n_taxa, fill = taxonomy3, group = taxonomy3)) +
  geom_bar(stat = "identity", colour = NA, linewidth = 0.25) +
  labs(x = "Site",
       y = "Number of OTUs",
       fill = "Class") +
  scale_fill_manual(values = c("grey90", "grey75", brewer.pal(11, "Paired")[11:1])) +
  scale_y_continuous(expand = c(0.01, 0.01)) +  
  theme_classic() +
  theme(axis.title.y = element_text(size = 12),
        axis.title.x = element_blank(),
        axis.text.y = element_text(size = 10),
        axis.text.x = element_text(size = 10, angle = 45, hjust = 1))
dev.off()



#### _Genotype ####
# Note: Genotype analysis runs much faster than site because only 10 levels vs. 15

#### __16S ####
# ASV level, n = 1 included
otu_16S <- otu_table(input_filt_rare_16S$data_loaded, taxa_are_rows = T)
tax_16S <- tax_table(as.matrix(input_filt_rare_16S$taxonomy_loaded))
map_16S <- sample_data(input_filt_rare_16S$map_loaded)
input_phy_16S <- phyloseq(otu_16S, tax_16S, map_16S)
venn_entry_16S <- ps_venn(ps = input_phy_16S,
                         group = "Entry",
                         fraction = 0.000000000000000000001,
                         weight = FALSE,
                         relative = FALSE,
                         plot = FALSE)
length(venn_entry_16S)
length(venn_entry_16S[[11]]) # 9330 ASV present once at all entries
names(venn_entry_16S)[1:11]
venn_entry_16S_sub <- venn_entry_16S[1:11]
entry_names <- c(levels(input_filt_rare_16S$map_loaded$Entry), "All")
for (i in 1:length(venn_entry_16S_sub)) {
  venn_entry_16S_sub[[i]] <- as.data.frame(venn_entry_16S_sub[[i]]) %>%
    set_names(c("ASV_ID")) %>%
    mutate(Entry = entry_names[i])
}
venn_entry_16S_toPlot <- venn_entry_16S_sub[[1]]
for (i in 2:length(venn_entry_16S_sub)) {
  venn_entry_16S_toPlot <- rbind(venn_entry_16S_toPlot, venn_entry_16S_sub[[i]])
}
venn_entry_16S_toPlot <- venn_entry_16S_toPlot %>%
  left_join(., input_filt_rare_16S$taxonomy_loaded, by = c("ASV_ID" = "taxonomy7"))
top_taxa_16S <- as.data.frame(table(venn_entry_16S_toPlot$taxonomy2)) %>%
  arrange(desc(Freq)) %>%
  slice_head(n = 12) %>%
  mutate(Var1 = as.character(Var1))
venn_entry_16S_toPlot <- venn_entry_16S_toPlot %>%
  group_by(Entry, taxonomy2) %>%
  summarise(n_taxa = n()) %>%
  drop_na(taxonomy2)
venn_entry_16S_toPlot$taxonomy2[venn_entry_16S_toPlot$taxonomy2 %notin% top_taxa$Var1] = 'Other'
venn_entry_16S_toPlot <- venn_entry_16S_toPlot %>%
  group_by(Entry, taxonomy2) %>%
  summarise(n_taxa = sum(n_taxa)) %>%
  mutate(taxonomy2 = factor(taxonomy2,
                            levels = c("Other", rev(top_taxa$Var1))))
pdf("InitialFigs/Venn_Entry_16S.pdf", width = 7, height = 5)
ggplot(venn_entry_16S_toPlot, aes(Entry, n_taxa, fill = taxonomy2, group = taxonomy2)) +
  geom_bar(stat = "identity", colour = NA, linewidth = 0.25) +
  labs(x = "Entry",
       y = "Number of ASVs",
       fill = "Phylum") +
  scale_fill_manual(values = c("grey75", brewer.pal(12, "Paired")[12:1])) +
  scale_y_continuous(expand = c(0.01, 0.01)) +  
  theme_classic() +
  theme(axis.title.y = element_text(size = 12),
        axis.title.x = element_blank(),
        axis.text.y = element_text(size = 10),
        axis.text.x = element_text(size = 10, angle = 45, hjust = 1))
dev.off()



# ASV level, n = 2 cutoff
otu_16S <- otu_table(input_filt_rare_16S_n2$data_loaded, taxa_are_rows = T)
tax_16S <- tax_table(as.matrix(input_filt_rare_16S_n2$taxonomy_loaded))
map_16S <- sample_data(input_filt_rare_16S_n2$map_loaded)
input_phy_16S <- phyloseq(otu_16S, tax_16S, map_16S)
venn_entry_16S <- ps_venn(ps = input_phy_16S,
                          group = "Entry",
                          fraction = 0.000000000000000000001,
                          weight = FALSE,
                          relative = FALSE,
                          plot = FALSE)
length(venn_entry_16S)
length(venn_entry_16S[[11]]) # 9330 ASV present once at all entries (~35%)
names(venn_entry_16S)[1:11]
venn_entry_16S_sub <- venn_entry_16S[1:11]
entry_names <- c(levels(input_filt_rare_16S_n2$map_loaded$Entry), "All")
for (i in 1:length(venn_entry_16S_sub)) {
  venn_entry_16S_sub[[i]] <- as.data.frame(venn_entry_16S_sub[[i]]) %>%
    set_names(c("ASV_ID")) %>%
    mutate(Entry = entry_names[i])
}
venn_entry_16S_toPlot <- venn_entry_16S_sub[[1]]
for (i in 2:length(venn_entry_16S_sub)) {
  venn_entry_16S_toPlot <- rbind(venn_entry_16S_toPlot, venn_entry_16S_sub[[i]])
}
venn_entry_16S_toPlot <- venn_entry_16S_toPlot %>%
  left_join(., input_filt_rare_16S_n2$taxonomy_loaded, by = c("ASV_ID" = "taxonomy7"))
top_taxa_16S <- as.data.frame(table(venn_entry_16S_toPlot$taxonomy2)) %>%
  arrange(desc(Freq)) %>%
  slice_head(n = 12) %>%
  mutate(Var1 = as.character(Var1))
venn_entry_16S_toPlot <- venn_entry_16S_toPlot %>%
  group_by(Entry, taxonomy2) %>%
  summarise(n_taxa = n()) %>%
  drop_na(taxonomy2)
venn_entry_16S_toPlot$taxonomy2[venn_entry_16S_toPlot$taxonomy2 %notin% top_taxa_16S$Var1] = 'Other'
venn_entry_16S_toPlot <- venn_entry_16S_toPlot %>%
  group_by(Entry, taxonomy2) %>%
  summarise(n_taxa = sum(n_taxa)) %>%
  mutate(taxonomy2 = factor(taxonomy2,
                            levels = c("Other", rev(top_taxa_16S$Var1))))
venn_text_sum_16S <- venn_entry_16S_toPlot %>%
  group_by(Entry) %>%
  summarise(tot = sum(n_taxa)) %>%
  mutate(y = c(9500, rep(200, 10)))
pdf("InitialFigs/Venn_Entry_16S_n2.pdf", width = 7, height = 5)
ggplot(venn_entry_16S_toPlot, aes(Entry, n_taxa, fill = taxonomy2, group = taxonomy2)) +
  geom_bar(stat = "identity", colour = NA, linewidth = 0.25) +
  geom_text(data = venn_text_sum_16S,
            aes(Entry, y, label = tot), inherit.aes = F) +
  labs(x = "Entry",
       y = "Number of ASVs",
       fill = "Phylum") +
  scale_fill_manual(values = c("grey75", brewer.pal(12, "Paired")[12:1])) +
  scale_y_continuous(expand = c(0.02, 0.02)) +  
  theme_classic() +
  theme(axis.title.y = element_text(size = 12),
        axis.title.x = element_blank(),
        axis.text.y = element_text(size = 10),
        axis.text.x = element_text(size = 10, angle = 45, hjust = 1))
dev.off()



# OTU level, n = 1 included
venn_entry_16S <- ps_venn(ps = input_phy_16S,
                          group = "Entry",
                          fraction = 0.000000000000000000001,
                          weight = FALSE,
                          relative = FALSE,
                          plot = FALSE)
length(venn_entry_16S)
length(venn_entry_16S[[11]]) # 5172 OTU present once at all entries
names(venn_entry_16S)[1:11]
venn_entry_16S_sub <- venn_entry_16S[1:11]
entry_names <- c(levels(input_filt_rare_16S$map_loaded$Entry), "All")
for (i in 1:length(venn_entry_16S_sub)) {
  venn_entry_16S_sub[[i]] <- as.data.frame(venn_entry_16S_sub[[i]]) %>%
    set_names(c("ASV_ID")) %>%
    mutate(Entry = entry_names[i])
}
venn_entry_16S_toPlot <- venn_entry_16S_sub[[1]]
for (i in 2:length(venn_entry_16S_sub)) {
  venn_entry_16S_toPlot <- rbind(venn_entry_16S_toPlot, venn_entry_16S_sub[[i]])
}
venn_entry_16S_toPlot <- venn_entry_16S_toPlot %>%
  left_join(., otu_tax_16S, by = c("ASV_ID" = "OTU_ID"))
top_taxa_16S <- as.data.frame(table(venn_entry_16S_toPlot$taxonomy2)) %>%
  arrange(desc(Freq)) %>%
  slice_head(n = 12) %>%
  mutate(Var1 = as.character(Var1))
venn_entry_16S_toPlot <- venn_entry_16S_toPlot %>%
  group_by(Entry, taxonomy2) %>%
  summarise(n_taxa = n()) %>%
  drop_na(taxonomy2)
venn_entry_16S_toPlot$taxonomy2[venn_entry_16S_toPlot$taxonomy2 %notin% top_taxa_16S$Var1] = 'Other'
venn_entry_16S_toPlot <- venn_entry_16S_toPlot %>%
  group_by(Entry, taxonomy2) %>%
  summarise(n_taxa = sum(n_taxa)) %>%
  mutate(taxonomy2 = factor(taxonomy2,
                            levels = c("Other", rev(top_taxa_16S$Var1))))
venn_text_sum_16S <- venn_entry_16S_toPlot %>%
  group_by(Entry) %>%
  summarise(tot = sum(n_taxa)) %>%
  mutate(y = c(5300, rep(200, 10)))
pdf("InitialFigs/Venn_Entry_OTU_16S.pdf", width = 7, height = 5)
ggplot(venn_entry_16S_toPlot, aes(Entry, n_taxa, fill = taxonomy2, group = taxonomy2)) +
  geom_bar(stat = "identity", colour = NA, linewidth = 0.25) +
  geom_text(data = venn_text_sum_16S,
            aes(Entry, y, label = tot), inherit.aes = F) +
  labs(x = "Entry",
       y = "Number of OTUs",
       fill = "Phylum") +
  scale_fill_manual(values = c("grey75", brewer.pal(12, "Paired")[12:1])) +
  scale_y_continuous(expand = c(0.02, 0.02)) +  
  theme_classic() +
  theme(axis.title.y = element_text(size = 12),
        axis.title.x = element_blank(),
        axis.text.y = element_text(size = 10),
        axis.text.x = element_text(size = 10, angle = 45, hjust = 1))
dev.off()



# OTU level, n = 2 cutoff
otu_16S <- otu_table(tax_sum_OTU_16S_n2, taxa_are_rows = T)
tax_16S <- tax_table(as.matrix(otu_tax_16S_n2))
map_16S <- sample_data(input_filt_rare_16S_n2$map_loaded)
input_phy_16S <- phyloseq(otu_16S, tax_16S, map_16S)
venn_entry_16S <- ps_venn(ps = input_phy_16S,
                          group = "Entry",
                          fraction = 0.000000000000000000001,
                          weight = FALSE,
                          relative = FALSE,
                          plot = FALSE)
length(venn_entry_16S)
length(venn_entry_16S[[11]]) # 5165 OTU present once at all entries
names(venn_entry_16S)[1:11]
venn_entry_16S_sub <- venn_entry_16S[1:11]
entry_names <- c(levels(input_filt_rare_16S_n2$map_loaded$Entry), "All")
for (i in 1:length(venn_entry_16S_sub)) {
  venn_entry_16S_sub[[i]] <- as.data.frame(venn_entry_16S_sub[[i]]) %>%
    set_names(c("ASV_ID")) %>%
    mutate(Entry = entry_names[i])
}
venn_entry_16S_toPlot <- venn_entry_16S_sub[[1]]
for (i in 2:length(venn_entry_16S_sub)) {
  venn_entry_16S_toPlot <- rbind(venn_entry_16S_toPlot, venn_entry_16S_sub[[i]])
}
venn_entry_16S_toPlot <- venn_entry_16S_toPlot %>%
  left_join(., otu_tax_16S_n2, by = c("ASV_ID" = "OTU_ID"))
top_taxa_16S <- as.data.frame(table(venn_entry_16S_toPlot$taxonomy2)) %>%
  arrange(desc(Freq)) %>%
  slice_head(n = 12) %>%
  mutate(Var1 = as.character(Var1))
venn_entry_16S_toPlot <- venn_entry_16S_toPlot %>%
  group_by(Entry, taxonomy2) %>%
  summarise(n_taxa = n()) %>%
  drop_na(taxonomy2)
venn_entry_16S_toPlot$taxonomy2[venn_entry_16S_toPlot$taxonomy2 %notin% top_taxa_16S$Var1] = 'Other'
venn_entry_16S_toPlot <- venn_entry_16S_toPlot %>%
  group_by(Entry, taxonomy2) %>%
  summarise(n_taxa = sum(n_taxa)) %>%
  mutate(taxonomy2 = factor(taxonomy2,
                            levels = c("Other", rev(top_taxa_16S$Var1))))
venn_text_sum_16S <- venn_entry_16S_toPlot %>%
  group_by(Entry) %>%
  summarise(tot = sum(n_taxa)) %>%
  mutate(y = c(5300, rep(200, 10)))
pdf("InitialFigs/Venn_Entry_OTU_16S_n2.pdf", width = 7, height = 5)
ggplot(venn_entry_16S_toPlot, aes(Entry, n_taxa, fill = taxonomy2, group = taxonomy2)) +
  geom_bar(stat = "identity", colour = NA, linewidth = 0.25) +
  geom_text(data = venn_text_sum_16S,
            aes(Entry, y, label = tot), inherit.aes = F) +
  labs(x = "Entry",
       y = "Number of OTUs",
       fill = "Phylum") +
  scale_fill_manual(values = c("grey75", brewer.pal(12, "Paired")[12:1])) +
  scale_y_continuous(expand = c(0.02, 0.02)) +  
  theme_classic() +
  theme(axis.title.y = element_text(size = 12),
        axis.title.x = element_blank(),
        axis.text.y = element_text(size = 10),
        axis.text.x = element_text(size = 10, angle = 45, hjust = 1))
dev.off()



#### __ITS ####
# ASV level, n = 1 included
venn_entry_ITS <- ps_venn(ps = input_phy_ITS,
                         group = "Entry",
                         fraction = 0.000000000000000000001,
                         weight = FALSE,
                         relative = FALSE,
                         plot = FALSE)
length(venn_entry_ITS)
length(venn_entry_ITS[[11]]) # 606 ASV present once at all entries
names(venn_entry_ITS)[1:11]
venn_entry_ITS_sub <- venn_entry_ITS[1:11]
entry_names <- c(levels(input_filt_rare_ITS$map_loaded$Entry), "All")
for (i in 1:length(venn_entry_ITS_sub)) {
  venn_entry_ITS_sub[[i]] <- as.data.frame(venn_entry_ITS_sub[[i]]) %>%
    set_names(c("ASV_ID")) %>%
    mutate(Entry = entry_names[i])
}
venn_entry_ITS_toPlot <- venn_entry_ITS_sub[[1]]
for (i in 2:length(venn_entry_ITS_sub)) {
  venn_entry_ITS_toPlot <- rbind(venn_entry_ITS_toPlot, venn_entry_ITS_sub[[i]])
}
venn_entry_ITS_toPlot <- venn_entry_ITS_toPlot %>%
  left_join(., input_filt_rare_ITS$taxonomy_loaded, by = c("ASV_ID" = "taxonomy8"))
top_taxa_ITS <- as.data.frame(table(venn_entry_ITS_toPlot$taxonomy3)) %>%
  arrange(desc(Freq)) %>%
  slice_head(n = 12) %>%
  mutate(Var1 = as.character(Var1))
venn_entry_ITS_toPlot <- venn_entry_ITS_toPlot %>%
  group_by(Entry, taxonomy3) %>%
  summarise(n_taxa = n()) %>%
  drop_na(taxonomy3)
venn_entry_ITS_toPlot$taxonomy3[venn_entry_ITS_toPlot$taxonomy3 %notin% top_taxa_ITS$Var1] = 'Other'
venn_entry_ITS_toPlot <- venn_entry_ITS_toPlot %>%
  group_by(Entry, taxonomy3) %>%
  summarise(n_taxa = sum(n_taxa)) %>%
  mutate(taxonomy3 = factor(taxonomy3,
                            levels = c("Other", "NA",
                                       rev(subset(top_taxa_ITS$Var1, top_taxa_ITS$Var1 != "NA")))))
pdf("InitialFigs/Venn_Entry_ITS.pdf", width = 7, height = 5)
ggplot(venn_entry_ITS_toPlot, aes(Entry, n_taxa, fill = taxonomy3, group = taxonomy3)) +
  geom_bar(stat = "identity", colour = NA, linewidth = 0.25) +
  labs(x = "Genotype",
       y = "Number of ASVs",
       fill = "Class") +
  scale_fill_manual(values = c("grey90", "grey75", brewer.pal(11, "Paired")[11:1])) +
  scale_y_continuous(expand = c(0.01, 0.01)) +  
  theme_classic() +
  theme(axis.title.y = element_text(size = 12),
        axis.title.x = element_blank(),
        axis.text.y = element_text(size = 10),
        axis.text.x = element_text(size = 10, angle = 45, hjust = 1))
dev.off()



# ASV level, n = 2 cutoff
otu_ITS <- otu_table(input_filt_rare_ITS_n2$data_loaded, taxa_are_rows = T)
tax_ITS <- tax_table(as.matrix(input_filt_rare_ITS_n2$taxonomy_loaded))
map_ITS <- sample_data(input_filt_rare_ITS_n2$map_loaded)
input_phy_ITS <- phyloseq(otu_ITS, tax_ITS, map_ITS)
venn_entry_ITS <- ps_venn(ps = input_phy_ITS,
                          group = "Entry",
                          fraction = 0.000000000000000000001,
                          weight = FALSE,
                          relative = FALSE,
                          plot = FALSE)
length(venn_entry_ITS)
length(venn_entry_ITS[[11]]) # 606 ASV present once at all entries
names(venn_entry_ITS)[1:11]
venn_entry_ITS_sub <- venn_entry_ITS[1:11]
entry_names <- c(levels(input_filt_rare_ITS_n2$map_loaded$Entry), "All")
for (i in 1:length(venn_entry_ITS_sub)) {
  venn_entry_ITS_sub[[i]] <- as.data.frame(venn_entry_ITS_sub[[i]]) %>%
    set_names(c("ASV_ID")) %>%
    mutate(Entry = entry_names[i])
}
venn_entry_ITS_toPlot <- venn_entry_ITS_sub[[1]]
for (i in 2:length(venn_entry_ITS_sub)) {
  venn_entry_ITS_toPlot <- rbind(venn_entry_ITS_toPlot, venn_entry_ITS_sub[[i]])
}
venn_entry_ITS_toPlot <- venn_entry_ITS_toPlot %>%
  left_join(., input_filt_rare_ITS_n2$taxonomy_loaded, by = c("ASV_ID" = "taxonomy8"))
top_taxa_ITS <- as.data.frame(table(venn_entry_ITS_toPlot$taxonomy3)) %>%
  arrange(desc(Freq)) %>%
  slice_head(n = 12) %>%
  mutate(Var1 = as.character(Var1))
venn_entry_ITS_toPlot <- venn_entry_ITS_toPlot %>%
  group_by(Entry, taxonomy3) %>%
  summarise(n_taxa = n()) %>%
  drop_na(taxonomy3)
venn_entry_ITS_toPlot$taxonomy3[venn_entry_ITS_toPlot$taxonomy3 %notin% top_taxa_ITS$Var1] = 'Other'
venn_entry_ITS_toPlot <- venn_entry_ITS_toPlot %>%
  group_by(Entry, taxonomy3) %>%
  summarise(n_taxa = sum(n_taxa)) %>%
  mutate(taxonomy3 = factor(taxonomy3,
                            levels = c("Other", "NA",
                                       rev(subset(top_taxa_ITS$Var1, top_taxa_ITS$Var1 != "NA")))))
venn_text_sum_ITS <- venn_entry_ITS_toPlot %>%
  group_by(Entry) %>%
  summarise(tot = sum(n_taxa)) %>%
  mutate(y = c(620, rep(30, 10)))
pdf("InitialFigs/Venn_Entry_ITS_n2.pdf", width = 7, height = 5)
ggplot(venn_entry_ITS_toPlot, aes(Entry, n_taxa, fill = taxonomy3, group = taxonomy3)) +
  geom_bar(stat = "identity", colour = NA, linewidth = 0.25) +
  geom_text(data = venn_text_sum_ITS,
            aes(Entry, y, label = tot), inherit.aes = F) +
  labs(x = "Genotype",
       y = "Number of ASVs",
       fill = "Class") +
  scale_fill_manual(values = c("grey90", "grey75", brewer.pal(11, "Paired")[11:1])) +
  scale_y_continuous(expand = c(0.02, 0.02)) +  
  theme_classic() +
  theme(axis.title.y = element_text(size = 12),
        axis.title.x = element_blank(),
        axis.text.y = element_text(size = 10),
        axis.text.x = element_text(size = 10, angle = 45, hjust = 1))
dev.off()



# OTU level, n = 1 included
otu_tax_ITS <- input_filt_rare_ITS$taxonomy_loaded %>%
  group_by(OTU_ID) %>%
  slice_head(n = 1) %>%
  as.data.frame()
rownames(otu_tax_ITS) <- otu_tax_ITS$OTU_ID
otu_ITS <- otu_table(tax_sum_OTU_ITS, taxa_are_rows = T)
tax_ITS <- tax_table(as.matrix(otu_tax_ITS))
map_ITS <- sample_data(input_filt_rare_ITS$map_loaded)
input_phy_ITS <- phyloseq(otu_ITS, tax_ITS, map_ITS)
venn_entry_ITS <- ps_venn(ps = input_phy_ITS,
                          group = "Entry",
                          fraction = 0.000000000000000000001,
                          weight = FALSE,
                          relative = FALSE,
                          plot = FALSE)
length(venn_entry_ITS)
length(venn_entry_ITS[[11]]) # 553 ASV present once at all entries
names(venn_entry_ITS)[1:11]
venn_entry_ITS_sub <- venn_entry_ITS[1:11]
entry_names <- c(levels(input_filt_rare_ITS$map_loaded$Entry), "All")
for (i in 1:length(venn_entry_ITS_sub)) {
  venn_entry_ITS_sub[[i]] <- as.data.frame(venn_entry_ITS_sub[[i]]) %>%
    set_names(c("ASV_ID")) %>%
    mutate(Entry = entry_names[i])
}
venn_entry_ITS_toPlot <- venn_entry_ITS_sub[[1]]
for (i in 2:length(venn_entry_ITS_sub)) {
  venn_entry_ITS_toPlot <- rbind(venn_entry_ITS_toPlot, venn_entry_ITS_sub[[i]])
}
venn_entry_ITS_toPlot <- venn_entry_ITS_toPlot %>%
  left_join(., otu_tax_ITS, by = c("ASV_ID" = "OTU_ID"))
top_taxa_ITS <- as.data.frame(table(venn_entry_ITS_toPlot$taxonomy3)) %>%
  arrange(desc(Freq)) %>%
  slice_head(n = 12) %>%
  mutate(Var1 = as.character(Var1))
venn_entry_ITS_toPlot <- venn_entry_ITS_toPlot %>%
  group_by(Entry, taxonomy3) %>%
  summarise(n_taxa = n()) %>%
  drop_na(taxonomy3)
venn_entry_ITS_toPlot$taxonomy3[venn_entry_ITS_toPlot$taxonomy3 %notin% top_taxa_ITS$Var1] = 'Other'
venn_entry_ITS_toPlot <- venn_entry_ITS_toPlot %>%
  group_by(Entry, taxonomy3) %>%
  summarise(n_taxa = sum(n_taxa)) %>%
  mutate(taxonomy3 = factor(taxonomy3,
                            levels = c("Other", "NA",
                                       rev(subset(top_taxa_ITS$Var1, top_taxa_ITS$Var1 != "NA")))))
venn_text_sum_ITS <- venn_entry_ITS_toPlot %>%
  group_by(Entry) %>%
  summarise(tot = sum(n_taxa)) %>%
  mutate(y = c(575, rep(100, 10)))
pdf("InitialFigs/Venn_Entry_OTU_ITS.pdf", width = 7, height = 5)
ggplot(venn_entry_ITS_toPlot, aes(Entry, n_taxa, fill = taxonomy3, group = taxonomy3)) +
  geom_bar(stat = "identity", colour = NA, linewidth = 0.25) +
  geom_text(data = venn_text_sum_ITS,
            aes(Entry, y, label = tot), inherit.aes = F) +
  labs(x = "Genotype",
       y = "Number of OTUs",
       fill = "Class") +
  scale_fill_manual(values = c("grey90", "grey75", brewer.pal(11, "Paired")[11:1])) +
  scale_y_continuous(expand = c(0.02, 0.02)) +  
  theme_classic() +
  theme(axis.title.y = element_text(size = 12),
        axis.title.x = element_blank(),
        axis.text.y = element_text(size = 10),
        axis.text.x = element_text(size = 10, angle = 45, hjust = 1))
dev.off()



# OTU level, n = 2 cutoff
otu_tax_ITS_n2 <- input_filt_rare_ITS_n2$taxonomy_loaded %>%
  group_by(OTU_ID) %>%
  slice_head(n = 1) %>%
  as.data.frame()
rownames(otu_tax_ITS_n2) <- otu_tax_ITS_n2$OTU_ID
otu_ITS <- otu_table(tax_sum_OTU_ITS_n2, taxa_are_rows = T)
tax_ITS <- tax_table(as.matrix(otu_tax_ITS_n2))
map_ITS <- sample_data(input_filt_rare_ITS_n2$map_loaded)
input_phy_ITS <- phyloseq(otu_ITS, tax_ITS, map_ITS)
venn_entry_ITS <- ps_venn(ps = input_phy_ITS,
                          group = "Entry",
                          fraction = 0.000000000000000000001,
                          weight = FALSE,
                          relative = FALSE,
                          plot = FALSE)
length(venn_entry_ITS)
ntax <- as.data.frame(matrix(NA, nrow = length(venn_entry_ITS), ncol = 2)) %>%
  set_names(c("Group", "nOTU"))
for (i in 1:length(venn_entry_ITS)) {
  ntax$Group[i] <- names(venn_entry_ITS[i])
  ntax$nOTU[i] <- length(venn_entry_ITS[[i]])
}
length(venn_entry_ITS[[10]]) # 542 OTUs present once at all entries. Note all is 10 not 11! 
names(venn_entry_ITS)[1:11] # PI 650758 missing from first 11. This must be because it had 0. Need to assign 0.
sum(names(venn_entry_ITS) == "PI 650758")
venn_entry_ITS_sub <- venn_entry_ITS[1:10]
entry_names <- c(names(venn_entry_ITS)[1:9], "All")
for (i in 1:length(venn_entry_ITS_sub)) {
  venn_entry_ITS_sub[[i]] <- as.data.frame(venn_entry_ITS_sub[[i]]) %>%
    set_names(c("ASV_ID")) %>%
    mutate(Entry = entry_names[i])
}
venn_entry_ITS_toPlot <- venn_entry_ITS_sub[[1]]
for (i in 2:length(venn_entry_ITS_sub)) {
  venn_entry_ITS_toPlot <- rbind(venn_entry_ITS_toPlot, venn_entry_ITS_sub[[i]])
}
venn_entry_ITS_toPlot <- venn_entry_ITS_toPlot %>%
  left_join(., otu_tax_ITS_n2, by = c("ASV_ID" = "OTU_ID"))
top_taxa_ITS <- as.data.frame(table(venn_entry_ITS_toPlot$taxonomy3)) %>%
  arrange(desc(Freq)) %>%
  slice_head(n = 12) %>%
  mutate(Var1 = as.character(Var1))
venn_entry_ITS_toPlot <- venn_entry_ITS_toPlot %>%
  group_by(Entry, taxonomy3) %>%
  summarise(n_taxa = n()) %>%
  drop_na(taxonomy3)
venn_entry_ITS_toPlot$taxonomy3[venn_entry_ITS_toPlot$taxonomy3 %notin% top_taxa_ITS$Var1] = 'Other'
venn_entry_ITS_toPlot <- venn_entry_ITS_toPlot %>%
  group_by(Entry, taxonomy3) %>%
  summarise(n_taxa = sum(n_taxa)) %>%
  mutate(taxonomy3 = factor(taxonomy3,
                            levels = c("Other", "NA",
                                       rev(subset(top_taxa_ITS$Var1, top_taxa_ITS$Var1 != "NA")))))
venn_entry_ITS_toPlot <- venn_entry_ITS_toPlot %>%
  ungroup() %>%
  add_row(Entry = "PI 650758", taxonomy3 = "Other", n_taxa = 0)
venn_text_sum_ITS <- venn_entry_ITS_toPlot %>%
  group_by(Entry) %>%
  summarise(tot = sum(n_taxa)) %>%
  mutate(y = c(575, rep(30, 10)))
pdf("InitialFigs/Venn_Entry_OTU_ITS_n2.pdf", width = 7, height = 5)
ggplot(venn_entry_ITS_toPlot, aes(Entry, n_taxa, fill = taxonomy3, group = taxonomy3)) +
  geom_bar(stat = "identity", colour = NA, linewidth = 0.25) +
  geom_text(data = venn_text_sum_ITS,
            aes(Entry, y, label = tot), inherit.aes = F) +
  labs(x = "Genotype",
       y = "Number of OTUs",
       fill = "Class") +
  scale_fill_manual(values = c("grey90", "grey75", brewer.pal(11, "Paired")[11:1])) +
  scale_y_continuous(expand = c(0.02, 0.02)) +  
  theme_classic() +
  theme(axis.title.y = element_text(size = 12),
        axis.title.x = element_blank(),
        axis.text.y = element_text(size = 10),
        axis.text.x = element_text(size = 10, angle = 45, hjust = 1))
dev.off()



#### 8. Site Subsets ####
# Make subsets of each site and test effect of genotype
# According to the significant interaction, there should be some w/ sig. effect and others not.
# Record R2 and p values from PERMANOVA. Plot multipanel PCoA
# Also try some combinations of sites with similar composition/few unique taxa
# For example: Brookings, Casselton, Grandin

# Make Site and Entry characters
input_filt_rare_16S$map_loaded$Site <- as.character(input_filt_rare_16S$map_loaded$Site)
input_filt_rare_ITS$map_loaded$Site <- as.character(input_filt_rare_ITS$map_loaded$Site)
input_filt_rare_16S$map_loaded$Entry <- as.character(input_filt_rare_16S$map_loaded$Entry)
input_filt_rare_ITS$map_loaded$Entry <- as.character(input_filt_rare_ITS$map_loaded$Entry)

# 16S
l_16S <- list()
for (i in 1:length(unique(input_filt_rare_16S$map_loaded$Site))) {
  l_16S[[i]] <- filter_data(input_filt_rare_16S, 
                            filter_cat = "Site", 
                            keep_vals = unique(input_filt_rare_16S$map_loaded$Site)[i])
}

l_bc_16S <- list()
for (i in 1:length(l_16S)) {
  l_bc_16S[[i]] <- calc_dm(l_16S[[i]]$data_loaded)
}

res_df_16S <- as.data.frame(matrix(NA, nrow = length(l_16S), ncol = 3)) %>%
  set_names(c("Site", "R2", "p"))
for (i in 1:length(l_16S)) {
  set.seed(1000)
  pm <- adonis2(l_bc_16S[[i]] ~ l_16S[[i]]$map_loaded$Entry)
  res_df_16S$Site[i] <- unique(input_filt_rare_16S$map_loaded$Site)[i]
  res_df_16S$R2[i] <- pm$R2[1]
  res_df_16S$p[i] <- pm$`Pr(>F)`[1]
}

for (i in 1:length(l_16S)) {
  pcoa <- cmdscale(l_bc_16S[[i]], k = nrow(l_16S[[i]]$map_loaded) - 1, eig = T)
  pcoaA1 <- paste("PC1: ", round((eigenvals(pcoa)/sum(eigenvals(pcoa)))[1]*100, 1), "%")
  pcoaA2 <- paste("PC2: ", round((eigenvals(pcoa)/sum(eigenvals(pcoa)))[2]*100, 1), "%")
  l_16S[[i]]$map_loaded$Axis01 <- vegan::scores(pcoa)[,1]
  l_16S[[i]]$map_loaded$Axis02 <- vegan::scores(pcoa)[,2]
}

site_ords_16S <- l_16S[[1]]$map_loaded
for (i in 2:length(l_16S)) {
  site_ords_16S <- rbind(site_ords_16S, l_16S[[i]]$map_loaded)
}

res_df_16S <- res_df_16S %>%
  mutate(Symbol = ifelse(p <= 0.001,
                         "***",
                         ifelse(p > 0.001 & p <= 0.01,
                                "**",
                                ifelse(p > 0.01 & p <= 0.05,
                                       "*",
                                       "N.S."))))
micro.hulls.16S <- ddply(site_ords_16S, c("Site", "Entry"), find_hull)
pdf("InitialFigs/Beta_Bray_Site_16S.pdf", width = 7, height = 7)
g14 <- ggplot(site_ords_16S, aes(Axis01, Axis02)) +
  geom_polygon(data = micro.hulls.16S, 
               aes(colour = Entry),
               alpha = 0, show.legend = F) +
  geom_point(size = 2, pch = 16, alpha = 0.75, aes(colour = Entry)) +
  geom_text(data = subset(res_df_16S, p > 0.05),
            aes(x = 0, y = 0.3, label = Symbol),
            size = 3, inherit.aes = F) +
  geom_text(data = subset(res_df_16S, p <= 0.05),
            aes(x = 0, y = 0.26, label = Symbol),
            size = 8, inherit.aes = F) +
  labs(x = "PC1", 
       y = "PC2",
       colour = "Genotype",
       title = "a) Archaea/Bacteria") +
  scale_colour_brewer(palette = "Paired") +
  scale_y_continuous(expand = c(0.1, 0.05),
                     breaks = c(-0.2, 0, 0.2)) +
  scale_x_continuous(breaks = c(-0.2, 0, 0.2),
                     labels = c("-0.2", "0", "0.2")) +
  facet_wrap(~ Site, ncol = 3) +
  guides(colour = guide_legend(override.aes = list(alpha = 1, size = 3))) +
  theme_bw() +  
  theme(legend.position = "right",
        axis.title = element_text(face = "bold", size = 12), 
        axis.text = element_text(size = 8),
        plot.title = element_text(vjust = -1, hjust = 0.5),
        strip.text = element_text(size = 8),
        strip.background = element_rect(fill = "white"))
g14
dev.off()
leg_gen <- get_legend(g14)
g14 <- g14 + theme(legend.position = "none")



# ITS
l_ITS <- list()
for (i in 1:length(unique(input_filt_rare_ITS$map_loaded$Site))) {
  l_ITS[[i]] <- filter_data(input_filt_rare_ITS, 
                            filter_cat = "Site", 
                            keep_vals = unique(input_filt_rare_ITS$map_loaded$Site)[i])
}

l_bc_ITS <- list()
for (i in 1:length(l_ITS)) {
  l_bc_ITS[[i]] <- calc_dm(l_ITS[[i]]$data_loaded)
}

res_df_ITS <- as.data.frame(matrix(NA, nrow = length(l_ITS), ncol = 3)) %>%
  set_names(c("Site", "R2", "p"))
for (i in 1:length(l_ITS)) {
  set.seed(1000)
  pm <- adonis2(l_bc_ITS[[i]] ~ l_ITS[[i]]$map_loaded$Entry)
  res_df_ITS$Site[i] <- unique(input_filt_rare_ITS$map_loaded$Site)[i]
  res_df_ITS$R2[i] <- pm$R2[1]
  res_df_ITS$p[i] <- pm$`Pr(>F)`[1]
}

for (i in 1:length(l_ITS)) {
  pcoa <- cmdscale(l_bc_ITS[[i]], k = nrow(l_ITS[[i]]$map_loaded) - 1, eig = T)
  pcoaA1 <- paste("PC1: ", round((eigenvals(pcoa)/sum(eigenvals(pcoa)))[1]*100, 1), "%")
  pcoaA2 <- paste("PC2: ", round((eigenvals(pcoa)/sum(eigenvals(pcoa)))[2]*100, 1), "%")
  l_ITS[[i]]$map_loaded$Axis01 <- vegan::scores(pcoa)[,1]
  l_ITS[[i]]$map_loaded$Axis02 <- vegan::scores(pcoa)[,2]
}

site_ords_ITS <- l_ITS[[1]]$map_loaded
for (i in 2:length(l_ITS)) {
  site_ords_ITS <- rbind(site_ords_ITS, l_ITS[[i]]$map_loaded)
}

res_df_ITS <- res_df_ITS %>%
  mutate(Symbol = ifelse(p <= 0.001,
                         "***",
                         ifelse(p > 0.001 & p <= 0.01,
                                "**",
                                ifelse(p > 0.01 & p <= 0.05,
                                       "*",
                                       "N.S."))))
micro.hulls.ITS <- ddply(site_ords_ITS, c("Site", "Entry"), find_hull)
pdf("InitialFigs/Beta_Bray_Site_ITS.pdf", width = 7, height = 7)
g15 <- ggplot(site_ords_ITS, aes(Axis01, Axis02)) +
  geom_polygon(data = micro.hulls.ITS, 
               aes(colour = Entry),
               alpha = 0, show.legend = F) +
  geom_point(size = 2, pch = 16, alpha = 0.75, aes(colour = Entry)) +
  geom_text(data = subset(res_df_ITS, p > 0.05),
            aes(x = 0, y = 0.3, label = Symbol),
            size = 3, inherit.aes = F) +
  geom_text(data = subset(res_df_ITS, p <= 0.05),
            aes(x = 0, y = 0.26, label = Symbol),
            size = 8, inherit.aes = F) +
  labs(x = "PC1", 
       y = "PC2",
       colour = "Genotype",
       title = "b) Fungi") +
  scale_colour_brewer(palette = "Paired") +
  scale_y_continuous(expand = c(0.1, 0.05),
                     breaks = c(-0.2, 0, 0.2)) +
  scale_x_continuous(breaks = c(-0.2, 0, 0.2),
                     labels = c("-0.2", "0", "0.2")) +
  facet_wrap(~ Site, ncol = 3) +
  guides(colour = guide_legend(override.aes = list(alpha = 1, size = 3))) +
  theme_bw() +  
  theme(legend.position = "none",
        axis.title = element_text(face = "bold", size = 12), 
        axis.text = element_text(size = 8),
        plot.title = element_text(vjust = -1, hjust = 0.5),
        strip.text = element_text(size = 8),
        strip.background = element_rect(fill = "white"))
g15
dev.off()

# Brookings, Casselton, Grandin
detach("package:fields", unload=TRUE)
detach("package:spam", unload=TRUE)
BroCasGra_16S <- filter_data(input_filt_rare_16S,
                             filter_cat = "Site",
                             keep_vals = c("Brookings", "Casselton", "Grandin"))
bc_bcg_16S <- calc_dm(BroCasGra_16S$data_loaded)
set.seed(1150)
ado19 <- adonis2(bc_bcg_16S ~ BroCasGra_16S$map_loaded$Site * BroCasGra_16S$map_loaded$Entry)
ado19 # Site***

BroCasGra_ITS <- filter_data(input_filt_rare_ITS,
                             filter_cat = "Site",
                             keep_vals = c("Brookings", "Casselton", "Grandin"))
bc_bcg_ITS <- calc_dm(BroCasGra_ITS$data_loaded)
set.seed(1150)
ado20 <- adonis2(bc_bcg_ITS ~ BroCasGra_ITS$map_loaded$Site * BroCasGra_ITS$map_loaded$Entry)
ado20 # Site***



#### 9. Genotype Subsets ####
# Do the same thing as above but switch Site and Entry
# 16S
l_16S <- list()
for (i in 1:length(unique(input_filt_rare_16S$map_loaded$Entry))) {
  l_16S[[i]] <- filter_data(input_filt_rare_16S, 
                            filter_cat = "Entry", 
                            keep_vals = unique(input_filt_rare_16S$map_loaded$Entry)[i])
}

l_bc_16S <- list()
for (i in 1:length(l_16S)) {
  l_bc_16S[[i]] <- calc_dm(l_16S[[i]]$data_loaded)
}

res_df_16S <- as.data.frame(matrix(NA, nrow = length(l_16S), ncol = 3)) %>%
  set_names(c("Site", "R2", "p"))
for (i in 1:length(l_16S)) {
  set.seed(1000)
  pm <- adonis2(l_bc_16S[[i]] ~ l_16S[[i]]$map_loaded$Site)
  res_df_16S$Site[i] <- unique(input_filt_rare_16S$map_loaded$Site)[i]
  res_df_16S$R2[i] <- pm$R2[1]
  res_df_16S$p[i] <- pm$`Pr(>F)`[1]
}

for (i in 1:length(l_16S)) {
  pcoa <- cmdscale(l_bc_16S[[i]], k = nrow(l_16S[[i]]$map_loaded) - 1, eig = T)
  pcoaA1 <- paste("PC1: ", round((eigenvals(pcoa)/sum(eigenvals(pcoa)))[1]*100, 1), "%")
  pcoaA2 <- paste("PC2: ", round((eigenvals(pcoa)/sum(eigenvals(pcoa)))[2]*100, 1), "%")
  l_16S[[i]]$map_loaded$Axis01 <- vegan::scores(pcoa)[,1]
  l_16S[[i]]$map_loaded$Axis02 <- vegan::scores(pcoa)[,2]
}

Entry_ords_16S <- l_16S[[1]]$map_loaded
for (i in 2:length(l_16S)) {
  Entry_ords_16S <- rbind(Entry_ords_16S, l_16S[[i]]$map_loaded)
}

res_df_16S <- res_df_16S %>%
  mutate(Symbol = ifelse(p <= 0.001,
                         "***",
                         ifelse(p > 0.001 & p <= 0.01,
                                "**",
                                ifelse(p > 0.01 & p <= 0.05,
                                       "*",
                                       "N.S."))))
micro.hulls.16Sg <- ddply(Entry_ords_16S, c("Entry", "Site"), find_hull)
pdf("InitialFigs/Beta_Bray_Genotype_16S.pdf", width = 7, height = 7)
g16 <- ggplot(Entry_ords_16S, aes(Axis01, Axis02)) +
  geom_polygon(data = micro.hulls.16Sg, 
               aes(colour = Site),
               alpha = 0, show.legend = F) +
  geom_point(size = 2, pch = 16, alpha = 0.75, aes(colour = Site)) +
  geom_text(data = subset(res_df_16S, p > 0.05),
            aes(x = 0, y = 0.3, label = Symbol),
            size = 3, inherit.aes = F) +
  geom_text(data = subset(res_df_16S, p <= 0.05),
            aes(x = 0, y = 0.34, label = Symbol),
            size = 8, inherit.aes = F) +
  labs(x = "PC1", 
       y = "PC2",
       colour = "Site") +
  scale_colour_viridis_d() +
  scale_y_continuous(expand = c(0.1, 0.1),
                     breaks = c(-0.2, 0, 0.2)) +
  scale_x_continuous(breaks = c(-0.2, 0, 0.2),
                     labels = c("-0.2", "0", "0.2")) +
  facet_wrap(~ Entry, ncol = 2) +
  guides(colour = guide_legend(override.aes = list(alpha = 1, size = 3))) +
  theme_bw() +  
  theme(legend.position = "right",
        axis.title = element_text(face = "bold", size = 12), 
        axis.text = element_text(size = 8),
        plot.title = element_text(vjust = -1),
        strip.text = element_text(size = 8),
        strip.background = element_rect(fill = "white"))
g16
dev.off()
leg_site <- get_legend(g16)
g16 <- g16 + theme(legend.position = "none")



# ITS
l_ITS <- list()
for (i in 1:length(unique(input_filt_rare_ITS$map_loaded$Entry))) {
  l_ITS[[i]] <- filter_data(input_filt_rare_ITS, 
                            filter_cat = "Entry", 
                            keep_vals = unique(input_filt_rare_ITS$map_loaded$Entry)[i])
}

l_bc_ITS <- list()
for (i in 1:length(l_ITS)) {
  l_bc_ITS[[i]] <- calc_dm(l_ITS[[i]]$data_loaded)
}

res_df_ITS <- as.data.frame(matrix(NA, nrow = length(l_ITS), ncol = 3)) %>%
  set_names(c("Site", "R2", "p"))
for (i in 1:length(l_ITS)) {
  set.seed(1000)
  pm <- adonis2(l_bc_ITS[[i]] ~ l_ITS[[i]]$map_loaded$Site)
  res_df_ITS$Site[i] <- unique(input_filt_rare_ITS$map_loaded$Site)[i]
  res_df_ITS$R2[i] <- pm$R2[1]
  res_df_ITS$p[i] <- pm$`Pr(>F)`[1]
}

for (i in 1:length(l_ITS)) {
  pcoa <- cmdscale(l_bc_ITS[[i]], k = nrow(l_ITS[[i]]$map_loaded) - 1, eig = T)
  pcoaA1 <- paste("PC1: ", round((eigenvals(pcoa)/sum(eigenvals(pcoa)))[1]*100, 1), "%")
  pcoaA2 <- paste("PC2: ", round((eigenvals(pcoa)/sum(eigenvals(pcoa)))[2]*100, 1), "%")
  l_ITS[[i]]$map_loaded$Axis01 <- vegan::scores(pcoa)[,1]
  l_ITS[[i]]$map_loaded$Axis02 <- vegan::scores(pcoa)[,2]
}

Entry_ords_ITS <- l_ITS[[1]]$map_loaded
for (i in 2:length(l_ITS)) {
  Entry_ords_ITS <- rbind(Entry_ords_ITS, l_ITS[[i]]$map_loaded)
}

res_df_ITS <- res_df_ITS %>%
  mutate(Symbol = ifelse(p <= 0.001,
                         "***",
                         ifelse(p > 0.001 & p <= 0.01,
                                "**",
                                ifelse(p > 0.01 & p <= 0.05,
                                       "*",
                                       "N.S."))))
micro.hulls.ITSg <- ddply(Entry_ords_ITS, c("Entry", "Site"), find_hull)
pdf("InitialFigs/Beta_Bray_Genotype_ITS.pdf", width = 7, height = 7)
g17 <- ggplot(Entry_ords_ITS, aes(Axis01, Axis02)) +
  geom_polygon(data = micro.hulls.ITSg, 
               aes(colour = Site),
               alpha = 0, show.legend = F) +
  geom_point(size = 2, pch = 16, alpha = 0.75, aes(colour = Site)) +
  geom_text(data = subset(res_df_ITS, p > 0.05),
            aes(x = 0, y = 0.3, label = Symbol),
            size = 3, inherit.aes = F) +
  geom_text(data = subset(res_df_ITS, p <= 0.05),
            aes(x = 0, y = 0.34, label = Symbol),
            size = 8, inherit.aes = F) +
  labs(x = "PC1", 
       y = "PC2",
       colour = "Genotype") +
  scale_colour_viridis_d() +
  scale_y_continuous(expand = c(0.1, 0.1),
                     breaks = c(-0.4, -0.2, 0, 0.2, 0.4)) +
  scale_x_continuous(breaks = c(-0.2, 0, 0.2),
                     labels = c("-0.2", "0", "0.2")) +
  facet_wrap(~ Entry, ncol = 2) +
  guides(colour = guide_legend(override.aes = list(alpha = 1, size = 3))) +
  theme_bw() +  
  theme(legend.position = "none",
        axis.title = element_text(face = "bold", size = 12), 
        axis.text = element_text(size = 8),
        plot.title = element_text(vjust = -1),
        strip.text = element_text(size = 8),
        strip.background = element_rect(fill = "white"))
g17
dev.off()


# Combined
pdf("FinalFigs/FigureS2.pdf", width = 9, height = 11)
plot_grid(g14, g15, leg_gen,
          g16, g17, leg_site,
          ncol = 3, rel_widths = c(0.42, 0.42, 0.16), rel_heights = c(0.55, 0.45))
dev.off()



#### 10. NCBI ####
#### _16S ####
repset_16S <- readFasta("data/repset_16S_filt.fasta") %>%
  separate(Header, into = c("Header", "Sequence_ID"), sep = " ")

# Need to get info for BioSamples
bs_16S <- input_filt_16S$map_loaded %>%
  arrange(Site) %>%
  mutate(geo_loc_name = paste("USA: ", State.x, sep = "")) %>%
  mutate(LONG_pos = abs(LONG)) %>%
  mutate(lat_lon = paste(LAT, "N", LONG_pos, "W", sep = " ")) %>%
  dplyr::select(sampleID, geo_loc_name, lat_lon)
write.csv(bs_16S, "data/mimarks_metadata_16S.csv", row.names = F)

# Now assign ASVs to samples (not comprehensive, just first sample)
info <- input_filt_16S$data_loaded
names(info) <- input_filt_16S$map_loaded$sampleID
for (i in 1:ncol(info)) {
  for (j in 1:nrow(info)) {
    ifelse(info[j, i] > 0, info[j, i] <- names(info)[i], info[j, i] <- "")
  }
}
info_cat <- info
info_cat <- info_cat %>%
  mutate_all(na_if, "") %>%
  mutate(unite(., "sample_name", c(names(info)), sep = ", ")) %>%
  mutate(sample_name = gsub("NA, ", "", sample_name)) %>%
  mutate(sample_name = gsub(", NA", "", sample_name)) %>%
  rownames_to_column(var = "Sequence_ID") %>%
  dplyr::select(Sequence_ID, sample_name)
info_first <- info_cat %>%
  separate(sample_name, into = c("sample_name", "Junk"), sep = ", ") %>%
  dplyr::select(Sequence_ID, sample_name)
info_first <- info_first %>%
  filter(Sequence_ID %in% repset_16S$Header)
#write_tsv(info_first, file = "~/Desktop/Sunflower/NCBI/biosample_assignment_16S.tsv")

# Filter NCBI flagged
flagged <- read.csv("data/NCBI_filter.csv") %>%
  separate(ASV, into = c("Junk", "ASV_ID"), sep = "\\|")
f <- readFasta("data/repset_16S_filt.fasta") %>%
  separate(Header, into = c("Header", "Sequence_ID"), sep = " ") %>%
  filter(Header %notin% flagged$ASV_ID)
microseq::writeFasta(f, "~/Desktop/Sunflower/NCBI/repset_used_filtered_16S.fasta")
sum(f$Header %in% flagged$ASV_ID)
info_first <- read_tsv("~/Desktop/Sunflower/NCBI/biosample_assignment_16S.tsv") %>%
  filter(Sequence_ID %notin% flagged$ASV_ID)
#write_tsv(info_first, file = "~/Desktop/Sunflower/NCBI/biosample_assignment_16S_filt.tsv")


#### _ITS ####
# Need to make "Source Modifiers" table
# Template was not downloading
# Go based off the example
repset_ITS <- readFasta("data/repset_ITS_filt.fasta") %>%
  separate(Header, into = c("Header", "Sequence_ID"), sep = " ")
sm <- data.frame("Sequence_ID" = repset_ITS$Header,
                 "Organism" = "uncultured fungus",
                 "Isolation-Source" = "sunflower rhizosphere soil",
                 "Country" = "USA",
                 "Collection-Date" = 2020,
                 "Strain" = repset_ITS$Header)
#write_tsv(sm, file = "~/Desktop/Sunflower/NCBI/source_modifiers.tsv")

# Filter NCBI flagged (35)
flagged <- read.csv("data/NCBI_filter_ITS.csv") %>%
  separate(ASV, into = c("Junk", "ASV_ID"), sep = "\\|")
f <- readFasta("data/repset_ITS_filt.fasta") %>%
  separate(Header, into = c("Header", "Sequence_ID"), sep = " ") %>%
  filter(Header %notin% flagged$ASV_ID)
microseq::writeFasta(f, "~/Desktop/Sunflower/NCBI/repset_used_filtered_ITS.fasta")
sum(f$Header %in% flagged$ASV_ID)
sm_filt <- sm %>%
  filter(Sequence_ID %in% f$Header)
#write_tsv(sm_filt, file = "~/Desktop/Sunflower/NCBI/source_modifiers_filt.tsv")

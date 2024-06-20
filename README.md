# Sunflower Rhizosphere Genotype x Environment

## Project Description (ABSTRACT)

>The rhizosphere microbial community may be important for mediating critical phenotypic attributes of cultivated sunflowers, namely resistance to fungal pathogen Sclerotinia sclerotiorum, which causes Sclerotinia basal stalk rot. The major goal of this study is to quantify the environmental variation in sunflower - rhizosphere associations. This study builds on previous work that identified microbial "taxa of interest" that potentially influence sclerotinia resistance. These analyses characterize plant-rhizosphere associations across the sunflower growing region to better understand how broadly relevant those initial findings were. These data will provide valuable information on the ubiquity of the taxa-of-interest, and how plant genotype-microbe associations vary among environments. 

### Experimental Questions
#### *Which ASVs associate strongly with sunflower genotype(s) (and/or Sclerotinia resistance)?* 
#### *Are there ASVs whose abundances are more strongly determined by genotype than by location?*
#### *How do plant genotype-microbe associations vary among environments?*

### Methods
- Common garden field experiment: 10 genotypes of sunflower each planted in 15 different locations spanning latitude of continental USA
- Harvested rhizosphere (soil around root) of sunflower plants, extracted DNA, characterized bacterial (16S) and fungal (ITS) communities via marker gene sequencing
- Sequence data processed using DADA2

## List of files

### R scripts

#### _ASV / OTU processing_
**ASV_OTU_scripts/dada2_sungxe_16S.Rmd** - bioinformatic processing of raw reads into 16S ASVs (Based on Fierer Lab DADA2 pipeline)  

**ASV_OTU_scripts/dada2_sungxe_ITS.Rmd** - bioinformatic processing of raw reads into ITS ASVs (Based on Fierer Lab DADA2 pipeline)  

**ASV_OTU_scripts/cluster97.Rmd** - bioinformatic processing of DADA2 16S ASVs into 97% OTU clusters using DECIPHER  

NOTE: these scripts used on Innes server

#### _16S Analyses_
**01_preprocess.Rmd** - initial processing of raw ASV table output from dada2 into cleaned tables for downstream analyses  

**02_explr.vis.Rmd** - exploration and visualization of bacterial community patterns across dataset  

**03_models.tmp.Rmd** - initial modeling attempts, incuding different approached for normalization and different models (e.g. glm, liner mixed effects, binomial)  

**04_RF.tmp.Rmd** - attempts at constructing Random Forest models to identify ASVs predictive of genotype. None of these models looked very good  

**05_ClimateData.Rmd** - incorporation of climate data for the locations, including building maps and adding new data to metadata file (script adapted from Kyle Keepers)  

**06_network.Rmd** - constructs a binomial network to visualize ASVs that associate with certain plant genotypes (need to identify better ASVs of interest before building this network)  

**07_MoreModels.Rmd** - further linear modeling attempts, including summarizing taxonomy at genera level  

**08_Cluster97.Rmd** - analysis and linear modeling of OTU 97% cluster  

#### _ITS Analyses_
**01_preprocess_ITS copy.Rmd** - initial processing of raw ASV table output from dada2 into cleaned tables for downstream analyses (for R1/R2 processed data)   

**01_preprocess_ITS_R1 copy.Rmd** - initial processing of raw ASV table output from dada2 into cleaned tables for downstream analyses (for R1 only data)   

**02_exploration_ITS copy.Rmd** - exploration and visualization of fungal community patterns across dataset


>[!NOTE]
>Large data files including sequence data will be deposited in (...link TBD, currently on Innes server)
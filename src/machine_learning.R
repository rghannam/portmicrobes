# =========================================================================================
# Script contents:
#
# Preprocess data, select amplicon sequence variant (ASV) features that will be used 
# throughout all downstream analysis
#
# Create machine learning models at all levels of taxonomic resolution
# =========================================================================================

rm(list=ls())
seed=81
set.seed(seed)

# Call count table, taxonomy table and metadata from GitHub
pm.counts<-'https://github.com/rghannam/pm_workflow/raw/master/data/counts_taxonomy/seqtabmergedNoC.rds'
pm.taxa<-'https://github.com/rghannam/pm_workflow/raw/master/data/counts_taxonomy/taxa.rds'
pm.metadata<-'https://github.com/rghannam/pm_workflow/raw/master/data/metadata/pm_metadata.csv'
pm.counts<-readRDS(url(pm.counts))
pm.taxa<-readRDS(url(pm.taxa))
pm.metadata<-read.csv(url(pm.metadata))
dim(pm.counts)
dim(pm.taxa)
dim(pm.metadata)

# Remove any arbitrary characters from sampleID's and match the metadata samples to the sequencing count table
pm.metadata$sample_name<-gsub("-", "", pm.metadata$sample_name)
row.names(pm.metadata)<-pm.metadata$sample_name
row.names(pm.counts)<-gsub("-", "", row.names(pm.counts))
rownames(pm.counts)
rownames(pm.metadata)

rowskeep<-rownames(pm.counts)
#pm.metadata<-pm.metadata[!pm.metadata$sample_name %in% rowskeep, ] # gives all missing, return 0, notrun
pm.metadata<-pm.metadata[pm.metadata$sample_name %in% rowskeep, ] # gives all missing

# Ensure observations are identical before workflow begins
if(!identical(rownames(pm.metadata), rownames(pm.counts))) stop("observations not identical")

# Manipulate counts,taxa,metadata that will make downstream analysis more efficient
colnames(pm.counts)<-NULL
dim(pm.counts)
options(max.print=2000)
pm.counts[,1514]

pm.taxa_t<-t(pm.taxa)
dim(pm.taxa_t)
colnames(pm.taxa_t) <- NULL
pm.taxa_re_t<-t(pm.taxa_t)
dim(pm.taxa_re_t)
dim(pm.metadata)
pm.metadata[, 1] <- as.factor(pm.metadata[, 1])

# Create phyloseq object
require(phyloseq)
sampledata = sample_data(pm.metadata)
ps.counts = otu_table(pm.counts, taxa_are_rows = FALSE) # leave false because we arbitrarily set column names with colnames()->NULL, we can always map back to respective ASV
ps.taxa = tax_table(pm.taxa_re_t)
ps.sampledata = sample_data(pm.metadata)
pm.ps = phyloseq(ps.counts, ps.taxa, ps.sampledata)
pm.ps
otu_table(pm.ps)[1:5, 1:5]

# Subset phyloseq object by open water, those that have reads >= 1000 and are from summer
pm.ps2<-subset_samples(pm.ps, sample_type == "open water" & reads_passed_filter_no_chimeras_overthousand == "yes" & season == 'Summer' )
pm.ps2

# Remove singletons and 0 sum features that would arise after subsetting
pm.ps2<-filter_taxa(pm.ps2, function(x) sum(x) > 1, TRUE)
pm.ps2

require(microbiome)
summarize_phyloseq(pm.ps2)
options(max.print=9999)
otu_table(pm.ps2)[1:400, 1:20]

# Filter out any ASVs that are non bacterial
pm.ps3 <- pm.ps2 %>%
  subset_taxa(
    Kingdom == "Bacteria" &
      Family  != "mitochondria" &
      Class   != "Chloroplast"
  )

pm.ps3

# Identify the total prevalence and abundance across all samples
# e.g. for sp1: 'Prevalence' means that this ASV has a value in 476 samples 
# TotalAbundance: summed total abundance of 376179 across all of the samples (476) that sp1 is in
prev.abund = apply(X = otu_table(pm.ps3),
                   MARGIN = ifelse(taxa_are_rows(pm.ps3), yes = 1, no = 2),
                   FUN = function(x){sum(x > 0)})
prev.abund = data.frame(Prevalence = prev.abund,
                        TotalAbundance = taxa_sums(pm.ps3),
                        tax_table(pm.ps3))
prev.abund
levels(prev.abund$Kingdom)
colSums(prev.abund[1])

# Prevalence filter by rows to keep only those ASVs within 15 or more samples
# Logical index of ASVs that are within >=15 samples
prev.abund2<-prev.abund[(prev.abund$Prevalence>=15) ,]


# Identify how many times an ASV that maps back to a phylum appears through all samples
# e.g. it's 18612 times that an ASV maps back to Actinobacteria were found in all samples (1218)
pm.totalavgprev.phylum <- plyr::ddply(prev.abund2, "Phylum", function(df1){cbind(mean(df1$Prevalence),sum(df1$Prevalence))})
as.data.frame(pm.totalavgprev.phylum) 
names(pm.totalavgprev.phylum) <- c("Phylum", "mean", "total")
colSums(pm.totalavgprev.phylum[3])

# Prevalence filter - remove everything that is not in 15 or more samples
prevdf1 = subset(prev.abund2, Phylum %in% get_taxa_unique(pm.ps3, "Phylum"))
prevdf1

# Execute prevalence filter, using `prune_taxa()` function
keepTaxa = rownames(prevdf1)
pm.ps4 = prune_taxa(keepTaxa, pm.ps3)
pm.ps4


# ==== Create input matrices and class labels ====
rm(list=ls()[! ls() %in% c("pm.ps4", "seed")])
set.seed(seed)
pkgs <- c("caretEnsemble", "caret", "randomForest", "plyr", "tidyverse", "data.table", "csv", "magrittr") # Load packages
invisible(sapply(pkgs,require, character = TRUE))


# Convert phyloseq counts to df, remove highly correlated variables, create local/regional class labels index
pstodf<-function(physeq, ...) {
  physeq<-filter_taxa(physeq, function(x) sum(x) > 1, TRUE) # rm 0sumfeat/singletons
  physeq <- microbiome::transform(physeq, transform = "log10", target = "OTU", shift=1) # log10(x+1 tform)
  #return(physeq) # we want it to a df, don't return physeq 
  #export as a df
  physeq.otu = as(otu_table(physeq), "matrix")
  if(taxa_are_rows(physeq)){physeq.otu <- t(physeq.otu)}
  physeq.otu.df = as.data.frame(physeq.otu)
  head(physeq.otu.df)[1:5]
  return(physeq.otu.df)
}

physeq.otu.df<-pstodf(pm.ps4)

dat.port<-physeq.otu.df


# Create local index for class labels
localind<-as.data.frame(rownames(dat.port))
names(localind)[1]<-"local"

localind$local<-gsub('[0-9]+', '', localind$local) # Remove any numerical after character string
library(stringr)
#  Remove any unwanted patterns in the strings from 'sampleid' factor
localind<-unlist(localind) # Convert to vector for stringr input
str(localind)
localind<-str_remove(localind, "post")
localind<-str_remove(localind, "point")
localind<-str_extract(localind, "^.{1,3}")
localind<-as.data.frame(localind)
names(localind)[1]<-"local"
dat.port<-cbind(dat.port, localind)
percentage <- prop.table(table(dat.port$local)) * 100
cbind(freq=table(dat.port$local), percentage=percentage)

# Randomize matrix
dat.port <- dat.port[sample(nrow(dat.port)),]
row.names(dat.port)[1:10]
dat.port$local[1:10]

(dat.port)[3215] # Sanity check for labels in proper order


# Create regional index for class labels
localind<-as.data.frame(rownames(physeq.otu.df))
names(localind)[1]<-"regionind"

localind$regionind<-gsub('[0-9]+', '', localind$regionind)
library(stringr)

localind<-unlist(localind)
str(localind)
localind<-str_remove(localind, "post")
localind<-str_remove(localind, "point")
localind<-str_extract(localind, "^.{1,3}")
regionind<-as.data.frame(localind)
names(regionind)[1]<-"regionind"

regionind$regionind <- str_replace_all(regionind$regionind, # Column we want to search
                                       c("BUS" = "ASIA","SIN" = "ASIA","HK" = "ASIA", # ASIA
                                         "EUM" = "EUR","EUN" = "EUR","EUR" = "EUR","EUV" = "EUR", "EUW" = "EUR", # EUROPE
                                         "BAL" = "EAST","CHA" = "EAST","GAL" = "EAST","NEW" = "EAST","NOF" = "EAST","NY" = "EAST",  # EAST
                                         "LA" = "WEST","OAK" = "WEST","SEA" = "WEST",
                                         "DS" = "LAKES", "GB" = "LAKES", "KEW" = "LAKES") # WEST
)

regionind$regionind<-as.factor(regionind$regionind)
levels(regionind$regionind)
names(regionind)[1]<-"region"


dat.region<-cbind(physeq.otu.df,regionind)
dat.region[3215] # Sanity

# Randomize matrix
dat.region <- dat.region[sample(nrow(dat.region)),]
row.names(dat.region)[1:10]
dat.region$region[1:10]

# ==== Run models ====
rm(list=ls()[! ls() %in% c("seed", "dat.port", "pm.ps4", "physeq.otu.df",
                           "trctrlbase", "modtypesbase", "metric")])

library("rlang")
library("Rcpp")
library("csv")
library("plyr")
library("dplyr")
require('caretEnsemble')
require('randomForest')
require('caret')
library('ModelMetrics')

set.seed(seed)

# Define train control, metric and base model ensembles
trctrlbase<- trainControl(method="repeatedcv", 
                          number=10, 
                          repeats=3, 
                          search="random", # No search parameters selected for global train control
                          summaryFunction = multiClassSummary, # Multiclass summary function, we're using 20 labels for port and 5 for region
                          classProbs=TRUE,
                          #savePredictions = "final",
                          preProc=c("center","scale")) # Remove mean value of each feature and divide non constant features by standard deviation
metric <- "Accuracy"

modtypesbase <- list(rf = caretModelSpec(method="rf", verbose=FALSE))


# ==== ASV local model =====
asv.local <-caretList(local~., 
                      data=dat.port, ntree=501,
                      trControl = trctrlbase, 
                      tuneList = modtypesbase)

asv.local
asv.local[["rf"]][["finalModel"]]


# ==== Regional model =====
asv.region <-caretList(region~., 
                       data=dat.region, ntree=501,
                       trControl = trctrlbase, 
                       tuneList = modtypesbase)

asv.region
asv.region[["rf"]][["finalModel"]]


# ==== Taxonomic models Phylum - Genus ====

# Extract taxonomic level dataframes from ps object
require(metagMisc)
require(microbiome)

pkgs <- c("hpgltools", "DESeq2", "plyr", "dplyr", "dbplyr","stringr", "phyloseq",
          "data.table", "tidyverse", "tibble", "caretEnsemble", "caret", "randomForest", "plyr", 
          "tidyverse", "data.table", "csv", "magrittr") 
invisible(sapply(pkgs,require, character = TRUE))

rm(list=ls()[! ls() %in% c("pm.ps4", "seed", "trctrlbase", "modtypesbase", "metric")])
seed=81
set.seed(seed)


# Agglomerate all taxonomic levels
pm.ps.phylum = tax_glom(pm.ps4, "Phylum", NArm = TRUE, bad_empty=c(NA, "", " ", "\t"))
pm.ps.phylum
pm.ps.class = tax_glom(pm.ps4, "Class", NArm = TRUE, bad_empty=c(NA, "", " ", "\t"))
pm.ps.class
pm.ps.order = tax_glom(pm.ps4, "Order", NArm = TRUE, bad_empty=c(NA, "", " ", "\t"))
pm.ps.order
pm.ps.family = tax_glom(pm.ps4, "Family", NArm = TRUE, bad_empty=c(NA, "", " ", "\t"))
pm.ps.family
pm.ps.genus = tax_glom(pm.ps4, "Genus", NArm = TRUE, bad_empty=c(NA, "", " ", "\t"))
pm.ps.genus

# Logarithmize counts: log10(x+1)
pm.ps.phylum <- microbiome::transform(pm.ps.phylum, transform = "log10", target = "OTU", shift=1) # log10(x+1 tform)
pm.ps.class <- microbiome::transform(pm.ps.class, transform = "log10", target = "OTU", shift=1) # log10(x+1 tform)
pm.ps.order <- microbiome::transform(pm.ps.order, transform = "log10", target = "OTU", shift=1) # log10(x+1 tform)
pm.ps.family <- microbiome::transform(pm.ps.family, transform = "log10", target = "OTU", shift=1) # log10(x+1 tform)
pm.ps.genus <- microbiome::transform(pm.ps.genus, transform = "log10", target = "OTU", shift=1) # log10(x+1 tform)


# Create a count matrix of Class feature for each taxonomic level
# This way, can get the actual feature (bacterial class name etc) rather than sp3
#Phylum
phy.melt <- psmelt(pm.ps.phylum)
phy.melt$Phylum <- as.character(phy.melt$Phylum)
phy.melt <- aggregate(Abundance~Sample+Phylum, phy.melt, FUN=sum)
phy.mat <- cast(phy.melt, Sample ~ Phylum)
rownames(phy.mat)<-phy.mat[,1]
phy.mat[1]<-NULL

#Class
class.melt <- psmelt(pm.ps.class)
class.melt$Class <- as.character(class.melt$Class)
class.melt <- aggregate(Abundance~Sample+Class, class.melt, FUN=sum)
class.mat <- cast(class.melt, Sample ~ Class)
rownames(class.mat)<-class.mat[,1]
class.mat[1]<-NULL

#Order
order.melt <- psmelt(pm.ps.order)
order.melt$Order <- as.character(order.melt$Order)
order.melt <- aggregate(Abundance~Sample+Order, order.melt, FUN=sum)
order.mat <- cast(order.melt, Sample ~ Order)
rownames(order.mat)<-order.mat[,1]
order.mat[1]<-NULL

#Family
family.melt <- psmelt(pm.ps.family)
family.melt$Family <- as.character(family.melt$Family)
family.melt <- aggregate(Abundance~Sample+Family, family.melt, FUN=sum)
family.mat <- cast(family.melt, Sample ~ Family)
rownames(family.mat)<-family.mat[,1]
family.mat[1]<-NULL

#Genus
genus.melt <- psmelt(pm.ps.genus)
genus.melt$Genus <- as.character(genus.melt$Genus)
genus.melt <- aggregate(Abundance~Sample+Genus, genus.melt, FUN=sum)
genus.mat <- cast(genus.melt, Sample ~ Genus)
rownames(genus.mat)<-genus.mat[,1]
genus.mat[1]<-NULL

# Convert to matrix
# Store taxa count matrices in a list
taxalist<-list(phy.mat,class.mat,order.mat,family.mat,genus.mat)
listNames<-c("Phylum", "Class", "Order","Family", "Genus")
names(taxalist)<-listNames

# Create class labels for all taxonomic level matrices
# Can't use same vectors we used for creating class labels for asv models, casting with reshape puts the samples
# in alphabetical order.

# Local labels
localind<-as.data.frame(rownames(taxalist[[1]]))
names(localind)[1]<-"local"

localind$local<-gsub('[0-9]+', '', localind$local)
library(stringr)
localind<-unlist(localind)
str(localind)
localind<-str_remove(localind, "post")
localind<-str_remove(localind, "point")
localind<-str_extract(localind, "^.{1,3}")
localind<-as.data.frame(localind)
names(localind)[1]<-"local"

# Append the local labels to each element of the list
taxalist.loc<-lapply(taxalist, function(x) 
  cbind(x, localind = localind$local))
options(max.print=2000)

# Sanity check for labels in proper order
rownames(taxalist.loc[["Order"]])
taxalist.loc[["Order"]][["localind"]]
phylabs<-taxalist.loc[["Phylum"]][["localind"]]
orderlabs<-taxalist.loc[["Order"]][["localind"]]
if(!identical(phylabs, orderlabs)) stop("not identical")

# Check class balance
percentage <- prop.table(table(taxalist.loc[[1]]$localind)) * 100
cbind(freq=table(taxalist.loc[[1]]$localind), percentage=percentage)



# Regional labels
localind<-as.data.frame(rownames(taxalist[[1]]))
names(localind)[1]<-"regionind"

localind$regionind<-gsub('[0-9]+', '', localind$regionind)
library(stringr)
localind<-unlist(localind)
str(localind)
localind<-str_remove(localind, "post")
localind<-str_remove(localind, "point")
localind<-str_extract(localind, "^.{1,3}")
regionind<-as.data.frame(localind)
names(regionind)[1]<-"regionind"

regionind$regionind <- str_replace_all(regionind$regionind, # Column we want to search
                                       c("BUS" = "ASIA","SIN" = "ASIA","HK" = "ASIA", #ASIA
                                         "EUM" = "EUR","EUN" = "EUR","EUR" = "EUR","EUV" = "EUR", "EUW" = "EUR", #EUROPE
                                         "BAL" = "EAST","CHA" = "EAST","GAL" = "EAST","NEW" = "EAST","NOF" = "EAST","NY" = "EAST",  #EAST
                                         "LA" = "WEST","OAK" = "WEST","SEA" = "WEST",
                                         "DS" = "LAKES", "GB" = "LAKES", "KEW" = "LAKES") #WEST
)

regionind$regionind<-as.factor(regionind$regionind)
levels(regionind$regionind)
names(regionind)[1]<-"region"

# Append the regional labels to each element of the list
taxalist.reg<-lapply(taxalist, function(x) 
  cbind(x, regionind = regionind$region))

# Sanity checks for order of proper class labels
rownames(taxalist.reg[["Order"]])
taxalist.reg[["Order"]][["regionind"]]
phylabs<-taxalist.reg[["Phylum"]][["regionind"]]
orderlabs<-taxalist.reg[["Order"]][["regionind"]]
if(!identical(phylabs, orderlabs)) stop("not identical")

# Randomize the observations (samples) across all levels of taxonomic resolution for region/ports
for (i in 1:length(taxalist.loc)){
  taxalist.loc[[i]] <- taxalist.loc[[i]][sample(nrow(taxalist.loc[[i]])),]
}

for (i in 1:length(taxalist.reg)){
  taxalist.reg[[i]] <- taxalist.reg[[i]][sample(nrow(taxalist.reg[[i]])),]
}


local.models<-taxalist.loc
region.models<-taxalist.reg

# Loop the caretlist models across all levels of taxonomic resolution

#Local
for(i in 1:length(local.models)){
  local.models[[i]]<-caretList(localind~.,
                               data=local.models[[i]], ntree=501,
                               trControl=trctrlbase, tuneList=modtypesbase)
}
local.models


#Region
for(i in 1:length(region.models)){
  region.models[[i]]<-caretList(regionind~.,
                                data=region.models[[i]], ntree=501,
                                trControl=trctrlbase, tuneList=modtypesbase)
}
region.models


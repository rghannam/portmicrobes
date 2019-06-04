# This will normalize count matrices at all levels of taxonomic resolution


detachAllPackages <- function() {
  
  basic.packages <- c("package:stats","package:graphics","package:grDevices","package:utils","package:datasets","package:methods","package:base")
  
  package.list <- search()[ifelse(unlist(gregexpr("package:",search()))==1,TRUE,FALSE)]
  
  package.list <- setdiff(package.list,basic.packages)
  
  if (length(package.list)>0)  for (package in package.list) detach(package, character.only=TRUE)
  
}

detachAllPackages()

library("rlang")
library("ggplot2")
library("Rcpp")
library("ggplot2")
library("vegan")
library("csv")
library("plotly")
library("plyr")
library("dplyr")
library("hpgltools")


#  Explore data through phyloseq and generate a taxonomy table across all levels of taxonomic resolution
seed=81
set.seed(seed)
setwd("~/pm_workflow") # Change to whatever directory

# Load data
pm.counts<-'https://github.com/rghannam/pm_workflow/raw/master/data/counts_taxonomy/pm_counts1218.rds'
pm.counts<-readRDS(url(pm.counts))
dim(pm.counts)

pm.taxa<-'https://github.com/rghannam/pm_workflow/raw/master/data/counts_taxonomy/pm_taxonomy.rds'
pm.taxa<-readRDS(url(pm.taxa))
dim(pm.taxa)

pm.metadata<-'https://github.com/rghannam/pm_workflow/raw/master/data/metadata/pm_metadata.csv'
pm.metadata<-read.csv(url(pm.metadata))
dim(pm.metadata)

# Data manipulation for phyloseq object format
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

# Samples must match from counts to metadata, lets format the metadata to match the count table
row.names(pm.metadata)<-pm.metadata$sample_name
dim(pm.counts)
samplenames<-rownames(pm.counts)
pm.metadata<-pm.metadata[match(samplenames, pm.metadata$sample_name),]

library(phyloseq)
# Create phyloseq object
sampledata = sample_data(pm.metadata)
counts = otu_table(pm.counts, taxa_are_rows = FALSE) # Leave false because we arbitrarily set column names with colnames()->NULL, we can always map back to respective ASV
taxa = tax_table(pm.taxa_re_t)
sampledata = sample_data(pm.metadata)
pm.ps = phyloseq(counts, taxa, sampledata)
pm.ps
otu_table(pm.ps)[1:5, 1:5]
#  Our working phyloseq object is now ready at a dimension of 119097 ASVs by 1218 samples
rm(list=ls()[! ls() %in% c("seed", "pm.ps")]) # Clear memory, our counts, taxa table & metadata are stored in phyloseq object



# Phyloseq analysis

# Explore taxa at various levels to see number of ASVs (features) per taxonomic rank 
options(max.print=2000) # n=1218
table(tax_table(pm.ps)[, "Phylum"], exclude = NULL) #3066
table(tax_table(pm.ps)[, "Class"], exclude = NULL) #5430
table(tax_table(pm.ps)[, "Order"], exclude = NULL) #14336
table(tax_table(pm.ps)[, "Family"], exclude = NULL) #31669
table(tax_table(pm.ps)[, "Genus"], exclude = NULL) #65072

# Transform to relative abundance. Save as new object.

# Exclude the uncharacterized or NA taxa
# Create new phyloseq object containing all ASVs that actually do characterize as a Phylum
p<-subset_taxa(pm.ps, !is.na(Phylum) & !Phylum %in% c("", "uncharacterized"))
c<-subset_taxa(pm.ps, !is.na(Class) & !Class %in% c("", "uncharacterized"))
o<-subset_taxa(pm.ps, !is.na(Order) & !Order %in% c("", "uncharacterized"))
f<-subset_taxa(pm.ps, !is.na(Family) & !Family %in% c("", "uncharacterized"))
g<-subset_taxa(pm.ps, !is.na(Genus) & !Genus %in% c("", "uncharacterized"))
x<-list()
alltaxa<-c(list(p,c,o,f,g), x)
rm(list=ls()[! ls() %in% c("alltaxa", "seed")]) # Clear memory
listNames <- c("phylum","class", "order", "family", "genus") # Renames to proper taxonomic level
names(alltaxa) <- listNames


# Compute prevalence of each feature, store as data.frame
#  This is a table containing information about prevalence/total abundance of all of our samples and ASVs
#  Prevalence = number samples that contain that ASV
#  e.g. sp4 maps to Bacteroidetes and sp4 was found in 695 samples
#  TotalAbundance = summed ASV counts across all samples mapping back to that ASV

#Phylum
prev.abund.phylum = apply(X = otu_table(alltaxa$phylum),
                          MARGIN = ifelse(taxa_are_rows(alltaxa$phylum), yes = 1, no = 2),
                          FUN = function(x){sum(x > 0)})
#  Add taxonomy and total read counts to this data.frame
prev.abund.phylum = data.frame(Prevalence = prev.abund.phylum,
                               TotalAbundance = taxa_sums(alltaxa$phylum),
                               tax_table(alltaxa$phylum))

#Class
prev.abund.class = apply(X = otu_table(alltaxa$class),
                         MARGIN = ifelse(taxa_are_rows(alltaxa$class), yes = 1, no = 2),
                         FUN = function(x){sum(x > 0)})
prev.abund.class = data.frame(Prevalence = prev.abund.class,
                              TotalAbundance = taxa_sums(alltaxa$class),
                              tax_table(alltaxa$class))

#Order
prev.abund.order = apply(X = otu_table(alltaxa$order),
                         MARGIN = ifelse(taxa_are_rows(alltaxa$order), yes = 1, no = 2),
                         FUN = function(x){sum(x > 0)})
prev.abund.order = data.frame(Prevalence = prev.abund.order,
                              TotalAbundance = taxa_sums(alltaxa$order),
                              tax_table(alltaxa$order))

#Family
prev.abund.family = apply(X = otu_table(alltaxa$family),
                          MARGIN = ifelse(taxa_are_rows(alltaxa$family), yes = 1, no = 2),
                          FUN = function(x){sum(x > 0)})
prev.abund.family = data.frame(Prevalence = prev.abund.family,
                               TotalAbundance = taxa_sums(alltaxa$family),
                               tax_table(alltaxa$family))

#Genus
prev.abund.genus = apply(X = otu_table(alltaxa$genus),
                         MARGIN = ifelse(taxa_are_rows(alltaxa$genus), yes = 1, no = 2),
                         FUN = function(x){sum(x > 0)})
prev.abund.genus = data.frame(Prevalence = prev.abund.genus,
                              TotalAbundance = taxa_sums(alltaxa$genus),
                              tax_table(alltaxa$genus))


# Which taxa/features consist of low prevalence ASVs?
#  e.g. Bacteroidetes was in a total prevalence of 105958
#  This is the total sum of how many times an ASV that maps to Bacteroidetes is contained in any number of samples
#  colSums on all 'Phylum=Bacteroidetes' should = 105958

#Phylum
pm.totalavgprev.phylum <- plyr::ddply(prev.abund.phylum, "Phylum", function(df1){cbind(mean(df1$Prevalence),sum(df1$Prevalence))})
as.data.frame(pm.totalavgprev.phylum)
names(pm.totalavgprev.phylum) <- c("Phylum", "mean", "total")

#Class
pm.totalavgprev.class <- plyr::ddply(prev.abund.class, "Class", function(df1){cbind(mean(df1$Prevalence),sum(df1$Prevalence))})
as.data.frame(pm.totalavgprev.class)
names(pm.totalavgprev.class) <- c("Class", "mean", "total")

#Order
pm.totalavgprev.order <- plyr::ddply(prev.abund.order, "Order", function(df1){cbind(mean(df1$Prevalence),sum(df1$Prevalence))})
as.data.frame(pm.totalavgprev.order)
names(pm.totalavgprev.order) <- c("Order", "mean", "total")

#Family
pm.totalavgprev.family <- plyr::ddply(prev.abund.family, "Family", function(df1){cbind(mean(df1$Prevalence),sum(df1$Prevalence))})
as.data.frame(pm.totalavgprev.family)
names(pm.totalavgprev.family) <- c("Family", "mean", "total")

#Genus
pm.totalavgprev.genus <- plyr::ddply(prev.abund.genus, "Genus", function(df1){cbind(mean(df1$Prevalence),sum(df1$Prevalence))})
as.data.frame(pm.totalavgprev.genus)
names(pm.totalavgprev.genus) <- c("Genus", "mean", "total")



# ---- Skip unsupervised prevalence filtering ----

# Here, we skip unsupervised prevalence filtering.
# We can define taxa to remove via 'filterPhyla' but lets choose a threshold instead
# Lets find which taxa are prevalent and in high abundance and use these for setting up filtering thresholds

#pm.gg.df.Phylum = subset(pm.prev.Phylum, Phylum %in% get_taxa_unique(pm.psPhylum, "Phylum"))
#pm.plot.Phylum<-ggplot(pm.gg.df.Phylum, aes(TotalAbundance, Prevalence / nsamples(pm.psPhylum),color=Phylum)) +
#  geom_hline(yintercept = 0.05, alpha = 0.5, linetype = 2) + geom_point(size = 0.3, alpha = 0.7) +
#  scale_x_log10() +  xlab("Total Abundance") + ylab("Prevalence [Frac. Samples]") +
#  facet_wrap(~Phylum) + theme(legend.position="none")
#pm.plot.Phylum # Each point is a different ASV

#prevThresh = 0.05 * nsamples(pm.psPhylum)
#prevThresh

# Execute prevalence filter, using `prune_taxa()` function
#keepTaxa = rownames(pm.gg.df.Phylum)[(pm.gg.df.Phylum$Prevalence >= prevThresh)] # Define which taxa to keep
#pm.psPhylum2 = prune_taxa(keepTaxa, pm.psPhylum) # Prune from phyloseq object
#pm.psPhylum2
#pm.prev.Phylum2 = apply(X = otu_table(pm.psPhylum2),
                        #MARGIN = ifelse(taxa_are_rows(pm.psPhylum2), yes = 1, no = 2),
                        #FUN = function(x){sum(x > 0)})
#  Add taxonomy and total read counts to this data.frame
#pm.prev.Phylum2 = data.frame(Prevalence = pm.prev.Phylum2,
                             #TotalAbundance = taxa_sums(pm.psPhylum2),
                             #tax_table(pm.psPhylum2))

#pm.totalavgprev.Phylum2 <- plyr::ddply(pm.prev.Phylum2, "Phylum", function(df1){cbind(mean(df1$Prevalence),sum(df1$Prevalence))})
#as.data.frame(pm.totalavgprev.Phylum2)
#pm.totalavgprev.Phylum2
#pm.totalavgprev.Phylum

# ---- End skip unsupervised prevalence filtering ----


#  Notice we're trimming # of phyla used in our phyloseq object and overall taxa (ASVs)
#  This in turn reduces abundance
#  What we're effetively doing is only keeping the ones that are within a certain % of samples.
#  Our dimensionality of data goes down too much just by filtering by 5% prev
#  This is not a good way of filtering since:
#    1) We want more than 11 phyla to make high throughput differential expression pairwise comparisons
#    2) We want machine learning models to have more features to work with if we use a low resolution matrix (e.g. phylum)
#  We don't want to flush out our signatures, prevalence probing the environment is unknown atm
#  If we can design a probe for low prevalence detection, it's fine to keep them in our dataset
#  Thinking about keeping this model robust enough to extend to field application

#  Define prevalence threshold as 1.23% of total samples (15/1218)
#  Effectively what we're doing is removing sparsity and low prevalent features for downstream analysis, machine learning etc

#  Setting a prevalence threshold of >=15 per taxonomic ranking irrespective of abundance
#  Reasoning, importance of abundance for our workflow is to be decided later and not to be governed by prevalence thresholding, although slightly correlated



# Prevalence filtering supervised
#  Phylum
pm.totalavgprev.phylum.subset <- subset(pm.totalavgprev.phylum, pm.totalavgprev.phylum$total <= 15) # Filter entries with <=15 prevalence
extract.phylum<-as.data.frame(pm.totalavgprev.phylum.subset[,1])
names(extract.phylum) <- ("Phylum") 
extract.phylum <- as.list(as.vector(t(extract.phylum))) # Create a list of the extracted so we can subset and subtract these

filterphyla = (extract.phylum)
pm.ps.phylum.sub = subset_taxa(alltaxa$phylum, !Phylum %in% filterphyla)
pm.ps.phylum.sub
alltaxa$phylum

length(get_taxa_unique(alltaxa$phylum, taxonomic.rank = "Phylum")) # 71
length(get_taxa_unique(pm.ps.phylum.sub, taxonomic.rank = "Phylum")) # 45

# Class
pm.totalavgprev.class.subset <- subset(pm.totalavgprev.class, pm.totalavgprev.class$total <= 15)
extract.class<-as.data.frame(pm.totalavgprev.class.subset[,1])
names(extract.class) <- ("Class") 
extract.class <- as.list(as.vector(t(extract.class)))

filterclass = (extract.class)
pm.ps.class.sub = subset_taxa(alltaxa$class, !Class %in% filterclass)
pm.ps.class.sub
alltaxa$class

length(get_taxa_unique(alltaxa$class, taxonomic.rank = "Class")) # 143
length(get_taxa_unique(pm.ps.class.sub, taxonomic.rank = "Class")) # 88

#Order
pm.totalavgprev.order.subset <- subset(pm.totalavgprev.order, pm.totalavgprev.order$total <= 15)
extract.order<-as.data.frame(pm.totalavgprev.order.subset[,1])
names(extract.order) <- ("Order") 
extract.order <- as.list(as.vector(t(extract.order)))

filterorder = (extract.order)
pm.ps.order.sub = subset_taxa(alltaxa$order, !Order %in% filterorder)
pm.ps.order.sub
alltaxa$order

length(get_taxa_unique(alltaxa$order, taxonomic.rank = "Order")) # 357
length(get_taxa_unique(pm.ps.order.sub, taxonomic.rank = "Order")) # 214

# Family
pm.totalavgprev.family.subset <- subset(pm.totalavgprev.family, pm.totalavgprev.family$total <= 15)
extract.family<-as.data.frame(pm.totalavgprev.family.subset[,1])
names(extract.family) <- ("Family") 
extract.family <- as.list(as.vector(t(extract.family)))

filterfamily = (extract.family)
pm.ps.family.sub = subset_taxa(alltaxa$family, !Family %in% filterfamily)
pm.ps.family.sub
alltaxa$family

length(get_taxa_unique(alltaxa$family, taxonomic.rank = "Family")) # 527
length(get_taxa_unique(pm.ps.family.sub, taxonomic.rank = "Family")) # 345

# Genus
pm.totalavgprev.genus.subset <- subset(pm.totalavgprev.genus, pm.totalavgprev.genus$total <= 15)
extract.genus<-as.data.frame(pm.totalavgprev.genus.subset[,1])
names(extract.genus) <- ("Genus") 
extract.genus <- as.list(as.vector(t(extract.genus)))

filtergenus = (extract.genus)
pm.ps.genus.sub = subset_taxa(alltaxa$genus, !Genus %in% filtergenus)
pm.ps.genus.sub
alltaxa$genus

length(get_taxa_unique(alltaxa$genus, taxonomic.rank = "Genus")) # 1945
length(get_taxa_unique(pm.ps.genus.sub, taxonomic.rank = "Genus")) # 809

# Agglomerate our prevalence filtered counts by taxa rankings
# e.g. combine all ASVs into respective taxonomy
rank_names(pm.ps.phylum.sub)
rm(list=ls()[! ls() %in% c("seed", "pm.ps.phylum.sub", "pm.ps.class.sub", "pm.ps.order.sub", "pm.ps.family.sub", "pm.ps.genus.sub", "alltaxa")]) # Clear memory

#Phylum
pm.ps.phylum = tax_glom(pm.ps.phylum.sub, "Phylum", NArm = TRUE, bad_empty=c(NA, "", " ", "\t"))
pm.ps.phylum # 45

#Class
pm.ps.class = tax_glom(pm.ps.class.sub, "Class", NArm = TRUE, bad_empty=c(NA, "", " ", "\t"))
pm.ps.class

#Order
pm.ps.order = tax_glom(pm.ps.order.sub, "Order", NArm = TRUE, bad_empty=c(NA, "", " ", "\t"))
pm.ps.order

#Family
pm.ps.family = tax_glom(pm.ps.family.sub, "Family", NArm = TRUE, bad_empty=c(NA, "", " ", "\t"))
pm.ps.family

#Genus
pm.ps.genus = tax_glom(pm.ps.genus.sub, "Genus", NArm = TRUE, bad_empty=c(NA, "", " ", "\t"))
pm.ps.genus


rm(list=ls()[! ls() %in% c("pm.ps.phylum", "pm.ps.class", "pm.ps.order", "pm.ps.family", "pm.ps.genus", "seed")]) # Clear memory
pkgs <- c("phyloseq", "reshape", "reshape2", "ggplot2")
invisible(sapply(pkgs,require, character = TRUE))

# Generate count matrices across all taxonomic ranks
#Phylum
pm.ps.glom.phylum <- psmelt(pm.ps.phylum)
pm.ps.glom.phylum$Phylum <- as.character(pm.ps.glom.phylum$Phylum)
pm.ps.glom.phylum
pm.ps.glom.phylum <- aggregate(Abundance~Sample+Phylum, pm.ps.glom.phylum, FUN=sum)
pm.phylumMatrix <- cast(pm.ps.glom.phylum, Sample ~ Phylum)
dim(pm.phylumMatrix)
head(pm.phylumMatrix)

#Class
pm.ps.glom.class <- psmelt(pm.ps.class)
pm.ps.glom.class$Class <- as.character(pm.ps.glom.class$Class)
pm.ps.glom.class
pm.ps.glom.class <- aggregate(Abundance~Sample+Class, pm.ps.glom.class, FUN=sum)
pm.classMatrix <- cast(pm.ps.glom.class, Sample ~ Class)
dim(pm.classMatrix)
head(pm.classMatrix)

#Order
pm.ps.glom.order <- psmelt(pm.ps.order)
pm.ps.glom.order$order <- as.character(pm.ps.glom.order$Order)
pm.ps.glom.order
pm.ps.glom.order <- aggregate(Abundance~Sample+order, pm.ps.glom.order, FUN=sum)
pm.orderMatrix <- cast(pm.ps.glom.order, Sample ~ order)
dim(pm.orderMatrix)
head(pm.orderMatrix)

#Family
pm.ps.glom.family <- psmelt(pm.ps.family)
pm.ps.glom.family$family <- as.character(pm.ps.glom.family$Family)
pm.ps.glom.family
pm.ps.glom.family <- aggregate(Abundance~Sample+family, pm.ps.glom.family, FUN=sum)
pm.familyMatrix <- cast(pm.ps.glom.family, Sample ~ family)
dim(pm.familyMatrix)
head(pm.familyMatrix)

#Genus
pm.ps.glom.genus <- psmelt(pm.ps.genus)
pm.ps.glom.genus$genus <- as.character(pm.ps.glom.genus$Genus)
pm.ps.glom.genus
pm.ps.glom.genus <- aggregate(Abundance~Sample+genus, pm.ps.glom.genus, FUN=sum)
pm.genusMatrix <- cast(pm.ps.glom.genus, Sample ~ genus)
dim(pm.genusMatrix)
head(pm.genusMatrix)

rm(list=ls()[! ls() %in% c("pm.phylumMatrix", "pm.classMatrix", "pm.orderMatrix", "pm.familyMatrix", "pm.genusMatrix", "seed")]) # Clear memory

pm.classMatrix<-as.data.frame(pm.classMatrix)
pm.familyMatrix<-as.data.frame(pm.familyMatrix)
pm.genusMatrix<-as.data.frame(pm.genusMatrix)
pm.orderMatrix<-as.data.frame(pm.orderMatrix)
pm.phylumMatrix<-as.data.frame(pm.phylumMatrix)

taxa_counts_list<-list(pm.classMatrix, pm.familyMatrix, pm.genusMatrix, pm.orderMatrix, pm.phylumMatrix)


# Normalize phyloseq count tables

path <- "~/pm_workflow" 
setwd("~/pm_workflow")

pkgs <- c("hpgltools", "DESeq2", "plyr", "dplyr", "dbplyr","stringr",
          "data.table", "tidyverse", "tibble") 
invisible(sapply(pkgs,require, character = TRUE))

# Import a list of all of the taxonomic count matrices
taxa_counts_list<-'https://github.com/rghannam/pm_workflow/raw/master/data/counts_taxonomy/taxa_counts_list.rds'
taxa_counts_list<-readRDS(url(taxa_counts_list))


listNames <- c("Class","Family","Genus", "Order", "Phylum") # Renames elements to proper taxonomic level
names(taxa_counts_list) <- listNames

for(i in 1:length(taxa_counts_list)){ # Drops first Sample column, since SampleID is stored as row.name, it's fine
  rownames(taxa_counts_list[[i]]) <- taxa_counts_list[[i]][,1]
  taxa_counts_list[[i]] <- taxa_counts_list[[i]][,-1,drop = FALSE]
}

# Transpose all for DA analysis later
tlist<-taxa_counts_list
tlist$Class<-as.data.frame(t(tlist$Class))
tlist$Family<-as.data.frame(t(tlist$Family))
tlist$Genus<-as.data.frame(t(tlist$Genus))
tlist$Order<-as.data.frame(t(tlist$Order))
tlist$Phylum<-as.data.frame(t(tlist$Phylum))

rm(list=ls()[! ls() %in% c("tlist", "seed")])

# Subset our metadata based on count table samples
pm.metadata<-'https://github.com/rghannam/pm_workflow/raw/master/data/metadata/pm_metadata.csv'
pm.metadata<-read.csv(url(pm.metadata))
dim(pm.metadata)

pm.counts<-'https://github.com/rghannam/pm_workflow/raw/master/data/counts_taxonomy/pm_counts1218.rds'
pm.counts<-readRDS(url(pm.counts))
dim(pm.counts)

rownames(pm.counts)
pm.metadata$sample_name
row.names(pm.counts)<-gsub("-", "", row.names(pm.counts)) # Remove '-" from sample names
pm.metadata$sample_name<-gsub("-", "", pm.metadata$sample_name)

# Arrange metadata
row.names(pm.metadata)<-pm.metadata$sample_name
dim(pm.counts)

# We need to match the samples from counts to our metadata
# Create vector of sample names from counts
samplenames<-rownames(pm.counts)
pm.metadata.arrange<-pm.metadata[match(samplenames, pm.metadata$sample_name),]

# Append metadata to subset our main count table per what samples we want (open water)
colnames(pm.metadata.arrange) # The sample variables we can subset our main count table by
names(pm.metadata.arrange)[names(pm.metadata.arrange) == 'sequence_run'] <- 'batch' # Rename 'sequence_run' to 'batches' for experimental class object (DA analysis)
pm.metadata.arrange <- pm.metadata.arrange[pm.metadata.arrange$sample_type == "open water" & 
                                             pm.metadata.arrange$reads_passed_filter_no_chimeras_overthousand == "yes" & 
                                             pm.metadata.arrange$season == 'Summer', ]
as.data.frame(rownames(pm.metadata.arrange))
length(unique(pm.metadata.arrange$sample_name))

allrownames<-as.data.frame(tlist[[1]])

tt<-t(tlist[[1]])
allrownames<-row.names(tt)
pm.metadata.arrange2=cbind(pm.metadata.arrange,allrownames)

library(dbplyr)
class(pm.metadata.arrange2)
pm.metadata.arrange2$sample_name<-as.factor(pm.metadata.arrange2$sample_name)
levels(pm.metadata.arrange2$sample_name)

metadata.subset.work<-pm.metadata.arrange2
library(dplyr)
metadata.subset.rearranged<-metadata.subset.work %>% arrange(sample_name, desc(allrownames))
metadata.subset.rearranged$allrownames<-NULL

#  Here's where we remove the taxa variables and rename them to X, we will apply their respective names back to them later on
dim(metadata.subset.rearranged) # Sanity check: 1218 observations
metadata.subset.rearranged$geo_port<-as.factor(metadata.subset.rearranged$geo_port)
metadata.subset.rearranged$batch<-as.factor(metadata.subset.rearranged$batch)
metadata.subset.rearranged$sample_name<-as.factor(metadata.subset.rearranged$sample_name)
names(metadata.subset.rearranged)[1]<-"sampleid" # had to rename to sampleid
names(metadata.subset.rearranged)[37]<-"condition" # had to rename to geo_port to condition

seed=81
set.seed(seed)

# Create experimental class object of all taxonomic levels (Phylum-Genus)
# We want to create an object where we can pass to a high throughput pairwise function for DESEq2 differential expression.

#Phylum
rm(list=ls()[! ls() %in% c("metadata.subset.rearranged","tlist", "seed")]) # Clear memory
save(list = ls(), file = "~/pm_workflow/normalize/pm_phylum_pairwise_env.RData") # Save these 2 objects as RData
pm.pairwise.phylum.env <- new.env()
load("~/pm_workflow/normalize/pm_phylum_pairwise_env.RData", envir=pm.pairwise.phylum.env) # Load the RData into the new environment

# Create experimental class object matching all of the count samples with metadata samples so we can perform high through pairwise between 20 levels (locations)
pm.exp.phylum <- create_expt(count_dataframe=pm.pairwise.phylum.env$tlist[[5]],
                             metadata=pm.pairwise.phylum.env$metadata.subset.rearranged)

# Perform high throughput pairwise comparisons and estimate dispersions
pm.pairwise.phylum<-deseq2_pairwise(input = pm.exp.phylum, batches = 'batch',
                                    model_batch = FALSE, model_cond = TRUE)

#Class
rm(list=ls()[! ls() %in% c("metadata.subset.rearranged","tlist", "seed")])
save(list = ls(), file = "~/pm_workflow/normalize/pm_class_pairwise_env.RData")
pm.pairwise.class.env <- new.env()
load("~/pm_workflow/normalize/pm_class_pairwise_env.RData", envir=pm.pairwise.class.env)

pm.exp.class <- create_expt(count_dataframe=pm.pairwise.class.env$tlist[[1]],
                            metadata=pm.pairwise.class.env$metadata.subset.rearranged)

pm.pairwise.class<-deseq2_pairwise(input = pm.exp.class, batches = 'batch',
                                   model_batch = FALSE, model_cond = TRUE)

#Order
rm(list=ls()[! ls() %in% c("metadata.subset.rearranged","tlist", "seed")])
save(list = ls(), file = "~/pm_workflow/normalize/pm_order_pairwise_env.RData")
pm.pairwise.order.env <- new.env()
load("~/pm_workflow/normalize/pm_order_pairwise_env.RData", envir=pm.pairwise.order.env)

pm.exp.order <- create_expt(count_dataframe=pm.pairwise.order.env$tlist[[4]],
                            metadata=pm.pairwise.order.env$metadata.subset.rearranged)

pm.pairwise.order<-deseq2_pairwise(input = pm.exp.order, batches = 'batch',
                                   model_batch = FALSE, model_cond = TRUE)

#Family
rm(list=ls()[! ls() %in% c("metadata.subset.rearranged","tlist", "seed")])
save(list = ls(), file = "~/pm_workflow/normalize/pm_family_pairwise_env.RData")
pm.pairwise.family.env <- new.env()
load("~/pm_workflow/normalize/pm_family_pairwise_env.RData", envir=pm.pairwise.family.env)

pm.exp.family <- create_expt(count_dataframe=pm.pairwise.family.env$tlist[[2]],
                             metadata=pm.pairwise.family.env$metadata.subset.rearranged)

pm.pairwise.family<-deseq2_pairwise(input = pm.exp.family, batches = 'batch',
                                    model_batch = FALSE, model_cond = TRUE)

#Genus
# Add 1 to every count so we can compute
gplus1<-pm.pairwise.genus.env$tlist[[3]] + 1
gplus1
rm(list=ls()[! ls() %in% c("metadata.subset.rearranged","tlist", "seed", "gplus1")])
save(list = ls(), file = "~/pm_workflow/normalize/pm_genus_pairwise_env.RData")
pm.pairwise.genus.env <- new.env()
load("~/pm_workflow/normalize/pm_genus_pairwise_env.RData", envir=pm.pairwise.genus.env)

pm.exp.genus <- create_expt(count_dataframe=gplus1,
                            metadata=pm.pairwise.genus.env$metadata.subset.rearranged)

pm.pairwise.genus<-deseq2_pairwise(input = pm.exp.genus, batches = 'batch',
                                   model_batch = FALSE, model_cond = TRUE)


# Extract normalized count table from each experimental class object
seed=81
set.seed(seed)
rm(list=ls()[! ls() %in% c("seed", "pm.pairwise.phylum", "pm.pairwise.class", 
                           "pm.pairwise.order", "pm.pairwise.family", "pm.pairwise.genus")])

normalized_phylum <- DESeq2::counts(pm.pairwise.phylum[["run"]], normalized=TRUE)
normalized_phylum_t<-t(normalized_phylum) # Transpose for machine learning ready matrix
dim(normalized_phylum_t)

normalized_class <- DESeq2::counts(pm.pairwise.class[["run"]], normalized=TRUE)
normalized_class_t<-t(normalized_class)
dim(normalized_class_t)

normalized_order <- DESeq2::counts(pm.pairwise.order[["run"]], normalized=TRUE)
normalized_order_t<-t(normalized_order)
dim(normalized_order_t)

normalized_family <- DESeq2::counts(pm.pairwise.family[["run"]], normalized=TRUE)
normalized_family_t<-t(normalized_family)
dim(normalized_family_t)

normalized_genus <- DESeq2::counts(pm.pairwise.genus[["run"]], normalized=TRUE)
normalized_genus_t<-t(normalized_genus)
dim(normalized_genus_t)

rm(list=ls()[! ls() %in% c("seed", "normalized_phylum_t", "normalized_class_t", 
                           "normalized_order_t", "normalized_family_t", "normalized_genus_t")])

taxa_counts_list_normed<-list(normalized_phylum_t, normalized_class_t, normalized_order_t, normalized_family_t, normalized_genus_t)
listNames <- c("Phylum","Class","Order", "Family", "Genus")
names(taxa_counts_list_normed) <- listNames


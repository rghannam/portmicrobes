# This will allow Enrichment Factor (EF) identification and is also a feature selection and dimensionality reduction method
# This is displayed at the taxonomic level of Class but applicable to all levels of resolutions

rm(list=ls())
seed=81
set.seed(seed)
detachAllPackages <- function() {
  
  basic.packages <- c("package:stats","package:graphics","package:grDevices","package:utils","package:datasets","package:methods","package:base")
  
  package.list <- search()[ifelse(unlist(gregexpr("package:",search()))==1,TRUE,FALSE)]
  
  package.list <- setdiff(package.list,basic.packages)
  
  if (length(package.list)>0)  for (package in package.list) detach(package, character.only=TRUE)
  
}

detachAllPackages()

library("rlang", lib.loc="~/R/x86_64-redhat-linux-gnu-library/3.5")
library("ggplot2", lib.loc="~/R/x86_64-redhat-linux-gnu-library/3.5")
library("Rcpp", lib.loc="/usr/lib64/R/library")
library("ggplot2", lib.loc="~/R/x86_64-redhat-linux-gnu-library/3.5")
library("vegan")
library("csv")
library("plotly")
library("plyr")
library("dplyr")

seed=81
set.seed(seed)

rm(list=ls())
setwd("~/pm_workflow")


# Generate a taxonomic class count matrix from all reads
pkgs <- c("phyloseq", "ggplot2", "rlang", "Rcpp", "vegan", "csv", "plotly", "plyr", "dplyr", "reshape", "reshape2") # load packages
invisible(sapply(pkgs,require, character = TRUE))

# Load the final count table from data_preprocessing.R
# These counts were subsetted pertaining to those used in analysis of this study (1218 surface water samples)
pm.counts<-'https://github.com/rghannam/pm_workflow/raw/master/data/counts_taxonomy/pm_counts1218.rds'
pm.counts<-readRDS(url(pm.counts))
dim(pm.counts)

pm.taxa<-'https://github.com/rghannam/pm_workflow/raw/master/data/counts_taxonomy/pm_taxonomy.rds'
pm.taxa<-readRDS(url(pm.taxa))
dim(pm.taxa)

pm.metadata<-'https://github.com/rghannam/pm_workflow/raw/master/data/metadata/pm_metadata.csv'
pm.metadata<-read.csv(url(pm.metadata))
dim(pm.metadata)

# Manipulate data so we can convert to a phyloseq object properly
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
counts = otu_table(pm.counts, taxa_are_rows = FALSE) # leave false because we arbitrarily set column names with colnames()->NULL, we can always map back to respective ASV
taxa = tax_table(pm.taxa_re_t)
sampledata = sample_data(pm.metadata)
pm.ps = phyloseq(counts, taxa, sampledata)
pm.ps
otu_table(pm.ps)[1:5, 1:5]

# Exclude the uncharacterized or NA Classes
c<-subset_taxa(pm.ps, !is.na(Class) & !Class %in% c("", "uncharacterized"))
x<-list()
alltaxa<-list(c, x)
listNames <- c("class")
names(alltaxa) <- listNames

# Calculate prevalence and abundance of each Class and how many ASVs were assigned to each Class
prev.abund.class = apply(X = otu_table(alltaxa$class),
                         MARGIN = ifelse(taxa_are_rows(alltaxa$class), yes = 1, no = 2),
                         FUN = function(x){sum(x > 0)})
prev.abund.class = data.frame(Prevalence = prev.abund.class,
                              TotalAbundance = taxa_sums(alltaxa$class),
                              tax_table(alltaxa$class))

pm.totalavgprev.class <- plyr::ddply(prev.abund.class, "Class", function(df1){cbind(mean(df1$Prevalence),sum(df1$Prevalence))})
as.data.frame(pm.totalavgprev.class)
names(pm.totalavgprev.class) <- c("Class", "mean", "total")


# Subset the counts containing only those that were within >=15 samples
pm.totalavgprev.class.subset <- subset(pm.totalavgprev.class, pm.totalavgprev.class$total <= 15) # Filter entries with <=15 prevalence
extract.class<-as.data.frame(pm.totalavgprev.class.subset[,1])
names(extract.class) <- ("Class") 
extract.class <- as.list(as.vector(t(extract.class)))

filterclass = (extract.class)
pm.ps.class.sub = subset_taxa(alltaxa$class, !Class %in% filterclass)
pm.ps.class.sub
alltaxa$class

length(get_taxa_unique(alltaxa$class, taxonomic.rank = "Class")) # 143
length(get_taxa_unique(pm.ps.class.sub, taxonomic.rank = "Class")) # 88
# Create a phyloseq object agglomerated by class
pm.ps.class = tax_glom(pm.ps.class.sub, "Class", NArm = TRUE, bad_empty=c(NA, "", " ", "\t"))
pm.ps.class

# Create a count matrix of Class feature
pm.ps.class.melt <- psmelt(pm.ps.class)
pm.ps.class.melt$Class <- as.character(pm.ps.class.melt$Class)
pm.ps.class.melt
pm.ps.class.melt <- aggregate(Abundance~Sample+Class, pm.ps.class.melt, FUN=sum)
pm.ClassMatrix <- cast(pm.ps.class.melt, Sample ~ Class)
dim(pm.ClassMatrix)
head(pm.ClassMatrix)
pm.ClassMatrix$Sample<-gsub("-", "", pm.ClassMatrix$Sample) # Remove - from sample_names

# pm.ClassMatrix is the object that will be used for DA analysis

seed=81
set.seed=seed
rm(list=ls()[! ls() %in% c("seed", "pm.ClassMatrix")])
pkgs <- c("hpgltools", "DESeq2", "plyr", "dplyr", "dbplyr","stringr",
          "data.table", "tidyverse", "tibble") 
invisible(sapply(pkgs,require, character = TRUE))

# This was for importing all taxonomic levels from source but here we will just use Class
taxa.import<-function(x, transpose = 1, ...){
  alist<-sapply(x, read.csv, simplify = FALSE, USE.NAMES = TRUE)
  stopifnot(is.list(alist))
  message("Importing taxonomy count table")
  listNames <- c("Class")
  names(alist) <- listNames
  for(i in 1:length(alist)){ # drops first Sample column, since SampleID is stored as row.name, it's fine.
    rownames(alist[[i]]) <- alist[[i]][,1]
    alist[[i]] <- alist[[i]][,-1,drop = FALSE]
  }
  colnames(alist$Class)<-NULL
  c<-c.mat.t<-t(alist$Class) 
  c<-as.data.frame(c.mat.t)
  alist<-list()
  alist<-c(list(c), alist)
  listNames <- c("Class")
  names(alist) <- listNames
  return(alist)
  return(c)
}

# Import pm.ClassMatrix
class_level<-taxa.import('https://github.com/rghannam/pm_workflow/raw/master/data/differential_abundance/pm.classmatrix.csv') # puts them in proper order?
row.names(class_level[[1]])
colnames(class_level[[1]])

# Import metadata
pm.metadata<-'https://github.com/rghannam/pm_workflow/raw/master/data/metadata/pm_metadata.csv'
pm.metadata<-read.csv(url(pm.metadata))
dim(pm.metadata)

# Import counts
pm.counts<-'https://github.com/rghannam/pm_workflow/raw/master/data/counts_taxonomy/pm_counts1218.rds'
pm.counts<-readRDS(url(pm.counts))
dim(pm.counts)

rownames(pm.counts)
pm.metadata$sample_name
row.names(pm.counts)<-gsub("-", "", row.names(pm.counts))
pm.metadata$sample_name<-gsub("-", "", pm.metadata$sample_name)

# Arrange metadata
row.names(pm.metadata)<-pm.metadata$sample_name
dim(pm.counts)

# Match the samples from counts to our metadata
samplenames<-rownames(pm.counts)
pm.metadata.arrange<-pm.metadata[match(samplenames, pm.metadata$sample_name),]

# Append metadata to subset our main count table per sample filtering criteria
colnames(pm.metadata.arrange) # The sample variables we can subset our main count table by
names(pm.metadata.arrange)[names(pm.metadata.arrange) == 'sequence_run'] <- 'batch' # Rename 'sequence_run' to 'batches' for experimental class object/DA analysis
pm.metadata.arrange <- pm.metadata.arrange[pm.metadata.arrange$sample_type == "open water" & 
                                             pm.metadata.arrange$reads_passed_filter_no_chimeras_overthousand == "yes" & 
                                             pm.metadata.arrange$season == 'Summer', ]
as.data.frame(rownames(pm.metadata.arrange))
length(unique(pm.metadata.arrange$sample_name))
allrownames<-as.data.frame(class_level[[1]])

tt<-t(class_level[[1]])
allrownames<-row.names(tt)
pm.metadata.arrange2=cbind(pm.metadata.arrange,allrownames) # Cbind metadata to the rownames of the class_level object

library(dbplyr)
class(pm.metadata.arrange2)
pm.metadata.arrange2$sample_name<-as.factor(pm.metadata.arrange2$sample_name)
levels(pm.metadata.arrange2$sample_name)

metadata.subset.work<-pm.metadata.arrange2
library(dplyr)
metadata.subset.rearranged<-metadata.subset.work %>% arrange(sample_name, desc(allrownames))
metadata.subset.rearranged$allrownames<-NULL

# Remove the taxonomic Class variables and rename them to X.., we will apply their respective names back to them later on
dim(metadata.subset.rearranged) # Sanity check: 1218 observations
metadata.subset.rearranged$geo_port<-as.factor(metadata.subset.rearranged$geo_port)
metadata.subset.rearranged$batch<-as.factor(metadata.subset.rearranged$batch)
metadata.subset.rearranged$sample_name<-as.factor(metadata.subset.rearranged$sample_name)
names(metadata.subset.rearranged)[1]<-"sampleid" # Rename for DA analysis
names(metadata.subset.rearranged)[37]<-"condition" # Rename for DA analysis

seed=81
set.seed(seed)
# Create an experimental class object at taxonomy: Class
rm(list=ls()[! ls() %in% c("metadata.subset.rearranged","class_level", "seed")]) # Clear memory
save(list = ls(), file = "~/pm_workflow/pm.class.pairwise.env.RData") # Save these 2 objects as RData
pm.pairwise.class.env <- new.env()
load("~/pm_workflow/pm.class.pairwise.env.RData", envir=pm.pairwise.class.env) # Load the RData into the new environment

require(hpgltools)
# Create experimental class object matching all of the count samples with metadata samples so we can perform high through pairwise between all 20 levels, or locations
pm.exp.class <- create_expt(count_dataframe=pm.pairwise.class.env$class_level[[1]],
                            metadata=pm.pairwise.class.env$metadata.subset.rearranged)

# Perform high throughput pairwise comparisons and estimate dispersions
pm.pairwise.class<-deseq2_pairwise(input = pm.exp.class, batches = 'batch',
                                   model_batch = FALSE, model_cond = TRUE)

rm(list=ls()[! ls() %in% c("pm.pairwise.class", "seed")]) # Clear memory
pkgs <- c("RG2", "hpgltools", "DESeq2", "dplyr", "plyr", "dbplyr","stringr",
          "data.table", "tidyverse", "tibble") 

# Run the above (lines 224-235) to obtain pm.pairwise.class object or just call from GitHub
pm.pairwise.class<-'https://github.com/rghannam/pm_workflow/raw/master/data/differential_abundance/pm.pairwise.class.RData'
pm.pairwise.class<-readRDS(url(pm.pairwise.class))

pwise.list.class<-as.list(pm.pairwise.class$all_tables)
pwise.list.class<- mapply(`[<-`, pwise.list.class, 'Pairwise', value = names(pwise.list.class), SIMPLIFY = FALSE)
pwise.list.class2<- mapply(`[<-`, pwise.list.class, 'Numerator', value = names(pwise.list.class), SIMPLIFY = FALSE)
pwise.list.class3<- mapply(`[<-`, pwise.list.class2, 'Denominator', value = names(pwise.list.class2), SIMPLIFY = FALSE)
#  Stores all the elements in the list in one dataframe, also appends rowname (Pairwise) to the first column
#  Our observations are currently a pairwise comparison for every location (20) vs every location per class (88) for 16720 observations.
pwise.class<-do.call(rbind.data.frame, pwise.list.class3)

colnames(pwise.class)
pwise.class.subset<-pwise.class
require(data.table)

#  Set to data.table for regular expression manipulation
setDT(pwise.class.subset, keep.rownames = TRUE)[]
class(pwise.class.subset)
colnames(pwise.class.subset)[1] <- "class.Var" # Rownames to columns
pwise.class.subset$class.Var<-gsub("^.*[A-Z]","", pwise.class.subset$class.Var) # Remove character string, keeps numerical which pertains to which Class it is (from our Class matrix)
pwise.class.subset$Numerator<-gsub("_.*", "", pwise.class.subset$Numerator)
pwise.class.subset$Denominator<-gsub(".*_", "", pwise.class.subset$Denominator)
#  Remove any string after '_vs' in the variable 'Numerator', we need 'Numerator' and 'Denominator' variables later on
#  To grep later to find prevalence of signatures.

#  We want to replace the numerical X1-X88 in pwise.class.subset$class.Var with the actual taxonomic Class they correspond to
class.matrix<-'https://github.com/rghannam/pm_workflow/raw/master/data/differential_abundance/pm.classmatrix.csv'
class.matrix<-read.csv(url(class.matrix))
dim(class.matrix)
row.names(class.matrix)<-class.matrix$Sample
class.matrix[1]<-NULL
class.matrix<-t(class.matrix)

#  Take index of the counts to append later
classInd <- arrayInd(1:88, dim(class.matrix)) # Index the class
classInd

classInd<-rownames(class.matrix)[classInd[,1]] # Extract rownames
classInd
classIndwork<-classInd
classIndwork<-as.data.frame(classIndwork)
colnames(classIndwork)[1] <- "class.index" # Change name of variable
classIndwork<-cbind(classed.variable = 1:88, classIndwork) # Class the count.index variable


#  Store all of the variables in a new object (x) so values are unperturbed moving forward
x <- pwise.class.subset[,2:10, drop=FALSE]
pwise.class.subset[] <- setNames(classIndwork$class.index, classIndwork$classed.variable)[as.character(unlist(pwise.class.subset))] # Replace numeric in 'class.Var' with actual corresponding Class from prior indexing
pwise.class.subset<-cbind(x, pwise.class.subset) # cbind all of the original values stored in x back to main dataframe

class.df <- pwise.class.subset[,1:10, drop=FALSE] # Returns data.frame with values converted to Class names
head(class.df, 10)
#  Returns a data.frame of all pairwise comparisons per location per class, along with statistics from differential abundances for subsequent feature selection
#  We can use this dataframe for further subsetting and downstream analytics

lengthclass<-length(unique(class.df$class.Var)) # Sanity check for all class
lengthclass # Should be 88

#  We need to subset dataframes based on Numerator and Denominator
#  This can be done based on +/- logFC values (+ = Num,  - = Denom) (based on the experimental design)

detachAllPackages <- function() {
  
  basic.packages <- c("package:stats","package:graphics","package:grDevices","package:utils","package:datasets","package:methods","package:base")
  
  package.list <- search()[ifelse(unlist(gregexpr("package:",search()))==1,TRUE,FALSE)]
  
  package.list <- setdiff(package.list,basic.packages)
  
  if (length(package.list)>0)  for (package in package.list) detach(package, character.only=TRUE)
  
}

detachAllPackages() # (unloads all but base R) because tidyverse and plyr mess with tibble functions
require(tidyverse)
class.df.subsetNUM <- class.df[class.df$logFC >0, ] # Subset 'Numerator'
#  Rearrange in alphabetical Class for downstream processing such as appending abundance sums
class.df.subsetNUM<-class.df.subsetNUM %>% arrange((Numerator)) # Arrange in alphabetical Class and the same way as downstream summary data.frames will be returned as
class.df.subset.num<-class.df.subsetNUM %>% 
  as.tibble() %>% 
  count(Numerator, class.Var)

#  class.df.subset.num has the Prevalence (n) of how many times a particular class was in a higher differential abundance against another location
#  e.g. class.df.subset.num, Acidobacteria was higher in Wilhelmshaven against 14 other locations

# Do the denominator now
class.df.subsetDENOM <- class.df[class.df$logFC <0, ] # Subset 'Denominator'
class.df.subsetDENOM<-class.df.subsetDENOM %>% arrange((Denominator))
class.df.subset.denom<-class.df.subsetDENOM %>% 
  as.tibble() %>% 
  count(Denominator, class.Var)

names(class.df.subset.num)[names(class.df.subset.num) == 'Numerator'] <-'Location'
names(class.df.subset.denom)[names(class.df.subset.denom) == 'Denominator'] <-'Location'

#  The logic behind this is that a +logFC == Numerator (location) has a higher differential abundance for that Class than it does its Denominator (location).
#  Likewise, -logFC == Denominator has a higher differential abundance than its Numerator
#  If we are separating by these values we then have two separate data.frames which effectively contain the location and the count it has of being + or - against another location, and for each class (signature)
#  In either data.frame, those counts (n) are accurate for true cases where that class is differentially expressed higher in that location
#  e.g. in either num/denom data.frames, the counts are basically all = +logFC.
# class.df.subsetNUM = all +logFC (+logFC remains +)
# class.df.subsetDENOM = all -logFC (since these are denominators, any -logFC is basically +) because of the fact we're treating them out as separate data

head(class.df.subset.denom)

#  Now we need to append abundances per Class
#  We need a vector of 'sampleid' from our class matrix so we can set conditions per location - create a location index
#  e.g. Baltimore1-Baltimore30=BALT
#  Then summarize all BALT for each Class, append those abundances to class.df.subsetNUM/class.df.subsetDENOM

class.matrix2<-'https://github.com/rghannam/pm_workflow/raw/master/data/differential_abundance/pm.classmatrix.csv'
class.matrix2<-read.csv(url(class.matrix2))
dim(class.matrix2)
row.names(class.matrix2)<-class.matrix2$Sample
class.matrix2[1]<-NULL

sampleidInd <- arrayInd(1:1218, dim(class.matrix2)) # Index the class
sampleidInd
sampleidInd<-rownames(class.matrix2)[sampleidInd[,1]] # Extract rownames
sampleidInd<-as.data.frame(sampleidInd)
sampleidInd$sampleidInd<-gsub('[0-9]+', '', sampleidInd$sampleidInd) # Regex to remove any numerical after character string

library(stringr)
#  Remove any unwanted patterns in the strings from 'sampleid' factor
sampleidInd<-unlist(sampleidInd) # Convert to vector for stringr input
str(sampleidInd)
sampleidInd<-str_remove(sampleidInd, "post")
sampleidInd<-str_remove(sampleidInd, "point")
sampleidInd<-str_extract(sampleidInd, "^.{1,3}")
#  Keep at maximum 3 characters in string
sampleidInd<-as.data.frame(sampleidInd)
str(sampleidInd) # 1218 observations at 20 levels, what we want
length(unique(sampleidInd$sampleidInd)) # Length matches proper levels for locations (20)

colnames(sampleidInd) <- c("location.index") # Rename column to location index
class.matrix2$location.index = sampleidInd$location.index # Send the index to the Class matrix
#  We've appended the location index to the Class matrix

require(reshape2)
require(dplyr)
require(tidyr)

#  Summarise all columns by groups within a factor (location.index)
#  This gives us all of the sums per location for each Class, now we can append to NUM/DENOM dfs.
#  Note: if there's a Class for location X in Numerator, that class won't be in location X in Denominator
summarise.class<-class.matrix2 %>% 
  group_by(location.index) %>% 
  summarise_all(funs(sum))
summarise.class<-as.data.frame(summarise.class)

#  Make a column rownames
rownames(summarise.class) <- summarise.class$location.index
summarise.class <- data.frame(summarise.class[,-1])
summarise.class.t<-t(summarise.class)
summarise.class.t<-as.data.frame(summarise.class.t)
#  summarise.class.t now has rownames of Class and variables as locations
#  We can rename the variables to match NUM/DENOM and separate into a NUM/DENOM df

#  Note: we haven't filtered by logFC, basemean or any of the stat variables yet
#  Find unique levels to Denominator factor
length(unique(class.df.subsetDENOM$Denominator)) # 19 locations in 'Denominator'
length(unique(class.df.subsetNUM$Numerator)) # 19 locations in 'Numerator'
length(unique(class.df.subsetDENOM$class.Var)) # 88 class in 'Denominator'
length(unique(class.df.subsetNUM$class.Var)) # 88 class in 'Numerator'

#  Remame column names from Class matrix to match how they are in class.df.subsetNUM/DENOM
#  Although Numerators don't have Baltimore and start with Busan, they're technically in different alphabetical Class, we will handle this downstream.
sum.phy<-summarise.class.t
colnames(sum.phy) <- c("Baltimore", "Busan", "Charleston", "Duluth", "Martigues", "Naples", "Rotterdam", "Venice", "Wilhelmshaven", "Galveston", "GreenBay", "HongKong", "Keweenaw", "LosAngeles", "NewOrleans", "Norfolk", "NewYork", "Oakland", "Seattle", "Singapore")
# ^ Rename columns to respective location

# Arrange df vars by position - landroni
# 'vars' must be a named vector, e.g. c("var.name"=1)
arrange.vars <- function(data, vars){
  # stop if not a data.frame (but should work for matrices as well)
  stopifnot(is.data.frame(data))
  
  # sort out inputs
  data.nms <- names(data)
  var.nr <- length(data.nms)
  var.nms <- names(vars)
  var.pos <- vars
  # sanity checks
  stopifnot( !any(duplicated(var.nms)), 
             !any(duplicated(var.pos)) )
  stopifnot( is.character(var.nms), 
             is.numeric(var.pos) )
  stopifnot( all(var.nms %in% data.nms) )
  stopifnot( all(var.pos > 0), 
             all(var.pos <= var.nr) )
  
  # prepare output
  out.vec <- character(var.nr)
  out.vec[var.pos] <- var.nms
  out.vec[-var.pos] <- data.nms[ !(data.nms %in% var.nms) ]
  stopifnot( length(out.vec)==var.nr )
  
  # re-arrange vars by position
  data <- data[ , out.vec]
  return(data)
}

# Rearrange variables of Class matrix to the same Class as they are in class.df.subsetNUM/DENOM
# Note: Denominator lacks Wilhelmshaven and Numerator lacks Baltimore, hence the starting order
# Rearrange the Denominator
sum.phy.rearranged_den<-arrange.vars(sum.phy, c("Baltimore"=1, "Busan"=2, "Charleston"=3, 
                                            "Duluth"=4,"Galveston"=5,"GreenBay"=6,"HongKong"=7,
                                            "Keweenaw"=8,"LosAngeles"=9,"Martigues"=10,"Naples"=11,
                                            "NewOrleans"=12,"NewYork"=13, "Norfolk"=14,"Oakland"=15,"Rotterdam"=16,
                                            "Seattle"=17,"Singapore"=18,"Venice"=19))

names(sum.phy.rearranged_den)
rownames(sum.phy.rearranged_den)

# Rearrange the Numerator
sum.phy.rearranged_num<-arrange.vars(sum.phy, c("Busan"=1, "Charleston"=2, 
                                            "Duluth"=3,"Galveston"=4,"GreenBay"=5,"HongKong"=6,
                                            "Keweenaw"=7,"LosAngeles"=8,"Martigues"=9,"Naples"=10,
                                            "NewOrleans"=11,"NewYork"=12, "Norfolk"=13,"Oakland"=14,"Rotterdam"=15,
                                            "Seattle"=16,"Singapore"=17,"Venice"=18,"Wilhelmshaven"=19))

names(sum.phy.rearranged_num)
rownames(sum.phy.rearranged_num)

num<-class.df.subsetNUM
denom<-class.df.subsetDENOM

# Extract out Denominator/class.Var only
num.sub<-as.data.frame(num[,c("Numerator","class.Var")])
denom.sub<-as.data.frame(denom[,c("Denominator","class.Var")])
num.names<-as.vector(num.sub$Numerator)
denom.names<-as.vector(denom.sub$Denominator)
class.names.num<-as.vector(num.sub$class.Var)
class.names.denom<-as.vector(denom.sub$class.Var)

# Makes rownames the locations from denominator,numerator
rownames(num.sub) <- make.names(num.sub[,1], unique = TRUE) 
rownames(denom.sub) <- make.names(denom.sub[,1], unique = TRUE) 
# We can regex the numericals out of the rownames, and now we can remove the column: Denominator
keep <- ("class.Var")

# Keep class.var
num.sub <- subset(num.sub, select = keep) 
num.sub <- num.sub[keep] 
denom.sub <- subset(denom.sub, select = keep) 
denom.sub <- denom.sub[keep] 

# Transpose the denominators from above code into column names
num.dat<-as.data.frame(num.sub)
denom.dat<-as.data.frame(denom.sub)
N<-num.dat
D<-denom.dat

# Remove numerical and . from column names
#Num
names(num.dat) <- gsub("[[:digit:]]", "", names(num.dat) )
names(num.dat) <- gsub("\\.", "", names(num.dat) )
names(num.dat)

#Denom
names(denom.dat) <- gsub("[[:digit:]]", "", names(denom.dat) )
names(denom.dat) <- gsub("\\.", "", names(denom.dat) )
names(denom.dat)

#  Rename variables and sanity checks
#Num
rownames(num.dat) <- make.names(num.dat[,1], unique = TRUE) # Renamed to classVar from class.Var
num.bind.2<-cbind(num.dat,num.names)
class(num.bind.2)
str(num.bind.2)
num.bind.2$classVar<-as.factor(num.bind.2$classVar)
class(num.bind.2$classVar)
levels(num.bind.2$classVar)# 88 levels (same as our class matrix)
main.num<-num.bind.2 # Working Numerator data.frame

#Denom
rownames(denom.dat) <- make.names(denom.dat[,1], unique = TRUE)
denom.bind.2<-cbind(denom.dat,denom.names)
class(denom.bind.2)
str(denom.bind.2)
denom.bind.2$classVar<-as.factor(denom.bind.2$classVar)
class(denom.bind.2$classVar)
levels(denom.bind.2$classVar)
main.denom<-denom.bind.2 # Working Denominator data.frame

# Obtain levels of 'Numerator'
levels.num<-levels(main.num$num.names)
uniq.num <- unique(main.num$num.names)

# Obtain levels of 'Denominator'
levels.denom<-levels(main.denom$denom.names)
uniq.denom <- unique(main.denom$denom.names)

#  What variables are present in our main Class matrix?
colnames.class<-colnames(sum.phy.rearranged_num)
colnames.class

# How does it match denom?
levels(main.num$num.names) # No Baltimore in Numerator
levels(main.denom$denom.names) # No Wilhelmshaven in Denominator this is how we know how to varrange.


# Lets separate out the values belonging to numerator from demominator
# Effectively will have each Class and its correpsonding location for which there was a pairwise comparison made
# These should sum to all pairwise comparisons 16720
# Store counts of num/denom in lists
vec_countsD<-list()
for (i in 1:20) {
  vec_countsD[[i]]<-sum.phy.rearranged_den[!is.na(sum.phy.rearranged_den[[i]]), i]
}

vec_countsN<-list()
for (i in 1:20) {
  vec_countsN[[i]]<-sum.phy.rearranged_num[!is.na(sum.phy.rearranged_num[[i]]), i]
}


# Take denom levels
denom_lev<-list()
for (i in 1:length(levels.denom)){
  denom_lev[[i]]<-subset(main.denom, denom.names == uniq.denom[i])
}

# Store values (vec_countsD) from denom levels.
dlist = list()
for (n in 1:length(levels.denom)){
  dlist[[n]]<-data.frame(sapply(denom_lev[[n]], factor, levels=classInd, labels = vec_countsD[[n]]))
}

detachAllPackages()

library(plyr)
dlist2<-rbind.fill(dlist)
denom_lev2<-rbind.fill(denom_lev) # 9800

# Take num levels
num_lev<-list()
for (i in 1:length(levels.num)){
  num_lev[[i]]<-subset(main.num, num.names == uniq.num[i])
}
nlist = list()
for (n in 1:length(levels.num)){
  nlist[[n]]<-data.frame(sapply(num_lev[[n]], factor, levels=classInd, labels = vec_countsN[[n]]))
}

library(plyr)
nlist2<-rbind.fill(nlist)
num_lev2<-rbind.fill(num_lev) #6920

colnames(nlist2)[colnames(nlist2)=="classVar"] <- "Abundance.summed"
nlist2$num.names<-NULL
colnames(dlist2)[colnames(dlist2)=="classVar"] <- "Abundance.summed"
dlist2$denom.names<-NULL

#  Still keeping Numerators/Denominators separate, append statistics from our differential expression table
#Num
class.df.num<-cbind(class.df.subsetNUM, nlist2)
class.df.num

#Denom
class.df.denom<-cbind(class.df.subsetDENOM, dlist2)
class.df.denom
str(class.df.denom)

#  Extract Numerator/Denominator vectors
#Num
numInd<-class.df.num$Numerator
numInd<-as.data.frame(numInd)
length(unique(numInd$numInd)) # 19 locations in num
#Denom
denomInd<-class.df.denom$Denominator
denomInd<-as.data.frame(denomInd)
length(unique(denomInd$denomInd)) # 19 locations in denom

length(unique(rownames(summarise.class))) # 20 locations
dim(summarise.class)

detachAllPackages() #  There are issues running tidyverse and plyr functions in sequence
require(tidyverse)

prev.num<-class.df %>% 
  as.tibble() %>% 
  count(Numerator, class.Var)

prev.denom<-class.df %>% 
  as.tibble() %>% 
  count(Denominator, class.Var)

library(plyr)

#  Summarise colums per location
sum.num<-ddply(prev.num, "Numerator", numcolwise(sum)) # e.g. How many times this Numerator appears in class.df
sum(sum.num$n) # Should = the length of obs in class.df of 16720

sum.denom<-ddply(prev.denom, "Denominator", numcolwise(sum))
sum(sum.denom$n)

# Filter by some desired parameters based on statistic variables
# Here, we subset by logFC and adjP
colnames(class.df)
detachAllPackages()
require(tidyverse)
#Num
class.df.num.2 <- class.df.num[class.df.num$logFC >= 2 &
                                 class.df.num$adj.P.Val <= 0.05, ]
class.df.num.3<-class.df.num.2 %>% 
  as.tibble() %>% 
  count(Numerator, class.Var, Abundance.summed)

#Denom
# Note: we used <= -2, but this is technically still 'logFC 2' for positive enrichment
class.df.denom.2 <- class.df.denom[class.df.denom$logFC <= -2 &
                                     class.df.denom$adj.P.Val <= 0.05, ]
class.df.denom.3<-class.df.denom.2 %>% 
  as.tibble() %>% 
  count(Denominator, class.Var, Abundance.summed)

#  Rename columns
colnames(class.df.num.3)[colnames(class.df.num.3)=="n"] <- "Prevalence"
colnames(class.df.denom.3)[colnames(class.df.denom.3)=="n"] <- "Prevalence"
class(class.df.denom.3)
#  How many levels passed our filters?
unique.levels <- length(unique(class.df.num.3$Numerator))
unique.levels # 19 locations passed
unique.levels <- length(unique(class.df.denom.3$Denominator))
unique.levels # 19 locations passed

unique.levels.d <- length(unique(class.df.denom.3$class.Var))
unique.levels.d # 56/88 passed filter.

unique.levels.n <- length(unique(class.df.num.3$class.Var))
unique.levels.n # 50/88 passed filter.

Numerator<-class.df.num.3
Denominator<-class.df.denom.3
names(Numerator)[names(Numerator) == 'Numerator'] <-'Location'
names(Denominator)[names(Denominator) == 'Denominator'] <-'Location'
class.df.num.denom <- rbind(Numerator, Denominator) # Append both num/denom dfs

library(data.table) 
class.df.num.denom <- data.table(class.df.num.denom) # Convert to data.table to group and summarize
detachAllPackages()
library(dplyr)
class.df.num.denom<-class.df.num.denom
class.df.num.denom[,3] <- as.numeric(as.matrix(class.df.num.denom[,3])) # Originally a factor

#  Summarize prevalence of all the conditions matching identical Location and Class variables.
#  e.g.) if Charleston was +logFC for Class: Firmicutes over 2 other Locations in Numerator
#  and was negative (which is technically +logFC) against 1 other location in Denominator, it sums to 3 while preserving 'Abunance.summed'
sum.sigs.class<-class.df.num.denom %>% group_by(Location, class.Var, Abundance.summed) %>% summarize(Prevalence = sum(Prevalence))

#  Create average abundance vector as an extra filtering parameter to be added to our final dataframe
colnames(sampleidInd)
percentage <- prop.table(table(sampleidInd$location.index)) * 100
getfreq<-cbind(freq=table(sampleidInd$location.index), percentage=percentage)
getfreq<-as.data.frame(getfreq)
row.names(getfreq) <- c("Baltimore", "Busan", "Charleston", "Duluth", "Martigues", "Naples", "Rotterdam", "Venice", "Wilhelmshaven", "Galveston", "GreenBay", "HongKong", "Keweenaw", "LosAngeles", "NewOrleans", "Norfolk", "NewYork", "Oakland", "Seattle", "Singapore")
freq <- c("freq")
getfreq<-getfreq[,freq]
getfreq<-as.data.frame(getfreq)
getfreq.labels.vec<-dplyr::pull(getfreq, getfreq)

sum.sigs.class2<-sum.sigs.class
keep.x <- c("Location")
sum.sigs.class2<-sum.sigs.class2[,keep.x]
sum.sigs.class2<-as.data.frame(sum.sigs.class2)
# sum.sigs2 is now a vector of all location in class of sum.sigs, where we need to apply 'nsamples' variable

sum.sigs.class3<-data.frame(sapply(sum.sigs.class2, factor, 
                                   levels=c("Baltimore", "Busan", "Charleston", "Duluth", 
                                            "Martigues", "Naples", "Rotterdam", "Venice", "Wilhelmshaven", 
                                            "Galveston", "GreenBay", "HongKong", "Keweenaw", 
                                            "LosAngeles", "NewOrleans", "Norfolk", "NewYork", "Oakland", 
                                            "Seattle", "Singapore"),                
                                   labels = c(getfreq.labels.vec)))
class(sum.sigs.class3)
sum.sigs.class.work<-sum.sigs.class
class(sum.sigs.class.work)
sum.sigs.class.work<-as.data.frame(sum.sigs.class.work)
class(sum.sigs.class.work)
sum.sigs.class4<-cbind(sum.sigs.class.work, sum.sigs.class3)
colnames(sum.sigs.class4)
colnames(sum.sigs.class4)[5] <- "nsamples"

head(sum.sigs.class4)
class(sum.sigs.class4$nsamples)
sum.sigs.class4[,5] <- as.integer(as.matrix(sum.sigs.class4[,5]))
head(sum.sigs.class4)
class(sum.sigs.class4$nsamples)
require(dplyr)
sum.sigs.class4<-tbl_df(sum.sigs.class4) # Make it into a tbl class
class(sum.sigs.class4)
# Define the new columns
sum.sigs.class4<-sum.sigs.class4%>%
  mutate(avgabund=Abundance.summed/nsamples)
head(sum.sigs.class4) # This has a working dataframe of prevalence/abundance.summed/prevalence/avg abundance.

length(unique(sum.sigs.class4$class.Var))

# Which didn't pass/make filter?
unique.sum.sigs.FINAL <- length(unique(sum.sigs.class4$class.Var))
unique.sum.sigs.FINAL # 56/88 passed filter.
enrichment.factors<-sum.sigs.class4
enrichment.factors$class.Var<-as.factor(enrichment.factors$class.Var)
levels(enrichment.factors$class.Var)
write.csv(enrichment.factors, "~/pm_workflow/enrichment.factors.csv")
# ^ enrichment.factors is used for Enrichment Factor (EF) assignment and dendrogram figure

class.passed.filter<-as.character(unique(enrichment.factors$class.Var))

# Main frame
summarise.class.t2<-summarise.class.t
rownamesallclass<-rownames(summarise.class.t2)

classdifference<-setdiff(rownamesallclass,class.passed.filter)

classdifference # These are the taxonomic Classes that did not make it past enrichment factor filtering criteria


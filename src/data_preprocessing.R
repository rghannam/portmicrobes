# ---- Data Preprocessing ----

seed=81
rm(list=ls()[! ls() %in% c("seed")])
set.seed(seed)
setwd("~/pm_workflow/") # CHANGE ME - this will be working directory

pm.counts<-'https://github.com/rghannam/pm_workflow/raw/master/data/counts_taxonomy/pm_counts.rds'
pm.counts<-readRDS(url(pm.counts))
dim(pm.counts)

pm.taxa<-'https://github.com/rghannam/pm_workflow/raw/master/data/counts_taxonomy/pm_taxonomy.rds'
pm.taxa<-readRDS(url(pm.taxa))
dim(pm.taxa)

pm.metadata<-'https://github.com/rghannam/pm_workflow/raw/master/data/metadata/pm_metadata.csv'
pm.metadata<-read.csv(url(pm.metadata))
dim(pm.metadata)

row.names(pm.metadata)<-pm.metadata$sample_name
dim(pm.counts)

# Need to match the samples from count table to our metadata
# Create vector of sample names from count table

samplenames<-rownames(pm.counts)
pm.metadata.arrange<-pm.metadata[match(samplenames, pm.metadata$sample_name),]

options(max.print=2500)
row.names(pm.counts)<-gsub("-", "", row.names(pm.counts)) # gsub to remove any - for downstream analysis
rownames(pm.counts)
pm.counts[,1514] # rename check

# Append metadata to subset our main count table as we please
pm.counts2=cbind(pm.metadata.arrange,pm.counts)
# Filter/subset data based on metadata criteria 
colnames(pm.metadata.arrange) # The sample variables we can subset our main count table by
colnames(pm.counts2[1:5])

# We are particular about using open water samples for our downstream analysis
# We also want samples with reads that passed the filter > 1000
pm.counts.subset <- pm.counts2[pm.counts2$sample_type == "open water" & 
                                 pm.counts2$reads_passed_filter_no_chimeras_overthousand == "yes" & 
                                 pm.counts2$season == 'Summer', ]
as.data.frame(rownames(pm.counts.subset))
length(unique(pm.counts.subset$sample_name)) # Ensure we have 20 levels to our location factor

# Remove metadata from main count table
x<-c(157987:158034) # Define range of columns to remove (keep everything but appended metadata)
y<-c(157987:158034)

index = which(x %in% y ) # Index of position of values 'x', 'y'
index
pm.counts.subset = pm.counts.subset[-index]

dim(pm.counts.subset) # Should have same column dimensions as original counts: 1218x157986
pm.counts.subset<-as.matrix(pm.counts.subset)
colnames(pm.counts.subset[1:2])

# Removal of singletons and zero sum features
# After subsetting from metadata, there will be ASVs from other samples that don't belong to open water, such as a swab of the object, for instance.
# Remove singletons
i <- (colSums(pm.counts.subset, na.rm=T) != 1) # T if colSum is not 0, F otherwise
pm.counts_nz <- pm.counts.subset[, i] # these are all the ASVs that are >1 and in >1 sample.
pm.counts_z <- pm.counts.subset[, !i] # this is a matrix of where an ASV is only in quantity "1" for only 1 sample
dim(pm.counts_nz)
dim(pm.counts_z)

# Remove all zero sum columns
i <- (colSums(pm.counts_nz, na.rm=F) != 0) # T if colSum is not 0, F otherwise
pm.counts_nz_F <- pm.counts_nz[, i] # these are all the ASVs that are >1 and in >1 sample.
pm.counts_z_F <- pm.counts_nz[, !i] # this is a matrix of where an ASV is only in quantity "1" for only 1 sample
dim(pm.counts_nz_F)
dim(pm.counts_z_F) # all zero features

# Sanity check for if any zero sum columns remain
dd1 = pm.counts_z_F[,colSums(pm.counts_z_F) > 1]
dd1
dd2 = pm.counts_nz_F[,colSums(pm.counts_nz_F) == 1] # zero singletons remain
dd2
# No col sums sum to 1 or 0.
rm(list=ls()[! ls() %in% c("pm.counts_nz_F", "seed")])
dim(pm.counts_nz_F) # We now have a 1218 x 119097 matrix for our downstream analysis

saveRDS(pm.counts_nz_F, file="~/pm_workflow/.../output/pm_counts1218.rds") # Main count table to be used from here on out


# Remove correlated variables (crude feature selection/dimensionality reduction)

dim(pm.counts_nz_F)
pm.counts_nz_F.del.column <- pm.counts_nz_F[, -c(35001:119097)]
# Find attributes that are highly correlated use cutoff 0.7 or .75 (standard)
set.seed(seed)
cutoff <- 0.70
correlations <- cor(pm.counts_nz_F.del.column[,1:35000])
library(caret)
highlyCorrelated <- findCorrelation(correlations, cutoff=cutoff)
for (value in highlyCorrelated) {
  print(names(pm.counts_nz_F.del.column)[value])
}
dim(pm.counts_nz_F)
dim(pm.counts_nz_F.del.column)

# Create a new dataset without highly correlated features
ASV.counts.noncorvars <- pm.counts_nz_F.del.column[,-highlyCorrelated]
dim(pm.counts_nz_F.del.column)
dim(ASV.counts.noncorvars) # 6916
ASV.counts.noncorvars[,1218]
# Column headers are actual DNA sequences

saveRDS(ASV.counts.noncorvars, file="~/pm_workflow/.../output/pm_counts_noncorvars.rds") # Non correlated vars count table


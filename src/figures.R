# ============================== Main Figures =============================================
# Script contents:
#
# Preprocess data, select amplicon sequence variant (ASV) features that will be used 
# throughout all downstream analysis
#
# Statistics of read counts of dominant taxa: Phyla/Class
# All figures
# =========================================================================================

# ==== Data preprocessing and statistics of dominant taxa ====
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

#length after prevalence filtering
length(get_taxa_unique(pm.ps3, taxonomic.rank = "Phylum")) # 39
length(get_taxa_unique(pm.ps4, taxonomic.rank = "Phylum")) # 24

length(get_taxa_unique(pm.ps3, taxonomic.rank = "Class")) # 68
length(get_taxa_unique(pm.ps4, taxonomic.rank = "Class")) # 38

length(get_taxa_unique(pm.ps3, taxonomic.rank = "Order")) # 190
length(get_taxa_unique(pm.ps4, taxonomic.rank = "Order")) # 114

length(get_taxa_unique(pm.ps3, taxonomic.rank = "Family")) # 468
length(get_taxa_unique(pm.ps4, taxonomic.rank = "Family")) # 223

length(get_taxa_unique(pm.ps3, taxonomic.rank = "Genus")) # 1579
length(get_taxa_unique(pm.ps4, taxonomic.rank = "Genus")) # 485

# Run some statistics of the dominant phyla/classes for the port microbiome biogeography paper

# Compute ASVs for each phyla/class 
ASV.phyla<-as.data.frame(table(tax_table(pm.ps4)[, "Phylum"], exclude = NULL))
colSums(ASV.phyla[2]) # here are the ASV sums per phyla.
dim(ASV.phyla)

ASV.class<-as.data.frame(table(tax_table(pm.ps4)[, "Class"], exclude = NULL))
colSums(ASV.class[2]) # here are the ASV sums per phyla.
dim(ASV.class)

# Find the % total of ASVs from the phyla/class of interest comprise
for(i in 1:length(ASV.class)){
  perctotal.class<-as.data.frame(ASV.class$Freq/3214*100)
}

for(i in 1:length(ASV.phyla)){
  perctotal.phylum<-as.data.frame(ASV.phyla$Freq/3214*100)
}

summarize_phyloseq(pm.ps4)

# Agglomerate phyloseq object by phylum and class
pm.ps.phylum = tax_glom(pm.ps4, "Phylum", NArm = TRUE, bad_empty=c(NA, "", " ", "\t"))
pm.ps.phylum # 24

pm.ps.class = tax_glom(pm.ps4, "Class", NArm = TRUE, bad_empty=c(NA, "", " ", "\t"))
pm.ps.class # 38

# Generate count matrices to sum abundance of top phyla and top classes.

# Phylum
pm.ps.phylum <- psmelt(pm.ps.phylum)
pm.ps.phylum$Phylum <- as.character(pm.ps.phylum$Phylum)
pm.ps.phylum
pm.ps.phylum <- aggregate(Abundance~Sample+Phylum, pm.ps.phylum, FUN=sum)
pm.ps.phylum.mat <- cast(pm.ps.phylum, Sample ~ Phylum)
dim(pm.ps.phylum.mat)
head(pm.ps.phylum.mat)

# Class
pm.ps.class <- psmelt(pm.ps.class)
pm.ps.class$Class <- as.character(pm.ps.class$Class)
pm.ps.class
pm.ps.class <- aggregate(Abundance~Sample+Class, pm.ps.class, FUN=sum)
pm.ps.class.mat <- cast(pm.ps.class, Sample ~ Class)
dim(pm.ps.class.mat)
head(pm.ps.class.mat)


colSums(pm.ps.phylum[3]) # 14085539, this is the total abundance of all taxa through all 1218 samples
colSums(pm.ps.class[3]) # 14085539, sanity check

# Identify sums of abundance of each phylum and class, and % total reads
readstats<-function(glommed.mat, ...){
  glommed.mat[1]<-NULL
  df.glommed.mat<-as.data.frame(colSums(glommed.mat[,]))
  glommed.col.name<-colnames(glommed.mat)
  df.glommed.mat<-cbind(df.glommed.mat, glommed.col.name)
  names(df.glommed.mat)[1]<-"summed.abundance"
  names(df.glommed.mat)[2]<-"taxa"
  colSums(df.glommed.mat[1])
  for(i in 1:length(df.glommed.mat)){
    freq.taxa<-as.data.frame(df.glommed.mat$summed.abundance/14085539*100)
    names(freq.taxa)[1]<-"percentreads"
  }
  df.glommed.mat<-cbind(df.glommed.mat, freq.taxa)
}

phylum.stats<-readstats(pm.ps.phylum.mat) # Returns a df of stats per phyla
class.stats<-readstats(pm.ps.class.mat) # Per class

# Identify the dominant phya and class based on percentage of reads, those within >=10% of total reads
# This returns the summed abundance (total reads), taxa and percent reads of total reads per taxa
dom.phyla<-phylum.stats[phylum.stats$percentreads >= 10, ]

# Identify those dominant classes that are >=5% of total reads
dom.class<-class.stats[class.stats$percentreads >= 5, ]
dom.phyla
dom.class


colSums(ASV.phyla[2]) ## 3214, this is asv count.
ASV.phyla # Find how many ASVs belong to each phyla
#2       Actinobacteria  243
#4        Bacteroidetes  944
#6        Cyanobacteria  148
#20      Proteobacteria 1384


colSums(ASV.class[2])
ASV.class # Find how many ASVs belong to each class

# Ensure that these classes comprise >=40% of their respective phyla's ASV content
#phyla: Actinobacteria
#1       Acidimicrobiia  101
101/243*100 # [1] 41.56379%
#3       Actinobacteria  129
129/243*100 # [1] 53.08642%

#phyla: Bacteroidetes
#9          Bacteroidia  927
927/944*100 # [1] 98.19915%

#phyla: Cyanobacteria
#31    Oxyphotobacteria  148
148/148*100 # [1] 100%

#phyla: Proteobacteria 1384
#4  Alphaproteobacteria  630
630/1384*100 # [1] 45.52023
#19 Gammaproteobacteria  671
671/1384*100 # [1] 48.48266

# % total asvs
ASV.phyla
phyla.x<-ASV.phyla$Var1
class.x<-ASV.class$Var1
perctotal.phylum<-cbind(perctotal.phylum, phyla.x) # Gives % total ASVs
perctotal.class<-cbind(perctotal.class, class.x)

# Here are the statistics for % total reads for both Phyla and Class, Table S1
phylum.stats
class.stats


rm(list=ls()[!ls() %in% c("seed", "pm.ps4",
                          "pm.taxa", "pm.metadata", "pm.counts")])

detachAllPackages <- function() {
  
  basic.packages <- c("package:stats","package:graphics","package:grDevices","package:utils","package:datasets","package:methods","package:base")
  
  package.list <- search()[ifelse(unlist(gregexpr("package:",search()))==1,TRUE,FALSE)]
  
  package.list <- setdiff(package.list,basic.packages)
  
  if (length(package.list)>0)  for (package in package.list) detach(package, character.only=TRUE)
  
}
detachAllPackages()

library("rlang", lib.loc="~/R/x86_64-redhat-linux-gnu-library/3.5")
library("ggplot2", lib.loc="~/R/x86_64-redhat-linux-gnu-library/3.5")
library("Rcpp", lib.loc="~/R/x86_64-redhat-linux-gnu-library/3.5")
library("ggplot2", lib.loc="~/R/x86_64-redhat-linux-gnu-library/3.5")
library("vegan")
library("csv")
library("plotly")
library("plyr")
library("dplyr")
library("phyloseq")
library("reshape")
library("reshape2")
library("plotly")

seed=81
set.seed(seed)



# ==== Figure 1 ====
pkgs <- c("ggplot2", "ggmap", "maps", "mapdata") #load packages
invisible(sapply(pkgs,require, character = TRUE))

world <- map_data("world")
ggplot() + geom_polygon(data = world, aes(x=long, y = lat, group = group)) + 
  coord_fixed(1.3)

gg1 <- ggplot() + 
  geom_polygon(data = world, aes(x=long, y = lat, group = group), fill = "#EDEDED", color = "darkgray") + 
  coord_fixed(1.3)
gg1

# These labels are coordinates latitude/longitude of the sampling location ports from this study
labs <- data.frame(
  long = c(114.171654, 103.822441, 129.106323, -89.9617, -76.27636, 
           -74.111938, -79.897705, -94.827751, -76.33596, 
           14.249008, 14.23758, 12.369708, 4.54956, 8.1069, 
           -88.650375, -92.180984, -87.944046, -118.276642, 
           -122.500671, -122.405815),
  
  lat = c(22.29073, 1.260632, 35.076935, 29.922594, 39.124435, 40.646053,
          32.857979, 29.411531, 36.981754, 40.823879, 40.79857,
          45.423679, 51.903225, 53.5323, 47.239178, 46.714516, 
          44.557785, 33.713993, 37.812164, 47.630363),
  names = c("SWFSC-FED", "NWFSC"),
  stringsAsFactors = FALSE
)  

gg1 + 
  geom_point(data = labs, aes(x = long, y = lat), color = "#000000", size = 8) +
  geom_point(data = labs, aes(x = long, y = lat), color = "#FCDF8A", size = 7) + theme_bw()



# ==== Figure 2 ====
rm(list=ls()[!ls() %in% c("seed", "pm.ps4",
                          "pm.taxa", "pm.metadata", "pm.counts")])
pm.ps4
pkgs <- c("phyloseq", "ggplot2", "rlang", "Rcpp", "vegan", "csv", "plotly", "plyr", "dplyr", "reshape", "reshape2") # Load packages
invisible(sapply(pkgs,require, character = TRUE))

# Merge phyloseq object by geo_port (20 locations/levels) variable from metadata
colnames(pm.ps4@sam_data)
merged.pm.ps.location = merge_samples(pm.ps4, "geo_port")
merged.pm.ps.location
otu_table(merged.pm.ps.location)[1:5, 1:5]

# Create a phyloseq object agglomerated by bacterial class
pm.ps.class.by.location = tax_glom(merged.pm.ps.location, "Class", NArm = TRUE, bad_empty=c(NA, "", " ", "\t"))
pm.ps.class.by.location

# Create a count matrix of class features by each sampling location
library(reshape)
pm.ps.class.melt <- psmelt(pm.ps.class.by.location)
pm.ps.class.melt$Class <- as.character(pm.ps.class.melt$Class)
pm.ps.class.melt
pm.ps.class.melt <- aggregate(Abundance~Sample+Class, pm.ps.class.melt, FUN=sum)
pm.ClassMatrix.loc <- cast(pm.ps.class.melt, Sample ~ Class)
dim(pm.ClassMatrix.loc)
head(pm.ClassMatrix.loc)
pm.ClassMatrix.loc

# Remove all bacterial classes but the dominant six used in this study, and agglomerate all "other" into "other" so we can display as a relative abundance taxa plot
pm.ClassMatrix.loc2<-pm.ClassMatrix.loc
keepclasscols<-c("Sample", "Acidimicrobiia", "Actinobacteria", 
                 "Alphaproteobacteria", "Bacteroidia", 
                 "Gammaproteobacteria", "Oxyphotobacteria")
pm.ClassMatrix.loc2<-pm.ClassMatrix.loc2[,keepclasscols]
dim(pm.ClassMatrix.loc2)

# Convert all other taxa into other variable.
# Remove the top six to sum rows without them.
pm.ClassMatrix.loc.notop6<-pm.ClassMatrix.loc[, ! names(pm.ClassMatrix.loc) %in% keepclasscols, drop = F]

# Sum all of the rows so that way they can collapse into "Other" category across all 20 locations.
Other<-rowSums(pm.ClassMatrix.loc.notop6)
cbindall<-cbind(pm.ClassMatrix.loc2, Other)
cbindall$Sample

row.names(cbindall) <- cbindall$Sample
row.names(cbindall)
cbindall[1] <- NULL

# Convert counts to relative abundance
cbindall.RA=sweep(cbindall,1,rowSums(cbindall),"/")
rowSums(cbindall.RA)

relabund.dom.classes<-as.matrix(cbindall.RA)
relabund.dom.classes<-as.data.frame(relabund.dom.classes)

# Reimport a formatted data that can be used for a taxa plot
plot.data<-'https://github.com/rghannam/pm_workflow/raw/master/data/graphing_inputs/fig2/abundance.plot.formatted.csv'
plot.data<-read.csv(url(plot.data))

class(plot.data$stack)
class(plot.data$xaxis)
class(plot.data$value)
plot.data$stack <- factor(plot.data$stack, levels = c('Acidimicrobiia', 'Actinobacteria', 'Alphaproteobacteria', 'Bacteroidia', 'Gammaproteobacteria', 'Oxyphotobacteria', 'Other'))
# Order as Figure 1 is: East West Europe Asia Lakes # order by region
plot.data$xaxis <- factor(plot.data$xaxis, levels = c('Keweenaw', 'Green Bay', 'Duluth', 'Singapore', 
                                                      'Hong Kong', 'Busan', 'Wilhelmshaven', 'Venice', 
                                                      'Rotterdam', 'Naples', 'Martigues', 'Seattle', 
                                                      'Oakland', 'Los Angeles', 'Norfolk', 'New York', 
                                                      'New Orleans', 'Galveston', 'Charleston', 'Baltimore'))


# Define plotting parameters
figyaxis <- list(title = "Port",
                 showline = TRUE,
                 showgrid = TRUE,
                 showticklabels = TRUE,
                 linecolor = 'black',
                 linewidth = 1,
                 tickformat = "%",
                 autotick = TRUE,
                 ticks = 'outside',
                 tickcolor = 'black',
                 tickwidth = 2,
                 ticklen = 5,
                 tickfont = list(family = 'Arial',
                                 size = 14,
                                 color = '#000000'))

figxaxis <- list(title = "Relative Abundance (% all reads)",
                 showline = TRUE,
                 showgrid = TRUE,
                 showticklabels = TRUE,
                 linecolor = 'black',
                 linewidth = 1,
                 autotick = TRUE,
                 tickformat = "%",
                 ticks = 'outside',
                 tickcolor = 'black',
                 tickwidth = 2,
                 ticklen = 5,
                 tickfont = list(family = 'Arial',
                                 size = 14,
                                 color = '#000000'))

marker_style <- list(line = list(width = 0.5,
                                 color = 'rgb(0, 0, 0)'));
t <- list(
  family = "Arial",
  size = 12,
  color = 'black')


# Establish colorset
library(ggsci)
pal_locuszoom()(7)
pal_ucscgb()(7)
pal_startrek()(7)

colors1<-c("#CC0C00FF" ,"#FF7F0EFF", "#FFCD00FF" ,"#00B5E2FF" ,"#00AF66FF" ,"#8755B5FF" ,"#FFFFFF")
plot.Fig2 <-plot.data %>% group_by(stack) %>% arrange(xaxis) %>%
  plot_ly( x = ~value, y = ~xaxis, color= ~stack, type = 'bar', orientation = 'h', mode='marker',
           marker=marker_style, textposition = 'auto',
           colors=colors1,
           #colorscale='Viridis',
           line = list(color = "#000000", width = 1)) %>%
  layout(title = "Figure 2", barmode = 'stack',
         xaxis = figxaxis, font=t,
         yaxis = figyaxis)

plot.Fig2



# ==== Figure 3 ====
rm(list=ls()[!ls() %in% c("seed", "pm.ps4",
                          "pm.taxa", "pm.metadata", "pm.counts")])

pkgs <- c("ggplot2", "rlang", "Rcpp", "vegan", "csv", "plotly", "plyr", "dplyr", "reshape", "reshape2", "tidyverse") # Load packages
invisible(sapply(pkgs,require, character = TRUE))
pm.ps4
length(get_taxa_unique(pm.ps4, taxonomic.rank = "Class")) 

# Agglomerate by class
pm.ps.class.1218x38 = tax_glom(pm.ps4, "Class", NArm = TRUE, bad_empty=c(NA, "", " ", "\t"))
pm.ps.class.1218x38

# Melt to dataframe
pm.ps.class.melt <- psmelt(pm.ps.class.1218x38)
pm.ps.class.melt$Class <- as.character(pm.ps.class.melt$Class)
pm.ps.class.melt
pm.ps.class.melt <- aggregate(Abundance~Sample+Class, pm.ps.class.melt, FUN=sum)
pm.ClassMatrix1218 <- cast(pm.ps.class.melt, Sample ~ Class)
dim(pm.ClassMatrix1218)

# Only keep those top 6 dominant
keepclasscols<-c("Sample", "Acidimicrobiia", "Actinobacteria", "Alphaproteobacteria", "Bacteroidia", "Gammaproteobacteria", "Oxyphotobacteria")
pm.ClassMatrix1218<-pm.ClassMatrix1218[,keepclasscols]
dim(pm.ClassMatrix1218)

pm.ClassMatrix1218$Sample
row.names(pm.ClassMatrix1218) <- pm.ClassMatrix1218$Sample
row.names(pm.ClassMatrix1218)
pm.ClassMatrix1218[1] <- NULL # Remove Sample column

# Convert to relative abundance
pm.ClassMatrix1218.RA=sweep(pm.ClassMatrix1218,1,rowSums(pm.ClassMatrix1218),"/")
rowSums(pm.ClassMatrix1218.RA)

pm.ClassMatrix1218.RA # This is used to format plotting data


plot.data2<-'https://github.com/rghannam/pm_workflow/raw/master/data/graphing_inputs/fig3/boxplot_formatted.csv'
plot.data2<-read.csv(url(plot.data2))

# Define plot parameters
marker_style <- list(line = list(width = 0.5,
                                 color = 'rgb(0, 0, 0)'));

figyaxis2 <- list(title = "Relative Abundance",
                  showline = TRUE,
                  showgrid = TRUE,
                  showticklabels = TRUE,
                  linecolor = 'black',
                  linewidth = 1,
                  autotick = TRUE,
                  tick0=0,
                  dtick=5,
                  ticks = 'outside',
                  tickcolor = 'black',
                  tickwidth = 2,
                  ticklen = 5,
                  tickfont = list(family = 'Arial',
                                  size = 12,
                                  color = '#000000'))

figxaxis2 <- list(title = "Region",
                  showline = TRUE,
                  showgrid = TRUE,
                  showticklabels = TRUE,
                  linecolor = 'black',
                  linewidth = 1,
                  autotick = TRUE,
                  tick0=0,
                  dtick=5,
                  ticks = 'outside',
                  tickcolor = 'black',
                  tickwidth = 2,
                  ticklen = 5,
                  tickfont = list(family = 'Arial',
                                  size = 16,
                                  color = '#000000'))

plot.data2$group <- factor(plot.data2$group, levels = c('Acidimicrobiia', 'Actinobacteria', 'Alphaproteobacteria', 'Bacteroidia', 'Gammaproteobacteria', 'Oxyphotobacteria'))

plot.data2$region <- factor(plot.data2$region, levels = c('East', 'West', 'Europe', 'Asia', 'Great Lakes'))

Fig3.jitters <- plot_ly(plot.data2, x = ~region, y = ~abundance, color = ~group,marker=marker_style,type = "box",
                        colors=c("#CC0C00FF" ,"#FF7F0EFF", "#FFCD00FF" ,"#00B5E2FF" ,"#00AF66FF" ,"#8755B5FF")) %>% 
  layout(boxmode = "group", title = 'X', showlegend=T, legend = list(orientation = 'h'),
         xaxis = figxaxis2,
         yaxis = figyaxis2)

Fig3.jitters # With jitters 

# Remove jitters
plot.Fig3 <- plotly_build(Fig3.jitters) 
plot.Fig3$x$data <- lapply(plot.Fig3$x$data, FUN = function(x){
  x$marker = list(opacity = 0)
  return(x)
})

plot.Fig3

# Analysis
# Find the mean of of each bacterial class at each region

# Create average abundance vector as an extra filtering parameter to be added to the final dataframe
colnames(plot.data2)
percentage <- prop.table(table(plot.data2$region)) * 100
percentage
getfreq<-cbind(freq=table(plot.data2$region), percentage=percentage)
getfreq<-as.data.frame(getfreq)
getfreq
row.names(getfreq) <- c("East", "West", "Europe", "Asia", "Great Lakes")
freq <- c("freq")
getfreq<-getfreq[,freq]
getfreq<-as.data.frame(getfreq)
getfreq.labels.vec<-dplyr::pull(getfreq, getfreq)
getfreq.labels.vec

sum.region<-ddply(plot.data2, "region", numcolwise(sum))
sum.region # Find the total samples (n) for each region (1218 total)

# Find total abundance
boxplotX2<-plot.data2
class(boxplotX2$abundance)
boxplotX3<-boxplotX2 %>% 
  group_by(region) %>% 
  do({
    sum_value = sum(distinct(., region, abundance)$abundance);
    mutate(., sum_value = sum_value)
  })
head(boxplotX3)
require(tibble)
tbl_df(boxplotX3) # Make it into a tbl class

boxplotX4<-boxplotX3 %>%
  mutate(avgabund=abundance/sum_value)
head(boxplotX4) # This has a working dataframe of prevalence/abundance.summed/prevalence/avg abundance.
boxplotX4

# This sums the classes by region (abundregion) to divide total samples per region
require(data.table) 
DT <- data.table(boxplotX4) 
head(DT)
DT<-DT[ , .(abundregion = sum(abundance)), by = .(group,region)] # This gives each region per class total abundance
DT

# Then divide these values by vec: sum.region
sum.region
require(tidyverse)
avgabund="avgabund"
DT[ , "avgabund"] <- avgabund
sum.region
DT

# Just to divide these totals by samples (n) can get an average relative abundance to compare to plot.Fig3
region_split <- split(DT, DT$region)
reg_names <- c("East", "West", "Europe", "Asia", "Great Lakes")
for (i in 1:length(region_split)) {
  assign(reg_names[i], region_split[[i]])
  East<-region_split[[1]]
  West<-region_split[[2]]
  Europe<-region_split[[3]]
  Asia<-region_split[[4]]
  Lakes<-region_split[[5]]
  East$avgabund<-region_split[[1]]$abundregion/355
  West$avgabund<-region_split[[2]]$abundregion/182
  Europe$avgabund<-region_split[[3]]$abundregion/294
  Asia$avgabund<-region_split[[4]]$abundregion/191
  Lakes$avgabund<-region_split[[5]]$abundregion/196
  region_split<-c(list(East, West, Europe, Asia, Lakes))
}
region_split # This gives average relative abund of each class at each region

# Sanity check
# Sum samples of each region 
# Sum the factor based on unique levels
library(dplyr)
ddply(DT,~region,summarise,number_of_distinct_orders=length(unique(region)))



# ==== Figure 4 ====
# These models come from machine_learning.R
pkgs <- c("tidyverse", "phyloseq", "caretEnsemble", "caret", "randomForest", "plyr", "tidyverse", "data.table", "csv", "magrittr") # Load packages
invisible(sapply(pkgs,require, character = TRUE))

# Load ASV local(port) and regional models
asv.local<-'https://github.com/rghannam/pm_workflow/raw/master/data/machine_learning/models/asv.local.rds'
asv.region<-'https://github.com/rghannam/pm_workflow/raw/master/data/machine_learning/models/asv.region.rds'
asv.local<-readRDS(url(asv.local))
asv.region<-readRDS(url(asv.region))

# Review variable importance from caret
getFromNamespace('varImp','caret')
function (object, ...) 
{
  UseMethod("varImp")
}

getS3method('varImp','randomForest')
function (object, ...) 
{
  code <- varImpDependencies("rf")
  code$varImp(object, ...)
}

code <- caret:::varImpDependencies('rf')
code$varImp

# vi = important variable
regvi<-varImp(asv.region$rf)
regvi
regvi[["importance"]]
regvidf<-as.data.frame(regvi[["importance"]])  # Coerce overall vi's to df, arranged in order

portvi<-caret::varImp(asv.local$rf)
portvi
portvidf<-as.data.frame(portvi[["importance"]])
head(portvidf)

# Replace the important variables under the speciesID (ASV ID) with what actual taxa they correspond to
portasvs<-rownames(portvidf)
portasvs<-str_remove(portasvs, "sp")
regionasvs<-rownames(regvidf)
regionasvs<-str_remove(regionasvs, "sp")

# Take column and convert it to a row in df.
regvidf$asvordered<-regionasvs # Appended an ASV ordered column so we can compare the two columns: asvordered/Overall
portvidf$asvordered<-portasvs

# Subset by overall variable importance >=1
regvidf.ss1 <- regvidf[regvidf$Overall >= "1", ]
portvidf.ss1 <- portvidf[portvidf$Overall >= "1", ]

# Find differences between port/regional models.
sd.reg<-as.data.frame(setdiff(regvidf.ss1$asvordered,portvidf.ss1$asvordered)) # Has the vars that are in region but not in port
sd.port<-as.data.frame(setdiff(portvidf.ss1$asvordered,regvidf.ss1$asvordered)) # Has the vars that are in port but not in region

regvidf.ss1$region <- TRUE
portvidf.ss1$port <- TRUE
portregionmerge <- merge(regvidf.ss1, portvidf.ss1, all=TRUE) # Just ordered the important vars in one df

library(compare)
library(tidyverse)
regvidf.ss1 %>% remove_rownames() %>% column_to_rownames(var = "asvordered")
portvidf.ss1 %>% remove_rownames() %>% column_to_rownames(var = "asvordered")

common <- intersect(regvidf.ss1$asvordered, portvidf.ss1$asvordered) # Find shared ASVs between port/regional models 
common<-as.data.frame(common)

library(plyr)

sd1<-as.data.frame(setdiff(regvidf.ss1$asvordered,portvidf.ss1$asvordered)) # Has the ones that are in region but not in port
sd2<-as.data.frame(setdiff(portvidf.ss1$asvordered,regvidf.ss1$asvordered)) # Has the ones that are in port but not in region

require(phyloseq)
# Try collapsing all ASVS into respective classes, e.g. assign a spX value (from pm.ps4@tax_table) to each bacterial Class, respectively
# So that the ASV can match with its variable importance
tax_table(pm.ps4)[1:20, 1:5]

# Extract taxonomy table
taxadf<-as.data.frame(tax_table(pm.ps4))
head(taxadf)
taxadt<-as.data.table(taxadf)
head(taxadt)
taxadf[,4:6]<-NULL
head(taxadf)
taxadt[,4:6]<-NULL
taxadf # Assign rownames to column
setDT(taxadf, keep.rownames = TRUE)[]
head(taxadf)

library(data.table)
# Collapses all classes and shows which variable of ASV exactly belongs to this class
# Can use this in linking the shared/different variables from our object: common or sd.port/sd.reg
taxadfcollapse<-setDT(taxadf)[, lapply(.SD, paste, collapse = "; "), by = Class]
head(taxadfcollapse)
taxadfcollapse[,3:4]<-NULL

head(taxadfcollapse)

# This shows the distribution of ASVs collapses by Class in 3214 ps objet, extract those that are important vars for the models
# e.g. 630/3214 ASVs belong to Alphaproteobacteria
library(magrittr)
res = taxadf %>% 
  melt(id.vars = "Class") %>% 
  dcast(Class ~ variable)

head(res)
colSums(res[,-1]) # Sanity check to see if sums to 3214

# This also helps explain why the model chose certain vars as imp, as the majority of ASVs are within this
# category of classes, they're the most used just by nature of the model.

# Numbers are formatted as strings but can convert

# This means: Acidimicrobiia 101/3214*100 accounted for X % of total ASV features.
# And so, an agglomeration is a feature reduction method that accumulates all of those ASVs into 1 variable, or composites them
# However - we lose accuracy dependent on the taxonomic resolution we choose to look at
head(res) # rn is the total asvs for that class.

# Take shared as a vector - and extract class by sp.
names(common)[1]<-"asvs"
spvector<-as.vector(common$asvs)
spvector<-as.numeric(spvector)
head(taxadf)

taxadf2<-taxadf
taxadf2$rn<-gsub('sp', '', taxadf2$rn)
class(taxadf2$rn)

extractedsharedasvs <- taxadf2[taxadf2$rn %in% spvector, ]
extractedsharedasvs[,1:3]<-NULL

# Append vector as column
extractedsharedasvs$ASV = spvector

length(unique(extractedsharedasvs$Class)) # There are only 8 unique classes of the 68 here

astibble<-as_tibble(extractedsharedasvs)
astibble[!duplicated(astibble$Class), ] # Gives the 8 unique classes that are within this dataframe

# What this shows is that of the 68 ASVS that are shared between the models that have varImp >=1
# Only 8 classes comprise the 68 ASVs - of these 8, the top dominant six from our analysis are inside of this

extractedsharedasvs$ASV
regvidf.ss1$asvordered
df.r<-regvidf.ss1[regvidf.ss1$asvordered %in% spvector,] # Extracted data for overall var IMP
portvidf.ss1$asvordered
df.p<-portvidf.ss1[portvidf.ss1$asvordered %in% spvector,]
# Based on the common vector - we now have the overall importance for region.

# Now - append to one dataframe
# Rename column so we know if its port of region
colnames(df.r)[1]<-"vimpRegion"
colnames(df.p)[1]<-"vimpPort"

# Extract column from extractedsharedasvs and so as to append to the asvordered column to this dataframe
Classordered<-as.vector(extractedsharedasvs$Class)
vimpRegion<-as.vector(df.r$vimpRegion)
vimpPort<-as.vector(df.p$vimpPort)
extractedsharedasvs$vimpRegion<-vimpRegion
extractedsharedasvs$vimpPort<-vimpPort

# Now we have a final dataframe with:
# Shared ASVs between model varimp>=1%
# Overall important variables from each model
# Which class the ASV belonged to

#detach("package:plyr", unload=TRUE) 
esasvsgpclass<-as.data.frame(extractedsharedasvs %>% group_by(Class) %>% summarize(count=n()))
esasvsgpclass
colSums(esasvsgpclass[2])
# This groups the unique entries and shows proportion of classes per total ASVs
# e.g. out of 68 shared ASVs, 62 were in the top 6 dominant classes, 91.17% of important asvs are within the dominant taxa

extractedsharedasvs # Has distribution of 68 shared ASVs between each models along with vImp region vImp port
esasvsgpclass # Shows the distribution of the 68 - bins 68 into respective ASV
portregionmerge # Has the merged vimp >1
dim(portvidf.ss1) # 250 vimp > 1
dim(regvidf.ss1) # 92 vimp > 1
res # Distribution of 3214 asvs

# Extract values from extractedsharedasvs for plotting important predictor variables.
# Order classes together and then apply a classcategory factor.
vImp.df<-extractedsharedasvs
vImp.df<-with(vImp.df, vImp.df[order(Class),]) # Takes multiple arguments, we just need to arrange by class
class(vImp.df)
vImp.df<-as.data.frame(vImp.df)
vImp.df
nclasses<-as.data.frame(vImp.df %>% group_by(Class) %>% summarize(count=n())) # Get the number of classes
nclasses

# There should be 8 that are outside of the 6 dominant, class these as "Other".
levels(vImp.df$Class)[levels(vImp.df$Class)=="Planctomycetacia"] <- "Other"
levels(vImp.df$Class)[levels(vImp.df$Class)=="Verrucomicrobiae"] <- "Other"
levels(vImp.df$Class)

# Plot these important predictors
# Define plotting parameters
library(plotly)
marker_style <- list(line = list(width = 0.5,
                                 color = 'rgb(0, 0, 0)'));
figyaxis <- list(title = "Variable Importance (Overall)",
                 showline = TRUE,
                 showgrid = TRUE,
                 showticklabels = TRUE,
                 linecolor = 'black',
                 linewidth = 1,
                 autotick = TRUE,
                 tick0=0,
                 dtick=5,
                 ticks = 'outside',
                 tickcolor = 'black',
                 tickwidth = 2,
                 ticklen = 5,
                 tickfont = list(family = 'Arial',
                                 size = 14,
                                 color = 'rgb(82, 82, 82)'))

figxaxis <- list(title = "Class",
                 showline = TRUE,
                 showgrid = TRUE,
                 showticklabels = TRUE,
                 linecolor = 'black',
                 linewidth = 1,
                 autotick = TRUE,
                 tick0=0,
                 dtick=5,
                 ticks = 'outside',
                 tickcolor = 'black',
                 tickwidth = 2,
                 ticklen = 5,
                 tickfont = list(family = 'Arial',
                                 size = 14,
                                 color = 'rgb(82, 82, 82)'))


plot.Fig4B.Jitters <- vImp.df %>%
  plot_ly(type = 'violin') %>%
  add_trace(
    x = ~Class[vImp.df$Class == 'Acidimicrobiia'],
    y = ~vimpPort[vImp.df$Class == 'Acidimicrobiia'],
    legendgroup = 'Port',
    scalegroup = 'Port',
    name = 'Port',
    box = list(
      visible = F
    ),
    meanline = list(
      visible = F
    ),
    line = list(
      color = 'blue'
    )
  ) %>%
  add_trace(
    x = ~Class[vImp.df$Class == 'Acidimicrobiia'],
    y = ~vimpRegion[vImp.df$Class == 'Acidimicrobiia'],
    legendgroup = 'Reg',
    scalegroup = 'Reg',
    name = 'Reg',
    box = list(
      visible = F
    ),
    meanline = list(
      visible = F
    ),
    line = list(
      color = 'pink'
    )
  ) %>% 
  
  #add class
  add_trace(
    x = ~Class[vImp.df$Class == 'Actinobacteria'],
    y = ~vimpPort[vImp.df$Class == 'Actinobacteria'],
    legendgroup = 'Port',
    scalegroup = 'Port',
    name = 'Port',
    box = list(
      visible = F
    ),
    meanline = list(
      visible = F
    ),
    line = list(
      color = 'blue'
    )
  ) %>%
  add_trace(
    x = ~Class[vImp.df$Class == 'Actinobacteria'],
    y = ~vimpRegion[vImp.df$Class == 'Actinobacteria'],
    legendgroup = 'Reg',
    scalegroup = 'Reg',
    name = 'Reg',
    box = list(
      visible = F
    ),
    meanline = list(
      visible = F
    ),
    line = list(
      color = 'pink'
    )
  ) %>% 
  
  #end add class
  
  #add class
  add_trace(
    x = ~Class[vImp.df$Class == 'Alphaproteobacteria'],
    y = ~vimpPort[vImp.df$Class == 'Alphaproteobacteria'],
    legendgroup = 'Port',
    scalegroup = 'Port',
    name = 'Port',
    box = list(
      visible = F
    ),
    meanline = list(
      visible = F
    ),
    line = list(
      color = 'blue'
    )
  ) %>%
  add_trace(
    x = ~Class[vImp.df$Class == 'Alphaproteobacteria'],
    y = ~vimpRegion[vImp.df$Class == 'Alphaproteobacteria'],
    legendgroup = 'Reg',
    scalegroup = 'Reg',
    name = 'Reg',
    box = list(
      visible = F
    ),
    meanline = list(
      visible = F
    ),
    line = list(
      color = 'pink'
    )
  ) %>% 
  
  #end add class
  
  #add class
  add_trace(
    x = ~Class[vImp.df$Class == 'Bacteroidia'],
    y = ~vimpPort[vImp.df$Class == 'Bacteroidia'],
    legendgroup = 'Port',
    scalegroup = 'Port',
    name = 'Port',
    box = list(
      visible = F
    ),
    meanline = list(
      visible = F
    ),
    line = list(
      color = 'blue'
    )
  ) %>%
  add_trace(
    x = ~Class[vImp.df$Class == 'Bacteroidia'],
    y = ~vimpRegion[vImp.df$Class == 'Bacteroidia'],
    legendgroup = 'Reg',
    scalegroup = 'Reg',
    name = 'Reg',
    box = list(
      visible = F
    ),
    meanline = list(
      visible = F
    ),
    line = list(
      color = 'pink'
    )
  ) %>% 
  
  #end add class
  
  #add class
  add_trace(
    x = ~Class[vImp.df$Class == 'Gammaproteobacteria'],
    y = ~vimpPort[vImp.df$Class == 'Gammaproteobacteria'],
    legendgroup = 'Port',
    scalegroup = 'Port',
    name = 'Port',
    box = list(
      visible = F
    ),
    meanline = list(
      visible = F
    ),
    line = list(
      color = 'blue'
    )
  ) %>%
  add_trace(
    x = ~Class[vImp.df$Class == 'Gammaproteobacteria'],
    y = ~vimpRegion[vImp.df$Class == 'Gammaproteobacteria'],
    legendgroup = 'Reg',
    scalegroup = 'Reg',
    name = 'Reg',
    box = list(
      visible = F
    ),
    meanline = list(
      visible = F
    ),
    line = list(
      color = 'pink'
    )
  ) %>% 
  
  #end add class
  
  #add class
  add_trace(
    x = ~Class[vImp.df$Class == 'Oxyphotobacteria'],
    y = ~vimpPort[vImp.df$Class == 'Oxyphotobacteria'],
    legendgroup = 'Port',
    scalegroup = 'Port',
    name = 'Port',
    box = list(
      visible = F
    ),
    meanline = list(
      visible = F
    ),
    line = list(
      color = 'blue'
    )
  ) %>%
  add_trace(
    x = ~Class[vImp.df$Class == 'Oxyphotobacteria'],
    y = ~vimpRegion[vImp.df$Class == 'Oxyphotobacteria'],
    legendgroup = 'Reg',
    scalegroup = 'Reg',
    name = 'Reg',
    box = list(
      visible = F
    ),
    meanline = list(
      visible = F
    ),
    line = list(
      color = 'pink'
    )
  ) %>% 
  
  #end add class
  
  #add class
  add_trace(
    x = ~Class[vImp.df$Class == 'Other'],
    y = ~vimpPort[vImp.df$Class == 'Other'],
    legendgroup = 'Port',
    scalegroup = 'Port',
    name = 'Port',
    box = list(
      visible = F
    ),
    meanline = list(
      visible = F
    ),
    line = list(
      color = 'blue'
    )
  ) %>%
  add_trace(
    x = ~Class[vImp.df$Class == 'Other'],
    y = ~vimpRegion[vImp.df$Class == 'Other'],
    legendgroup = 'Reg',
    scalegroup = 'Reg',
    name = 'Reg',
    box = list(
      visible = F
    ),
    meanline = list(
      visible = F
    ),
    line = list(
      color = 'pink'
    )
  ) %>% 
  
  #end add class
  
  layout(boxmode = "group", title = 'X', showlegend=TRUE,
         xaxis = figxaxis,
         yaxis = figyaxis,
         violingap = 0,
         violingroupgap = 0,
         violinmode = 'overlay',
         legend = list(
           tracegroupgap = 0
         )
  )
plot.Fig4B.Jitters

# Remove points
plot.Fig4B <- plotly_build(plot.Fig4B.Jitters)
plot.Fig4B$x$data <- lapply(plot.Fig4B$x$data, FUN = function(x){
  x$marker = list(opacity = 0)
  return(x)
})

plot.Fig4B # Remove points, set opacity to markers=0.

# Make a pie chart of the ASV distribution now
# Get %'s of each class
nclasses
nclassvec<-as.vector(nclasses$count)
nclasses.perc <- (nclassvec/68)*100
nclasses.perc
category = c('Acidimicrobiia', 'Actinobacteria', 'Alphaproteobacteria', 'Bacteroidia',
             'Gammaproteobacteria', 'Oxyphotobacteria', 'Other') 
percentASV = c("4.41", "4.41", "17.64", "33.82", "26.47", "4.41", "8.82") # other is 6/68
legend = c('Acidimicrobiia', 'Actinobacteria', 'Alphaproteobacteria', 'Bacteroidia', 'Gammaproteobacteria', 'Oxyphotobacteria', 'Other')
piedf = data.frame(category, percentASV, legend)
piedf
piedf$category <- factor(piedf$category, levels = c('Acidimicrobiia', 'Actinobacteria', 'Alphaproteobacteria', 'Bacteroidia', 'Gammaproteobacteria', 'Oxyphotobacteria', 'Other'))
class(piedf$category)
piedf$category

# Legend attributes
l <- list(
  font = list(
    family = "sans-serif",
    size = 12,
    color = "#000000"),
  #bgcolor = "#E2E2E2",
  bordercolor = "#000000",
  borderwidth = 2,
  x= 1,
  y= 0.5)

plot.Fig4A <- plot_ly(piedf, labels = ~legend, values = ~percentASV, type = 'pie', sort=FALSE, # Orders it like df...don't have to use colors=~category, although not possible here in pitchart.
                      textposition = 'inside',
                      textinfo = 'percent',
                      hole=0.6,
                      insidetextfont = list(color = '#FFFFFF'),
                      hoverinfo = 'text',
                      text = ~paste('%', percentASV, ' percent'),
                      marker = list(colors = c("#CC0C00FF" ,"#FF7F0EFF", "#FFCD00FF" ,"#00B5E2FF" ,"#00AF66FF" ,"#8755B5FF" ,"#000000"),
                                    line = list(color = '#FFFFFF', width = 2)),
                      #The 'pull' attribute can also be used to create space between the sectors
                      showlegend = TRUE) %>%
  layout(title = 'Total ASV distribution', legend=l,autosize = T, width = 250, height = 250)
plot.Fig4A



# ==== Figure 5 ====

# Refer to differential_abundance.R for the plotting data was acquired
# Also stored under enrichment.factors.csv

# ==== Figure 6 ====
seed=81
set.seed(seed)
local.models<-'https://github.com/rghannam/pm_workflow/raw/master/data/machine_learning/models/taxarank.local.rds'
region.models<-'https://github.com/rghannam/pm_workflow/raw/master/data/machine_learning/models/taxarank.region.rds'
asv.local<-'https://github.com/rghannam/pm_workflow/raw/master/data/machine_learning/models/asv.local.rds'
asv.region<-'https://github.com/rghannam/pm_workflow/raw/master/data/machine_learning/models/asv.region.rds'
local.models<-readRDS(url(local.models))
region.models<-readRDS(url(region.models))
asv.local<-readRDS(url(asv.local))
asv.region<-readRDS(url(asv.region))

local.models
region.models
asv.local

pkgs <- c("caretEnsemble", "caret", "randomForest", "plyr", "tidyverse", "data.table", "csv", "magrittr") # Load packages
invisible(sapply(pkgs,require, character = TRUE))

# Create dataframe of model stats
local.models[[1]][["rf"]][["finalModel"]] # conf matrix of phylum local model
region.models[[1]][["rf"]][["finalModel"]] # conf matrix of phylum local model

# Store stats of each model in a list for local/regional
locstats<-local.models
for(i in 1:length(locstats)){
  locstats[[i]]<-list(local.models[[i]][["rf"]][["results"]])
}

regstats<-region.models
for(i in 1:length(regstats)){
  regstats[[i]]<-list(region.models[[i]][["rf"]][["results"]])
}


# Find mtry with the lowest logLoss
#TODO loop or function this
loc.min.phy<-locstats[["Phylum"]][[1]][which.min(locstats[["Phylum"]][[1]]$logLoss),]
loc.min.class<-locstats[["Class"]][[1]][which.min(locstats[["Class"]][[1]]$logLoss),]
loc.min.order<-locstats[["Order"]][[1]][which.min(locstats[["Order"]][[1]]$logLoss),]
loc.min.family<-locstats[["Family"]][[1]][which.min(locstats[["Family"]][[1]]$logLoss),]
loc.min.genus<-locstats[["Genus"]][[1]][which.min(locstats[["Genus"]][[1]]$logLoss),]

# Regional
reg.min.phy<-regstats[["Phylum"]][[1]][which.min(regstats[["Phylum"]][[1]]$logLoss),]
reg.min.class<-regstats[["Class"]][[1]][which.min(regstats[["Class"]][[1]]$logLoss),]
reg.min.order<-regstats[["Order"]][[1]][which.min(regstats[["Order"]][[1]]$logLoss),]
reg.min.family<-regstats[["Family"]][[1]][which.min(regstats[["Family"]][[1]]$logLoss),]
reg.min.genus<-regstats[["Genus"]][[1]][which.min(regstats[["Genus"]][[1]]$logLoss),]

# ASV models
loc.min.asv<-asv.local[["rf"]][["results"]][which.min(asv.local[["rf"]][["results"]]$logLoss),]
reg.min.asv<-asv.region[["rf"]][["results"]][which.min(asv.region[["rf"]][["results"]]$logLoss),]

locstats.df<-rbind(loc.min.phy,loc.min.class,loc.min.order,loc.min.family,loc.min.genus, loc.min.asv)
regstats.df<-rbind(reg.min.phy,reg.min.class,reg.min.order,reg.min.family,reg.min.genus, reg.min.asv)


# Create a dataframe formatted for plotting
resolution<-c("Phylum", "Class", "Order", "Family", "Genus", "ASV")
# Extract local stats
logLoss.loc<-locstats.df$logLoss
acc.loc<-locstats.df$Accuracy
# Extract regional stats
logLoss.reg<-regstats.df$logLoss
acc.reg<-regstats.df$Accuracy

# Bind to plottable dataframe
plot.loc.stats<-data.frame(resolution,logLoss.loc,acc.loc)
plot.reg.stats<-data.frame(resolution,logLoss.reg,acc.reg)
plot.loc.stats
locstats
loc.min.asv

colnames(locstats.df)
keeps<-c("Accuracy", "logLoss", "logLossSD", "Mean_Precision", "Mean_Recall", "mtry")
local.table2<-locstats.df[(names(locstats.df) %in% keeps)]
regional.table2<-regstats.df[(names(regstats.df) %in% keeps)]

localstats.df<-local.table2
regionalstats.df<-regional.table2

rownames(localstats.df)<-c("Phylum", "Class", "Order", "Family", "Genus", "ASV")
rownames(regionalstats.df)<-c("Phylum", "Class", "Order", "Family", "Genus", "ASV")

# Plot
require(plotly)
plot.loc.stats$resolution <- factor(plot.loc.stats$resolution, levels = c('Phylum', 'Class', 'Order', 'Family', "Genus", "ASV")) # this shows port then region
plot.reg.stats$resolution <- factor(plot.reg.stats$resolution, levels = c('Phylum', 'Class', 'Order', 'Family', "Genus", "ASV")) # this shows port then region

figyaxis <- list(title = "Relative Abundance",
                 showline = TRUE,
                 showgrid = TRUE,
                 showticklabels = TRUE,
                 linecolor = 'black',
                 linewidth = 1,
                 autotick = TRUE,
                 tick0=0,
                 dtick=5,
                 ticks = 'outside',
                 tickcolor = 'black',
                 tickwidth = 2,
                 ticklen = 5,
                 tickfont = list(family = 'Arial',
                                 size = 12,
                                 color = '#000000'))

figxaxis <- list(title = "Region",
                 showline = TRUE,
                 showgrid = TRUE,
                 showticklabels = TRUE,
                 linecolor = 'black',
                 linewidth = 1,
                 autotick = TRUE,
                 tick0=0,
                 dtick=5,
                 ticks = 'outside',
                 tickcolor = 'black',
                 tickwidth = 2,
                 ticklen = 5,
                 tickfont = list(family = 'Arial',
                                 size = 16,
                                 color = '#000000'))


# Region
plot.Fig6.region <- plot_ly(plot.reg.stats, x = ~resolution, y = ~acc.reg, type = 'scatter', mode ='markers', fillcolor='#F68C5B1A', marker = list(
  color = '#393B45',
  size = 8,
  line = list(
    color = '#393B45',
    width = 1
  )
),
line = list(color = '#393B45'),
marker=list(color='#393B45'),
showlegend = TRUE, name = 'Accuracy') %>%
  add_trace(y = ~logLoss.reg, type = 'scatter', mode = "markers",
            fill = 'tonexty', fillcolor='#F6A5384D', marker = list(
              color = '#FFFFFF',
              size = 8,
              line = list(
                color = '#393B45'
              )
            ),
            showlegend = TRUE, name = 'logLoss') %>%
  layout(title = "Region",
         paper_bgcolor='#FFFFFF', plot_bgcolor='#FFFFFF',
         xaxis = figxaxis,
         yaxis = figyaxis,
         legend=list(orientation='h'),
         title = "ML model metrics",
         xaxis = list(title = "Model Metric (%)"),
         margin = list(l = 65)
  )

plot.Fig6.region

# Local
plot.Fig6.port <- plot_ly(plot.loc.stats, x = ~resolution, y = ~acc.loc, type = 'scatter', mode ='markers', fillcolor='#F68C5B1A', marker = list(
  color = '#393B45',
  size = 8,
  line = list(
    color = '#393B45',
    width = 1
  )
),
line = list(color = '#393B45'),
marker=list(color='#393B45'),
showlegend = TRUE, name = 'Accuracy') %>%
  add_trace(y = ~logLoss.loc, type = 'scatter', mode = "markers",
            fill = 'tonexty', fillcolor='#D249514D', marker = list(
              color = '#FFFFFF',
              size = 8,
              line = list(
                color = '#393B45'
              )
            ),
            showlegend = TRUE, name = 'logLoss') %>%
  layout(title = "Port",
         paper_bgcolor='#FFFFFF', plot_bgcolor='#FFFFFF',
         xaxis = figxaxis,
         yaxis = figyaxis,
         legend=list(orientation='h'),
         title = "ML model metrics",
         xaxis = list(title = "Model Metric (%)"),
         margin = list(l = 65)
  )

plot.Fig6.port



# ==== Figure 7 ====

# Have to re-model (Machine learning) with taxonomic Class level at region with importance=TRUE
# gives the overall variable importance per region
seed=81
set.seed(seed)
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

pkgs <- c("hpgltools", "DESeq2", "plyr", "dplyr", "dbplyr","stringr", "phyloseq",
          "data.table", "tidyverse", "tibble", "caretEnsemble", "caret", "randomForest", "plyr", 
          "tidyverse", "data.table", "csv", "magrittr") 
invisible(sapply(pkgs,require, character = TRUE))

rm(list=ls()[! ls() %in% c("pm.ps4", "seed", "trctrlbase", "modtypesbase", "metric")])
seed=81
set.seed(seed)

pm.ps.class = tax_glom(pm.ps4, "Class", NArm = TRUE, bad_empty=c(NA, "", " ", "\t"))
pm.ps.class
pm.ps.class <- microbiome::transform(pm.ps.class, transform = "log10", target = "OTU", shift=1) # log10(x+1 tform)

require(reshape)
#Class
class.melt <- psmelt(pm.ps.class)
class.melt$Class <- as.character(class.melt$Class)
class.melt <- aggregate(Abundance~Sample+Class, class.melt, FUN=sum)
class.mat <- cast(class.melt, Sample ~ Class)
rownames(class.mat)<-class.mat[,1]
class.mat[1]<-NULL

# Create regional class labels
localind<-as.data.frame(rownames(class.mat))
names(localind)[1]<-"regionind"

localind$regionind<-gsub('[0-9]+', '', localind$regionind) # remove any numerical after character string
library(stringr)
#  Remove any unwanted patterns in the strings from 'sampleid' factor
localind<-unlist(localind) # convert to vector for stringr input
str(localind)
localind<-str_remove(localind, "post")
localind<-str_remove(localind, "point")
localind<-str_extract(localind, "^.{1,3}")
regionind<-as.data.frame(localind)
names(regionind)[1]<-"regionind"

regionind$regionind <- str_replace_all(regionind$regionind, # column to to search
                                       c("BUS" = "ASIA","SIN" = "ASIA","HK" = "ASIA", # ASIA
                                         "EUM" = "EUR","EUN" = "EUR","EUR" = "EUR","EUV" = "EUR", "EUW" = "EUR", # EUROPE
                                         "BAL" = "EAST","CHA" = "EAST","GAL" = "EAST","NEW" = "EAST","NOF" = "EAST","NY" = "EAST",  # EAST
                                         "LA" = "WEST","OAK" = "WEST","SEA" = "WEST",
                                         "DS" = "LAKES", "GB" = "LAKES", "KEW" = "LAKES") # WEST
)

regionind$regionind<-as.factor(regionind$regionind)
levels(regionind$regionind)
names(regionind)[1]<-"region"

# Append regional labels
class.mat2<-cbind(class.mat,regionind)

# Check class balance
percentage <- prop.table(table(class.mat2$region)) * 100
cbind(freq=table(class.mat2$region), percentage=percentage)

# Sanity check for append properly
class.mat2[39]

# Randomize observations
class.mat2 <- class.mat2[sample(nrow(class.mat2)),]
class.mat2[39]

# Model
region.class.vartrue <-caretList(region~.,
                                 data=class.mat2, ntree=501, importance=TRUE,
                                 trControl=trctrlbase, tuneList=modtypesbase)

region.class.vartrue
impvars.class.region<-varImp(region.class.vartrue$rf)
impvars.class.region


# Dataframe only those top 10
x<-c("Acidimicrobiia", "Bacteroidia", "Verrucomicrobiae", "Rhodothermia", "Actinobacteria",
     "Planctomycetacia", "Gammaproteobacteria", "Mollicutes", "Oxyphotobacteria", "Lentisphaeria")
y<-c(100.00, 97.88, 74.77, 74.54, 71.17, 61.90, 59.14, 58.63, 51.94, 51.94)
vImp.df<-as.data.frame(cbind(x, y))
vImp.df[,2] <- as.numeric(as.character(vImp.df[,2])) # Have to use as.char before as.num, factors are stored internally as integers, need to apply factor level labels (actual values) or it will convert to integer codes, non numeric chars such as . ruins this.
class(vImp.df$x) # Factor
class(vImp.df$y)

vImp.df=vImp.df %>%
  arrange(y) %>%
  mutate(x=factor(x,x))

# Plot variable importance of these top 10
plot.Fig7.varImp = ggplot(vImp.df, aes(x=x, y=y)) +
  geom_segment( aes(x=x, xend=x, y=0, yend=y ), color=ifelse(vImp.df$x %in% c("Acidimicrobiia","Bacteroidia","Actinobacteria","Oxyphotobacteria", "Gammaproteobacteria"), "#F7756D", "#00BEC4"), size=ifelse(vImp.df$x %in% c("A","D"), .08, .5) ) +
  geom_point( color=ifelse(vImp.df$x %in% c("Acidimicrobiia","Bacteroidia","Actinobacteria","Oxyphotobacteria", "Gammaproteobacteria" ), "#F7756D", "#00BEC4"), size=ifelse(vImp.df$x %in% c("A","D"), 5, 3) ) +
  theme_light() +
  coord_flip() +
  theme(
    legend.position="none",
    panel.grid.major.y = element_blank(),
    panel.border = element_blank(),
    axis.ticks.y = element_blank()
  ) +
  xlab("Class Predictor") +
  ylab("Overall Variable Importance (%)")
plot.Fig7.varImp


# Provide heatmap to attach to varImp class plot
impvars.class.region
# Order variable importance by by cl.varImp: consistent with previous graphs: east coast/west coast/europe/asia/lakes
varImp.region<-rbind(c(44.78, 59.58, 100.00, 33.91, 70.71),c(57.30, 84.48,   97.88, 42.31, 56.89), #Acidimicrobiia,Bacteroidia
                     c(74.65, 74.77,  40.04, 26.31, 63.95), c(50.08, 53.31,  37.18, 74.54, 49.19), #Verrucomicrobiae,Rhodothermia
                     c(71.17, 43.10,  39.89, 40.83, 70.83), c(49.67, 57.69,  55.17, 25.59, 61.90), #Actinobacteria, Planctomycetacia
                     c(59.14, 46.31,  39.45, 33.12, 55.44), c(58.63, 25.64,  35.17, 20.66, 35.67), #Gammaproteobacteria, Mollicutes
                     c(37.84, 43.13,  51.94, 19.41, 49.11), c(51.94, 31.46,  28.36, 16.39, 23.38)) #Oxyphotobacteria, Lentisphaeria

class(varImp.region)

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

varImp.region<-as.data.frame(varImp.region)
#varImp.region<-as.matrix(varImp.region)
#change in order east west europe asia lakes
varImp.region.arr<-arrange.vars(varImp.region, c("V2"=1, "V5"=2, "V3"=3, "V1"=4, "V4"=5))
varImp.region.arr<-as.matrix(varImp.region.arr)


X<-c("East Coast U.S", "West Coast U.S", "Europe", "Asia", "Great Lakes")

plot.Fig7.heatmap.small <- plot_ly(
  x = ~X,
  z = varImp.region.arr, color=varImp.region.arr, type = "heatmap",   showticklabels = TRUE, autosize=F, width=250, height=700)

plot.Fig7.heatmap.small # Use this for actual block size
# * Needs to be reflected/rotated 180 to append to a varImp plot.



# ==== Figure 8 =====
top10taxa<-c("Acidimicrobiia", "Bacteroidia", "Verrucomicrobiae", "Rhodothermia", "Actinobacteria",
             "Planctomycetacia", "Gammaproteobacteria", "Mollicutes", "Oxyphotobacteria", "Lentisphaeria")
top10taxa<-as.vector(top10taxa)

# Create dataframe with Enrichment Factor (EF) data
# These are the top 10 taxa
df<-data.frame("Acidimicrobiia"=c(12,4,12,4,9),
               "Bacteroidia"=c(0,0,0,0,0),
               "Verrucomicrobiae"=c(7,7,8,1,7),
               "Rhodothermia"=c(15,12,12,7,0),
               "Actinobacteria"=c(14,1,13,3,14),
               "Planctomycetacia"=c(13,5,6,4,10),
               "Gammaproteobacteria"=c(2,0,1,10,8),
               "Mollicutes"=c(15,10,16,17,0),
               "Oxyphotobacteria"=c(16,3,2,11,13),
               "Lentisphaeria"=c(0,10,6,19,0))
rows<-c("East", "West", "Europe", "Asia", "Lakes")
row.names(df)<-rows

df <-rbind(rep(20,10) , rep(0,10) , df)
colnames(df) <- c("1" , "2" , "3" , "4" , "5", "6", "7", "8", "9", "10")


# New colors
colors1<-c("#CC0C00FF" ,"#FF7F0EFF", "#FFCD00FF" ,"#00B5E2FF" ,"#00AF66FF" ,"#8755B5FF" ,"#FFFFFF")

colors_border=c("#CC0C00FF", "#FF7F0EFF", "#FFCD00FF", "#00B5E2FF","#00AF66FF")
colors_in=c("#CC0C0066", "#FF7F0E66", "#FFCD0066", "#00B5E266", "#00AF6666") #Add 40% opacity to inside (66)

# Prepare title
mytitle <- c("East Coast U.S", "West Coast U.S", "Europe", "Asia", "Great Lakes")

# Split the screen in 6 parts
par(mar=rep(0.8,4))
par(mfrow=c(2,3))

library(fmsb)
# Loop for each plot
for(i in 1:5){
  # Customize the polar chart here
  radarchart(df[c(1,2,i+2),], axistype=1, 
             # Custom polygon
             pcol=colors_border[i] , pfcol=colors_in[i] , plwd=2, plty=1 , 
             # Custom the grid
             cglcol="grey", cglty=1, axislabcol="black", caxislabels=seq(0,20,5), cglwd=0.8,
             # Custom labels
             vlcex=0.8,
             # Title
             title=mytitle[i]
  )
}




# ============================== Supplementary Figures ====================================

# ==== Figure S1 ====
# Refer to Figure 3 of how  generate plot.data2 dataframe was generated
plot.data2<-'https://github.com/rghannam/pm_workflow/raw/master/data/graphing_inputs/fig3/boxplot_formatted.csv'
plot.data2<-read.csv(url(plot.data2))

seed=81
set.seed(seed)
pkgs <- c("phyloseq", "ggplot2", "rlang", "Rcpp", "vegan", "csv", "plotly", "plyr", "dplyr", "reshape", "reshape2") # Load packages
invisible(sapply(pkgs,require, character = TRUE))
levels(plot.data2$group)

# Order the variables
plot.data2$group <- factor(plot.data2$group, levels = c('Acidimicrobiia', 'Actinobacteria', 'Alphaproteobacteria', 'Bacteroidia', 'Gammaproteobacteria', 'Oxyphotobacteria'))
plot.data2$region <- factor(plot.data2$region, levels = c('East', 'West', 'Europe', 'Asia', 'Great Lakes'))

# Plot
require(plotly)
figyaxis <- list(title = "Port",
                 showline = TRUE,
                 showgrid = TRUE,
                 showticklabels = TRUE,
                 linecolor = 'black',
                 linewidth = 1,
                 tickformat = "%",
                 autotick = TRUE,
                 ticks = 'outside',
                 tickcolor = 'black',
                 tickwidth = 2,
                 ticklen = 5,
                 tickfont = list(family = 'Arial',
                                 size = 14,
                                 color = 'rgb(82, 82, 82)'))

figxaxis <- list(title = "Relative Abundance (% all reads)",
                 showline = TRUE,
                 showgrid = TRUE,
                 showticklabels = TRUE,
                 linecolor = 'black',
                 linewidth = 1,
                 autotick = TRUE,
                 tickformat = "%",
                 ticks = 'outside',
                 tickcolor = 'black',
                 tickwidth = 2,
                 ticklen = 5,
                 tickfont = list(family = 'Arial',
                                 size = 14,
                                 color = 'rgb(82, 82, 82)'))


# Use main data.frame to subset and plot each class and its dispersion
bp.Acidimicrobiia<- plot.data2[plot.data2$group == "Acidimicrobiia", ] #1
as.data.frame(rownames(bp.Acidimicrobiia))

bp.Actinobacteria<- plot.data2[plot.data2$group == "Actinobacteria", ] #2
as.data.frame(rownames(bp.Actinobacteria))

bp.Alphaproteobacteria<- plot.data2[plot.data2$group == "Alphaproteobacteria", ] #3
as.data.frame(rownames(bp.Alphaproteobacteria))

bp.Gammaproteobacteria<- plot.data2[plot.data2$group == "Gammaproteobacteria", ] #4
as.data.frame(rownames(bp.Gammaproteobacteria))

bp.Bacteroidia<- plot.data2[plot.data2$group == "Bacteroidia", ] #5
as.data.frame(rownames(bp.Bacteroidia))

bp.Oxyphotobacteria<- plot.data2[plot.data2$group == "Oxyphotobacteria", ] #6
as.data.frame(rownames(bp.Oxyphotobacteria))

# Extract out one graph for each class
bp.1 <- plot_ly(bp.Acidimicrobiia, y = ~abundance, x=~region, type="violin", color = I("#9D1B44")) %>%
  layout(title = "Acidimicrobiia", xaxis = figxaxis, yaxis=figyaxis, showlegend=F, autosize=F, width=250, height=250)
bp.1

bp.2 <- plot_ly(bp.Actinobacteria, y = ~abundance, x=~region, type="violin", color = I("#F46D43")) %>%
  layout(title = "Actinobacteria", xaxis = figxaxis, yaxis=figyaxis, showlegend=F, autosize=F, width=250, height=250)
bp.2

bp.3 <- plot_ly(bp.Alphaproteobacteria, y = ~abundance, x=~region, type="violin", color = I("#FEE08B")) %>%
  layout(title = "Alphaproteobacteria", xaxis = figxaxis, yaxis=figyaxis, showlegend=F, autosize=F, width=250, height=250)
bp.3

bp.4 <- plot_ly(bp.Gammaproteobacteria, y = ~abundance, x=~region, type="violin", color = I("#66C2A5")) %>%
  layout(title = "Gammaproteobacteria", xaxis = figxaxis, yaxis=figyaxis, showlegend=F, autosize=F, width=250, height=250)
bp.4

bp.5 <- plot_ly(bp.Bacteroidia, y = ~abundance, x=~region, type="violin", color = I("#3288BD")) %>%
  layout(title = "Bacteroidia", xaxis = figxaxis, yaxis=figyaxis, showlegend=F, autosize=F, width=250, height=250)
bp.5

bp.6 <- plot_ly(bp.Oxyphotobacteria, y = ~abundance, x=~region, type="violin", color = I("#5E4FA2")) %>%
  layout(title = "Oxyphotobacteria", xaxis = figxaxis, yaxis=figyaxis, showlegend=F, autosize=F, width=250, height=250)
bp.6



# ==== Figure S2 ====
rm(list=ls()[!ls() %in% c("seed", "pm.ps4")])
set.seed(seed)

colnames(pm.ps4@sam_data)

# Extract out counts from phyloseq object as a dataframe
physeq.otu = as(otu_table(pm.ps4), "matrix")
if(taxa_are_rows(pm.ps4)){physeq.otu <- t(physeq.otu)}
physeq.otu.df = as.data.frame(physeq.otu)
head(physeq.otu.df)[1:5]
dim(physeq.otu.df)

# Find distances
require(vegan)
abund_dist <- vegdist(physeq.otu.df) # Bray
abund_dist

# Run analysis of similarity (ANOSIM) at both spatial scales
abund_ano.Local <- anosim(abund_dist, pm.ps4@sam_data$geo_port) # Use geo_port factor for local scales
abund_ano.Local
plot(abund_ano.Local)

abund_ano.Region <- anosim(abund_dist, pm.ps4@sam_data$geo_region) # Use geo_region factor for regional scales
abund_ano.Region
plot(abund_ano.Region)

samdata<-pm.ps4@sam_data
samdata$geo_region<-as.factor(samdata$geo_region)
samdata$geo_port<-as.factor(samdata$geo_port)

levels(samdata$geo_port)
levels(samdata$geo_region)



# ==== Figure S3 ====
pm.ps4@sam_data$geo_region
pm.ps4@sam_data$geo_region <- factor(pm.ps4@sam_data$geo_region, levels = c('East Coast U.S', 'West Coast U.S', 'Europe', 'Asia', 'Great Lakes'))

pm.ps.dist = ordinate(pm.ps4, method = "PCoA", distance = "jaccard")
pm.ps.dist

pcoa = plot_ordination(pm.ps4, pm.ps.dist,
                       type="sample_name", 
                       color ="geo_region") +
  theme_bw() +
  stat_ellipse(geom = "polygon", 
               alpha = 0.2, 
               aes(fill = geo_region))


pcoa

#pcoa[["plot_env"]][["DF"]] # apply axis values to each sample.



# ==== Figure S4 ====
# Refer to Figure 4 to create this object (vImp.df)
# The frequency polygon will be created from this

# Plot all of the 68 important vars, respectively
vimpregion<-as.vector(vImp.df$vimpRegion)
vimpport<-as.vector(vImp.df$vimpPort)
classes<-as.vector(vImp.df$Class)
color_pallete<-c("#CC0C00FF" ,"#FF7F0EFF", "#FFCD00FF" ,"#00B5E2FF" ,"#00AF66FF" ,"#8755B5FF", "#000000")
vimpregion
vimpport

df<-data.frame(vimpregion,vimpport, classes )

vimpport<-ggplot(df, aes(x=vimpport, color=classes)) +
  geom_freqpoly(binwidth = 2) +
  scale_color_manual(labels = c('Acidimicrobiia', 'Actinobacteria', 'Alphaproteobacteria', 'Bacteroidia', 'Gammaproteobacteria', 'Oxyphotobacteria', 'Other'), 
                     values = color_pallete) +
  theme(legend.position="top") +   theme_minimal()
vimpport

vimpregion<-ggplot(df, aes(x=vimpregion, color=classes)) +
  geom_freqpoly(binwidth = 2) +
  scale_color_manual(labels = c('Acidimicrobiia', 'Actinobacteria', 'Alphaproteobacteria', 'Bacteroidia', 'Gammaproteobacteria', 'Oxyphotobacteria', 'Other'), 
                     values = color_pallete) +
  theme(legend.position="top") +   theme_minimal()
vimpregion



# ==== Figure S5 ====
rm(list=ls()[!ls() %in% c("seed", "pm.ps4")])

head(otu_table(pm.ps4)[1:5,1:10])
physeq.otu<-as(otu_table(pm.ps4), "matrix")
physeq.otu = as.data.frame(physeq.otu)
seed=81
set.seed(seed)

pm.taxa<-'https://github.com/rghannam/pm_workflow/raw/master/data/counts_taxonomy/taxa.rds'
pm.metadata<-'https://github.com/rghannam/pm_workflow/raw/master/data/metadata/pm_metadata.csv'
pm.taxa<-readRDS(url(pm.taxa))
pm.metadata<-read.csv(url(pm.metadata))

row.names(pm.metadata)<-pm.metadata$sample_name
row.names(pm.metadata)<-gsub("-", "", row.names(pm.metadata))
pm.metadata<-pm.metadata[-1]
library(dplyr)
library(tibble)
pm.metadata <- pm.metadata %>% rownames_to_column("sample_name")
row.names(pm.metadata)<-pm.metadata$sample_name


# Perform string manipulations on dataframe to rename the DNA string from which the ASV belongs to
# Original ASV table
asvtable<-'https://github.com/rghannam/pm_workflow/raw/master/data/counts_taxonomy/seqtabmergedNoC.rds'
asvtable<-readRDS(url(asvtable))
dim(asvtable)

colnames(physeq.otu)
# Extract the features used in machine learning regional model
library(tidyverse)
spvars<-colnames(physeq.otu)
spvars<-str_remove(spvars, "sp") # Remove the sp

# Based on position, match the sp to nucleotide string from asvtable
spvars<-as.numeric(spvars) # cCnvert to numeric so indexing remains in bounds
dat<-asvtable[,spvars] # Extract those columns based on number (position) of the spvars vector
dim(dat) # Same dim as spvars

# Sanity checks
dat2<-colnames(dat)
dat.feat<-dat2[10]

dnastring<-colnames(asvtable)  # Get colnames from original count table
og.feat<-dnastring[12] # 3rd position in original is 1st in dat2..since the regional model had sp3 as our first column/feature etc..., then 14..., 14 from original should match to 2 of dat
dat.feat
og.feat

# dat2 has the nucleotide strings of sp.

if(!identical(dat.feat, og.feat)) stop("string not identical")

# Create regional dataframe with nucleotide strings as features/column headers
physeq.otu.df<-as.data.frame(dat)
dim(physeq.otu.df)

# Match rownames too.
rows<-rownames(physeq.otu.df)
physeq.otu.df2<-physeq.otu.df[rows,]
dim(physeq.otu.df2)
# SANITY check number 2
# 14 of original should match to 2 of physeq.otu.df2
og<-dnastring[11]
subset<-colnames(physeq.otu.df2)
subset<-subset[9]
if(!identical(subset, og)) stop("string not identical")
dim(physeq.otu.df2)

row.names(physeq.otu.df2)<-gsub("-", "", row.names(physeq.otu.df2)) # Remove - from sample ids so they match to metadata or rowvectors of samples containing environmental variables later on
# End perform string manipulation, now the sp's are renamed to actual DNA strings

# Only keep those samples with env data.
env_samples_vec<-read.csv('https://github.com/rghannam/pm_workflow/raw/master/data/graphing_inputs/fig5S-6S/env_samples_vec.csv', header=FALSE)

row.names(physeq.otu.df2)[1:200]
rows1218<-row.names(physeq.otu)

pm.metadata <- pm.metadata[rownames(pm.metadata) %in% rows1218, ] # Match metadata samples with those in count table
physeq.otu.df2<-physeq.otu.df2[rownames(physeq.otu.df2) %in% rows1218, ]

pm.metadata2 <- pm.metadata[rownames(pm.metadata) %in% env_samples_vec$V1, ]
physeq.otu.df3 <- physeq.otu.df2[rownames(physeq.otu.df2) %in% env_samples_vec$V1, ]
row.names(pm.metadata2)[1:200] # Match
row.names(physeq.otu.df3)[1:200] # Match

colnames(pm.metadata2)
# Since we're plotting physiochemical conditions, define columns to keep in metadata
keeps <- c("sample_name", "conduc","salinity", "total_dissolved_solids",
           "surf_temp", "ph", "diss_oxygen", "geo_region")

pm.metadata3<-pm.metadata2[(names(pm.metadata2) %in% keeps)]

# Remove the units that were required for attributes in metadata
head(pm.metadata3)[1:7]
pm.metadata3$conduc<-gsub(" uS/cm", "", as.character(pm.metadata3$conduc) , n)
pm.metadata3$salinity<-gsub(" psu", "", as.character(pm.metadata3$salinity) , n)
pm.metadata3$total_dissolved_solids<-gsub(" mg/L", "", as.character(pm.metadata3$total_dissolved_solids) , n)
pm.metadata3$surf_temp<-gsub(" deg C", "", as.character(pm.metadata3$surf_temp) , n)
pm.metadata3$diss_oxygen<-gsub(" mg/L", "", as.character(pm.metadata3$diss_oxygen) , n)

# Rename levels in geo_region
library(plyr)
pm.metadata3$geo_region <- revalue(pm.metadata3$geo_region, c("West Coast U.S"="West"))
pm.metadata3$geo_region <- revalue(pm.metadata3$geo_region, c("East Coast U.S"="East"))
pm.metadata3$geo_region <- revalue(pm.metadata3$geo_region, c("Great Lakes"="Lakes"))

# Create phyloseq object of samples that we have environmental conditions for
levels(pm.metadata3$geo_region)
row.names(physeq.otu.df3)[1:1176]
row.names(pm.metadata3)[1:1176]
head(physeq.otu.df3)[1:10]
head(pm.metadata3)[1:8]

pm.taxa_t<-t(pm.taxa)
dim(pm.taxa_t)
pm.taxa_re_t<-t(pm.taxa_t)
dim(pm.taxa_re_t)
dim(pm.metadata3) # [1] 1176   8
pm.metadata3[, 1] <- as.factor(pm.metadata3[, 1])

sampledata = sample_data(pm.metadata3)
counts = otu_table(physeq.otu.df3, taxa_are_rows = FALSE) # leave false because we arbitrarily set column names with colnames()->NULL, we can always map back to respective ASV
taxa = tax_table(pm.taxa_re_t)
sampledata = sample_data(pm.metadata)
pm.ps.env = phyloseq(counts, taxa, sampledata)
pm.ps.env
otu_table(pm.ps.env)[1:5, 1:10]
head(otu_table(pm.ps4)[1:5,1:10])

length(get_taxa_unique(pm.ps4, taxonomic.rank = "Class")) # 38
length(get_taxa_unique(pm.ps.env, taxonomic.rank = "Class")) # 38
uniqcl<-get_taxa_unique(pm.ps.env, taxonomic.rank=rank_names(pm.ps.env)[3], errorIfNULL=TRUE) # Get unique classes
uniqcl

pm.ps.env.normed = tax_glom(pm.ps.env, "Class", NArm = TRUE, bad_empty=c(NA, "", " ", "\t"))
pm.ps.env.normed
pm.ps.env.normed <- microbiome::transform(pm.ps.env.normed, transform = "log10", target = "OTU", shift=1) # log10(x+1 tform)
pm.ps.env.normed


# Use package microbiomeSeq for this
seed=81
set.seed(seed)
library(phyloseq)
library(ggplot2)
library(microbiomeSeq)

colnames(pm.ps.env.normed@sam_data)[colnames(pm.ps.env.normed@sam_data)=="geo_region"]<-"region"

physeq<-pm.ps.env.normed

# Arrange environmental vars
physeq@sam_data
arrange.vars <- function(data, vars){
  # Stop if not a data.frame (but should work for matrices as well)
  stopifnot(is.data.frame(data))
  
  # Sort out inputs
  data.nms <- names(data)
  var.nr <- length(data.nms)
  var.nms <- names(vars)
  var.pos <- vars
  # Sanity checks
  stopifnot( !any(duplicated(var.nms)),
             !any(duplicated(var.pos)) )
  stopifnot( is.character(var.nms),
             is.numeric(var.pos) )
  stopifnot( all(var.nms %in% data.nms) )
  stopifnot( all(var.pos > 0),
             all(var.pos <= var.nr) )
  
  # Prepare output
  out.vec <- character(var.nr)
  out.vec[var.pos] <- var.nms
  out.vec[-var.pos] <- data.nms[ !(data.nms %in% var.nms) ]
  stopifnot( length(out.vec)==var.nr )
  
  # Re-arrange vars by position
  data <- data[ , out.vec]
  return(data)
}

# Remove all character strings - sanity
colnames(physeq@sam_data)
physeq@sam_data$ph
physeq@sam_data$conduc<-gsub("[^0-9\\.]", "", physeq@sam_data$conduc)
physeq@sam_data$diss_oxygen<-gsub("[^0-9\\.]", "", physeq@sam_data$diss_oxygen)
physeq@sam_data$salinity<-gsub("[^0-9\\.]", "", physeq@sam_data$salinity)
physeq@sam_data$total_dissolved_solids<-gsub("[^0-9\\.]", "", physeq@sam_data$total_dissolved_solids)
physeq@sam_data$surf_temp<-gsub("[^0-9\\.]", "", physeq@sam_data$surf_temp)

# Keep only those environmental variables we're interested in
keeps<-c("conduc","diss_oxygen","ph","salinity","total_dissolved_solids","surf_temp","region")
physeq@sam_data<-physeq@sam_data[,keeps]

# Convert to numeric
physeq@sam_data$conduc<-as.numeric(physeq@sam_data$conduc)
physeq@sam_data$diss_oxygen<-as.numeric(physeq@sam_data$diss_oxygen)
physeq@sam_data$ph<-as.numeric(physeq@sam_data$ph)
physeq@sam_data$salinity<-as.numeric(physeq@sam_data$salinity)
physeq@sam_data$total_dissolved_solids<-as.numeric(physeq@sam_data$total_dissolved_solids)
physeq@sam_data$surf_temp<-as.numeric(physeq@sam_data$surf_temp)

# Arrange vars
phy.arranged<-arrange.vars(physeq@sam_data, c("conduc"=1, "diss_oxygen"=2, "ph"=3,
                                              "salinity"=4,"total_dissolved_solids"=5,"surf_temp"=6, "region"=7))

# Plot
phy.arranged$region<-as.factor(phy.arranged$region)
plot.FigS5 <- plot_anova_env(phy.arranged, grouping_column = "region", pValueCutoff = 0.001, 
                             select.variables = c("conduc", "diss_oxygen", 
                                                  "ph", "salinity", "total_dissolved_solids", "surf_temp"))


print(plot.FigS5)



# ==== Figure S6 ====
physeq <- taxa_level(physeq, "Class")
env.taxa.cor <- taxa.env.correlation(physeq, grouping_column = "region", method = "pearson", 
                                     pvalue.threshold = 0.001, padjust.method = "BH", adjustment = 5, num.taxa = 43, 
                                     select.variables = NULL)


plot.FigS6 <- plot_taxa_env(env.taxa.cor)
print(plot.FigS6)



# ==== Figure S7 ====
# Find significant environmental variables

# ADONIS

seed=81
set.seed=seed

# logarithmize the counts: log10(x+1)
pm.ps.env.log10 <- microbiome::transform(pm.ps.env, transform = "log10", target = "OTU", shift=1) # log10(x+1 tform)
pm.ps.env.log10
head(otu_table(pm.ps.env.log10)[1:5,1:5])

physeq.otu.log10<-as(otu_table(pm.ps.env.log10), "matrix")
physeq.otu.log10 = as.data.frame(physeq.otu.log10)

env.counts<-physeq.otu.log10
env.metadata<-pm.metadata3

# Row names to column in env.metadata
env.metadata$sample_name<-rownames(env.metadata)

samples.ordered<-rownames(env.counts)
samples.ordered2<-as.data.frame(samples.ordered)
env.metadata<-env.metadata[match(samples.ordered, env.metadata$sample_name),] # Arranged rows in metadata by the rows in count table

# Remove unwanted columns
colnames(env.metadata)
env.metadata[1]<-NULL # sample_name
env.metadata[7]<-NULL # geo_region

head(rownames(env.metadata)[1:5])
head(rownames(env.counts)[1:5])

# Just a check to ensure that the samples in meta_table are in the same order as in abund_table
env.metadata2<-env.metadata[rownames(env.counts),]

# Filter out any samples taxas that have zero entries 
env.counts2<-subset(env.counts,rowSums(env.counts)!=0)

dim(env.counts2)
dim(env.metadata2)
head(env.counts2)[1:3]
head(env.metadata2)[1:6]

# Can convert to relative frequencies but notrun
# abund_table<-abund_table/rowSums(abund_table)

# Use PERMANOVA to find significant environmental variables
library(vegan)
table.adonis <- adonis(env.counts2 ~ ., data=env.metadata2, method = "bray")
table.adonis

(0.83333+0.00286+0.00216)*100

# Calculate bray curtis distance matrix just to test statistics on regions 
# Use pm.ps.env.normed.justenvironment - normalized phyloseq object
rm(list=ls()[!ls() %in% c("seed", "pm.ps.env.normed", "pm.ps4", "physeq",
                          "pm.ps.env", "pm.metadata3", "pm.ps.env.log10")])
set.seed(seed)
distanceMethodList

# Manipulate ps object
colnames(pm.ps.env.log10@sam_data)
physeq<-pm.ps.env.log10
colnames(physeq@sam_data)[colnames(physeq@sam_data)=="geo_region"]<-"region"

# Remove all character strings - sanity
physeq@sam_data$conduc<-gsub("[^0-9\\.]", "", physeq@sam_data$conduc)
physeq@sam_data$diss_oxygen<-gsub("[^0-9\\.]", "", physeq@sam_data$diss_oxygen)
physeq@sam_data$salinity<-gsub("[^0-9\\.]", "", physeq@sam_data$salinity)
physeq@sam_data$total_dissolved_solids<-gsub("[^0-9\\.]", "", physeq@sam_data$total_dissolved_solids)
physeq@sam_data$surf_temp<-gsub("[^0-9\\.]", "", physeq@sam_data$surf_temp)

# Keep only those environmental variables we're interested in
keeps<-c("conduc","diss_oxygen","ph","salinity","total_dissolved_solids","surf_temp","region")
physeq@sam_data<-physeq@sam_data[,keeps]

# Convert to numeric
physeq@sam_data$conduc<-as.numeric(physeq@sam_data$conduc)
physeq@sam_data$diss_oxygen<-as.numeric(physeq@sam_data$diss_oxygen)
physeq@sam_data$ph<-as.numeric(physeq@sam_data$ph)
physeq@sam_data$salinity<-as.numeric(physeq@sam_data$salinity)
physeq@sam_data$total_dissolved_solids<-as.numeric(physeq@sam_data$total_dissolved_solids)
physeq@sam_data$surf_temp<-as.numeric(physeq@sam_data$surf_temp)

# Calculate bray distances of physloeq object
pmps_bray <- phyloseq::distance(physeq, method = "bray")
pmps_bray

# Make a data frame from the sample_data
sampledf <- data.frame(sample_data(physeq))

require(vegan)
# Adonis test
# Do the 5 regions we collected samples from have different centroids?
adregion<-adonis(pmps_bray ~ region, data = sampledf) # Using the geo_region (renamed region) var from metadata
adregion

# Check dispersions on region
beta <- betadisper(pmps_bray, sampledf$region)
beta
plot(beta)
permutest(beta)

anova(beta)

physeq@sam_data
physeq@sam_data$region
physeq@sam_data$region <- factor(physeq@sam_data$region, levels = c('East Coast U.S', 'West Coast U.S', 'Europe', 'Asia', 'Great Lakes'))
physeq

# Constrained Analysis of Principal Coordinates (CAP)
colnames(physeq@sam_data)
# Remove data points with missing metadata
bray_not_na <- physeq %>%
  subset_samples(
    !is.na(conduc) & 
      !is.na(salinity) &
      !is.na(total_dissolved_solids) & 
      !is.na(surf_temp) & 
      !is.na(ph) & 
      !is.na(diss_oxygen)
  )
bray_not_na # Contains all environmental variables

# CAP ordinate
cap_ord <- ordinate(
  physeq = physeq, 
  method = "CAP",
  distance = pmps_bray,
  formula = ~ conduc + salinity + total_dissolved_solids + surf_temp + ph + diss_oxygen
)
cap_ord
# CAP plot
cap_plot <- plot_ordination(
  physeq = physeq, 
  ordination = cap_ord, 
  color = "region", 
  axes = c(1,2)
) + 
  #aes(shape = Station) + 
  #colorsspectral=c('#9E0142', '#F46D43', '#FEE08B','#66C2A5', '#3288BD', '#5E4FA2')
  geom_point(aes(colour = region), alpha = 1, size = 4) + 
  geom_point(colour = "white", size = 1.5) + 
  theme_bw() +
  scale_color_manual(values = c("#F8766D", "#DACE00", "#00BF7D", "#00B0F6", 
                                "#E76BF3")
  )
cap_plot

# Customize the points
cap_plot$layers <- cap_plot$layers[-1]
cap_plot
cap_plot2<- cap_plot + geom_point(size=3)
cap_plot2
cap_plot2<- cap_plot + geom_point(size=3, alpha=0.8)
cap_plot2

# Now add the environmental variables as arrows
arrowmat <- vegan::scores(cap_ord, display = "bp")
arrowmat
# Add labels, make a data.frame
arrowdf <- data.frame(labels = rownames(arrowmat), arrowmat)

# Define the arrow aesthetic mapping
arrow_map <- aes(xend = CAP1, 
                 yend = CAP2, 
                 x = 0, 
                 y = 0, 
                 shape = NULL, 
                 color = NULL, 
                 label = labels)

label_map <- aes(x = 1.3 * CAP1, 
                 y = 1.3 * CAP2, 
                 shape = NULL, 
                 color = NULL, 
                 label = labels)

arrowhead = arrow(length = unit(0.02, "npc"))

# Make a new graphic
cap_plot2 + 
  geom_segment(
    mapping = arrow_map, 
    size = .5, 
    data = arrowdf, 
    color = "black", 
    arrow = arrowhead
  ) + 
  geom_text(
    mapping = label_map, 
    size = 4,  
    data = arrowdf, 
    show.legend = FALSE
  )


anovacapord<-anova(cap_ord)
anovacapord
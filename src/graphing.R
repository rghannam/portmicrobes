rm(list=ls())

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
library("phyloseq")

seed=81
set.seed(seed)


# ---- Generate Figures for the Port Microbes workflow ----



# ---- Figure 1 ----

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



# ---- Figure 2 ----

pkgs <- c("phyloseq", "ggplot2", "rlang", "Rcpp", "vegan", "csv", "plotly", "plyr", "dplyr", "reshape", "reshape2") # Load packages
invisible(sapply(pkgs,require, character = TRUE))

# --- Read in data ---

# We'll be using the 1218x119097 count table from data_preprocessing.R

# Call counts/taxonomy/metadata from GitHub
pm.counts<-'https://github.com/rghannam/pm_workflow/raw/master/data/counts_taxonomy/pm_counts1218.rds'
pm.counts<-readRDS(url(pm.counts))
dim(pm.counts)

pm.taxa<-'https://github.com/rghannam/pm_workflow/raw/master/data/counts_taxonomy/pm_taxonomy.rds'
pm.taxa<-readRDS(url(pm.taxa))
dim(pm.taxa)

pm.metadata<-'https://github.com/rghannam/pm_workflow/raw/master/data/metadata/pm_metadata.csv'
pm.metadata<-read.csv(url(pm.metadata))
dim(pm.metadata)

# Manipulate data so we can convert to a phyloseq object properly and variable names match downstream.
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

# Merge phyloseq object by geo_port (20 locations/levels) variable from metadata
colnames(sampledata)
merged.pm.ps.location = merge_samples(pm.ps, "geo_port")
merged.pm.ps.location
otu_table(merged.pm.ps.location)[1:5, 1:5]

c<-subset_taxa(merged.pm.ps.location, !is.na(Class) & !Class %in% c("", "uncharacterized"))
x<-list() # Store in list as habit to do all levels of taxonomy, here we just do class
alltaxa<-list(c, x)
listNames <- c("class") # Renames to proper taxonomic level
names(alltaxa) <- listNames


# --- Compute prevalence of each bacterial feature ---
prev.abund.class = apply(X = otu_table(alltaxa$class),
                         MARGIN = ifelse(taxa_are_rows(alltaxa$class), yes = 1, no = 2),
                         FUN = function(x){sum(x > 0)})
prev.abund.class = data.frame(Prevalence = prev.abund.class,
                              TotalAbundance = taxa_sums(alltaxa$class),
                              tax_table(alltaxa$class))
#write.csv(prev.abund.class, '~/pm_workflow/.../output/prev.abund.class.csv', row.names = TRUE)


# --- Which bacterial class consist of low prevalence ASVs? ---
pm.totalavgprev.class <- plyr::ddply(prev.abund.class, "Class", function(df1){cbind(mean(df1$Prevalence),sum(df1$Prevalence))})
as.data.frame(pm.totalavgprev.class)
names(pm.totalavgprev.class) <- c("Class", "mean", "total")
#write.csv(pm.totalavgprev.class, '~/pm_workflow/.../output/pm.totalavgprev.class', row.names = TRUE)

# Filter entries with <=15 prevalence
pm.totalavgprev.class.subset <- subset(pm.totalavgprev.class, pm.totalavgprev.class$total <= 15)
extract.class<-as.data.frame(pm.totalavgprev.class.subset[,1])
names(extract.class) <- ("Class") 
extract.class <- as.list(as.vector(t(extract.class)))

filterclass = (extract.class)
pm.ps.class.sub = subset_taxa(alltaxa$class, !Class %in% filterclass)
pm.ps.class.sub
alltaxa$class

length(get_taxa_unique(alltaxa$class, taxonomic.rank = "Class")) # 143 unique classes before filtering
length(get_taxa_unique(pm.ps.class.sub, taxonomic.rank = "Class")) # 87 unique classes after filtering

# Create a phyloseq object agglomerated by bacterial class
pm.ps.class = tax_glom(pm.ps.class.sub, "Class", NArm = TRUE, bad_empty=c(NA, "", " ", "\t"))
pm.ps.class

#saveRDS(pm.ps.class, file = "~/pm_workflow/.../output/pm.ps.class.rds")


# Create a count matrix of class features by each sampling location
library(reshape)
pm.ps.class.melt <- psmelt(pm.ps.class)
pm.ps.class.melt$Class <- as.character(pm.ps.class.melt$Class)
pm.ps.class.melt
pm.ps.class.melt <- aggregate(Abundance~Sample+Class, pm.ps.class.melt, FUN=sum)
pm.ClassMatrix.loc <- cast(pm.ps.class.melt, Sample ~ Class)
dim(pm.ClassMatrix.loc)
head(pm.ClassMatrix.loc)
write.csv(pm.ClassMatrix.loc, '~/PM_workflow/.../output/pm.ClassMatrix.loc.csv', row.names = FALSE)

# Remove all bacterial classes but the dominant six used in this study, and agglomerate all "other" into "other" so we can display as a relative abundance taxa plot
pm.ClassMatrix.loc2<-pm.ClassMatrix.loc
keepclasscols<-c("Sample", "Acidimicrobiia", "Actinobacteria", "Alphaproteobacteria", "Bacteroidia", "Gammaproteobacteria", "Oxyphotobacteria")
pm.ClassMatrix.loc2<-pm.ClassMatrix.loc2[,keepclasscols]
dim(pm.ClassMatrix.loc2)

# Convert all other into other variable.
# Remove the top six, so we can sum rows without them.
pm.ClassMatrix.loc.notop6<-pm.ClassMatrix.loc[, ! names(pm.ClassMatrix.loc) %in% keepclasscols, drop = F]

# Sum all of the rows so that way we can collapse into "Other" category across all 20 locations.
Other<-rowSums(pm.ClassMatrix.loc.notop6)
cbindall<-cbind(pm.ClassMatrix.loc2, Other)
cbindall$Sample

row.names(cbindall) <- cbindall$Sample
row.names(cbindall)
cbindall[1] <- NULL

# Convert counts to relative abundance
cbindall.RA=sweep(cbindall,1,rowSums(cbindall),"/")
rowSums(cbindall.RA)

# Write this to a csv and we will use this to format data for graphing a taxa plot
relabund.dom.classes<-as.matrix(cbindall.RA)
write.csv(relabund.dom.classes, "~/pm_workflow/.../output/relabund.dom.classes.csv")

# Call formatted data from relabund.dom.classes.csv ordered by bacterial class from GitHub
plot.data<-'https://github.com/rghannam/pm_workflow/raw/master/data/graphing_inputs/fig2/abundance.plot.formatted.csv'
plot.data<-read.csv(url(plot.data))
colnames(plot.data)
plot.data$xaxis

# Relevel the stack var
plot.data$stack <- factor(plot.data$stack, levels = c('Acidimicrobiia', 'Actinobacteria', 'Alphaproteobacteria', 'Gammaproteobacteria', 'Bacteroidia', 'Oxyphotobacteria', 'Other'))

# Order as Figure 1 is: East West Europe Asia Lakes
plot.data$xaxis <- factor(plot.data$xaxis, levels = c('Keweenaw', 'Green Bay', 'Duluth', 'Singapore', 'Hong Kong', 'Busan', 'Wilhelmshaven', 'Venice', 'Rotterdam', 'Naples', 'Martigues', 'Seattle', 'Oakland', 'Los Angeles', 'Norfolk', 'New York', 'New Orleans', 'Galveston', 'Charleston', 'Baltimore'))

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

marker_style <- list(line = list(width = 0.5,
                                 color = 'rgb(0, 0, 0)'));
t <- list(
  family = "Avenir",
  size = 12,
  color = 'black')

colorsspectral=c('#9E0142', '#F46D43', '#FEE08B','#66C2A5', '#3288BD', '#5E4FA2', "#FFFFFF")

# Plot
plot.Fig2 <-plot.data %>% group_by(stack) %>% arrange(xaxis) %>%
  plot_ly( x = ~value, y = ~xaxis, color= ~stack, type = 'bar', orientation = 'h', mode='marker',
           marker=marker_style, textposition = 'auto',
           colors=colorsspectral,
           line = list(color = "#000000", width = 1)) %>%
  layout(title = "Figure 2", barmode = 'stack',
         xaxis = figxaxis, font=t,
         yaxis = figyaxis)

plot.Fig2



# ---- Figure 3 ----
rm(list=ls()[! ls() %in% c("pm.ps", "pm.counts", "pm.metadata", "pm.taxa", "seed")])
pkgs <- c("ggplot2", "rlang", "Rcpp", "vegan", "csv", "plotly", "plyr", "dplyr", "reshape", "reshape2", "tidyverse") # Load packages
invisible(sapply(pkgs,require, character = TRUE))

c<-subset_taxa(pm.ps, !is.na(Class) & !Class %in% c("", "uncharacterized"))
x<-list()
alltaxa<-list(c, x)
listNames <- c("class")
names(alltaxa) <- listNames

prev.abund.class = apply(X = otu_table(alltaxa$class),
                         MARGIN = ifelse(taxa_are_rows(alltaxa$class), yes = 1, no = 2),
                         FUN = function(x){sum(x > 0)})
prev.abund.class = data.frame(Prevalence = prev.abund.class,
                              TotalAbundance = taxa_sums(alltaxa$class),
                              tax_table(alltaxa$class))
#write.csv(prev.abund.class, '~/pm_workflow/.../output/prev.abund.class.csv', row.names = TRUE)

pm.totalavgprev.class <- plyr::ddply(prev.abund.class, "Class", function(df1){cbind(mean(df1$Prevalence),sum(df1$Prevalence))})
as.data.frame(pm.totalavgprev.class)
names(pm.totalavgprev.class) <- c("Class", "mean", "total")
#write.csv(pm.totalavgprev.class, '~/pm_workflow/.../output/pm.totalavgprev.class', row.names = TRUE)

pm.totalavgprev.class.subset <- subset(pm.totalavgprev.class, pm.totalavgprev.class$total <= 15)
extract.class<-as.data.frame(pm.totalavgprev.class.subset[,1])
names(extract.class) <- ("Class") 
extract.class <- as.list(as.vector(t(extract.class)))

filterclass = (extract.class)
pm.ps.class.sub = subset_taxa(alltaxa$class, !Class %in% filterclass)
pm.ps.class.sub
alltaxa$class

length(get_taxa_unique(alltaxa$class, taxonomic.rank = "Class")) # 143 unique classes before filtering
length(get_taxa_unique(pm.ps.class.sub, taxonomic.rank = "Class")) # 88 unique classes after filtering

pm.ps.class.1218x88 = tax_glom(pm.ps.class.sub, "Class", NArm = TRUE, bad_empty=c(NA, "", " ", "\t"))
pm.ps.class.1218x88

# Here, we need to use our pm.ps.class agglomerated phyloseq object and only keep the dominant classes
pm.ps.class.1218x88

pm.ps.class.melt <- psmelt(pm.ps.class.1218x88)
pm.ps.class.melt$Class <- as.character(pm.ps.class.melt$Class)
pm.ps.class.melt
pm.ps.class.melt <- aggregate(Abundance~Sample+Class, pm.ps.class.melt, FUN=sum)
pm.ClassMatrix1218 <- cast(pm.ps.class.melt, Sample ~ Class)
dim(pm.ClassMatrix1218)

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

pm.ClassMatrix1218.RA # This is usedto format plot.data2

# We will then use all samples (n = 1218) by 6 classes for a total of 7,308 rows of plot data

# This file is a formatted version of pm.ClassMatrix1218.RA object
plot.data2<-'https://github.com/rghannam/pm_workflow/raw/master/data/graphing_inputs/fig3/box_plot_formatted.csv'
plot.data2<-read.csv(url(plot.data2))
colnames(plot.data2)

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
                  tickfont = list(family = 'Avenir',
                                  size = 16,
                                  color = '#000000'))

plot.data2$group <- factor(plot.data2$group, levels = c('Acidimicrobiia', 'Actinobacteria', 'Alphaproteobacteria', 'Gammaproteobacteria', 'Bacteroidia', 'Oxyphotobacteria'))

plot.data2$region <- factor(plot.data2$region, levels = c('East', 'West', 'Europe', 'Asia', 'Lakes'))

Fig3.jitters <- plot_ly(plot.data2, x = ~region, y = ~abundance, color = ~group,marker=marker_style,type = "box",
                        colors=c('#9E0142', '#F46D43', '#FEE08B','#66C2A5', '#3288BD', '#5E4FA2')) %>% 
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
# *Lets find the mean of of each bacterial class at each region region

# Create average abundance vector as an extra filtering parameter to be added to our final dataframe
colnames(plot.data2)
percentage <- prop.table(table(plot.data2$region)) * 100
percentage
getfreq<-cbind(freq=table(plot.data2$region), percentage=percentage)
getfreq<-as.data.frame(getfreq)
getfreq
row.names(getfreq) <- c("East", "West", "Europe", "Asia", "Lakes")
freq <- c("freq")
getfreq<-getfreq[,freq]
getfreq<-as.data.frame(getfreq)
getfreq.labels.vec<-dplyr::pull(getfreq, getfreq)
getfreq.labels.vec

sum.region<-ddply(plot.data2, "region", numcolwise(sum))
sum.region # Find the total samples (n) for each region (1218 total)

# Now we have a vector of the total samples from each region, now lets find total abundance.
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


# This sums the classes by region (abundregion) so we can divide total samples per region
require(data.table) 
DT <- data.table(boxplotX4) 
head(DT)
DT<-DT[ , .(abundregion = sum(abundance)), by = .(group,region)] # This gives us each region per class total abundance
DT

# Then we have to divide these values by vec: sum.region
sum.region
require(tidyverse)
avgabund="avgabund"
DT[ , "avgabund"] <- avgabund
sum.region
DT

# Just to divide these totals by samples (n) so we can get an average relative abundance to compare to plot.Fig3
region_split <- split(DT, DT$region)
reg_names <- c("East", "West", "Europe", "Asia", "Lakes")
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
region_split # This gives us avg rel abund of each class at each region.

# Sanity check
# Sum samples of each region 
# Sum the factor based on unique levels
library(dplyr)
ddply(DT,~region,summarise,number_of_distinct_orders=length(unique(region)))



# ---- Figure 4 ----
rm(list=ls()[! ls() %in% c("pm.counts", "pm.metadata", "pm.taxa", "seed")])

pkgs <- c("caretEnsemble", "caret", "randomForest", "plyr", "tidyverse", "data.table", "csv", "magrittr") # Load packages
invisible(sapply(pkgs,require, character = TRUE))

# Here, call the ASV machine learning models from GitHub
# These were generated through machine_learning.R
ASV.region.model<-'https://github.com/rghannam/pm_workflow/raw/master/data/machine_learning/models/region/ASV.region.model.rds'
ASV.region.model<-readRDS(url(ASV.region.model))
ASV.region.model

ASV.port.model<-'https://github.com/rghannam/pm_workflow/raw/master/data/machine_learning/models/port/ASV.port.model.rds'
ASV.port.model<-readRDS(url(ASV.port.model))
ASV.port.model

# These are counts after removing highly correlated variables to reduce features from 119097 to 6916
ASV.counts.6916vars<-'https://github.com/rghannam/pm_workflow/raw/master/data/counts_taxonomy/pm_counts_noncorvars.rds'
ASV.counts.6916vars<-readRDS(url(ASV.counts.6916vars))


pm.taxa<-'https://github.com/rghannam/pm_workflow/raw/master/data/counts_taxonomy/pm_taxonomy.rds'
pm.taxa<-readRDS(url(pm.taxa))
dim(pm.taxa)

# Now lets read in our taxonomy table and replace values with species.
head(pm.taxa[1:5])
pm.counts<-ASV.counts.6916vars
colnames(pm.counts)<-NULL
pm.taxa_t<-t(pm.taxa)
dim(pm.taxa_t)
colnames(pm.taxa_t) <- NULL
pm.taxa_re_t<-t(pm.taxa_t)
dim(pm.taxa_re_t)

# Phyloseq object creation - note that we're using 6916 features, not 119097 like from earlier figures
counts = otu_table(pm.counts, taxa_are_rows = FALSE) # leave false because we arbitrarily set column names with colnames()->NULL, we can always map back to respective ASV
taxa = tax_table(pm.taxa_re_t)
pm.ps = phyloseq(counts, taxa)
pm.ps
otu_table(pm.ps)[1:5, 1:5]
tax_table(pm.ps)[1:5, 1:5]

# Find shared impportant vars from ASV port/region models

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

# vi = important variable
code <- caret:::varImpDependencies('rf')
code$varImp

regvi<-varImp(ASV.region.model$rf)
regvi
regvi[["importance"]]
regvidf<-as.data.frame(regvi[["importance"]])  # Coerce vi's to df, arranged in order
regvidf

portvi<-caret::varImp(ASV.port.model$rf)
portvi
portvidf<-as.data.frame(portvi[["importance"]]) # Coerce vi's to df, arranged in order
portvidf 
head(portvidf)


# So now lets take our vi dfs and replace the column names with sp (species, or actual taxonomy).
# Take column and convert it to a row in df.
asvs<-as.vector(1:6916)
regvidf$asvordered<-asvs # Appended an ASV ordered column so we can compare the 2 columns: asvordered/Overall
portvidf$asvordered<-asvs

# Subset by overall variable importance >=1
regvidf.ss1 <- regvidf[regvidf$Overall >= "1", ]
portvidf.ss1 <- portvidf[portvidf$Overall >= "1", ]

# Find differences between port/regional models.
sd.reg<-as.data.frame(setdiff(regvidf.ss1$asvordered,portvidf.ss1$asvordered)) # Has the vars that are in region but not in port
sd.port<-as.data.frame(setdiff(portvidf.ss1$asvordered,regvidf.ss1$asvordered)) # Has the vars that are in port but not in region

regvidf.ss1$region <- TRUE
portvidf.ss1$port <- TRUE
portregionmerge <- merge(regvidf.ss1, portvidf.ss1, all=TRUE) # Just ordered the important vars in one df.

library(compare)
library(tidyverse)
regvidf.ss1 %>% remove_rownames() %>% column_to_rownames(var = "asvordered")
portvidf.ss1 %>% remove_rownames() %>% column_to_rownames(var = "asvordered")

common <- intersect(regvidf.ss1$asvordered, portvidf.ss1$asvordered) # Find shared ASVs between port/regional models 
common<-as.data.frame(common)

library(plyr)
#negate_match_df <- function (x, y, on = NULL) 
#{
#  if (is.null(on)) {
#    on <- intersect(names(x), names(y))
#    message("Matching on: ", paste(on, collapse = ", "))
#  }
#  keys <- join.keys(x, y, on)
#  x[!(keys$x %in% keys$y), , drop = FALSE]
#}
#diff <- negate_match_df(sd.port,sd.reg) # same as sd.port

sd1<-as.data.frame(setdiff(regvidf.ss1$asvordered,portvidf.ss1$asvordered)) # Has the ones that are in region but not in port
sd2<-as.data.frame(setdiff(portvidf.ss1$asvordered,regvidf.ss1$asvordered)) # Has the ones that are in port but not in region

# Try collapsing all ASVS into respective classes, e.g. assign a spX value (from pm.ps@tax_table) to each bacterial Class, respectively
# So that we can match this value with the variable importance
tax_table(pm.ps)[1:20, 1:5]

taxadf<-as.data.frame(tax_table(pm.ps))
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
# Collapses all classes and shows which variable of ASV exactly belongs to this class.
# We can use this in linking the shared/different variables from our object: common or sd.port/sd.reg
taxadfcollapse<-setDT(taxadf)[, lapply(.SD, paste, collapse = "; "), by = Class]
head(taxadfcollapse)
taxadfcollapse[,3:4]<-NULL

head(taxadfcollapse)

# This shows the distribution of ASVs collapses by Class
library(magrittr)
res = taxadf %>% 
  melt(id.vars = "Class") %>% 
  dcast(Class ~ variable)

head(res)
colSums(res[,-1]) # Sanity check to see if sums to 6916

# This also helps explain why the model chose certain vars as imp, as the majority of ASVs are within this
# Category of classes, they're the most used just by nature of the model.
# Numbers are formatted as strings but can convert

# e.g this means: bacteroideia 1464/6916*100 accounted for X % of asv features.
# And so, an agglomeration is a feature reduciton method that accumulates all of those ASVs into 1 variable, or composites them
# However - we lose accuracy dependent on the taxonomic resolution we choose to look at.
head(res) # rn is the total asvs for that class.

# Take shared as a vector - and extract class by sp.
names(common)[1]<-"asvs"
spvector<-as.vector(common$asvs)
head(taxadf)

extractedsharedasvs<-taxadf[spvector,"Class"]

# Now, append vector as column
extractedsharedasvs$ASV = spvector

length(unique(extractedsharedasvs$Class)) # There's only 10 unique classes of the 134 here.

astibble<-as.tibble(extractedsharedasvs)
astibble[!duplicated(astibble$Class), ] # Gives the 10 unique, 

# What this shows is that of the 134 ASVS that are shared between the models that have varImp >=1
# Only 10 classes comprise the 134 ASVs - of these 10, the top dominant six from our analysis are inside of this.

extractedsharedasvs$ASV
regvidf.ss1$asvordered
df.r<-regvidf.ss1[regvidf.ss1$asvordered %in% spvector,] # Extracted data for overall var IMP
portvidf.ss1$asvordered
df.p<-portvidf.ss1[portvidf.ss1$asvordered %in% spvector,]
# Based on the common vector - we now have the overall importance for region.

# Now - append to 1 df.
# Rename column so we know if its port of region
colnames(df.r)[1]<-"vimpRegion"
colnames(df.p)[1]<-"vimpPort"

# Extract column from extractedsharedasvs and so ww can append to the asvordered column to this dataframe
Classordered<-as.vector(extractedsharedasvs$Class)
vimpRegion<-as.vector(df.r$vimpRegion)
vimpPort<-as.vector(df.p$vimpPort)
extractedsharedasvs$vimpRegion<-vimpRegion
extractedsharedasvs$vimpPort<-vimpPort

# Now we have a final dataframe with:
# Shared ASVs between model varimp>=1%
# Overall important variables from each model
# Which class the ASV belonged to

library("dplyr")
esasvsgpclass<-as.data.frame(extractedsharedasvs %>% group_by(Class) %>% summarize(count=n()))
esasvsgpclass

# This groups the unique entries and shows proportion of classes per total ASVs
# e.g) out of 134, 126 were in the top 6 classes.

extractedsharedasvs # Has distribution of 134 shared asvs between each models along with vImp region vImp port.
esasvsgpclass # Shows the distribution of the 134 - bins 134 into respective ASV
portregionmerge # Has the merged vimp >1
portvidf.ss1 # 361 vimp >1
regvidf.ss1 # 160 vimp > 1
res # Distribution of 6916 asvs

# Extract values from extractedsharedasvs for plotting important predictor variables.
# Order classes together and then apply a classcategory factor.
vImp.df<-extractedsharedasvs
vImp.df<-with(vImp.df, vImp.df[order(Class),]) # Takes multiple arguments, we just need to arrange by class
class(vImp.df)
vImp.df<-as.data.frame(vImp.df)
vImp.df
nclasses<-as.data.frame(vImp.df %>% group_by(Class) %>% summarize(count=n())) #Get the number of classes
nclasses
levels(vImp.df$Class)[levels(vImp.df$Class)=="Planctomycetacia"] <- "Other"
levels(vImp.df$Class)[levels(vImp.df$Class)=="Rhodothermia"] <- "Other"
levels(vImp.df$Class)[levels(vImp.df$Class)=="SL56_marine_group"] <- "Other"
levels(vImp.df$Class)[levels(vImp.df$Class)=="Verrucomicrobiae"] <- "Other"

# Now lets plot these important predictors.


# Define plotting parameters
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

nclasses
# Lets plot a pie chart of ASV distribution now

category = c('Acidimicrobiia', 'Actinobacteria', 'Alphaproteobacteria', 'Bacteroidia',
             'Gammaproteobacteria', 'Oxyphotobacteria', 'Other') 
percentASV = c("2.98", "8.95", "19.40", "24.62", "19.40", "18.65", "5.94")
legend = c('Acidimicrobiia', 'Actinobacteria', 'Alphaproteobacteria', 'Gammaproteobacteria', 'Bacteroidia', 'Oxyphotobacteria', 'Other')
piedf = data.frame(category, percentASV, legend)
piedf
piedf$category <- factor(piedf$category, levels = c('Acidimicrobiia', 'Actinobacteria', 'Alphaproteobacteria', 'Gammaproteobacteria', 'Bacteroidia', 'Oxyphotobacteria', 'Other'))
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
                    insidetextfont = list(color = '#FFFFFF'),
                    hoverinfo = 'text',
                    text = ~paste('%', percentASV, ' percent'),
                    marker = list(colors = c("#9E0142", '#F46D43', '#FEE08B', '#66C2A5', '#3288BD', '#5E4FA2', '#000000'),
                                  line = list(color = '#FFFFFF', width = 2)),
                    #The 'pull' attribute can also be used to create space between the sectors
                    showlegend = TRUE) %>%
  layout(title = 'Total ASV distribution', legend=l,autosize = T, width = 250, height = 250)
plot.Fig4A



# ---- Figure 5 ----

# This figure was generated through rawgraphs.io using sum.sigs.class5 from differential_abundance.R

# Used only EF > 10 and unique, these are the 24 classes in  Fig5A, and Fig5B was made in adobe illustrator



# ---- Figure 6 ----
rm(list=ls())
seed=81
set.seed(seed)
# These results come from machine_learning.R script, I aggregated the macro-averaged results together into a csv for both port/region. Lets call them.
pkgs <- c("caretEnsemble", "caret", "randomForest", "plyr", "tidyverse", "data.table", "csv", "magrittr") # Load packages
invisible(sapply(pkgs,require, character = TRUE))

# Load in regional and port machine learning model stats by each taxonomic level
region.ml.stats<-'https://github.com/rghannam/pm_workflow/raw/master/data/machine_learning/model_stats/region/region_alltaxa_stats.csv'
region.ml.stats<-read.csv(url(region.ml.stats))
colnames(region.ml.stats)

port.ml.stats<-'https://github.com/rghannam/pm_workflow/raw/master/data/machine_learning/model_stats/port/port_alltaxa_stats.csv'
port.ml.stats<-read.csv(url(port.ml.stats))
colnames(port.ml.stats)

region.ml.stats$ResolutionRegion <- factor(region.ml.stats$ResolutionRegion, levels = c('Phylum', 'Class', 'Order', 'Family', "Genus", "ASV")) # this shows port then region
port.ml.stats$ResolutionPort <- factor(port.ml.stats$ResolutionPort, levels = c('Phylum', 'Class', 'Order', 'Family', "Genus", "ASV")) # this shows port then region

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
                 tickfont = list(family = 'Avenir',
                                 size = 16,
                                 color = '#000000'))

region.ml.stats
port.ml.stats

plot.Fig6.region <- plot_ly(region.ml.stats, x = ~ResolutionRegion, y = ~AccuracyRegion, type = 'scatter', mode ='markers', fillcolor='#F68C5B1A', marker = list(
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
  add_trace(y = ~logLossRegion, type = 'scatter', mode = "markers",
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


plot.Fig6.port <- plot_ly(port.ml.stats, x = ~ResolutionPort, y = ~AccuracyPort, type = 'scatter', mode ='markers', fillcolor='#F68C5B1A', marker = list(
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
  add_trace(y = ~logLossPort, type = 'scatter', mode = "markers",
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



# ---- Figure 7 ----
rm(list=ls())
seed=81
set.seed(seed)
pkgs <- c("caretEnsemble", "caret", "randomForest", "csv") # Load packages
invisible(sapply(pkgs,require, character = TRUE))

# Load regional class model
class.region.ml<-'https://github.com/rghannam/pm_workflow/raw/master/data/machine_learning/models/region/class.region.model.rds'
class.region.ml<-readRDS(url(class.region.ml))
class.region.ml

cl.varImp<-varImp(class.region.ml$rf)
cl.varImp
class.region.ml

x<-c("Acidimicrobiia", "Actinobacteria", "Holophagae", "Oxyphotobacteria", "Gammaproteobacteria",
     "Phycisphaerae", "Thermoplasmata", "Rhodothermia", "Pla3_lineage", "Armatimonadia")
y<-c(100.00, 97.13, 76.83, 69.92, 68.811, 59.83, 44.572, 44.50, 39.52, 39.28)
vImp.df<-as.data.frame(cbind(x, y))
vImp.df[,2] <- as.numeric(as.character(vImp.df[,2])) # Have to use as.char before as.num, factors are stored internally as integers, need to apply factor level labels (actual values) or it will convert to integer codes, non numeric chars such as . ruins this.
class(vImp.df$x) # Factor
class(vImp.df$y)

vImp.df=vImp.df %>%
  arrange(y) %>%
  mutate(x=factor(x,x))

plot.Fig7.varImp = ggplot(vImp.df, aes(x=x, y=y)) +
  geom_segment( aes(x=x, xend=x, y=0, yend=y ), color=ifelse(vImp.df$x %in% c("Acidimicrobiia","Actinobacteria","Oxyphotobacteria", "Gammaproteobacteria"), "#F7756D", "#00BEC4"), size=ifelse(vImp.df$x %in% c("A","D"), .08, .5) ) +
  geom_point( color=ifelse(vImp.df$x %in% c("Acidimicrobiia","Actinobacteria","Oxyphotobacteria", "Gammaproteobacteria" ), "#F7756D", "#00BEC4"), size=ifelse(vImp.df$x %in% c("A","D"), 5, 3) ) +
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

cl.varImp
# Order variable importance by by cl.varImp: consistent with previous graphs: east coast/west coast/europe/asia/lakes
varImp.region<-rbind(c(100.00,37.523,81.57,40.30,45.09),c(59.23,58.629,56.11,74.98,97.13), #Acidimicrobiia,Actinobacteria
                     c(55.27,30.131,61.86,76.83,19.38), c(69.92,30.569,34.04,31.01,42.18), #Holophagae,Oxyphotobacteria
                     c(55.93,68.811,55.31,64.77,57.29), c(59.83,19.756,56.43,57.05,32.21), #Gammaproteobacteria,Phycisphaerae
                     c(37.87,44.572,30.60,36.80,34.53), c(30.43,22.407,40.58,44.50,24.45), #Thermoplasmata,Rhodothermia
                     c(28.72,17.346,39.52,31.81,23.05), c(39.28,24.451,31.81,29.01,19.53)) #Pla3_lineage,Armatimonadia

#re-order so aesthetically they are arranged in same order.
Y<-c("Acidimicrobiia", "Actinobacteria", "Holophagae", "Oxyphotobacteria", "Gammaproteobacteria",
     "Phycisphaerae", "Thermoplasmata", "Rhodothermia", "Pla3_lineage", "Armatimonadia")
Y <- factor(Y, levels = c("Armatimonadia", "Pla3_lineage", "Rhodothermia", "Thermoplasmata", 
                          "Phycisphaerae", "Gammaproteobacteria", "Oxyphotobacteria", "Holophagae", "Actinobacteria", "Acidimicrobiia"))

X<-c("East U.S", "West U.S", "Europe", "Asia", "Lakes")

#plot.Fig7.heatmap <- plot_ly(
#  x = ~X, y=~Y,
#  z = varImp.region, type = "heatmap"
#)
#plot.Fig7.heatmap # don't use, have to change legend and not formatted properly

plot.Fig7.heatmap.small <- plot_ly(
  x = ~X,
  z = varImp.region, color=varImp.region, type = "heatmap",   showticklabels = TRUE, autosize=F, width=250, height=700)

plot.Fig7.heatmap.small # Use this for actual block size



# ---- Figure 8 ----

# topef.oftop10varIMP is data from differential_abundance.R - from sum.sigs.class5

# Create a df of the top 10 important predictor variables from Fig8A and give highest EF for each class for each region

topef.oftop10varIMP<-'https://github.com/rghannam/pm_workflow/raw/master/data/graphing_inputs/fig8/EF_for_top10_predictors.csv'
topef.oftop10varIMP<-read.csv(url(topef.oftop10varIMP))
row.names(topef.oftop10varIMP)<-topef.oftop10varIMP$X
topef.oftop10varIMP[1]<-NULL

topef.oftop10varIMP <-rbind(rep(20,10) , rep(0,10) , topef.oftop10varIMP)
colnames(topef.oftop10varIMP) <- c("1" , "2" , "3" , "4" , "5", "6", "7", "8", "9", "10")

# New colors
colors_border=c("#9d1b44", "#f26c44", "#69c3a5", "#3388bd","#fcdf8a")
colors_in=c("#9d1b4466", "#f26c4466", "#69c3a566", "#3388bd66", "#fcdf8a66") #Add 40% opacity to inside (66)

# Prepare title
mytitle <- c("East U.S", "West U.S", "Europe", "Asia", "Great Lakes")

# Split the screen in 6 parts
par(mar=rep(0.8,4))
par(mfrow=c(2,3))

library(fmsb)
# Loop for each plot
for(i in 1:5){
  
  # Customize the polar chart here
  radarchart(topef.oftop10varIMP[c(1,2,i+2),], axistype=1, 
          
             #custom polygon
             pcol=colors_border[i] , pfcol=colors_in[i] , plwd=2, plty=1 , 
             
             #custom the grid
             cglcol="grey", cglty=1, axislabcol="black", caxislabels=seq(0,20,5), cglwd=0.8,
             
             #custom labels
             vlcex=0.8,
             
             #title
             title=mytitle[i]
  )
}



# ---- Supplemental Figures ----

# ---- Figure S1 ----
rm(list=ls())
seed=81
set.seed(seed)
pkgs <- c("phyloseq", "ggplot2", "rlang", "Rcpp", "vegan", "csv", "plotly", "plyr", "dplyr", "reshape", "reshape2") # Load packages
invisible(sapply(pkgs,require, character = TRUE))

# Refer to Figure 2 on plot.data was generated (line 72)
plot.data2<-'https://github.com/rghannam/pm_workflow/raw/master/data/graphing_inputs/fig3/box_plot_formatted.csv'
plot.data2<-read.csv(url(plot.data2))
colnames(plot.data2)

plot.data2$group <- factor(plot.data2$group, levels = c('Acidimicrobiia', 'Actinobacteria', 'Alphaproteobacteria', 'Gammaproteobacteria', 'Bacteroidia', 'Oxyphotobacteria'))
plot.data2$region <- factor(plot.data2$region, levels = c('East', 'West', 'Europe', 'Asia', 'Lakes'))


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


# Lets use that main data.frame to subset and plot each class and its dispersion
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

# Would need to make 1 for each class and just put three by three in one figure.
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



# ---- Figure S2 ----
rm(list=ls())
seed=81
set.seed(seed)

ASV.counts.6916vars<-'https://github.com/rghannam/pm_workflow/raw/master/data/counts_taxonomy/pm_counts_noncorvars.rds'
ASV.counts.6916vars<-readRDS(url(ASV.counts.6916vars))

pm.counts<-ASV.counts.6916vars

pm.taxa<-'https://github.com/rghannam/pm_workflow/raw/master/data/counts_taxonomy/pm_taxonomy.rds'
pm.taxa<-readRDS(url(pm.taxa))
dim(pm.taxa)

pm.metadata<-'https://github.com/rghannam/pm_workflow/raw/master/data/metadata/pm_metadata.csv'
pm.metadata<-read.csv(url(pm.metadata))
row.names(pm.metadata)<-pm.metadata$sample_name
row.names(pm.metadata)<-gsub("-", "", row.names(pm.metadata)) # Used "-" character some samples so some extra work needs to be done to account for this for downstream analysis
pm.metadata<-pm.metadata[-1]
library(dplyr)
library(tibble)
pm.metadata <- pm.metadata %>% rownames_to_column("sample_name")
row.names(pm.metadata)<-pm.metadata$sample_name

colnames(pm.counts)<-NULL
pm.taxa_t<-t(pm.taxa)
dim(pm.taxa_t)
colnames(pm.taxa_t) <- NULL
pm.taxa_re_t<-t(pm.taxa_t)
dim(pm.taxa_re_t)
dim(pm.metadata)

sampledata = sample_data(pm.metadata)
counts = otu_table(pm.counts, taxa_are_rows = FALSE) # leave false because we arbitrarily set column names with colnames()->NULL, we can always map back to respective ASV
taxa = tax_table(pm.taxa_re_t)
sampledata = sample_data(pm.metadata)
pm.ps.6916 = phyloseq(counts, taxa, sampledata)
pm.ps.6916
otu_table(pm.ps.6916)[1:5, 1:5]


# ANOSIM plot

# Use the 6916x1218 matrix.
dim(pm.ps.6916@sam_data)

abund_dist.6916x1218 <- vegdist(ASV.counts.6916vars)
abund_dist.6916x1218

colnames(pm.ps.6916@sam_data)
abund_ano.Local <- anosim(abund_dist.6916x1218, pm.ps.6916@sam_data$geo_port) # Use geo_port factor for local scales
abund_ano.Local
plot(abund_ano.Local)

abund_ano.Region <- anosim(abund_dist.6916x1218, pm.ps.6916@sam_data$geo_region) # Use geo_region factor for regional scales
abund_ano.Region
plot(abund_ano.Region)


# ---- Figure S3 ----

# Continuation from Figure S2 using the phyloseq object pm.ps.6916

pm.ps.dist.jaccard <- distance(pm.ps.6916,  method = "jaccard")
pm.ps.dist.jaccard

library(ape)

pcoa.lingoes<-pcoa(pm.ps.dist.jaccard, correction="lingoes", rn=NULL)

pcoa = plot_ordination(pm.ps.6916, pcoa.lingoes,
                              type="sample_name", 
                              color ="geo_region") +
  theme_bw() +
  stat_ellipse(geom = "polygon", 
               alpha = 0.2, 
               aes(fill = geo_region))

pcoa

#pcoa[["plot_env"]][["DF"]] # apply axis values to each sample.



# ---- Figure S4 ----

# Refer to Figure 4 - line 490, to create this object (vImp.df)
vImp.df<-extractedsharedasvs
vImp.df<-with(vImp.df, vImp.df[order(Class),]) # Takes multiple arguments, we just need to arrange by class
class(vImp.df)
vImp.df<-as.data.frame(vImp.df)
vImp.df

nclasses<-as.data.frame(vImp.df %>% group_by(Class) %>% summarize(count=n())) #Get the number of classes
nclasses

levels(vImp.df$Class)[levels(vImp.df$Class)=="Planctomycetacia"] <- "Other"
levels(vImp.df$Class)[levels(vImp.df$Class)=="Rhodothermia"] <- "Other"
levels(vImp.df$Class)[levels(vImp.df$Class)=="SL56_marine_group"] <- "Other"
levels(vImp.df$Class)[levels(vImp.df$Class)=="Verrucomicrobiae"] <- "Other"


# Refer to Figure 4 (line X) for vImp.df object, we will base our interleaved histograms from this object

# Lets plot all of the 134 important vars, respectively ie) show 7 little plots, with 4 acido etc...
vimpregion<-as.vector(vImp.df$vimpRegion)
vimpport<-as.vector(vImp.df$vimpPort)
classes<-as.vector(vImp.df$Class)
color_pallete<-c("#9D1B44", '#F46D43', '#FEE08B', '#66C2A5', '#3288BD', '#5E4FA2', '#000000')
vimpregion
vimpport

df<-data.frame(vimpregion,vimpport, classes )

# Interleaved histograms

vimpport<-ggplot(df, aes(x=vimpport, color=classes)) + 
  geom_histogram(fill="white", position="dodge2")+
  scale_color_manual(labels = c('Acidimicrobiia', 'Actinobacteria', 'Alphaproteobacteria', 'Gammaproteobacteria', 'Bacteroidia', 'Oxyphotobacteria', 'Other'), 
                     values = color_pallete) +
  theme(legend.position="top") +   theme_minimal() + geom_freqpoly()
vimpport

vimpregion<-ggplot(df, aes(x=vimpregion, color=classes)) + 
  geom_histogram(fill="white", position="dodge2")+
  scale_color_manual(labels = c('Acidimicrobiia', 'Actinobacteria', 'Alphaproteobacteria', 'Gammaproteobacteria', 'Bacteroidia', 'Oxyphotobacteria', 'Other'), 
                     values = color_pallete) +
  theme(legend.position="top") +   theme_minimal() + geom_freqpoly()
vimpregion


# ---- Figure S5/S6 ----
rm(list=ls())
seed=81
set.seed(seed)
ASV.counts.6916vars<-'https://github.com/rghannam/pm_workflow/raw/master/data/counts_taxonomy/pm_counts_noncorvars.rds'
ASV.counts.6916vars<-readRDS(url(ASV.counts.6916vars))

pm.counts<-ASV.counts.6916vars

pm.taxa<-'https://github.com/rghannam/pm_workflow/raw/master/data/counts_taxonomy/pm_taxonomy.rds'
pm.taxa<-readRDS(url(pm.taxa))
dim(pm.taxa)

pm.metadata<-'https://github.com/rghannam/pm_workflow/raw/master/data/metadata/pm_metadata.csv'
pm.metadata<-read.csv(url(pm.metadata))
row.names(pm.metadata)<-pm.metadata$sample_name
row.names(pm.metadata)<-gsub("-", "", row.names(pm.metadata))
pm.metadata<-pm.metadata[-1]
library(dplyr)
library(tibble)
pm.metadata <- pm.metadata %>% rownames_to_column("sample_name")
row.names(pm.metadata)<-pm.metadata$sample_name

# Import vector containing only those samples for which we have environmental conditions collected
env_samples_vec<-read.csv('https://github.com/rghannam/pm_workflow/raw/master/data/graphing_inputs/fig5S-6S/env_samples_vec.csv', header=FALSE)

row.names(pm.counts)[1:200]
rows1218<-row.names(pm.counts)


pm.metadata <- pm.metadata[rownames(pm.metadata) %in% rows1218, ] # Match metadata samples with those in count table
pm.metadata2 <- pm.metadata[rownames(pm.metadata) %in% env_samples_vec$V1, ]
pm.counts2 <- pm.counts[rownames(pm.counts) %in% env_samples_vec$V1, ]
row.names(pm.metadata2)[1:200] # Match
row.names(pm.counts2)[1:200] # Match

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
row.names(pm.counts2)[1:1176]
row.names(pm.metadata3)[1:1176]
head(pm.counts2)[1:10]
head(pm.metadata3)[1:8]
colnames(pm.counts2)<-NULL

pm.taxa_t<-t(pm.taxa)
dim(pm.taxa_t)
colnames(pm.taxa_t) <- NULL
pm.taxa_re_t<-t(pm.taxa_t)
dim(pm.taxa_re_t)
dim(pm.metadata3) # [1] 1514   16
pm.metadata3[, 1] <- as.factor(pm.metadata3[, 1])

sampledata = sample_data(pm.metadata3)
counts = otu_table(pm.counts2, taxa_are_rows = FALSE) # leave false because we arbitrarily set column names with colnames()->NULL, we can always map back to respective ASV
taxa = tax_table(pm.taxa_re_t)
sampledata = sample_data(pm.metadata)
pm.ps.env = phyloseq(counts, taxa, sampledata)
pm.ps.env
otu_table(pm.ps.env)[1:5, 1:5]

# Extract abundance matrix from the phyloseq object
pm.ps.env_abund = as(otu_table(pm.ps.env), "matrix")
# Transpose if necessary
#if(taxa_are_rows(pm.ps.env)){pm.ps.env_abund <- t(pm.ps.env_abund)}
# Coerce to data.frame
pm.ps.env_abund.df = as.data.frame(pm.ps.env_abund)

head(pm.ps.env_abund.df[1:5]) # named by taxa table
rm(list=ls()[! ls() %in% c("pm.ps.env_abund.df","seed", "pm.ps.env")]) # Clear memory
pm.ps.env

# Convert to DESeq2 and normalize counts
pm.ps.env.deseq = phyloseq_to_deseq2(pm.ps.env, ~ geo_region)
pm.ps.env.deseq

# Calculate geometric means prior to estimate size factors

gm_mean = function(x, na.rm=TRUE){
  exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))
}

library(DESeq2)

# To just reproduce this plot, skip this step if you lack the processing power, and just use pm.ps.env
geoMeans = apply(counts(pm.ps.env.deseq), 1, gm_mean)
pm.ps.env.deseqsf = estimateSizeFactors(pm.ps.env.deseq, geoMeans = geoMeans)
pm.ps.env.deseq.norm = DESeq(pm.ps.env.deseqsf, fitType="local")

ncts <- counts(pm.ps.env.deseq.norm, normalized=TRUE)
ncts<-t(ncts)
head(ncts[1:5])

saveRDS(ncts,"~/PMP_workflow/Figures/Fig5-6S/outputs//1176/ncts1176.rds" )
saveRDS(pm.ps.env.deseq.norm,"~/PMP_workflow/Figures/Fig5-6S/outputs/1176/pm.ps.env.deseq.norm.rds")

# Just call the normalized object from GitHub
ncts<-'https://github.com/rghannam/pm_workflow/raw/master/data/graphing_inputs/fig5S-6S/ncts1176.rds'
ncts<-readRDS(url(ncts))

ncts2<-ncts

pm.taxa<-'https://github.com/rghannam/pm_workflow/raw/master/data/counts_taxonomy/pm_taxonomy.rds'
pm.taxa<-readRDS(url(pm.taxa))
dim(pm.taxa)

pm.metadata<-'https://github.com/rghannam/pm_workflow/raw/master/data/metadata/pm_metadata.csv'
pm.metadata<-read.csv(url(pm.metadata))
row.names(pm.metadata)<-pm.metadata$sample_name
row.names(pm.metadata)<-gsub("-", "", row.names(pm.metadata))
pm.metadata<-pm.metadata[-1]
library(dplyr)
library(tibble)
pm.metadata <- pm.metadata %>% rownames_to_column("sample_name")
row.names(pm.metadata)<-pm.metadata$sample_name

row.names(pm.counts)[1:200]
rows1176<-row.names(ncts2)
env_samples_vec<-read.csv('https://github.com/rghannam/pm_workflow/raw/master/data/graphing_inputs/fig5S-6S/env_samples_vec.csv', header=FALSE)

pm.metadata <- pm.metadata[rownames(pm.metadata) %in% rows1176, ] # match count table
pm.metadata2 <- pm.metadata[rownames(pm.metadata) %in% env_samples_vec$V1, ]
ncts3 <- ncts2[rownames(ncts2) %in% env_samples_vec$V1, ]
row.names(pm.metadata2)[1:200]
row.names(ncts2)[1:200]

colnames(pm.metadata2)
keeps <- c("sample_name", "conduc","salinity", "total_dissolved_solids",
           "surf_temp", "ph", "diss_oxygen", "geo_region")

pm.metadata3<-pm.metadata2[(names(pm.metadata2) %in% keeps)]

head(pm.metadata3)[1:7]
pm.metadata3$conduc<-gsub(" uS/cm", "", as.character(pm.metadata3$conduc) , n)
pm.metadata3$salinity<-gsub(" psu", "", as.character(pm.metadata3$salinity) , n)
pm.metadata3$total_dissolved_solids<-gsub(" mg/L", "", as.character(pm.metadata3$total_dissolved_solids) , n)
pm.metadata3$surf_temp<-gsub(" deg C", "", as.character(pm.metadata3$surf_temp) , n)
pm.metadata3$diss_oxygen<-gsub(" mg/L", "", as.character(pm.metadata3$diss_oxygen) , n)
pm.metadata3$ph<-as.character(pm.metadata3$ph)

library(plyr)
pm.metadata3$geo_region <- revalue(pm.metadata3$geo_region, c("West Coast U.S"="West"))
pm.metadata3$geo_region <- revalue(pm.metadata3$geo_region, c("East Coast U.S"="East"))
pm.metadata3$geo_region <- revalue(pm.metadata3$geo_region, c("Great Lakes"="Lakes"))

# Rename columns
colnames(pm.metadata3)[colnames(pm.metadata3)=="conduc"] <- "Conductivity"
colnames(pm.metadata3)[colnames(pm.metadata3)=="salinity"] <- "Salinity"
colnames(pm.metadata3)[colnames(pm.metadata3)=="total_dissolved_solids"] <- "TDS"
colnames(pm.metadata3)[colnames(pm.metadata3)=="surf_temp"] <- "Temp"
colnames(pm.metadata3)[colnames(pm.metadata3)=="ph"] <- "pH"
colnames(pm.metadata3)[colnames(pm.metadata3)=="diss_oxygen"] <- "ODO"

pm.metadata3$Conductivity<-as.numeric(pm.metadata3$Conductivity)
class(pm.metadata3$Conductivity)
pm.metadata3$Salinity<-as.numeric(pm.metadata3$Salinity)
pm.metadata3$TDS<-as.numeric(pm.metadata3$TDS)
pm.metadata3$Temp<-as.numeric(pm.metadata3$Temp)
pm.metadata3$pH<-as.numeric(pm.metadata3$pH)
pm.metadata3$ODO<-as.numeric(pm.metadata3$ODO)

pm.metadata3[1]<-NULL
head(pm.metadata3)[1:7]

levels(pm.metadata3$geo_region)
row.names(ncts3)[1:1176]
row.names(pm.metadata3)[1:1176]
head(ncts3)[1:10]
head(pm.metadata3)[1:7]
colnames(ncts3)<-NULL

# Recreate phyloseq object with normalized data
pm.taxa_t<-t(pm.taxa)
dim(pm.taxa_t)
colnames(pm.taxa_t) <- NULL
pm.taxa_re_t<-t(pm.taxa_t)
dim(pm.taxa_re_t)
#pm.metadata3[, 1] <- as.factor(pm.metadata3[, 1])

sampledata = sample_data(pm.metadata3)
counts = otu_table(ncts3, taxa_are_rows = FALSE)
taxa = tax_table(pm.taxa_re_t)
sampledata = sample_data(pm.metadata3)
pm.ps.env.normed.justenvironment = phyloseq(counts, taxa, sampledata)
pm.ps.env.normed.justenvironment
otu_table(pm.ps.env.normed.justenvironment)[1:5, 1:5]

head(pm.ps.env.normed.justenvironment@sam_data)[1:5]
rownames(pm.ps.env.normed.justenvironment@otu_table)[1:5]

# Use package microbiomeSeq for this
seed=81
set.seed(seed)
library(phyloseq)
library(ggplot2)
library(microbiomeSeq)

colnames(pm.ps.env.normed.environment.sample.data.only@sam_data)[colnames(pm.ps.env.normed.environment.sample.data.only@sam_data)=="geo_region"]<-"region"

physeq <- taxa_level(pm.ps.env.normed.environment.sample.data.only, "Class")
physeq
#names(pm.ps.env.normed.environment.sample.data.only@sam_data) <- letters[1:7];

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


# ---- Figure S5 ----
phy.arranged<-arrange.vars(physeq@sam_data, c("Conductivity"=1, "ODO"=2, "pH"=3,
                                              "Salinity"=4,"TDS"=5,"Temp"=6))


plot.FigS5 <- plot_anova_env(phy.arranged, grouping_column = "region", pValueCutoff = 0.001, 
                             select.variables = c("Conductivity", "ODO", 
                                                  "pH", "Salinity", "TDS", "Temp"))


print(plot.FigS5)

# ---- Figure S6 ----
env.taxa.cor <- taxa.env.correlation(physeq, grouping_column = "region", method = "pearson", 
                                     pvalue.threshold = 0.001, padjust.method = "BH", adjustment = 5, num.taxa = 43, 
                                     select.variables = NULL)


plot.FigS6 <- plot_taxa_env(env.taxa.cor)
print(plot.FigS6)

# ---- Identify significant environmental variables  ----
rm(list=ls())
seed=81
set.seed(seed)
ncts<-'https://github.com/rghannam/pm_workflow/raw/master/data/graphing_inputs/fig5S-6S/ncts1176.rds'
ncts<-readRDS(url(ncts))
ncts2<-ncts

pm.taxa<-'https://github.com/rghannam/pm_workflow/raw/master/data/counts_taxonomy/pm_taxonomy.rds'
pm.taxa<-readRDS(url(pm.taxa))
dim(pm.taxa)

pm.metadata<-'https://github.com/rghannam/pm_workflow/raw/master/data/metadata/pm_metadata.csv'
pm.metadata<-read.csv(url(pm.metadata))
row.names(pm.metadata)<-pm.metadata$sample_name
row.names(pm.metadata)<-gsub("-", "", row.names(pm.metadata))
pm.metadata<-pm.metadata[-1]
library(plyr)
library(dplyr)
library(tibble)
pm.metadata <- pm.metadata %>% rownames_to_column("sample_name")
row.names(pm.metadata)<-pm.metadata$sample_name

row.names(ncts2)[1:200]
rows1176<-row.names(ncts2)
env_samples_vec<-read.csv('https://github.com/rghannam/pm_workflow/raw/master/data/graphing_inputs/fig5S-6S/env_samples_vec.csv', header=FALSE)

pm.metadata <- pm.metadata[rownames(pm.metadata) %in% rows1176, ]
pm.metadata2 <- pm.metadata[rownames(pm.metadata) %in% env_samples_vec$V1, ]
ncts3 <- ncts2[rownames(ncts2) %in% env_samples_vec$V1, ]
row.names(pm.metadata2)[1:200]
row.names(ncts2)[1:200]

colnames(pm.metadata2)
keeps <- c("sample_name", "conduc","salinity", "total_dissolved_solids",
           "surf_temp", "ph", "diss_oxygen", "geo_region")

pm.metadata3<-pm.metadata2[(names(pm.metadata2) %in% keeps)]

head(pm.metadata3)[1:7]
pm.metadata3$conduc<-gsub(" uS/cm", "", as.character(pm.metadata3$conduc) , n)
pm.metadata3$salinity<-gsub(" psu", "", as.character(pm.metadata3$salinity) , n)
pm.metadata3$total_dissolved_solids<-gsub(" mg/L", "", as.character(pm.metadata3$total_dissolved_solids) , n)
pm.metadata3$surf_temp<-gsub(" deg C", "", as.character(pm.metadata3$surf_temp) , n)
pm.metadata3$diss_oxygen<-gsub(" mg/L", "", as.character(pm.metadata3$diss_oxygen) , n)
pm.metadata3$ph<-as.character(pm.metadata3$ph)

library(plyr)
pm.metadata3$geo_region <- revalue(pm.metadata3$geo_region, c("West Coast U.S"="West"))
pm.metadata3$geo_region <- revalue(pm.metadata3$geo_region, c("East Coast U.S"="East"))
pm.metadata3$geo_region <- revalue(pm.metadata3$geo_region, c("Great Lakes"="Lakes"))

colnames(pm.metadata3)[colnames(pm.metadata3)=="conduc"] <- "Conductivity"
colnames(pm.metadata3)[colnames(pm.metadata3)=="salinity"] <- "Salinity"
colnames(pm.metadata3)[colnames(pm.metadata3)=="total_dissolved_solids"] <- "TDS"
colnames(pm.metadata3)[colnames(pm.metadata3)=="surf_temp"] <- "Temp"
colnames(pm.metadata3)[colnames(pm.metadata3)=="ph"] <- "pH"
colnames(pm.metadata3)[colnames(pm.metadata3)=="diss_oxygen"] <- "ODO"

pm.metadata3$Conductivity<-as.numeric(pm.metadata3$Conductivity)
class(pm.metadata3$Conductivity)
pm.metadata3$Salinity<-as.numeric(pm.metadata3$Salinity)
pm.metadata3$TDS<-as.numeric(pm.metadata3$TDS)
pm.metadata3$Temp<-as.numeric(pm.metadata3$Temp)
pm.metadata3$pH<-as.numeric(pm.metadata3$pH)
pm.metadata3$ODO<-as.numeric(pm.metadata3$ODO)

pm.metadata3[1]<-NULL
head(pm.metadata3)[1:7]

levels(pm.metadata3$geo_region)
row.names(ncts3)[1:1176]
row.names(pm.metadata3)[1:1176]
head(ncts3)[1:10]
head(pm.metadata3)[1:7]

pm.taxa_t<-t(pm.taxa)
dim(pm.taxa_t)
colnames(pm.taxa_t) <- NULL
pm.taxa_re_t<-t(pm.taxa_t)
dim(pm.taxa_re_t)
#pm.metadata3[, 1] <- as.factor(pm.metadata3[, 1])

sampledata = sample_data(pm.metadata3)
counts = otu_table(ncts3, taxa_are_rows = FALSE) # leave false because we arbitrarily set column names with colnames()->NULL, we can always map back to respective ASV
taxa = tax_table(pm.taxa_re_t)
sampledata = sample_data(pm.metadata3)
pm.ps.env.normed.justenvironment = phyloseq(counts, taxa, sampledata)
pm.ps.env.normed.justenvironment
otu_table(pm.ps.env.normed.justenvironment)[1:5, 1:5]

head(pm.ps.env.normed.justenvironment@sam_data)[1:5]
rownames(pm.ps.env.normed.justenvironment@otu_table)[1:5]


# Make a count matrix agglomerated to taxonomic Class
pm.ps.env.normed.class = tax_glom(pm.ps.env.normed.justenvironment, "Class", NArm = TRUE, bad_empty=c(NA, "", " ", "\t"))
pm.ps.env.normed.class
library(reshape)
library(reshape2)
pm.ps.env.normed.class.melt <- psmelt(pm.ps.env.normed.class)
pm.ps.env.normed.class.melt$Class <- as.character(pm.ps.env.normed.class.melt$Class)
pm.ps.env.normed.class.melt
pm.ps.env.normed.class.melt <- aggregate(Abundance~Sample+Class, pm.ps.env.normed.class.melt, FUN=sum)
pm.ps.env.normed.class.df <- cast(pm.ps.env.normed.class.melt, Sample ~ Class)
dim(pm.ps.env.normed.class.df)
head(pm.ps.env.normed.class.df)
pm.ps.env.normed.class.df2<-pm.ps.env.normed.class.df

# Use pm_metadata3 for metadata and pm.ps.env.normed.class.df2 for count table
# Lets format it properly.
env.counts<-pm.ps.env.normed.class.df2
env.metadata<-pm.metadata3

# Row names to column in env.metadata
env.metadata$sample_name<-rownames(env.metadata)

# Get a vector of row names to arrange by.
rownames(env.counts)<-env.counts[,1] # put Sample column to rownames

samples.ordered<-rownames(env.counts)
samples.ordered2<-as.data.frame(samples.ordered)
env.metadata<-env.metadata[match(samples.ordered, env.metadata$sample_name),] # Arranged rows in metadata by the rows in count table

# Remove unwanted columns
colnames(env.metadata)
env.metadata[8]<-NULL
env.metadata[7]<-NULL
colnames(env.counts)
env.counts[1]<-NULL

# Just a check to ensure that the samples in meta_table are in the same order as in abund_table
env.metadata2<-env.metadata[rownames(env.counts),]

# Filter out any samples taxas that have zero entries 
env.counts2<-subset(env.counts,rowSums(env.counts)!=0)
env.counts2

dim(env.counts2)
dim(env.metadata2)
head(env.counts2)[1:3]
head(env.metadata2)[1:6]

# Can convert to relative frequencies but skip this
#abund_table<-abund_table/rowSums(abund_table)

# Use PERMANOVA to find significant environmental variables
library(vegan)
table.adonis <- adonis(env.counts2 ~ ., data=env.metadata2, method = "bray")
table.adonis


# ---- Figure S7 ----

# Calculate bray curtis distance matrix just to test statistics on regions 

# Use pm.ps.env.normed.justenvironment - normalized phyloseq object
distanceMethodList
pmps_bray <- phyloseq::distance(pm.ps.env.normed.justenvironment, method = "bray")
pmps_bray

# Make a data frame from the sample_data
sampledf <- data.frame(sample_data(pm.ps.env.normed.justenvironment))

# Adonis test # do the 5 regions we collected samples from have different centroids?
adregion<-adonis(pmps_bray ~ geo_region, data = sampledf) # Using the geo_region var from metadata
adregion

beta <- betadisper(pmps_bray, sampledf$geo_region)
beta
plot(beta)
permutest(beta)

anova(beta)

pm.ps.env.normed.justenvironment@sam_data
pm.ps.env.normed.justenvironment@sam_data$geo_region <- factor(pm.ps.env.normed.justenvironment@sam_data$geo_region, levels = c('East', 'West', 'Europe', 'Asia', 'Lakes'))

# Constrained Analysis of Principal Coordinates (CAP)
colnames(pm.ps.env.normed.justenvironment@sam_data)
# Remove data points with missing metadata
bray_not_na <- pm.ps.env.normed.justenvironment %>%
  subset_samples(
    !is.na(Conductivity) & 
      !is.na(Salinity) &
      !is.na(TDS) & 
      !is.na(Temp) & 
      !is.na(pH) & 
      !is.na(ODO)
  )
bray_not_na # Contains all environmental variables

# CAP ordinate
cap_ord <- ordinate(
  physeq = pm.ps.env.normed.justenvironment, 
  method = "CAP",
  distance = pmps_bray,
  formula = ~ Conductivity + Salinity + TDS + Temp + pH + ODO
)
cap_ord
# CAP plot
cap_plot <- plot_ordination(
  physeq = pm.ps.env.normed.justenvironment, 
  ordination = cap_ord, 
  color = "geo_region", 
  axes = c(1,2)
) + 
  #aes(shape = Station) + 
  #colorsspectral=c('#9E0142', '#F46D43', '#FEE08B','#66C2A5', '#3288BD', '#5E4FA2')
  geom_point(aes(colour = geo_region), alpha = 1, size = 4) + 
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

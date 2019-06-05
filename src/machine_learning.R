# Machine learning models across all levels of taxonomic resolution

library("rlang")
library("Rcpp")
library("csv")
library("plyr")
library("dplyr")
require('caretEnsemble')
require('randomForest')
require('caret')
library('ModelMetrics')

# Define global train control / metric / modeling type
seed=81
set.seed(seed)
trctrlbase<- trainControl(method="repeatedcv", 
                          number=10, 
                          repeats=3, 
                          search="random", # No search parameters selected for global train control
                          summaryFunction = multiClassSummary, # Multiclass summary function, we're using 20 labels for port and 5 for region
                          classProbs=TRUE,
                          preProc=c("center","scale")) # Remove mean value of each feature and divide non constant features by standard deviation
metric <- "Accuracy"

modtypesbase <- list(rf       = caretModelSpec(method="rf", verbose=FALSE))

# ---- Model local locations (20 ports) ----
# Load in port/local labels
localInd<-read.csv('https://github.com/rghannam/pm_workflow/raw/master/data/machine_learning/indexing/local_labels.csv', col.names=c("localInd"), header=FALSE)
levels(localInd$localInd)

# Load in matrices from normalize_counts.R

ncts<-'https://github.com/rghannam/pm_workflow/raw/master/data/counts_taxonomy/taxa_counts_list_normed.rds'
ncts<-readRDS(url(ncts))
dim(ncts[[1]])

# Data manipulations on our normalized count matrices
# Convert elements of this list to data.frames
for (i in 1:length(ncts)){
  ncts[[i]] <- as.data.frame(ncts[[i]])
}

# Append local labels to each element of the list
ncts<-lapply(ncts, function(x) 
  cbind(x, localInd = localInd$localInd))
class(ncts[[1]]$localInd)

# Randomize the observations (samples) across all levels of taxonomic resolution
for (i in 1:length(ncts)){
  ncts[[i]] <- ncts[[i]][sample(nrow(ncts[[i]])),]
}

# Sanity checks
row.names(ncts[[2]])[1:5]
ncts[[2]]$localInd[1:5]


#  Retrieve class frequencies to account for class imbalance and to see if we need to weight some locations
percentage <- prop.table(table(ncts[[1]]$localInd)) * 100
cbind(freq=table(ncts[[1]]$localInd), percentage=percentage)

seed=81
set.seed(seed)

# Models
# Phylum
phylum.port.model <-caretList(localInd~., 
                                    data=ncts[[1]], ntree=2000, importance=TRUE,
                                    trControl = trctrlbase, 
                                    tuneList = modtypesbase)

phylum.port.model
phylum.port.model[["rf"]][["finalModel"]]

# Class
class.port.model <-caretList(localInd~., 
                              data=ncts[[2]], ntree=2000, importance=TRUE,
                              trControl = trctrlbase, 
                              tuneList = modtypesbase)

class.port.model
class.port.model[["rf"]][["finalModel"]]

# Order
order.port.model <-caretList(localInd~., 
                             data=ncts[[3]], ntree=2000, importance=TRUE,
                             trControl = trctrlbase, 
                             tuneList = modtypesbase)

order.port.model
order.port.model[["rf"]][["finalModel"]]

# Family
family.port.model <-caretList(localInd~., 
                             data=ncts[[4]], ntree=2000, importance=TRUE,
                             trControl = trctrlbase, 
                             tuneList = modtypesbase)

family.port.model
family.port.model[["rf"]][["finalModel"]]

# Genus
genus.port.model <-caretList(localInd~., 
                              data=ncts[[5]], ntree=2000, importance=TRUE,
                              trControl = trctrlbase, 
                              tuneList = modtypesbase)

genus.port.model
genus.port.model[["rf"]][["finalModel"]]

# Script not finished....


# ---- Model regional locations (5 regions) ----
# Load in region labels
regionInd<-read.csv('https://github.com/rghannam/pm_workflow/raw/master/data/machine_learning/indexing/region_labels.csv', col.names=c("regionInd"), header=FALSE)
dim(regionInd)
levels(regionInd$regionInd)

# Apply regional class labels
ncts<-'https://github.com/rghannam/pm_workflow/raw/master/data/counts_taxonomy/taxa_counts_list_normed.rds'
ncts<-readRDS(url(ncts))
dim(ncts[[1]])

rm(list=ls()[! ls() %in% c("ncts", "regionInd", "localInd", "modtypesbase", "trctrlbase", "seed")])

# Data manipulations on our normalized count matrices
# Convert elements of this list to data.frames
for (i in 1:length(ncts)){
  ncts[[i]] <- as.data.frame(ncts[[i]])
}

# Append the labels to each element of the list
ncts<-lapply(ncts, function(x) 
  cbind(x, regionInd = regionInd$regionInd))
class(ncts[[1]]$regionInd)

# Randomize the observations (samples) across all levels of taxonomic resolution
for (i in 1:length(ncts)){
  ncts[[i]] <- ncts[[i]][sample(nrow(ncts[[i]])),]
}

# Sanity checks
row.names(ncts[[2]])[1:5]
ncts[[2]]$regionInd[1:5]


#  Retrieve class frequencies to account for class imbalance and to see if we need to weight some locations
percentage <- prop.table(table(ncts[[1]]$regionInd)) * 100
cbind(freq=table(ncts[[1]]$regionInd), percentage=percentage)

seed=81
set.seed(seed)

# Models
# Phylum
phylum.region.model <-caretList(regionInd~., 
                              data=ncts[[1]], ntree=2000, importance=TRUE,
                              trControl = trctrlbase, 
                              tuneList = modtypesbase)

phylum.region.model
phylum.region.model[["rf"]][["finalModel"]]

# Class
class.region.model <-caretList(regionInd~., 
                                data=ncts[[2]], ntree=2000, importance=TRUE,
                                trControl = trctrlbase, 
                                tuneList = modtypesbase)

class.region.model
class.region.model[["rf"]][["finalModel"]]

# Order
order.region.model <-caretList(regionInd~., 
                               data=ncts[[3]], ntree=2000, importance=TRUE,
                               trControl = trctrlbase, 
                               tuneList = modtypesbase)

order.region.model
order.region.model[["rf"]][["finalModel"]]

# Family
family.region.model <-caretList(regionInd~., 
                               data=ncts[[4]], ntree=2000, importance=TRUE,
                               trControl = trctrlbase, 
                               tuneList = modtypesbase)

family.region.model
family.region.model[["rf"]][["finalModel"]]

# Genus
genus.region.model <-caretList(regionInd~., 
                                data=ncts[[5]], ntree=2000, importance=TRUE,
                                trControl = trctrlbase, 
                                tuneList = modtypesbase)

genus.region.model
genus.region.model[["rf"]][["finalModel"]]


# ---- ASV models ----
rm(list=ls()[! ls() %in% c("modtypesbase", "trctrlbase", "seed")])

# Load in ASV count matrix
ASV.counts.6916vars<-'https://github.com/rghannam/pm_workflow/raw/master/data/counts_taxonomy/pm_counts_noncorvars.rds'
ASV.counts.6916vars<-readRDS(url(ASV.counts.6916vars))

# Load in local ASV labels
localIndASV<-read.csv('https://github.com/rghannam/pm_workflow/raw/master/data/machine_learning/indexing/local_labelsASV.csv', col.names=c("localIndASV"), header=FALSE)
levels(localIndASV$localIndASV)

# Load in regional ASV labels
regionIndASV<-read.csv('https://github.com/rghannam/pm_workflow/raw/master/data/machine_learning/indexing/region_labelsASV.csv', col.names=c("regionIndASV"), header=FALSE)
levels(regionIndASV$regionIndASV)

# Apply local vector and randomize matrix
ASV.local<-cbind(ASV.counts.6916vars, localIndASV)
ASV.local <- ASV.local[sample(nrow(ASV.local)),]
row.names(ASV.local)[1:10]
ASV.local$localIndASV[1:10]

# Apply regional vector and randomize matrix
ASV.region<-cbind(ASV.counts.6916vars, regionIndASV)
ASV.region <- ASV.region[sample(nrow(ASV.region)),]
row.names(ASV.region)[1:10]
ASV.region$regionIndASV[1:10]

# ASV local model
ASV.port.model <-caretList(localIndASV~., 
                                data=ASV.local, ntree=2000, importance=TRUE,
                                trControl = trctrlbase, 
                                tuneList = modtypesbase)


# ASV regional model
ASV.region.model <-caretList(regionIndASV~., 
                           data=ASV.region, ntree=2000, importance=TRUE,
                           trControl = trctrlbase, 
                           tuneList = modtypesbase)

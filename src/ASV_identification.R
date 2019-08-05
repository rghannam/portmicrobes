# Identify Amplicon Sequence Variants (ASVs)
# https://benjjneb.github.io/dada2/

# Sequencing files are on accessions PRJNA542890 / PRJNA542685

# Reproducible workflow for run-per-run DNA sequencing read processing.

# ==============================================================
#   Run 1 - pertains to sequencing run 1/3
#   refer to metadata to find which samples belong to which run
# ==============================================================

rm(list=ls())
seed=81
set.seed(seed)

setwd("~/pm_workflow") # CHANGE ME TO THE PROPER DIRECTORY (where you want this run's outputs to go)
library(dada2); packageVersion("dada2")

# The only thing you need to change is this 'path' object
path <- "~/pm_workflow/data/sequencing/inputs-fastqs/run1" # CHANGE ME TO THE PROPER DIRECTORY (containing the fastq files)
list.files(path)

# Forward and reverse fastq filenames have format: SAMPLENAME_R1_001.fastq and SAMPLENAME_R2_001.fastq
fnFs <- sort(list.files(path, pattern="_R1_001.fastq", full.names = TRUE))
fnRs <- sort(list.files(path, pattern="_R2_001.fastq", full.names = TRUE))
# Extract sample names, assuming filenames have format: SAMPLENAME_XXX.fastq
sample.names <- sapply(strsplit(basename(fnFs), "_"), `[`, 1)

# Place filtered files in filtered/ subdirectory
filtFs <- file.path(path, "filtered", paste0(sample.names, "_F_filt.fastq.gz"))
filtRs <- file.path(path, "filtered", paste0(sample.names, "_R_filt.fastq.gz"))

out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, trimLeft=c(19,20), truncLen=c(280,190),
                     maxN=0, maxEE=c(2,5), truncQ=2, rm.phix=TRUE,
                     compress=TRUE, multithread=TRUE)

keep <- out[,"reads.out"] > 50 # Keep only samples that have 50 or more reads after filtering (out)
filtFs1 <- file.path(filtFs)[keep]
filtRs1 <- file.path(filtRs)[keep]

sample.namesF <- sapply(strsplit(basename(filtFs1), "_"), `[`, 1) # Assumes filename = samplename_XXX.fastq.gz
sample.namesR <- sapply(strsplit(basename(filtRs1), "_"), `[`, 1) # Assumes filename = samplename_XXX.fastq.gz
if(!identical(sample.namesF, sample.namesR)) stop("Forward and reverse files do not match.") # Check if F/R match.

sample.names.kept <- sapply(strsplit(basename(filtFs1), "_"), `[`, 1)
names(filtFs1) <- sample.names.kept
names(filtRs1) <- sample.names.kept
set.seed(seed)

# Learn forward error rates from kept > 50 read samples
errF <- learnErrors(filtFs1, nbases=1e8, multithread=TRUE) # 100790190 total bases in 373297 reads from 16 samples will be used for learning the error rates 6 rounds
errR <- learnErrors(filtRs1, nbases=1e8, multithread=TRUE) # 102830400 total bases in 571280 reads from 25 samples will be used for learning the error rates 9 rounds

# Plot profiles
plotQualityProfile(fnFs[1:2])
plotQualityProfile(fnRs[1:2])
plotErrors(errF, nominalQ=TRUE)
plotErrors(errR, nominalQ=TRUE)

# Merge should only keep what you specified with 'keep' object that filterd low count reads.
mergers<- vector("list", length(sample.names.kept))
names(mergers) <- sample.names.kept

# Simplified clump for dereplicating and inferring sequence variants while sparing your console.
for(sam in sample.names.kept) {
  cat("Processing:", sam, "\n")
  derepF <- derepFastq(filtFs1[[sam]])
  ddF <- dada(derepF, err=errF, multithread=TRUE)
  derepR <- derepFastq(filtRs1[[sam]])
  ddR <- dada(derepR, err=errR, multithread=TRUE)
  merger <- mergePairs(ddF, derepF, ddR, derepR)
  mergers[[sam]] <- merger # Use mergers as it's a list of accepted SampleID's
}


seqtab <- makeSequenceTable(mergers)
dada2:::pfasta(head(colnames(seqtab)))
# Export this sequence table as a serialized R object that can be easily imported with low memory costs later on
# Not recommended that subsequent workflow use CSV or any delimited format.
dim(seqtab)

# It's generally wise to have fully tracked reads through the pipeline, so even though we're not going to merge seqtabnochim before we assign taxonomy, I want to have the N reads passed after removing chimeras.
seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE, verbose=TRUE)
dim(seqtab.nochim)
# Identified 144543 bimeras out of 234064 input sequences - dim:649x89521

# The invdividual dim from seqtab.nochim for each run should be consistent with pooling runs and removing chimeras dim.
# That way we know if our track object is telling the truth on where our reads went in chimera removal; run by run basis or all together
getN <- function(x) sum(getUniques(x))
# This is input/filterd/merged/tabled in seqtab, not seqtab.nochim - self note - disregard
# Do not need to track denoised reads with base dada2 algorithm, it doesn't remove any reads through denoising.

track <- cbind(out, sapply(mergers, getN), rowSums(seqtab), rowSums(seqtab.nochim))
# If processing a single sample, remove the sapply calls: e.g. replace sapply(dadaFs, getN) with getN(dadaFs)
colnames(track) <- c("input", "filtered", "merged", "seqtab", "seqtabnochim")
rownames(track) <- sample.names
head(track)

# We will assign taxonomy after all runs are done - after merging all 3 ASV tables from each run.

# save rds per your workflow directory
saveRDS(seqtab, file="~/pm_workflow/run1/seqtabRun1.rds")
saveRDS(seqtab.nochim, file="~/pm_workflow/run1/seqtab.nochimRun1.rds")
write.csv(track, "~/pm_workflow/run1/trackedreads.csv")
# Save the environment for sequencing run 1
save(list = ls(), file = "~/pm_workflow/run1/dada2run1.RData")

#Sequencing run 2/3
# =========================================================
#       Run 2 - pertains to sequencing run 2/3
# =========================================================

rm(list=ls())
seed=81
set.seed(seed)

setwd("~/pm_workflow") # CHANGE ME TO THE PROPER DIRECTORY (where you want this run's outputs to go)
library(dada2); packageVersion("dada2")

path <- "~/pm_workflow/data/sequencing/inputs-fastqs/run2" # CHANGE ME TO THE PROPER DIRECTORY (containing the fastq files)
list.files(path)

fnFs <- sort(list.files(path, pattern="_R1_001.fastq", full.names = TRUE))
fnRs <- sort(list.files(path, pattern="_R2_001.fastq", full.names = TRUE))
sample.names <- sapply(strsplit(basename(fnFs), "_"), `[`, 1)


filtFs <- file.path(path, "filtered", paste0(sample.names, "_F_filt.fastq.gz"))
filtRs <- file.path(path, "filtered", paste0(sample.names, "_R_filt.fastq.gz"))

out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, trimLeft=c(19,20), truncLen=c(280,190),
                     maxN=0, maxEE=c(2,5), truncQ=2, rm.phix=TRUE,
                     compress=TRUE, multithread=TRUE)

keep <- out[,"reads.out"] > 50
filtFs1 <- file.path(filtFs)[keep]
filtRs1 <- file.path(filtRs)[keep]

sample.namesF <- sapply(strsplit(basename(filtFs1), "_"), `[`, 1) 
sample.namesR <- sapply(strsplit(basename(filtRs1), "_"), `[`, 1) 
if(!identical(sample.namesF, sample.namesR)) stop("Forward and reverse files do not match.")

sample.names.kept <- sapply(strsplit(basename(filtFs1), "_"), `[`, 1)
names(filtFs1) <- sample.names.kept
names(filtRs1) <- sample.names.kept
set.seed(seed)


errF <- learnErrors(filtFs1, nbases=1e8, multithread=TRUE) #105837840 total bases in 391992 reads from 15 samples
errR <- learnErrors(filtRs1, nbases=1e8, multithread=TRUE) #102640680 total bases in 570226 reads from 22 samples

plotQualityProfile(fnFs[1:2])
plotQualityProfile(fnRs[1:2])
plotErrors(errF, nominalQ=TRUE)
plotErrors(errR, nominalQ=TRUE)

mergers<- vector("list", length(sample.names.kept))
names(mergers) <- sample.names.kept

for(sam in sample.names.kept) {
  cat("Processing:", sam, "\n")
  derepF <- derepFastq(filtFs1[[sam]])
  ddF <- dada(derepF, err=errF, multithread=TRUE)
  derepR <- derepFastq(filtRs1[[sam]])
  ddR <- dada(derepR, err=errR, multithread=TRUE)
  merger <- mergePairs(ddF, derepF, ddR, derepR)
  mergers[[sam]] <- merger 
}


seqtab <- makeSequenceTable(mergers)
dada2:::pfasta(head(colnames(seqtab)))

dim(seqtab)

seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE, verbose=TRUE)
dim(seqtab.nochim)

getN <- function(x) sum(getUniques(x))

track <- cbind(out, sapply(mergers, getN), rowSums(seqtab), rowSums(seqtab.nochim))
colnames(track) <- c("input", "filtered", "merged", "seqtab", "seqtabnochim")
rownames(track) <- sample.names
head(track)

# save rds per your workflow directory
saveRDS(seqtab, file="~/pm_workflow/run2/seqtabRun2.rds")
saveRDS(seqtab.nochim, file="~/pm_workflow/run2/seqtab.nochimRun2.rds")
write.csv(track, "~/pm_workflow/run2/trackedreads.csv")
# Save the environment for sequencing run 1
save(list = ls(), file = "~/pm_workflow/run2/dada2run2.RData")


#Sequencing run 3/3
# =========================================================
#       Run 3 - pertains to sequencing run 3/3
# =========================================================

rm(list=ls())
seed=81
set.seed(seed)

setwd("~/pm_workflow")  #CHANGE ME TO THE PROPER DIRECTORY contaning sequencing run
library(dada2); packageVersion("dada2")

path <- "~/pm_workflow/data/sequencing/inputs-fastqs/run3" # CHANGE ME to the directory containing the fastq files after unzipping.
list.files(path)

fnFs <- sort(list.files(path, pattern="_R1_001.fastq", full.names = TRUE))
fnRs <- sort(list.files(path, pattern="_R2_001.fastq", full.names = TRUE))
sample.names <- sapply(strsplit(basename(fnFs), "_"), `[`, 1)


filtFs <- file.path(path, "filtered", paste0(sample.names, "_F_filt.fastq.gz"))
filtRs <- file.path(path, "filtered", paste0(sample.names, "_R_filt.fastq.gz"))

out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, trimLeft=c(19,20), truncLen=c(280,190),
                     maxN=0, maxEE=c(2,5), truncQ=2, rm.phix=TRUE,
                     compress=TRUE, multithread=TRUE)

keep <- out[,"reads.out"] > 50 
filtFs1 <- file.path(filtFs)[keep]
filtRs1 <- file.path(filtRs)[keep]

sample.namesF <- sapply(strsplit(basename(filtFs1), "_"), `[`, 1) 
sample.namesR <- sapply(strsplit(basename(filtRs1), "_"), `[`, 1) 
if(!identical(sample.namesF, sample.namesR)) stop("Forward and reverse files do not match.")

sample.names.kept <- sapply(strsplit(basename(filtFs1), "_"), `[`, 1)
names(filtFs1) <- sample.names.kept
names(filtRs1) <- sample.names.kept
set.seed(seed)

errF <- learnErrors(filtFs1, nbases=1e8, multithread=TRUE) # 107697330 total bases in 398879 reads from 11 samples will be used for learning the error rates.
errR <- learnErrors(filtRs1, nbases=1e8, multithread=TRUE) # 103126320 total bases in 572924 reads from 18 samples will be used for learning the error rates.

plotQualityProfile(fnFs[1:2])
plotQualityProfile(fnRs[1:2])
plotErrors(errF, nominalQ=TRUE)
plotErrors(errR, nominalQ=TRUE)


mergers<- vector("list", length(sample.names.kept))
names(mergers) <- sample.names.kept

for(sam in sample.names.kept) {
  cat("Processing:", sam, "\n")
  derepF <- derepFastq(filtFs1[[sam]])
  ddF <- dada(derepF, err=errF, multithread=TRUE)
  derepR <- derepFastq(filtRs1[[sam]])
  ddR <- dada(derepR, err=errR, multithread=TRUE)
  merger <- mergePairs(ddF, derepF, ddR, derepR)
  mergers[[sam]] <- merger 
}


seqtab <- makeSequenceTable(mergers)
dada2:::pfasta(head(colnames(seqtab)))
dim(seqtab)

seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE, verbose=TRUE)
dim(seqtab.nochim)

getN <- function(x) sum(getUniques(x))

track <- cbind(out, sapply(mergers, getN), rowSums(seqtab), rowSums(seqtab.nochim))
colnames(track) <- c("input", "filtered", "merged", "seqtab", "seqtabnochim")
rownames(track) <- sample.names
head(track)

saveRDS(seqtab, file="~/pm_workflow/run3/seqtabRun3.rds")
saveRDS(seqtab.nochim, file="~/pm_workflow/run3/seqtab.nochimRun3.rds")
write.csv(track, "~/pm_workflow/run3/trackedreads.csv")
# Save the environment for sequencing run 1
save(list = ls(), file = "~/pm_workflow/run3/dada2run3.RData")


# =========================================================
#       Merging Runs
# =========================================================

# Reproducible workflow for run-per-run dada2 processing.
seed=81
set.seed(seed)
rm(list=ls())
setwd("~/pm_workflow")  #CHANGE ME TO THE PROPER DIRECTORY contaning where the SILVA training set is
library(dada2); packageVersion("dada2")

# Merge multiple runs (if necessary) - import sequence tables - change to proper directory
seqtab1 <- readRDS("~/pm_workflow/run1/seqtabRun1.rds")
seqtab2 <- readRDS("~/pm_workflow/run2/seqtabRun2.rds")
seqtab3 <- readRDS("~/pm_workflow/run3/seqtabRun3.rds")
seqtab.all <- mergeSequenceTables(seqtab1, seqtab2, seqtab3)
dim(seqtab.all) # 1514 x 514359

# Remove chimeras
seqtabmergedNoC <- removeBimeraDenovo(seqtab.all, method="consensus", multithread=TRUE)
dim(seqtabmergedNoC) #1514 x 117397
# Assign taxonomy
taxa <- assignTaxonomy(seqtabmergedNoC, "~/pm_workflow/data/sequencing/training-species-assignment/silva_nr_v132_train_set.fa", multithread=TRUE)
head(taxa)
dim(taxa)

taxasp <- addSpecies(taxa, "~/pm_workflow/.../input/Training/silva_species_assignment_v132.fa") # Directory where SILVA training set is

getN <- function(x) sum(getUniques(x)) # Reference for tracking seqtab

# Track reads after chimera removal
trackmerged <- cbind(rowSums(seqtab.all), rowSums(seqtabmergedNoC))
colnames(trackmerged) <- c("seqtaball", "seqtabmergednoc")
head(trackmerged)


write.csv(trackmerged, file = "~/pm_workflow/mergedruns/trackmerged.csv")
saveRDS(seqtab.all, file= "~/pm_workflow/mergedruns/seqtaball.rds")
saveRDS(seqtabmergedNoC, file="~/pm_workflow/mergedruns/seqtabmergedNoC.rds")

saveRDS(taxa, file="~/pm_workflow/mergedruns/taxa.rds")
saveRDS(taxasp, file="~/pm_workflow/mergedruns/taxasp.rds")

# ===========================================================================================
# trackmerged: tracks the reads through the pipeline from raw to non-chimeric final counts.
# seqtab.all: sequence count table pre chimera removal.
# seqtabmergedNoC: sequence count table after removing chimeras - this is what we will use for subsequent analysis
# taxa: taxonomic table corresponding to each ASV in seqtabmergedNoC
# taxasp: taxonomic table corresponding to each ASV in seqtabmergedNoC with species assignment
# ===========================================================================================
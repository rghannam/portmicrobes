# =========================================================
# Initial setup for analysis
# Required packages for port microbes reproducible workflow
# =========================================================

# ==== Detach all packages ====
detachAllPackages <- function() {
  
  basic.packages <- c("package:stats","package:graphics","package:grDevices","package:utils","package:datasets","package:methods","package:base")
  
  package.list <- search()[ifelse(unlist(gregexpr("package:",search()))==1,TRUE,FALSE)]
  
  package.list <- setdiff(package.list,basic.packages)
  
  if (length(package.list)>0)  for (package in package.list) detach(package, character.only=TRUE)
  
}
detachAllPackages()

# Required packages
.cran_packages  <-  c("caret", "caretEnsemble", "e1071", "ggplot2",
                      "ggmap", "maps", "mapsdata", "randomForest",
                      "vegan", "plyr", "dplyr", "reshape", "stringr", "tidyr", "tidyverse",
                      "reshape2", "devtools", "rlang", "Rcpp", "csv", "tibble", "plotly",
                      "magrittr", "data.table")

.bioc_packages <- c("phyloseq", "DESeq2", "dada2")

.github_packages <- c("umerijaz/microbiomeSeq", "elsayed-lab/hpgltools")


# Install CRAN packages (if not already installed)
.inst <- .cran_packages %in% installed.packages()
if (any(!.inst)) {
  install.packages(.cran_packages[!.inst], repos = "http://cran.rstudio.com/")
}


# Install Bioconductor packages
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install()

BiocManager::install(c("phyloseq", "DESeq2", "dada2", "GenomicFeatures", "hpgltools", "microbiome"))


# Alternative (not run)
.inst <- .bioc_packages %in% installed.packages()
if (any(!.inst)) {
  source("http://bioconductor.org/biocLite.R")
  biocLite(.bioc_packages[!.inst])
}


# Install GitHub packages
.inst <- .github_packages %in% installed.packages()
if (any(!.inst)) {
  devtools::install_github(.github_packages[!.inst])
}


# Load packages
pkgs <- c("caret", "caretEnsemble", "e1071", "ggplot2", "ggmap", "maps",
          "mapsdata", "randomForest","vegan", "plyr", "dplyr", "reshape", 
          "stringr", "tidyr", "tidyverse",
          "reshape2", "devtools", "rlang", "Rcpp", "csv", "tibble", "plotly",
          "phyloseq", "DESeq2", "dada2", "microbiomeSeq", "GenomicFeatures", "hpgltools", "microbiome") # Load the PM package where all data is stored compressed for easy accessibility.
invisible(sapply(pkgs,require, character = TRUE))



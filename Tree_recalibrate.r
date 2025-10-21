library(PRROC)
library(phangorn)
library(plyr)
library(reshape2)
library(RERconverge)

# Get command-line arguments
args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 2) {
  stop("Usage: Rscript estimate_tree.R <alignment_directory_name> <location>")
}

alndir_name <- args[1]
alndir_location <- args[2]
setwd(alndir_location)

# Construct output filename automatically
output_file <- paste0(alndir_name, ".recal.tre")

estimatePhangornTreeAll(
  alndir      = alndir_name,
  treefile    = "species_tree.tre",
  output.file = output_file
)
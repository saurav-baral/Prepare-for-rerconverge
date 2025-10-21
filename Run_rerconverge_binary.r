library(RERconverge)

args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 2) {
  stop("Usage: Rscript estimate_tree.R <alignment_directory_name> <location>")
}

rerconverge_location <- args[2]

setwd(rerconverge_location)

Trees=readTrees(paste('recalibrated_trees.tree',sep=""))
mamRERw=RERconverge::getAllResiduals(Trees, transform = "sqrt",  n.pcs = 0, use.weights = T,
weights=NULL,norm="scale")

phenotype_foreground <- strsplit(args[1], ",")[[1]]

phenvBin=foreground2Paths(phenotype_foreground, Trees, clade="all")
res=correlateWithBinaryPhenotype(mamRERw, phenvBin, min.sp=10, min.pos=2, weighted="auto", winsorizeRER=NULL, winsorizetrait=NULL,
bootstrap=TRUE, bootn=100, sort = T)




write.csv(res, 'correlations.csv')


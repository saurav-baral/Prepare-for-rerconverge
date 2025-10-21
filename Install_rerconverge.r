required_pkgs <- c("usethis", "devtools", "PRROC", "phangorn", "plyr", "reshape2")
for (pkg in required_pkgs) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg, repos = "https://cloud.r-project.org")
  }
}

if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager", repos = "https://cloud.r-project.org")
}

bioc_pkgs <- c("impute", "ggtree")

# Find which are missing
missing_pkgs <- bioc_pkgs[!sapply(bioc_pkgs, requireNamespace, quietly = TRUE)]

# Install only if missing
if (length(missing_pkgs) > 0) {
  BiocManager::install(missing_pkgs, ask = FALSE, update = FALSE)
}


library(devtools)

if (!requireNamespace("RERconverge")) {
  devtools::install_github("nclark-lab/RERconverge")
}
library(RERconverge)


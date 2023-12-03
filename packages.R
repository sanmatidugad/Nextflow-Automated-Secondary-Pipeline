if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install(version = "3.16")

rm(list = ls())
Sys.setenv(TZ='EDT')

ReQ_packages = c("dplyr", "optparse","data.table", "tools", "colorspace")
BioCpackages = c("subSeq", "edgeR")

for (pack in 1:length(ReQ_packages)) {
  if(ReQ_packages[pack] %in% rownames(installed.packages()) == FALSE) {
    install.packages(ReQ_packages[pack], repos='http://cran.us.r-project.org')
  }
}

for (pack in 1:length(BioCpackages)){
  if (BioCpackages[pack] %in% rownames(installed.packages()) == FALSE){
    BiocManager::install(BioCpackages[pack])
  }
}



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

suppressPackageStartupMessages(suppressWarnings(library(optparse)))
suppressPackageStartupMessages(suppressWarnings(library(subSeq)))
suppressPackageStartupMessages(suppressWarnings(library(data.table)))
suppressPackageStartupMessages(suppressWarnings(library(dplyr)))
suppressPackageStartupMessages(suppressWarnings(library(tools)))
suppressPackageStartupMessages(suppressWarnings(library(colorspace)))

args=commandArgs(trailingOnly = T)

data = NULL    #count_matrix_dir
target = NULL   #normal target 
output = NULL
normal_count = NULL    #number of normal samples
project = NULL

## all these go in the command
#data = "gtex_tcga_ovary_100.csv"
#target = 1000000
#normal_count = 100
#project = "test1"

if ("-data" %in% args && "-normal_count" %in% args && "-target" %in% args && "-project" %in% args) {    
  data_index = which(args == "-data")    #Getting index of the "--input" option
  normal_count_index = which(args == '-normal_count')
  target_index = which(args == "-target")
  project_index = which(args == "-project")
  
  data = args[data_index + 1]
  normal_count = as.numeric(args[normal_count_index + 1]) 
  target = as.numeric(args[target_index + 1])
  project <- args[project_index + 1]
}

gene_count_data = read.csv(data, row.names = 1, header = TRUE, stringsAsFactors = FALSE, strip.white = TRUE, check.names = FALSE)
data = data.frame(gene_count_data)

options(scipen = 20)

main <- function(counts_matrix_dir, output_dir, target){
  df = data.frame(counts_matrix_dir)
  #df = read.csv(counts_matrix_dir, row.names = 1, header = TRUE, as.is = TRUE, stringsAsFactors = FALSE, strip.white = TRUE, check.names=FALSE)
  #print(colnames(df))
  check_df_input(df)
  num_cols <- unlist(lapply(df, is.numeric))
  
  if (length(num_cols) != ncol(df)){
    stop("ERROR: Type Error. Please check that columns contain only numeric values")
  }
  
  gc_totals <- apply(df[, num_cols], 2, sum)
  
  subsampled_df = df
  #seed = gen_seed(df[,4:ncol(df)], condition_list)
  subsampled_df[,num_cols] <- subsample_seq(subsampled_df[,num_cols], target)
  
  #Sum total number of gene counts for each sample
  new_totals <- apply(subsampled_df[, num_cols], 2, sum)
  gc_totals <-cbind(gc_totals, new_totals)
  
  #Strip extension from filename
  filename <- file_path_sans_ext(basename(project))
  
  #Add a subsampling tag to output file name
  fn <- paste(filename, "-", sprintf("%.0f", target), ".csv", sep="")
  
  #Use rownames as first column and add the label the column as "SYMBOL"
  write.csv(subsampled_df, paste(output_dir), row.names=TRUE)
  print(paste("Writing subsampled data to", sprintf("%.0f", target), "to", file.path(output_dir)))
  
  fn_gc <- paste(file_path_sans_ext(output_dir), ".gc_totals.csv", sep="")
  print(paste("Writing data gc total counts to", sprintf("%.0f", target), "to", fn_gc))
  write.csv(gc_totals,fn_gc)
  return(subsampled_df)
  
}

subsample_seq <- function(df,target, condition_list){
  output <- data.frame(matrix(0, ncol = ncol(df), nrow = nrow(df)))
  colnames(output) <- colnames(df)
  condition_list <- create_condition_list(df[,1:ncol(df)])
  seed = gen_seed(df[,1:ncol(df)], condition_list)
  
  for(i in 1:ncol(df)){
    total_gc <- sum(df[,i])
    if(sum(total_gc)>target){
      reduce_by <- target/total_gc    #proportion
      ssdata <- generateSubsampledMatrix(df, reduce_by, seed, replication = 1) 
      output[,i] <-ssdata[,i]
    } else{
      output[,i] <- df[,i]}
  }
  return(output)
}

create_condition_list <- function(df){
  #Make a list of samples
  samples <- names(df)
  samples <- as.matrix(samples)
  
  #Construct a matrix assigning ~1/2 samples to condition "A" and the rest to condition "B"
  #This is a hack needed to utilize the subsampling function from the subSeq library that was created for differential expression 
  sample_type <- 1:length(samples)
  num <- round(length(samples)/2)
  sample_type[1:num] <- "A"
  sample_type[(num+1):length(samples)] <- "B"
  
  #Convert to factorized list
  condition_list <- factor(c(sample_type))
  return(condition_list)
}

gen_seed <- function(df, condition_list){
  #Get ready for subsampling: this first operation is just performed to get the seed. Using a subset of data speeds it up. 
  #To generate the seed, a list of samples need to be provided with at least 2 different conditions, i.e. "A" and "B"
  ### There are issues creating a seed when the sample is too small--need to check for a minimum of 8 samples
  num = ncol(df)/2
  ss = subsample(df[,num:(num+4)], c(.01, .5), method="edgeR", treatment=condition_list[num:(num+4)])
  #Get the seed and start subsampling
  seed = getSeed(ss)
  print(paste("seed value used for subsampling:", seed))
  return(seed)
}

check_df_input <- function(df){
  if (ncol(df) <= 8 ){
    stop("Data must have > 8 columns")
  }
}

output = main(data, output_dir = paste(project, "-subsampled-", target,".csv",sep = ""), target=target)

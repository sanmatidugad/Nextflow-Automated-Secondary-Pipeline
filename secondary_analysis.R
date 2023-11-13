rm(list=ls())    #clear the environment in R. Remove all objects (variables, functions, etc.) from the current workspace.
Sys.setenv(TZ='EDT')

ReQ_Modules = c("plyr","dplyr", "optparse", "matrixStats", "ggplot2", "reshape2")#, "tidyverse")


for (each in 1:length(ReQ_Modules)) {
  if(ReQ_Modules[each] %in% rownames(installed.packages()) == FALSE) {
    install.packages(ReQ_Modules[each])
  }
}

#invisible(lapply(ReQ_Modules, library, character.only = TRUE))

suppressPackageStartupMessages(suppressWarnings(library(plyr)))
suppressPackageStartupMessages(suppressWarnings(library(dplyr)))
suppressPackageStartupMessages(suppressWarnings(library(optparse)))
suppressPackageStartupMessages(suppressWarnings(library(matrixStats)))
suppressPackageStartupMessages(suppressWarnings(library(ggplot2)))
suppressPackageStartupMessages(suppressWarnings(library(reshape2)))

args=commandArgs(trailingOnly = T)

input_file = NULL
reference = NULL
project = NULL
normal_count = NULL
tumor_count = NULL

if ("-i" %in% args && "-r" %in% args && "-p" %in% args && "-n" %in% args){ # && "-tumor_columns" %in% args) {    
  input_index = which(args == "-i")    #Getting index of the "--input" option
  reference_index = which(args == "-r")
  project_index = which(args == "-p")
  normal_columns_index = which(args == "-n")
  #tumor_columns_index = which(args == "--tumor_columns")
  
  input_file <-  args[input_index + 1] #"union_normalized.txt"
  reference <- args[reference_index + 1] #"ID_SYMBOL.csv"
  project <- args[project_index + 1]
  normal_count = as.numeric(args[normal_columns_index + 1])
  
  #tumor_count = as.numeric(tumor_columns_index + 1)
}

print(paste("your input file is:", input_file))


#This function calculates mean of third quartile of all  normalized expression profiles
get_Q3_mean <- function(df){
  df.matrix <- data.matrix(df)    #convert data frame to data.matrix
  Quartile <- colQuantiles(df.matrix)    #calculate Quartiles for each normalized expression profile
  Q3.mean <- mean(Quartile[,"75%"])    #take mean of the third quartile
  return(Q3.mean)
} 

#input_file = "test1_normalized.csv"
#reference = "ID_SYMBOL.csv"
#normal_count = 10

main = function(input_file, reference, project){
  #### PART A
  df1 = read.csv(input_file, row.names = 1)    #reading input file
  print(paste("Dimension of the gene count table is", dim(df1)))
  
  gene_symbols = read.csv(reference, sep = ",", row.names = 2)    #reading reference file
  rownames(gene_symbols) = gsub("\\.\\d+$", "", rownames(gene_symbols))
  head(gene_symbols)
  
  df = df1 %>% group_by(rowname = row.names(.)) %>% summarise_all(sum) %>% ungroup()
  #summing up the duplicates.
  df = as.data.frame(df)
  rownames(df) = df$rowname  #replace row names with first column
  df = df[,c(-1)]   #remove gene id column
  
  Q3.mean = get_Q3_mean(df)    #calculate the mean of third quartile of each normalized expression profile
  df = df/Q3.mean*1000+1    #normalize the data by diving by third quartile, multiple by 1000 and add one to prevent log from going to Inf
  df = log2(df)    #take log2 of normalized data
  is.na(df) <- sapply(df, is.infinite)
  df[is.na(df)] <- 0
  df = round(df, 6)
  
  df = merge(gene_symbols, df, by = "row.names", all.y = T)    #merge gene names using gene ids.
  
  write.csv(df, file = paste(project, "-", "scaled.csv", sep = ""), row.names = F)
  print(paste("Scaled data file was successfully written to", paste(project, "-", "scaled.csv", sep = "")))
  
  #### PART B
  rownames(df) = df[,1]    #after merge make gene id as row names
  df = df[, c(-1)]    #remove the gene id column 
  df = df[!duplicated(rownames(df)),]    #remove duplicate row names
  normal_samples = df[,c(2:(normal_count+1))]    #the first column is SYMBOL and normals start from the second column
  #print(colnames(normal_samples))
  normal_samples = t(normal_samples)
  normal_samples = as.data.frame(normal_samples)
  norms.m <- melt(normal_samples, measure.vars = 1:ncol(normal_samples),  value.name = "Expression", variable.name="SYMBOL")
  A <- boxplot(Expression~SYMBOL, data=norms.m, plot = F, range = 1.5)
  boxstats <- A$stats
  colnames(boxstats)<-A$names
  rownames(boxstats)<-c('Min','Q1','Median','Q3','Max')
  boxstats <-data.frame(t(boxstats))
  boxstats$IQR <-as.matrix(boxstats[1:nrow(boxstats),"Max"]-boxstats[1:nrow(boxstats), "Min"])
  boxstats$HIGH <- boxstats[1:nrow(boxstats), "Max"]
  boxstats$LOW <- boxstats[1:nrow(boxstats), "Min"]
  boxstats$VERY_HIGH <- boxstats[1:nrow(boxstats), "Median"]+2*boxstats[1:nrow(boxstats),"IQR"]
  boxstats$VERY_LOW <- boxstats[1:nrow(boxstats), "Median"]-2*boxstats[1:nrow(boxstats),"IQR"]
  boxstats1 = merge(gene_symbols, boxstats, by = "row.names", all.y = T)
  write.csv(boxstats1, paste(project, "-" , "boxstats_threshold.csv", sep = ""), row.names = F)
  print(paste("boxplot threshold data file was successfully written to", paste(project, "-" , "boxstats_threshold.csv", sep = "")))
  
  df = df[,c(-1)]
  num_genes <- ncol(t(df))
  genes.over <-as.matrix(df[1:num_genes, ]  > boxstats[1:num_genes, "Max"])
  genes.under <-as.matrix(df[1:num_genes, ]  < boxstats[1:num_genes, "Min"])
  IQR <-as.matrix(boxstats[1:num_genes, "Max"]-boxstats[1:num_genes, "Min"])
  genes.high <-as.matrix(df[1:num_genes, ]  > boxstats[1:num_genes, "Median"]+2*IQR[1:num_genes])
  genes.low <-as.matrix(df[1:num_genes, ]  < boxstats[1:num_genes, "Median"]-2*IQR[1:num_genes])
  
  genes.states <- genes.high
  genes.states[] <- "NORMAL"
  genes.states[genes.over] <- "HIGH"
  genes.states[genes.under] <- "LOW"
  genes.states[genes.high] <- "VERY HIGH"
  genes.states[genes.low] <- "VERY LOW"
  
  ### melt gene states so that for each sample, we know aberrant expressed genes
  genes.states.m <- melt(genes.states, varnames = c('SYMBOL', 'SAMPLE_ID'))
  genes.stats.m.aberrant <- genes.states.m[genes.states.m$value != 'NORMAL',]
  
  genes.states1 = as.data.frame(genes.states)
  dim(genes.states1)
  
  genes.states1 = merge(gene_symbols, genes.states1, by = "row.names", all = T)
  
  write.csv(genes.states1, paste(project, "-" , "gene_states.csv", sep = ""), row.names = F)    # # #Save data to files
  print(paste("gene states data file was successfully written to", paste(project, "-" , "gene_states", sep = "")))
  write.csv(genes.stats.m.aberrant, paste(project, "-" , "aberant_gene_states.csv", sep = ""), row.names = F)
  print(paste("aberant gene states data file was successfully written to", paste(project, "-" , "aberant_gene_states.csv", sep = "")))
  
  ## PART C
  exclude_normal_samples = genes.stats.m.aberrant[!genes.stats.m.aberrant$SAMPLE_ID %in% rownames(normal_samples),] 
  exclude_normal_samples
  gene_stats = table(exclude_normal_samples$SYMBOL, exclude_normal_samples$value)
  #write.csv(gene_stats, "curated_gene_state_statistics.csv", row.names = T)
  
  format_gene_stats = dcast(data.frame(gene_stats), Var1~Var2, value.var = "Freq")
  rownames(format_gene_stats) = format_gene_stats[,1]
  format_gene_stats = format_gene_stats[,-1]
  
  tumor_count = dim(df)[2] - normal_count   #number of tumor samples = total columns - normal normal
  print(paste("Tumor sample count is:", tumor_count))
  
  format_gene_stats$aberant_sum = rowSums(format_gene_stats)
  format_gene_stats$NORMAL = tumor_count - format_gene_stats$aberant_sum
  format_gene_stats$pc_HIGH = round((format_gene_stats$HIGH/tumor_count)*100, 2)
  format_gene_stats$pc_LOW = round((format_gene_stats$LOW/tumor_count)*100, 2)
  format_gene_stats$pc_VERY_LOW = round((format_gene_stats$`VERY LOW`/tumor_count)*100,2)
  format_gene_stats$pc_VERY_HIGH = round((format_gene_stats$`VERY HIGH`/tumor_count)*100,2)
  
  format_gene_stats = format_gene_stats[, -which(names(format_gene_stats) == "aberant_sum")]
  
  format_gene_stats = merge(gene_symbols, format_gene_stats, by = "row.names", all.y = T)
  write.csv(format_gene_stats, paste(project, "-", "gene_statistics.csv", sep = ""), row.names = F )
  print(paste("gene wise statistics was successfully written to", paste(project, "-", "gene_statistics.csv", sep = "")))
  
  ## PART 4
  rownames(genes.states1) = genes.states1[,1]    #converting first column to rownames
  genes.states1 = genes.states1[,c(-1)]    #excluding the first  column
  genes.states1 = genes.states1[!duplicated(rownames(genes.states1)),]
  genes.states1 = genes.states1[,c(-1)]    #excluding symbol column
  
  #count_gene_states  = genes.states1 %>% map(function(x) table(x))
  count_gene_states = lapply(data.frame(genes.states1), table)
  split_data = ldply(count_gene_states, data.frame)
  
  #reshape_data = data.frame(pivot_wider(split_data, id_cols = ".id", names_from = "Var1", values_from = "Freq"))    #from tidyverse
  reshape_data = reshape(split_data, idvar = ".id", timevar = "Var1", direction = "wide")
  reshape_data <- replace(reshape_data, is.na(reshape_data), 0)
  
  write.csv(reshape_data, paste(project, "-", "sample_statistics.csv", sep = ""), row.names = F)
  print(paste("statistics for each samples was successfully written to: ", paste(project, "-", "sample_statistics.csv", sep = "")))
  
}

main(input_file, reference, project)

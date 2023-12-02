## MAINTAINER: Sanmati Dugad

## Command Line:
# Rscript z_analysis.v1.R \
# -i <test>_normalized.txt \
# -s <test>-gene_states.csv \
# -gc <test>_subsampled_<reads>.gc_totals.csv \
# -n normal_count {numeric value}
# -a actionable_data {gene export file}

rm(list=ls())    #clear the environment in R. Remove all objects (variables, functions, etc.) from the current workspace.
Sys.setenv(TZ='EDT')

ReQ_Modules = c("dplyr", "tidyr", "data.table", "matrixStats")

for (each in 1:length(ReQ_Modules)) {
  if(ReQ_Modules[each] %in% rownames(installed.packages()) == FALSE) {
    install.packages(ReQ_Modules[each])
  }
} 

suppressPackageStartupMessages(suppressWarnings(library(dplyr)))
suppressPackageStartupMessages(suppressWarnings(library(tidyr)))
suppressPackageStartupMessages(suppressWarnings(library(data.table)))
suppressPackageStartupMessages(suppressWarnings(library(matrixStats)))
suppressPackageStartupMessages(suppressWarnings(library(utils)))

args=commandArgs(trailingOnly = T)

normalized_data = NULL    #normalized data
states_data = NULL    #gene states
gene_counts = NULL    #gene counts file from subsampled output
actionable_data = NULL    #actionable gene list
normal_count = NULL    #number of normal samples

if ("-i" %in% args ){ #&& "-s" %in% args && "-gc" %in% args && "-n" %in% args && "-a" %in% args ) {    
  normalized_data_index = which(args == "-i")
  states_data_index = which(args == "-s")
  gene_counts_index = which(args == "-gc")
  normal_count_index = which(args == "-n")
  actionable_data_index = which(args == "-a")
  
  normalized_data = args[normalized_data_index + 1] 
  states_data =  args[states_data_index + 1] 
  gene_counts = args[gene_counts_index + 1]
  normal_count = args[normal_count_index + 1]
  actionable_data = args[actionable_data_index + 1]
}

## Assign Varaibles:
# setwd("/home/data/Nextflow_Secondary_Pipeline/")
# normalized_data = "20231117_Lymph_Node_Z-scores/lt1_normalized.csv"


normalized_input = read.csv(normalized_data, row.names = 1)
normal_count = as.numeric(normal_count)
tumor_sample = colnames(normalized_input[,(normal_count+1):(ncol(normalized_input)), drop = FALSE])

main = function(normalized_data, states_input, gene_counts, actionable_data, i){
  print(paste(i))
  normalized_input = read.csv(normalized_data, row.names = 1)
  states_input = read.csv(states_data, row.names = 1)
  counts_input = read.csv(gene_counts, row.names = 1)
  actionable_input = read.csv(actionable_data)
  normal_count = as.numeric(normal_count)
  
  normalized_input[, c(i)] = round(normalized_input[, c(i)], 3)    # round up tumor data point to 
  normalized_input$Normalized_Per_Million = round((normalized_input[,c(i)] * 1000000) / counts_input[c(i), c("new_totals")] , 3)    #Per Million Normalized Data
  normalized_input$Normal_Mean = round(rowMeans(normalized_input[,1:normal_count]), 3)  # calculating Normal Means
  normalized_input= normalized_input %>% mutate(Normal_Mean = ifelse(i != 0 & Normal_Mean == 0, 1, Normal_Mean))
  normalized_input$Normal_Mean_per_Million = round(rowMeans((normalized_input[,1:normal_count]*1000000)/ counts_input[1:normal_count, c("new_totals")]) , 3)
  normalized_input$Normal_Median = round(apply(normalized_input[,1:normal_count], MARGIN = 1, median), 4)    #calculating Normal Medians
  normalized_input$Normal_SD = round(apply(normalized_input[,1:normal_count], 1, sd) , 3)    #calculating Standard Deviation 
  normalized_input$Normal_SD[normalized_input$Normal_SD == 0] = 1    # replacing standard deviation of zero by 1
  normalized_input$Z_score = round((normalized_input[, c(i)] - normalized_input$Normal_Mean) / normalized_input$Normal_SD , 3)   #calculating Z-scores
  normalized_input$Foldchange = round(normalized_input[, c(i)] / normalized_input$Normal_Mean, 3)    # calculate foldchange values
  
  normalized_input$Z_score[is.na(normalized_input$Z_score)] = 0    #performing for a single vector at a time. cost of is.na function is n^2
  normalized_input$Z_score[normalized_input$Z_score == Inf] = 0
  normalized_input$Z_score[normalized_input$Z_score == -Inf] = 0
  
  normalized_input$Foldchange[is.na(normalized_input$Foldchange)] = 0
  normalized_input$Foldchange[normalized_input$Foldchange == Inf] = 0
  normalized_input$Foldchange[normalized_input$Foldchange == -Inf] = 0
  
  states_tumor = states_input[, c("SYMBOL", i), drop = FALSE]    #selecting SYMBOL column and tumor sample column in gene state data
  colnames(states_tumor) = c("SYMBOL", "Gene_State")    #renaming tumor sample in gene state data    
  normalized_input.1 = merge(normalized_input, states_tumor, by = "row.names", all.x = TRUE)    #merging normalized data table with gene state tumor data
  
  #Adding False Positive and False Negative Check
  normalized_input.1$Potential_FP_or_FN = NA
  #FALSE POSITIVE: median normal < 10, normalized expression < 10, Aberrantly Called (Very High, High, Low, Very Low)
  normalized_input.1$Potential_FP_or_FN <- ifelse(normalized_input.1[, c(i)] < 10 & normalized_input.1[, c("Normal_Median")] < 10 & normalized_input.1[, c("Gene_State")] != "NORMAL", "+", normalized_input.1$Potential_FP_or_FN)
  #FALSE NEGATIVE: median normal < 10, normalized expression < 10, Not Aberrantly Called (Normal)
  normalized_input.1$Potential_FP_or_FN <- ifelse(normalized_input.1[, c(i)] < 10 & normalized_input.1[, c("Normal_Median")] < 10 & normalized_input.1[, c("Gene_State")] == "NORMAL", "-", normalized_input.1$Potential_FP_or_FN)
  normalized_input.1$Potential_FP_or_FN[is.na(normalized_input.1$Potential_FP_or_FN)] = 0
  rownames(normalized_input.1) = normalized_input.1[,1]
  
  normalized_input.2 = normalized_input.1[, c("SYMBOL", i, "Normalized_Per_Million","Normal_Mean","Normal_Mean_per_Million","Normal_Median", "Normal_SD", "Z_score", "Foldchange", "Gene_State", "Potential_FP_or_FN")]    #selecting columns for final output
  
  write.csv(normalized_input.2, file = paste(i, "-z_stats.csv", sep = ""), row.names = T)  
  
  # Only Actionable Genes.
  a = actionable_input[trimws(actionable_input$gene) %in% trimws(normalized_input.2$SYMBOL),]
   normalized_input.3 = normalized_input.2[normalized_input.2$SYMBOL %in% a$gene, ]
   write.csv(normalized_input.3, paste(i, "-z_stats-actionable.csv", sep = ""))

}

for(i in tumor_sample){    # i is each tumor sample
  main(normalized_data, states_input, gene_counts, actionable_data, i)
}

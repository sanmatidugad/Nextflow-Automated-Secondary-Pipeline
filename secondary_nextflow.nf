#!/usr/bin/env/nextflow

params.input = ""
params.normal_count = ""
params.normal_target = ""
params.tumor_target = ""
params.project = ""

params.outdir = "/home/data/Nextflow_Secondary_Pipeline/"
date = new java.util.Date()

process SUBSAMPLE {

    input:
    path subsample_run
    path subsample_data

    output:
    file "${params.project}_normal_${params.normal_target}.csv"
    file "${params.project}_tumor_${params.tumor_target}.csv"
    file "${params.project}_normal_${params.normal_target}.gc_totals.csv"
    file "${params.project}_tumor_${params.tumor_target}.gc_totals.csv"
    file "${params.project}_data_subsampled.csv"

    """
    Rscript ${subsample_run} \
    -data ${subsample_data} \
    -normal_count ${params.normal_count} \
    -normal_target ${params.normal_target} \
    -tumor_target ${params.tumor_target} \
    -project ${params.project}
    """

}

process NORMALIZATION {
    memory 8.GB
    
    input:
    path normalization_run
    path normalization_data
    path node_modules
    
    output:
    file "${params.project}_normalized.csv"

    """
    node ${normalization_run} ${normalization_data} ${params.project}_normalized.csv
    """
}

process SECONDARY_ANALYSIS {
   memory 8.GB
   
   input:
   path secondary_run
   path normalization_input
   path reference
   
   output:
   file "${params.project}-scaled.csv"
   file "${params.project}-boxstats_threshold.csv"
   file "${params.project}-aberant_gene_states.csv"
   file "${params.project}-gene_states.csv"
   file "${params.project}-gene_statistics.csv"
   file "${params.project}-sample_statistics.csv"
   
   """
   Rscript ${secondary_run} -i ${normalization_input} -r ${reference} -p ${params.project} -n ${params.normal_count}
   """

}

workflow {

SUBSAMPLE("${params.outdir}v2_subsample.R", "${params.outdir}${params.input}")
NORMALIZATION("${params.outdir}normalizer_09292023.js", SUBSAMPLE.out[4], "${params.outdir}node_modules/")
SECONDARY_ANALYSIS("${params.outdir}secondary_analysis.R", NORMALIZATION.out[0], "${params.outdir}ID_SYMBOL.csv")
}



#!/usr/bin/env/nextflow

params.input = ""
params.normal_count = ""
params.target = ""
params.project = ""

params.outdir = "/home/data/Nextflow_Secondary_Pipeline/"
date = new java.util.Date()


process SUBSAMPLE {

    input:
    path subsample_run
    path subsample_data

    output:
    file "${params.project}-subsampled-${params.target}.csv"
    file "${params.project}-subsampled-${params.target}.gc_totals.csv"

    """
    Rscript ${subsample_run} \
    -data ${subsample_data} \
    -normal_count ${params.normal_count} \
    -target ${params.target} \
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

process ZSTATS {
   memory 4.GB
   
   input:
   path zstats_run
   path normalization_input
   path gene_states
   path gene_counts
   path actionable_genes
   
//   output:
//   file "${sample_id}-z_stats.csv"
//   file "${sample_id}-z_stats-actionable.csv"
   
   """
   Rscript ${zstats_run} \
   -i ${normalization_input} \
   -s ${gene_states} \
   -gc ${gene_counts} \
   -n ${params.normal_count} \
   -a ${actionable_genes}
   """

}

workflow {
SUBSAMPLE("${params.outdir}subsample.v1.R", "${params.outdir}${params.input}")
NORMALIZATION("${params.outdir}normalizer.js", SUBSAMPLE.out[0], "${params.outdir}node_modules/")
SECONDARY_ANALYSIS("${params.outdir}secondary_analysis.R", NORMALIZATION.out[0], "${params.outdir}ID_SYMBOL.csv")
ZSTATS("${params.outdir}z_analysis.v1.R", NORMALIZATION.out[0],  SECONDARY_ANALYSIS.out[3], SUBSAMPLE.out[1], "${params.outdir}20230922_genes_edited.csv")
}


// Command-Line
// nextflow secondary_nextflow.v1.nf --input lymph_nodes.csv --normal_count 19 --target 10000000 --project lymph

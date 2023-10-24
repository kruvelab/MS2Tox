#!/bin/bash nextflow


process TRAIN_DOWN_MODEL {
    container = "quay.io/ida_rahu/toxicity_model_r:latest"

    //cpus 8
    //memory '256 GB'
    //time '48h'

    publishDir "${params.outdir}/$assay/", mode: 'move', pattern: "*.rda"
    publishDir "${params.resultsdir}/", mode: 'move', pattern: "*.tsv"


    input:
    tuple val(assay), path(data), path(test_data), val(model), val(correlation)

    output:
    path "*.rda"
    path "*.tsv"

    shell:
    '''
    Rscript !{baseDir}/bin/down_training.R \
     --assay !{assay} \
     --data !{data} \
     --correlation !{correlation} \
     --test_data !{test_data} \
     --model !{model}
    '''
}



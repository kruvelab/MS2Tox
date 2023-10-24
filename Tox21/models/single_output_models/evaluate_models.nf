#!/bin/bash nextflow
nextflow.enable.dsl=2

process EVALUATE_MODELS {
    container = "quay.io/ida_rahu/toxicity_model_r:latest"

    cpus 16
    memory '128 GB'
    time '48h'
        
    publishDir "${params.outdir}/", mode: 'move', pattern: '*.tsv'

    input:
    tuple val(assay), path(test_data), val(model), path(model_file), val(correlation)

    output:
    path '*.tsv'

    shell:
    '''
    Rscript !{baseDir}/bin/evaluate_models.R \
     --assay !{assay} \
     --test_data !{test_data} \
     --model !{model} \
     --model_file !{model_file} \
     --correlation !{correlation} 
    '''
}

workflow {
    test_data_assay_ch = Channel.fromPath(params.test_data).flatten().map {file -> tuple(file.name.split('.tsv')[0].split('_')[1], file) }
    model_files_folder_ch = Channel.fromPath(params.model_files_folder).flatten().map {file -> tuple(file.name.split('_')[-2], file.name.split('_')[-3], file)}
    corr_ch = Channel.from(params.correlation)
    test_data_and_models_ch = test_data_assay_ch.combine(model_files_folder_ch, by: 0)
    final_ch = test_data_and_models_ch.combine(corr_ch)
    EVALUATE_MODELS(final_ch)
}

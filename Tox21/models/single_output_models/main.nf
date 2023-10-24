#!/bin/bash nextflow
nextflow.enable.dsl=2

// Importing modules
include { TRAIN_MODEL } from './modules/train_model'
include {TRAIN_UP_MODEL} from './modules/train_up_model'
include {TRAIN_DOWN_MODEL} from './modules/train_down_model'
include {TRAIN_SMOTE_MODEL} from './modules/train_smote_model'
include {TRAIN_ROSE_MODEL} from './modules/train_rose_model'

workflow {
    
    data_assay_ch = Channel.fromPath(params.data).flatten().map {file -> tuple(file.name.split('.tsv')[0], file) }
    test_data_assay_ch = Channel.fromPath(params.test_data).flatten().map {file -> tuple(file.name.split('.tsv')[0].split('_')[1], file) }
    data_ch = data_assay_ch.join(test_data_assay_ch)
    corr_ch = Channel.from(params.correlation)
    models_ch = Channel.fromList(params.models_to_test)

    combined_ch = data_ch.combine(models_ch)
    combined_ch_final = combined_ch.combine(corr_ch)

    TRAIN_UP_MODEL(combined_ch_final)
    TRAIN_DOWN_MODEL(combined_ch_final)        
    TRAIN_SMOTE_MODEL(combined_ch_final)
    TRAIN_ROSE_MODEL(combined_ch_final)
    TRAIN_MODEL(combined_ch_final)
}

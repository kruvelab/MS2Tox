#!/bin/bash nextflow
nextflow.enable.dsl=2

process TEST_MODEL {
    container = "quay.io/ida_rahu/toxicity_model_r:latest"

    cpus 16
    memory '128 GB'
    time '48h'
        
    publishDir "${params.outdir}/", mode: 'move', pattern: '*.tsv'

    input:
    tuple path(data), path(model), val(cutoff), val(output_file)

    output:
    path '*.tsv'

    shell:
    '''
    Rscript !{baseDir}/test_models.R \
     --data !{data} \
     --model !{model} \
     --cutoff !{cutoff} \
     --output_file !{output_file}
    '''
}

workflow {

    models_ch = Channel.from([
    ["./data/data07/real_life_test_set_07_nr.ahr.tsv", "./trained_models/twoClassSummary/data07/nr.ahr/down_xgbTree_nr.ahr_0.7.rda", 0.524823, 'nr_ahr_xgbTree_07_down.tsv'], 
    ["./data/data07/real_life_test_set_07_nr.ar.lbd.tsv", "./trained_models/twoClassSummary/data07/nr.ar.lbd/smote_Rborist_nr.ar.lbd_0.7.rda", 0.012000, 'nr_ar_lbd_Rborist_07_smote.tsv'],
    ["./NeuralNets/data/data08/real_life_test_set_08_nr.ar.tsv", "./trained_models/twoClassSummary/data08/nr.ar/gbm_nr.ar_0.8.rda", 0.018136, 'nr_ar_gbm_08_original.tsv'],
    ["./data/data07/real_life_test_set_07_nr.aromatase.tsv", "./trained_models/twoClassSummary/data07/nr.aromatase/up_gbm_nr.aromatase_0.7.rda", 0.409162, 'nr_romatase_gbm_07_up.tsv'],
    ["./data/data09/real_life_test_set_09_nr.er.lbd.tsv", "./trained_models/twoClassSummary/data09/nr.er.lbd/rf_nr.er.lbd_0.9.rda", 0.012000, 'nr_er_lbd_rf_09_original.tsv'],
    ["./data/data09/real_life_test_set_09_nr.er.tsv", "./trained_models/twoClassSummary/data09/nr.er/up_rf_nr.er_0.9.rda", 0.062000, 'nr_er_rf_09_up.tsv'],
    ["./data/data07/real_life_test_set_07_nr.ppar.gamma.tsv", "./trained_models/twoClassSummary/data07/nr.ppar.gamma/smote_ranger_nr.ppar.gamma_0.7.rda", 0.076000, 'nr_ppar_gamma_ranger_07_smote.tsv'],
    ["./data/data08/real_life_test_set_08_sr.are.tsv", "./trained_models/twoClassSummary/data08/sr.are/down_xgbDART_sr.are_0.8.rda", 0.317319, 'sr_are_xgbDART_08_down.tsv'],
    ["./data/data08/real_life_test_set_08_sr.atad5.tsv", "./trained_models/twoClassSummary/data08/sr.atad5/up_gbm_sr.atad5_0.8.rda", 0.438668, 'sr_atad5_gbm_08_up.tsv'],
    ["./data/data07/real_life_test_set_07_sr.hse.tsv", "./trained_models/twoClassSummary/data07/sr.hse/rose_gbm_sr.hse_0.7.rda", 0.692024, 'sr_hse_gbm_07_rose.tsv'],
    ["./data/data08/real_life_test_set_08_sr.mmp.tsv", "./trained_models/twoClassSummary/data08/sr.mmp/up_rf_sr.mmp_0.8.rda", 0.208000, 'sr_mmp_rf_08_up.tsv'],
    ["./data/data08/real_life_test_set_08_sr.p53.tsv", "./trained_models/twoClassSummary/data08/sr.p53/down_gbm_sr.p53_0.8.rda", 0.453718, 'sr_p53_gbm_08_down.tsv']])
    

    TEST_MODEL(models_ch)
}

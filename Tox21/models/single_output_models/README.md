This folder contains a Nextflow workflow that was used for training (`main.nf`) and testing (`evaluate_models.nf`) the single-output models. The folder with data generated in `Split_data_train_and_intermediate_test_sets.Rmd` is needed to run the code.

All the trained models are saved as `.rda` files into specified folders (see `nextflow.config`).

While running the evaluation pipeline (for testing models on an intermediate test set), all the relevant statistics about the trained models are written into separate `.tsv` files.

The code in `Select_single_output_models.ipynb` was used to concatenate the relevant statistics about models and select the single-output models for final evaluation (models that showed the lowest false positive rate at 90% of recall on the intermediate test set).

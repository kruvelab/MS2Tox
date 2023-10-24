This folder contains the code that was used for testing the models on a real-life test set.

For testing the single-output models, the Nextflow workflow was used (`main.nf`).
The for-loop was used similarly for testing the multi-output model (see `test_models.R`), and all the predictions were written into one file.

The file `Final_evaluation.ipynb` shows one example of how the different sampling iterations affected the model's outcome

The file `SHAP_analysis.Rmd` contains the code that was used for Shapley Additive exPlanations (SHAP) analysis.

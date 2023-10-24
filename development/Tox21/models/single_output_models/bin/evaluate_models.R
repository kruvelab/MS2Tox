library(optparse)
library(dplyr)
library(data.table)
library(stringr)
library(caret)
library(keras)
library(PRROC)
library(MLmetrics)
library(purrr)
library(pROC)
library(ROSE)
library(performanceEstimation)
library(ROCR)
library(ROCit)

option_list <- list( 
    make_option('--assay', type='character',
    help='One of the twelve Tox21 assay.'),
    make_option('--test_data', type='character',
    help='Data that is used for testing the performance of the model.'),
    make_option('--model', type='character',
    help='Model that was trained.'),
    make_option('--model_file', type='character',
    help='Trained model file .rda'),
    make_option('--correlation', type='character',
    help='Maximum correlation between fingerprint features'));

option_parser <- OptionParser(option_list=option_list)
options <- parse_args(option_parser)
assay <- options$assay
correlation <- options$correlation
test_data <- fread(options$test_data, header=T)
model <- options$model
file_name <- basename(options$model_file)
splitted_length <- length(strsplit(file_name, '_')[[1]])
trained_model <- readRDS(options$model_file)

test_data  <- as.data.frame(test_data)[test_data[, 1] != 99999, ]
test_data[, 1] <- relevel(as.factor(as.factor(test_data[, 1]) %>% make.names()), ref='X1')
colnames(test_data)[1] <- 'tox_class'


predictions <- predict(trained_model, dplyr::select(test_data, -1))
predictions_prob <- predict(trained_model, dplyr::select(test_data, -1), type='prob')

statistics_df <- data.frame(matrix(ncol = 44, nrow = 1))
colnames(statistics_df) <- c('model', 'assay', 'correlation', 'sampling', 'accuracy', 'balanced_accuracy', 'kappa',
                             'F1', 'precision', 'recall', 'PRAUC', 'sensitivity', 'specificity', 'ROCAUC', 
                             'cutoff_09', 'FPR_09', 'accuracy_09', 'balanced_accuracy_09', 'kappa_09', 'F1_09', 
                             'precision_09', 'recall_09', 'sensitivity_09', 'specificity_09',
                             'cutoff_095', 'FPR_095', 'accuracy_095', 'balanced_accuracy_095', 'kappa_095', 'F1_095', 
                             'precision_095', 'recall_095', 'sensitivity_095', 'specificity_095',
                             'cutoff_099', 'FPR_099', 'accuracy_099', 'balanced_accuracy_099', 'kappa_099', 'F1_099', 
                             'precision_099', 'recall_099', 'sensitivity_099', 'specificity_099')

conf_mat <- confusionMatrix(predictions, test_data$tox_class, mode='everything')
pred <- ROCR::prediction(predictions_prob$X1, test_data$tox_class)
perf1 <- ROCR::performance(pred, 'auc')
perf2 <- ROCR::performance(pred, 'aucpr')


roc_emp <- rocit(predictions_prob$X1, relevel(test_data$tox_class, ref='X0'))
threshold1 <- 0.9
cutoff1 <- roc_emp$Cutoff[min(which(roc_emp$TPR >= threshold1))]
FPR1 <- roc_emp$FPR[min(which(roc_emp$TPR >= threshold1))]
pred_labels1 <- ifelse(predictions_prob$X1 >= cutoff1, 'X1', 'X0')
conf_mat1 <- confusionMatrix(relevel(as.factor(pred_labels1), ref='X1'), test_data$tox_class, mode='everything')

statistics_df$cutoff_09 <- cutoff1
statistics_df$FPR_09 <- FPR1
statistics_df$accuracy_09 <- conf_mat1$overall['Accuracy']
statistics_df$balanced_accuracy_09 <- conf_mat1$byClass['Balanced Accuracy']
statistics_df$kappa_09 <- conf_mat1$overall['Kappa']
statistics_df$F1_09 <- conf_mat1$byClass['F1']
statistics_df$precision_09 <- conf_mat1$byClass['Precision']
statistics_df$recall_09 <- conf_mat1$byClass['Recall']
statistics_df$sensitivity_09 <- conf_mat1$byClass['Sensitivity']
statistics_df$specificity_09 <- conf_mat1$byClass['Specificity']

threshold2 <- 0.95
cutoff2 <- roc_emp$Cutoff[min(which(roc_emp$TPR >= threshold2))]
FPR2 <- roc_emp$FPR[min(which(roc_emp$TPR >= threshold2))]
pred_labels2 <- ifelse(predictions_prob$X1 >= cutoff2, 'X1', 'X0')
conf_mat2 <- confusionMatrix(relevel(as.factor(pred_labels2), ref='X1'), test_data$tox_class, mode='everything')

statistics_df$cutoff_095 <- cutoff2
statistics_df$FPR_095 <- FPR2
statistics_df$accuracy_095 <- conf_mat2$overall['Accuracy']
statistics_df$balanced_accuracy_095 <- conf_mat2$byClass['Balanced Accuracy']
statistics_df$kappa_095 <- conf_mat2$overall['Kappa']
statistics_df$F1_095 <- conf_mat2$byClass['F1']
statistics_df$precision_095 <- conf_mat2$byClass['Precision']
statistics_df$recall_095 <- conf_mat2$byClass['Recall']
statistics_df$sensitivity_095 <- conf_mat2$byClass['Sensitivity']
statistics_df$specificity_095 <- conf_mat2$byClass['Specificity']

threshold3 <- 0.99
cutoff3 <- roc_emp$Cutoff[min(which(roc_emp$TPR >= threshold3))]
FPR3 <- roc_emp$FPR[min(which(roc_emp$TPR >= threshold3))]
pred_labels3 <- ifelse(predictions_prob$X1 >= cutoff3, 'X1', 'X0')
conf_mat3 <- confusionMatrix(relevel(as.factor(pred_labels3), ref='X1'), test_data$tox_class, mode='everything')

statistics_df$cutoff_099 <- cutoff3
statistics_df$FPR_099 <- FPR3
statistics_df$accuracy_099 <- conf_mat3$overall['Accuracy']
statistics_df$balanced_accuracy_099 <- conf_mat3$byClass['Balanced Accuracy']
statistics_df$kappa_099 <- conf_mat3$overall['Kappa']
statistics_df$F1_099 <- conf_mat3$byClass['F1']
statistics_df$precision_099 <- conf_mat3$byClass['Precision']
statistics_df$recall_099 <- conf_mat3$byClass['Recall']
statistics_df$sensitivity_099 <- conf_mat3$byClass['Sensitivity']
statistics_df$specificity_099 <- conf_mat3$byClass['Specificity']

statistics_df$model <- model
statistics_df$assay <- assay
statistics_df$correlation <- correlation
statistics_df$sampling <- ifelse(splitted_length ==3, 'NA', strsplit(file_name, '_')[[1]][1]) 
statistics_df$accuracy <- conf_mat$overall['Accuracy']
statistics_df$balanced_accuracy <- conf_mat$byClass['Balanced Accuracy']
statistics_df$kappa <- conf_mat$overall['Kappa']
statistics_df$F1 <- conf_mat$byClass['F1']
statistics_df$precision <- conf_mat$byClass['Precision']
statistics_df$recall <- conf_mat$byClass['Recall']
statistics_df$sensitivity <- conf_mat$byClass['Sensitivity']
statistics_df$specificity <- conf_mat$byClass['Specificity']
statistics_df$PRAUC <- perf2@y.values[[1]]
statistics_df$ROCAUC <- perf1@y.values[[1]]

sampling_to_name <- ifelse(splitted_length ==3, '', strsplit(file_name, '_')[[1]][1]) 
write.table(statistics_df, paste0(sampling_to_name, '_', model, '_', assay, '_', correlation, '.tsv'), row.names=F, col.names=T, sep='\t', quote=F)

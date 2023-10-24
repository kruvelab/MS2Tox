library(dplyr)
library(caret)
library(data.table)
library(ROCit)
library(ROCR)
library(optparse)

option_list <- list( 
  make_option('--data', type='character',
              help='test data'),
  make_option('--model', type='character',
              help='model that is going to be tested'),
  make_option('--cutoff', type='numeric',
              help='cutoff value of the model'),
  make_option('--output_file', type='character',
              help='output file name'))

option_parser <- OptionParser(option_list=option_list)
options <- parse_args(option_parser)

test_data <- as.data.frame(fread(options$data, header=T, sep='\t'))
test_data <- test_data[complete.cases(test_data[, 1]), ]
test_data[, 1] <- relevel(as.factor(as.factor(test_data[, 1]) %>% make.names()), ref='X1')
colnames(test_data)[1] <- 'tox_class'

model_cutoff <- options$cutoff
model <- readRDS(options$model)


sample_func <- function(x) {
  sample(c(0, 1), 1, prob=c(1-x, x))
}

df_predictions_10000 <- data.frame(matrix(ncol=0, nrow=nrow(test_data)))
df_predictions_prob_10000 <- data.frame(matrix(ncol=0, nrow=nrow(test_data)))

for (i in 1:10000) {
  df_sampled <- apply(dplyr::select(test_data, -c(1)), c(1, 2), sample_func)
  predictions <- predict(model, df_sampled)
  df_predictions_10000 <- cbind(df_predictions_10000, as.data.frame(predictions))
  predictions_prob <- predict(model, df_sampled, type='prob')
  df_predictions_prob_10000 <- cbind(df_predictions_prob_10000, as.data.frame(predictions_prob$X1))
}

test_data$average_probs <- rowMeans(df_predictions_prob_10000, na.rm=TRUE)
test_data$pred_labels_default <- ifelse(test_data$average_probs >= 0.5, 'X1', 'X0')
test_data$pred_labels_custom <- ifelse(test_data$average_probs >= model_cutoff, 'X1', 'X0')

write.table(df_predictions_prob_10000, paste0('full_', options$output_file), row.names=F, col.names=T, sep='\t', quote=F)

statistics_df <- data.frame(matrix(ncol = 30, nrow = 1))
colnames(statistics_df) <- c('accuracy', 'balanced_accuracy', 'kappa', 'F1', 'precision', 'recall', 
                             'PRAUC', 'sensitivity', 'specificity', 'ROCAUC', 
                             'accuracy_custom', 'balanced_accuracy_custom', 'kappa_custom', 'F1_custom', 'precision_custom', 'recall_custom', 
                             'PRAUC_custom', 'sensitivity_custom', 'specificity_custom', 'ROCAUC_custom', 
                             'cutoff_09', 'FPR_09', 'accuracy_09', 'balanced_accuracy_09', 'kappa_09', 'F1_09', 
                             'precision_09', 'recall_09', 'sensitivity_09', 'specificity_09')


conf_mat <- confusionMatrix(relevel(as.factor(test_data$pred_labels_default), ref='X1'), test_data$tox_class, mode='everything')
pred <- ROCR::prediction(test_data$average_probs, test_data$tox_class)
perf1 <- ROCR::performance(pred, 'auc')
perf2 <- ROCR::performance(pred, 'aucpr')

roc_emp <- rocit(test_data$average_probs, relevel(test_data$tox_class, ref='X0'))
threshold1 <- 0.9
cutoff1 <- roc_emp$Cutoff[min(which(roc_emp$TPR >= threshold1))]
FPR1 <- roc_emp$FPR[min(which(roc_emp$TPR >= threshold1))]
pred_labels1 <- ifelse(test_data$average_probs >= cutoff1, 'X1', 'X0')
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


conf_mat_custom <- confusionMatrix(relevel(as.factor(test_data$pred_labels_custom), ref='X1'), test_data$tox_class, mode='everything')
statistics_df$accuracy_custom <- conf_mat_custom$overall['Accuracy']
statistics_df$balanced_accuracy_custom <- conf_mat_custom$byClass['Balanced Accuracy']
statistics_df$kappa_custom <- conf_mat_custom$overall['Kappa']
statistics_df$F1_custom <- conf_mat_custom$byClass['F1']
statistics_df$precision_custom <- conf_mat_custom$byClass['Precision']
statistics_df$recall_custom <- conf_mat_custom$byClass['Recall']
statistics_df$sensitivity_custom <- conf_mat_custom$byClass['Sensitivity']
statistics_df$specificity_custom <- conf_mat_custom$byClass['Specificity']

write.table(statistics_df, options$output_file, row.names=F, col.names=T, sep='\t', quote=F)
library(optparse)
library(rcdk)
library(ChemmineR)
library(ChemmineOB)
library(data.table)
library(dplyr)
library(fingerprint)
library(tidyverse)
library(rJava)
library(tidyr)
library(reshape2)
library(stringr)
library(anticlust)
library(caret)
library(naivebayes)
library(caTools)
library(bnclassify)
library(kknn)
library(LiblineaR)
library(ada)
library(plyr)
library(xgboost)
library(LogicReg)
library(e1071)
library(ranger)
library(ordinalForest)
library(adabag)
library(MASS)
library(klaR)
library(gbm)
library(ranger)
library(Rborist)
library(randomForest)
library(ipred)
library(C50)
library(logicFS)
library(neuralnet)
library(nnet)
library(keras)
library(PRROC)
library(MLmetrics)
library(purrr)
library(pROC)
library(ROSE)
library(performanceEstimation)
library(ROCR)

option_list <- list(
    make_option('--assay', type='character',
    help='One of the twelve Tox21 assay.'),
    make_option('--data', type='character',
    help='Data that is used for training the model.'),
    make_option('--correlation', type='character',
    help='Maximum correlation between fingerprint features'),
    make_option('--test_data', type='character',
    help='Data that is used for testing the performance of the model.'),
    make_option('--model', type='character',
    help='Model that is going to be trained.'));

option_parser <- OptionParser(option_list=option_list)
options <- parse_args(option_parser)
assay <- options$assay
data <- fread(options$data, header=T)
correlation <- options$correlation
test_data <- fread(options$test_data, header=T)
model <- options$model

data <- as.data.frame(data)[data[, 1] != 99999, ]
data[, 1] <- relevel(as.factor(as.factor(data[, 1]) %>% make.names()), ref='X1')
colnames(data)[1] <- 'Class'


set.seed(9560)
rose_train <- ROSE(Class ~ ., data=data)$data
rose_train$Class <- relevel(rose_train$Class, ref='X1')

test_data  <- as.data.frame(test_data)[test_data[, 1] != 99999, ]
test_data[, 1] <- relevel(as.factor(as.factor(test_data[, 1]) %>% make.names()), ref='X1')
colnames(test_data)[1] <- 'tox_class'

set.seed(123)

train_control <- trainControl(method='repeatedcv', repeats=10, savePredictions=T, classProbs=T,
                              summaryFunction=twoClassSummary)
 
unlink(".RData")
set.seed(5627)
t <- try(
trained_model <- train(Class ~ ., data = rose_train,
                       method=model, metric='ROC',
                       trControl=train_control,
                       tuneLength=30))

if (inherits(t, 'try-error')) {
  model_NULL <- NULL
  saveRDS(model_NULL, paste0('rose_', model, '_', assay, '_', correlation, '_NULL.rda'))
  fwrite(data.table(), paste0('rose_', model, '_', assay, '_', correlation, '_NULL.tsv'))
  } else {
saveRDS(trained_model, paste0('rose_', model, '_', assay, '_', correlation, '.rda'))

predictions <- predict(trained_model, dplyr::select(test_data, -1))
predictions_prob <- predict(trained_model, dplyr::select(test_data, -1), type='prob')

statistics_df <- data.frame(matrix(ncol = 14, nrow = 1))
colnames(statistics_df) <- c('model', 'assay', 'correlation', 'sampling', 'accuracy', 'balanced_accuracy', 'kappa',
                             'F1', 'precision', 'recall', 'PRAUC', 'sensitivity', 'specificity', 'ROCAUC')

conf_mat <- confusionMatrix(predictions, test_data$tox_class, mode='everything')
pred <- ROCR::prediction(predictions_prob$X1, test_data$tox_class)
perf1 <- ROCR::performance(pred, 'auc')
perf2 <- ROCR::performance(pred, 'aucpr')

statistics_df$model <- model
statistics_df$assay <- assay
statistics_df$correlation <- correlation
statistics_df$sampling <- 'ROSE'
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

write.table(statistics_df, paste0('rose_', model, '_', assay, '_', correlation, '.tsv'), row.names=F, col.names=T, sep='\t', quote=F)}


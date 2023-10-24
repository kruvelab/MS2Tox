FROM bioconductor/bioconductor_docker:RELEASE_3_15
LABEL authors="Ida Rahu" \
      description="Docker image for Tox21 model training in R"

RUN R -e "BiocManager::install(c('optparse', 'rcdk', 'ChemmineR', 'ChemmineOB', 'data.table', 'dplyr', 'fingerprint', 'tidyverse', 'rJava', 'tidyr', 'reshape2', 'stringr', 'anticlust', 'caret', 'naivebayes', 'caTools', 'bnclassify', 'kknn', 'LiblineaR', 'fastAdaboost', 'ada', 'plyr', 'deepboost', 'xgboost', 'LogicReg', 'e1071', 'ranger', 'ordinalForest', 'adabag', 'MASS', 'klaR', 'gbm', 'ranger', 'Rborist', 'randomForest', 'extraTrees', 'ipred', 'C50', 'logicFS', 'neuralnet', 'nnet', 'keras', 'PRROC', 'MLmetrics', 'purrr', 'pROC', 'ROSE', 'performanceEstimation', 'ROCR', 'ROCit', 'tensorflow'))"
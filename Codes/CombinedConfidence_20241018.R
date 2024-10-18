library(tidyverse)
library(MS2Tox)


# #Installing MS2Tox and using unzipping function
# remotes::install_github("kruvelab/MS2Tox")
# 
# MS2Tox::UnZip_SIRIUS5(folderwithSIRIUSdata)


#ConfidenceLevels

# folderwithSIRIUSdata <- " ~/SIRIUSfiles"
# 
# conf_level <- CombinedConfidence(folderwithSIRIUSdata, maxrank = 1)

CombinedConfidence <- function(folderwithSIRIUSdata, maxrank = 1){
  final_scores <- SiriusScores1(folderwithSIRIUSdata) %>%
    left_join(SiriusScores2(folderwithSIRIUSdata),
              by = c("id", "adduct", "molecularFormula")) %>%
    mutate(ConfidenceScore = as.numeric(ConfidenceScore)) %>%
    mutate(topCSIscore = as.numeric(topCSIscore)) %>%
    drop_na(rank)

  final_scores["ConfidenceScore"][final_scores["ConfidenceScore"] == "-Inf"] <- NA

  maxConfScore <- 0

  maxCSIscore <- min((final_scores %>%
                        select(topCSIscore) %>%
                        na.omit())$topCSIscore)

  final_scores$ConfidenceScore[is.na(final_scores$ConfidenceScore)] <- maxConfScore
  final_scores$topCSIscore[is.na(final_scores$topCSIscore)] <- maxCSIscore


  new_confidenceScore <- final_scores %>%
    mutate(ConfidenceScore = as.numeric(ConfidenceScore)) %>%
    mutate(topCSIscore = as.numeric(topCSIscore)) %>%
    mutate(absMassErrorPrecursor = abs(`massErrorPrecursor(ppm)`))

  new_confidenceScore <- new_confidenceScore %>%
    mutate(MassScoreScaled = ((absMassErrorPrecursor-min(absMassErrorPrecursor)))/(max(absMassErrorPrecursor)-min(absMassErrorPrecursor))*4+1) %>%
    mutate(CSIScoreScaled = ((topCSIscore-min(topCSIscore)))/(max(topCSIscore)-min(topCSIscore))*4+1) %>%
    mutate(ConfidenceScoreScaled = ((ConfidenceScore-min(ConfidenceScore)))/(max(ConfidenceScore)-min(ConfidenceScore))*4+1) %>%
    mutate(MassScoreScaled = 6-MassScoreScaled) %>%
    mutate(NewConfidenceScore = (MassScoreScaled+CSIScoreScaled+ConfidenceScoreScaled)/3) %>%
    filter(rank <= maxrank) %>%
    select(id, adduct, molecularFormula, NewConfidenceScore, rank)

  return(new_confidenceScore)

}

SiriusScores1 <- function(folderwithSIRIUSdata) {
  setwd(folderwithSIRIUSdata)
  subfolders_confidenceScore <- dir(folderwithSIRIUSdata, all.files = TRUE, recursive = TRUE, pattern = ".info")
  subfolders_confidenceScore <- subfolders_confidenceScore[grepl("scores", subfolders_confidenceScore)]


  ScoreTable <- tibble()
  for(subfolder in subfolders_confidenceScore){
    subfolder_name_split <- str_split(subfolder, "/")
    #subfolder_name_split <- str_split(subfolder, "_")
    formula_split <- str_split(subfolder_name_split[[1]][3], "_")
    adduct_split <- str_split(formula_split[[1]][2], ".i")

    adduct <- adduct_split[[1]][1]
    #foldernr <- subfolder_name_split[[1]][1]
    id <- subfolder_name_split[[1]][1]
    predform <- formula_split[[1]][1]


    all_scores <- read_delim(subfolder, delim = "\t", col_names = FALSE, show_col_types = FALSE)

    siriusscore <- all_scores %>%
      filter(X1 == "sirius.scores.SiriusScore")

    if (dim(siriusscore)[1] != 0){
      siriusscore <- siriusscore$X2
    } else {
      siriusscore = "missing"
    }

    confidencescore <- all_scores %>%
      filter(X1 == "fingerid.ConfidenceScore")

    if (dim(confidencescore)[1] != 0){
      confidencescore <- confidencescore$X2
    } else{
      confidencescore = "missing"
    }

    isotopescore <- all_scores %>%
      filter(X1 == "sirius.scores.IsotopeScore")

    if (dim(isotopescore)[1] != 0){
      isotopescore <- isotopescore$X2
    } else{
      isotopescore = "missing"
    }

    topCSIscore <- all_scores %>%
      filter(X1 == "fingerid.blast.TopCSIScore")

    if (dim(topCSIscore)[1] != 0){
      topCSIscore <- topCSIscore$X2
    } else{
      topCSIscore = "missing"
    }

    scoretable1 <- tibble(id) %>%
      #mutate(foldernumber = foldernr) %>%
      mutate(adduct = adduct) %>%
      mutate(molecularFormula = predform) %>%
      mutate(SiriusScore = siriusscore) %>%
      mutate(ConfidenceScore = confidencescore) %>%
      mutate(IsotopeScore = isotopescore) %>%
      mutate(topCSIscore = topCSIscore)


    names(scoretable1)[names(scoretable1) == "value"] <- "id"

    ScoreTable <- ScoreTable %>%
      rbind(scoretable1) %>%
      unique()
  }

  ScoreTable1 <- ScoreTable #%>%
    #mutate(id = as.double(id))

  return(ScoreTable1)
}


SiriusScores2 <- function(folderwithSIRIUSdata){
  setwd(folderwithSIRIUSdata)
  subfolders_peaksexplained <- dir(folderwithSIRIUSdata, all.files = TRUE, recursive = TRUE, pattern = "formula_candidates.tsv")
  
  ScoreTable <- tibble()
  for(subfolder in subfolders_peaksexplained){
    subfolder_name_split <- str_split(subfolder, "/")
    
    #foldernr <- subfolder_name_split[[1]][1]
    id <- subfolder_name_split[[1]][1]
    
    
    all_scores <- read_delim(subfolder, delim = "\t", show_col_types = FALSE) 
    
    
    if (length(all_scores) != 0){
      all_scores_notempty <- all_scores %>%
        mutate(id = id) %>%
        #mutate(foldernumber = foldernr) %>%
        select(id, molecularFormula, adduct, SiriusScore, `massErrorPrecursor(ppm)`, rank)
      
      sel_H <- all_scores_notempty %>%
        filter(adduct == "[M + H]+") %>%
        mutate(adduct = "[M+H]+")
      sel_Na <- all_scores_notempty %>%
        filter(adduct == "[M + Na]+") %>%
        mutate(adduct = "[M+Na]+")
      sel_Hneg <- all_scores_notempty %>%
        filter(adduct == "[M - H]-") %>%
        mutate(adduct = "[M-H]-")
      
      all_scores_notempty <- rbind(sel_H, sel_Na, sel_Hneg)
      
      ScoreTable <- ScoreTable %>%
        rbind(all_scores_notempty) %>%
        unique()
      
    } else {
      print(paste("File", subfolder, "missing data"))
    }
  }
  ScoreTable <- ScoreTable #%>%
    #mutate(id = as.double(id))
  
  return(ScoreTable)
  
}



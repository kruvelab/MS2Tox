# library(dplyr)
# library(rlist)
# library(Rdisop)
# library(tibble)
# library(readr)
# library(magrittr)
# library(xgboost)
# library(tidyselect)
# library(stringr)
# library(rcdklibs)
# library(rcdk)


#' @import dplyr
#' @import rlist
#' @import Rdisop
#' @import tibble
#' @import readr
#' @import stringr
#' @import xgboost
#' @import tidyselect
#' @importFrom magrittr %>%



#results_static <- FishLC50Prediction(folderwithSIRIUSfiles, "static")
#Predictions can be made either for static or flow-through mode

#' @export
FishLC50Prediction <- function(folderwithSIRIUSfiles, LC50mode = "static") {
  mode <- LC50mode
  error_message = ""
  if(mode == "static"){
    #FishModel <- readRDS("R/20211118_fish_stat_finalapplic.rds")
    FishModel <- readRDS(system.file("models", "20211118_fish_stat_finalapplic.rds", package = "MS2Tox"))
  } else if(mode == "flow"){
    #FishModel <- readRDS("R/20211126_fish_flow_finalapplic.rds")
    FishModel <- readRDS(system.file("models","20211126_fish_flow_finalapplic.rds", package = "MS2Tox"))
  } else{
    error_message = 'LC50mode must me either "static" or "flow"'
  }
  if (error_message != 'LC50mode must me either "static" or "flow"') {
    inputtable <- FpTableForPredictions(folderwithSIRIUSfiles)
    predictions <- inputtable %>%
      mutate(LC50_predicted = predict(FishModel, newdata = inputtable)) %>%
      select(id, foldernumber, LC50_predicted, predform, predion) %>%
      unique()
    #return(predictions)

  } else {
    print('LC50mode must me either static or flow')
    predictions = "LC50mode must me either static or flow"
  }
  return(predictions)
}

#' @export
FpTableForPredictions <- function(folderwithSIRIUSfiles){
  subfolder <- dir(folderwithSIRIUSfiles, all.files = TRUE, recursive = TRUE, pattern = ".fpt")
  subfolder_score <- dir(folderwithSIRIUSfiles, all.files = TRUE, recursive = TRUE, pattern = ".info")

  fp_names_pos <- paste("Un", read_delim(paste(folderwithSIRIUSfiles,"/csi_fingerid.tsv", sep = ""), delim = "\t")$absoluteIndex, sep = "")
  fp_names_neg <- paste("Un", read_delim(paste(folderwithSIRIUSfiles,"/csi_fingerid_neg.tsv", sep = "" ), delim = "\t")$absoluteIndex, sep = "")
  fp_names_common <- as.data.frame(fp_names_pos)  %>%
    mutate(fp_names_neg = fp_names_pos) %>%
    inner_join(as.data.frame(fp_names_neg))
  fp_names_common <- as.vector(fp_names_common$fp_names_pos)

  selected_rank1_table <- FingerPrintTable(subfolder, fp_names_pos, fp_names_neg, fp_names_common, folderwithSIRIUSfiles) %>%
    inner_join(SiriusScoreRank1(subfolder_score, folderwithSIRIUSfiles), by = c("id", "foldernumber", "predion"))

  final_table_mass <- selected_rank1_table %>%
    mutate(exactMass = "mass")

  for(n in 1:length(final_table_mass$predform)){
    molecule = getMolecule(final_table_mass$predform[[n]])
    final_table_mass$exactMass[[n]] = getMass(molecule)
  }

  final_table_mass <- final_table_mass %>%
    mutate(exactMass = as.numeric(exactMass)) %>%
    select(id, foldernumber, predform, predion, exactMass, everything())

  return(final_table_mass)
}

#' @export
FingerPrintTable <- function(subfolder, fp_names_pos, fp_names_neg, fp_names_common, folderwithSIRIUSfiles){

  subfolder_pos <- c()
  subfolder_neg <- c()
  for(subfold in subfolder){
    if ((grepl("]2+", subfold, fixed=TRUE) | grepl("]+", subfold, fixed=TRUE)) & grepl("/fingerprint", subfold, fixed=TRUE)){
      subfolder_pos <- subfolder_pos %>%
        list.append(subfold)
    }

    if ((grepl("]2-", subfold, fixed=TRUE) | grepl("]-", subfold, fixed=TRUE)) & grepl("/fingerprint", subfold, fixed=TRUE)){
      subfolder_neg <- subfolder_neg %>%
        list.append(subfold)
    }
  }

  fp_pos <- FingerPrintTablePOS(subfolder_pos, folderwithSIRIUSfiles)

  if(nrow(fp_pos) != 0){
    colnames(fp_pos) <- c(fp_names_pos, "predion", "id", "foldernumber", "predform")
  }

  fp_neg <- FingerPrintTableNEG(subfolder_neg, folderwithSIRIUSfiles)

  if(nrow(fp_neg) != 0){
    colnames(fp_neg) <- c(fp_names_neg, "predion", "id", "foldernumber", "predform")
  }

  if (!is.null(subfolder_neg) & !is.null(subfolder_pos)){

  fp_all <- fp_pos %>%
    select(id, foldernumber, predform, predion, all_of(fp_names_common)) %>%
    bind_rows(fp_neg %>% select(id, foldernumber, predform, predion, all_of(fp_names_common))) %>%
    unique()
  }


  if (is.null(subfolder_neg) & !is.null(subfolder_pos)){

    fp_all <- fp_pos %>%
      select(id, foldernumber, predform, predion, all_of(fp_names_common)) %>%
      unique()
  }

  if (!is.null(subfolder_neg) & is.null(subfolder_pos)){

    fp_all <- fp_neg %>%
      select(id, foldernumber, predform, predion, all_of(fp_names_common)) %>%
      unique()
  }

  if (is.null(subfolder_neg) & is.null(subfolder_pos)){
    print("NO data found in subfolders")
  }

  return(fp_all)
}

#' @export
FingerPrintTablePOS <- function(subfolder, folderwithSIRIUSfiles){
  fingerprint_data <- tibble()
  for(direct in subfolder){
# subfolder with data must be on a form where id is on second place after nr_ (0_Einc270001_Einc270001)
    file_name <- str_split(direct, "_")
    comp_name <- str_split(direct, "/")

    sir_fold <- file_name[[1]][1]
    id_this <- file_name[[1]][2] # if id in some other position, this can be changed
    pred_ion <- comp_name[[1]][3]

    filedata <- read_delim(paste(folderwithSIRIUSfiles, direct, sep = "/"), delim = " ", col_names = FALSE)
    filedata <- as_tibble(t(filedata))
    filedata <- filedata %>%
      mutate(predion = pred_ion) %>%
      mutate(predion = sub("\\..*", "", predion)) %>%
      mutate(id = id_this) %>%
      mutate(sir_fol_nr = sir_fold) %>%
      mutate(predform = sub("\\_.*", "", predion))
    fingerprint_data <- fingerprint_data %>%
      bind_rows(filedata)
  }
  return(fingerprint_data)
}

#' @export
FingerPrintTableNEG <- function(subfolder, folderwithSIRIUSfiles){
  fingerprint_data <- tibble()
  for(direct in subfolder){
# subfolder with data must be on a form where id is on second place after nr_ (0_Einc270001_Einc270001)
      file_name <- str_split(direct, "_")
      comp_name <- str_split(direct, "/")

      sir_fold <- file_name[[1]][1]
      id_this <- file_name[[1]][2]
      pred_ion <- comp_name[[1]][3]

      filedata <- read_delim(paste(folderwithSIRIUSfiles, direct, sep = "/"), delim = " ", col_names = FALSE)
      filedata <- as_tibble(t(filedata))
      filedata <- filedata %>%
        mutate(predion = pred_ion) %>%
        mutate(predion = sub("\\..*", "", predion)) %>%
        mutate(id = id_this) %>%
        mutate(sir_fol_nr = sir_fold) %>%
        mutate(predform = sub("\\_.*", "", predion))
      fingerprint_data <- fingerprint_data %>%
        bind_rows(filedata)

  }
  return(fingerprint_data)
}

#' @export
SiriusScoreRank1 <- function(subfolder_score, folderwithSIRIUSfiles){
  scores_table <- tibble()
  scores_table <- scores_table %>%
    mutate(id = "id") %>%
    mutate(foldernumber = "sirfolnr") %>%
    mutate(siriusscore = "siriusscore")

  for (filename in subfolder_score) {
    if (grepl("/scores", filename, fixed=TRUE)){

      file_name_score <- str_split(filename, "_")
      comp_name_score <- str_split(filename, "/")

      foldernumber <- file_name_score[[1]][1]
      id <- file_name_score[[1]][2] #if id is not in second place, [2] must be changed
      pred_st <- comp_name_score[[1]][3]

      fileConnection <- file(paste(folderwithSIRIUSfiles, filename, sep = "/"))
      record <- readLines(fileConnection)
      close(fileConnection)

      SiriusScore <- substring(grep('sirius.scores.SiriusScore', record, value = TRUE, fixed = TRUE),27)

      filedata <- data.frame(id , foldernumber)

      filedata <- filedata %>%
        mutate(predion = pred_st) %>%
        mutate(predion = sub("\\..*", "", predion)) %>%
        mutate(siriusscore = SiriusScore)

      scores_table <- scores_table%>%
        bind_rows(filedata)
    }
  }
  data_scores <- scores_table %>%
    unique() %>%
    select(-predion) %>% #if two compound have same siriusscore they get same rank
    unique() %>%
    mutate(siriusscore = as.numeric(siriusscore)) %>%
    group_by(id, foldernumber) %>%
    arrange(desc(siriusscore)) %>%
    mutate(rank = row_number()) %>%
    ungroup() %>%
    left_join(scores_table %>% unique() %>%  mutate(siriusscore = as.numeric(siriusscore)), by = c("id", "foldernumber", "siriusscore")) %>%
    filter(rank == 1) %>%
    select(-rank, -siriusscore) %>%
    unique()

  return(data_scores)
}






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
#' @import rcdklibs
#' @import rcdk
#' @import stats
#' @import utils
#' @importFrom magrittr %>%



#results_static <- FishLC50Prediction(folderwithSIRIUSfiles, "static")
#Predictions can be made either for static or flow-through mode

#' @export
FishLC50Prediction <- function(fingerprints, LC50mode = "static") {

  mode <- LC50mode
  error_message = ""
  if(mode == "static"){
    #FishModel <- readRDS("R/20211118_fish_stat_finalapplic.rds")
    FishModel <- readRDS(system.file("models", "20211118_fish_stat_finalapplic.rds", package = "MS2Tox"))
  } else if(mode == "flow"){
    #FishModel <- readRDS("R/20211126_fish_flow_finalapplic.rds")
    FishModel <- readRDS(system.file("models","20211126_fish_flow_finalapplic.rds", package = "MS2Tox"))
  } else{
    error_message = 'LC50mode must be either "static" or "flow"'
  }
  if (error_message != 'LC50mode must be either "static" or "flow"') {
    if (is.character(fingerprints))
      inputtable <- FpTableForPredictions(fingerprints)
    else
      inputtable <- fingerprints
    predictions <- inputtable %>%
      mutate(LC50_predicted = predict(FishModel, newdata = inputtable)) %>%
      select(id, foldernumber, LC50_predicted, predform, predion) %>%
      unique()
    #return(predictions)

  } else {
    print('LC50mode must be either static or flow')
    predictions = "LC50mode must be either static or flow"
  }
  return(predictions)
}


#' @export
FpTableForPredictions <- function(folderwithSIRIUSfiles){
  subfolder <- dir(folderwithSIRIUSfiles, all.files = TRUE, recursive = TRUE, pattern = ".fpt")
  subfolder_score <- dir(folderwithSIRIUSfiles, all.files = TRUE, recursive = TRUE, pattern = ".info")

  fp_names_pos <- paste("Un", read_delim(paste(folderwithSIRIUSfiles,"/csi_fingerid.tsv", sep = ""), delim = "\t", show_col_types = FALSE)$absoluteIndex, sep = "")
  fp_names_neg <- paste("Un", read_delim(paste(folderwithSIRIUSfiles,"/csi_fingerid_neg.tsv", sep = "" ), delim = "\t", show_col_types = FALSE)$absoluteIndex, sep = "")
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

  # options for progress bar
  options(width = 80)
  n <- length(subfolder)
  
  for(direct in subfolder){
# subfolder with data must be on a form where id is on second place after nr_ (0_Einc270001_Einc270001)
    file_name <- str_split(direct, "_")
    comp_name <- str_split(direct, "/")

    sir_fold <- file_name[[1]][1]

    # detect patRoon featureID and parse to output
    if (str_detect(file_name[[1]][2], "M") & str_detect(file_name[[1]][3], "R")) {
      id_this <- str_c(file_name[[1]][2], file_name[[1]][3], file_name[[1]][4], sep = "_")
    } else id_this <- file_name[[1]][2]
    
    pred_ion <- comp_name[[1]][3]

    filedata <- read_delim(paste(folderwithSIRIUSfiles, direct, sep = "/"), delim = " ", col_names = FALSE, show_col_types = FALSE)
    filedata <- as_tibble(t(filedata))
    filedata <- filedata %>%
      mutate(predion = pred_ion) %>%
      mutate(predion = sub("\\..*", "", predion)) %>%
      mutate(id = id_this) %>%
      mutate(sir_fol_nr = sir_fold) %>%
      mutate(predform = sub("\\_.*", "", predion))
    fingerprint_data <- fingerprint_data %>%
      bind_rows(filedata)

    # add progress bar
    extra <- nchar('||100%')
    width <- options()$width
    step <- round(direct / n * (width - extra))
    text <- sprintf('|%s%s|% 3s%%', strrep('=', step),
                    strrep(' ', width - step - extra), round(direct / n * 100))
    cat(text)
    Sys.sleep(0.05)
    cat(if (direct == n) '\n' else '\014')
      
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
    
      # detect patRoon featureID and parse to output
      if (str_detect(file_name[[1]][2], "M") & str_detect(file_name[[1]][3], "R")) {
        id_this <- str_c(file_name[[1]][2], file_name[[1]][3], file_name[[1]][4], sep = "_")
      } else id_this <- file_name[[1]][2]
    
      pred_ion <- comp_name[[1]][3]

      filedata <- read_delim(paste(folderwithSIRIUSfiles, direct, sep = "/"), delim = " ", col_names = FALSE, show_col_types = FALSE)
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
    dplyr::filter(rank == 1) %>%
    select(-rank, -siriusscore) %>%
    unique()

  return(data_scores)
}





#' @export
LC50fromSMILES <- function(compoundslistwithSMILES, LC50mode = "static"){

  mode <- LC50mode
  error_message = ""
  if(mode == "static"){
    #FishModel <- readRDS("R/20211118_fish_stat_finalapplic.rds")
    FishModel <- readRDS(system.file("models", "20211118_fish_stat_finalapplic.rds", package = "MS2Tox"))
  } else if(mode == "flow"){
    #FishModel <- readRDS("R/20211126_fish_flow_finalapplic.rds")
    FishModel <- readRDS(system.file("models","20211126_fish_flow_finalapplic.rds", package = "MS2Tox"))
  } else{
    error_message = 'LC50mode must be either "static" or "flow"'
  }
  if (error_message != 'LC50mode must be either "static" or "flow"') {

    inputtable <- SMILESFingerprints(compoundslistwithSMILES)
    predictions <- inputtable %>%
      mutate(LC50_predicted = predict(FishModel, newdata = inputtable)) %>%
      select(SMILES, LC50_predicted) %>%
      unique()

  } else {
    print('LC50mode must be either static or flow')
    predictions = "LC50mode must be either static or flow"
  }
  return(predictions)
}




#' @export
SMILESFingerprints <- function(compoundslistwithSMILES){

  #substructure fingerprints
  row <- c(1:308) #first one has mass, 307 fingerprints actually
  list <- as.data.frame(row)

  col_names_substr <- c(paste("Un", c(55:361), sep = ""), "exactMass")

  substr_table <- as.data.frame(t(list)) %>%
    rownames_to_column()

  error_comp <- tibble()

  for (n in 1:length(compoundslistwithSMILES$SMILES)){

    SMILES <- compoundslistwithSMILES$SMILES[[n]]
    exactmass <- as.numeric(compoundslistwithSMILES$exactMass[[n]] )

    mol2 <- parse.smiles(SMILES)[[1]]

    if(!is.null(mol2)){

      substr_fingerprints <- get.fingerprint(mol2,
                                             type = "substructure")

      table <- as.data.frame(substr_fingerprints@bits) %>%
        mutate(fp = 1) %>%
        mutate(fp = as.numeric(fp))
      colnames(table) <- c("row", SMILES)

      table <- table %>%
        mutate(row = as.numeric(row))

      datarow <- list %>%
        left_join(table) %>%
        select(-row)

      newrow <- as.data.frame(t(datarow)) %>%
        rownames_to_column() %>%
        mutate(V308 = exactmass)

      substr_table <- substr_table %>%
        add_row(newrow)
    }

    else{
      wrong_smile <- as.data.frame(SMILES)
      error_comp <- error_comp %>%
        rbind(wrong_smile)
      print(SMILES)
    }
  }

  substr_table[is.na(substr_table)] <- 0

  substr_table <- substr_table %>%
    dplyr::filter(rowname != "row")

  names(substr_table) <- c("SMILES", col_names_substr)

  #MACCS fingerprints
  row <- c(1:166)

  list <- as.data.frame(row)

  col_names_maccs <- paste("Un", c(362:527), sep = "")

  maccs_table <- as.data.frame(t(list)) %>%
    rownames_to_column()

  for (n in 1:length(compoundslistwithSMILES$SMILES)){

    SMILES <- compoundslistwithSMILES$SMILES[[n]]

    mol2 <- parse.smiles(SMILES)[[1]]

    if(!is.null(mol2)){

      maccs_fingerprints <- get.fingerprint(mol2,
                                            type = "maccs")

      table <- as.data.frame(maccs_fingerprints@bits) %>%
        mutate(fp = 1) %>%
        mutate(fp = as.numeric(fp))
      colnames(table) <- c("row", SMILES)

      table <- table %>%
        mutate(row = as.numeric(row))

      datarow <- list %>%
        left_join(table) %>%
        select(-row)

      newrow <- as.data.frame(t(datarow)) %>%
        rownames_to_column()

      maccs_table <- maccs_table %>%
        add_row(newrow)
    }

    else{
      print(SMILES)
    }
  }

  maccs_table[is.na(maccs_table)] <- 0

  maccs_table <- maccs_table %>%
    dplyr::filter(rowname != "row")

  names(maccs_table) <- c("SMILES", col_names_maccs)


  #pubchem fingerprints
  row <- c(1:881)

  list <- as.data.frame(row)

  col_names_pubchem <- paste("Un", c(528:1408), sep = "")

  pubchem_table <- as.data.frame(t(list)) %>%
    rownames_to_column()

  for (n in 1:length(compoundslistwithSMILES$SMILES)){

    SMILES <- compoundslistwithSMILES$SMILES[[n]]

    mol2 <- parse.smiles(SMILES)[[1]]

    if(!is.null(mol2)){

      pubchem_fingerprints <- get.fingerprint(mol2,
                                              type = "pubchem")

      table <- as.data.frame(pubchem_fingerprints@bits) %>%
        mutate(fp = 1) %>%
        mutate(fp = as.numeric(fp))
      colnames(table) <- c("row", SMILES)

      table <- table %>%
        mutate(row = as.numeric(row))

      datarow <- list %>%
        left_join(table) %>%
        select(-row)

      newrow <- as.data.frame(t(datarow)) %>%
        rownames_to_column()

      pubchem_table <- pubchem_table %>%
        add_row(newrow)
    }

    else{
      print(SMILES)
    }
  }

  pubchem_table[is.na(pubchem_table)] <- 0

  pubchem_table <- pubchem_table %>%
    dplyr::filter(rowname != "row")

  names(pubchem_table) <- c("SMILES", col_names_pubchem)

  #KlekRoth fingerprints

  row <- c(as.numeric(1:4860))

  list <- as.data.frame(row)

  col_names_KlekotaRoth <- paste("Un", c(1409:6268), sep = "")


  KlekotaRoth_table <- as.data.frame(t(list)) %>%
    rownames_to_column()

  for (n in 1:length(compoundslistwithSMILES$SMILES)){

    SMILES <- compoundslistwithSMILES$SMILES[[n]]

    mol2 <- parse.smiles(SMILES)[[1]]

    if(!is.null(mol2)){

      KlekotaRoth_fingerprints <- get.fingerprint(mol2,
                                                  type = "kr")

      table <- as.data.frame(KlekotaRoth_fingerprints@bits) %>%
        mutate(fp = 1) %>%
        mutate(fp = as.numeric(fp))
      colnames(table) <- c("row", SMILES)

      table <- table %>%
        mutate(row = as.numeric(row))

      datarow <- list %>%
        left_join(table) %>%
        select(-row)

      newrow <- as.data.frame(t(datarow)) %>%
        rownames_to_column()

      KlekotaRoth_table <- KlekotaRoth_table %>%
        add_row(newrow)
    }

    else{
      print(SMILES)
    }
  }

  KlekotaRoth_table[is.na(KlekotaRoth_table)] <- 0

  KlekotaRoth_table <- KlekotaRoth_table %>%
    dplyr::filter(rowname != "row")

  names(KlekotaRoth_table) <- c("SMILES", col_names_KlekotaRoth)


  #custommadeSMARTS
  row <- c(1:202)
  list <- as.data.frame(row)

  smartslist <- c("C1(CCC3C2CCCC3)CCCCC12","CC(C)CC(C(~O)[O,N])N","CC(C(C(~O)[O,N])N)O",
                  "[OH0]=[CH0](O)[CH1][CH2][OH1]","[OH0]=[CH0](N)[CH1](N)[CH2][OH1]","[OH0]=[CH0](N)[CH1]N",
                  "[CH2]([CH1]([CH2][OH0][PH0](~O)(~O)~O)~O)~O","[CH2]([CH1]([CH1]([CH1]([CH1]1~[#7])~O)~O)[OH0]1)[OH0][PH0](~O)(~O)~O",
                  "[#6]~C([CH0]([NH1][CH1](~[#6])[CH0](~[!#1])~O)~[OH0D1])~[!#1]",
                  "[CH1](~[!#1])[CH1]([CH1]([CH1]([CH1]([CH2]~O)~O)~O)~O)[NH1][CH0](~[CH3D1])~[OH0D1]","c(:c:c:n1[CH2][CH2]~C):c1",
                  "[CH2](c(:c:n:c1):n1)[CH1]([CH0](~[!#1])~O)~N","[CH2]([CH2]~[!#1])[NH0](~[#6])~[#6]","[#6]~C([NH1][CH1]([CH2]c(:c:n:c1:c:c:c:c2):c12)[CH0](~[!#1])~O)~[OH0D1]",
                  "[#6]~S[CH2][CH2][CH1]([CH0](~[!#1])~O)~N","[#6]~C([NH1][CH1]([CH2][CH2][CH0](~[!#1])~O)[CH0](~[!#1])~O)~[OH0D1]",
                  "[#6]~Cc(:c:c:c(:c1~[!#1])~[!#1]):c1","[#6]~C([CH0](~[!#1])~[!#1])[NH0](~[#6])[CH0](~[!#1])~O","c(:c(~[!#1]):n:c(:n1)~[!#1])(:c1~[!#1])~[!#1]",
                  "[#6]~C[CH0](~[!#1])[OH0][CH1](~[#6])~[#6]","[#6]~C[CH2][CH2][CH2][CH0]([NH1]~[#6])~[OH0D1]",
                  "[#6]~C[OH0][CH1]([CH1]([CH1]([CH1]([CH2]1)~O)~O)~O)[OH0]1","[CH1]([CH0](~C)(~[#6])~[!#1])(~C)[OH0][CH1](~C)~[!#1]",
                  "[#6]~C(~[!#1])[CH2][CH1]([OH0]~[#6])[OH0]~[#6]","[#6]~C([CH0](~[!#1])~O)[NH1][CH0]([CH2]~[!#1])~[OH0D1]","[#6]~n(:c:c:n1):c1",
                  "[#6]~Oc(:c:c:c:c1~[!#1]):c1~[!#1]","[#6]~C(~C)[OH0][PH0](~O)(~O)~O","c(:c:c:c(:c1)[NH0](~[!#1])~[!#1]):c1",
                  "[#6]~C([CH1](~[#6])~[!#1])[NH1][CH0](~[#6])~[OH0D1]","[CH2](~[#6])[NH1][CH0](~[!#1])~[!#1]","[#6]~O[CH1]([CH0](~C)(~[#6])~[!#1])~[#6]",
                  "[#6]~O[CH1]([CH2]~[!#1])[CH2]~[!#1]","[#6]~C(~[!#1])[NH0](~[#6])~[#6]","[CH2]([CH1]([CH1]([CH1]([CH1]1~[!#1])~O)~O)[OH0]1)~O",
                  "[#6]~Cc(:c:c:c(:c1)~[!#1]):c1~[!#1]","[#6]~Cc(:c:c:c:c1~[!#1]):c1~[!#1]","[#6]~C([CH1]([CH1]([CH1]([CH1]1[OH0]c(:c:c:c:c2):c2)~O)~O)~O)[OH0]1",
                  "[CH2]([CH0](~[!#1])~[!#1])[OH0]~[#6]","[#6]~C(~[#6])[CH0](~[!#1])[OH0]~[#6]","c(:c:c:c(:c1)[CH0](~[!#1])(~[!#1])~[!#1]):c1",
                  "[#6]~C[OH0]c(:c:c:c(:c1)~[!#1]):c1","[#6]~C=[CH1][CH2][CH2][CH2][CH0](~[!#1])~O","[#6]~Oc(:c:c(:c:c1)[CH2]~[!#1]):c1~[!#1]",
                  "[#6]~N[CH1]([CH2]~[#6])[CH0](~[!#1])~O","[#6]~C[CH2][CH0](~[!#1])[OH0][CH1](~[#6])~[#6]","[#6]~C[CH2][OH0][CH0](~[!#1])~[#6]",
                  "[#6]~C([CH2][OH0][PH0](~O)(~O)~O)~[!#1]","[#6]~C([CH0](~[!#1])~O)[NH1][CH0]([CH1]([CH1](~[#6])~[!#1])~[!#1])~[OH0D1]",
                  "c(:c:c(:c([CH0]~[!#1]):c1~[!#1])~[!#1]):c1","[#6]~Cc(:c:c:c(:c1[CH2][CH1]=[CH0](~[#6])~[#6])~[#8]):c1","c(:c:c(:c:c1)~[!#1]):c1[CH1]~[!#1]",
                  "[#6]~C[CH2][CH2][CH2][CH2][CH2][CH2][CH2][CH0](~[OH0D1])[OH0][CH1](~[#6])~[#6]","[#6]~O[CH0](~[!#1])[CH2]~[!#1]",
                  "[#6]~O[CH1]([CH1]([CH1]([CH1]([CH1]1[CH2]~O)~[OH1D1])~[OH1D1])~O)[OH0]1","[#6]~N[CH2][CH2][CH2][CH1]~[!#1]",
                  "c(:c:c(:c:c1)~[!#1]):c1[CH2]~[!#1]","[#6]~C([CH1]([CH0](~C)~[!#1])~[!#1])[OH0]~[#6]","[CH2]([CH2][NH1][CH0](~[!#1])~[OH0D1])[CH0](~[!#1])~O",
                  "[#6]~C([CH2][CH2][CH1]([CH0](~[#6])(~[#6])[CH1]([CH2][CH2]1)~O)[CH0]12~[#6])[CH1]2[CH2]~[#6]","[#6]~Cc(:c:c(:c:c1)~[!#1]):c1~[!#1]",
                  "[CH2]([CH0]([CH2][OH0][CH1]1~[OH0D2])([CH1]1~[OH1D1])~[OH1D1])~[OH1D1]","[#6]~Oc(:c:c:c(:c1)~[!#1]):c1~[!#1]",
                  "[#6]~N[CH2][CH1]([CH0](~[!#1])~[!#1])~[#6]","[CH2]([CH2][CH1]([CH0]1(~[CH3D1])~C)~C)[CH1]1~[!#1]","[#6]~O[CH0](~[!#1])[CH1](~[#6])~[!#1]",
                  "[#6]~C[OH0][PH0](~O)(~[!#1])[OH0][CH2]~[#6]","[#6]~c(:c(:c:c(:c1)[OH0]~[#6])~[!#1]):c1~[!#1]","[#6]~C([CH0](~[!#1])~[!#1])[OH0][CH0](~[!#1])~[#6]",
                  "[#6]~C[CH2][NH0](~[#6])~[!#1]","[CH1](~[!#1])[CH1]([CH1]([CH1]([CH1]([CH2]~[!#1])~O)~O)~O)~[!#1]","[#6]~C[CH1]([CH2]~[#6])[OH0]~[#6]",
                  "[#6]~C([CH2]~[#6])[NH0](~[#6])~[#6]","[#6]~C[CH1](~[#6])[CH0](~[#6])[CH0](~[!#1])[OH0]~[#6]","[#6]~C[OH0][CH1](~[#6])[OH0][CH2]~[#6]",
                  "[#6]~C[CH2][CH0](~O)[OH0][CH2][CH1]~[#6]","[#6]~c(:c:c(:c(:c1)~[!#1])~[!#1]):c1~[!#1]","[#6]~c(:c:c:c:c1):c1[OH0]~[#6]",
                  "[#6]~c(:c:c(:c:c1~[!#1])~[!#1]):c1","[#6]~C[CH1]=[CH1][CH2][CH2][CH0](~[!#1])~O","[#6]~C(~[#6])(~[#6])[OH0]~[#6]",
                  "[#6]~N[CH0](~[!#1])[NH1]~[#6]","c(:c:c(:c:c1)~[!#1]):c1[CH0](~[!#1])~[!#1]",
                  "[#6]~C([CH1]([CH1]([CH1]([CH1]1~[OH0D2])~O)[OH0]~[#6])~O)[OH0]1","[#6]~C[OH0][CH1](~[#6])[OH0][CH1](~[#6])[CH0](~[!#1])~[!#1]",
                  "[#6]~c(:c:c:c(:c1[OH0]~[#6])~[!#1]):c1","[#6]~c(:c:c:c:c1):c1[OH0]c(:c:c:c:c2):c2","[#6]~c(:c:c:c:c1[CH2]~[!#1]):c1",
                  "c(:c(:c:c(:c1[CH1]~[!#1])~[!#1])~[!#1]):c1~[!#1]","[#6]~C([CH0]([NH1][CH1]([CH2]c(:c:c:c:c1):c1)[CH0](~[!#1])~O)~[OH0D1])~[!#1]",
                  "c(:c:c:c(:c1)[CH1](~[#6])~[!#1]):c1","c(:c:c(:c:c1~[!#1])[CH0]~[!#1]):c1","[#6]~C(~[!#1])[CH1](~[#6])[OH0][CH0](~[!#1])~[#6]",
                  "[#6]~Sc(:c:c:c:c1):c1","[CH2]([CH2][OH0][PH0](~O)~O)~[!#1]","[#6]~C[OH0][CH1]([CH0](~C)(~[#6])~[!#1])~[!#1]",
                  "[#6]~C(c(:c:c:c:c1~[!#1]):c1)~[!#1]","c(:c:c(:c:c1~[!#1])[CH2]~[!#1]):c1","[#6]~C[CH2][CH2][CH2][CH2][CH2][CH2][CH2][CH2][CH0](~O)[OH0]~[#6]",
                  "[#6]~C[CH2][CH1]([CH0](~[!#1])~O)[NH1][CH0](~[#6])~[OH0D1]","[#6]~c(:c:c:c(:c1~[!#1])~[!#1]):c1~[!#1]","[#6]~C[OH0][CH0](~[!#1])[CH1]~[!#1]",
                  "c(:c:c(:c(:c1)[NH0](~[!#1])~[!#1])~[!#1]):c1~[!#1]","[CH2]([CH1]([CH0](~[!#1])~O)[NH1][CH0]([CH1]([CH2]~[#6])~[#7])~[OH0D1])~[!#1]",
                  "[#6]~c(:c:c:c:n1~[!#1]):c1","c(:c([CH2]~[!#1]):c(:c:c1~[!#1])~[!#1]):c1~[!#1]","[#6]~Nc(:c:c:c:c1):c1",
                  "[#6]~C[CH0](~[#6])([CH0](~[!#1])[OH0]~[#6])~[!#1]","[#6]~O[CH0](~[!#1])[CH2][CH2]~[#6]",
                  "c(:c:c:c(:c1[CH1]~[!#1])~[!#1]):c1","[#6]~C(=[CH1]~[!#1])[CH0](~[!#1])[OH0]~[#6]",
                  "[#6]~c(:c(:c:c:c1~[!#1])~[!#1]):c1[OH0]~[#6]","[#6]~C([CH1]([CH1]([CH1]([CH1]~[!#1])~[!#1])~[!#1])[OH0]~[#6])~[!#1]",
                  "[CH2]([CH1]([CH0]([NH1]~[#6])~[OH0D1])~[!#1])~[#6]","[#6]~C[NH1][CH2][CH2]~[!#1]","[#6]~C[CH1]([CH1](~[#6])[CH0](~[!#1])[OH0]~[#6])~[!#1]",
                  "[#6]~C([NH1][CH1]([CH2]~[!#1])[CH0](~[!#1])~O)~O","[#6]~c(:c:c(:c:c1)[OH0]~[#6]):c1~[!#1]",
                  "[#6]~C(c(:c:c:c(:c1)~[!#1]):c1)~[!#1]","[#6]~C([CH0](~[!#1])~O)[NH1][CH0]([CH1]([CH2]c(:c:c:c:c1):c1)~[#7])~[OH0D1]",
                  "[#6]~c(:c:c:c(:c1)[CH2]~[!#1]):c1~[!#1]","[#6]~c(:c:c:c:c1):c1[CH2]~[!#1]","c(:c:c:c(:c1)[OH0][CH1]~[!#1]):c1",
                  "c(:c(:c:c(:c1~[!#1])~[!#1])[CH0]~[!#1]):c1~[!#1]","[#6]~C(~[#6])[CH0](~[!#1])[OH0][CH1](~[#6])~[#6]","[#6]~C([CH1]([CH1]([CH1]([CH1]1~O)[OH0]~[#6])~O)~O)[OH0]1",
                  "[#6]~Oc(:c:c(:c:c1~[!#1])~[!#1]):c1","[#6]~C[OH0][CH0](~[!#1])[CH2][CH2]~[#6]","[#6]~C[OH0][CH1](~[#6])[OH0][CH0](~[!#1])~[#6]",
                  "[#6]~C([CH0]([NH1][CH1]([CH2][CH0](~[!#1])~O)[CH0](~[!#1])~O)~O)~N","[#6]~C(~[#6])[NH0](~C)~[!#1]","[CH1][NH1][CH0](=[OH0])[CH1]([CH2]c)[NH2]",
                  "[CH3][CH2][CH1]=[CH1][CH2][CH1]=[CH1][CH2][CH1]=[CH1][CH2][CH1]","c[OH0][CH1]([CH1])c","[CH1][OH0]c",
                  "[CH3][CH1]([CH1])[OH1]","[CH3][CH2][CH2][CH2][CH2][CH1]","[CH3][CH2][CH2][CH2][CH2][CH2][CH2][CH2][CH1]","[CH1][CH2]c",
                  "[CH1][CH1]=[CH1][CH1]","[CH0]([CH3])[CH2][CH2][CH0]","[CH1][CH0](=[OH0])[NH1][CH1]([CH2]c)[CH0](=[OH0])[OH1]","[CH1]n",
                  "[CH3][OH0][CH1]","[CH1][CH2][CH1]([CH1])[OH1]","[CH1][CH1]([CH1]([CH2][OH1])[OH1])[OH1]","c[CH3]",
                  "[CH1][CH2][CH2]c","[CH1][CH2][CH2][CH2][CH0](=[OH0])[OH1]","[CH1][CH2][CH0](=[OH0])[NH2]",
                  "[CH1][CH2][CH0][CH3]","[CH0][CH0](=[OH0])[OH1]","[CH0]([CH1])[CH0](=[OH0])[OH0][CH3]","[CH1][CH2][CH2][CH0][OH1]",
                  "[CH0]([CH2][OH1])[OH0][CH1]","[CH1][CH2][CH1]","[CH1][CH2][CH2][CH1]","[CH1]c","[CH1][CH2][OH0][CH1]","[CH1][CH1]([CH2][CH1])[CH1]",
                  "[CH2]([CH2][CH0](=[OH0])[OH1])[CH1]","c[CH0][OH1]","[CH3][OH0][CH0](=[OH0])c","[CH1][CH0](=[OH0])[OH0][CH1]","[CH0][CH2][OH1]",
                  "[NH0][CH0](=[OH0])[CH1]","[CH1][CH2][CH0](=[OH0])[OH1]","[CH0]([CH1])c","[CH1][CH2][CH2][CH2][CH2][NH2]","[CH0][CH0](=[OH0])[OH0][CH1]",
                  "c[CH0]","[CH1][CH2][CH1]=[CH0]","[CH1][CH1]([CH3])[OH0][CH1]","c[CH0](=[OH0])[OH0][CH1]","c[CH0](=[OH0])c","[CH1][CH2][CH2][CH0]",
                  "[CH0][CH1]","[CH1][NH0]","[CH1][CH2][CH0](=[OH0])c","[CH0][CH0]","[CH1][CH2][CH0]","c[CH0](=[OH0])[CH1]","c[CH1]=[CH1][CH0]",
                  "[NH0]c","[CH1][CH2][OH0][PH0](=[OH0])([OH1])[OH1]","c[OH0][CH0]","[CH0][CH0][CH3]","[CH3][CH1]([CH3])[CH1]","c[CH1]=[OH0]",
                  "[CH3][CH2][CH2][CH2][CH2][CH1]=[CH1][CH2][CH1]","c[OH0][CH2][CH0]","[CH0]([CH1])[OH1]","c[CH1]=[CH1][CH0](=[OH0])[OH0][CH1]",
                  "[CH3][CH0](=[OH0])[NH1][CH1]","[CH0](c)c","[CH0][NH2]","[CH0]([CH0])[OH1]","[CH0][NH1]c","[CH0]([CH3])[OH0][CH1]",
                  "[CH3][OH0][CH0][CH0]","c[CH0](=[OH0])[OH1]","[CH1][OH0][CH0](=[OH0])[CH3]")

  col_names_custom <- c("Un8179","Un8181","Un8182","Un8183","Un8184","Un8185","Un8187","Un8190","Un8191","Un8192","Un8193","Un8194",
                        "Un8197","Un8198","Un8199","Un8200","Un8202","Un8203","Un8205","Un8206","Un8207","Un8210","Un8213","Un8215",
                        "Un8217","Un8218","Un8219","Un8220","Un8221","Un8222","Un8224","Un8225","Un8226","Un8228","Un8229","Un8230",
                        "Un8232","Un8233","Un8234","Un8235","Un8236","Un8237","Un8239","Un8240","Un8241","Un8242","Un8243","Un8244",
                        "Un8245","Un8246","Un8247","Un8248","Un8250","Un8251","Un8252","Un8253","Un8254","Un8255","Un8256","Un8258",
                        "Un8259","Un8260","Un8261","Un8262","Un8263","Un8264","Un8266","Un8267","Un8268","Un8269","Un8270","Un8271",
                        "Un8272","Un8273","Un8274","Un8275","Un8276","Un8277","Un8280","Un8282","Un8283","Un8284","Un8285","Un8286",
                        "Un8287","Un8290","Un8292","Un8294","Un8295","Un8296","Un8297","Un8299","Un8303","Un8304","Un8308","Un8309",
                        "Un8310","Un8311","Un8313","Un8314","Un8315","Un8316","Un8318","Un8319","Un8320","Un8321","Un8323","Un8324",
                        "Un8325","Un8327","Un8329","Un8332","Un8333","Un8336","Un8337","Un8338","Un8339","Un8340","Un8341","Un8342",
                        "Un8345","Un8347","Un8348","Un8349","Un8351","Un8353","Un8354","Un8355","Un8356","Un8359","Un8360","Un8363",
                        "Un8367","Un8374","Un8376","Un8378","Un8379","Un8380","Un8381","Un8382","Un8383","Un8384","Un8385","Un8386",
                        "Un8387","Un8388","Un8390","Un8394","Un8395","Un8396","Un8397","Un8399","Un8400","Un8401","Un8402","Un8403",
                        "Un8405","Un8406","Un8407","Un8408","Un8409","Un8410","Un8411","Un8412","Un8413","Un8414","Un8415","Un8416",
                        "Un8417","Un8420","Un8421","Un8422","Un8423","Un8424","Un8426","Un8427","Un8428","Un8431","Un8432","Un8433",
                        "Un8434","Un8436","Un8437","Un8438","Un8439","Un8440","Un8441","Un8443","Un8444","Un8445","Un8446","Un8447",
                        "Un8448","Un8450","Un8451","Un8453","Un8454","Un8455","Un8456","Un8457","Un8460","Un8461")


  custom_table <- as.data.frame(t(list)) %>%
    rownames_to_column()

  for (n in 1:length(compoundslistwithSMILES$SMILES)){

    SMILES <- compoundslistwithSMILES$SMILES[[n]]

    mol2 <- parse.smiles(SMILES)[[1]]

    if(!is.null(mol2)){

      custom_fingerprints <- get.fingerprint(mol2,
                                             type = "substructure",
                                             substructure.pattern = smartslist)

      table <- as.data.frame(custom_fingerprints@bits) %>%
        mutate(fp = 1) %>%
        mutate(fp = as.numeric(fp))
      colnames(table) <- c("row", SMILES)

      table <- table %>%
        mutate(row = as.numeric(row))

      datarow <- list %>%
        left_join(table) %>%
        select(-row)

      newrow <- as.data.frame(t(datarow)) %>%
        rownames_to_column()

      custom_table <- custom_table %>%
        add_row(newrow)
    }

    else{
      print(SMILES)
    }
  }

  custom_table[is.na(custom_table)] <- 0

  custom_table <- custom_table %>%
    dplyr::filter(rowname != "row")

  names(custom_table) <- c("SMILES", col_names_custom)

  row <- c(1:113)
  list <- as.data.frame(row)

  smartslistring <- c("[CH1]([CH1][OH0][CH1]-1~[#7])([CH1]-1~[#8])~[#8]","c(:c(:c(:c:c:1):o:c:c:2):c:2):c:1","C(C(C(C(C-1~[#8])~[#8])~[#8])~[#8])(C-1~[#8])~[#8]",
                      "c(:c:o:c(:c:c(:c:c:1~[#8])~[#8]):c:1:2):c:2~[#8]","c(:c(:c:c:c:1~[#17])~[#8]):c:1~[#8]","c(:c(:n:c(:n:1)~[#8])~[#7])(:c:1~[#8])~[#7]",
                      "c(:c(:c:c:c:1)~[#8]):c:1","C(CC(OC=1)~[#8])C=1","c(:c:c:c:c:1:c:c:n:2):c:1:2","c(:c(:c:c:c:1)~[#16]):c:1",
                      "[CH1]([CH1][CH1]([CH2]-1)~[#8])[CH1]-1~[#8]","c(:c:o:c:c:1):c:1~[#8]","c(:c(:n:c(:n:1)~[#7])~[#7]):c:1~[#8]",
                      "C(CCC(C-1CCC-2)C-2)(C-1CC=C-3C(C(CC-4)CCC-5)C-5)C-3-4",  "c(:c:c(:c:n:n:1):c:1:c:2):c:2","c(:c(:c:c:c:1CCCO-2)~[#8]):c:1-2",
                      "C(CCOC-1~[#8])C-1","c(:c:c:c:c:1):c:1~[#7]","c(:c:c:c(:c:1C(C-2)~[#8])N-2):c:1","c(:c(:c(:c:c:1)~[#8])~[#8])(:c:1)~[#8]",
                      "C(CCC-1)(C-1)~[#8]","c(:c(:c:c(:c:1)~[#8])~[#8])(:c:1~[#8])~[#8]","c(:n:c:c:1):n:1","c(:n:c:n:1):n:1",
                      "C(C(C(OC-1)~[#8])~[#7])(C-1~[#8])~[#8]","[CH1]([CH2][CH2][CH0][CH1]-1)[CH1]-1","c(:c(:c:c:c:1:c:c:c(:o:2)~[#8])~[#8]):c:1:2",
                      "c(:n:c:n:c:1~[#7])(:c:1)~[#7]","n(:c:c:c:c:1):c:1","C(COC-1)C-1","C(CCCC-1CC(CC-2)~[#8])C-1-2","C(C(COC-1)~[#8])(C-1~[#8])~[#8]",
                      "c(:c:o:c(:c:c(:c:c:1~[#8])~[#8]):c:1:2)(:c:2)~[#8]","[CH1]([CH0][CH2][CH2][CH2]-1)[CH0]-1","c(:c(:c:c(:c:1[CH0]([CH1][CH1][OH0]-2)~[#8])~[#8])~[#8]):c:1-2",
                      "C(CCC(C-1CCC-2CCC-3)C-2-3)C-1","c(:c:c:n:c:1~[#8]):n:1","n(:c:c:c(:c:1:c:c:c:c:2)~[#8]):c:1:2",   "[CH1]([CH1][CH1]([CH2]-1)~[#8])[CH0]-1~[#8]",
                      "c(:c(:c:c(:c:1)~[#8])~[#8]):c:1~[#8]","c(:c:c:c(CCCC-1):c-1:2):c:2","C(CCCC-1~[#8])(C-1)~[#8]",
                      "[CH1]([CH2][CH1][OH0][CH1]-1)[CH1]-1","c(:c:c(:c(:c:c:c:o:1):c:1:2)~[#8])(:c:2)~[#8]",
                      "[CH0]([CH1]([CH1]([CH1]=[CH1][OH0]-1)[CH1]-2~[#8])[CH1]-1~[#8])([CH1]-2-3)[OH0]-3",
                      "C(CCC(C-1CC(C-2~[#8])~[#8])C-2)(C-1CC=C-3C(C(CC-4)CCC-5)C-5)C-3-4","[CH0]([CH2][CH2][CH1][CH0]-1)[CH1]-1",
                      "C(CC(C-1)~[#8])C-1~[#8]","c(:c:c:c(:c:1)~[#9])(:c:1)~[#7]","C(COC-1)(C-1)~[#8]","c(:n:c:n:c:1~[#7])(:c:1)~[#8]",
                      "c(:c:c:n:c:1~[#7]):n:1","c(:c(:c(:c(:c:1)~[#8])[CH0]([CH2][CH1]-2)~[#8])[OH0]-2):c:1~[#8]","C(CCCC-1~[#8])C-1",
                      "n(:c:c:c:1):c:1","c(:c:c(:n:c:1)~[#8]):n:1","c(:c:o:c(:c:c:c:c:1):c:1:2)(:c:2~[#8])~[#8]",
                      "[CH1]([CH2][CH2][CH0][CH1]-1)[CH2]-1","[CH0]([CH1][CH2][CH2]-1)[CH0]-1~[#8]","c(:n:c:n:c:1)(:c:1)~[#7]",
                      "C(C(COc:1:c:c(:c:c:2)~[#8])~[#8])c:1:2","c(:c:c(:c(:c:1:c:c:c:2~[#8]):c:2)~[#8]):o:1","c(:c(:c:c:c:1~[#8])~[#8])(:c:1)~[#8]",
                      "C(C=CCC-1CCCC-2)C-1-2","C(CCC(C-1CC(C-2CCC-3)~[#8])C-2-3)C-1","c(:c:c(:c(:c:1:c:c:c:2~[#8]):c:2~[#8])~[#8]):o:1",
                      "C(CCCC-1CCCC-2)C-1-2","[CH0]([CH2][CH2][CH1]([CH0]-1)~[#8])([CH1]-1[CH2][CH2][CH0]-2)[CH1]-2","C(CCC-1CCCC-2)C-1-2",
                      "C(CC(CC-1~[#8])~[#8])(C-1~[#8])~[#8]","C(CC(C-1)~[#8])(O-1)~[#7]","C(CC(CC-1~[#8])~[#8])O-1","c(:c:c:c:c:1OCCC-2~[#8])(:c:1-2)~[#8]",
                      "C(CCC-1C(C(CC-2)C(C(C-3)CC(C-4)~[#8])C-4)C-3)C-1-2",
                      "c(:n:c(:n:c:1)~[#8])(:c:1)~[#7]","C(C(C(C-1)~[#8])~[#8])O-1",
                      "[CH0]([CH1]([CH2][CH2][CH0]([CH1]-1[CH2][CH1]=[CH0]-2[CH1]([CH0]([CH2][CH2]-3)[CH2][CH2][CH0]-4)[CH2]-4)[CH0]-2-3)[CH0]-1[CH2][CH2]-5)[CH1]-5~[#8]",
                      "[CH1]([CH0]-1)[OH0]-1","C(=COCC-1CCC-2)C-1-2","n(:c:n:c:1:c(:n:c:n:2)~[#7]):c:1:2","C(CCc(:c:c:c:c:1):c:1-2)O-2",
                      "c(:c:o:c(:c:1:c:c:c:2~[#8]):c:2)(:c:1~[#8])~[#8]","[CH0](=[CH1][OH0][CH1][CH1]-1)[CH1]-1","c(:c:c(:n:c:1~[#7])~[#7]):n:1",
                      "c(:c:c(:c:c:1)~[#8]):n:1","C(C=CC-1~[#8])O-1","c(:n:n:c:1):c:1","c(:c(:c(:c(:c(:c:c(:c:1)~[#8])~[#8]):c:1:2)~[#8])~[#8]):o:2",
                      "[CH1]([CH1][CH1][CH1]-1)[OH0]-1","[CH0]([CH1][CH0][CH2][CH2]-1)[CH0]-1","c(:c:c:c:c:1:n:c:n:2):c:1:2","n(:c:n:c:1:c(:n:c:n:2)~[#8]):c:1:2",
                      "c(:n:c:n:c:1:n:c:n:2):c:1:2","c(:n:c:c:c:1~[#8])(:n:1)~[#7]",  "[CH0]([CH2][CH2][CH0][CH1]-1)[CH2]-1","[CH0](=[CH1][CH0]([OH0]-1)~[#8])[CH2]-1","c(:c:n:c(:n:1)~[#8]):c:1~[#8]",
                      "C(COc(:c:c:c:c:1):c:1-2)(C-2~[#8])~[#8]","c(:c(:c:c(:c:c:c(:c:1)~[#8]):c:1:2)~[#8]):o:2","[CH2]([CH2][CH1][CH0][CH1]-1)[CH1]-1",
                      "c(:c(:c:c(:c:1:o:c:c:2):c:2~[#8])~[#8])(:c:1)~[#8]","C(COc(:c:c:c:c:1):c:1-2)C-2~[#8]","C(CC=CO-1)C-1",
                      "[CH0]([CH2][CH2][CH1]([CH0]-1)~[#8])[CH1]-1","C(COC-1~[#8])C-1","C(COC-1~[#8])(C-1~[#8])~[#8]","n(:c:c:n:1):c:1~[#16]",
                      "[CH1]([CH1][CH1][CH1][CH2]-1)[OH0]-1","c(:n:c:n:c:1:c:c:c:c:2)(:c:1:2)~[#8]","[CH1]([CH1][CH1]([CH2][CH1]-1~[#8])~[#8])[OH0]-1",
                      "c1ccccc1-c","[cD2]1ccc[cD2]c1-c(c)c",
                      "[C](-,=[C]-,=1)-,=[C](-,=[C]-,=2)-[C](-,=[C]-,=[C]-,=[C]-,=2-,=[O])-,=[C](-,=[C]-,=3)-,=[C]-,=1-,=[C](-,=[C]-,=4)-,=[C](-,=[C]-,=3)-,=[C]-,=[C]-,=4")


  col_names_ring <- c("Un8468","Un8480","Un8483","Un8484","Un8486","Un8489","Un8491","Un8494","Un8495","Un8496","Un8497","Un8498",
                      "Un8501","Un8502","Un8505","Un8516","Un8521","Un8524","Un8529","Un8539","Un8540","Un8541","Un8542","Un8546",
                      "Un8547","Un8553","Un8556","Un8558","Un8559","Un8563","Un8566","Un8568","Un8575","Un8581","Un8582","Un8584",
                      "Un8585","Un8587","Un8588","Un8589","Un8606","Un8607","Un8610","Un8611","Un8618","Un8619","Un8620","Un8632",
                      "Un8635","Un8637","Un8638","Un8641","Un8647","Un8648","Un8651","Un8655","Un8663","Un8664","Un8670","Un8673",
                      "Un8677","Un8678","Un8682","Un8686","Un8691","Un8692","Un8698","Un8706","Un8710","Un8723","Un8724","Un8730",
                      "Un8736","Un8740","Un8743","Un8749","Un8750","Un8751","Un8758","Un8760","Un8767","Un8770","Un8772","Un8777",
                      "Un8779","Un8782","Un8783","Un8788","Un8789","Un8797","Un8802","Un8808","Un8814","Un8820","Un8834","Un8843",
                      "Un8849","Un8852","Un8854","Un8858","Un8860","Un8864","Un8871","Un8874","Un8881","Un8888","Un8889","Un8901",
                      "Un8903","Un8908","Un8922","Un8923","Un8924")

  ring_table <- as.data.frame(t(list)) %>%
    rownames_to_column()

  for (n in 1:length(compoundslistwithSMILES$SMILES)){
    SMILES <- compoundslistwithSMILES$SMILES[[n]]

    mol2 <- parse.smiles(SMILES)[[1]]

    if(!is.null(mol2)){

      ring_fingerprints <- get.fingerprint(mol2,
                                           type = "substructure",
                                           substructure.pattern = smartslistring)

      table <- as.data.frame(ring_fingerprints@bits) %>%
        mutate(fp = 1) %>%
        mutate(fp = as.numeric(fp))
      colnames(table) <- c("row", SMILES)

      table <- table %>%
        mutate(row = as.numeric(row))

      datarow <- list %>%
        left_join(table) %>%
        select(-row)

      newrow <- as.data.frame(t(datarow)) %>%
        rownames_to_column()

      ring_table <- ring_table %>%
        add_row(newrow)
    }

    else{
      print(SMILES)
    }
  }

  ring_table[is.na(ring_table)] <- 0

  ring_table <- ring_table %>%
    dplyr::filter(rowname != "row")

  names(ring_table) <- c("SMILES", col_names_ring)


  #Combining all FP together
  final_table <- substr_table %>%
    left_join(maccs_table) %>%
    left_join(pubchem_table) %>%
    left_join(KlekotaRoth_table) %>%
    left_join(custom_table) %>%
    left_join(ring_table)

  return(final_table)
}




# opening zip files for SIRIUS 5

#' @export
UnZip_SIRIUS5 <- function(folderwithSIRIUSfiles) {
  subfolder_fp_zip <- dir(folderwithSIRIUSfiles, all.files = TRUE, recursive = TRUE, pattern = "fingerprints")

  for (zipF in subfolder_fp_zip){
    outfolder <- str_split(zipF, "/")[[1]][1]
    outDir <- paste(folderwithSIRIUSfiles, outfolder, "fingerprints1",
                    sep = "/")
    zipFile <- paste(folderwithSIRIUSfiles, zipF, sep = "/")
    unzip(zipFile, exdir = outDir)
  }


  subfolder_scores_zip <- dir(folderwithSIRIUSfiles, all.files = TRUE, recursive = TRUE, pattern = "scores")

  for (zipF in subfolder_scores_zip) {
    outfolder <- str_split(zipF, "/")[[1]][1]
    outDir <- paste(folderwithSIRIUSfiles, outfolder, "scores1",
                    sep = "/")
    zipFile <- paste(folderwithSIRIUSfiles, zipF, sep = "/")
    unzip(zipFile, exdir = outDir)
  }
}

# purpose: activity annotation of nodes by a molecular netowrk
# date: 2024-11-19
# author: Yvonne Kreutzer <yvonne.kreutzer@mmk.su.se>
#
#_______________________________________________________________________________#
#

#_______________________________________________________________________________#

# libraries and functions
library(tidyverse)
library(janitor)



# 1) loading data---------------------------------------------------------------

# 1.1) best performing network settings______________________________________________

setwd(paste( "C:/Users/", admin, "OneDrive - Kruvelab/Yvonne/01_Projects/Paper_1_MS2Network/03_output/1_analyzed/Molecular network finals", sep = "/"))
df_best_performance <- read.delim("20240924_MolNetwork_test_performance.tsv") %>% 
  mutate(assay = ifelse(assay == "nr_ah_r", "nr_ahr", assay))

# 1.2) similarities between spectra__________________________________________________
# load similarities, calculated between measured WWTP spectra and spectra from the MS2 dataset
# MS2DeepScore similarities were separately calculated according to the tuturial 
# available at https://github.com/matchms/ms2deepscore/tree/main

df_spectra_similarites <- read.delim(r"(C:\Users\yvkr1259\OneDrive - Kruvelab\Yvonne\01_Projects\Paper_1_MS2Network\03_output\2_publish\20241127_example_similarities.tsv)",
                                     sep = "\t") 
# only keep the connections between measured spectra to spectra in the MS2 dataset
df_spectra_similarites <- df_spectra_similarites %>%
  mutate(similarity = ifelse(similarity < 0, 0, similarity)) %>%
  mutate(connection_type = case_when(
    (grepl("^[0-9]+$", ID_A) & grepl("^[0-9]+$", ID_B)) ~ "only_measured_spectra",
    (grepl("^[0-9]+$", ID_A) & !(grepl("^[0-9]+$", ID_B))) ~ "measured_to_tox",
    (grepl("^[0-9]+$", ID_B) & !(grepl("^[0-9]+$", ID_A))) ~ "measured_to_tox",
    TRUE ~ "only_tox"
  )) %>%
  filter(connection_type == "measured_to_tox")

# 1.3) tox data______________________________________________________________________

df_tox21 <- read.delim(r"(C:\Users\yvkr1259\OneDrive - Kruvelab\Yvonne\01_Projects\Paper_1_MS2Network\03_output\2_publish\other SI\InChIKey14_split_tox21_data.tsv)") %>% 
  clean_names() %>% 
  rename(inchikey14 = in_ch_i_key14)

unique_assays <- colnames(df_tox21[, 3:8])

# 2) annotation of measured WWTP spectra----------------------------------------

# get relevant spectra IDs
spectra_IDs = df_spectra_similarites %>%
  select(ID_A, ID_B) %>%
  unlist() %>%
  .[grepl("^[0-9]+$", .)] %>%
  unique()

# 2.1) annotation loop_______________________________________________________________

# getting propability lablels 
df_spectra_annotations <- tibble()

for (current_assay in unique_assays) {
  
  df_intermediate_results <- tibble()
  cat(sprintf("%s \t %s\n", Sys.time(), current_assay))
  
  # get best network parameters for selected assay
  current_assay_best_performance = df_best_performance %>% 
    distinct(assay, similarity_threshold) %>% 
    filter(assay == current_assay)
  similarity_thresshold = unlist(current_assay_best_performance$similarity_threshold)
  
  # apply similarity threshold
  current_spectra_similarites = df_spectra_similarites %>% 
    filter(similarity >= similarity_thresshold)
  
  
  for (target_compound in spectra_IDs) {
    # get all compounds that are directly connected to the target compound
    target_compound_connections = current_spectra_similarites %>% 
      filter(ID_A == target_compound|ID_B == target_compound)
    
    connected_compounds <- target_compound_connections %>%
      select(ID_A, ID_B) %>%
      pivot_longer(everything(), values_to = "compound") %>%
      pull(compound) %>%
      unique()
    connected_compounds <- setdiff(connected_compounds, target_compound)
    
    propability_label <- df_tox21 %>%
      filter(inchikey14 %in% connected_compounds) %>%
      select(!!sym(current_assay)) %>%
      drop_na() %>%
      summarise(propability_label = mean(!!sym(current_assay), na.rm = TRUE)) %>%
      pull(propability_label)
    
    
    current_results = data.table(
      inchikey14 = target_compound,
      assay = current_assay,
      count_connected_compound = length(connected_compounds),
      prob_label = propability_label)
    
    df_intermediate_results = rbind(df_intermediate_results, current_results)  
  }
  
  df_spectra_annotations = rbind(df_spectra_annotations, df_intermediate_results)
}

df_spectra_annotations = df_spectra_annotations %>% 
  rename(ID = inchikey14)

# 2.2) apply threshold to probabilities and get results__________________________

# add prediction thresholds and make predictions
df_spectra_annotations = df_spectra_annotations %>% 
  left_join(df_best_performance %>%
              select(assay, metrics_set, prediction_threshold),
            relationship = "many-to-many")

df_spectra_annotations = df_spectra_annotations %>% 
  mutate(predicted_label = ifelse(prob_label>= prediction_threshold, 1, 0))

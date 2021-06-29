# create metadata tables for TCGA lung cancer samples


# load R packages
library(tidyverse)


# load clinical information from the TCGA pan-cancer atlas publication
# https://gdc.cancer.gov/about-data/publications/pancanatlas
# https://doi.org/10.1016/j.cell.2018.03.022
tcga_pancan_metadata <- data.table::fread(file = "./data/raw/PancanAtlas_publication/clinical_PANCAN_patient_with_followup.tsv") %>%
  as_tibble()

tcga_luad_metadata <- tcga_pancan_metadata %>%
  filter(acronym == "LUAD") %>%
  select(sample = bcr_patient_barcode, tcga_cohort = acronym, gender, year_of_diagnosis = year_of_initial_pathologic_diagnosis, age_at_diagnosis = age_at_initial_pathologic_diagnosis, race, ethnicity, tissue_source_site, pathologic_stage, tobacco_smoking_history, stopped_smoking_year, year_of_tobacco_smoking_onset, histological_type, number_pack_years_smoked) %>%
  mutate(across(.cols = everything(), .fns = ~ replace(.x, .x %in% c("[Discrepancy]", "[Unknown]", "[Not Available]", "[Not Evaluated]"), NA))) %>%
  mutate(across(.cols = c("year_of_diagnosis", "age_at_diagnosis", "stopped_smoking_year", "year_of_tobacco_smoking_onset", "number_pack_years_smoked"), .fns = as.numeric)) %>%
  #mutate_at(.vars = c("year_of_diagnosis", "age_at_diagnosis", "stopped_smoking_year", "year_of_tobacco_smoking_onset", "number_pack_years_smoked"), .funs = as.numeric) %>%
  #mutate_if(.predicate = colnames(.) %in% c("year_of_diagnosis", "age_at_diagnosis", "stopped_smoking_year", "year_of_tobacco_smoking_onset", "number_pack_years_smoked"), .funs = as.numeric) %>%
  mutate(sample = str_replace_all(sample, "-", ".")) %>%
  mutate(date_of_birth = year_of_diagnosis-age_at_diagnosis) %>%
  mutate(smoker = if_else(tobacco_smoking_history == "Lifelong Non-smoker", "NO", "YES")) %>%
  select(sample, tcga_cohort, gender, date_of_birth, everything())

tcga_lusc_metadata <- tcga_pancan_metadata %>%
  filter(acronym == "LUSC") %>%
  select(sample = bcr_patient_barcode, tcga_cohort = acronym, gender, year_of_diagnosis = year_of_initial_pathologic_diagnosis, age_at_diagnosis = age_at_initial_pathologic_diagnosis, race, ethnicity, tissue_source_site, pathologic_stage, tobacco_smoking_history, stopped_smoking_year, year_of_tobacco_smoking_onset, histological_type, number_pack_years_smoked) %>%
  mutate(across(.cols = everything(), .fns = ~ replace(.x, .x %in% c("[Discrepancy]", "[Unknown]", "[Not Available]", "[Not Evaluated]"), NA))) %>%
  mutate(across(.cols = c("year_of_diagnosis", "age_at_diagnosis", "stopped_smoking_year", "year_of_tobacco_smoking_onset", "number_pack_years_smoked"), .fns = as.numeric)) %>%
  #mutate_at(.vars = c("year_of_diagnosis", "age_at_diagnosis", "stopped_smoking_year", "year_of_tobacco_smoking_onset", "number_pack_years_smoked"), .funs = as.numeric) %>%
  #mutate_if(.predicate = colnames(.) %in% c("year_of_diagnosis", "age_at_diagnosis", "stopped_smoking_year", "year_of_tobacco_smoking_onset", "number_pack_years_smoked"), .funs = as.numeric) %>%
  mutate(sample = str_replace_all(sample, "-", ".")) %>%
  mutate(date_of_birth = year_of_diagnosis-age_at_diagnosis) %>%
  mutate(smoker = if_else(tobacco_smoking_history == "Lifelong Non-smoker", "NO", "YES")) %>%
  select(sample, tcga_cohort, gender, date_of_birth, everything())


write_tsv(x = tcga_luad_metadata, file = "./data/processed/metadata/tcga_luad_metadata_pancancer_atlas.txt")
write_tsv(x = tcga_lusc_metadata, file = "./data/processed/metadata/tcga_lusc_metadata_pancancer_atlas.txt")

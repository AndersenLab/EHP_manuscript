#!/usr/bin/env Rscript
library(tidyverse)

# set working dir to base directory of repository
setwd(paste0(dirname(rstudioapi::getActiveDocumentContext()$path),"/.."))

#===================================================================#
# Step 1: Read data files, clean, join, filter
#===================================================================#
# read in raw data from SG
raw_dat3 <- readxl::read_excel("data/raw/Data_Andersen_All.xlsx", na = c("NA", ""))

# clean it
proc_dat3 <- raw_dat3 %>%
  dplyr::mutate(filter = case_when(latin_name == "Daphnia ambigua" & is.na(duration_d) ~ "remove",
                                   TRUE ~ "keep")) %>%
  dplyr::filter(filter == "keep") %>%
  dplyr::select(-filter) %>%
  dplyr::mutate(effect_unit = ifelse(effect_unit == "mg/l", "mg/L", effect_unit)) # fix mg/l

# get the cas numbers from the Andersen set - why are there only 18 cas numbers here? 
cas_keep <- proc_dat3 %>%
  dplyr::rename(chem_name2 = chem_name) %>%
  dplyr::distinct(cas, chem_name2) 

cas_keepv <- cas_keep %>%
  dplyr::pull(cas)

chem_keepv <- cas_keep %>%
  dplyr::pull(chem_name2)

# read in windy-boyd data set
raw_dat4 <- readxl::read_excel("data/raw/Data_Boyd_EnviroTox.xlsx", na = c("NA", "")) 

# filter the data to cas in andersen experiments
proc_dat4 <- raw_dat4 %>%
  dplyr::filter(cas %in% cas_keepv) %>%
  dplyr::left_join(cas_keep) %>%
  dplyr::mutate(chem_name = chem_name2) %>%
  dplyr::select(-chem_name2) %>%
  dplyr::mutate(source = ifelse(source == "EnviroToxDB", "EnviroTox DB", source))

# join the processed 3 and 4 data
join_dat <- dplyr::bind_rows(proc_dat3, proc_dat4) %>%
  dplyr::mutate(id = paste0(cas, chem_name, group, latin_name, strain, test_type, test_statistic, duration_d, effect_value,
                            effect_unit, endpoint, source)) %>%
  dplyr::distinct(id, .keep_all = T) %>%
  dplyr::select(-id)

# look for duplications since the row numbers above don't make sense. 262 rows duplicated in proc_dat3, some up to 4 times.
# Let's just trust it for now.
dup.test <- proc_dat3 %>%
  dplyr::mutate(id = paste0(cas, chem_name, group, latin_name, strain, test_type, test_statistic, duration_d, effect_value,
                            effect_unit, endpoint, source)) %>%
  dplyr::group_by(id) %>%
  dplyr::mutate(n = n()) %>%
  dplyr::filter(n > 1) %>%
  dplyr::distinct(id, .keep_all = T)

# save the joined data
readr::write_csv(join_dat, file = "data/processed/03_clean.csv")

# # OLD DATA READ IN
# #===================================================================#
# # Step X: Read second data file and clean if needed
# #===================================================================#
# 
# # load effect data taken from 'master' tab in Master_File.xlxs
# raw_dat <- readxl::read_excel("data/raw/Master_File.xlsx", sheet = 'master', na = c("NA", ""))
# 
# # clean it
# proc_dat <- raw_dat %>%
#   dplyr::select(-Standard.Error) %>% # not needed b/c all NA
#   janitor::clean_names(case = "snake") %>% # fix var names to all be snake case
#   dplyr::select(cas, chem_name, group, latin_name, strain, test_type, test_statistic, duration_h, effect_value) # reorder for clarity
# 
# # save it
# readr::write_csv(proc_dat, file = "data/processed/01_clean.csv")
# 
# #===================================================================#
# # Step X: Read second data file and clean if needed
# #===================================================================#
# # read in raw data
# raw_dat2 <- readxl::read_excel("data/raw/02_clean_SG.xlsx", na = c("NA", ""))
# 
# # the only issue is that Daphnia ambigua has a single duration as "instar"
# # save it
# readr::write_csv(raw_dat2, file = "data/processed/02_clean.csv")
# 

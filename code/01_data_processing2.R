#!/usr/bin/env Rscript
library(tidyverse)

# set working dir to base directory of repository
setwd(paste0(dirname(rstudioapi::getActiveDocumentContext()$path),"/.."))

# output the CAS#s tested in Widmayer 2022, these are the wrong CAS numbers!
# https://github.com/AndersenLab/toxin_dose_responses/blob/master/manuscript_tables/supp.table.1.csv file is not correct!
# DO NOT USE
# 1910425 = paraquat ce - RIGHT
# 75365730 = paraquat all_cas - WRONG
# 7646857 = Zinc chloride ce - RIGHT
# 646857 = Zinc dichloride all_cas
# all_cas <- data.table::fread("data/raw/Widmayer_2022_cas.csv") %>%
#  dplyr::mutate(cas = stringr::str_replace_all(CAS, pattern = "-", replacement = "")) %>% # fix cas to match ce
#  dplyr::mutate(cas = case_when(Toxicant == "Methylmercury dichloride" ~ 115093,
#                                Toxicant == "Aldicarb" ~ 116063,
#                                Toxicant == "Paraquat" ~ 1910425,
#                                Toxicant == "Zinc dichloride" ~ 7646857)) 
#  dplyr::select(chem_name = Toxicant, cas)
# rio::export(all_cas, file = "data/processed/archive/toxicant_cas.csv")

# pull in CAS numbers and molecular weights from the data/processed/toxcast_mw.xlsx file
# that was exported from toxcast. Note, this file has a missing MW for Mancozeb.

#==============================================================================#
# Testing: Get original cleaned data and match coded data prcessing
#==============================================================================#
orig <- data.table::fread("data/processed/03_clean.csv") 
names(orig)  
unique(orig$source)

#===================================================================#
# Step 1: Read in Widmayer data and make list of 18 chemicals
#===================================================================#
# read in raw C. elegans data from SG: 18 chemicals from 22 total raw
ce <- readxl::read_excel("data/raw/Data_Andersen_All.xlsx", na = c("NA", "")) %>%
  dplyr::filter(source == "Widmayer 2022")  %>%
  dplyr::mutate(effect_unit = ifelse(effect_unit == "mg/l", "mg/L", effect_unit)) %>% # fix mg/l
  dplyr::mutate(duration_d = ifelse(source == "Widmayer 2022", 2, duration_d)) %>%
  dplyr::relocate(endpoint, .after = source) # flip to match output
length(unique(ce$chem_name))

# get chemical names, cas, and molecular weight data from toxcast_mw and edit mancozeb
# Also, paraquat CAS is for paraquat dichloride with MW of ~257, which might not match paraquat used in other studies.
cas_df <- readxl::read_excel("data/raw/toxcast_mw.xlsx") %>%
  dplyr::mutate(INPUT = as.character(INPUT)) %>%
  dplyr::distinct(.keep_all = T) %>%
  dplyr::mutate(AVERAGE_MASS = ifelse(CASRN == "8018-01-7", 541.08, AVERAGE_MASS)) %>% # Add Mancozeb
  dplyr::select(CASRN, AVERAGE_MASS, PREFERRED_NAME) %>%
  dplyr::mutate(cas = stringr::str_replace_all(CASRN, pattern = "-", replacement = ""),
                cas = as.numeric(cas)) %>%
  dplyr::left_join(dplyr::select(ce, cas, chem_name)) %>%
  dplyr::distinct(chem_name, .keep_all = T) %>%
  dplyr::filter(!is.na(chem_name)) %>% # remove paraquat and keep only paraquat dichloride cas
  dplyr::select(cas, chem_name, AVERAGE_MASS)

# test names match perfectly to orig
identical(names(orig), names(ce))

#===================================================================#
# Step 2: Pull the Boyd data from the EHP file and keep all good chems
#===================================================================#
# read in original boyd data from Scott: 763 raw observations, 
boyd_orig <- readxl::read_excel("data/raw/Data_Boyd_EnviroTox.xlsx", na = c("NA", "")) %>%
  dplyr::filter(source == "Boyd et al")
  length(unique(boyd_orig$cas))

# read in boyd data, convert AC50 to mg/L, filter out extrapolated AC50s.
boyd <- data.table::fread("data/raw/Boyd/ehp.1409645.s002.csv") %>%
  dplyr::slice(-1) %>% # drop first row
  janitor::row_to_names(1) %>% # set names
  dplyr::mutate(CAS = stringr::str_replace_all(CASRN, pattern = "-", replacement = ""),
                CAS = as.numeric(CAS)) %>% # fix CAS for joining convert NOCAS to NA
  dplyr::left_join(cas_df, by = c("CAS" = "cas")) %>%
  dplyr::mutate(AC50_uM = as.numeric(`C.elegans AC50 (Hill function)`),
                AC50 = (AC50_uM / 1e6) * AVERAGE_MASS * 1000) %>%
  dplyr::mutate(QC = ifelse(AC50_uM < 200.5, "PASS", "FAIL")) #AC50s at or below the maximum tested concentration of 200uM

# get list of boyd chems 467 pass filter, 459 not in Widmayer (8 overlap)
boyd_chem <- boyd %>%
  dplyr::distinct(CAS, .keep_all = T) %>% # remvoe dups
  dplyr::filter(!is.na(CAS)) %>% # remove NOCAS
  dplyr::select(cas = CAS, chem_name = Name) %>%
  dplyr::filter(!(cas %in% cas_df$cas)) # remove chems that overlap with widmayer

# ouput the cas numbers
rio::export(boyd_chem, "data/processed/boyd_chem.csv")

# add to cas list
all_ce_cas <- readxl::read_excel("data/raw/toxcast_boyd_mw.xlsx") %>%
  dplyr::select(CASRN = INPUT, AVERAGE_MASS, PREFERRED_NAME) %>%
  dplyr::mutate(cas = stringr::str_replace_all(CASRN, pattern = "-", replacement = ""),
                cas = as.numeric(cas)) %>%
  dplyr::filter(!is.na(AVERAGE_MASS)) %>% # remove NO MW
  dplyr::mutate(AVERAGE_MASS = ifelse(cas == 27176870, 326.49, AVERAGE_MASS)) %>%
  dplyr::select(cas, chem_name = PREFERRED_NAME, AVERAGE_MASS) %>%
  dplyr::filter(!(cas %in% cas_df$cas)) %>%
  dplyr::bind_rows(cas_df)

# process boyd data to join by adding other variables
proc_boyd <- boyd %>%
  dplyr::left_join(all_ce_cas)
  dplyr::mutate(group = "NEMATODE_COPAS1",
                latin_name = "Caenorhabditis elegans",
                strain = NA_character_,
                test_type = NA_character_,
                test_statistic = "AC50",
                duration_d = NA_real_, # do we have a duration of exposure? yes right?
                effect_unit = "mg/L",
                source = "Boyd et al",
                endpoint = "Growth") %>%
  select(cas = CAS, chem_name, group:duration_d, effect_value = AC50, effect_unit:endpoint, QC)

# test names match perfectly to orig
identical(names(orig), names(proc_boyd))

#===================================================================#
# Step 2: Read Karmaus 2022 folder (RAT_EMPIRICAL)
#===================================================================#
# read in from Karmaus folder - 1885 chems filtered to 14 after matching ce cas
# Pull the data from toxsci-21-0357-File010.xlsx only.
kar <- readxl::read_excel("data/raw/Karmaus/toxsci-21-0357-File010.xlsx", na = c("NA", "")) %>%
  dplyr::mutate(cas = stringr::str_replace_all(CASRN, pattern = "-", replacement = ""),
                cas = as.numeric(cas)) %>% # fix cas to match ce
  dplyr::filter(cas %in% all_ce_cas$cas) %>%
  dplyr::left_join(all_ce_cas)
length(unique(kar$CASRN))

# process the karmaus data for joining
names(orig)
glimpse(orig)
glimpse(kar)
proc_kar <- kar %>%
  dplyr::mutate(cas = as.integer(cas),
                group = "RAT_EMPIRICAL",
                latin_name = "Rat",
                strain = NA_character_,
                test_type = NA_character_,
                test_statistic = "LD50",
                duration_d = NA_integer_,
                effect_value = LD50_mgkg,
                effect_unit = "mg/kg",
                source = "Karmaus 2022",
                endpoint = "Mortality") %>%
  dplyr::select(cas, chem_name, group:endpoint) 
  
# test names match perfectly to orig
identical(names(orig), names(proc_kar))
#===================================================================#
# Step 3: Read and process Comptox folder - RAT test model predictions
#===================================================================#
# read in the output from comptox: 30891 raw observations, 14 pass TEST filter
comptox <- readxl::read_excel("data/raw/Comptox/20240717_comptox_download.xlsx", na = c("NA", ""),
                              sheet = "Toxval Details") %>%
  dplyr::filter(SOURCE == "TEST") # filter to acute oral LD50 TEST data

# process the data for joining, 9 chems pass ce filter
names(orig)
glimpse(orig)
glimpse(comptox)

proc_comptox <- comptox %>%
  dplyr::mutate(cas = stringr::str_replace_all(CASRN, pattern = "-", replacement = ""),
                cas = as.numeric(cas)) %>% # fix cas to match ce
  dplyr::filter(cas %in% all_ce_cas$cas) %>% # filter to ce chems
  dplyr::left_join(all_ce_cas) %>% # join chem names
  dplyr::mutate(cas = as.integer(cas),
                group = "RAT_EMPIRICAL",
                latin_name = "Rat",
                strain = NA_character_,
                test_type = NA_character_,
                test_statistic = "LD50",
                duration_d = NA_integer_,
                effect_value = TOXVAL_NUMERIC,
                effect_unit = "mg/kg",
                source = "EPA CompTox",
                endpoint = "Mortality",
                QC = NA_character_) %>%
  dplyr::select(cas, chem_name, group:endpoint, QC) 

# test names match perfectly to orig
identical(names(orig), names(proc_comptox))

#===================================================================#
# Step 5: Load zebrafish data from toxcast Rdata file
#===================================================================#
# load the mc5 toxcast data
load("data/processed/mc5-6_winning_model_fits-flags_invitrodb_v4_1_SEPT2023.Rdata")

# Get assay description data so we can use biological_process_target (BPT) to filter to relevant assays
assays <- readxl::read_excel("data/raw/assay_annotations_invitrodb_v4_1_SEPT2023.xlsx", sheet = "annotations_combined")

# view all the assay biological targets
sort(unique(assays$biological_process_target))

# filter assays to desired zebrafish assays then join mc5 data
zf1 <- assays %>%
  dplyr::filter(organism == "zebrafish" & assay_design_type_sub %in% c("embryo development","embryonic mortality")) %>% # The Padilla et al. 2012 data are here.
  dplyr::left_join(., mc5) %>% # join in the toxcast data
  dplyr::mutate(cas = stringr::str_replace_all(casn, pattern = "-", replacement = ""),
                cas = as.numeric(cas)) %>% #remove hyphens
  dplyr::filter(cas %in% all_ce_cas$cas) %>% # filter to worm cas # cas_numbers
  dplyr::left_join(all_ce_cas) %>% # join MW data
  mutate(AC50_mg_L = (ac50 / 1e6) * AVERAGE_MASS * 1000, # Convert µM to mg/L
         AC20_mg_L = (ac20 / 1e6) * AVERAGE_MASS * 1000,
         AC10_mg_L = (ac10 / 1e6) * AVERAGE_MASS * 1000,
         AC5_mg_L = (ac5 / 1e6) * AVERAGE_MASS * 1000,
         BMD_mg_L = (bmd / 1e6) * AVERAGE_MASS * 1000,
         ACC_mg_L = (acc / 1e6) * AVERAGE_MASS * 1000) %>%
  dplyr::mutate(flag.length = ifelse(is.na(flag.length), 0, flag.length)) %>%
  dplyr::mutate(QC = ifelse(fitc %in% c(37, 38, 41, 42) & (flag.length < 5 | is.na(flag.length)), "PASS", "FAIL")) %>% # QC filter 
  dplyr::group_by(aeid, casn) %>%
  dplyr::arrange(desc(fitc)) %>% # arrange the data so best fits are at the top of each assay, chem
  dplyr::mutate(aeid_chem_pass = ifelse(sum(QC == "PASS") >= 1, "PASS", "FAIL"),
                aeid_chem_pass_n = sum(QC == "PASS"),
                aeid_chem_fail_n = sum(QC == "FAIL")) %>%
  dplyr::group_by(aeid) %>%
  dplyr::mutate(total_n_compounds = length(unique(chnm)),
                pass_temp = paste0(casn,"_",aeid_chem_pass),
                total_n_compounds_pass = sum(stringr::str_count(unique(pass_temp), pattern = "_PASS"))) %>%
  dplyr::select(aeid, casn, chnm, spid, fitc, QC, aeid_chem_pass, aeid_chem_pass_n, aeid_chem_fail_n, total_n_compounds, pass_temp, total_n_compounds_pass, AC50_mg_L:ACC_mg_L, biological_process_target, everything()) %>%
  dplyr::ungroup() %>%
  dplyr::arrange(aeid, chnm, desc(QC), desc(fitc)) %>%
  dplyr::mutate(n_assays = length(unique(aeid)))

# process the data for joining, 9 chems pass ce filter
names(orig)
glimpse(orig)
glimpse(zf1)

proc_zf1 <- zf1 %>%
  dplyr::mutate(cas = as.integer(stringr::str_replace_all(casn, pattern = "-", replacement = "")),
                group = ifelse(aeid == 1372, "TC_ZF_Tanguay", "TC_ZF_Padilla"), # updated group name to reflect toxcast_zebrafish_source
                latin_name = "Danio rerio",
                strain = NA_character_,
                effect_unit = "mg/L",
                source = "TOXCAST",
                chem_name = case_when(chnm == "2,4-Dichlorophenoxyacetic acid" ~ "2,4-D",
                                      chnm == "Methylmercuric(II) chloride" ~ "Methylmercury chloride",
                                      TRUE ~ chnm),
                test_type = NA,
                duration_d = timepoint_hr / 24,
                endpoint = ifelse(assay_component_endpoint_name == "Tanguay_ZF_120hpf_MORT", "Mortality", "TERATOSCORE")) %>% # add variables to join
  dplyr::select(cas = casn, chem_name, group, latin_name, strain, test_type,
                duration_d, effect_unit, source, endpoint,
                AC50_mg_L, AC20_mg_L, AC10_mg_L, AC10_mg_L, AC5_mg_L, BMD_mg_L, ACC_mg_L,
                hitc, fitc, QC, flag.length) %>% # select what's needed
  tidyr::pivot_longer(cols = AC50_mg_L:ACC_mg_L,
                      names_to = "test_statistic", values_to = "effect_value") %>% # reshape long
    dplyr::mutate(test_statistic = stringr::str_replace(test_statistic, pattern = "_mg_L", replacement = "")) %>% # make clean test_stat
  dplyr::select(cas, chem_name, group, latin_name, strain, test_type,
                test_statistic, duration_d, effect_value, effect_unit, source, endpoint, QC) %>% # select proper order
  #dplyr::filter(!is.na(effect_value), QC == "PASS") %>% # filter by QC
  dplyr::filter(!is.na(effect_value)) %>%
  dplyr::select(cas, chem_name, group:endpoint, QC) # set proper order of variables

# test names match perfectly to orig
identical(names(orig), names(proc_zf1))

#===================================================================#
# Step 4: Read Zebrafish folder to get Su and Scholz data
#===================================================================#
# read in su data: 427 raw observations, 164 chem pass ce filter
# Leave duration as NA, and include all
zf_su <- readxl::read_excel("data/raw/Zebrafish/Su et al 2021/1-s2.0-S0048969721027765-mmc2.xls",
                            na = c("NA", ""),
                            sheet = "FET-LC50") %>%
  dplyr::mutate(cas = as.numeric(stringr::str_replace_all(`CAS number`, pattern = "-", replacement = ""))) %>% # fix cas to match ce
  dplyr::group_by(`Duration (h)`) %>%
  dplyr::mutate(n_chem_timepoint = length(unique(`Chemical name`))) %>%
  dplyr::ungroup() 

# process the data for joining, 9 chems pass ce filter
names(orig)
glimpse(orig)
glimpse(all_ce_cas)
glimpse(zf_su)

proc_zf_su <- zf_su %>%
  dplyr::filter(cas %in% all_ce_cas$cas, `Duration (h)` %in% c("120", "96", "48", "72")) %>% # filter to ce chems and useful time points
  dplyr::left_join(all_ce_cas) %>% # join chem names
  dplyr::mutate(group = "ZF_EMBRYO_SU",
                latin_name = "Danio rerio",
                strain = NA_character_,
                test_type = `Test type`,
                test_statistic = "LC50",
                duration_d = as.numeric(`Duration (h)`) / 24,
                effect_value = as.numeric(`LC50 (µg/L)`) / 1000, #convert to mg/L
                effect_unit = "mg/L",
                source = "Su et al 2021",
                endpoint = "Mortality",
                QC = NA_character_) %>%
  dplyr::select(cas, chem_name, group:endpoint, QC)

# test names match perfectly to orig
identical(names(orig), names(proc_zf_su))

#-------------------------------------------------------------------------------
# TIM WORKGING HERE to read in and process below - could use all timepoints and leave duration as NA? SG
###########
############
# read in Scholz data: 156 raw observation, 7 chem pass ce filter
# only 5 chem at LC50 96 hours, and fewer at other timepoints, we can try and pair with su data for 48 hours and 96 hours.
zf_scholz <- readxl::read_excel("data/raw/Zebrafish/Scholz et al 2016/annex2_fet_en.xlsx",
                                na = c("NA", ""),
                                sheet = "ZFET final for AFT comparison") %>%
  dplyr::mutate(cas = stringr::str_replace_all(CAS, pattern = "-", replacement = "")) %>% # fix cas to match ce
  dplyr::filter(cas %in% all_ce_cas$cas) %>%
  dplyr::mutate(source = "scholz et al. 2016") #%>%
  dplyr::select(cas = 5, `LC50 0-24 hpf (mg/L)`:`LC50 24-120 hpf (mg/L)`, source, source2 = 92) %>%
  dplyr::select(cas, source, source2, `LC50 0-48 hpf (mg/L)`, `LC50 0-96 hpf (mg/L)`) #%>% 
  #dplyr::filter(!is.na(`LC50 0-96 hpf (mg/L)`))
length(unique(zf_scholz$cas))

# Can we pair su and scholz  data to increase the number of compounds with LC50s at 96 hours?
# need to convert units too!

#===================================================================#
# Step 6: Get EnviroTox DB from source
#===================================================================#
# getting latest export of database from SG
enviroTox <- readxl::read_excel("data/raw/envirotox_20240729124104.xlsx",
                   na = c("NA", ""),
                   sheet = "test") %>%
  dplyr::mutate(strain = NA_character_)

# process the data for joining
names(orig)
glimpse(orig)
glimpse(enviroTox)

# Need to clean up the source2, source columns
proc_enviroTox <- enviroTox %>%
  dplyr::filter(CAS %in% cas_df$cas) %>% # filter to Widmayer2022 cas # cas_numbers
  dplyr::mutate(source = "EnviroTox DB",
                duration_d = `Duration (hours)` / 24) %>%
  dplyr::left_join(., cas_df, by = c("CAS" = "cas")) %>%
  dplyr::select(cas = CAS, chem_name, group = `Trophic Level`, latin_name = `Latin name`,
                strain, test_type = `Test type`, test_statistic = `Test statistic`, duration_d,
                effect_value = `Effect value`, effect_unit = Unit, source, endpoint = Effect)

# test names match perfectly to orig
identical(names(orig), names(proc_enviroTox))

#===================================================================#
# Step 7: Pull the Boyd data from the EHP file and keep all good chems
#===================================================================#
# read in original boyd data from Scott: 763 raw observations, 
boyd_orig <- readxl::read_excel("data/raw/Data_Boyd_EnviroTox.xlsx", na = c("NA", "")) %>%
  dplyr::filter(source == "Boyd et al")
  length(unique(boyd_orig$cas))

# read in boyd data, convert AC50 to mg/L, filter out extrapolated AC50s.
boyd <- data.table::fread("data/raw/Boyd/ehp.1409645.s002.csv") %>%
  dplyr::slice(-1) %>% # drop first row
  janitor::row_to_names(1) %>% # set names
  dplyr::mutate(CAS = stringr::str_replace_all(CASRN, pattern = "-", replacement = ""),
                CAS = as.numeric(CAS)) %>% # fix CAS for joining convert NOCAS to NA
  dplyr::left_join(cas_df, by = c("CAS" = "cas")) %>%
  dplyr::mutate(AC50_uM = as.numeric(`C.elegans AC50 (Hill function)`),
                AC50 = (AC50_uM / 1e6) * AVERAGE_MASS * 1000) %>%
  dplyr::filter(AC50_uM < 200.5) # filter to AC50s at or below the maximum tested concentration of 200uM

# process boyd data to join by adding other variables
proc_boyd <- boyd %>%
  dplyr::mutate(group = "NEMATODE_COPAS1",
                latin_name = "Caenorhabditis elegans",
                strain = NA_character_,
                test_type = NA_character_,
                test_statistic = "AC50",
                duration_d = NA_real_, # do we have a duration of exposure? yes right?
                effect_unit = "mg/L",
                source = "Boyd et al",
                endpoint = "Growth") %>%
  select(cas = CAS, chem_name, group:duration_d, effect_value = AC50, effect_unit:endpoint)

# test names match perfectly to orig
identical(names(orig), names(proc_boyd))
#===================================================================#
# Step 8: Read NIEHS_ICE folder
#===================================================================#
# read in NIEHS_ICE Folder rat derm data: 922 raw observations, 10 chems match ce cas, 4 chems pass response modifier filter
ice_derm <- readxl::read_excel("data/raw/NIEHS_ICE/ICE_Acute_Dermal_Toxicity.xlsx", na = c("NA", "")) %>%
  dplyr::mutate(cas = stringr::str_replace_all(CASRN, pattern = "-", replacement = "")) %>% # fix cas to match ce
  dplyr::filter(cas %in% cas_df$cas) %>% 
  dplyr::filter(is.na(`Response Modifier`))
length(unique(ice_derm$CASRN))

# read in NIEHS_ICE Folder rat oral data: 12958 ra observations, 13 chem pass ce filter
ice_oral <- readxl::read_excel("data/raw/NIEHS_ICE Folder/ICE_Acute_Oral_Toxicity.xlsx", na = c("NA", "")) %>%
  dplyr::mutate(cas = stringr::str_replace_all(CASRN, pattern = "-", replacement = "")) %>% # fix cas to match ce
  dplyr::filter(cas %in% cas_df$cas) %>% # 13 chem pass
  dplyr::filter(is.na(`Response Modifier`)) %>%# 11 chem pass
  dplyr::filter(Endpoint == "LD50") # 11 chem pass
length(unique(ice_oral$CASRN))

# read in NIEHS_ICE Folder rat dart data: 3629 raw observations, 8 chems match ce cas.
ice_dart <- readxl::read_excel("data/raw/NIEHS_ICE Folder/ICE_DART.xlsx", na = c("NA", "")) %>%
  dplyr::mutate(cas = stringr::str_replace_all(CASRN, pattern = "-", replacement = "")) %>% # fix cas to match ce
  dplyr::filter(cas %in% cas_df$cas) # 8 chem pass filter
length(unique(ice_dart$CASRN))

#===================================================================#
# Step 9: Shape to original data format
#===================================================================#
orig <- data.table::fread("data/processed/03_clean.csv") 
names(orig)  
unique(orig$source)

# get molecular weight data
# Get comptox MW data to convert ToxCast to w/v unit
mw_data <- readxl::read_excel("data/processed/toxcast_mw.xlsx") %>%
  dplyr::mutate(INPUT = as.character(INPUT)) %>%
  dplyr::distinct(.keep_all = T) %>%
  dplyr::mutate(AVERAGE_MASS = ifelse(CASRN == "8018-01-7", 541.08, AVERAGE_MASS),  # fix Mancozeb
                cas = stringr::str_replace_all(CASRN, pattern = "-", replacement = "")) %>%
  dplyr::select(cas, AVERAGE_MASS)

# 1 joining widmayer ce data
ce_join <- ce %>%
  dplyr::select(names(orig))
# 2 joining boyd ce data
boyd_join <- boyd %>%
  dplyr::select(names(orig))
# 3 joingin Karmaus empirical rat data
kar_join <- kar %>%
  dplyr::left_join(., all_cas) %>%
  dplyr::left_join(., mw_data) %>%
  dplyr::mutate(group = "RAT_EMPIRICAL",
                latin_name = "Rat",
                strain = NA_character_,
                test_type = NA_character_,
                test_statistic = "LD50",
                duration_d = NA_integer_,
                effect_value = LD50_mgkg,
                effect_unit = "mg/kg",
                source = "Karmaus 2022",
                endpoint = "Mortality") %>%
  dplyr::select(names(orig))
# 4 joining NIEHS_ICE data from Rat and Rabbit
names(ice_derm)
glimpse(ice_derm)
ice_derm_join <-  ice_derm 

names(ice_oral)
glimpse(ice_oral)
ice_oral_join <-  ice_oral %>%
  
  
  
  
  names(ice_oral)
names(ice_dart)

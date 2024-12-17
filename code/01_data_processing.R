#!/usr/bin/env Rscript
library(tidyverse)

# set working dir to base directory of repository
setwd(paste0(dirname(rstudioapi::getActiveDocumentContext()$path),"/.."))

#==============================================================================#
# Testing: Get original cleaned data and match coded data processing
#==============================================================================#
orig <- data.table::fread("data/processed/03_clean.csv") %>%
  dplyr::mutate(orig_source = NA_character_, .after = source) %>%
  dplyr::mutate(QC = NA_character_)

#===================================================================#
# Step 1: Read in Widmayer data and make list of 18 chemicals
#===================================================================#
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

# read in raw C. elegans data from SG: 18 chemicals from 22 total raw, could add anthelmintics
ce <- readxl::read_excel("data/raw/Data_Andersen_All.xlsx", na = c("NA", "")) %>%
  dplyr::filter(source == "Widmayer 2022")  %>%
  dplyr::mutate(effect_unit = ifelse(effect_unit == "mg/l", "mg/L", effect_unit)) %>% # fix mg/l
  dplyr::mutate(duration_d = ifelse(source == "Widmayer 2022", 2, duration_d),
                source = "Widmayer et al. 2022") %>%
  dplyr::relocate(endpoint, .after = source) %>% # flip to match output
  dplyr::mutate(orig_source = "Widmayer et al. 2022", .after = source) %>%
  dplyr::mutate(QC = NA_character_) # add QC to match
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
#===================================================================#
# Step 2: Pull the Boyd data from the EHP file and keep all good chems
#===================================================================#
# # read in original boyd data from Scott: 763 raw observations, 
# boyd_orig <- readxl::read_excel("data/raw/Data_Boyd_EnviroTox.xlsx", na = c("NA", "")) %>%
#   dplyr::filter(source == "Boyd et al")
# length(unique(boyd_orig$cas))
# # get list of boyd chems 467 pass filter, 459 not in Widmayer (8 overlap)
# boyd_chem <- boyd %>%
#   dplyr::distinct(CAS, .keep_all = T) %>% # remvoe dups
#   dplyr::filter(!is.na(CAS)) %>% # remove NOCAS
#   dplyr::select(cas = CAS, chem_name = Name) %>%
#   dplyr::filter(!(cas %in% cas_df$cas)) # remove chems that overlap with widmayer for now
# 
# # output the cas numbers, these are loaded into toxcast dashboard to output molecular weights
# rio::export(boyd_chem, "data/processed/boyd_chem.csv")

# Molecular weights are added to complete cas list from boyd and adding widmayer
all_ce_cas <- readxl::read_excel("data/raw/toxcast_boyd_mw.xlsx") %>%
  dplyr::select(CASRN = INPUT, AVERAGE_MASS, PREFERRED_NAME) %>%
  dplyr::mutate(cas = stringr::str_replace_all(CASRN, pattern = "-", replacement = ""),
                cas = as.numeric(cas)) %>%
  dplyr::filter(!is.na(AVERAGE_MASS)) %>% # remove NO MW
  dplyr::mutate(AVERAGE_MASS = ifelse(cas == 27176870, 326.49, AVERAGE_MASS)) %>%
  dplyr::select(cas, chem_name = PREFERRED_NAME, AVERAGE_MASS) %>%
  dplyr::filter(!(cas %in% cas_df$cas)) %>%
  dplyr::bind_rows(cas_df) # add back the Widmayer 2022 compounds

# read in boyd data, convert AC50 to mg/L, filter out extrapolated AC50s.
boyd <- #data.table::fread("data/raw/Boyd/ehp.1409645.s002.csv") %>%
  readxl::read_excel("data/raw/Boyd/ehp.1409645.s002.xlsx", sheet = "Excel Table S2") %>%
  janitor::row_to_names(1) %>% # set names
  dplyr::mutate(CAS = stringr::str_replace_all(CASRN, pattern = "-", replacement = "")) %>%
  dplyr::filter(!(grepl(CAS, pattern = "NOCAS"))) %>%
  dplyr::mutate(CAS = as.numeric(CAS)) %>% # fix CAS for joining convert NOCAS to NA
  dplyr::left_join(all_ce_cas, by = c("CAS" = "cas")) %>%
  dplyr::mutate(AC50_uM = as.numeric(`C.elegans AC50 (Hill function)`),
                AC50 = (AC50_uM / 1e6) * AVERAGE_MASS * 1000) %>%
  dplyr::mutate(QC = ifelse(AC50_uM < 200.5, "PASS", "FAIL")) #AC50s at or below the maximum tested concentration of 200uM

# NO MOLECULAR WEIGHT - COULD RECOVER MANUALLY
boyd_noCAS <- boyd %>%
  dplyr::filter(is.na(AVERAGE_MASS))

# process boyd data to join by adding other variables
proc_boyd <- boyd %>%
  dplyr::filter(!is.na(chem_name)) %>% #drop no molecular weight chems
  dplyr::mutate(group = "NEMATODE_COPAS1",
                latin_name = "Caenorhabditis elegans",
                strain = NA_character_,
                test_type = NA_character_,
                test_statistic = "AC50",
                duration_d = NA_real_, # do we have a duration of exposure? yes right?
                effect_unit = "mg/L",
                source = "Boyd et al. 2016",
                orig_source = "Boyd et al. 2016",
                endpoint = "Growth") %>%
  select(cas = CAS, chem_name, group:duration_d, effect_value = AC50, effect_unit:endpoint, QC)
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

# process the karmaus data for joining
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
                source = "Karmaus et al. 2022",
                orig_source = NA_character_,
                endpoint = "Mortality",
                QC = NA_character_) %>%
  dplyr::select(cas, chem_name, group:endpoint, QC) 
#===================================================================#
# Step 3: Read and process Comptox folder - RAT test model predictions
#===================================================================#
# read in the output from comptox: 30891 raw observations, 14 pass TEST filter
comptox <- readxl::read_excel("data/raw/Comptox/20240717_comptox_download.xlsx", na = c("NA", ""),
                              sheet = "Toxval Details") %>%
  dplyr::filter(SOURCE == "TEST") # filter to acute oral LD50 TEST data

# process the data for joining, 9 chems pass ce filter
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
                orig_source = LONG_REF,
                endpoint = "Mortality",
                QC = NA_character_) %>%
  dplyr::select(cas, chem_name, group:endpoint, QC) 
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
  #dplyr::filter(chnm == "Methyl isothiocyanate" | chnm == "Butafenacil" | casn == "NOCAS_34742") %>% # test milbamectin no AC50 enven though it is published in padilla 2012? OK, no AVERAGE_MASS so lost
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
proc_zf1 <- zf1 %>%
  dplyr::mutate(cas = as.integer(stringr::str_replace_all(casn, pattern = "-", replacement = "")),
                group = ifelse(aeid == 1372, "TC_ZF_Tanguay", "TC_ZF_Padilla"), # updated group name to reflect toxcast_zebrafish_source
                latin_name = "Danio rerio",
                strain = NA_character_,
                effect_unit = "mg/L",
                source = "TOXCAST",
                orig_source = ifelse(group == "TC_ZF_Tanguay", "Tanguay_lab_OSU", "Padilla_lab_EPA"),
                chem_name = case_when(chnm == "2,4-Dichlorophenoxyacetic acid" ~ "2,4-D",
                                      chnm == "Methylmercuric(II) chloride" ~ "Methylmercury chloride",
                                      TRUE ~ chnm),
                test_type = NA,
                duration_d = timepoint_hr / 24,
                endpoint = ifelse(assay_component_endpoint_name == "Tanguay_ZF_120hpf_MORT", "Mortality", "TERATOSCORE")) %>% # add variables to join
  dplyr::select(cas, chem_name, group, latin_name, strain, test_type,
                duration_d, effect_unit, source, orig_source, endpoint,
                AC50_mg_L, AC20_mg_L, AC10_mg_L, AC10_mg_L, AC5_mg_L, BMD_mg_L, ACC_mg_L,
                hitc, fitc, QC, flag.length) %>% # select what's needed
  tidyr::pivot_longer(cols = AC50_mg_L:ACC_mg_L,
                      names_to = "test_statistic", values_to = "effect_value") %>% # reshape long
    dplyr::mutate(test_statistic = stringr::str_replace(test_statistic, pattern = "_mg_L", replacement = "")) %>% # make clean test_stat
  dplyr::select(cas, chem_name, group, latin_name, strain, test_type,
                test_statistic, duration_d, effect_value, effect_unit, source, orig_source, endpoint, QC) %>% # select proper order
  #dplyr::filter(!is.na(effect_value), QC == "PASS") %>% # filter by QC
  dplyr::filter(!is.na(effect_value)) %>%
  dplyr::select(cas, chem_name, group:endpoint, QC) # set proper order of variables

# chemical counts
proc_zf1 %>%
  dplyr::group_by(group) %>%
  dplyr::summarise(chem_n = length(unique(cas)))
#===================================================================#
# Step 4: Read Zebrafish folder to get Su data
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

# process the data for joining
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
                orig_source = Reference,
                endpoint = "Mortality",
                QC = NA_character_) %>%
  dplyr::select(cas, chem_name, group:endpoint, QC)
#==================================================================#
# Step 5: Read Zebrafish folder to get Scholz et al. 2016 data
#===================================================================
# read in Scholz data: 11 time points, 4 time points with > 20 chems
zf_scholz <- readxl::read_excel("data/raw/Zebrafish/Scholz et al 2016/annex2_fet_en.xlsx",
                                na = c("NA", ""),
                                sheet = "ZFET final for AFT comparison") %>%
  dplyr::mutate(cas = as.integer(stringr::str_replace_all(CAS, pattern = "-", replacement = ""))) %>% # fix cas to match ce
  dplyr::filter(cas %in% all_ce_cas$cas) %>%
  dplyr::mutate(source = "scholz et al. 2016")

# count numbers of chems for each timepoint
zf_scholz %>%
  dplyr::select(cas, `LC50 0-24 hpf (mg/L)`:`LC50 24-120 hpf (mg/L)`, source = 92) %>%
  tidyr::pivot_longer(cols = `LC50 0-24 hpf (mg/L)`:`LC50 24-120 hpf (mg/L)`) %>%
  dplyr::filter(!is.na(value)) %>%
  dplyr::group_by(name) %>%
  dplyr::summarise(nchem = length(unique(cas)))
  
# All the LC50 0-120 (140) hpf (mg/L) values are from Padilla et al. 2012...
# I can't recreate know how scholz calculated the padilla et al. 2012 data into an LC50, best left as separate data source
# The scholz data are not identical to toxcast data, which makes sense because toxcast has TERATOSCORE for Padilla not mortality, still they are close
scholz_v_toxcast <- zf_scholz %>%
  dplyr::select(`Common name`, cas, `LC50 0-120 (140) hpf (mg/L)`, source = 92) %>%
  dplyr::filter(source == "Padilla et al. 2012") %>%
  dplyr::filter(!is.na(`LC50 0-120 (140) hpf (mg/L)`)) %>%
  dplyr::left_join(dplyr::filter(proc_zf1, group == "TC_ZF_Padilla", test_statistic == "AC50") %>% dplyr::select(cas, effect_value, duration_d)) %>%
  dplyr::group_by(cas) %>%
  dplyr::mutate(mean_ev = mean(effect_value),
                geometric_mean_ev = (prod(effect_value)^(1/n()))) %>% # datasets are often different and neither matches the table 1 in padilla et al. 2012
  dplyr::arrange(cas)

# are the LC50 0-120 (140) hpf (mg/L) and LC50 0-120 hpf (mg/L) values similar?
# No not very often, we can flag the padilla with a QC flag b/c I don't know what 120 hpf (140) means anyway.
scholz_v_120_test <- zf_scholz %>%
  dplyr::select(cas, chnm = `Common name`, `LC50 0-120 (140) hpf (mg/L)`, `LC50 0-120 hpf (mg/L)`, orig_source = 92) %>%
  dplyr::mutate(value = ifelse(orig_source == "Padilla et al. 2012", `LC50 0-120 (140) hpf (mg/L)`, `LC50 0-120 hpf (mg/L)`),
                group = ifelse(orig_source == "Padilla et al. 2012", "Padilla et al. 2012", "other")) %>%
  dplyr::filter(!is.na(value)) %>%
  dplyr::group_by(chnm, group) %>%
  dplyr::mutate(mean = mean(value)) %>%
  dplyr::distinct(chnm, group, mean) %>%
  dplyr::ungroup() %>%
  tidyr::pivot_wider(names_from = group, values_from = mean) %>%
  dplyr::mutate(diff = `Padilla et al. 2012` - other) %>%
  dplyr::filter(!is.na(diff))

ggplot(scholz_v_120_test) +
  aes(x = `Padilla et al. 2012`, y = other) +
  ggplot2::geom_point()

# process the scholz data
proc_zf_scholz <- zf_scholz %>%
  dplyr::select(`Common name`, cas,
                `LC50 0-120 (140) hpf (mg/L)`,
                `LC50 0-120 hpf (mg/L)`,
                `LC50 0-48 hpf (mg/L)`,
                `LC50 0-96 hpf (mg/L)`,
                orig_source = 92) %>%
  dplyr::filter(cas %in% all_ce_cas$cas) %>%
  tidyr::pivot_longer(col = `LC50 0-120 (140) hpf (mg/L)`:`LC50 0-96 hpf (mg/L)`, values_to = "effect_value") %>%
  dplyr::filter(!is.na(effect_value)) %>%
  dplyr::mutate(group = "SCHOLZ_ZF",
                latin_name = "Danio rerio",
                strain = NA_character_,
                test_type = NA_character_,
                test_statistic = "LC50",
                duration_d = dplyr::case_when(name == "LC50 0-120 (140) hpf (mg/L)" ~ 5,
                                              name == "LC50 0-120 hpf (mg/L)" ~ 5,
                                              name == "LC50 0-48 hpf (mg/L)" ~ 2,
                                              name == "LC50 0-96 hpf (mg/L)" ~ 4,
                                              TRUE ~ NA_real_),
                effect_unit = "mg/L",
                source = "Scholz et al. 2016",
                endpoint = "Mortality",
                QC = dplyr::case_when(orig_source == "Padilla et al. 2012" ~ "FAIL", # set Padilla LC50s to FAIL b/c they are suspect
                                      TRUE ~ NA_character_)) %>%
  dplyr::left_join(all_ce_cas) %>%
  dplyr::select(cas, chem_name, group:duration_d, effect_value, effect_unit:source, orig_source, endpoint, QC)
#===================================================================#
# Step 6: Get EnviroTox DB from source
#===================================================================#
# getting latest export of database from SG
enviroTox <- readxl::read_excel("data/raw/envirotox_20240729124104.xlsx",
                   na = c("NA", ""),
                   sheet = "test") %>%
  dplyr::mutate(strain = NA_character_)

# Need to clean up the source2, source columns
proc_enviroTox <- enviroTox %>%
  dplyr::filter(CAS %in% all_ce_cas$cas) %>% # filter to cas_numbers
  dplyr::mutate(source = "EnviroTox DB",
                duration_d = `Duration (hours)` / 24,
                QC = NA_character_) %>%
  dplyr::left_join(., all_ce_cas, by = c("CAS" = "cas")) %>%
  dplyr::select(cas = CAS, chem_name, group = `Trophic Level`, latin_name = `Latin name`,
                strain, test_type = `Test type`, test_statistic = `Test statistic`, duration_d,
                effect_value = `Effect value`, effect_unit = Unit, source, orig_source = Source, endpoint = Effect, QC) 

# could clean endpoints a bit - commenting out for now, but an example would be:  
# dplyr::mutate(dplyr::case_when(endpoint %in% c("Mortality/Growth",
#                                                  "Mortality, Mortality",
#                                                  "Mortality, Survival") ~ "Mortality",
#                                  endpoint == "Immobilization: Change in the failure to respond or lack of movement after mechanical stimulation." ~ "Immobility",
#                                  endpoint == "Intoxication, Immobile" ~ "Immobility",
#                                  TRUE ~ endpoint)

# look at many endpoints
enviroTox_eps <- proc_enviroTox %>%
  dplyr::group_by(endpoint) %>%
  dplyr::summarize(n_obs = n()) %>%
  dplyr::arrange(endpoint) %>%
  dplyr::filter(n_obs > 30) # filter to endpoints with many observations.
nrow(enviroTox_eps)
enviroTox_eps$endpoint

# look at original endpoints - ok lots there too
orig_eps <- orig %>%
  dplyr::filter(source == "EnviroTox DB") %>%
  dplyr::distinct(endpoint) %>%
  dplyr::arrange(endpoint) %>%
  dplyr::pull()
length(orig_eps)

# shorten the endpoint list
proc_enviroTox2 <- proc_enviroTox %>%
  dplyr::filter(endpoint %in% enviroTox_eps$endpoint)
#===================================================================#
# Step 7: Read NIEHS_ICE folder
#===================================================================#
# read in NIEHS_ICE Folder rat derm data: 6 chems match ce cas and pass filters - dropped from study
ice_derm <- readxl::read_excel("data/raw/NIEHS_ICE/Acute_Dermal_Toxicity.xlsx", na = c("NA", "")) %>%
  dplyr::mutate(cas = stringr::str_replace_all(CASRN, pattern = "-", replacement = "")) %>% # fix cas to match ce
  dplyr::filter(cas %in% all_ce_cas$cas) %>% 
  dplyr::filter(Response != "Not available" & Endpoint == "LD50") %>%
  dplyr::filter(is.na(`Response Modifier`))
length(unique(ice_derm$CASRN))

# read in NIEHS_ICE Folder rat oral data: 16167 raw observations, 4283 in ce chems, 2295 LD50s, 1414 pass filters
ice_oral <- readxl::read_excel("data/raw/NIEHS_ICE/Acute_Oral_Toxicity.xlsx", na = c("NA", "")) %>%
  dplyr::mutate(cas = stringr::str_replace_all(CASRN, pattern = "-", replacement = "")) %>% # fix cas to match ce
  dplyr::filter(cas %in% all_ce_cas$cas) %>% 
  dplyr::filter(Endpoint == "LD50") %>% # 11 chem pass
  dplyr::mutate(QC = ifelse(is.na(`Response Modifier`), "PASS", "FAIL"))
length(unique(ice_oral$CASRN))

# process it
proc_ice_oral <- ice_oral %>%
  dplyr::select(cas, test_type = Assay, test_statistic = Endpoint, effect_value = Response,
                effect_unit = `Response Unit`, orig_source = Reference, QC) %>%
  dplyr::mutate(cas = as.integer(cas)) %>%
  dplyr::left_join(dplyr::select(all_ce_cas, cas, chem_name)) %>%
  dplyr::mutate(effect_value = as.numeric(effect_value),
                group = "NIEHS_RAT",
                latin_name = "Rat",
                strain = NA_character_,
                duration_d = NA_real_,
                source = "NIEHS_ICE",
                endpoint = "Acute Oral Toxicity") %>%
  dplyr::select(names(proc_enviroTox))

# read in NIEHS_ICE Folder rat dart data: 138326 raw observations, 77302 match ce cas, 2168 match "LOEL" endpoint
ice_dart <- readxl::read_excel("data/raw/NIEHS_ICE/DART.xlsx", na = c("NA", "")) %>%
  dplyr::mutate(cas = stringr::str_replace_all(CASRN, pattern = "-", replacement = "")) %>% # fix cas to match ce
  dplyr::filter(cas %in% all_ce_cas$cas) %>%# 8 chem pass filter
  dplyr::filter(Endpoint == "LOEL" & Route == "Oral" & Lifestage == "Adult",
                Assay %in% c("DART, Developmental malformation",
                                                   "DART, In life observation",
                                                   "DART, Reproductive performance")) # males and females not distinguished here

# process it
proc_ice_dart <- ice_dart %>%
  dplyr::select(cas, Species, strain = Strain, endpoint = Assay, test_statistic = Endpoint, effect_value = Response,
                effect_unit = `Response Unit`, orig_source = Reference) %>%
  dplyr::mutate(cas = as.integer(cas)) %>%
  dplyr::left_join(dplyr::select(all_ce_cas, cas, chem_name)) %>%
  dplyr::mutate(effect_value = as.numeric(effect_value),
                group = dplyr::case_when(Species == "Rat" ~ "NIHEHS_RAT",
                                         Species == "Rabbit" ~ "NIHEHS_RABBIT",
                                         Species == "Mouse" ~ "NIHEHS_MOUSE"),
                latin_name = Species,
                strain = NA_character_,
                duration_d = NA_real_,
                test_type = NA_character_,
                source = "NIEHS_ICE",
                QC = NA_character_) %>%
  dplyr::select(names(proc_enviroTox))
#===================================================================#
# Step 9: Shape to original data format
#===================================================================#
# bind the data sources together
bind_dat <- dplyr::bind_rows(ce, proc_boyd, proc_comptox, proc_enviroTox2,
                             proc_ice_dart, proc_ice_oral, proc_kar,
                             proc_zf_scholz, proc_zf_su, proc_zf1)

# save data
readr::write_csv(bind_dat, file = "data/processed/00_data.csv")

# summarize the data
n.sp <- length(unique(bind_dat$latin_name))
n.drugs <- length(unique(bind_dat$cas))
n.drugs.per.sp <- bind_dat %>%
  dplyr::group_by(latin_name) %>%
  dplyr::mutate(n.drugs.in.sp = length(unique(cas))) %>%
  dplyr::ungroup() %>%
  dplyr::distinct(latin_name, .keep_all = T) %>%
  dplyr::mutate(avg.n.drugs.sp = mean(n.drugs.in.sp),
                sd.n.drugs.sp = sd(n.drugs.in.sp)) #%>%
p <- ggplot(n.drugs.per.sp) +
  aes(x = n.drugs.in.sp) +
  geom_histogram(bins = 1000) +
  theme_bw()
p
n.groups <- length(unique(bind_dat$group))
n.endpoints <- length(unique(bind_dat$endpoint))
n.test_statistics <- length(unique(bind_dat$test_statistic))

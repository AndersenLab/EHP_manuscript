#!/usr/bin/env Rscript
library(tidyverse)

# set working dir to base directory of repository
setwd(paste0(dirname(rstudioapi::getActiveDocumentContext()$path),"/.."))

#------------------------------------------------------------------------------#
# part 1: # read in Karmaus data read in the Point Estimates Only sheet
#------------------------------------------------------------------------------#
Karmaus <-  readxl::read_excel(path = "data/raw/raw2/Karmaus Folder/toxsci-21-0357-File010.xlsx", sheet = "Point Estimates Only") %>%
  dplyr::select(2:3)

#------------------------------------------------------------------------------#
# part 2: # read in the ICE data
#------------------------------------------------------------------------------#
ice_oral <- readxl::read_excel(path = "data/raw/raw2/NIEHS_ICE Folder/Acute_Oral_Toxicity.xlsx", sheet = "Data") %>%
  dplyr::select(CASRN, Endpoint, Response, res_mod = `Response Modifier`, res_unit = `Response Unit`, Mixture) %>%
  dplyr::filter(res_unit == "mg/kg" & is.na(res_mod) & Mixture == "Chemical", Endpoint == "LD50")

ice_dart <- readxl::read_excel(path = "data/raw/raw2/NIEHS_ICE Folder/DART.xlsx", sheet = "Data") %>%
  dplyr::select(CASRN, Endpoint, Response, res_unit = `Response Unit`, Assay) %>%
  dplyr::filter(Endpoint == "LOEL")
unique(ice_dart$Assay)
# We don't want to mix the assays together when estimating the chemical endpoint mean so we can setup a unique ID
# for the assay and endpoint.
table(ice_dart$Assay)
# use all assays for now, most will be filtered out anyway

#------------------------------------------------------------------------------#
# part 3: # read in the Rat data
#------------------------------------------------------------------------------#
# get the CAS from the nematode data use this to 
cas <- readxl::read_excel(path = "~/repos/EHP_manuscript/data/raw/Data_Andersen_All.xlsx") %>%
  dplyr::distinct(cas) 
rio::export(cas, file = "data/raw/cas.txt")

# check to see if manual export from web is the same as old data
old_rat <- readxl::read_excel(path = "~/repos/EHP_manuscript/data/raw/Data_Andersen_All.xlsx") %>%
  dplyr::filter(group == "RAT_TEST")

# get the new data
new_rat_raw <- readxl::read_excel(path = "data/raw/raw2/rat/CCD-Batch-Search_2024-06-04_07_04_22.xlsx", sheet = "Toxval Details") 
# process it
new_rat <- new_rat_raw %>%
  dplyr::filter(SOURCE == "TEST") %>%
  dplyr::select(chem_name = NAME, casrn = CASRN, source = SOURCE, Response = TOXVAL_NUMERIC_ORIGINAL, res_unit = TOXVAL_UNITS,
                type = TOXVAL_TYPE, LONG_REF) %>%
  dplyr::mutate(cas = stringr::str_replace_all(casrn, pattern = "-", replacement = ""),
                cas = as.double(cas))

# join the old and new to see what's up
test_join <- full_join(old_rat, new_rat, by = c("cas" = "cas"))
missing <- test_join %>%
  dplyr::filter(is.na(chem_name.y)) 

missing.vec <- missing$chem_name.x
missing.vec
# output these for scott
old_out <- old_rat
rio::export(old_out, file = "data/raw/old_rat_test.csv")
new_out <- new_rat
rio::export(new_out, file = "data/raw/new_rat_test.csv")

# nickel is in the raw data but not in the TEST data. 
TRUE %in% (grepl(new_rat_raw$CASRN, pattern = "7718-54-9"))
grepl(new_rat_raw$s, pattern = "7718549")
glimpse(new_rat)


#------------------------------------------------------------------------------#
# part 4: # read in the Zebra data (tanguay pedilla, other)
#------------------------------------------------------------------------------#



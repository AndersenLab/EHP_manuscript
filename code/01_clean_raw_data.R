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
  dplyr::mutate(effect_unit = ifelse(effect_unit == "mg/l", "mg/L", effect_unit)) %>% # fix mg/l
  dplyr::mutate(duration_d = ifelse(source == "Widmayer 2022", 2, duration_d))

# look for duplication since the row numbers above don't make sense. 262 rows duplicated in proc_dat3, some up to 4 times.
# Let's just trust it for now.
dup.test <- proc_dat3 %>%
  dplyr::mutate(id = paste0(cas, chem_name, group, latin_name, strain, test_type, test_statistic, duration_d, effect_value,
                            effect_unit, endpoint, source)) %>%
  dplyr::group_by(id) %>%
  dplyr::mutate(n = n()) %>%
  dplyr::filter(n > 1) %>%
  dplyr::distinct(id, .keep_all = T)

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

# clean up the endpoints - Let's add the changes here so we can export them cleanly for manuscript
eps <- tibble::tibble(endpoint = unique(join_dat$endpoint), proc_ep = NA_character_) %>%
  dplyr::mutate(proc_ep = dplyr::case_when(endpoint %in% c("Mortality/Growth",
                                                         "Mortality, Mortality",
                                                         "Mortality, Survival") ~ "Mortality",
                                           endpoint == "Immobilization: Change in the failure to respond or lack of movement after mechanical stimulation." ~ "Immobility",
                                           endpoint == "Intoxication, Immobile" ~ "Immobility",
                                           TRUE ~ endpoint))


# add the endpoint changes to joined data here
join_dat_proc <- join_dat %>%
  dplyr::left_join(eps) %>%
  dplyr::select(-endpoint) %>%
  dplyr::rename(endpoint = proc_ep)

# # SG asked about mortality and intoxication for inverts - can we merge 2day mortality and 2day intoxication?
# merge1 <- join_dat_proc %>%
#   dplyr::filter(group == "INVERT") %>%
#   dplyr::group_by(endpoint, duration_d) %>%
#   dplyr::mutate(n.in.group = n()) %>%
#   dplyr::distinct(latin_name, endpoint, duration_d, n.in.group) %>%
#   dplyr::filter(endpoint == "Intoxication")
  
# save the joined data
readr::write_csv(join_dat_proc, file = "data/processed/03_clean.csv")

# summarize the data
n.sp <- length(unique(join_dat_proc$latin_name))
n.drugs <- length(unique(join_dat_proc$cas))
n.drugs.per.sp <- join_dat_proc %>%
  dplyr::group_by(latin_name) %>%
  dplyr::mutate(n.drugs.in.sp = length(unique(cas))) %>%
  dplyr::ungroup() %>%
  dplyr::distinct(latin_name, .keep_all = T) %>%
  dplyr::mutate(avg.n.drugs.sp = mean(n.drugs.in.sp),
                sd.n.drugs.sp = sd(n.drugs.in.sp)) #%>%
p <- ggplot(n.drugs.per.sp) +
  aes(x = n.drugs.in.sp) +
  geom_histogram() +
  theme_bw()
p
n.groups <- length(unique(join_dat_proc$group))
n.endpoints <- length(unique(join_dat_proc$endpoint))
n.test_statistics <- length(unique(join_dat_proc$test_statistic))
#===================================================================#
# Step 2: Read in data from ToxCast with tcpl package
# https://cran.r-project.org/web/packages/tcpl/tcpl.pdf
#===================================================================#
## Store the current config settings, so they can be reloaded at the end
## of the examples
conf_store <- tcpl::tcplConfList()
tcpl::tcplConfExample()

## Load all of level 0 for multiple-concentration data, note 'mc' is the
## default value for type
test <- tcpl::tcplLoadData(lvl = 0)
## Load all of level 1 for single-concentration
tcplLoadData(lvl = 1, type = "sc")

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

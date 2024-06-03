#!/usr/bin/env Rscript
library(tidyverse)
library(readxl)
library(expss)

# set working dir to base directory of repository
setwd(paste0(dirname(rstudioapi::getActiveDocumentContext()$path),"/.."))

# load functions
source("code/functions.R")

# load ToxCast data from Glaberman
load("data/processed/mc5-6_winning_model_fits-flags_invitrodb_v4_1_SEPT2023.Rdata")

# read in the data dictionary to help join toxcast
dat.dic <- data.table::fread("data/raw/toxcast_data_dictionaryv1.csv")

#------------------------------------------------------------------------------#
# Get ToxCast data that match the nematode toxicants and apply
# simple filters based on discussion with Katie
#------------------------------------------------------------------------------#
# Define CAS numbers of interest
cas_numbers <- c("10108-64-2", "114-26-1", "115-09-3", "115-86-6", "116-06-3",
                 "121-75-5", "16752-77-5", "175013-18-0", "1912-24-9", "2921-88-2",
                 "4685-14-7", "5234-68-4", "63-25-2", "7646-85-7", "7718-54-9",
                 "7761-88-8", "8018-01-7", "94-75-7")

# Get comptox MW data to convert ToxCast to w/v unit
mw_data <- read_excel("data/processed/toxcast_mw.xlsx") %>%
  dplyr::mutate(INPUT = as.character(INPUT)) %>%
  dplyr::distinct(.keep_all = T) %>%
  dplyr::mutate(AVERAGE_MASS = ifelse(CASRN == "8018-01-7", 541.08, AVERAGE_MASS)) %>% # fix Manxozeb
  dplyr::select(CASRN, AVERAGE_MASS)

# Get assay description data so we can use biological_process_target (BPT) to filter to relevant assays
assays <- readxl::read_excel("data/raw/assay_annotations_invitrodb_v4_1_SEPT2023.xlsx", sheet = "annotations_combined")

# view all the assay biological targets
sort(unique(assays$biological_process_target))

# create BPT list for filtering
bpts <- c("cell cycle", "cell death", "cell proliferation", "cell stress", "cell viability",
          "cytotoxicity", "mitochondrial depolarization", "regulation of development",
          "viability")

# filter to assays with proper biological process targets bpts
assays_bpt_filter <- assays %>%
  dplyr::filter(biological_process_target %in% bpts & organism == "human") 

# Filter to worm cas, filter to relevant assays, set QC filters and count compounds within assays
mc5_worm_bpt_filter <- assays_bpt_filter %>%
  dplyr::left_join(., mc5) %>%
  dplyr::filter(casn %in% cas_numbers) %>% # filter to worm cas
  dplyr::left_join(mw_data, by = c("casn" = "CASRN")) %>% # join MW data
  mutate(AC50_mg_L = (ac50 / 1e6) * AVERAGE_MASS * 1000, # Convert µM to mg/L {orig. ac50 * MW * 1000} yikes!
         AC20_mg_L = (ac20 / 1e6) * AVERAGE_MASS * 1000,
         AC10_mg_L = (ac10 / 1e6) * AVERAGE_MASS * 1000,
         AC5_mg_L = (ac5 / 1e6) * AVERAGE_MASS * 1000,
         BMD_mg_L = (bmd / 1e6) * AVERAGE_MASS * 1000,
         ACC_mg_L = (acc / 1e6) * AVERAGE_MASS * 1000) %>%
  dplyr::mutate(QC = ifelse(fitc > 35 & (flag.length < 5 | is.na(flag.length)), "PASS", "FAIL")) %>%
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

# make a table of assay counts by BPT and intended_target_family_sub
dist_assays <- mc5_worm_bpt_filter %>%
  dplyr::distinct(aeid, .keep_all = T)
table(dist_assays$biological_process_target)
table(dist_assays$intended_target_family_sub) # new filter is not much different from filtering to cytotoxicity in intended_target_family_sub

# make diagnositc plots for QC like Katie Paul-Friedman
mc5_2 <- mc5 %>%
  dplyr::mutate(flag.length = ifelse(is.na(flag.length), 0, flag.length))
qc1_all <- ggplot(mc5_2) +
  aes(x = flag.length, fill = as.factor(fitc)) +
  scale_fill_viridis_d() +
  geom_bar(position = 'dodge') +
  theme_bw()+
  labs(x = "Flag count", y = "Count of curve fitc", fill = "Fit category",
       subtitle = "Full ToxCast datas\nn = 3.2 million rows")
qc1_all

# filter to informative like like Katie Paul-Friedman
mc5_3 <- mc5_2 %>%
  dplyr::filter(fitc %in% c(1, 36:42))
qc1_all2 <- ggplot(mc5_3) +
  aes(x = flag.length, fill = as.factor(fitc)) +
  scale_fill_viridis_d() +
  geom_bar(position = 'dodge') +
  theme_bw()+
  labs(x = "Flag count", y = "Count of curve fitc", fill = "Fit category",
       subtitle = "Filtered ToxCast data\nn = 300 thousand rows")
qc1_all2

# Show the worm data like Katie Paul-Friedman
qc1_worms1 <- ggplot(mc5_worm_bpt_filter) +
  aes(x = flag.length, fill = as.factor(fitc)) +
  geom_bar(position = 'dodge') +
  scale_fill_viridis_d() +
  geom_bar(position = 'dodge') +
  theme_bw()+
  labs(x = "Flag count", y = "Count of curve fitc", fill = "Fit category",
       subtitle = "Worm compound filtered ToxCast data\nn = 2858 rows")
qc1_worms1

# Show the worm data like Katie Paul-Friedman
qc2_worms2 <- ggplot(mc5_worm_bpt_filter %>% dplyr::filter(QC == "PASS")) +
  aes(x = flag.length, fill = as.factor(fitc)) +
  geom_bar(position = 'dodge') +
  scale_fill_viridis_d() +
  geom_bar(position = 'dodge') +
  theme_bw()+
  labs(x = "Flag count", y = "Count of curve fitc", fill = "Fit category",
       subtitle = "Worm compound|QC filtered ToxCast data\nn = 349 rows")
qc2_worms2

# put these together
qc_grid <- cowplot::plot_grid(qc1_all, qc1_all2, qc1_worms1, qc2_worms2, labels = c("a", "b", "c", "d"), align = "hv", axis = "tblr", ncol = 4)
cowplot::ggsave2(qc_grid, file = "plots/toxcast_qc_filtering.png", width = 16, height = 4)

# plot the effect of filtering
qc_plot_bar <- ggplot(mc5_worm_bpt_filter %>% dplyr::distinct(aeid, casn, .keep_all = T)) +
  aes(x = as.character(aeid), fill = aeid_chem_pass) +
  geom_bar() +
  geom_hline(yintercept = 16, linetype = 2, linewidth = 0.25) +
  geom_hline(yintercept = 10, linetype = 2, linewidth = 0.25, color = "red") +
  labs(x = "Toxcast assay id", y = "Number of compounds shared with Widmayer data",
       subtitle = glue::glue("ToxCast relevant assays based on biological process target n={mc5_worm_bpt_filter %>% distinct(n_assays) %>% pull(n_assays)}
                             QC PASS = fitc > 36 & flag.length <= 4"),
       fill = "compound QC") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90),
        panel.grid = element_blank())
qc_plot_bar
cowplot::ggsave2(qc_plot_bar, file = "plots/toxcast_qc_by_assay_id.png", width = 12, height = 6)

#==============================================================================#
# Summary of analysis - The ToxCast data has few high quality estimates
# of curve parameters for endpoints related to cyotoxicity. Among the 
# 156 cytotoxicity related assays, none have more than 6 compounds with high
# quality curve fits. This is not enough to perform regression.
#==============================================================================#

#==============================================================================#
# Testing deming regression with CIs
#==============================================================================#
set.seed(33)
n <- 15
x <- rnorm(n, 0, 1)
y <- 3.1 + 6 * x + rnorm(n, 0, 2)
x2 <- x + rnorm(n, 0, 2)

set.seed(33)
library(deming)
fit1 <- lm(y ~ x2)
fit2 <- deming(y ~ x2)
Tfit <- lm(y ~ x)
summary(fit1)
fit2
?deming()

# try with SimplyAgree package
library(SimplyAgree)
?dem_reg()
datSA <- tibble::tibble(x = x2, y = y)
fit3 <- dem_reg(x = "x",
                y = "y",
                data = datSA)
check(fit3)

# #==============================================================================#
# # DO NOT DELETE - If needed, the code to join the Toxcast data for orthreg 
# # function is below in the section "Get The other dataset to join".
# #==============================================================================#
# #------------------------------------------------------------------------------#
# # Attempt 1: Get ToxCast data that match the nematode toxicants
# #------------------------------------------------------------------------------#
# # Define CAS numbers of interest
# cas_numbers <- c("10108-64-2", "114-26-1", "115-09-3", "115-86-6", "116-06-3",
#                  "121-75-5", "16752-77-5", "175013-18-0", "1912-24-9", "2921-88-2",
#                  "4685-14-7", "5234-68-4", "63-25-2", "7646-85-7", "7718-54-9",
#                  "7761-88-8", "8018-01-7", "94-75-7")
# 
# ## Get comptox MW data to convert ToxCast to w/v unit
# mw_data <- read_excel("data/processed/toxcast_mw.xlsx") %>%
#   dplyr::mutate(INPUT = as.character(INPUT)) %>%
#   dplyr::distinct(.keep_all = T) %>%
#   dplyr::mutate(AVERAGE_MASS = ifelse(CASRN == "8018-01-7", 541.08, AVERAGE_MASS)) %>% # fix Manxozeb
#   dplyr::select(CASRN, AVERAGE_MASS)
# 
# ## Filter the dataset for those CAS numbers, join mw data, transform to uM  
# mc5_worm <- mc5 %>% 
#   dplyr::filter(casn %in% cas_numbers) %>% # filter to worm cas
#   dplyr::left_join(mw_data, by = c("casn" = "CASRN")) %>% # join MW data
#   mutate(AC50_mg_L = (ac50 / 1e6) * AVERAGE_MASS * 1000, # Convert µM to mg/L {Why not ac50 * MW? TAC - this is mg/L}
#          AC20_mg_L = (ac20 / 1e6) * AVERAGE_MASS * 1000,
#          AC10_mg_L = (ac10 / 1e6) * AVERAGE_MASS * 1000,
#          AC5_mg_L = (ac5 / 1e6) * AVERAGE_MASS * 1000,
#          BMD_mg_L = (bmd / 1e6) * AVERAGE_MASS * 1000,
#          ACC_mg_L = (acc / 1e6) * AVERAGE_MASS * 1000) %>%
#   dplyr::group_by(spid) %>%
#   dplyr::mutate(n_compounds = length(unique(chnm))) %>%
#   dplyr::ungroup()
# 
# # look at table of fitc for filtering quality data - wow, lots of poor quality fits. Ideally we would filter to 
# # just 37,38 or 41,42 and < 3 flags
# expss::cross_cases(mc5_worm, flag.length, fitc)
# 
# # filter using fitc and flags - used 3+ flags as cuttoff (taken from toxcast documentation as )
# mc5_proc <- mc5_worm %>%
#   dplyr::mutate(QC = case_when(fitc %in% c(37, 38, 41, 42) & (flag.length < 3 | is.na(flag.length)) ~ "PASS",
#                                fitc %in% c(37, 38, 41, 42) & flag.length >= 3 ~ "FLAG",
#                                !(fitc %in% c(37, 38, 41, 42)) & (flag.length < 3 | is.na(flag.length)) ~ "CURVE",
#                                !(fitc %in% c(37, 38, 41, 42)) & flag.length >= 3 ~ "FLAG+CURVE",
#                                TRUE ~ "wtf?"))
# 
# # QC breakdown
# expss::cross_cases(mc5_proc, QC)
# 
# # filtered
# mc5_final <- mc5_proc %>%
#   dplyr::filter(QC == "PASS") %>% # apply QC filter
#   dplyr::mutate(group = "CELL", # make a group for these data
#                 latin_name = "unknown") %>% 
#   dplyr::select(cas = casn, chem_name = chnm, group, latin_name,
#                 spid, aenm, ac50, ac20, fitc, QC) %>%
#   dplyr::group_by(chem_name) %>%
#   dplyr::mutate(n_aenm = length(unique(aenm))) %>%
#   dplyr::ungroup()
# 
# # get a table of unique endpoints per group
# table(mc5_final$chem_name)
# 
# #------------------------------------------------------------------------------#
# # Attempt 2: Try getting all the viability endpoints then apply filters
# #------------------------------------------------------------------------------#
# # use the assay annotations to filter to the cytotoxicity assays
# assays <- readxl::read_excel("data/raw/assay_annotations_invitrodb_v4_1_SEPT2023.xlsx", sheet = "annotations_combined")
# 
# # filter to cytotoxicity
# filtered_assays <- assays %>%
#   dplyr::filter(intended_target_family_sub == "cytotoxicity" & organism == "human") 
# 
# # join the full mc5 data, filter to our compounds, add quality filters to see quality issues
# mc5_v2 <- filtered_assays %>%
#   dplyr::left_join(., mc5) %>%
#   dplyr::filter(casn %in% cas_numbers) %>% # filter to worm cas
#   dplyr::left_join(mw_data, by = c("casn" = "CASRN")) %>% # join MW data
#   mutate(AC50_mg_L = (ac50 / 1e6) * AVERAGE_MASS * 1000, # Convert µM to mg/L {orig. ac50 * MW * 1000} yikes!
#           AC20_mg_L = (ac20 / 1e6) * AVERAGE_MASS * 1000,
#           AC10_mg_L = (ac10 / 1e6) * AVERAGE_MASS * 1000,
#           AC5_mg_L = (ac5 / 1e6) * AVERAGE_MASS * 1000,
#           BMD_mg_L = (bmd / 1e6) * AVERAGE_MASS * 1000,
#           ACC_mg_L = (acc / 1e6) * AVERAGE_MASS * 1000) %>%
#   dplyr::mutate(QC = case_when(fitc %in% c(37, 38, 41, 42) & (flag.length < 3 | is.na(flag.length)) ~ "PASS",
#                                fitc %in% c(37, 38, 41, 42) & flag.length >= 3 ~ "FLAG",
#                                !(fitc %in% c(37, 38, 41, 42)) & (flag.length < 3 | is.na(flag.length)) ~ "CURVE",
#                                !(fitc %in% c(37, 38, 41, 42)) & flag.length >= 3 ~ "FLAG+CURVE",
#                                TRUE ~ "wtf?")) %>%
#   dplyr::mutate(QC2 = case_when((flag.length < 4 | is.na(flag.length)) & hitc >= 0.5 & ac50 > 10^logc_min ~ "PASS",
#                                TRUE ~ "FAIL")) %>% # scott filters
#   dplyr::group_by(aeid) %>%
#   dplyr::mutate(total_n_compounds = length(unique(chnm))) %>%
#   dplyr::group_by(aeid, QC) %>%
#   dplyr::mutate(n_compounds_QC = length(unique(chnm))) %>%
#   dplyr::group_by(aeid, QC2) %>%
#   dplyr::mutate(n_compounds_QC2 = length(unique(chnm))) %>%
#   dplyr::group_by(aeid, casn) %>%
#   dplyr::mutate(n_spid_in_aeid_casn = length(unique(spid))) %>%
#   dplyr::ungroup() %>%
#   dplyr::mutate(n_assays = length(unique(aeid)))
# 
# # plot the effect of filtering
# tc_qc_plot <- ggplot(mc5_v2 %>% dplyr::distinct(aeid, QC, .keep_all = T)) +
#   aes(x = as.character(aeid), y = n_compounds_QC, fill = QC) +
#   geom_bar(stat = "identity") +
#   geom_hline(yintercept = 16, linetype = 2, linewidth = 0.25) +
#   labs(x = "", y = "Number of compounds",
#        subtitle = glue::glue("ToxCast QC1 n={mc5_v2 %>% distinct(n_assays) %>% pull(n_assays)} cytotoxicity assays\nvalues > 16 b/c replicate samples in assay")) +
#   theme_bw() +
#   theme(axis.text.x = element_text(angle = 90))
# tc_qc_plot
# 
# # try now with Scott G. filtering only QC2 - def better.
# tc_qc2_plot <- ggplot(mc5_v2 %>% dplyr::distinct(aeid, QC2, .keep_all = T)) +
#   aes(x = as.character(aeid), y = n_compounds_QC2, fill = QC2) +
#   geom_bar(stat = "identity") +
#   geom_hline(yintercept = 16, linetype = 2, linewidth = 0.25) +
#   labs(x = "", y = "Number of compounds",
#        subtitle = glue::glue("ToxCast QC2 n={mc5_v2 %>% distinct(n_assays) %>% pull(n_assays)} cytotoxicity assays\nvalues > 16 b/c replicate samples in assay")) +
#   theme_bw() +
#   #facet_wrap(~timepoint_hr) +
#   theme(axis.text.x = element_text(angle = 90))
# tc_qc2_plot
# 
# # lets look at the distribution of timepoints for the assays
# tc_time_p <- ggplot(mc5_v2 %>% dplyr::distinct(aeid, .keep_all = T)) +
#   aes(x = as.character(aeid), y = timepoint_hr) +
#   geom_bar(stat = "identity") +
#   labs(x = "", y = "Timepoint (hours)",
#        subtitle = "ToxCast assay timepoints") +
#   theme_bw() +
#   theme(axis.text.x = element_text(angle = 90))
# tc_time_p
# 
# # look at the celltypes as replicates?
# # look at agonist and antagonists to remove those assasy?
# 
# 
# # OK, filter by QC2 
# # then group by hours and celltype, then show the most relevant assays? UGh, so many weird endpopints
# mc5_v2_select <- mc5_v2 %>%
#   #dplyr::filter(QC2 == "PASS") %>% # leave out filter for now
#   dplyr::arrange(aeid, chnm, AC50_mg_L) %>%
#   dplyr::select(AC50_mg_L, AC20_mg_L, AC10_mg_L, AC10_mg_L, AC5_mg_L, BMD_mg_L, ACC_mg_L,
#                 flag.length, fitc, hitc, ac50, logc_min, QC, QC2, chnm, timepoint_hr, spid, aeid, assay_component_endpoint_name, tissue, cell_format, cell_short_name,
#                 assay_source_desc, assay_component_desc, assay_component_target_desc, assay_component_endpoint_desc)
# 
# # look at the distribution of AC50s for passing assays for each compound - focus on 24 hours
# tp_plot <- ggplot(mc5_v2_select %>% dplyr::filter(QC2 == "PASS")) +
#   aes(x = chnm, y = AC50_mg_L, fill = fitc) +
#   geom_jitter(shape = 21) +
#   geom_boxplot(outlier.shape = NA, alpha = 0.5) +
#   labs(title = "cytotoxicity QC2 pass AC50 distributions by timepoint and drug") +
#   facet_wrap(~timepoint_hr) + 
#   theme_bw() +
#   theme(axis.text.x = element_text(angle = 90))
# 
# # filter to 24 hour timepoint and manually review assays to choose best dataset
# mc5_v2_final_1 <- mc5_v2_select %>%
#   dplyr::filter(QC2 == "PASS" & timepoint_hr == 24)
# assays_to_review <- mc5_v2_final_1 %>%
#   dplyr::distinct(aeid, .keep_all = T) %>%
#   dplyr::select(17:25)
# rio::export(assays_to_review, file = "data/processed/assays_to_review.csv")
# 
# # try filtering to these manually currated assays and see how many drugs are present
# man1 <- mc5_v2_select %>%
#   dplyr::filter(assay_component_endpoint_name %in% c("TOX21_RT_HEPG2_FLO_24hr_viability",
#                              "TOX21_RT_HEK293_GLO_24hr_viability",
#                              "CCTE_Simmons_CellTiterGLO_HEK293T",
#                              "BSK_SAg_PBMCCytotoxicity"))
# 
# man_plot1 <- ggplot(man1) +
#   aes(x = chnm, y = ac50, fill = as.character(fitc)) +
#   geom_jitter(shape = 21) +
#   geom_boxplot(outlier.shape = NA, alpha = 0.5) +
#   facet_wrap(~assay_component_endpoint_name) + 
#   theme_bw() +
#   theme(axis.text.x = element_text(angle = 90))
# 
# man_plot2 <- ggplot(man1) +
#   aes(x = chnm, y = AC50_mg_L, fill = as.character(fitc)) +
#   geom_jitter(shape = 21) +
#   geom_boxplot(outlier.shape = NA, alpha = 0.5) +
#   facet_wrap(~assay_component_endpoint_name) + 
#   theme_bw() +
#   theme(axis.text.x = element_text(angle = 90))
# 
# ggsave(tc_qc_plot, filename = "plots/toxcast_QC_cytotoxicity.png", width = 7, height = 5)
# ggsave(tc_qc2_plot, filename = "plots/toxcast_QC2_cytotoxicity.png", width = 7, height = 5)
# ggsave(tc_time_p, filename = "plots/toxcast_timepoints_cytotoxicity.png", width = 7, height = 5)
# ggsave(man_plot1, filename = "plots/toxcast_manual_cur1.png", width = 7, height = 7)
# ggsave(man_plot2, filename = "plots/toxcast_manual_cur2.png", width = 7, height = 7)
# #------------------------------------------------------------------------------#
# # Get The other dataset to join
# #------------------------------------------------------------------------------#
# # load cleaned data and drop group for now - causes eror with pwOrthReg function when assigning latin_name to group
# dat <- data.table::fread("data/processed/03_clean.csv")
# 
# # use data.dic to guide the join
# glimpse(dat)
# 
# # Setup the TOXCAST data to work with the ploting functions
# mc5_join <- mc5_v2 %>%
#   dplyr::mutate(casn = as.integer(stringr::str_replace_all(casn, pattern = "-", replacement = "")),
#                 group = paste("TOXCAST", organism, sep = "_"),
#                 latin_name = NA_character_,
#                 effect_unit = "mg/L",
#                 source = "TOXCAST",
#                 chem_name = case_when(chnm == "2,4-Dichlorophenoxyacetic acid" ~ "2,4-D",
#                                       chnm == "Methylmercuric(II) chloride" ~ "Methylmercury chloride",
#                                       TRUE ~ chnm),
#                 test_type = NA,
#                 duration_d = timepoint_hr / 24,
#                 endpoint = assay_component_endpoint_name) %>%
#   dplyr::select(cas = casn, chem_name, group, latin_name, strain = spid, test_type,
#                 duration_d, effect_unit, source, endpoint,
#                 AC50_mg_L, AC20_mg_L, AC10_mg_L, AC10_mg_L, AC5_mg_L, BMD_mg_L, ACC_mg_L,
#                 hitc, fitc, QC, QC2, flag.length) %>%
#   tidyr::pivot_longer(cols = AC50_mg_L:ACC_mg_L,
#                       names_to = "test_statistic", values_to = "effect_value") %>%
#   dplyr::select(cas, chem_name, group, latin_name, strain, test_type,
#                 test_statistic, duration_d, effect_value, effect_unit, source, endpoint, QC2)
# 
# # Bind them together
# bind <- dplyr::bind_rows(dat, mc5_join %>% dplyr::select(-QC2))
# 
# # look at the distributions of AC50s
# effect_hist <- ggplot(bind) +
#   aes(x = effect_value) +
#   geom_histogram() +
#   facet_wrap(~group, scales = "free")
# effect_hist
# 
# # ugh, these values are often completely insane
# effect_hist2 <- ggplot(mc5_join) +
#   aes(x = effect_value, fill = QC2) +
#   geom_histogram() +
#   facet_wrap(~endpoint, scales = "free")
# effect_hist2
# #===========================================================#
# # Figure 3 - C. elegans toxicity data compared with rats,
# #            fish, invertebrates, and algae.
# #===========================================================#
# # 3a
# # Widmayer vs RAT
# test <- bind %>%
#   dplyr::filter(group == "NEMATODE" | group == "TOXCAST_human")
# 
# # run the orth regressions across all pairs
# or_test <- pwOrthReg(data = test, group = "group", limit.comp = "NEMATODE", min.n = 5, message = T, plot = T)
# or_df <- data.table::rbindlist(or_test$orthregs)
# or_test$plots[341] # random test
# or_test$plots[587]# random test
# or_test$plots[75]# random test
# 
# # lets just look at these select assay endpoints
# test_small <- bind %>%
#   dplyr::filter(group == "NEMATODE" | group == "TOXCAST_human") %>%
#   dplyr::filter(endpoint %in% c("Growth", "TOX21_RT_HEPG2_FLO_24hr_viability",
#                                 "TOX21_RT_HEK293_GLO_24hr_viability", "CCTE_Simmons_CellTiterGLO_HEK293T",
#                                 "BSK_SAg_PBMCCytotoxicity"))
# 
# # run the orth regressions across all pairs - ugh, all ugly. 
# or_test_small <- pwOrthReg(data = test_small, group = "group", limit.comp = "NEMATODE", min.n = 5, message = T, plot = T)
# or_df_small <- data.table::rbindlist(or_test_small$orthregs)
# or_test_small$plots
# 

# # Filter ToxCast for criteria based on Friedman 2020 and updated with new release note
# mc5_filtered <- mc5 %>%
#   dplyr::filter(flag.length < 4 & hitc >= 0.5 & ac50 > 10^logc_min) # Scott, do we want these OR my filters or below?
# 
# # The fifth percentile of that distribution was used to identify a minimum bioactive concentration for each chemical in ToxCast,
# # regardless of the specific biological pathways involved.
# ac50_5th_percentile <- mc5 %>%
#   group_by(casn) %>%
#   summarise(fifth_percentile_ac50 = quantile(ac50, probs = 0.05, na.rm = TRUE))

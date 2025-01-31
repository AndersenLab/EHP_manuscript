#!/usr/bin/env Rscript
library(tidyverse)

# set working dir to base directory of repository
setwd(paste0(dirname(rstudioapi::getActiveDocumentContext()$path),"/.."))

# load functions
source("code/functions.R")

# load cleaned data and drop group for now - causes eror with pwOrthReg function when assigning latin_name to group
dat <- data.table::fread("data/processed/00_data.csv")

# filter to zebrafish studies
zf <- dat %>%
  dplyr::filter(latin_name == "Danio rerio")

# find groups - 4 groups, note Scholz_ZF has Padilla mortality data in it.
unique(zf$group)
table(zf$group)

#------------------------------------------------------------------------------#
# Part 1: look at all Padilla data
#------------------------------------------------------------------------------#
TC_padilla <- dat %>%
  dplyr::filter(group == "NEMATODE" | group == "TC_ZF_Padilla" | group == "NEMATODE_COPAS1") 

# run the orth regressions across all pairs to NEMATODE with QC filter
or1 <- pwOrthReg(data = TC_padilla, group = "group",  limit.comp = "NEMATODE", min.n = 5, message = T, QC = "ignore", plot = T)
or1df <- data.table::rbindlist(or1$orthregs)
# highest r2
print(glue::glue("slope = {round(or1df[1]$orth.reg.slope, digits = 2)}, y-intercept = {round(or1df[1]$orth.reg.intercept, digits = 2)}, r^2 = {round(or1df[1]$orth.reg.r.squared, digits = 2)}"))
Padilla <- cowplot::plot_grid(plotlist = or1$plots, ncol = 4)
or1$plots[[1]]

cowplot::ggsave2(Padilla, filename = "plots/TC_padilla.png", width = 15, height = 15)

#------------------------------------------------------------------------------#
# Part 1b: look at all high quality Padilla data
#------------------------------------------------------------------------------#
# run the orth regressions across all pairs to NEMATODE with QC filter
or1_qc <- pwOrthReg(data = TC_padilla, group = "group",  limit.comp = "NEMATODE", min.n = 5, message = T, QC = "filter", plot = T)
or1df_qc <- data.table::rbindlist(or1_qc$orthregs)
# highest r2
print(glue::glue("slope = {round(or1df_qc[1]$orth.reg.slope, digits = 2)}, y-intercept = {round(or1df_qc[1]$orth.reg.intercept, digits = 2)}, r^2 = {round(or1df_qc[1]$orth.reg.r.squared, digits = 2)}"))
Padilla_qc <- cowplot::plot_grid(plotlist = or1_qc$plots, ncol = 4)
or1_qc$plots[[1]]

cowplot::ggsave2(Padilla_qc, filename = "plots/TC_padilla_qc_filter.png", width = 15, height = 15)

# look at original data for Padilla
orig_padilla <- data.table::fread("data/processed/03_clean.csv") %>%
  dplyr::filter(grepl(source, pattern = "Padilla")) %>%
  dplyr::arrange(chem_name)

new_padilla <- TC_padilla %>%
  dplyr::filter(group == "TC_ZF_Padilla" & test_statistic == "AC50") %>%
  dplyr::arrange(chem_name) %>%
  dplyr::filter(chem_name %in% orig_padilla$chem_name) %>%
  dplyr::select(chem_name, cas, effect_value_new = effect_value)

join_padilla <- left_join(orig_padilla, new_padilla)

comp <- ggplot(join_padilla) +
  aes(x = effect_value_new, y = effect_value) +
  geom_abline(slope = 1, intercept = 0, linetype = 2, color = "black") +
  geom_point() +
  theme_bw() +
  ggrepel::geom_label_repel(aes(label = chem_name),
                  box.padding   = 0.35, 
                  point.padding = 0.5,
                  segment.color = 'grey50') +
  labs(x = "TC_ZF_Padilla new AC50 (mg/L)", y = "ZF_EMBRYO1 original AC50 (mg/L)")

cowplot::ggsave2(comp, filename = glue::glue("plots/{today}_TC_padilla_new_vs_ZF_EMBRYO1_orig.png"), width = 6, height = 6)
rio::export(join_padilla, file = glue::glue("data/processed/{today}_TC_padilla_new_vs_ZF_EMBRYO1_orig.csv"))
#------------------------------------------------------------------------------#
# Part 2: look at all Tanguay data
#------------------------------------------------------------------------------#
TC_tang <- dat %>%
  dplyr::filter(group == "NEMATODE" | group == "TC_ZF_Tanguay" | group == "NEMATODE_COPAS1") 

# run the orth regressions across all pairs to NEMATODE with QC filter
or2 <- pwOrthReg(data = TC_tang, group = "group",  limit.comp = "NEMATODE", min.n = 5, message = T, QC = "ignore", plot = T)
or2df <- data.table::rbindlist(or2$orthregs)
# highest r2
print(glue::glue("slope = {round(or2df[1]$orth.reg.slope, digits = 2)}, y-intercept = {round(or2df[1]$orth.reg.intercept, digits = 2)}, r^2 = {round(or2df[1]$orth.reg.r.squared, digits = 2)}"))
tang <- cowplot::plot_grid(plotlist = or2$plots, ncol = 4)
cowplot::ggsave2(tang, filename = "plots/Tang.png", width = 15, height = 15)

#------------------------------------------------------------------------------#
# Part 3: look at all Scholz data
#------------------------------------------------------------------------------#
scholz <- dat %>%
  dplyr::filter(group == "NEMATODE" | group == "SCHOLZ_ZF" | group == "NEMATODE_COPAS1") 

# run the orth regressions across all pairs to NEMATODE with QC filter
or3 <- pwOrthReg(data = scholz, group = "group",  limit.comp = "NEMATODE", min.n = 5, message = T, QC = "ignore", plot = T)
or3df <- data.table::rbindlist(or3$orthregs)
# highest r2
print(glue::glue("slope = {round(or3df[1]$orth.reg.slope, digits = 2)}, y-intercept = {round(or3df[1]$orth.reg.intercept, digits = 2)}, r^2 = {round(or3df[1]$orth.reg.r.squared, digits = 2)}"))
scholz_plots <- cowplot::plot_grid(plotlist = or3$plots, ncol = 2)
cowplot::ggsave2(scholz_plots, filename = "plots/scholz.png", width = 8, height = 8)
# not a lot of data overlap

#------------------------------------------------------------------------------#
# Part 4: Envirotox fish data alone
#------------------------------------------------------------------------------#
etox <- dat %>%
  dplyr::filter(group == "NEMATODE" | group == "FISH" | group == "NEMATODE_COPAS1") %>%
  dplyr::filter(case_when(group %in% c("NEMATODE",  "NEMATODE_COPAS1") ~ T,
                          group == "FISH" & latin_name == "Danio rerio" ~ T,
                          TRUE ~ FALSE)) %>%
  #dplyr::filter(case_when(group == "NEMATODE" ~ T,
  #                        group == "NEMATODE_COPAS1" ~ T,
  #                        source == "EnviroTox DB" & duration_d == 4 & test_statistic %in% c("LC50") ~ T, # removing EC50s has a big filtering effect, even though some are coded as mortality
  #                        TRUE ~ F)) %>%
  dplyr::mutate(endpoint2 = dplyr::case_when(endpoint %in% c("Mortality/Growth",
                                                             "Mortality, Mortality",
                                                             "Mortality, Survival") ~ "Mortality",
                                             endpoint == "Immobilization: Change in the failure to respond or lack of movement after mechanical stimulation." ~ "Immobility",
                                             endpoint == "Intoxication, Immobile" ~ "Immobility",
                                             TRUE ~ endpoint)) %>%
  dplyr::mutate(endpoint = endpoint2) %>%
  dplyr::select(-endpoint2)

# run the orth regressions across all pairs to NEMATODE with QC filter
or4 <- pwOrthReg(data = etox, group = "group",  limit.comp = "NEMATODE", min.n = 5, message = T, QC = "ignore", plot = T)
or4df <- data.table::rbindlist(or4$orthregs)
# highest r2
print(glue::glue("slope = {round(or4df[1]$orth.reg.slope, digits = 2)},
                 y-intercept = {round(or4df[1]$orth.reg.intercept, digits = 2)},
                 r^2 = {round(or4df[1]$orth.reg.r.squared, digits = 2)}"))
FISH_plots <- cowplot::plot_grid(plotlist = or4$plots, ncol = 4)
cowplot::ggsave2(FISH_plots, filename = "plots/zf_envirotox.png", width = 15, height = 15)
# not a lot of data overlap
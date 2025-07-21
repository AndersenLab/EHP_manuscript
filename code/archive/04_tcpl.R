#!/usr/bin/env Rscript
#devtools::install_github("USEPA/CompTox-ToxCast-tcpl")
library(tidyverse)
library(data.table)
library(plotly)
library(tcplfit2)
library(tcpl)

# set working dir to base directory of repository
setwd(paste0(dirname(rstudioapi::getActiveDocumentContext()$path),"/.."))

# load functions
source("code/functions.R")

# configure tcpl
tcplConfPrompt(db   = "invitrodb_v2",
               user = "root",
               host = "localhost",
               drvr = "MySQL")

# load cleaned data and drop group for now - causes eror with pwOrthReg function when assigning latin_name to group
dat <- data.table::fread("data/processed/03_clean.csv") %>%
  dplyr::mutate(cas2 = paste0(stringr::str_sub(cas, end = -4), "-",
                             stringr::str_sub(cas, start = -3, end = -2), "-",
                             stringr::str_sub(cas, start = -1))) # add proper cas id with hyphens

#===========================================================================####
# Part 1: Review database
#===========================================================================####
# look across Assay Registration Elements - Hierarchical structure.
# Assay source, assay, assay component, and assay endpoint are registered via tcpl scripting into a collection of tables.
# These assay element tables broadly describe who conducted the assay, what platform was used, what was being measured (raw readout), and how the measurement was interpreted (normalized component data).
# A hierarchical structure of the assay elements is as follows: assay source > assay > assay component > assay component endpoint.
# From a database v3.4 snapshot taken in December 2020, InvitroDB supported 32 assay sources, 763 assays, 2074 components, and 2780 endpoints.
# asid (assay source ID) - only 26 on 20230515, must be older database, ugh.
asids <- tcplLoadAsid() 
# aid (assay ID)
aids <- tcplLoadAid()
# acid (assay component ID)
acids <- tcplLoadAcid()
# aeid (assay endpoint ID)
aeids <- tcplLoadAeid()

#===========================================================================####
# Part 2: Pull all multi concentration level 5 data for our toxicants
#===========================================================================####
# pull all mc5 data
mc5 <- tcplPrepOtpt(
  tcplLoadData(lvl=5, # data level
               fld="aeid", # fields to query on
               val=aeids$aeid, # value for each field
               # values should match their corresponding 'fld'
               type = "mc") # return additional parameters from mc5_param 
)

# filter to our 22 chemicals - no silver nitrate
mc5.chem <- mc5 %>%
  dplyr::filter(casn %in% dat$cas2)

# what's missing? Just silver nitrate, the cas2 is correct though, 7761-88-8
missing <- tibble::tibble(missing = (unique(dat$cas2) %in% unique(mc5.chem$casn)),
                          cas2 = unique(dat$cas2),
                          chem_name = unique(dat$chem_name))

# export the rda file
#save(mc5.chem, file = "data/processed/mc5.chem.issue.rda")
#rio::export(as.data.frame(missing$cas2), file = "data/processed/chems.csv")
#rio::export(tibble::tibble(vars = names(mc5.chem)), file = "data/processed/vars.csv")
#===========================================================================####
# Part 3: Filter to rows with a model fit and endpoints we care about
#===========================================================================####
# with fits
mc5.chem.fit <- mc5.chem %>%
  dplyr::filter(nconc >= 5 & !is.na(modl_ga))

# with vialbility, cytotoxicity, mortality
mc5.chem.fit.end <- mc5.chem.fit %>%
  dplyr::filter(grepl(aenm, pattern = "viability|cytotoxicity|MORT"))



# export these
save(mc5.chem.fit.end, file = "data/processed/mc5.chem.fit.end.rda")


############################################################################ OLD
# # Provide the chemical name and assign to 'chnm'.
# chnm <- 'Chlorothalonil'
# # Load the chemical data from the database.
# chem <- tcplLoadChem(field = 'chnm',val = chnm)
# # Load mc5 data from the database for the specified chemical.
# Chlorothalonil.mc5 <- tcplLoadData(lvl = 5, # data level 
#                                    fld = 'spid', # field to query on
#                                    val = chem[,spid], # value for each field (fld)
#                                    type = 'mc',
#                                    add.fld=TRUE) # data type - MC
# 
# test <- Chlorothalonil.mc5 %>%
#   dplyr::filter(aeid == 5)


# # Create table of all assay endpoint ids (aeids) per assay source
# aeids <- tcplLoadAeid(fld="asid", # field to query on
#                       val=14, # value for each field
#                       # values should match their corresponding 'fld'
#                       add.fld = c("aid", "anm", "acid", "acnm")) # additional fields to return
# 
# # load multi concentration data level 5
# mc5 <- tcplPrepOtpt(
#   tcplLoadData(lvl=5, # data level
#                fld="aeid", # fields to query on
#                val=aeids$aeid, # value for each field
#                # values should match their corresponding 'fld'
#                type = "mc") # data type - MC
# )
# 
# # testing this chunk for methods
# # Create a function to list all available methods function (SC & MC).
# method_list <- function() {
#   # Single Concentration
#   ## Level 1
#   sc1 <- tcplMthdList(1, 'sc')
#   sc1[, lvl := "sc1"]
#   setnames(sc1, c("sc1_mthd", "sc1_mthd_id"), c("mthd", "mthd_id"))
#   ## Level 2
#   sc2 <- tcplMthdList(2, 'sc')
#   sc2[, lvl := "sc2"]
#   setnames(sc2, c("sc2_mthd", "sc2_mthd_id"), c("mthd", "mthd_id"))
#   
#   # Multiple Concentration
#   ## Level 2
#   mc2 <- tcplMthdList(2, 'mc')
#   mc2[, lvl := "mc2"]
#   setnames(mc2, c("mc2_mthd", "mc2_mthd_id"), c("mthd", "mthd_id"))
#   ## Level 3
#   mc3 <- tcplMthdList(3, 'mc')
#   mc3[, lvl := "mc3"]
#   setnames(mc3, c("mc3_mthd", "mc3_mthd_id"), c("mthd", "mthd_id"))
#   ## Level 4
#   mc4 <- tcplMthdList(4, 'mc')
#   mc4[, lvl := "mc4"]
#   setnames(mc4, c("mc4_mthd", "mc4_mthd_id"), c("mthd", "mthd_id"))
#   ## Level 5
#   mc5 <- tcplMthdList(5, 'mc')
#   mc5[, lvl := "mc5"]
#   setnames(mc5, c("mc5_mthd", "mc5_mthd_id"), c("mthd", "mthd_id"))
#   # Compile the Output
#   mthd.list <- rbind(sc1, sc2, mc2, mc3, mc4, mc5)
#   mthd.list <- mthd.list[, c("lvl", "mthd_id", "mthd", "desc")]
#   # Return the Results
#   return(mthd.list)
# }
# 
# # Run the 'method_list' functions and store output.
# amthds <- method_list()
# # Print the available methods list.
# amthds

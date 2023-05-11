#!/usr/bin/env Rscript
library(tidyverse)
library(data.table)
library(plotly)
library(tcplfit2)
library(tcpl)

# set working dir to base directory of repository
setwd(paste0(dirname(rstudioapi::getActiveDocumentContext()$path),"/.."))

# Testing tcpl
tcplConfList()
?tcplConf()

# setup Assay Elements
# List all assay source IDs
tcplLoadAsid() 
# Create table of all assay endpoint ids (aeids) per assay source
aeids <- tcplLoadAeid(fld="asid", # field to query on
                      val=14, # value for each field
                      # values should match their corresponding 'fld'
                      add.fld = c("aid", "anm", "acid", "acnm")) # additional fields to return

# load multi concentration data level 5
mc5 <- tcplPrepOtpt(
  tcplLoadData(lvl=5, # data level
               fld="aeid", # fields to query on
               val=aeids$aeid, # value for each field
               # values should match their corresponding 'fld'
               type = "mc") # data type - MC
)

# testing this chunk for methods
# Create a function to list all available methods function (SC & MC).
method_list <- function() {
  # Single Concentration
  ## Level 1
  sc1 <- tcplMthdList(1, 'sc')
  sc1[, lvl := "sc1"]
  setnames(sc1, c("sc1_mthd", "sc1_mthd_id"), c("mthd", "mthd_id"))
  ## Level 2
  sc2 <- tcplMthdList(2, 'sc')
  sc2[, lvl := "sc2"]
  setnames(sc2, c("sc2_mthd", "sc2_mthd_id"), c("mthd", "mthd_id"))
  
  # Multiple Concentration
  ## Level 2
  mc2 <- tcplMthdList(2, 'mc')
  mc2[, lvl := "mc2"]
  setnames(mc2, c("mc2_mthd", "mc2_mthd_id"), c("mthd", "mthd_id"))
  ## Level 3
  mc3 <- tcplMthdList(3, 'mc')
  mc3[, lvl := "mc3"]
  setnames(mc3, c("mc3_mthd", "mc3_mthd_id"), c("mthd", "mthd_id"))
  ## Level 4
  mc4 <- tcplMthdList(4, 'mc')
  mc4[, lvl := "mc4"]
  setnames(mc4, c("mc4_mthd", "mc4_mthd_id"), c("mthd", "mthd_id"))
  ## Level 5
  mc5 <- tcplMthdList(5, 'mc')
  mc5[, lvl := "mc5"]
  setnames(mc5, c("mc5_mthd", "mc5_mthd_id"), c("mthd", "mthd_id"))
  # Compile the Output
  mthd.list <- rbind(sc1, sc2, mc2, mc3, mc4, mc5)
  mthd.list <- mthd.list[, c("lvl", "mthd_id", "mthd", "desc")]
  # Return the Results
  return(mthd.list)
}

# Run the 'method_list' functions and store output.
amthds <- method_list()
# Print the available methods list.
amthds

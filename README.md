# EHP_manuscript
This repository holds data, code, and plots for the Environmental  Health Perspectives manuscript written with Dr. Scott Glaberman

## Usage
* Open the `02_orthogonal_regression.R` script in Rstudio.
* Run the code to see how the the `orthReg` function handles different comparisons.
    * The `orthReg` function takes three parameters, `data`, `x`, `y`
    * The `data` parameter should be set to the cleaned data `data/processed/01_clean.csv`, which was sourced from the master sheet in `data/raw/Master_File.xlsx`.
    * The `x` parameter is a vector with taxa, test_statistic, and duration_h for the taxa you want to be on the x-axis.
    * The `y` parameter is another vector of the same form as `x`, but it corresponds to what will fall on the y-axis.

## Examples
```
# set working dir to base directory of repository
setwd(paste0(dirname(rstudioapi::getActiveDocumentContext()$path),"/.."))

# load functions
source("code/functions.R")

# load cleaned data
dat <- data.table::fread("data/processed/01_clean.csv")

# run the orthReg function on Rainbow trout and D. magna with specific test_statistics and durations. 
test1 <- orthReg(data = dat,
                x = c("Oncorhynchus mykiss", "LC50", "96"),
                y = c("Daphnia magna", "EC50", "48"))
```
library(dplyr)
library(readxl)


## Loads ToxCast data
load("/Users/scott/Library/CloudStorage/Dropbox/Work Sync/Current Analysis/Pesticide ERA Analysis/ToxCast Resources/INVITRODB_V4_1_SUMMARY/mc5-6_winning_model_fits-flags_invitrodb_v4_1_SEPT2023.Rdata")
## OR
load("mc5-6_winning_model_fits-flags_invitrodb_v4_1_SEPT2023.Rdata")

## Filter ToxCast for criteria based on Friedman 2020 and updated with new release note
mc5_filtered <- mc5[mc5$flag.length < 4 & mc5$hitc >= 0.5 & mc5$ac50 > 10^mc5$logc_min, ]
head(mc5_filtered)


## The fifth percentile of that distribution was used to identify a minimum bioactive concentration for each chemical in ToxCast, regardless of the specific biological pathways involved.
ac50_5th_percentile <- mc5 %>%
  group_by(casn) %>%
  summarise(fifth_percentile_ac50 = quantile(ac50, probs = 0.05, na.rm = TRUE))

###############################################################################################################

## Now I want to work with the nematode data

## Define CAS numbers of interest
cas_numbers <- c("10108-64-2", "114-26-1", "115-09-3", "115-86-6", "116-06-3",
                 "121-75-5", "16752-77-5", "175013-18-0", "1912-24-9", "2921-88-2",
                 "4685-14-7", "5234-68-4", "63-25-2", "7646-85-7", "7718-54-9",
                 "7761-88-8", "8018-01-7", "94-75-7")

## Filter the dataset for those CAS numbers
mc5_filtered_worm <- mc5[mc5$casn %in% cas_numbers, ]


## Convert ToxCast to w/v unit
## Load in Comptox MW data
setwd("/Users/scott/Library/CloudStorage/Dropbox/Work Sync/Current Analysis/Nematode Analysis/new/toxcast comparison")
mw_data <- read_excel("toxcast_mw.xlsx")
## Remove dashes from ToxCast data
mc5_filtered_worm$casn <- gsub("-", "", mc5_filtered_worm$casn)
mw_data$INPUT <- as.character(mw_data$INPUT)  # Ensure mw_data CAS numbers are strings
# Perform the merge allowing Cartesian join
merged_data_worm <- merge(mc5_filtered_worm, mw_data, by.x = "casn", by.y = "INPUT", all.x = TRUE, allow.cartesian = TRUE)

# Calculate AC50 and ACC in µg/L using molecular weight
merged_data_worm <- merged_data_worm %>%
  mutate(AC50_ug_L = ac50 * AVERAGE_MASS * 1000, # Convert µM to µg/L
         ACC_ug_L = acc * AVERAGE_MASS * 1000)  # Convert µM to µg/L












## Convert ToxCast to w/v unit
## Load in Comptox MW data
setwd("/Users/scott/Library/CloudStorage/Dropbox/Work Sync/Current Analysis/Nematode Analysis/new/toxcast comparison")
mw_data <- read_excel("toxcast_mw.xlsx")

# Define CAS numbers of interest for nematode data
cas_numbers <- c("10108-64-2", "114-26-1", "115-09-3", "115-86-6", "116-06-3",
                 "121-75-5", "16752-77-5", "175013-18-0", "1912-24-9", "2921-88-2",
                 "4685-14-7", "5234-68-4", "63-25-2", "7646-85-7", "7718-54-9",
                 "7761-88-8", "8018-01-7", "94-75-7")

# Filter ToxCast data for nematode CAS numbers and criteria based on Friedman 2022
mc5_filtered_worm <- mc5 %>%
  filter(casn %in% cas_numbers, flag.length < 3, hitc >= 0.5, ac50 > 10^logc_min)

# Load molecular weight data
mw_data <- read_excel("path_to_your_toxcast_mw.xlsx")

# Prepare data for merging
mc5_filtered_worm$casn <- gsub("-", "", mc5_filtered_worm$casn) # Remove dashes
mw_data$INPUT <- as.character(mw_data$INPUT) # Ensure CAS numbers are strings


# Calculate AC50 and ACC in µg/L using molecular weight for each CAS number in mc5_filtered_worm
mc5_filtered_worm$AC50_ug_L <- NA  # Initialize the column for AC50 in µg/L
mc5_filtered_worm$ACC_ug_L <- NA  # Initialize the column for ACC in µg/L

# Loop through each row in mc5_filtered_worm to calculate the new values
for(i in 1:nrow(mc5_filtered_worm)){
  # Get the molecular weight for the current CAS number
  mw <- mw_data$AVERAGE_MASS[mw_data$INPUT == mc5_filtered_worm$casn[i]]
  if(length(mw) > 0){
    # Calculate AC50 and ACC in µg/L
    mc5_filtered_worm$AC50_ug_L[i] <- mc5_filtered_worm$ac50[i] * mw * 1000
    mc5_filtered_worm$ACC_ug_L[i] <- mc5_filtered_worm$acc[i] * mw * 1000
  }
}

# Now mc5_filtered_worm includes the AC50 and ACC values in µg/L
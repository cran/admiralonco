# Name: ADTTE
#
# Label: Time to Event Analysis Dataset
#
# Input: adsl, adrs, tte_source objects
library(admiral)
library(admiralonco)
library(admiral.test) # Contains example datasets from the CDISC pilot project
library(dplyr)
library(lubridate)

# ---- Load source datasets ----

# Use e.g. haven::read_sas to read in .sas7bdat, or other suitable functions
# as needed and assign to the variables below.
# For illustration purposes read in admiral test data

data("admiral_adsl")
data("admiral_adrs")

adsl <- admiral_adsl
adrs <- admiral_adrs

# ---- Derivations ----

# Add response date to ADSL / ADRS for duration of response calculation
adsl <- adsl %>%
  derive_vars_merged(
    dataset_add = adrs,
    filter_add = PARAMCD == "RSP" & AVALC == "Y" & ANL01FL == "Y",
    by_vars = vars(STUDYID, USUBJID),
    new_vars = vars(TEMP_RESPDT = ADT)
  )

adrs <- adrs %>%
  derive_vars_merged(
    dataset_add = adsl,
    by_vars = vars(STUDYID, USUBJID),
    new_vars = vars(TEMP_RESPDT)
  )

# Use pre-defined tte_source objects to derive Overall Survival, Progression
# Free Survival and Duration of Response
adtte <- derive_param_tte(
  dataset_adsl = adsl,
  start_date = RANDDT,
  event_conditions = list(death_event),
  censor_conditions = list(lastalive_censor, rand_censor),
  source_datasets = list(adsl = adsl, adrs = adrs),
  set_values_to = vars(PARAMCD = "OS", PARAM = "Overall Survival")
) %>%
  derive_param_tte(
    dataset_adsl = adsl,
    start_date = RANDDT,
    event_conditions = list(pd_event, death_event),
    censor_conditions = list(lasta_censor, rand_censor),
    source_datasets = list(adsl = adsl, adrs = adrs),
    set_values_to = vars(PARAMCD = "PFS", PARAM = "Progression Free Survival")
  ) %>%
  derive_param_tte(
    dataset_adsl = adsl,
    start_date = TEMP_RESPDT,
    event_conditions = list(pd_event, death_event),
    censor_conditions = list(lasta_censor),
    source_datasets = list(
      adsl = filter(adsl, !is.na(TEMP_RESPDT)),
      adrs = filter(adrs, !is.na(TEMP_RESPDT))
    ),
    set_values_to = vars(PARAMCD = "RSD", PARAM = "Duration of Response")
  )

# Derive analysis value and sequence
adtte <- adtte %>%
  derive_vars_duration(
    new_var = AVAL,
    start_date = STARTDT,
    end_date = ADT
  ) %>%
  derive_var_obs_number(
    by_vars = vars(STUDYID, USUBJID),
    order = vars(PARAMCD),
    check_type = "error"
  )

# Join any required ADSL variables
adtte <- adtte %>%
  derive_vars_merged(
    dataset_add = adsl,
    new_vars = vars(ARMCD, ARM, ACTARMCD, ACTARM, AGE, SEX),
    by_vars = vars(STUDYID, USUBJID)
  )

# ---- Save output ----

dir <- tempdir() # Change to whichever directory you want to save the dataset in
saveRDS(adtte, file = file.path(dir, "adtte.rds"), compress = "bzip2")

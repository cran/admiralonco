## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)
link <- function(text, url) {
  return(
    paste0(
      "[", text, "]",
      "(", url, ")"
    )
  )
}
dyn_link <- function(text,
                     base_url,
                     relative_url = "",
                     # Change to TRUE when admiral adopts multiversion docs
                     is_multiversion = FALSE,
                     multiversion_default_ref = "main") {
  url <- paste(base_url, relative_url, sep = "/")
  if (is_multiversion) {
    url <- paste(
      base_url,
      Sys.getenv("BRANCH_NAME", multiversion_default_ref),
      relative_url,
      sep = "/"
    )
  }
  return(link(text, url))
}
# Other variables
admiral_homepage <- "https://pharmaverse.github.io/admiral"

library(admiraldev)

## ---- warning=FALSE, message=FALSE--------------------------------------------
library(admiral)
library(admiralonco)
library(dplyr)
library(admiral.test)
library(lubridate)

## -----------------------------------------------------------------------------
data("admiral_adsl")
data("admiral_adrs")
adsl <- admiral_adsl
adrs <- admiral_adrs

## ----echo=FALSE---------------------------------------------------------------
knitr::kable(list_tte_source_objects(package = "admiralonco"))

## ---- eval=FALSE--------------------------------------------------------------
#  adsl_death_event <- event_source(
#    dataset_name = "adsl",
#    date = DTHDT,
#    set_values_to = vars(
#      EVNTDESC = "STUDY DEATH",
#      SRCDOM = "ADSL",
#      SRCVAR = "DTHDT"
#    )
#  )

## -----------------------------------------------------------------------------
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

## -----------------------------------------------------------------------------
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
    source_datasets = list(adsl = filter(adsl, !is.na(TEMP_RESPDT)), adrs = filter(adrs, !is.na(TEMP_RESPDT))),
    set_values_to = vars(PARAMCD = "RSD", PARAM = "Duration of Response")
  )

## ---- echo=FALSE--------------------------------------------------------------
dataset_vignette(
  adtte,
  display_vars = vars(USUBJID, PARAMCD, PARAM, STARTDT, ADT, CNSR)
)

## ---- eval=FALSE--------------------------------------------------------------
#  pd_nact_event <- event_source(
#    dataset_name = "adsl",
#    filter = PDDT < NACTDT | is.na(NACTDT),
#    date = PDDT,
#    set_values_to = vars(
#      EVNTDESC = "Disease Progression prior to NACT",
#      SRCDOM = "ADSL",
#      SRCVAR = "PDDT"
#    )
#  )
#  
#  death_nact_event <- event_source(
#    dataset_name = "adsl",
#    filter = DTHDT < NACTDT | is.na(NACTDT),
#    date = DTHDT,
#    set_values_to = vars(
#      EVNTDESC = "Death prior to NACT",
#      SRCDOM = "ADSL",
#      SRCVAR = "DTHDT"
#    )
#  )
#  
#  lasta_nact_censor <- censor_source(
#    dataset_name = "adsl",
#    date = LASTANDT,
#    set_values_to = vars(
#      EVNTDESC = "Last Tumor Assessment prior to NACT",
#      CNSDTDSC = "Last Tumor Assessment prior to NACT",
#      SRCDOM = "ADSL",
#      SRCVAR = "LASTANDT"
#    )
#  )
#  
#  adtte <- derive_param_tte(
#    dataset_adsl = adsl,
#    start_date = RANDDT,
#    event_conditions = list(pd_nact_event, death_nact_event),
#    censor_conditions = list(lasta_nact_censor, rand_censor),
#    source_datasets = list(adsl = adsl),
#    set_values_to = vars(PARAMCD = "PFSNACT", PARAM = "Progression Free Survival prior to NACT")
#  )

## ----eval=TRUE----------------------------------------------------------------
adtte <- adtte %>%
  derive_vars_duration(
    new_var = AVAL,
    start_date = STARTDT,
    end_date = ADT
  )

## ---- echo=FALSE--------------------------------------------------------------
dataset_vignette(
  adtte
)

## ---- eval=FALSE--------------------------------------------------------------
#  adtte_months <- adtte %>%
#    derive_vars_duration(
#      new_var = AVAL,
#      start_date = STARTDT,
#      end_date = ADT,
#      out_unit = "months"
#    )

## ----eval=TRUE----------------------------------------------------------------
adtte <- adtte %>%
  derive_var_obs_number(
    by_vars = vars(STUDYID, USUBJID),
    order = vars(PARAMCD),
    check_type = "error"
  )

## ---- echo=FALSE--------------------------------------------------------------
dataset_vignette(adtte)

## ----eval=TRUE----------------------------------------------------------------
adtte <- adtte %>%
  derive_vars_merged(
    dataset_add = adsl,
    new_vars = vars(ARMCD, ARM, ACTARMCD, ACTARM, AGE, SEX),
    by_vars = vars(STUDYID, USUBJID)
  )

## ---- echo=FALSE--------------------------------------------------------------
dataset_vignette(
  adtte,
  display_vars = vars(USUBJID, PARAMCD, CNSR, AVAL, ARMCD, AGE, SEX)
)


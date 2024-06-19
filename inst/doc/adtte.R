## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)
library(admiraldev)

## ----warning=FALSE, message=FALSE---------------------------------------------
library(admiral)
library(admiralonco)
library(pharmaverseadam)
library(dplyr)
library(lubridate)

## -----------------------------------------------------------------------------
data("adsl")
data("adrs_onco")
adrs <- adrs_onco

## ----echo=FALSE---------------------------------------------------------------
knitr::kable(list_tte_source_objects(package = "admiralonco"))

## ----eval=FALSE---------------------------------------------------------------
#  adsl_death_event <- event_source(
#    dataset_name = "adsl",
#    date = DTHDT,
#    set_values_to = exprs(
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
    by_vars = get_admiral_option("subject_keys"),
    new_vars = exprs(TEMP_RESPDT = ADT)
  )

## -----------------------------------------------------------------------------
adtte <- derive_param_tte(
  dataset_adsl = adsl,
  start_date = RANDDT,
  event_conditions = list(death_event),
  censor_conditions = list(lastalive_censor, rand_censor),
  source_datasets = list(adsl = adsl, adrs = adrs),
  set_values_to = exprs(PARAMCD = "OS", PARAM = "Overall Survival")
) %>%
  derive_param_tte(
    dataset_adsl = adsl,
    start_date = RANDDT,
    event_conditions = list(pd_event, death_event),
    censor_conditions = list(lasta_censor, rand_censor),
    source_datasets = list(adsl = adsl, adrs = adrs),
    set_values_to = exprs(PARAMCD = "PFS", PARAM = "Progression Free Survival")
  ) %>%
  derive_param_tte(
    dataset_adsl = filter(adsl, !is.na(TEMP_RESPDT)),
    start_date = TEMP_RESPDT,
    event_conditions = list(pd_event, death_event),
    censor_conditions = list(lasta_censor),
    source_datasets = list(adsl = adsl, adrs = adrs),
    set_values_to = exprs(PARAMCD = "RSD", PARAM = "Duration of Response")
  )

## ----echo=FALSE---------------------------------------------------------------
dataset_vignette(
  adtte,
  display_vars = exprs(USUBJID, PARAMCD, PARAM, STARTDT, ADT, CNSR)
)

## ----eval=FALSE---------------------------------------------------------------
#  pd_nact_event <- event_source(
#    dataset_name = "adsl",
#    filter = PDDT < NACTDT | is.na(NACTDT),
#    date = PDDT,
#    set_values_to = exprs(
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
#    set_values_to = exprs(
#      EVNTDESC = "Death prior to NACT",
#      SRCDOM = "ADSL",
#      SRCVAR = "DTHDT"
#    )
#  )
#  
#  lasta_nact_censor <- censor_source(
#    dataset_name = "adsl",
#    date = LASTANDT,
#    set_values_to = exprs(
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
#    set_values_to = exprs(PARAMCD = "PFSNACT", PARAM = "Progression Free Survival prior to NACT")
#  )

## ----eval=TRUE----------------------------------------------------------------
adtte <- adtte %>%
  derive_vars_duration(
    new_var = AVAL,
    start_date = STARTDT,
    end_date = ADT
  )

## ----echo=FALSE---------------------------------------------------------------
dataset_vignette(
  adtte
)

## ----eval=FALSE---------------------------------------------------------------
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
    by_vars = get_admiral_option("subject_keys"),
    order = exprs(PARAMCD),
    check_type = "error"
  )

## ----echo=FALSE---------------------------------------------------------------
dataset_vignette(adtte)

## ----eval=TRUE----------------------------------------------------------------
adtte <- adtte %>%
  derive_vars_merged(
    dataset_add = adsl,
    new_vars = exprs(ARMCD, ARM, ACTARMCD, ACTARM, AGE, SEX),
    by_vars = get_admiral_option("subject_keys")
  )

## ----echo=FALSE---------------------------------------------------------------
dataset_vignette(
  adtte,
  display_vars = exprs(USUBJID, PARAMCD, CNSR, AVAL, ARMCD, AGE, SEX)
)


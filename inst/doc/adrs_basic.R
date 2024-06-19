## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

library(admiraldev)

## ----message=FALSE------------------------------------------------------------
library(admiral)
library(admiralonco)
library(dplyr)
library(pharmaversesdtm)
library(pharmaverseadam)
library(lubridate)
library(stringr)
data("adsl")
data("rs_onco_recist")
data("tu_onco_recist")

rs <- rs_onco_recist
tu <- tu_onco_recist

rs <- convert_blanks_to_na(rs)
tu <- convert_blanks_to_na(tu)

## ----echo=FALSE---------------------------------------------------------------
# select subjects from adsl such that there is one subject without RS data
rs_subjects <- unique(rs$USUBJID)
adsl_subjects <- unique(adsl$USUBJID)
adsl <- filter(
  adsl,
  USUBJID %in% union(rs_subjects, setdiff(adsl_subjects, rs_subjects)[1])
)

## ----eval=TRUE----------------------------------------------------------------
adsl_vars <- exprs(RANDDT)
adrs <- derive_vars_merged(
  rs,
  dataset_add = adsl,
  new_vars = adsl_vars,
  by_vars = get_admiral_option("subject_keys")
)

## ----eval=TRUE, echo=FALSE----------------------------------------------------
dataset_vignette(
  adrs,
  display_vars = exprs(USUBJID, RSTESTCD, RSDTC, VISIT, RANDDT),
  filter = RSTESTCD == "OVRLRESP"
)

## -----------------------------------------------------------------------------
adrs <- adrs %>%
  filter(RSEVAL == "INVESTIGATOR" & RSTESTCD == "OVRLRESP") %>%
  mutate(
    PARAMCD = "OVR",
    PARAM = "Overall Response by Investigator",
    PARCAT1 = "Tumor Response",
    PARCAT2 = "Investigator",
    PARCAT3 = "RECIST 1.1"
  )

## ----echo=FALSE---------------------------------------------------------------
dataset_vignette(
  adrs,
  display_vars = exprs(USUBJID, VISIT, RSTESTCD, RSEVAL, PARAMCD, PARAM, PARCAT1, PARCAT2, PARCAT3)
)

## -----------------------------------------------------------------------------
adrs <- adrs %>%
  derive_vars_dt(
    dtc = RSDTC,
    new_vars_prefix = "A",
    highest_imputation = "D",
    date_imputation = "last"
  ) %>%
  mutate(AVISIT = VISIT)

## ----echo=FALSE---------------------------------------------------------------
dataset_vignette(
  adrs,
  display_vars = exprs(USUBJID, AVISIT, PARAMCD, PARAM, RSSTRESC, RSDTC, ADT, ADTF)
)

## -----------------------------------------------------------------------------
adrs <- adrs %>%
  mutate(
    AVALC = RSSTRESC,
    AVAL = aval_resp(AVALC)
  )

## ----echo=FALSE---------------------------------------------------------------
dataset_vignette(
  adrs,
  display_vars = exprs(USUBJID, AVISIT, PARAMCD, PARAM, RSSTRESC, AVALC, AVAL)
)

## -----------------------------------------------------------------------------
worst_resp <- function(arg) {
  case_when(
    arg == "NE" ~ 1,
    arg == "CR" ~ 2,
    arg == "PR" ~ 3,
    arg == "SD" ~ 4,
    arg == "NON-CR/NON-PD" ~ 5,
    arg == "PD" ~ 6,
    TRUE ~ 0
  )
}

adrs <- adrs %>%
  restrict_derivation(
    derivation = derive_var_extreme_flag,
    args = params(
      by_vars = c(get_admiral_option("subject_keys"), exprs(ADT)),
      order = exprs(worst_resp(AVALC), RSSEQ),
      new_var = ANL01FL,
      mode = "last"
    ),
    filter = !is.na(AVAL) & ADT >= RANDDT
  )

## ----echo=FALSE---------------------------------------------------------------
dataset_vignette(
  adrs,
  display_vars = exprs(USUBJID, AVISIT, PARAMCD, PARAM, AVALC, ADT, RANDDT, ANL01FL)
)

## ----eval=FALSE---------------------------------------------------------------
#  adrs <- adrs %>%
#    mutate(
#      ANL01FL = case_when(
#        !is.na(AVAL) & ADT >= RANDDT & ADT < NACTDT ~ "Y",
#        TRUE ~ NA_character_
#      )
#    )

## ----eval=FALSE---------------------------------------------------------------
#  adrs <- adrs %>%
#    derive_var_relative_flag(
#      by_vars = get_admiral_option("subject_keys"),
#      order = exprs(ADT, RSSEQ),
#      new_var = ANL02FL,
#      condition = AVALC == "PD",
#      mode = "first",
#      selection = "before",
#      inclusive = TRUE
#    )

## -----------------------------------------------------------------------------
adrs <- adrs %>%
  derive_extreme_records(
    dataset_ref = adsl,
    dataset_add = adrs,
    by_vars = get_admiral_option("subject_keys"),
    filter_add = PARAMCD == "OVR" & AVALC == "PD" & ANL01FL == "Y",
    order = exprs(ADT, RSSEQ),
    mode = "first",
    exist_flag = AVALC,
    false_value = "N",
    set_values_to = exprs(
      PARAMCD = "PD",
      PARAM = "Disease Progression by Investigator",
      PARCAT1 = "Tumor Response",
      PARCAT2 = "Investigator",
      PARCAT3 = "RECIST 1.1",
      AVAL = yn_to_numeric(AVALC),
      ANL01FL = "Y"
    )
  )

## ----echo=FALSE---------------------------------------------------------------
dataset_vignette(
  adrs,
  display_vars = exprs(USUBJID, AVISIT, PARAMCD, PARAM, AVALC, ADT, ANL01FL),
  filter = PARAMCD == "PD"
)

## -----------------------------------------------------------------------------
pd <- date_source(
  dataset_name = "adrs",
  date = ADT,
  filter = PARAMCD == "PD" & AVALC == "Y"
)

## ----eval=FALSE---------------------------------------------------------------
#  pd <- date_source(
#    dataset_name = "adsl",
#    date = PDDT
#  )

## -----------------------------------------------------------------------------
adrs <- adrs %>%
  derive_param_response(
    dataset_adsl = adsl,
    filter_source = PARAMCD == "OVR" & AVALC %in% c("CR", "PR") & ANL01FL == "Y",
    source_pd = pd,
    source_datasets = list(adrs = adrs),
    set_values_to = exprs(
      PARAMCD = "RSP",
      PARAM = "Response by Investigator (confirmation not required)",
      PARCAT1 = "Tumor Response",
      PARCAT2 = "Investigator",
      PARCAT3 = "RECIST 1.1",
      AVAL = yn_to_numeric(AVALC),
      ANL01FL = "Y"
    )
  )

## ----echo=FALSE---------------------------------------------------------------
dataset_vignette(
  adrs,
  display_vars = exprs(USUBJID, AVISIT, PARAMCD, PARAM, AVALC, ADT, ANL01FL),
  filter = PARAMCD == "RSP"
)

## -----------------------------------------------------------------------------
resp <- date_source(
  dataset_name = "adrs",
  date = ADT,
  filter = PARAMCD == "RSP" & AVALC == "Y"
)

## -----------------------------------------------------------------------------
adrs <- adrs %>%
  derive_param_clinbenefit(
    dataset_adsl = adsl,
    filter_source = PARAMCD == "OVR" & ANL01FL == "Y",
    source_resp = resp,
    source_pd = pd,
    source_datasets = list(adrs = adrs),
    reference_date = RANDDT,
    ref_start_window = 42,
    set_values_to = exprs(
      PARAMCD = "CB",
      PARAM = "Clinical Benefit by Investigator (confirmation for response not required)",
      PARCAT1 = "Tumor Response",
      PARCAT2 = "Investigator",
      PARCAT3 = "RECIST 1.1",
      AVAL = yn_to_numeric(AVALC),
      ANL01FL = "Y"
    )
  )

## ----echo=FALSE---------------------------------------------------------------
dataset_vignette(
  adrs,
  display_vars = exprs(USUBJID, AVISIT, PARAMCD, PARAM, AVALC, ADT, RANDDT, ANL01FL),
  filter = PARAMCD == "CB"
)

## -----------------------------------------------------------------------------
adrs <- adrs %>%
  derive_param_bor(
    dataset_adsl = adsl,
    filter_source = PARAMCD == "OVR" & ANL01FL == "Y",
    source_pd = pd,
    source_datasets = list(adrs = adrs),
    reference_date = RANDDT,
    ref_start_window = 42,
    set_values_to = exprs(
      PARAMCD = "BOR",
      PARAM = "Best Overall Response by Investigator (confirmation not required)",
      PARCAT1 = "Tumor Response",
      PARCAT2 = "Investigator",
      PARCAT3 = "RECIST 1.1",
      AVAL = aval_resp(AVALC),
      ANL01FL = "Y"
    )
  )

## ----echo=FALSE---------------------------------------------------------------
dataset_vignette(
  adrs,
  display_vars = exprs(USUBJID, AVISIT, PARAMCD, PARAM, AVALC, ADT, RANDDT, ANL01FL),
  filter = PARAMCD == "BOR"
)

## ----eval=FALSE---------------------------------------------------------------
#  aval_resp_new <- function(arg) {
#    case_when(
#      arg == "CR" ~ 7,
#      arg == "PR" ~ 6,
#      arg == "SD" ~ 5,
#      arg == "NON-CR/NON-PD" ~ 4,
#      arg == "PD" ~ 3,
#      arg == "NE" ~ 2,
#      arg == "MISSING" ~ 1,
#      TRUE ~ NA_real_
#    )
#  }

## -----------------------------------------------------------------------------
adrs <- adrs %>%
  derive_extreme_records(
    dataset_ref = adsl,
    dataset_add = adrs,
    by_vars = get_admiral_option("subject_keys"),
    filter_add = PARAMCD == "BOR" & AVALC %in% c("CR", "PR"),
    order = exprs(ADT, RSSEQ),
    mode = "first",
    exist_flag = AVALC,
    false_value = "N",
    set_values_to = exprs(
      PARAMCD = "BCP",
      PARAM = "Best Overall Response of CR/PR by Investigator (confirmation not required)",
      PARCAT1 = "Tumor Response",
      PARCAT2 = "Investigator",
      PARCAT3 = "RECIST 1.1",
      AVAL = yn_to_numeric(AVALC),
      ANL01FL = "Y"
    )
  )

## ----echo=FALSE---------------------------------------------------------------
dataset_vignette(
  adrs,
  display_vars = exprs(USUBJID, AVISIT, PARAMCD, PARAM, AVALC, ADT, ANL01FL),
  filter = PARAMCD == "BCP"
)

## -----------------------------------------------------------------------------
adrs <- adrs %>%
  derive_param_confirmed_resp(
    dataset_adsl = adsl,
    filter_source = PARAMCD == "OVR" & ANL01FL == "Y",
    source_pd = pd,
    source_datasets = list(adrs = adrs),
    ref_confirm = 28,
    set_values_to = exprs(
      PARAMCD = "CRSP",
      PARAM = "Confirmed Response by Investigator",
      PARCAT1 = "Tumor Response",
      PARCAT2 = "Investigator",
      PARCAT3 = "RECIST 1.1",
      AVAL = yn_to_numeric(AVALC),
      ANL01FL = "Y"
    )
  )

confirmed_resp <- date_source(
  dataset_name = "adrs",
  date = ADT,
  filter = PARAMCD == "CRSP" & AVALC == "Y"
)

adrs <- adrs %>%
  derive_param_clinbenefit(
    dataset_adsl = adsl,
    filter_source = PARAMCD == "OVR" & ANL01FL == "Y",
    source_resp = confirmed_resp,
    source_pd = pd,
    source_datasets = list(adrs = adrs),
    reference_date = RANDDT,
    ref_start_window = 42,
    set_values_to = exprs(
      PARAMCD = "CCB",
      PARAM = "Confirmed Clinical Benefit by Investigator",
      PARCAT1 = "Tumor Response",
      PARCAT2 = "Investigator",
      PARCAT3 = "RECIST 1.1",
      AVAL = yn_to_numeric(AVALC),
      ANL01FL = "Y"
    )
  ) %>%
  derive_param_confirmed_bor(
    dataset_adsl = adsl,
    filter_source = PARAMCD == "OVR" & ANL01FL == "Y",
    source_pd = pd,
    source_datasets = list(adrs = adrs),
    reference_date = RANDDT,
    ref_start_window = 42,
    ref_confirm = 28,
    set_values_to = exprs(
      PARAMCD = "CBOR",
      PARAM = "Best Confirmed Overall Response by Investigator",
      PARCAT1 = "Tumor Response",
      PARCAT2 = "Investigator",
      PARCAT3 = "RECIST 1.1",
      AVAL = aval_resp(AVALC),
      ANL01FL = "Y"
    )
  ) %>%
  derive_extreme_records(
    dataset_ref = adsl,
    dataset_add = adrs,
    by_vars = get_admiral_option("subject_keys"),
    filter_add = PARAMCD == "CBOR" & AVALC %in% c("CR", "PR"),
    order = exprs(ADT, RSSEQ),
    mode = "first",
    exist_flag = AVALC,
    false_value = "N",
    set_values_to = exprs(
      PARAMCD = "CBCP",
      PARAM = "Best Confirmed Overall Response of CR/PR by Investigator",
      PARCAT1 = "Tumor Response",
      PARCAT2 = "Investigator",
      PARCAT3 = "RECIST 1.1",
      AVAL = yn_to_numeric(AVALC),
      ANL01FL = "Y"
    )
  )

## ----echo=FALSE---------------------------------------------------------------
dataset_vignette(
  adrs,
  display_vars = exprs(USUBJID, AVISIT, PARAMCD, PARAM, AVALC, ADT, RANDDT, ANL01FL),
  filter = PARAMCD %in% c("CRSP", "CCB", "CBOR", "CBCP")
)

## -----------------------------------------------------------------------------
adrs_bicr <- rs %>%
  filter(
    RSEVAL == "INDEPENDENT ASSESSOR" & RSACPTFL == "Y" & RSTESTCD == "OVRLRESP"
  ) %>%
  mutate(
    PARAMCD = "OVRB",
    PARAM = "Overall Response by BICR",
    PARCAT1 = "Tumor Response",
    PARCAT2 = "Blinded Independent Central Review",
    PARCAT3 = "RECIST 1.1"
  )

## ----echo=FALSE---------------------------------------------------------------
dataset_vignette(
  adrs_bicr,
  display_vars = exprs(USUBJID, VISIT, RSTESTCD, RSEVAL, PARAMCD, PARAM, PARCAT1, PARCAT2, PARCAT3),
  filter = PARAMCD == "OVRR1"
)

## -----------------------------------------------------------------------------
adsldth <- adsl %>%
  select(!!!get_admiral_option("subject_keys"), DTHDT, !!!adsl_vars)

adrs <- adrs %>%
  derive_extreme_records(
    dataset_ref = adsldth,
    dataset_add = adsldth,
    by_vars = get_admiral_option("subject_keys"),
    filter_add = !is.na(DTHDT),
    exist_flag = AVALC,
    false_value = "N",
    set_values_to = exprs(
      PARAMCD = "DEATH",
      PARAM = "Death",
      PARCAT1 = "Reference Event",
      AVAL = yn_to_numeric(AVALC),
      ANL01FL = "Y",
      ADT = DTHDT
    )
  ) %>%
  select(-DTHDT)

## ----echo=FALSE---------------------------------------------------------------
dataset_vignette(
  adrs,
  display_vars = exprs(USUBJID, AVISIT, PARAMCD, PARAM, AVALC, ADT, ANL01FL),
  filter = PARAMCD == "DEATH"
)

## -----------------------------------------------------------------------------
adrs <- adrs %>%
  derive_extreme_records(
    dataset_ref = adsl,
    dataset_add = adrs,
    by_vars = get_admiral_option("subject_keys"),
    filter_add = PARAMCD == "OVR" & ANL01FL == "Y",
    order = exprs(ADT, RSSEQ),
    mode = "last",
    set_values_to = exprs(
      PARAMCD = "LSTA",
      PARAM = "Last Disease Assessment by Investigator",
      PARCAT1 = "Tumor Response",
      PARCAT2 = "Investigator",
      PARCAT3 = "RECIST 1.1",
      ANL01FL = "Y"
    )
  )

## ----echo=FALSE---------------------------------------------------------------
dataset_vignette(
  adrs,
  display_vars = exprs(USUBJID, AVISIT, PARAMCD, PARAM, AVALC, ADT, ANL01FL),
  filter = PARAMCD == "LSTA"
)

## -----------------------------------------------------------------------------
adslmdis <- adsl %>%
  select(!!!get_admiral_option("subject_keys"), !!!adsl_vars)

adrs <- adrs %>%
  derive_param_exist_flag(
    dataset_ref = adslmdis,
    dataset_add = tu,
    condition = TUEVAL == "INVESTIGATOR" & TUSTRESC == "TARGET" & VISIT == "SCREENING",
    false_value = "N",
    missing_value = "N",
    set_values_to = exprs(
      PARAMCD = "MDIS",
      PARAM = "Measurable Disease at Baseline by Investigator",
      PARCAT2 = "Investigator",
      PARCAT3 = "RECIST 1.1",
      AVAL = yn_to_numeric(AVALC),
      ANL01FL = "Y"
    )
  )

## ----echo=FALSE---------------------------------------------------------------
dataset_vignette(
  adrs,
  display_vars = exprs(USUBJID, AVISIT, PARAMCD, PARAM, AVALC, ADT, ANL01FL),
  filter = PARAMCD == "MDIS"
)

## ----eval=TRUE----------------------------------------------------------------
adrs <- adrs %>%
  derive_var_obs_number(
    by_vars = get_admiral_option("subject_keys"),
    order = exprs(PARAMCD, ADT, VISITNUM, RSSEQ),
    check_type = "error"
  )

## ----eval=TRUE, echo=FALSE----------------------------------------------------
dataset_vignette(
  adrs,
  display_vars = exprs(USUBJID, PARAMCD, ADT, VISITNUM, AVISIT, ASEQ),
  filter = USUBJID == "01-701-1015"
)

## ----eval=TRUE----------------------------------------------------------------
adrs <- adrs %>%
  derive_vars_merged(
    dataset_add = select(adsl, !!!negate_vars(adsl_vars)),
    by_vars = get_admiral_option("subject_keys")
  )

## ----eval=TRUE, echo=FALSE----------------------------------------------------
dataset_vignette(
  adrs,
  display_vars = exprs(USUBJID, RFSTDTC, RFENDTC, DTHDTC, DTHFL, AGE, AGEU),
  filter = USUBJID == "01-701-1015"
)

